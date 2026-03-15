#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: TP53 somatic variants database
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_tp53.py
Parse the NCI TP53 database (GRCh38) somatic and germline variant files.

Required files:
  SomaticVariants_GRCh38.csv   (primary)
  GermlineVariants_GRCh38.csv  (optional, if include_germline=True)

Download: https://tp53.cancer.gov

Key column usage:
  - Position  : g_description_GRCh38 (e.g. g.7675184A>G)
  - REF/ALT   : parsed from g_description_GRCh38
  - Consequence: Effect (case-insensitive)
  - Filter    : DNE_LOFclass != 'NA' (excludes variants with no functional evidence)
  - Cancer type: Morphology (histological diagnosis)
  - HGVSc     : c_description
  - HGVSp     : ProtDescription
  - Count     : TCGA_ICGC_GENIE_count if available, else 1 per row
"""

import os
import re

import pandas as pd

from parsers.common import (
    STANDARD_COLS,
    clean_allele,
    empty_standard_df,
    is_valid_allele,
    map_consequence,
    normalise_chrom,
    setup_logger,
)

log = setup_logger('TP53')

# Regex to parse g_description_GRCh38: e.g. "g.7674220C>T"
_GDESC_SNV_RE = re.compile(r'g\.(\d+)([ACGTacgt])>([ACGTacgt])')

# Consequence mapping — case-insensitive lookup applied at parse time
_TP53_CONSEQUENCE_MAP = {
    'missense':   'missense',
    'nonsense':   'nonsense',
    'silent':     'synonymous',
    'splice':     'splice_site',
    'frameshift': 'frameshift',
    'in-frame':   'inframe_indel',
    'in frame':   'inframe_indel',
    'deletion':   'frameshift',
    'insertion':  'frameshift',
    'complex':    'frameshift',
    'intronic':   'intronic',
    'other':      'unknown',
    'unknown':    'unknown',
    '':           'unknown',
}

# DNE_LOFclass values to include (exclude 'NA')
_DNE_LOF_INCLUDE = {'DNE_LOF', 'DNE', 'LOF'}


def parse_tp53(somatic_tsv: str,
               germline_tsv: str | None = None,
               include_germline: bool = False) -> pd.DataFrame:
    frames = []

    if os.path.exists(somatic_tsv):
        log.info('Parsing TP53 somatic file: %s', somatic_tsv)
        df = _parse_tp53_tsv(somatic_tsv, source_label='TP53_somatic')
        if df is not None and not df.empty:
            frames.append(df)
    else:
        log.warning('TP53 somatic file not found: %s', somatic_tsv)

    if include_germline and germline_tsv and os.path.exists(germline_tsv):
        log.info('Parsing TP53 germline file: %s', germline_tsv)
        df = _parse_tp53_tsv(germline_tsv, source_label='TP53_germline')
        if df is not None and not df.empty:
            frames.append(df)

    if not frames:
        return empty_standard_df()

    result = pd.concat(frames, ignore_index=True)
    log.info('TP53: %d total rows', len(result))
    return result


def _parse_tp53_tsv(path: str, source_label: str) -> pd.DataFrame | None:
    sep = ',' if path.endswith('.csv') else '\t'
    try:
        df_raw = pd.read_csv(path, sep=sep, dtype=str, low_memory=False,
                             na_filter=False)
        if len(df_raw.columns) == 1:
            sep = ',' if sep == '\t' else '\t'
            df_raw = pd.read_csv(path, sep=sep, dtype=str, low_memory=False,
                                 na_filter=False)
    except Exception as e:
        log.error('Failed to read TP53 file %s: %s', path, e)
        return None

    log.info('TP53 columns found: %s', list(df_raw.columns))
    df_raw = df_raw.fillna('')

    # Locate required columns
    pos_col      = _find_col(df_raw, ['g_description_GRCh38'])
    effect_col   = _find_col(df_raw, ['Effect', 'Mutation_type', 'mutation_type'])
    dnelof_col   = _find_col(df_raw, ['DNE_LOFclass'])
    morpho_col   = _find_col(df_raw, ['Morphology', 'Short_topo', 'Topography'])
    count_col    = _find_col(df_raw, ['TCGA_ICGC_GENIE_count', 'Somatic_count', 'Count'])
    hgvsp_col    = _find_col(df_raw, ['ProtDescription', 'HGVSp', 'AAchange'])
    hgvsc_col    = _find_col(df_raw, ['c_description', 'HGVSc'])

    if not pos_col:
        log.error('TP53: g_description_GRCh38 column not found — cannot parse positions')
        return None

    n_skipped_dnelof = 0
    n_skipped_pos    = 0
    rows = []

    for _, row in df_raw.iterrows():
        try:
            # --- DNE_LOFclass filter ---
            if dnelof_col:
                dnelof = str(row[dnelof_col]).strip()
                if dnelof not in _DNE_LOF_INCLUDE:
                    n_skipped_dnelof += 1
                    continue

            # --- Position and alleles from g_description_GRCh38 ---
            gdesc = str(row[pos_col]).strip()
            m = _GDESC_SNV_RE.search(gdesc)
            if not m:
                n_skipped_pos += 1
                continue

            pos = int(m.group(1))
            ref = m.group(2).upper()
            alt = m.group(3).upper()

            if not is_valid_allele(ref) or not is_valid_allele(alt):
                continue

            # TP53 is always chr17
            chrom = 'chr17'

            # --- Consequence (case-insensitive) ---
            effect_raw  = str(row[effect_col]).strip().lower() if effect_col else ''
            consequence = _TP53_CONSEQUENCE_MAP.get(effect_raw,
                          map_consequence(effect_raw))

            # --- Cancer type from Morphology ---
            cancer_type = str(row[morpho_col]).strip().lower() if morpho_col else ''
            if not cancer_type or cancer_type in ('', 'na', 'nan'):
                cancer_type = 'unspecified'

            # --- Sample count ---
            n_samples = 1
            if count_col:
                count_val = str(row[count_col]).strip()
                if count_val.isdigit():
                    n_samples = int(count_val)

            # --- HGVSc and HGVSp ---
            hgvsc = str(row[hgvsc_col]).strip() if hgvsc_col else ''
            hgvsp = str(row[hgvsp_col]).strip() if hgvsp_col else ''

            rows.append({
                'chrom':       chrom,
                'pos':         pos,
                'ref':         ref,
                'alt':         alt,
                'gene':        'TP53',
                'hgvsc':       hgvsc,
                'hgvsp':       hgvsp,
                'consequence': consequence,
                'cancer_type': cancer_type,
                'n_samples':   n_samples,
                'source':      source_label,
            })
        except (ValueError, TypeError, KeyError):
            continue

    log.info(
        'TP53 (%s): %d rows kept | %d excluded by DNE_LOFclass | %d excluded by unparseable position',
        source_label, len(rows), n_skipped_dnelof, n_skipped_pos
    )

    if not rows:
        log.warning('TP53: no rows produced from %s', path)
        return None

    df = pd.DataFrame(rows, columns=STANDARD_COLS)
    return df


def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None
