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
Parse the IARC TP53 database (GRCh38) somatic and germline variant files.

Required files:
  SomaticVariants_GRCh38.tsv   (primary)
  GermlineVariants_GRCh38.tsv  (optional, if include_germline=True)

Download:
  https://tp53.isb-cgc.org/  (ISB-CGC mirror — recommended, TSV downloads)
  or
  https://tp53.fr/           (IARC primary — may require navigation to Download)

The ISB-CGC mirror provides direct TSV downloads. Columns vary slightly
between releases. This parser handles both the ISB-CGC and IARC TSV formats
and will log which columns it finds.

Key column mappings attempted (in priority order):
  - Genomic position : 'g_description_GRCh38', 'Chr38', 'Chr_GRCh38'
  - REF/ALT          : 'WTcodon', 'Mutcodon', or inferred from g_description
  - Consequence      : 'Effect', 'Mutation_type'
  - Cancer type      : 'Topography', 'Primary_site', 'cancer_type'
  - Sample ID        : 'Sample_ID', 'sample_id', 'Sample_Name'
  - Occurrence count : 'Somatic_count', 'Count', 'somatic_count'
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

# Regex to parse g.description_GRCh38: e.g. "g.7674220C>T" or "g.7673802_7673803del"
_GDESC_SNV_RE  = re.compile(r'g\.(\d+)([ACGTacgt])>([ACGTacgt])')
_GDESC_CHR_RE  = re.compile(r'^chr(\w+):')

# Consequence mapping for TP53 database-specific terms
_TP53_CONSEQUENCE_MAP = {
    'Missense':                  'missense',
    'Nonsense':                  'nonsense',
    'Silent':                    'synonymous',
    'Splice':                    'splice_site',
    'Frameshift':                'frameshift',
    'In-frame':                  'inframe_indel',
    'In frame':                  'inframe_indel',
    'Deletion':                  'frameshift',
    'Insertion':                 'frameshift',
    'Complex':                   'frameshift',
    'Intronic':                  'intronic',
    'Other':                     'unknown',
    'Unknown':                   'unknown',
    '':                          'unknown',
}


def parse_tp53(somatic_tsv: str,
               germline_tsv: str | None = None,
               include_germline: bool = False) -> pd.DataFrame:
    """
    Parse IARC TP53 database TSV files.

    Parameters
    ----------
    somatic_tsv      : path to SomaticVariants_GRCh38.tsv
    germline_tsv     : path to GermlineVariants_GRCh38.tsv (optional)
    include_germline : whether to include germline variants

    Returns
    -------
    DataFrame with STANDARD_COLS schema.
    """
    frames = []

    if os.path.exists(somatic_tsv):
        log.info('Parsing TP53 somatic TSV: %s', somatic_tsv)
        df = _parse_tp53_tsv(somatic_tsv, source_label='TP53_somatic')
        if df is not None and not df.empty:
            frames.append(df)
    else:
        log.warning('TP53 somatic TSV not found: %s', somatic_tsv)

    if include_germline and germline_tsv and os.path.exists(germline_tsv):
        log.info('Parsing TP53 germline TSV: %s', germline_tsv)
        df = _parse_tp53_tsv(germline_tsv, source_label='TP53_germline')
        if df is not None and not df.empty:
            frames.append(df)

    if not frames:
        return empty_standard_df()

    result = pd.concat(frames, ignore_index=True)
    log.info('TP53: %d total rows', len(result))
    return result


def _parse_tp53_tsv(path: str, source_label: str) -> pd.DataFrame | None:
    # Auto-detect delimiter: CSV from tp53.cancer.gov or TSV from older releases
    sep = ',' if path.endswith('.csv') else '\t'
    try:
        df_raw = pd.read_csv(path, sep=sep, dtype=str, low_memory=False,
                             na_filter=False)
        # If only one column was parsed, the delimiter guess was wrong — retry
        if len(df_raw.columns) == 1:
            sep = ',' if sep == '\t' else '\t'
            df_raw = pd.read_csv(path, sep=sep, dtype=str, low_memory=False,
                                 na_filter=False)
    except Exception as e:
        log.error('Failed to read TP53 file %s: %s', path, e)
        return None

    log.debug('TP53 columns in %s: %s', path, list(df_raw.columns))
    df_raw = df_raw.fillna('')

    # Detect column names for this release
    chrom_col  = _find_col(df_raw, ['Chr38', 'Chr_GRCh38', 'chr38'])
    pos_col    = _find_col(df_raw, ['g_description_GRCh38', 'Genome_GRCh38_position'])
    effect_col = _find_col(df_raw, ['Effect', 'Mutation_type', 'mutation_type'])
    ctype_col  = _find_col(df_raw, ['Topography', 'Primary_site', 'cancer_type', 'primary_site'])
    sample_col = _find_col(df_raw, ['Sample_ID', 'sample_id', 'Sample_Name'])
    count_col  = _find_col(df_raw, ['Somatic_count', 'Count', 'somatic_count', 'count'])
    hgvsp_col  = _find_col(df_raw, ['c_description', 'HGVSp', 'AAchange'])
    hgvsc_col  = _find_col(df_raw, ['c_description', 'HGVSc'])
    ref_col    = _find_col(df_raw, ['WTbase', 'Ref_base', 'ref'])
    alt_col    = _find_col(df_raw, ['Mutbase', 'Alt_base', 'alt'])

    rows = []
    for _, row in df_raw.iterrows():
        try:
            # --- Chromosome ---
            chrom = ''
            if chrom_col:
                chrom = normalise_chrom(row[chrom_col])

            # --- Position and alleles ---
            pos = None
            ref = ''
            alt = ''

            gdesc = str(row.get(pos_col or '', '')).strip() if pos_col else ''

            # Try to parse genomic description (e.g. "g.7674220C>T")
            m = _GDESC_SNV_RE.search(gdesc)
            if m:
                pos = int(m.group(1))
                ref = m.group(2).upper()
                alt = m.group(3).upper()
            elif ref_col and alt_col:
                # Use separate columns if present
                ref = clean_allele(row[ref_col])
                alt = clean_allele(row[alt_col])
                # Position from numeric column if gdesc unparseable
                if pos_col and row[pos_col].isdigit():
                    pos = int(row[pos_col])

            if pos is None or not ref or not alt:
                continue
            if not is_valid_allele(ref) or not is_valid_allele(alt):
                continue

            if not chrom:
                # TP53 is always chr17
                chrom = 'chr17'

            # --- Consequence ---
            effect_raw  = str(row[effect_col]).strip() if effect_col else ''
            consequence = _TP53_CONSEQUENCE_MAP.get(effect_raw,
                          map_consequence(effect_raw))

            # --- Cancer type ---
            cancer_type = str(row[ctype_col]).strip().lower() if ctype_col else 'unspecified'
            if not cancer_type:
                cancer_type = 'unspecified'

            # --- Sample count ---
            n_samples = 1
            if count_col and row[count_col].isdigit():
                n_samples = int(row[count_col])
            elif sample_col and row[sample_col]:
                n_samples = 1

            hgvsp = str(row[hgvsp_col]).strip() if hgvsp_col else ''
            hgvsc = str(row[hgvsc_col]).strip() if hgvsc_col else ''

            rows.append({
                'chrom':        chrom,
                'pos':          pos,
                'ref':          ref,
                'alt':          alt,
                'gene':         'TP53',
                'hgvsc':        hgvsc,
                'hgvsp':        hgvsp,
                'consequence':  consequence,
                'cancer_type':  cancer_type,
                'n_samples':    n_samples,
                'source':       source_label,
            })
        except (ValueError, TypeError, KeyError):
            continue

    if not rows:
        log.warning('TP53: no rows produced from %s', path)
        return None

    df = pd.DataFrame(rows, columns=STANDARD_COLS)
    log.info('TP53 (%s): %d rows', source_label, len(df))
    return df


def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Return the first candidate column name present in df, or None."""
    for c in candidates:
        if c in df.columns:
            return c
    return None
