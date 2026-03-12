#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: AACR GENIE MAF
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_genie.py
Parse AACR Project GENIE data_mutations_extended.txt (MAF format, GRCh38).

Required file:
  data_mutations_extended.txt
  Download: https://www.synapse.org/#!Synapse:syn7222066
  Requires a free Synapse account and GENIE data access agreement.
  Select the latest public release (e.g. GENIE 15.0-public).

The MAF file follows the GDC MAF format specification. Key columns used:
  Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2,
  Hugo_Symbol, HGVSc, HGVSp_Short, Variant_Classification,
  Tumor_Sample_Barcode, Center (maps to cancer centre / panel)

GENIE does not annotate cancer type per mutation row directly. Cancer type
is retrieved from the companion clinical file: data_clinical_sample.txt
using SAMPLE_ID -> ONCOTREE_CODE.
"""

import os

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

log = setup_logger('GENIE')

_MAF_COLS = {
    'chrom':       'Chromosome',
    'pos':         'Start_Position',
    'ref':         'Reference_Allele',
    'alt':         'Tumor_Seq_Allele2',
    'gene':        'Hugo_Symbol',
    'hgvsc':       'HGVSc',
    'hgvsp':       'HGVSp_Short',
    'consequence': 'Variant_Classification',
    'sample_id':   'Tumor_Sample_Barcode',
}

_CLINICAL_SAMPLE_ID_COL    = 'SAMPLE_ID'
_CLINICAL_ONCOTREE_COL     = 'ONCOTREE_CODE'
_CLINICAL_CANCER_TYPE_COL  = 'CANCER_TYPE'


def parse_genie(maf_path: str,
                clinical_sample_path: str | None = None) -> pd.DataFrame:
    """
    Parse GENIE MAF file.

    Parameters
    ----------
    maf_path              : path to data_mutations_extended.txt
    clinical_sample_path  : path to data_clinical_sample.txt (optional;
                            enables cancer type annotation per sample)

    Returns
    -------
    DataFrame with STANDARD_COLS schema.
    One row per mutation per sample. Aggregation to unique samples per variant
    is performed at the merge step.
    """
    if not os.path.exists(maf_path):
        log.warning('GENIE MAF not found: %s  — skipping', maf_path)
        return empty_standard_df()

    log.info('Parsing GENIE MAF: %s', maf_path)

    # Build cancer type lookup from clinical file if available
    cancer_type_map: dict[str, str] = {}
    if clinical_sample_path and os.path.exists(clinical_sample_path):
        cancer_type_map = _load_clinical_sample(clinical_sample_path)
        log.info('GENIE clinical lookup: %d samples', len(cancer_type_map))

    # Read MAF; skip metadata lines starting with '#'
    try:
        df_maf = pd.read_csv(
            maf_path, sep='\t', comment='#', dtype=str,
            low_memory=False, na_filter=False,
        )
    except Exception as e:
        log.error('Failed to read GENIE MAF: %s', e)
        return empty_standard_df()

    # Validate required columns
    required = list(_MAF_COLS.values())
    missing  = [c for c in required if c not in df_maf.columns]
    if missing:
        log.error('GENIE MAF missing columns: %s', missing)
        return empty_standard_df()

    rows = []
    for _, row in df_maf.iterrows():
        try:
            chrom   = normalise_chrom(row[_MAF_COLS['chrom']])
            pos_str = row[_MAF_COLS['pos']]
            ref     = clean_allele(row[_MAF_COLS['ref']])
            alt     = clean_allele(row[_MAF_COLS['alt']])

            if not is_valid_allele(ref) or not is_valid_allele(alt):
                continue
            if alt in ('', '-', '.'):
                continue

            pos     = int(float(pos_str))
            gene    = str(row[_MAF_COLS['gene']]).strip()
            hgvsc   = str(row[_MAF_COLS['hgvsc']]).strip()
            hgvsp   = str(row[_MAF_COLS['hgvsp']]).strip()
            var_cls = str(row[_MAF_COLS['consequence']]).strip()
            sample  = str(row[_MAF_COLS['sample_id']]).strip()

            consequence = map_consequence(var_cls)

            # Cancer type: from clinical lookup or fall back to 'unspecified'
            cancer_type = cancer_type_map.get(sample, 'unspecified').lower()

            rows.append({
                'chrom':        chrom,
                'pos':          pos,
                'ref':          ref,
                'alt':          alt,
                'gene':         gene,
                'hgvsc':        hgvsc,
                'hgvsp':        hgvsp,
                'consequence':  consequence,
                'cancer_type':  cancer_type,
                'n_samples':    1,
                'source':       'GENIE',
            })
        except (ValueError, TypeError):
            continue

    if not rows:
        log.warning('GENIE: no rows produced')
        return empty_standard_df()

    df = pd.DataFrame(rows, columns=STANDARD_COLS)
    log.info('GENIE: %d rows', len(df))
    return df


def _load_clinical_sample(clinical_path: str) -> dict[str, str]:
    """
    Load GENIE data_clinical_sample.txt and return a dict of
    SAMPLE_ID -> cancer type label (CANCER_TYPE preferred, else ONCOTREE_CODE).
    """
    try:
        df = pd.read_csv(
            clinical_path, sep='\t', comment='#',
            dtype=str, low_memory=False, na_filter=False,
        )
    except Exception as e:
        log.warning('Could not load clinical sample file: %s', e)
        return {}

    lookup: dict[str, str] = {}
    for _, row in df.iterrows():
        sid = str(row.get(_CLINICAL_SAMPLE_ID_COL, '')).strip()
        if not sid:
            continue
        ctype = str(row.get(_CLINICAL_CANCER_TYPE_COL,
                            row.get(_CLINICAL_ONCOTREE_COL, 'unspecified'))).strip()
        lookup[sid] = ctype.lower()
    return lookup
