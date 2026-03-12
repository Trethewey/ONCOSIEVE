#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: DoCM API
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_docm.py
Parse the Database of Curated Mutations (DoCM).

DoCM curates high-confidence cancer driver mutations from the literature.
All variants have published functional evidence.

Data source:
  http://docm.info/api/v1/variants.tsv   (public, no authentication)
  The API returns GRCh37 coordinates. This parser requests GRCh38 if available
  or flags variants as requiring liftover.

Alternatively, provide a locally downloaded TSV.

NOTE: As of 2024, DoCM is no longer actively updated (last release ~2016)
but remains a useful curated set of high-confidence driver mutations.
All DoCM variants are included regardless of sample count (treated as
expert-curated evidence). They are promoted to at least Tier 2.
"""

import os
import time
from io import StringIO
from typing import Optional

import pandas as pd
import requests

from parsers.common import (
    STANDARD_COLS,
    clean_allele,
    empty_standard_df,
    is_valid_allele,
    map_consequence,
    normalise_chrom,
    setup_logger,
)

log = setup_logger('DoCM')

_API_URL = 'http://docm.info/api/v1/variants.tsv'
_TIMEOUT = 60


def parse_docm(tsv: Optional[str] = None,
               use_api: bool = True) -> pd.DataFrame:
    """
    Parse DoCM variants.

    Parameters
    ----------
    tsv     : path to locally downloaded DoCM TSV
    use_api : if True, fetch from API; falls back to tsv if API fails

    Returns
    -------
    DataFrame with STANDARD_COLS schema.
    Coordinates in DoCM are hg19 — flagged in 'source' column.
    Run bcftools liftover or CrossMap post-pipeline if mixing with hg38.
    """
    raw_text = None

    if use_api:
        log.info('Fetching DoCM from API: %s', _API_URL)
        try:
            resp = requests.get(_API_URL, timeout=_TIMEOUT)
            resp.raise_for_status()
            raw_text = resp.text
        except Exception as e:
            log.warning('DoCM API failed: %s — trying local file', e)

    if raw_text is None:
        if not tsv or not os.path.exists(tsv):
            log.warning('DoCM TSV not found and API unavailable — skipping')
            return empty_standard_df()
        log.info('Reading DoCM TSV: %s', tsv)
        with open(tsv) as fh:
            raw_text = fh.read()

    try:
        df_raw = pd.read_csv(StringIO(raw_text), sep='\t', dtype=str,
                             low_memory=False, na_filter=False)
    except Exception as e:
        log.error('Failed to parse DoCM TSV: %s', e)
        return empty_standard_df()

    log.debug('DoCM columns: %s', list(df_raw.columns))
    df_raw = df_raw.fillna('')

    # Column detection
    chrom_col  = _find_col(df_raw, ['hg19_chr', 'chr', 'chromosome'])
    pos_col    = _find_col(df_raw, ['hg19_start', 'start', 'position'])
    ref_col    = _find_col(df_raw, ['ref', 'reference', 'reference_allele'])
    alt_col    = _find_col(df_raw, ['alt', 'variant', 'alt_allele'])
    gene_col   = _find_col(df_raw, ['gene', 'Hugo_Symbol'])
    aa_col     = _find_col(df_raw, ['aa_change', 'amino_acid_change', 'HGVSp'])
    type_col   = _find_col(df_raw, ['variant_type', 'type', 'Variant_Classification'])
    ctype_col  = _find_col(df_raw, ['diseases', 'disease', 'cancer_type', 'Cancer_Type'])
    pmid_col   = _find_col(df_raw, ['pubmed_sources', 'pmid', 'source'])

    rows = []
    for _, row in df_raw.iterrows():
        try:
            chrom = normalise_chrom(row[chrom_col]) if chrom_col else ''
            pos   = int(row[pos_col]) if pos_col and row[pos_col].isdigit() else None
            ref   = clean_allele(row[ref_col]) if ref_col else ''
            alt   = clean_allele(row[alt_col]) if alt_col else ''

            if pos is None or not ref or not alt:
                continue
            if not is_valid_allele(ref) or not is_valid_allele(alt):
                continue

            gene        = str(row[gene_col]).strip() if gene_col else ''
            hgvsp       = str(row[aa_col]).strip() if aa_col else ''
            var_type    = str(row[type_col]).strip() if type_col else ''
            cancer_type = str(row[ctype_col]).strip().lower() if ctype_col else 'unspecified'
            consequence = map_consequence(var_type)

            rows.append({
                'chrom':        chrom,
                'pos':          pos,
                'ref':          ref,
                'alt':          alt,
                'gene':         gene,
                'hgvsc':        '',
                'hgvsp':        hgvsp,
                'consequence':  consequence,
                'cancer_type':  cancer_type,
                'n_samples':    1,   # DoCM = expert-curated; count not meaningful
                'source':       'DoCM_hg19',   # FLAG: coordinates need liftover
            })
        except (ValueError, TypeError, KeyError):
            continue

    if not rows:
        log.warning('DoCM: no rows produced')
        return empty_standard_df()

    df = pd.DataFrame(rows, columns=STANDARD_COLS)
    log.warning('DoCM: coordinates are hg19 — liftover required before '
                'merging with hg38 sources. Run run_liftover.sh first.')
    log.info('DoCM: %d rows', len(df))
    return df


def _find_col(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None
