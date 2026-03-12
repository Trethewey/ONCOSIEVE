#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: Cancer Hotspots API
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_hotspots.py
Parse MSK Cancer Hotspots v2 data (GRCh38).

Data source:
  https://www.cancerhotspots.org/files/hotspots_v2.xls   (local file)
  https://www.cancerhotspots.org/api/hotspots/single      (API — used if use_api=True)

Hotspot data provides statistically significant recurrence positions.
Each hotspot position may span multiple amino acid residues. This parser
produces one row per (gene, residue, variant type) combination.

Genomic coordinates in hotspots_v2.xls are hg19. The API returns hg38.
Use use_api=True (default) for GRCh38 coordinates.
"""

import time
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

log = setup_logger('CancerHotspots')

_API_BASE = 'https://www.cancerhotspots.org/api'
_TIMEOUT  = 30


def _fetch_hotspots_api(api_base: str) -> list[dict]:
    """Fetch all single-residue hotspots from the Cancer Hotspots API."""
    url = f'{api_base}/hotspots/single'
    try:
        resp = requests.get(url, timeout=_TIMEOUT)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        log.error('Cancer Hotspots API error: %s', e)
        return []


def parse_hotspots(xls: Optional[str] = None,
                   api_base: str = _API_BASE,
                   use_api: bool = True,
                   max_qvalue: float = 0.05) -> pd.DataFrame:
    """
    Parse Cancer Hotspots v2.

    Parameters
    ----------
    xls        : path to hotspots_v2.xls (used only if use_api=False)
    api_base   : Cancer Hotspots API base URL
    use_api    : if True, pull from API (returns GRCh38 coords)
    max_qvalue : q-value significance threshold

    Returns
    -------
    DataFrame with STANDARD_COLS schema.
    chrom/pos/ref/alt may be empty if the API does not return genomic coords;
    gene + HGVSp are used as fallback join keys at the merge step.
    """
    if use_api:
        log.info('Fetching Cancer Hotspots from API: %s', api_base)
        data = _fetch_hotspots_api(api_base)
    else:
        if not xls or not __import__('os').path.exists(xls):
            log.warning('Cancer Hotspots XLS not found: %s', xls)
            return empty_standard_df()
        log.info('Reading Cancer Hotspots XLS: %s', xls)
        try:
            df_xls = pd.read_excel(xls, dtype=str)
            data = df_xls.to_dict(orient='records')
        except Exception as e:
            log.error('Failed to read hotspots XLS: %s', e)
            return empty_standard_df()

    rows = []
    for hs in data:
        try:
            gene       = str(hs.get('gene', '') or '').strip()
            residue    = str(hs.get('residue', '') or '').strip()
            var_type   = str(hs.get('type', '') or '').strip().lower()
            qvalue     = float(hs.get('qValue', 1.0) or 1.0)
            n_samples  = int(hs.get('tumorCount', 0) or 0)

            if qvalue > max_qvalue:
                continue

            # Genomic coordinates (available in API response as genomicLocation)
            chrom = ''
            pos   = pd.NA
            ref   = ''
            alt   = ''
            genomic = hs.get('genomicLocation', {}) or {}
            if genomic:
                chrom = normalise_chrom(str(genomic.get('chromosome', '') or ''))
                pos   = genomic.get('start') or pd.NA
                ref   = clean_allele(str(genomic.get('referenceAllele', '') or ''))
                alt   = clean_allele(str(genomic.get('variantAllele', '') or ''))

            # Consequence from variant type
            consequence = _map_hs_type(var_type)

            # HGVSp-like from residue
            hgvsp = f'p.{residue}' if residue else ''

            rows.append({
                'chrom':        chrom,
                'pos':          int(pos) if pd.notna(pos) else pd.NA,
                'ref':          ref,
                'alt':          alt,
                'gene':         gene,
                'hgvsc':        '',
                'hgvsp':        hgvsp,
                'consequence':  consequence,
                'cancer_type':  'pan_cancer',
                'n_samples':    n_samples,
                'source':       'CancerHotspots',
            })
        except (ValueError, TypeError, AttributeError):
            continue

    if not rows:
        log.warning('CancerHotspots: no rows produced')
        return empty_standard_df()

    df = pd.DataFrame(rows, columns=STANDARD_COLS)
    log.info('CancerHotspots: %d hotspot entries', len(df))
    return df


def _map_hs_type(var_type: str) -> str:
    mapping = {
        'missense': 'missense',
        'truncating': 'nonsense',
        'indel': 'inframe_indel',
        'splice': 'splice_site',
        'silent': 'synonymous',
        'other': 'unknown',
    }
    return mapping.get(var_type.lower(), 'unknown')
