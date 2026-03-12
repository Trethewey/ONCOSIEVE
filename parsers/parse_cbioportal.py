#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: cBioPortal REST API
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_cbioportal.py
Pull mutation data from the cBioPortal public REST API (GRCh38 studies).

No data download required. Data is fetched via:
  https://www.cbioportal.org/api

This parser:
  - Enumerates all public studies (or TCGA PanCancer studies only)
  - For each study, retrieves all somatic mutations via /mutations endpoint
  - Maps Variant_Classification -> standardised consequence
  - Returns one row per (variant, cancer_type, study) with n_samples = 1
    (samples are counted per-variant at the merge step)

Rate limiting: the cBioPortal API is public but rate-limited.
Set config.cbioportal.request_delay_s >= 0.5 to avoid 429 responses.

NOTE: cBioPortal returns hg19 coordinates for older studies and hg38 for newer
ones. The 'referenceGenomeVersion' field in each study is checked; hg19 studies
are flagged and can be optionally lifted over (liftover not implemented here —
hg19 studies are skipped by default to avoid coordinate mixing).
Set SKIP_HG19 = False below if you want to include them (requires liftover).
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

log = setup_logger('cBioPortal')

API_BASE    = 'https://www.cbioportal.org/api'
SKIP_HG19   = True       # skip hg19 studies to avoid coordinate mixing
BATCH_SIZE  = 10_000     # max mutations per API request (cBioPortal default 10k)
SESSION_TIMEOUT = 30     # seconds


def _get(url: str, params: dict | None = None, delay: float = 0.5) -> Optional[dict | list]:
    """GET with simple retry logic."""
    for attempt in range(3):
        try:
            resp = requests.get(url, params=params, timeout=SESSION_TIMEOUT)
            if resp.status_code == 200:
                time.sleep(delay)
                return resp.json()
            if resp.status_code == 429:
                wait = int(resp.headers.get('Retry-After', 60))
                log.warning('Rate limited; waiting %ds', wait)
                time.sleep(wait)
            else:
                log.warning('HTTP %d for %s', resp.status_code, url)
                break
        except requests.RequestException as e:
            log.warning('Request error (attempt %d): %s', attempt + 1, e)
            time.sleep(5 * (attempt + 1))
    return None


def _get_studies(api_base: str, tcga_only: bool, delay: float) -> list[dict]:
    """Return list of public cBioPortal studies."""
    data = _get(f'{api_base}/studies', params={'projection': 'SUMMARY'}, delay=delay)
    if not data:
        return []
    studies = [s for s in data if s.get('publicStudy', False)]
    if tcga_only:
        studies = [s for s in studies if 'tcga' in s['studyId'].lower()]
    log.info('cBioPortal: %d studies to process', len(studies))
    return studies


def _get_mutations_for_study(study_id: str, api_base: str,
                              delay: float) -> list[dict]:
    """Fetch all mutations for a given study ID."""
    url = f'{api_base}/studies/{study_id}/mutations'
    mutations = []
    page = 0
    while True:
        params = {
            'projection': 'DETAILED',
            'pageSize':   BATCH_SIZE,
            'pageNumber': page,
        }
        data = _get(url, params=params, delay=delay)
        if not data:
            break
        mutations.extend(data)
        if len(data) < BATCH_SIZE:
            break
        page += 1
    return mutations


def parse_cbioportal(api_base: str = API_BASE,
                     use_tcga_pancancer_only: bool = False,
                     max_studies: int = 0,
                     request_delay_s: float = 0.5) -> pd.DataFrame:
    """
    Pull mutation data from cBioPortal public API.

    Parameters
    ----------
    api_base                  : cBioPortal API base URL
    use_tcga_pancancer_only   : restrict to TCGA studies only
    max_studies               : cap on number of studies (0 = no limit)
    request_delay_s           : delay between API calls in seconds

    Returns
    -------
    DataFrame with STANDARD_COLS schema.
    """
    studies = _get_studies(api_base, use_tcga_pancancer_only, request_delay_s)
    if not studies:
        log.error('cBioPortal: failed to retrieve study list')
        return empty_standard_df()

    if max_studies > 0:
        studies = studies[:max_studies]

    all_rows: list[dict] = []

    for i, study in enumerate(studies):
        study_id   = study['studyId']
        cancer_type = study.get('cancerType', {}).get('name', study_id).lower()
        ref_genome  = study.get('referenceGenome', 'hg19').lower()

        if SKIP_HG19 and ref_genome in ('hg19', 'grch37', 'grch37/hg19'):
            log.debug('Skipping hg19 study: %s', study_id)
            continue

        log.info('[%d/%d] Fetching %s (%s)', i + 1, len(studies),
                 study_id, ref_genome)

        mutations = _get_mutations_for_study(study_id, api_base, request_delay_s)
        if not mutations:
            continue

        for mut in mutations:
            try:
                # Genomic coordinates
                chrom = normalise_chrom(str(mut.get('chr', '') or ''))
                pos   = mut.get('startPosition')
                ref   = clean_allele(str(mut.get('referenceAllele', '') or ''))
                alt   = clean_allele(str(mut.get('variantAllele', '') or ''))

                if not chrom or pos is None:
                    continue
                if not is_valid_allele(ref) or not is_valid_allele(alt):
                    continue

                gene        = str(mut.get('gene', {}).get('hugoGeneSymbol', '') or '')
                hgvsc       = str(mut.get('annotationSummary', {}).get(
                                  'transcriptConsequenceSummary', {}).get(
                                  'hgvscShort', '') or '')
                hgvsp       = str(mut.get('annotationSummary', {}).get(
                                  'transcriptConsequenceSummary', {}).get(
                                  'hgvspShort', '') or '')
                var_class   = str(mut.get('mutationType', '') or '')
                sample_id   = str(mut.get('uniqueSampleKey', '') or '')
                consequence = map_consequence(var_class)

                all_rows.append({
                    'chrom':        chrom,
                    'pos':          int(pos),
                    'ref':          ref,
                    'alt':          alt,
                    'gene':         gene,
                    'hgvsc':        hgvsc,
                    'hgvsp':        hgvsp,
                    'consequence':  consequence,
                    'cancer_type':  cancer_type,
                    'n_samples':    1,
                    'source':       'cBioPortal',
                })

            except (KeyError, TypeError, ValueError):
                continue

        log.info('  -> %d mutations added (total so far: %d)',
                 len(mutations), len(all_rows))

    if not all_rows:
        log.warning('cBioPortal: no rows produced')
        return empty_standard_df()

    df = pd.DataFrame(all_rows, columns=STANDARD_COLS)
    log.info('cBioPortal: %d raw mutation rows', len(df))
    return df
