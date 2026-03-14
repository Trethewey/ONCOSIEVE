#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Standalone cBioPortal snapshot fetcher
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
fetch_cbioportal.py
Fetch all hg38 mutation data from the cBioPortal public REST API and save
as a versioned snapshot TSV.gz.

Run this script once to produce a static data file. Subsequent pipeline runs
load from the snapshot via --from-intermediates, eliminating live API
dependency and 503/rate-limit variability.

Usage:
    python tools/fetch_cbioportal.py
    python tools/fetch_cbioportal.py --output data/cbioportal/cbioportal_YYYYMMDD.tsv.gz
    python tools/fetch_cbioportal.py --max-studies 50   # test run
    python tools/fetch_cbioportal.py --max-retries 10 --retry-delay 30

Output:
    A TSV.gz file in STANDARD_COLS schema, compatible with the ONCOSIEVE
    intermediate format. Copy or symlink to intermediate/cbioportal.tsv.gz
    to use with --from-intermediates.

Retry behaviour:
    503 (Service Unavailable): retries up to --max-retries times with
        exponential backoff starting at --retry-delay seconds.
    429 (Rate Limited):        honours Retry-After header, then retries.
    Network errors:            retries up to --max-retries times.
    404 (No profile):          skips the study silently (expected for some).
    Other 4xx/5xx:             logs a warning and skips the study.
"""

import argparse
import gzip
import os
import sys
import time
from datetime import datetime

import pandas as pd
import requests

# Add project root to path so parsers.common is importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from parsers.common import (
    STANDARD_COLS,
    clean_allele,
    empty_standard_df,
    is_valid_allele,
    map_consequence,
    normalise_chrom,
    setup_logger,
)

log = setup_logger('fetch_cbioportal')

API_BASE        = 'https://www.cbioportal.org/api'
BATCH_SIZE      = 10_000
SESSION_TIMEOUT = 60     # longer than parser default — large studies can be slow
SKIP_HG19       = True


# =============================================================================
# HTTP helpers
# =============================================================================

def _post_with_retry(url: str,
                     params: dict,
                     body: dict,
                     delay: float,
                     max_retries: int,
                     retry_delay: float) -> list | None:
    """
    POST with retry logic for 503, 429, and network errors.

    Returns parsed JSON list on success, None if all retries exhausted
    or a non-retryable error is encountered.
    """
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.post(
                url,
                params=params,
                json=body,
                timeout=SESSION_TIMEOUT,
            )
            time.sleep(delay)

            if resp.status_code == 200:
                return resp.json()

            if resp.status_code == 404:
                # No mutation profile for this study — normal, do not retry.
                return None

            if resp.status_code == 429:
                wait = int(resp.headers.get('Retry-After', 60))
                log.warning('Rate limited (429); waiting %ds before retry', wait)
                time.sleep(wait)
                continue

            if resp.status_code == 503:
                wait = retry_delay * (2 ** (attempt - 1))  # exponential backoff
                log.warning(
                    'HTTP 503 for %s (attempt %d/%d) — retrying in %.0fs',
                    url, attempt, max_retries, wait
                )
                time.sleep(wait)
                continue

            # Other error — log and give up on this study.
            log.warning('HTTP %d for %s — skipping study', resp.status_code, url)
            return None

        except requests.RequestException as e:
            wait = retry_delay * (2 ** (attempt - 1))
            log.warning(
                'Network error (attempt %d/%d): %s — retrying in %.0fs',
                attempt, max_retries, e, wait
            )
            time.sleep(wait)

    log.warning('All %d retries exhausted for %s — skipping study', max_retries, url)
    return None


def _get_studies(api_base: str, delay: float,
                 max_retries: int, retry_delay: float) -> list[dict]:
    """Retrieve list of public hg38 studies."""
    log.info('Fetching study list from cBioPortal...')
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(
                f'{api_base}/studies',
                params={'projection': 'SUMMARY'},
                timeout=SESSION_TIMEOUT,
            )
            time.sleep(delay)
            if resp.status_code == 200:
                studies = [s for s in resp.json() if s.get('publicStudy', False)]
                log.info('Retrieved %d public studies', len(studies))
                return studies
            if resp.status_code == 503:
                wait = retry_delay * (2 ** (attempt - 1))
                log.warning('HTTP 503 fetching study list — retrying in %.0fs', wait)
                time.sleep(wait)
                continue
            log.error('HTTP %d fetching study list', resp.status_code)
            return []
        except requests.RequestException as e:
            wait = retry_delay * (2 ** (attempt - 1))
            log.warning('Network error fetching study list: %s — retrying in %.0fs', e, wait)
            time.sleep(wait)
    log.error('Failed to retrieve study list after %d attempts', max_retries)
    return []


def _fetch_study_mutations(study_id: str,
                            api_base: str,
                            delay: float,
                            max_retries: int,
                            retry_delay: float) -> list[dict]:
    """Fetch all mutations for a study, paging through BATCH_SIZE at a time."""
    url  = f'{api_base}/molecular-profiles/{study_id}_mutations/mutations/fetch'
    body = {'sampleListId': f'{study_id}_all'}
    all_mutations = []
    page = 0

    while True:
        params = {
            'projection': 'DETAILED',
            'pageSize':   BATCH_SIZE,
            'pageNumber': page,
        }
        data = _post_with_retry(url, params, body, delay, max_retries, retry_delay)

        if data is None:
            # 404 = no profile (silent skip) or retries exhausted (already logged)
            break

        all_mutations.extend(data)

        if len(data) < BATCH_SIZE:
            break   # last page

        page += 1

    return all_mutations


# =============================================================================
# Row builder
# =============================================================================

def _mutations_to_rows(mutations: list[dict], cancer_type: str) -> list[dict]:
    """Convert raw cBioPortal mutation dicts to STANDARD_COLS rows."""
    rows = []
    for mut in mutations:
        try:
            chrom = normalise_chrom(str(mut.get('chr', '') or ''))
            pos   = mut.get('startPosition')
            ref   = clean_allele(str(mut.get('referenceAllele', '') or ''))
            alt   = clean_allele(str(mut.get('variantAllele', '') or ''))

            if not chrom or pos is None:
                continue
            if not is_valid_allele(ref) or not is_valid_allele(alt):
                continue

            gene      = str(mut.get('gene', {}).get('hugoGeneSymbol', '') or '')
            hgvsc     = str(mut.get('annotationSummary', {}).get(
                            'transcriptConsequenceSummary', {}).get(
                            'hgvscShort', '') or '')
            hgvsp     = str(mut.get('annotationSummary', {}).get(
                            'transcriptConsequenceSummary', {}).get(
                            'hgvspShort', '') or '')
            var_class = str(mut.get('mutationType', '') or '')

            rows.append({
                'chrom':       chrom,
                'pos':         int(pos),
                'ref':         ref,
                'alt':         alt,
                'gene':        gene,
                'hgvsc':       hgvsc,
                'hgvsp':       hgvsp,
                'consequence': map_consequence(var_class),
                'cancer_type': cancer_type,
                'n_samples':   1,
                'source':      'cBioPortal',
            })
        except (KeyError, TypeError, ValueError):
            continue
    return rows


# =============================================================================
# Main
# =============================================================================

def fetch(api_base: str,
          output_path: str,
          max_studies: int,
          delay: float,
          max_retries: int,
          retry_delay: float) -> None:

    studies = _get_studies(api_base, delay, max_retries, retry_delay)
    if not studies:
        log.error('No studies retrieved — aborting.')
        sys.exit(1)

    if max_studies > 0:
        log.info('Limiting to first %d studies (--max-studies)', max_studies)
        studies = studies[:max_studies]

    # Filter hg19 studies
    if SKIP_HG19:
        hg38_studies = [
            s for s in studies
            if s.get('referenceGenome', 'hg19').lower()
            not in ('hg19', 'grch37', 'grch37/hg19')
        ]
        log.info(
            'Skipping %d hg19 studies; processing %d hg38 studies',
            len(studies) - len(hg38_studies), len(hg38_studies)
        )
        studies = hg38_studies

    n_total   = len(studies)
    all_rows  = []
    n_skipped = 0
    n_failed  = 0

    for i, study in enumerate(studies):
        study_id    = study['studyId']
        cancer_type = study.get('cancerType', {}).get('name', study_id).lower()
        ref_genome  = study.get('referenceGenome', 'hg38').lower()

        log.info('[%d/%d] Fetching %s (%s)', i + 1, n_total, study_id, ref_genome)

        mutations = _fetch_study_mutations(
            study_id, api_base, delay, max_retries, retry_delay
        )

        if not mutations:
            n_skipped += 1
            continue

        rows = _mutations_to_rows(mutations, cancer_type)
        if not rows:
            n_skipped += 1
            continue

        all_rows.extend(rows)
        log.info(
            '  -> %d mutations added (total so far: %d)',
            len(rows), len(all_rows)
        )

    log.info('Fetch complete: %d rows from %d studies (%d skipped/empty, %d failed)',
             len(all_rows), n_total - n_skipped, n_skipped, n_failed)

    if not all_rows:
        log.error('No rows produced — output not written.')
        sys.exit(1)

    # Write output
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    df = pd.DataFrame(all_rows, columns=STANDARD_COLS)
    df.to_csv(output_path, sep='\t', index=False, compression='gzip')
    log.info('Snapshot written: %s  (%d rows)', output_path, len(df))
    log.info('')
    log.info('To use this snapshot in the pipeline:')
    log.info('  cp %s intermediate/cbioportal.tsv.gz', output_path)
    log.info('  bash run_pipeline.sh --from-intermediates')


def main():
    date_str = datetime.now().strftime('%Y%m%d')
    default_output = f'data/cbioportal/cbioportal_{date_str}.tsv.gz'

    ap = argparse.ArgumentParser(
        description='Fetch cBioPortal mutation data and save as a versioned snapshot.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument('--output', default=default_output,
                    help=f'Output path (default: {default_output})')
    ap.add_argument('--api-base', default=API_BASE,
                    help=f'cBioPortal API base URL (default: {API_BASE})')
    ap.add_argument('--max-studies', type=int, default=0,
                    help='Limit number of studies (0 = all). Useful for testing.')
    ap.add_argument('--delay', type=float, default=0.5,
                    help='Delay between requests in seconds (default: 0.5)')
    ap.add_argument('--max-retries', type=int, default=8,
                    help='Max retries per request on 503/network error (default: 8)')
    ap.add_argument('--retry-delay', type=float, default=10.0,
                    help='Initial retry delay in seconds, doubles each attempt (default: 10)')
    args = ap.parse_args()

    log.info('cBioPortal snapshot fetcher')
    log.info('  Output      : %s', args.output)
    log.info('  Max retries : %d', args.max_retries)
    log.info('  Retry delay : %.0fs (exponential backoff)', args.retry_delay)
    log.info('  Request delay: %.1fs', args.delay)
    if args.max_studies:
        log.info('  Max studies : %d (test mode)', args.max_studies)
    log.info('')

    fetch(
        api_base    = args.api_base,
        output_path = args.output,
        max_studies = args.max_studies,
        delay       = args.delay,
        max_retries = args.max_retries,
        retry_delay = args.retry_delay,
    )


if __name__ == '__main__':
    main()
