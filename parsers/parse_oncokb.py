#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: OncoKB annotation API
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_oncokb.py
Query the OncoKB public API to annotate variants with oncogenicity.

No API token is required for the public endpoint. The endpoint returns
the 'oncogenic' field (Oncogenic, Likely Oncogenic, Predicted Oncogenic,
Likely Neutral, Unknown) for all submitted variants.

If a local allAnnotatedVariants.txt file exists (from a token-based download),
it is used instead to avoid API calls.
"""

import os
import re
import time

import pandas as pd
import requests

from parsers.common import (
    STANDARD_COLS,
    empty_standard_df,
    setup_logger,
)

log = setup_logger('OncoKB')

_API_BASE       = 'https://www.oncokb.org/api/v1'
_API_TOKEN      = None  # set via api_token argument — never hardcode here
_BATCH_ENDPOINT = f'{_API_BASE}/annotate/mutations/byProteinChange'
_BATCH_SIZE     = 50
_RETRY_DELAY    = 2.0
_REQUEST_DELAY  = 0.5
_MAX_RETRIES    = 3

_ONCOGENICITY_INCLUDE = {
    'Oncogenic',
    'Likely Oncogenic',
    'Predicted Oncogenic',
}

_COL_GENE   = 'Hugo Symbol'
_COL_ALT    = 'Alteration'
_COL_EFFECT = 'Mutation Effect'
_COL_ONCO   = 'Oncogenicity'
_COL_CTYPE  = 'Cancer Type'


def parse_oncokb(variants_file: str,
                 merged_df=None,
                 include_oncogenicity: list | None = None,
                 api_token: str | None = None) -> pd.DataFrame:
    """
    Annotate variants with OncoKB oncogenicity.

    If variants_file exists on disk, it is parsed directly.
    Otherwise, the public API is queried using (gene, hgvsp) pairs from merged_df.
    """
    if include_oncogenicity is None:
        include_oncogenicity = list(_ONCOGENICITY_INCLUDE)

    if os.path.exists(variants_file):
        log.info('OncoKB: using local file %s', variants_file)
        return _parse_file(variants_file, include_oncogenicity)

    if merged_df is None or merged_df.empty:
        log.warning('OncoKB: no local file and no merged_df provided — skipping')
        return empty_standard_df()

    if not api_token:
        log.warning('OncoKB: no api_token provided — skipping API query')
        return empty_standard_df()
    log.info('OncoKB: querying public API for %d unique variants', len(merged_df))
    return _query_api(merged_df, include_oncogenicity, api_token)


_COL_ALIASES = {
    _COL_GENE:   ['Hugo Symbol', 'hugoSymbol', 'Gene', 'gene'],
    _COL_ALT:    ['Alteration', 'alteration', 'Variant', 'variant', 'Protein Change'],
    _COL_EFFECT: ['Mutation Effect', 'mutationEffect', 'Mutation_Effect', 'mutation_effect'],
    _COL_ONCO:   ['Oncogenicity', 'oncogenicity', 'oncogenic'],
    _COL_CTYPE:  ['Cancer Type', 'cancerType', 'Cancer_Type', 'cancer_type'],
}


def _normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to canonical names using _COL_ALIASES."""
    rename = {}
    for canonical, aliases in _COL_ALIASES.items():
        if canonical in df.columns:
            continue
        for alias in aliases:
            if alias in df.columns:
                rename[alias] = canonical
                break
    if rename:
        log.info('OncoKB: renaming columns: %s', rename)
        df = df.rename(columns=rename)
    return df


def _parse_file(variants_file: str, include_oncogenicity: list) -> pd.DataFrame:
    try:
        df_raw = pd.read_csv(variants_file, sep='\t', dtype=str, low_memory=False)
    except Exception as e:
        log.error('Failed to read OncoKB file: %s', e)
        return empty_standard_df()

    df_raw = _normalise_columns(df_raw)

    required = [_COL_GENE, _COL_ALT, _COL_EFFECT, _COL_ONCO]
    missing = [c for c in required if c not in df_raw.columns]
    if missing:
        log.error('OncoKB file missing columns: %s  (available: %s)', missing, list(df_raw.columns))
        return empty_standard_df()

    df_raw = df_raw.fillna('')
    df_f = df_raw[df_raw[_COL_ONCO].isin(include_oncogenicity)].copy()
    log.info('OncoKB file: %d / %d entries pass filter', len(df_f), len(df_raw))

    return _build_output(
        genes=[str(r[_COL_GENE]).strip() for _, r in df_f.iterrows()],
        alts=[str(r[_COL_ALT]).strip() for _, r in df_f.iterrows()],
        effects=[str(r[_COL_EFFECT]).strip() for _, r in df_f.iterrows()],
        oncogenicities=df_f[_COL_ONCO].tolist(),
        cancer_types=[str(r.get(_COL_CTYPE, '')).strip() for _, r in df_f.iterrows()],
    )


def _query_api(merged_df: pd.DataFrame, include_oncogenicity: list, api_token: str) -> pd.DataFrame:
    pairs = (
        merged_df[['gene', 'hgvsp']]
        .dropna()
        .drop_duplicates()
        .query("gene != '' and hgvsp != ''")
    )

    if pairs.empty:
        log.warning('OncoKB: no valid (gene, hgvsp) pairs to query')
        return empty_standard_df()

    all_results = []
    total = len(pairs)
    n_batches = (total + _BATCH_SIZE - 1) // _BATCH_SIZE
    log.info('OncoKB API: %d variants in %d batches', total, n_batches)

    for batch_idx in range(n_batches):
        batch = pairs.iloc[batch_idx * _BATCH_SIZE:(batch_idx + 1) * _BATCH_SIZE]
        payload = [
            {
                'id': str(i),
                'hugoSymbol': str(row['gene']).strip(),
                'alteration': (str(row['hgvsp'])[2:]
                               if str(row['hgvsp']).startswith('p.')
                               else str(row['hgvsp'])),
            }
            for i, (_, row) in enumerate(batch.iterrows())
        ]

        for attempt in range(_MAX_RETRIES):
            try:
                resp = requests.post(
                    _BATCH_ENDPOINT,
                    json=payload,
                    headers={'Content-Type': 'application/json',
                             'Authorization': f'Bearer {api_token}',
                             'Accept': 'application/json'},
                    timeout=30,
                )
                if resp.status_code == 200:
                    all_results.extend(resp.json())
                    break
                elif resp.status_code == 429:
                    wait = _RETRY_DELAY * (attempt + 2)
                    log.warning('OncoKB rate limited — waiting %.1fs', wait)
                    time.sleep(wait)
                else:
                    log.warning('OncoKB batch %d: HTTP %d — %s',
                                batch_idx, resp.status_code, resp.text[:200])
                    break
            except requests.RequestException as e:
                log.warning('OncoKB batch %d attempt %d: %s', batch_idx, attempt + 1, e)
                time.sleep(_RETRY_DELAY)

        if batch_idx % 10 == 0:
            log.info('OncoKB API: processed %d / %d',
                     min((batch_idx + 1) * _BATCH_SIZE, total), total)
        time.sleep(_REQUEST_DELAY)

    if not all_results:
        log.warning('OncoKB API: no results returned')
        return empty_standard_df()

    return _parse_api_results(all_results, include_oncogenicity)


def _parse_api_results(results: list, include_oncogenicity: list) -> pd.DataFrame:
    genes, alts, effects, oncogenicities, cancer_types = [], [], [], [], []

    for item in results:
        onco = item.get('oncogenic', '') or ''
        if onco not in include_oncogenicity:
            continue
        query  = item.get('query', {})
        gene   = query.get('hugoSymbol', '') or ''
        alt    = query.get('alteration', '') or ''
        effect = (item.get('mutationEffect') or {}).get('knownEffect', '') or ''
        genes.append(gene)
        alts.append(alt)
        effects.append(effect)
        oncogenicities.append(onco)
        cancer_types.append('pan_cancer')

    log.info('OncoKB API: %d / %d pass oncogenicity filter', len(genes), len(results))
    if not genes:
        return empty_standard_df()

    return _build_output(genes, alts, effects, oncogenicities, cancer_types)


def _build_output(genes, alts, effects, oncogenicities, cancer_types) -> pd.DataFrame:
    rows = []
    for gene, alt, effect, onco, ctype in zip(
            genes, alts, effects, oncogenicities, cancer_types):
        rows.append({
            'chrom':       '',
            'pos':         pd.NA,
            'ref':         '',
            'alt':         '',
            'gene':        gene,
            'hgvsc':       '',
            'hgvsp':       _normalise_hgvsp(alt),
            'consequence': _map_effect(effect),
            'cancer_type': (ctype or 'pan_cancer').lower(),
            'n_samples':   0,
            'source':      f'OncoKB:{onco}',
        })
    df_out = pd.DataFrame(rows, columns=STANDARD_COLS)
    df_out['oncokb_oncogenicity'] = oncogenicities
    log.info('OncoKB: %d entries retained', len(df_out))
    return df_out


def _map_effect(effect: str) -> str:
    e = effect.lower()
    if 'gain' in e:
        return 'missense'
    if 'loss' in e:
        return 'missense'
    if 'neutral' in e:
        return 'synonymous'
    return 'unknown'


def _normalise_hgvsp(alteration: str) -> str:
    alt = alteration.strip()
    if alt.startswith('p.'):
        return alt
    if re.match(r'^[A-Z]\d+[A-Z*]$', alt):
        return f'p.{alt}'
    return alt
