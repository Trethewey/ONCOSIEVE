#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Post-processor: MANE Select transcript remapping
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================
"""
mane_remap.py

Remaps HGVSc annotations in the ONCOSIEVE whitelist to MANE Select transcripts.

The whitelist aggregates HGVSc strings from multiple sources (COSMIC, GENIE,
TP53, cBioPortal etc.), each of which may have used a different transcript for
the same gene. This script:

  1. Downloads the NCBI MANE Select summary table (GRCh38).
  2. Builds a gene -> MANE Select RefSeq NM_ accession lookup.
  3. For each variant in the whitelist, checks whether the HGVSc transcript
     prefix matches the MANE Select transcript for that gene.
  4. For mismatches, queries the Ensembl VEP REST API using genomic coordinates
     (chrom, pos, ref, alt) to obtain the MANE Select HGVSc annotation.
  5. Writes an updated TSV with a new column `hgvsc_mane` and a flag column
     `mane_remapped` (True/False).

NOTE: This script is experimental. VEP API queries are rate-limited and
slow for large variant sets. Run with --max-variants to test on a subset.
Variants where remapping fails retain their original hgvsc value.

Requirements:
  pip install requests pandas tqdm --break-system-packages

Usage:
  python mane_remap.py \\
      --whitelist output/pan_cancer_whitelist_GRCh38.tsv.gz \\
      --output    output/pan_cancer_whitelist_GRCh38.mane.tsv.gz

  # Test on first 500 variants only
  python mane_remap.py \\
      --whitelist output/pan_cancer_whitelist_GRCh38.tsv.gz \\
      --output    output/pan_cancer_whitelist_GRCh38.mane.tsv.gz \\
      --max-variants 500
"""

import argparse
import gzip
import io
import logging
import os
import sys
import time
import urllib.request
from typing import Optional

import pandas as pd
import requests

# =============================================================================
# Configuration
# =============================================================================

# NCBI MANE Select summary — GRCh38, updated with each MANE release
_MANE_URL = (
    'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/'
    'MANE.GRCh38.v1.5.summary.txt.gz'
)

# Ensembl VEP REST API — GRCh38
_VEP_BASE = 'https://rest.ensembl.org/vep/human/region'
_VEP_HEADERS = {
    'Content-Type': 'application/json',
    'Accept':       'application/json',
}

# Batch size for VEP API (max 200 per request)
_VEP_BATCH = 100

# Delay between VEP API requests (seconds) — respect rate limit
_VEP_DELAY = 0.5

# Max retries per VEP batch
_VEP_RETRIES = 3

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [mane_remap] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
log = logging.getLogger('mane_remap')


# =============================================================================
# MANE Select lookup
# =============================================================================

def load_mane_table(cache_path: str = 'mane_select.tsv.gz') -> dict[str, str]:
    """
    Download and parse NCBI MANE Select summary.
    Returns dict: gene_symbol -> RefSeq NM_ accession (versioned, e.g. NM_000546.6)
    Caches locally to avoid repeated downloads.
    """
    if not os.path.exists(cache_path):
        log.info('Downloading MANE Select table from NCBI...')
        urllib.request.urlretrieve(_MANE_URL, cache_path)
        log.info('Saved to %s', cache_path)
    else:
        log.info('Using cached MANE table: %s', cache_path)

    with gzip.open(cache_path, 'rt') as fh:
        df = pd.read_csv(fh, sep='\t', dtype=str, low_memory=False)

    log.info('MANE table: %d rows, columns: %s', len(df), list(df.columns[:8]))

    # Column names vary slightly between MANE releases — find them flexibly
    gene_col  = _find_col(df, ['symbol', 'gene_symbol', 'HGNC_symbol', 'name'])
    nm_col    = _find_col(df, ['RefSeq_nuc', 'refseq_nuc', 'NM_accession'])
    type_col  = _find_col(df, ['MANE_status', 'mane_status', 'status'])

    if not gene_col or not nm_col:
        raise RuntimeError(
            f'Cannot find gene/NM columns in MANE table. '
            f'Available: {list(df.columns)}'
        )

    # Keep MANE Select only (not MANE Plus Clinical)
    if type_col:
        df = df[df[type_col].str.contains('Select', case=False, na=False)]

    lookup = {}
    for _, row in df.iterrows():
        gene = str(row[gene_col]).strip()
        nm   = str(row[nm_col]).strip()
        if gene and nm and nm.startswith('NM_'):
            lookup[gene] = nm

    log.info('MANE Select lookup: %d genes', len(lookup))
    return lookup


def _find_col(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    """Return the first candidate column name that exists in df."""
    for c in candidates:
        if c in df.columns:
            return c
    # Case-insensitive fallback
    lower = {col.lower(): col for col in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def hgvsc_transcript(hgvsc: str) -> str:
    """Extract the transcript accession from a HGVSc string e.g. NM_000546.6:c.817C>T -> NM_000546.6"""
    if ':' in hgvsc:
        return hgvsc.split(':')[0].strip()
    return ''


def transcript_base(accession: str) -> str:
    """Strip version from accession: NM_000546.6 -> NM_000546"""
    return accession.split('.')[0] if '.' in accession else accession


def needs_remap(hgvsc: str, mane_nm: str) -> bool:
    """Return True if hgvsc does not use the MANE Select transcript."""
    if not hgvsc or not mane_nm:
        return False
    tx = transcript_base(hgvsc_transcript(hgvsc))
    mane_base = transcript_base(mane_nm)
    return bool(tx) and tx != mane_base


# =============================================================================
# VEP API
# =============================================================================

def _vep_region_string(chrom: str, pos: int, ref: str, alt: str) -> str:
    """
    Format a variant as a VEP region string.
    VEP uses 1-based coordinates and strand (+).
    For SNVs: '1 12345 12345 A/T +'
    For insertions/deletions: coordinates adjusted per VEP convention.
    """
    chrom = str(chrom).lstrip('chr')
    pos   = int(pos)

    if len(ref) == len(alt):
        # SNV or MNV
        end = pos + len(ref) - 1
        return f'{chrom} {pos} {end} {ref}/{alt} +'
    elif len(ref) > len(alt):
        # Deletion
        start = pos + 1
        end   = pos + len(ref) - 1
        return f'{chrom} {start} {end} {ref[1:]}/- +'
    else:
        # Insertion
        return f'{chrom} {pos} {pos} -/{alt[1:]} +'


def _parse_vep_response(results: list, mane_lookup: dict[str, str]) -> dict[str, str]:
    """
    Parse VEP batch response, extracting MANE Select HGVSc for each variant.
    Returns dict: vep_region_string -> hgvsc_mane
    """
    out = {}
    for result in results:
        region = result.get('input', '')
        transcript_consequences = result.get('transcript_consequences', [])

        mane_hgvsc = ''
        for tc in transcript_consequences:
            # VEP flags MANE Select transcripts directly
            if tc.get('mane_select') or tc.get('mane_plus_clinical') == 'NM':
                mane_hgvsc = tc.get('hgvsc', '')
                if mane_hgvsc:
                    break

            # Fallback: match against our MANE lookup by gene
            gene = tc.get('gene_symbol', '')
            if gene in mane_lookup:
                nm = mane_lookup[gene]
                tx = tc.get('transcript_id', '')
                if transcript_base(tx) == transcript_base(nm):
                    mane_hgvsc = tc.get('hgvsc', '')
                    if mane_hgvsc:
                        break

        if mane_hgvsc:
            out[region] = mane_hgvsc

    return out


def query_vep_batch(regions: list[str], mane_lookup: dict[str, str]) -> dict[str, str]:
    """
    Query Ensembl VEP REST API for a batch of region strings.
    Returns dict: region_string -> mane_hgvsc
    """
    payload = {'variants': regions}

    for attempt in range(_VEP_RETRIES):
        try:
            resp = requests.post(
                _VEP_BASE,
                headers=_VEP_HEADERS,
                json=payload,
                timeout=60,
            )
            if resp.status_code == 200:
                return _parse_vep_response(resp.json(), mane_lookup)
            elif resp.status_code == 429:
                wait = int(resp.headers.get('Retry-After', 5))
                log.warning('VEP rate limit hit — waiting %ds', wait)
                time.sleep(wait)
            else:
                log.warning('VEP returned %d: %s', resp.status_code, resp.text[:200])
                break
        except requests.RequestException as e:
            log.warning('VEP request error (attempt %d): %s', attempt + 1, e)
            time.sleep(2 ** attempt)

    return {}


# =============================================================================
# Main remapping logic
# =============================================================================

def remap_whitelist(
    whitelist_path: str,
    output_path: str,
    mane_cache: str = 'mane_select.tsv.gz',
    max_variants: int = 0,
) -> None:

    # Load whitelist
    log.info('Loading whitelist: %s', whitelist_path)
    df = pd.read_csv(whitelist_path, sep='\t', dtype=str, low_memory=False)
    log.info('Loaded %d variants', len(df))

    if max_variants > 0:
        df = df.head(max_variants)
        log.info('Testing on first %d variants', len(df))

    # Load MANE lookup
    mane_lookup = load_mane_table(mane_cache)

    # Initialise output columns
    df['hgvsc_mane']   = df['hgvsc']
    df['mane_remapped'] = False

    # Identify variants that need remapping
    remap_mask = df.apply(
        lambda row: needs_remap(
            str(row.get('hgvsc', '')),
            mane_lookup.get(str(row.get('gene', '')), '')
        ),
        axis=1,
    )

    n_needs_remap = remap_mask.sum()
    n_already_mane = (~remap_mask & df['hgvsc'].str.startswith('NM_', na=False)).sum()
    log.info('Variants already on MANE Select transcript: %d', n_already_mane)
    log.info('Variants requiring remapping: %d', n_needs_remap)

    if n_needs_remap == 0:
        log.info('No remapping required. Writing output.')
    else:
        # Build VEP region strings for variants needing remap
        remap_df = df[remap_mask].copy()
        remap_df['_region'] = remap_df.apply(
            lambda row: _vep_region_string(
                row['chrom'], row['pos'], row['ref'], row['alt']
            ),
            axis=1,
        )

        # Batch query VEP
        regions = remap_df['_region'].tolist()
        vep_results: dict[str, str] = {}
        n_batches = (len(regions) + _VEP_BATCH - 1) // _VEP_BATCH

        log.info('Querying Ensembl VEP in %d batches of %d...', n_batches, _VEP_BATCH)

        for i in range(0, len(regions), _VEP_BATCH):
            batch = regions[i:i + _VEP_BATCH]
            batch_num = i // _VEP_BATCH + 1
            log.info('  Batch %d / %d', batch_num, n_batches)
            results = query_vep_batch(batch, mane_lookup)
            vep_results.update(results)
            time.sleep(_VEP_DELAY)

        # Apply remapped HGVSc back to dataframe
        n_remapped = 0
        for idx, row in remap_df.iterrows():
            region = row['_region']
            if region in vep_results:
                df.at[idx, 'hgvsc_mane']    = vep_results[region]
                df.at[idx, 'mane_remapped']  = True
                n_remapped += 1

        n_failed = n_needs_remap - n_remapped
        log.info('Successfully remapped: %d', n_remapped)
        if n_failed > 0:
            log.warning(
                '%d variants could not be remapped — original hgvsc retained', n_failed
            )

    # Write output
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    df.drop(columns=['_region'], errors='ignore').to_csv(
        output_path, sep='\t', index=False,
        compression='gzip' if output_path.endswith('.gz') else None,
    )
    log.info('Written: %s  (%d rows)', output_path, len(df))

    # Summary
    n_success = df['mane_remapped'].astype(str).str.lower().eq('true').sum()
    log.info(
        'Summary: %d total | %d already MANE | %d remapped | %d unchanged',
        len(df), n_already_mane, n_success, len(df) - n_already_mane - n_success,
    )


# =============================================================================
# Entry point
# =============================================================================

def main() -> None:
    ap = argparse.ArgumentParser(
        description='Remap ONCOSIEVE whitelist HGVSc annotations to MANE Select transcripts.'
    )
    ap.add_argument(
        '--whitelist', required=True,
        help='Input whitelist TSV (output/pan_cancer_whitelist_GRCh38.tsv.gz)',
    )
    ap.add_argument(
        '--output', required=True,
        help='Output TSV path (.tsv.gz recommended)',
    )
    ap.add_argument(
        '--mane-cache', default='mane_select.tsv.gz',
        help='Local cache path for MANE Select table (default: mane_select.tsv.gz)',
    )
    ap.add_argument(
        '--max-variants', type=int, default=0,
        help='Limit to first N variants for testing (0 = no limit)',
    )
    args = ap.parse_args()

    remap_whitelist(
        whitelist_path = args.whitelist,
        output_path    = args.output,
        mane_cache     = args.mane_cache,
        max_variants   = args.max_variants,
    )


if __name__ == '__main__':
    main()
