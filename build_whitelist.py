#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Main pipeline: aggregate multi-source somatic variant data and build the whitelist
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
build_whitelist.py
Main pipeline script. Runs all parsers, merges outputs, applies filters,
assigns tiers, and writes the final whitelist in TSV and VCF formats.

Usage:
    python build_whitelist.py --config config.yaml [--skip-sources cosmic,genie]
"""

import argparse
import gzip
import logging
import os
import sys
from datetime import datetime

import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import yaml

from parsers.common import INCLUDED_CONSEQUENCES, STANDARD_COLS, setup_logger

log = setup_logger('build_whitelist')

# Sources that contribute to the sample count threshold
COUNT_SOURCES = {'COSMIC', 'cBioPortal', 'GENIE', 'TP53_somatic', 'TP53_germline'}

# Sources that are annotation/curated (not count-based)
CURATED_SOURCES = {'OncoKB', 'ClinVar', 'CancerHotspots'}

# VCF header template
_VCF_HEADER = """\
##fileformat=VCFv4.2
##fileDate={date}
##source=pan_cancer_whitelist_pipeline
##reference=GRCh38/hg38
##contig=<ID=1,length=248956422,assembly=GRCh38>
##contig=<ID=2,length=242193529,assembly=GRCh38>
##contig=<ID=3,length=198295559,assembly=GRCh38>
##contig=<ID=4,length=190214555,assembly=GRCh38>
##contig=<ID=5,length=181538259,assembly=GRCh38>
##contig=<ID=6,length=170805979,assembly=GRCh38>
##contig=<ID=7,length=159345973,assembly=GRCh38>
##contig=<ID=8,length=145138636,assembly=GRCh38>
##contig=<ID=9,length=138394717,assembly=GRCh38>
##contig=<ID=10,length=133797422,assembly=GRCh38>
##contig=<ID=11,length=135086622,assembly=GRCh38>
##contig=<ID=12,length=133275309,assembly=GRCh38>
##contig=<ID=13,length=114364328,assembly=GRCh38>
##contig=<ID=14,length=107043718,assembly=GRCh38>
##contig=<ID=15,length=101991189,assembly=GRCh38>
##contig=<ID=16,length=90338345,assembly=GRCh38>
##contig=<ID=17,length=83257441,assembly=GRCh38>
##contig=<ID=18,length=80373285,assembly=GRCh38>
##contig=<ID=19,length=58617616,assembly=GRCh38>
##contig=<ID=20,length=64444167,assembly=GRCh38>
##contig=<ID=21,length=46709983,assembly=GRCh38>
##contig=<ID=22,length=50818468,assembly=GRCh38>
##contig=<ID=X,length=156040895,assembly=GRCh38>
##contig=<ID=Y,length=57227415,assembly=GRCh38>
##contig=<ID=MT,length=16569,assembly=GRCh38>
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=CSQ,Number=1,Type=String,Description="Standardised consequence term">
##INFO=<ID=HGVSc,Number=1,Type=String,Description="HGVSc notation">
##INFO=<ID=HGVSp,Number=1,Type=String,Description="HGVSp notation">
##INFO=<ID=N_SAMPLES,Number=1,Type=Integer,Description="Total samples across count-based sources">
##INFO=<ID=N_CANCER_TYPES,Number=1,Type=Integer,Description="Number of distinct cancer types">
##INFO=<ID=CANCER_TYPES,Number=.,Type=String,Description="Pipe-delimited list of cancer types">
##INFO=<ID=SOURCES,Number=.,Type=String,Description="Pipe-delimited list of sources">
##INFO=<ID=WL_TIER,Number=1,Type=Integer,Description="Whitelist tier: 1=highest, 3=minimum threshold">
##INFO=<ID=ONCOKB,Number=1,Type=String,Description="OncoKB oncogenicity classification (if available)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def _apply_data_dir(cfg: dict, data_dir: str) -> dict:
    """
    Replace the leading path component of all relative file paths in
    data_sources with data_dir. The config.yaml uses paths like
    'data/cosmic/file.tsv' — this replaces 'data/' with the supplied
    data_dir so paths are never doubled.
    """
    if not data_dir:
        return cfg
    data_dir = data_dir.rstrip('/')
    path_keys = {'tsv', 'vcf', 'maf', 'clinical_sample', 'somatic_tsv',
                 'germline_tsv', 'variants_file', 'dir', 'chain'}

    def _repath(v: str) -> str:
        if not v or os.path.isabs(v):
            return v
        # Strip any leading path component (e.g. 'data/') then join with data_dir
        parts = v.replace('\\', '/').split('/', 1)
        remainder = parts[1] if len(parts) > 1 else parts[0]
        return os.path.join(data_dir, remainder)

    for source in cfg.get('data_sources', {}).values():
        if not isinstance(source, dict):
            continue
        for k, v in source.items():
            if k in path_keys and isinstance(v, str):
                source[k] = _repath(v)
    for k, v in cfg.get('reference', {}).items():
        if isinstance(v, str):
            cfg['reference'][k] = _repath(v)
    return cfg


def load_config(path: str, data_dir: str = '') -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)

    # Load and merge settings.yaml if specified
    settings_file = cfg.pop('settings_file', 'settings.yaml')
    if os.path.exists(settings_file):
        with open(settings_file) as fh:
            settings = yaml.safe_load(fh) or {}
        cfg['thresholds']            = settings.get('thresholds', {})
        cfg['tiering']               = settings.get('tiering', {})
        cfg['vaf_rescue']            = settings.get('vaf_rescue', {})
        cfg['included_consequences'] = settings.get('included_consequences', [])
        cfg['log_level']             = settings.get('log_level', 'INFO')
        for src in ('oncokb', 'clinvar', 'cancer_hotspots', 'cbioportal'):
            if src in settings and src in cfg.get('data_sources', {}):
                cfg['data_sources'][src].update(settings[src])
        log.info('Settings loaded from: %s', settings_file)
    else:
        log.warning('settings.yaml not found — using defaults')

    return cfg


def _finalise_config(cfg: dict, data_dir: str = '') -> dict:
    """Apply data_dir prefix after config is fully loaded."""
    if data_dir:
        cfg = _apply_data_dir(cfg, data_dir)
        log.info('Data directory: %s', data_dir)
    return cfg


def run_parsers(cfg: dict, skip: set, inter_dir: str = 'intermediate') -> dict[str, pd.DataFrame]:
    """Run all enabled parsers and save each intermediate immediately on completion."""
    from parsers.parse_cosmic      import parse_cosmic
    from parsers.parse_oncokb      import parse_oncokb
    from parsers.parse_cbioportal  import parse_cbioportal
    from parsers.parse_clinvar     import parse_clinvar
    from parsers.parse_genie       import parse_genie
    from parsers.parse_tp53        import parse_tp53
    from parsers.parse_hotspots    import parse_hotspots

    frames: dict[str, pd.DataFrame] = {}
    ds = cfg.get('data_sources', {})
    n_threads = cfg.get('performance', {}).get('threads', 4)
    os.makedirs(inter_dir, exist_ok=True)

    def _save(name: str, df: pd.DataFrame) -> None:
        if not df.empty:
            ipath = os.path.join(inter_dir, f'{name.lower()}.tsv.gz')
            df.to_csv(ipath, sep='\t', index=False, compression='gzip')
            log.info('Saved intermediate: %s  (%d rows)', ipath, len(df))

    def _should_run(name: str) -> bool:
        return name not in skip and ds.get(name, {}).get('enabled', True)

    if _should_run('cosmic'):
        log.info('=== Running COSMIC parser ===')
        frames['COSMIC'] = parse_cosmic(
            tsv_path   = ds['cosmic']['tsv'],
            vcf_path   = ds['cosmic'].get('vcf'),
            chunk_size = ds['cosmic'].get('chunk_size', 200_000),
            n_threads  = n_threads,
        )
        _save('COSMIC', frames['COSMIC'])

    if _should_run('cbioportal'):
        log.info('=== Running cBioPortal parser ===')
        frames['cBioPortal'] = parse_cbioportal(
            api_base                = ds['cbioportal'].get('api_base', 'https://www.cbioportal.org/api'),
            use_tcga_pancancer_only = ds['cbioportal'].get('use_tcga_pancancer_only', False),
            max_studies             = ds['cbioportal'].get('max_studies', 0),
            request_delay_s         = ds['cbioportal'].get('request_delay_s', 0.5),
        )

    if _should_run('clinvar'):
        log.info('=== Running ClinVar parser ===')
        frames['ClinVar'] = parse_clinvar(
            vcf_path       = ds['clinvar']['vcf'],
            include_clinsig = ds['clinvar'].get('include_clinsig'),
            somatic_only    = ds['clinvar'].get('somatic_only', True),
        )
        _save('ClinVar', frames['ClinVar'])

    if _should_run('genie'):
        log.info('=== Running GENIE parser ===')
        genie_maf = ds['genie']['maf']
        genie_dir = os.path.dirname(genie_maf)
        clinical  = os.path.join(genie_dir, 'data_clinical_sample.txt')
        frames['GENIE'] = parse_genie(
            maf_path              = genie_maf,
            clinical_sample_path  = clinical if os.path.exists(clinical) else None,
        )
        _save('GENIE', frames['GENIE'])

    # Build task list for parallel execution
    tasks = {}

    if _should_run('tp53'):
        tasks['TP53'] = lambda: parse_tp53(
            somatic_tsv      = ds['tp53']['somatic_tsv'],
            germline_tsv     = ds['tp53'].get('germline_tsv'),
            include_germline = ds['tp53'].get('include_germline', False),
        )

    if _should_run('cancer_hotspots'):
        tasks['CancerHotspots'] = lambda: parse_hotspots(
            tsv_path  = ds['cancer_hotspots'].get('tsv'),
            max_qvalue= ds['cancer_hotspots'].get('max_qvalue', 0.05),
        )


    if tasks:
        log.info('Running %d parsers in parallel (threads=%d)...', len(tasks), n_threads)
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            future_to_name = {executor.submit(fn): name for name, fn in tasks.items()}
            for future in as_completed(future_to_name):
                name = future_to_name[future]
                try:
                    frames[name] = future.result()
                    log.info('=== %s parser complete ===', name)
                    _save(name, frames[name])
                except Exception as exc:
                    log.error('%s parser failed: %s', name, exc)

    if _should_run('oncokb'):
        log.info('=== Running OncoKB parser ===')
        all_so_far = [df for df in frames.values() if not df.empty]
        merged_so_far = pd.concat(all_so_far, ignore_index=True) if all_so_far else pd.DataFrame()
        frames['OncoKB'] = parse_oncokb(
            variants_file        = ds['oncokb']['variants_file'],
            merged_df            = merged_so_far,
            include_oncogenicity = ds['oncokb'].get('include_oncogenicity'),
            api_token            = ds['oncokb'].get('api_token'),
        )
        _save('OncoKB', frames['OncoKB'])
    return frames



def _first_nonempty(series) -> str:
    """Return the first non-empty, non-null value from a pandas Series, or ''."""
    for val in series:
        if val and str(val).strip() and str(val).strip().lower() not in ('nan', 'none', 'na'):
            return str(val).strip()
    return ''


def merge_and_aggregate(frames: dict[str, pd.DataFrame],
                         included_consequences: set[str]) -> pd.DataFrame:
    """
    Merge all parser outputs and aggregate by variant key.

    Uses Polars for multi-threaded groupby, falling back to pandas if Polars
    is not installed. On a 40-core server, Polars completes in 2-3 minutes
    vs 6+ hours for the pandas single-threaded version.
    """
    all_dfs = [df for df in frames.values() if not df.empty]
    if not all_dfs:
        log.error('No parser produced any output')
        sys.exit(1)

    raw = pd.concat(all_dfs, ignore_index=True)
    log.info('Total raw rows across all sources: %d', len(raw))

    # Filter consequence types
    raw = raw[raw['consequence'].isin(included_consequences)].copy()
    log.info('After consequence filter: %d rows', len(raw))

    # Separate coordinate-resolved and coordinate-missing rows
    has_coords = raw['chrom'].notna() & (raw['chrom'] != '') & \
                 raw['pos'].notna() & (raw['ref'] != '') & (raw['alt'] != '')
    df_coord   = raw[has_coords].copy()
    df_nocoord = raw[~has_coords].copy()

    log.info('Rows with coordinates: %d', len(df_coord))
    log.info('Rows without coordinates (annotation-only): %d', len(df_nocoord))

    # Normalise chrom — strip chr prefix for consistency
    df_coord['chrom'] = df_coord['chrom'].astype(str).str.replace(r'^chr', '', regex=True)
    df_coord['pos']   = pd.to_numeric(df_coord['pos'], errors='coerce')
    df_coord['n_samples'] = pd.to_numeric(df_coord['n_samples'], errors='coerce').fillna(0)

    # Tag count-source rows
    count_pattern = '|'.join(COUNT_SOURCES)
    df_coord['_is_count_source'] = df_coord['source'].str.contains(count_pattern, na=False)

    log.info('Aggregating by variant key...')

    try:
        import polars as pl
        log.info('Using Polars for multi-threaded aggregation')
        df_agg = _aggregate_polars(df_coord)
    except ImportError:
        log.warning('Polars not installed — falling back to pandas (slow). Run: pip install polars')
        df_agg = _aggregate_pandas(df_coord)

    log.info('Aggregated: %d unique variants', len(df_agg))
    return df_agg


def _clean_cancer_type(val: str) -> str:
    """Sanitise a single cancer type string."""
    return str(val).replace('\t', ' ').replace('\n', ' ').replace('|', '-').strip()


def _aggregate_polars(df_coord: pd.DataFrame) -> pd.DataFrame:
    """Multi-threaded aggregation using Polars."""
    import polars as pl

    # Ensure oncokb_oncogenicity column exists
    if 'oncokb_oncogenicity' not in df_coord.columns:
        df_coord['oncokb_oncogenicity'] = ''

    # Fill nulls for string columns before converting
    str_cols = ['chrom', 'ref', 'alt', 'source', 'cancer_type', 'gene',
                'consequence', 'hgvsc', 'hgvsp', 'oncokb_oncogenicity']
    for c in str_cols:
        if c in df_coord.columns:
            df_coord[c] = df_coord[c].fillna('').astype(str)

    lf = pl.from_pandas(df_coord).lazy()

    key_cols = ['chrom', 'pos', 'ref', 'alt']

    def first_nonempty(col):
        """Return first non-empty, non-null, non-'nan' value."""
        return (
            pl.col(col)
            .filter(
                pl.col(col).is_not_null() &
                (pl.col(col).str.strip_chars() != '') &
                (pl.col(col).str.to_lowercase() != 'nan') &
                (pl.col(col).str.to_lowercase() != 'none') &
                (pl.col(col).str.to_lowercase() != 'na')
            )
            .first()
            .fill_null('')
            .alias(col)
        )

    # Sources: sorted unique pipe-joined
    sources_expr = (
        pl.col('source')
        .filter(pl.col('source').is_not_null() & (pl.col('source') != ''))
        .unique()
        .sort()
        .str.join('|')
        .alias('sources')
    )

    # Cancer types: sanitise, sorted unique pipe-joined
    cancer_expr = (
        pl.col('cancer_type')
        .filter(
            pl.col('cancer_type').is_not_null() &
            (pl.col('cancer_type') != '') &
            (pl.col('cancer_type').str.to_lowercase() != 'unspecified')
        )
        .map_elements(lambda v: _clean_cancer_type(v), return_dtype=pl.String)
        .unique()
        .sort()
        .str.join('|')
        .fill_null('unspecified')
        .alias('cancer_types')
    )

    n_cancer_expr = (
        pl.col('cancer_type')
        .filter(
            pl.col('cancer_type').is_not_null() &
            (pl.col('cancer_type') != '') &
            (pl.col('cancer_type').str.to_lowercase() != 'unspecified')
        )
        .n_unique()
        .alias('n_cancer_types')
    )

    # n_samples from count sources only
    n_samples_expr = (
        pl.when(pl.col('_is_count_source'))
        .then(pl.col('n_samples'))
        .otherwise(0.0)
        .sum()
        .cast(pl.Int64)
        .alias('n_samples')
    )

    df_agg = (
        lf.group_by(key_cols)
        .agg([
            sources_expr,
            cancer_expr,
            n_cancer_expr,
            n_samples_expr,
            first_nonempty('gene'),
            first_nonempty('consequence'),
            first_nonempty('hgvsc'),
            first_nonempty('hgvsp'),
            first_nonempty('oncokb_oncogenicity'),
        ])
        .collect(streaming=True)
    ).to_pandas()

    df_agg['n_cancer_types'] = df_agg['n_cancer_types'].fillna(1).astype(int)
    df_agg['n_samples']      = df_agg['n_samples'].fillna(0).astype(int)
    df_agg['pos']            = df_agg['pos'].astype(int)
    return df_agg


def _aggregate_pandas(df_coord: pd.DataFrame) -> pd.DataFrame:
    """Single-threaded pandas fallback aggregation."""
    key_cols = ['chrom', 'pos', 'ref', 'alt']

    if 'oncokb_oncogenicity' not in df_coord.columns:
        df_coord['oncokb_oncogenicity'] = ''

    def first_nonempty_agg(s):
        vals = s.dropna()
        vals = vals[~vals.str.strip().str.lower().isin(['', 'nan', 'none', 'na'])]
        return vals.iloc[0] if len(vals) else ''

    def clean_ct_agg(x):
        cleaned = []
        for t in x.dropna().unique():
            t = _clean_cancer_type(str(t))
            if t and t.lower() != 'unspecified':
                cleaned.append(t)
        return '|'.join(sorted(cleaned)) or 'unspecified'

    sources_agg        = df_coord.groupby(key_cols)['source'].agg(lambda x: '|'.join(sorted(x.dropna().unique())))
    cancer_type_agg    = df_coord.groupby(key_cols)['cancer_type'].agg(clean_ct_agg)
    n_cancer_types_agg = df_coord.groupby(key_cols)['cancer_type'].agg(
        lambda x: len(set(t for t in x.dropna().unique() if t and t != 'unspecified')) or 1)
    n_samples_agg      = df_coord[df_coord['_is_count_source']].groupby(key_cols)['n_samples'].sum()
    gene_agg           = df_coord.groupby(key_cols)['gene'].agg(first_nonempty_agg)
    consequence_agg    = df_coord.groupby(key_cols)['consequence'].agg(first_nonempty_agg)
    hgvsc_agg          = df_coord.groupby(key_cols)['hgvsc'].agg(first_nonempty_agg)
    hgvsp_agg          = df_coord.groupby(key_cols)['hgvsp'].agg(first_nonempty_agg)
    oncokb_agg         = df_coord.groupby(key_cols)['oncokb_oncogenicity'].agg(first_nonempty_agg)

    df_agg = pd.DataFrame({
        'sources':             sources_agg,
        'cancer_types':        cancer_type_agg,
        'n_cancer_types':      n_cancer_types_agg,
        'n_samples':           n_samples_agg,
        'gene':                gene_agg,
        'consequence':         consequence_agg,
        'hgvsc':               hgvsc_agg,
        'hgvsp':               hgvsp_agg,
        'oncokb_oncogenicity': oncokb_agg,
    }).reset_index()

    df_agg['n_samples']   = df_agg['n_samples'].fillna(0).astype(int)
    df_agg['pos']         = df_agg['pos'].astype(int)
    return df_agg


def apply_filters(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    Apply inclusion filters.

    A variant passes if ANY of the following conditions is true:
      1. n_samples >= min_samples AND n_cancer_types >= min_cancer_types
         (count-based threshold)
      2. OncoKB oncogenicity = Oncogenic | Likely Oncogenic | Predicted Oncogenic
         (expert-curated override)
      3. ClinVar is a source (pathogenic somatic evidence)
      4. CancerHotspots is a source (statistically significant hotspot)
    """
    thr = cfg.get('thresholds', {})
    min_samples      = thr.get('min_samples_total', 10)
    min_cancer_types = thr.get('min_cancer_types', 1)

    passes_count   = (df['n_samples'] >= min_samples) & \
                     (df['n_cancer_types'] >= min_cancer_types)
    passes_oncokb  = df['oncokb_oncogenicity'].isin(
                         ['Oncogenic', 'Likely Oncogenic', 'Predicted Oncogenic'])
    passes_clinvar = df['sources'].str.contains('ClinVar', na=False)
    passes_hotspot = df['sources'].str.contains('CancerHotspots', na=False)

    mask = passes_count | passes_oncokb | passes_clinvar | passes_hotspot

    df_pass = df[mask].copy()
    log.info('After filters: %d / %d variants retained', len(df_pass), len(df))
    log.info('  Count-based:   %d', passes_count.sum())
    log.info('  OncoKB:        %d', passes_oncokb.sum())
    log.info('  ClinVar:       %d', passes_clinvar.sum())
    log.info('  CancerHotspot: %d', passes_hotspot.sum())
    return df_pass


def assign_tiers(df: pd.DataFrame, cfg: dict | None = None) -> pd.DataFrame:
    """
    Assign whitelist tier to each variant.

    Tier 1: OncoKB Oncogenic/Likely Oncogenic
            OR (n_samples >= tier1_min_samples AND n_cancer_types >= tier1_min_cancer_types)
    Tier 2: OncoKB Predicted Oncogenic OR ClinVar/CancerHotspots
            OR (n_samples >= tier2_min_samples AND n_cancer_types >= tier2_min_cancer_types)
    Tier 3: passes base thresholds (minimum confidence)

    Thresholds are read from settings.yaml via cfg['tiering'].
    """
    tiering = (cfg or {}).get('tiering', {})
    t1_s  = tiering.get('tier1_min_samples',       50)
    t1_ct = tiering.get('tier1_min_cancer_types',   3)
    t2_s  = tiering.get('tier2_min_samples',        25)
    t2_ct = tiering.get('tier2_min_cancer_types',    2)

    # Vectorised tiering — avoids iterrows on large dataframes
    onco = df['oncokb_oncogenicity'].fillna('').astype(str)
    n_s  = df['n_samples'].astype(int)
    n_ct = df['n_cancer_types'].astype(int)
    srcs = df['sources'].fillna('').astype(str)

    is_tier1 = (
        onco.isin(['Oncogenic', 'Likely Oncogenic'])
        | (n_s >= t1_s) & (n_ct >= t1_ct)
    )
    is_tier2 = (
        (onco == 'Predicted Oncogenic')
        | srcs.str.contains('ClinVar|CancerHotspots', na=False)
        | (n_s >= t2_s) & (n_ct >= t2_ct)
    )

    df['wl_tier'] = 3
    df.loc[is_tier2 & ~is_tier1, 'wl_tier'] = 2
    df.loc[is_tier1, 'wl_tier'] = 1
    log.info('Tier distribution: Tier1=%d  Tier2=%d  Tier3=%d',
             (df['wl_tier'] == 1).sum(),
             (df['wl_tier'] == 2).sum(),
             (df['wl_tier'] == 3).sum())
    return df


def write_tsv(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    out_cols = [
        'chrom', 'pos', 'ref', 'alt', 'gene', 'hgvsc', 'hgvsp',
        'consequence', 'n_cancer_types', 'cancer_types', 'n_samples',
        'sources', 'oncokb_oncogenicity', 'wl_tier',
    ]
    df_out = df[out_cols].sort_values(['chrom', 'pos'])
    df_out.to_csv(path, sep='\t', index=False,
                  compression='gzip' if path.endswith('.gz') else None)
    log.info('Wrote TSV: %s  (%d rows)', path, len(df_out))


def write_vcf(df: pd.DataFrame, path: str) -> None:
    """Write whitelist as a sites-only VCF (no sample columns)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)

    # Sort by chromosome order then position
    chrom_order = {str(i): i for i in range(1, 23)}
    chrom_order.update({'X': 23, 'Y': 24, 'MT': 25, 'M': 25})
    df = df.copy()
    df['chrom'] = df['chrom'].astype(str).str.replace(r'^chr', '', regex=True).str.replace('chrM', 'MT')
    df['_chrom_sort'] = df['chrom'].map(lambda c: chrom_order.get(c, 99))
    df = df.sort_values(['_chrom_sort', 'pos']).drop(columns='_chrom_sort')

    header = _VCF_HEADER.format(date=datetime.today().strftime('%Y%m%d'))

    import subprocess, tempfile

    # Write to temp uncompressed file then bgzip for tabix compatibility
    tmp_path = path.replace('.gz', '') if path.endswith('.gz') else path + '.tmp'
    with open(tmp_path, 'wt') as fh:
        fh.write(header)
        for row in df.itertuples(index=False):
            info_parts = [
                f'GENE={row.gene}',
                f'CSQ={row.consequence}',
                f'HGVSc={row.hgvsc}',
                f'HGVSp={row.hgvsp}',
                f'N_SAMPLES={row.n_samples}',
                f'N_CANCER_TYPES={row.n_cancer_types}',
                f'CANCER_TYPES={row.cancer_types}',
                f'SOURCES={row.sources}',
                f'WL_TIER={row.wl_tier}',
            ]
            if getattr(row, 'oncokb_oncogenicity', None):
                info_parts.append(f'ONCOKB={row.oncokb_oncogenicity}')

            info = ';'.join(p for p in info_parts if not p.endswith('='))
            fh.write('\t'.join([
                str(row.chrom),
                str(row.pos),
                '.',
                str(row.ref),
                str(row.alt),
                '.',
                'PASS',
                info,
            ]) + '\n')

    if path.endswith('.gz'):
        subprocess.run(['bgzip', '-f', tmp_path], check=True)
        # bgzip writes to tmp_path + '.gz' by default, rename if needed
        bgz_path = tmp_path + '.gz'
        if bgz_path != path:
            os.rename(bgz_path, path)

    log.info('Wrote VCF: %s  (%d variants)', path, len(df))


# ── Entry point ───────────────────────────────────────────────────────────────


def log_database_versions(cfg: dict, settings_file: str = 'settings.yaml') -> None:
    """
    Log the version/date of each database used in the pipeline.
    Writes a human-readable summary to the log and to output/database_versions.txt.
    """
    import re
    import requests
    from datetime import datetime, timezone

    ds      = cfg.get('data_sources', {})
    out_dir = cfg.get('output', {}).get('dir', 'output')
    os.makedirs(out_dir, exist_ok=True)
    lines   = []

    def _add(source, version, path_or_url=''):
        msg = f'  {source:<25} {version}'
        if path_or_url:
            msg += f'  ({path_or_url})'
        log.info(msg)
        lines.append(msg)

    log.info('=' * 60)
    log.info('ONCOSIEVE — pan-cancer variant curation and rescue tool')
    log.info('Author   : Dr Christopher Trethewey')
    log.info('Email    : christopher.trethewey@nhs.net')
    log.info('=' * 60)
    log.info('DATABASE VERSIONS')
    log.info('=' * 60)

    # COSMIC — extract version from filename
    if ds.get('cosmic', {}).get('enabled'):
        tsv = ds['cosmic'].get('tsv', '')
        m = re.search(r'v(\d+)', tsv)
        version = f'v{m.group(1)}' if m else 'unknown'
        _add('COSMIC', version, tsv)

    # GENIE — hardcoded v19
    if ds.get('genie', {}).get('enabled'):
        maf = ds['genie'].get('maf', '')
        _add('GENIE', 'v19', maf)

    # ClinVar — read date from VCF header
    if ds.get('clinvar', {}).get('enabled'):
        vcf = ds['clinvar'].get('vcf', '')
        version = 'unknown'
        if os.path.exists(vcf):
            try:
                import pysam
                v = pysam.VariantFile(vcf)
                for rec in v.header.records:
                    if 'fileDate' in str(rec) or 'dbSNP' in str(rec):
                        version = str(rec).strip()
                        break
                v.close()
            except Exception:
                pass
        _add('ClinVar', version, vcf)

    # OncoKB — query API for data version
    if ds.get('oncokb', {}).get('enabled'):
        version = 'unknown'
        token = cfg.get('data_sources', {}).get('oncokb', {}).get('api_token', '')
        try:
            headers = {'Accept': 'application/json'}
            if token:
                headers['Authorization'] = f'Bearer {token}'
            resp = requests.get(
                'https://www.oncokb.org/api/v1/info',
                headers=headers, timeout=10
            )
            if resp.status_code == 200:
                info = resp.json()
                dv = info.get('dataVersion', {})
                version = f"{dv.get('version','?')} ({dv.get('date','?')})"
        except Exception:
            pass
        _add('OncoKB', version, 'https://www.oncokb.org/api/v1')

    # TP53
    if ds.get('tp53', {}).get('enabled'):
        tsv = ds['tp53'].get('somatic_tsv', '')
        version = 'unknown'
        if os.path.exists(tsv):
            mtime = os.path.getmtime(tsv)
            version = f'downloaded {datetime.fromtimestamp(mtime).strftime("%Y-%m-%d")}'
        _add('TP53 database', version, tsv)

    # cBioPortal — live API, log URL
    if ds.get('cbioportal', {}).get('enabled'):
        version = 'live API'
        try:
            resp = requests.get(
                f"{ds['cbioportal'].get('api_base','https://www.cbioportal.org/api')}/info",
                timeout=10
            )
            if resp.status_code == 200:
                info = resp.json()
                version = f"live API ({info.get('portalVersion', '?')})"
        except Exception:
            pass
        _add('cBioPortal', version, ds['cbioportal'].get('api_base', ''))

    # CancerHotspots — live API
    if ds.get('cancer_hotspots', {}).get('enabled'):
        _add('CancerHotspots', 'live API v2', 'https://www.cancerhotspots.org/api')

    log.info('=' * 60)

    # Write to file
    run_date = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    out_path = os.path.join(out_dir, 'database_versions.txt')
    with open(out_path, 'w') as fh:
        fh.write('ONCOSIEVE — pan-cancer variant curation and rescue tool\n')
        fh.write('Author   : Dr Christopher Trethewey\n')
        fh.write('Email    : christopher.trethewey@nhs.net\n')
        fh.write('=' * 60 + '\n')
        fh.write(f'Run date : {run_date}\n')
        fh.write(f'Settings : {settings_file}\n')
        fh.write('=' * 60 + '\n')
        for line in lines:
            fh.write(line.strip() + '\n')
    log.info('Database version log written to: %s', out_path)

def main():
    parser = argparse.ArgumentParser(
        description='Build pan-cancer whitelist from multiple somatic variant sources.'
    )
    parser.add_argument('--config', default='config.yaml',
                        help='Path to config.yaml (default: config.yaml)')
    parser.add_argument('--skip-sources', default='',
                        help='Comma-separated source names to skip '
                             '(e.g. cosmic,genie)')
    parser.add_argument('--from-intermediates', action='store_true',
                        help='Skip parsers and load from saved intermediate files')
    parser.add_argument('--rerun-oncokb', action='store_true',
                        help='Load all intermediates then re-run OncoKB against merged variants')
    parser.add_argument('--intermediate-only', action='store_true',
                        help='Run parsers and save intermediates; do not merge')
    parser.add_argument('--data-dir', default='',
                        help='Path to reference data directory. Overrides relative '
                             'paths in config.yaml (e.g. /srv/data/reference/)')
    args = parser.parse_args()

    cfg  = load_config(args.config)
    cfg  = _finalise_config(cfg, args.data_dir)
    skip = set(s.strip().lower() for s in args.skip_sources.split(',') if s.strip())

    log.setLevel(getattr(logging, cfg.get('log_level', 'INFO')))

    # Log database versions
    log_database_versions(cfg, settings_file=cfg.get('_settings_file', 'settings.yaml'))

    inter_dir = cfg.get('intermediate_dir', 'intermediate')
    out_dir   = cfg.get('output', {}).get('dir', 'output')
    prefix    = cfg.get('output', {}).get('prefix', 'pan_cancer_whitelist_GRCh38')
    os.makedirs(inter_dir, exist_ok=True)
    os.makedirs(out_dir,   exist_ok=True)

    # Run all parsers (or load from intermediates)
    if args.rerun_oncokb:
        log.info('--rerun-oncokb: loading all intermediates except OncoKB...')
        frames = {}
        for fname in os.listdir(inter_dir):
            if fname.endswith('.tsv.gz') and 'oncokb' not in fname.lower():
                name = fname.replace('.tsv.gz', '')
                path = os.path.join(inter_dir, fname)
                df = pd.read_csv(path, sep='\t', dtype=str, low_memory=False)
                for col in ('n_samples', 'n_cancer_types', 'pos'):
                    if col in df.columns:
                        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
                frames[name] = df
                log.info('Loaded intermediate: %s  (%d rows)', path, len(df))
        merged_so_far = pd.concat([df for df in frames.values() if not df.empty], ignore_index=True)
        log.info('OncoKB re-run: querying against %d merged variants', len(merged_so_far))
        from parsers.parse_oncokb import parse_oncokb
        ds = cfg['data_sources']
        frames['OncoKB'] = parse_oncokb(
            variants_file        = ds['oncokb']['variants_file'],
            merged_df            = merged_so_far,
            include_oncogenicity = ds['oncokb'].get('include_oncogenicity'),
            api_token            = ds['oncokb'].get('api_token'),
        )
        ipath = os.path.join(inter_dir, 'oncokb.tsv.gz')
        frames['OncoKB'].to_csv(ipath, sep='\t', index=False)
        log.info('Saved intermediate: %s  (%d rows)', ipath, len(frames['OncoKB']))
    elif args.from_intermediates:
        log.info('Loading from intermediate files in: %s', inter_dir)
        frames = {}
        for fname in os.listdir(inter_dir):
            if fname.endswith('.tsv.gz'):
                name = fname.replace('.tsv.gz', '')
                path = os.path.join(inter_dir, fname)
                df = pd.read_csv(path, sep='\t', dtype=str, low_memory=False)
                # Restore numeric dtypes lost when saving/loading as TSV
                for col in ('n_samples', 'n_cancer_types', 'pos'):
                    if col in df.columns:
                        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
                frames[name] = df
                log.info('Loaded intermediate: %s  (%d rows)', path, len(df))
    else:
        frames = run_parsers(cfg, skip, inter_dir=inter_dir)

    if args.intermediate_only:
        log.info('--intermediate-only flag set; stopping before merge')
        return

    # Merge and aggregate
    included_csq = set(cfg.get('included_consequences', list(INCLUDED_CONSEQUENCES)))
    df_merged = merge_and_aggregate(frames, included_csq)

    # Apply filters
    df_filtered = apply_filters(df_merged, cfg)

    # Assign tiers
    df_tiered = assign_tiers(df_filtered, cfg=cfg)

    # Write outputs
    tsv_path = os.path.join(out_dir, f'{prefix}.tsv.gz')
    vcf_path = os.path.join(out_dir, f'{prefix}.vcf.gz')
    write_tsv(df_tiered, tsv_path)
    write_vcf(df_tiered, vcf_path)

    log.info('Pipeline complete. Output in: %s', out_dir)
    log.info('  TSV: %s', tsv_path)
    log.info('  VCF: %s', vcf_path)
    log.info('  Index VCF with: tabix -p vcf %s', vcf_path)


if __name__ == '__main__':
    main()
