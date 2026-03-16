#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — validate_transcript_fix.py
#
# Comprehensive validation of transcript priority fix.
# Checks pipeline outputs against intermediate source files and known
# reference variants to confirm:
#
#   1. Transcript priority ordering is correct
#   2. hgvsc / hgvsp are fully populated (full ENST...:c. format)
#   3. transcript_source and is_mane_select columns are consistent
#   4. Key clinical variants have correct MANE Select transcripts
#   5. Row counts and column schema are correct
#   6. Determinism — no non-deterministic transcript selection
#   7. No regressions in variant counts, tiers, or source coverage
#
# Usage:
#   python validate_transcript_fix.py \
#       --output output/pan_cancer_whitelist_GRCh38_full.tsv.gz \
#       --base-tsv output/pan_cancer_whitelist_GRCh38.tsv.gz \
#       --intermediate-dir intermediate \
#       [--reference-output output/excellent/pan_cancer_whitelist_GRCh38_full.tsv.gz]
#
# Author : Dr Christopher Trethewey
# =============================================================================

import argparse
import os
import re
import sys
from collections import Counter

import pandas as pd

# ── Constants ────────────────────────────────────────────────────────────────

TRANSCRIPT_SOURCE_PRIORITY = {
    'CancerHotspots': 0,
    'GENIE':          1,
    'ClinVar':        2,
    'COSMIC':         3,
    'TCGA':           4,
    'TP53_somatic':   5,
    'TP53_germline':  5,
}

# Key clinical variants that MUST have correct transcripts
# These are the most clinically important variants in oncology
# Format: (chrom, pos, ref, alt, gene, expected_protein_residue_keyword)
KEY_VARIANTS = [
    ('12', 25245350, 'C', 'T', 'KRAS',   'Gly12'),    # KRAS G12D
    ('12', 25245350, 'C', 'A', 'KRAS',   'Gly12'),    # KRAS G12V
    ('12', 25245351, 'C', 'A', 'KRAS',   'Gly12'),    # KRAS G12C
    ('7',  140753336,'A', 'T', 'BRAF',   'Val600'),    # BRAF V600E
    ('3',  179218303,'G', 'A', 'PIK3CA', 'Glu545'),    # PIK3CA E545K
    ('3',  179234297,'A', 'G', 'PIK3CA', 'His1047'),   # PIK3CA H1047R
    ('17', 7675088,  'C', 'T', 'TP53',   'Arg175'),    # TP53 R175H
    ('17', 7674220,  'C', 'T', 'TP53',   'Arg248'),    # TP53 R248Q
    ('17', 7673803,  'G', 'A', 'TP53',   'Arg273'),    # TP53 R273C
    ('14', 104780214,'C', 'T', 'AKT1',   'Glu17'),     # AKT1 E17K
]

# Expected output columns for base TSV (build_whitelist.py output)
EXPECTED_BASE_COLS = [
    'chrom', 'pos', 'ref', 'alt', 'gene', 'hgvsc', 'hgvsp',
    'consequence', 'n_cancer_types', 'cancer_types', 'n_samples',
    'sources', 'oncokb_oncogenicity', 'clinvar_clinical_significance',
    'transcript_source', 'is_mane_select',
    'tp53_class', 'wl_tier',
]

# Expected output columns for full TSV (after post_process + post_pipeline)
EXPECTED_FULL_COLS = [
    'chrom', 'pos', 'ref', 'alt', 'gene',
    'transcript_id', 'is_mane_select', 'hgvsc', 'hgvsp',
    'protein_change', 'consequence',
    'n_cancer_types', 'cancer_types', 'n_samples', 'sources',
    'oncokb_oncogenicity', 'clinvar_clinical_significance',
    'transcript_source', 'tp53_class', 'wl_tier',
    'genome_version', 'revel_score',
]

INTERMEDIATE_SOURCES = {
    'cancerhotspots': 'CancerHotspots',
    'genie':          'GENIE',
    'clinvar':        'ClinVar',
    'cosmic':         'COSMIC',
    'tcga':           'TCGA',
    'tp53':           'TP53_somatic',
    'oncokb':         'OncoKB',
}

_RE_ENST = re.compile(r'ENST\d+', re.I)

# ── Test helpers ─────────────────────────────────────────────────────────────

class ValidationResult:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.warnings = 0
        self.messages = []

    def ok(self, msg):
        self.passed += 1
        self.messages.append(f'  PASS  {msg}')

    def fail(self, msg):
        self.failed += 1
        self.messages.append(f'  FAIL  {msg}')

    def warn(self, msg):
        self.warnings += 1
        self.messages.append(f'  WARN  {msg}')

    def info(self, msg):
        self.messages.append(f'  INFO  {msg}')

    def section(self, title):
        self.messages.append(f'\n{"=" * 70}')
        self.messages.append(f'  {title}')
        self.messages.append(f'{"=" * 70}')

    def report(self):
        for m in self.messages:
            print(m)
        print()
        print('=' * 70)
        print(f'  SUMMARY: {self.passed} passed, {self.failed} failed, '
              f'{self.warnings} warnings')
        print('=' * 70)
        return self.failed == 0


# ── Test 1: Column schema ───────────────────────────────────────────────────

def test_column_schema(r: ValidationResult, df: pd.DataFrame, label: str,
                       expected_cols: list):
    r.section(f'Column schema — {label}')

    for col in expected_cols:
        if col in df.columns:
            r.ok(f'Column present: {col}')
        else:
            r.fail(f'Missing column: {col}')

    # Check for unexpected columns (not necessarily a failure)
    extra = [c for c in df.columns if c not in expected_cols]
    if extra:
        r.info(f'Extra columns (not in expected list): {extra}')


# ── Test 2: New columns populated ────────────────────────────────────────────

def test_new_columns_populated(r: ValidationResult, df: pd.DataFrame):
    r.section('New columns populated')

    # transcript_source
    if 'transcript_source' in df.columns:
        ts = df['transcript_source'].fillna('')
        n_populated = (ts != '').sum()
        n_total = len(df)
        pct = 100 * n_populated / n_total if n_total else 0
        r.info(f'transcript_source populated: {n_populated:,} / {n_total:,} '
               f'({pct:.1f}%)')
        if n_populated > 0:
            r.ok('transcript_source has values')
        else:
            r.fail('transcript_source is entirely empty')

        # Check values are valid source names
        valid_sources = set(TRANSCRIPT_SOURCE_PRIORITY.keys()) | {''}
        invalid = set(ts.unique()) - valid_sources
        if invalid:
            r.fail(f'Invalid transcript_source values: {invalid}')
        else:
            r.ok('All transcript_source values are valid source names')

        # Distribution
        dist = ts[ts != ''].value_counts()
        for src, count in dist.items():
            r.info(f'  transcript_source = {src}: {count:,}')
    else:
        r.fail('transcript_source column missing')

    # is_mane_select
    if 'is_mane_select' in df.columns:
        mane = df['is_mane_select']
        n_true = mane.astype(str).str.lower().isin(['true', '1']).sum()
        r.info(f'is_mane_select = True: {n_true:,}')
        if n_true > 0:
            r.ok(f'is_mane_select has {n_true:,} True values')
        else:
            r.warn('is_mane_select has no True values '
                    '(expected if CancerHotspots intermediate is present)')
    else:
        r.fail('is_mane_select column missing')


# ── Test 3: is_mane_select consistency ───────────────────────────────────────

def test_mane_consistency(r: ValidationResult, df: pd.DataFrame):
    r.section('is_mane_select consistency with transcript_source')

    if 'transcript_source' not in df.columns or 'is_mane_select' not in df.columns:
        r.fail('Required columns missing for MANE consistency check')
        return

    ts = df['transcript_source'].fillna('')
    mane = df['is_mane_select'].astype(str).str.lower()

    # is_mane_select should be True iff transcript_source == CancerHotspots
    is_hotspot = (ts == 'CancerHotspots')
    is_mane_true = mane.isin(['true', '1'])

    # Hotspot source but not flagged as MANE
    mismatch_a = is_hotspot & ~is_mane_true
    if mismatch_a.sum() > 0:
        r.fail(f'{mismatch_a.sum()} rows have transcript_source=CancerHotspots '
               f'but is_mane_select is not True')
    else:
        r.ok('All CancerHotspots-sourced transcripts are flagged MANE Select')

    # Flagged as MANE but not from CancerHotspots
    mismatch_b = ~is_hotspot & is_mane_true
    if mismatch_b.sum() > 0:
        r.fail(f'{mismatch_b.sum()} rows have is_mane_select=True '
               f'but transcript_source is not CancerHotspots')
    else:
        r.ok('No non-CancerHotspots rows falsely flagged as MANE Select')


# ── Test 4: hgvsc fully populated ────────────────────────────────────────────

def test_hgvsc_format(r: ValidationResult, df: pd.DataFrame):
    r.section('hgvsc format — full transcript:c. notation')

    hgvsc = df['hgvsc'].fillna('')
    has_hgvsc = hgvsc[hgvsc != '']
    n_total = len(df)
    n_has = len(has_hgvsc)
    r.info(f'Rows with hgvsc: {n_has:,} / {n_total:,}')

    if n_has == 0:
        r.fail('No rows have hgvsc values')
        return

    # Check how many have full ENST prefix
    has_enst = has_hgvsc.str.contains(r'ENST\d+', case=False, na=False)
    n_enst = has_enst.sum()
    n_bare_c = (~has_enst).sum()
    pct_enst = 100 * n_enst / n_has if n_has else 0
    r.info(f'  With ENST prefix:  {n_enst:,} ({pct_enst:.1f}%)')
    r.info(f'  Bare c. notation:  {n_bare_c:,} '
           f'(expected for TCGA/TP53 sourced variants)')

    if pct_enst >= 80:
        r.ok(f'{pct_enst:.1f}% of hgvsc values have full transcript prefix')
    elif pct_enst >= 50:
        r.warn(f'Only {pct_enst:.1f}% of hgvsc values have transcript prefix')
    else:
        r.fail(f'Only {pct_enst:.1f}% of hgvsc values have transcript prefix')

    # Check hgvsc is NOT stripped to bare c. when ENST is available
    # (Regression check: post_pipeline.py used to strip ENST prefix)
    stripped = has_hgvsc.str.match(r'^c\.', na=False)
    n_stripped = stripped.sum()
    r.info(f'  Starting with bare "c.":  {n_stripped:,}')

    # Verify colon-separated format where ENST present
    has_colon = has_hgvsc[has_enst].str.contains(':', na=False)
    n_colon = has_colon.sum()
    if n_enst > 0:
        if n_colon == n_enst:
            r.ok('All ENST-prefixed hgvsc values use ENST:c. format')
        else:
            r.warn(f'{n_enst - n_colon} ENST-prefixed values missing colon separator')


# ── Test 5: hgvsp fully populated ────────────────────────────────────────────

def test_hgvsp_format(r: ValidationResult, df: pd.DataFrame):
    r.section('hgvsp format')

    hgvsp = df['hgvsp'].fillna('')
    has_hgvsp = hgvsp[hgvsp != '']
    n_total = len(df)
    n_has = len(has_hgvsp)
    r.info(f'Rows with hgvsp: {n_has:,} / {n_total:,}')

    if n_has == 0:
        r.warn('No rows have hgvsp values')
        return

    # Check for ENSP prefix (full format)
    has_ensp = has_hgvsp.str.contains(r'ENSP\d+', case=False, na=False)
    n_ensp = has_ensp.sum()
    pct_ensp = 100 * n_ensp / n_has if n_has else 0
    r.info(f'  With ENSP prefix: {n_ensp:,} ({pct_ensp:.1f}%)')

    # Check for p. prefix
    has_p = has_hgvsp.str.contains(r'p\.', na=False)
    n_p = has_p.sum()
    r.info(f'  With p. notation: {n_p:,}')


# ── Test 6: Transcript priority ordering ─────────────────────────────────────

def test_transcript_priority(r: ValidationResult, df: pd.DataFrame,
                              inter_dir: str):
    r.section('Transcript priority ordering')

    if not os.path.isdir(inter_dir):
        r.warn(f'Intermediate directory not found: {inter_dir}')
        return

    # Load all intermediate files
    sources = {}
    for fname in os.listdir(inter_dir):
        if not fname.endswith('.tsv.gz'):
            continue
        key = fname.replace('.tsv.gz', '')
        source_name = INTERMEDIATE_SOURCES.get(key, key)
        try:
            idf = pd.read_csv(os.path.join(inter_dir, fname), sep='\t',
                               dtype=str, low_memory=False)
            sources[source_name] = idf
        except Exception as e:
            r.warn(f'Could not load {fname}: {e}')

    if not sources:
        r.fail('No intermediate files loaded')
        return

    r.info(f'Loaded {len(sources)} intermediate files: {list(sources.keys())}')

    # For each variant in the output that has a transcript_source,
    # verify that no higher-priority source had a valid hgvsc for that variant
    if 'transcript_source' not in df.columns:
        r.fail('transcript_source column missing — cannot verify priority')
        return

    # Build lookup: (chrom, pos, ref, alt) -> {source: hgvsc}
    # This is expensive, so sample a subset
    sample = df[
        df['transcript_source'].fillna('').ne('') &
        df['hgvsc'].fillna('').ne('')
    ].head(5000)

    violations = 0
    checked = 0

    for _, row in sample.iterrows():
        chrom = str(row['chrom'])
        pos = str(row['pos'])
        ref = str(row['ref'])
        alt = str(row['alt'])
        winning_source = row['transcript_source']
        winning_priority = TRANSCRIPT_SOURCE_PRIORITY.get(winning_source, 99)

        # Check if any higher-priority source had hgvsc for this variant
        for src_name, idf in sources.items():
            src_priority = TRANSCRIPT_SOURCE_PRIORITY.get(src_name, 99)
            if src_priority >= winning_priority:
                continue  # Same or lower priority, skip

            match = idf[
                (idf['chrom'].astype(str) == chrom) &
                (idf['pos'].astype(str) == pos) &
                (idf['ref'].astype(str) == ref) &
                (idf['alt'].astype(str) == alt)
            ]
            if match.empty:
                continue

            # Check if this higher-priority source had a valid hgvsc
            hgvsc_vals = match['hgvsc'].dropna()
            hgvsc_vals = hgvsc_vals[
                ~hgvsc_vals.str.strip().str.lower().isin(['', 'nan', 'none', 'na'])
            ]
            if not hgvsc_vals.empty:
                violations += 1
                if violations <= 5:
                    r.fail(
                        f'Priority violation at {chrom}:{pos} {ref}>{alt}: '
                        f'winning source={winning_source} (priority={winning_priority}) '
                        f'but {src_name} (priority={src_priority}) also had hgvsc='
                        f'{hgvsc_vals.iloc[0]}'
                    )
                break

        checked += 1

    if violations == 0:
        r.ok(f'No priority violations found in {checked:,} sampled variants')
    else:
        r.fail(f'{violations} priority violations in {checked:,} sampled variants')


# ── Test 7: Key clinical variants ────────────────────────────────────────────

def test_key_variants(r: ValidationResult, df: pd.DataFrame):
    r.section('Key clinical variant annotations')

    df_str = df.copy()
    df_str['chrom'] = df_str['chrom'].astype(str)
    df_str['pos'] = pd.to_numeric(df_str['pos'], errors='coerce').astype('Int64')

    for chrom, pos, ref, alt, gene, expected_residue in KEY_VARIANTS:
        match = df_str[
            (df_str['chrom'] == chrom) &
            (df_str['pos'] == pos) &
            (df_str['ref'] == ref) &
            (df_str['alt'] == alt)
        ]

        label = f'{gene} {chrom}:{pos} {ref}>{alt}'

        if match.empty:
            r.fail(f'{label}: NOT FOUND in output')
            continue

        row = match.iloc[0]
        hgvsc = str(row.get('hgvsc', ''))
        hgvsp = str(row.get('hgvsp', ''))
        protein_change = str(row.get('protein_change', ''))
        transcript_source = str(row.get('transcript_source', ''))
        is_mane = str(row.get('is_mane_select', ''))
        tid = str(row.get('transcript_id', ''))

        r.info(f'{label}:')
        r.info(f'  transcript_id={tid}  transcript_source={transcript_source}  '
               f'is_mane={is_mane}')
        r.info(f'  hgvsc={hgvsc}')
        r.info(f'  hgvsp={hgvsp}  protein_change={protein_change}')

        # Check protein change contains expected residue
        prot_check = protein_change if protein_change not in ('', 'nan') else hgvsp
        if expected_residue in str(prot_check):
            r.ok(f'{label}: protein annotation contains {expected_residue}')
        else:
            r.fail(f'{label}: expected {expected_residue} in protein annotation, '
                   f'got {prot_check}')

        # Check hgvsc is not empty (these are well-known variants, should all have it)
        if hgvsc and hgvsc != 'nan':
            r.ok(f'{label}: hgvsc populated')
        else:
            r.warn(f'{label}: hgvsc is empty')


# ── Test 8: Row count and tier distribution ──────────────────────────────────

def test_row_counts(r: ValidationResult, df: pd.DataFrame, ref_df=None):
    r.section('Row counts and tier distribution')

    n = len(df)
    r.info(f'Total variants: {n:,}')

    if 'wl_tier' in df.columns:
        tiers = df['wl_tier'].astype(str).value_counts().sort_index()
        for tier, count in tiers.items():
            r.info(f'  Tier {tier}: {count:,}')

    if ref_df is not None:
        n_ref = len(ref_df)
        r.info(f'Reference output variants: {n_ref:,}')

        if n == n_ref:
            r.ok(f'Row count matches reference: {n:,}')
        elif abs(n - n_ref) <= 10:
            r.warn(f'Row count differs slightly from reference: '
                   f'{n:,} vs {n_ref:,} (diff={n - n_ref})')
        else:
            r.fail(f'Row count differs significantly from reference: '
                   f'{n:,} vs {n_ref:,} (diff={n - n_ref})')


# ── Test 9: Source coverage ──────────────────────────────────────────────────

def test_source_coverage(r: ValidationResult, df: pd.DataFrame, ref_df=None):
    r.section('Source coverage')

    if 'sources' not in df.columns:
        r.fail('sources column missing')
        return

    all_sources = df['sources'].fillna('').str.split('|').explode().str.strip()
    source_counts = all_sources[all_sources != ''].value_counts()

    for src, count in source_counts.items():
        r.info(f'  {src}: {count:,} variant-source pairs')

    expected_sources = {'COSMIC', 'GENIE', 'TCGA', 'ClinVar', 'OncoKB',
                        'CancerHotspots', 'TP53_somatic'}
    present = set(source_counts.index)
    missing = expected_sources - present
    if missing:
        r.warn(f'Expected sources not found in output: {missing}')
    else:
        r.ok('All expected sources present')

    if ref_df is not None and 'sources' in ref_df.columns:
        ref_sources = ref_df['sources'].fillna('').str.split('|').explode().str.strip()
        ref_counts = ref_sources[ref_sources != ''].value_counts()
        for src in expected_sources & present:
            new_count = source_counts.get(src, 0)
            ref_count = ref_counts.get(src, 0)
            if new_count == ref_count:
                r.ok(f'{src}: count matches reference ({new_count:,})')
            elif abs(new_count - ref_count) <= 5:
                r.ok(f'{src}: count within tolerance '
                     f'({new_count:,} vs ref {ref_count:,})')
            else:
                r.warn(f'{src}: count differs from reference: '
                       f'{new_count:,} vs {ref_count:,}')


# ── Test 10: Transcript ID consistency with hgvsc ────────────────────────────

def test_transcript_id_consistency(r: ValidationResult, df: pd.DataFrame):
    r.section('transcript_id consistency with hgvsc')

    if 'transcript_id' not in df.columns:
        r.warn('transcript_id column not present (may be added by post-processing)')
        return

    # Where hgvsc has ENST prefix, transcript_id should match
    hgvsc = df['hgvsc'].fillna('')
    tid = df['transcript_id'].fillna('')

    has_enst_in_hgvsc = hgvsc.str.contains(r'ENST\d+', case=False, na=False)
    subset = df[has_enst_in_hgvsc].copy()

    if subset.empty:
        r.warn('No hgvsc values with ENST prefix to check')
        return

    # Extract ENST from hgvsc
    extracted = subset['hgvsc'].str.extract(r'(ENST\d+)', flags=re.I)[0].str.upper()
    tid_sub = subset['transcript_id'].fillna('').str.upper()

    matches = (extracted == tid_sub)
    n_match = matches.sum()
    n_total = len(subset)

    if n_match == n_total:
        r.ok(f'transcript_id matches hgvsc ENST in all {n_total:,} rows')
    else:
        n_mismatch = n_total - n_match
        r.fail(f'{n_mismatch} rows have transcript_id/hgvsc ENST mismatch')
        mismatches = subset[~matches][['gene', 'hgvsc', 'transcript_id']].head(5)
        for _, m in mismatches.iterrows():
            r.info(f'  Mismatch: hgvsc={m["hgvsc"]}, transcript_id={m["transcript_id"]}')


# ── Test 11: Cross-reference against intermediates ───────────────────────────

def test_cross_reference_intermediates(r: ValidationResult, df: pd.DataFrame,
                                        inter_dir: str):
    r.section('Cross-reference: CancerHotspots transcripts in output')

    hotspots_path = os.path.join(inter_dir, 'cancerhotspots.tsv.gz')
    if not os.path.exists(hotspots_path):
        r.warn('CancerHotspots intermediate not found — skipping')
        return

    hs = pd.read_csv(hotspots_path, sep='\t', dtype=str, low_memory=False)
    hs = hs[hs['hgvsc'].notna() & (hs['hgvsc'].str.strip() != '')]
    r.info(f'CancerHotspots variants with hgvsc: {len(hs):,}')

    # For each CancerHotspots variant with hgvsc, check if the output uses it
    hs['chrom'] = hs['chrom'].astype(str).str.replace(r'^chr', '', regex=True)
    df_check = df.copy()
    df_check['chrom'] = df_check['chrom'].astype(str)

    merged = hs.merge(
        df_check[['chrom', 'pos', 'ref', 'alt', 'hgvsc', 'transcript_source']],
        on=['chrom', 'pos', 'ref', 'alt'],
        how='inner',
        suffixes=('_hs', '_out')
    )

    if merged.empty:
        r.warn('No CancerHotspots variants found in output (coordinate match)')
        return

    r.info(f'CancerHotspots variants found in output: {len(merged):,}')

    # Check how many use CancerHotspots as transcript_source
    uses_hs = (merged['transcript_source'] == 'CancerHotspots').sum()
    r.info(f'  Using CancerHotspots transcript: {uses_hs:,}')

    if uses_hs == len(merged):
        r.ok('All CancerHotspots variants use CancerHotspots transcript (priority 0)')
    elif uses_hs > 0:
        not_using = len(merged) - uses_hs
        r.warn(f'{not_using} CancerHotspots variants NOT using CancerHotspots '
               f'transcript (may indicate hgvsc was empty in CancerHotspots)')
        # Show examples
        examples = merged[merged['transcript_source'] != 'CancerHotspots'].head(3)
        for _, ex in examples.iterrows():
            r.info(f'  {ex["gene"]} {ex["chrom"]}:{ex["pos"]}: '
                   f'hs_hgvsc={ex["hgvsc_hs"]}, out_hgvsc={ex["hgvsc_out"]}, '
                   f'source={ex["transcript_source"]}')
    else:
        r.fail('No CancerHotspots variants use CancerHotspots transcript')


# ── Test 12: Determinism check ───────────────────────────────────────────────

def test_determinism_spot_check(r: ValidationResult, df: pd.DataFrame):
    r.section('Determinism spot check')

    # Check that within the output, same gene doesn't randomly get different
    # transcript_ids for same genomic position (which would indicate
    # non-determinism at the gene level)
    if 'transcript_source' not in df.columns:
        r.warn('transcript_source not available for determinism check')
        return

    # Group by (chrom, pos, ref, alt) — should have exactly 1 row each
    dups = df.groupby(['chrom', 'pos', 'ref', 'alt']).size()
    n_dups = (dups > 1).sum()
    if n_dups > 0:
        r.warn(f'{n_dups} coordinate keys appear more than once')
    else:
        r.ok('All coordinate keys are unique (no duplicate variants)')

    # Check that KRAS variants at same locus all use same transcript
    kras = df[df['gene'] == 'KRAS'].copy()
    if not kras.empty and 'transcript_id' in df.columns:
        kras_tids = kras.groupby('pos')['transcript_id'].nunique()
        multi_tid = kras_tids[kras_tids > 1]
        if len(multi_tid) > 0:
            # This is expected — different alt alleles at same pos can get
            # different transcripts. Only flag if transcript_source differs.
            r.info(f'KRAS: {len(multi_tid)} positions with multiple transcript_ids '
                   f'(expected if different alts)')
        else:
            r.ok('KRAS: consistent transcript_id per position')


# ── Test 13: Comparison with reference output ────────────────────────────────

def test_reference_comparison(r: ValidationResult, df: pd.DataFrame,
                               ref_df: pd.DataFrame):
    r.section('Comparison with reference (pre-fix) output')

    # Variants should be the same set (same chrom/pos/ref/alt)
    df_keys = set(zip(df['chrom'].astype(str), df['pos'].astype(str),
                      df['ref'].astype(str), df['alt'].astype(str)))
    ref_keys = set(zip(ref_df['chrom'].astype(str), ref_df['pos'].astype(str),
                       ref_df['ref'].astype(str), ref_df['alt'].astype(str)))

    in_both = df_keys & ref_keys
    only_new = df_keys - ref_keys
    only_ref = ref_keys - df_keys

    r.info(f'Variants in both:     {len(in_both):,}')
    r.info(f'Only in new output:   {len(only_new):,}')
    r.info(f'Only in reference:    {len(only_ref):,}')

    if len(only_new) == 0 and len(only_ref) == 0:
        r.ok('Variant sets are identical')
    elif len(only_new) <= 5 and len(only_ref) <= 5:
        r.warn('Minor differences in variant sets')
    else:
        r.fail(f'Significant difference in variant sets: '
               f'+{len(only_new)} / -{len(only_ref)}')

    # Compare transcript changes for key genes
    r.info('')
    r.info('Transcript changes for key genes:')
    for gene in ['KRAS', 'BRAF', 'PIK3CA', 'TP53', 'AKT1']:
        new_gene = df[df['gene'] == gene].head(5)
        ref_gene = ref_df[ref_df['gene'] == gene].head(5)

        if new_gene.empty:
            continue

        for _, nrow in new_gene.iterrows():
            key = (str(nrow['chrom']), str(nrow['pos']),
                   str(nrow['ref']), str(nrow['alt']))
            rmatches = ref_df[
                (ref_df['chrom'].astype(str) == key[0]) &
                (ref_df['pos'].astype(str) == key[1]) &
                (ref_df['ref'].astype(str) == key[2]) &
                (ref_df['alt'].astype(str) == key[3])
            ]
            if rmatches.empty:
                continue

            rrow = rmatches.iloc[0]
            new_hgvsc = str(nrow.get('hgvsc', ''))
            ref_hgvsc = str(rrow.get('hgvsc', ''))
            new_tsrc = str(nrow.get('transcript_source', ''))

            if new_hgvsc != ref_hgvsc:
                r.info(f'  {gene} {key[0]}:{key[1]} {key[2]}>{key[3]}: '
                       f'CHANGED from {ref_hgvsc[:50]} -> {new_hgvsc[:50]} '
                       f'(source={new_tsrc})')
            else:
                r.info(f'  {gene} {key[0]}:{key[1]}: unchanged ({new_hgvsc[:50]})')


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='ONCOSIEVE — validate transcript priority fix'
    )
    parser.add_argument('--output', required=True,
                        help='Path to output TSV (full or base)')
    parser.add_argument('--base-tsv', default='',
                        help='Path to base TSV from build_whitelist.py '
                             '(pre-post_process, for schema check)')
    parser.add_argument('--intermediate-dir', default='intermediate',
                        help='Path to intermediate directory')
    parser.add_argument('--reference-output', default='',
                        help='Path to reference (pre-fix) output TSV '
                             'for comparison')
    args = parser.parse_args()

    r = ValidationResult()

    print()
    print('=' * 70)
    print('  ONCOSIEVE — Transcript Priority Fix Validation')
    print('=' * 70)
    print(f'  Output:         {args.output}')
    print(f'  Intermediates:  {args.intermediate_dir}')
    if args.base_tsv:
        print(f'  Base TSV:       {args.base_tsv}')
    if args.reference_output:
        print(f'  Reference:      {args.reference_output}')

    # Load output
    if not os.path.exists(args.output):
        print(f'\nERROR: Output file not found: {args.output}')
        sys.exit(1)

    print(f'\nLoading output: {args.output}')
    df = pd.read_csv(args.output, sep='\t', dtype=str, low_memory=False)
    print(f'  {len(df):,} rows, {len(df.columns)} columns')

    # Load base TSV if provided
    base_df = None
    if args.base_tsv and os.path.exists(args.base_tsv):
        print(f'Loading base TSV: {args.base_tsv}')
        base_df = pd.read_csv(args.base_tsv, sep='\t', dtype=str,
                               low_memory=False)
        print(f'  {len(base_df):,} rows, {len(base_df.columns)} columns')

    # Load reference output if provided
    ref_df = None
    if args.reference_output and os.path.exists(args.reference_output):
        print(f'Loading reference: {args.reference_output}')
        ref_df = pd.read_csv(args.reference_output, sep='\t', dtype=str,
                              low_memory=False)
        print(f'  {len(ref_df):,} rows, {len(ref_df.columns)} columns')

    # Run all tests
    # Test 1: Column schema
    is_full = 'revel_score' in df.columns  # full TSV has revel_score
    if is_full:
        test_column_schema(r, df, 'full TSV', EXPECTED_FULL_COLS)
    if base_df is not None:
        test_column_schema(r, base_df, 'base TSV', EXPECTED_BASE_COLS)

    # Test 2: New columns populated
    # Use whichever has transcript_source
    check_df = base_df if (base_df is not None and
                            'transcript_source' in base_df.columns) else df
    test_new_columns_populated(r, check_df)

    # Test 3: MANE consistency
    test_mane_consistency(r, check_df)

    # Test 4: hgvsc format
    test_hgvsc_format(r, df)

    # Test 5: hgvsp format
    test_hgvsp_format(r, df)

    # Test 6: Transcript priority ordering
    test_transcript_priority(r, check_df, args.intermediate_dir)

    # Test 7: Key clinical variants
    test_key_variants(r, df)

    # Test 8: Row counts
    test_row_counts(r, df, ref_df)

    # Test 9: Source coverage
    test_source_coverage(r, df, ref_df)

    # Test 10: transcript_id consistency
    test_transcript_id_consistency(r, df)

    # Test 11: Cross-reference intermediates
    test_cross_reference_intermediates(r, check_df, args.intermediate_dir)

    # Test 12: Determinism
    test_determinism_spot_check(r, df)

    # Test 13: Reference comparison
    if ref_df is not None:
        test_reference_comparison(r, df, ref_df)

    # Final report
    print()
    success = r.report()
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
