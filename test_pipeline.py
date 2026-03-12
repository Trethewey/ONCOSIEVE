#!/usr/bin/env python3
"""
ONCOSIEVE — pipeline output test script
Author : Dr Christopher Trethewey
Email  : christopher.trethewey@nhs.net

Tests intermediate files and final output for common issues.
Run from the oncosieve directory:
    python3 test_pipeline.py
"""

import os
import sys
import gzip
import pandas as pd

INTERMEDIATE_DIR = 'intermediate'
OUTPUT_DIR       = 'output'
PREFIX           = 'pan_cancer_whitelist_GRCh38'

PASS = '\033[92mPASS\033[0m'
FAIL = '\033[91mFAIL\033[0m'
WARN = '\033[93mWARN\033[0m'

failures = 0


def check(label: str, condition: bool, detail: str = '', warn_only: bool = False):
    global failures
    if condition:
        print(f'  {PASS}  {label}')
    else:
        tag = WARN if warn_only else FAIL
        print(f'  {tag}  {label}')
        if detail:
            print(f'        {detail}')
        if not warn_only:
            failures += 1


# =============================================================================
# 1. Intermediate files
# =============================================================================
print('\n=== Intermediate files ===')

EXPECTED_INTERMEDIATES = ['cosmic', 'clinvar', 'genie', 'cancerhotspots', 'tp53']

for name in EXPECTED_INTERMEDIATES:
    path = os.path.join(INTERMEDIATE_DIR, f'{name}.tsv.gz')
    check(f'{name}.tsv.gz exists', os.path.exists(path))

# =============================================================================
# 2. Per-intermediate content checks
# =============================================================================
print('\n=== Intermediate content ===')

def load_intermediate(name):
    path = os.path.join(INTERMEDIATE_DIR, f'{name}.tsv.gz')
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, sep='\t', dtype=str, low_memory=False)


# ClinVar — check for tabs in cancer_type
df_clinvar = load_intermediate('clinvar')
if df_clinvar is not None:
    check(f'ClinVar: {len(df_clinvar)} rows', len(df_clinvar) > 0)
    tabs = df_clinvar['cancer_type'].dropna().str.contains('\t').sum() if 'cancer_type' in df_clinvar.columns else 0
    check('ClinVar: no tabs in cancer_type', tabs == 0, detail=f'{tabs} rows contain tabs')

# COSMIC — check cancer types are readable (not coso IDs)
df_cosmic = load_intermediate('cosmic')
if df_cosmic is not None:
    check(f'COSMIC: {len(df_cosmic)} rows', len(df_cosmic) > 0)
    if 'cancer_type' in df_cosmic.columns:
        sample = df_cosmic['cancer_type'].dropna().unique()[:5].tolist()
        coso_count = sum(1 for v in df_cosmic['cancer_type'].dropna() if str(v).startswith('coso'))
        check('COSMIC: cancer_type uses readable names (not coso IDs)',
              coso_count == 0,
              detail=f'{coso_count} rows still have coso IDs. Sample: {sample}')

# GENIE
df_genie = load_intermediate('genie')
if df_genie is not None:
    check(f'GENIE: {len(df_genie)} rows', len(df_genie) > 0)

# TP53
df_tp53 = load_intermediate('tp53')
if df_tp53 is not None:
    check(f'TP53: {len(df_tp53)} rows', len(df_tp53) > 0)

# CancerHotspots
df_hs = load_intermediate('cancerhotspots')
if df_hs is not None:
    check(f'CancerHotspots: {len(df_hs)} rows', len(df_hs) > 0)

# =============================================================================
# 3. Final output files
# =============================================================================
print('\n=== Output files ===')

tsv_path = os.path.join(OUTPUT_DIR, f'{PREFIX}.tsv.gz')
vcf_path = os.path.join(OUTPUT_DIR, f'{PREFIX}.vcf.gz')
tbi_path = vcf_path + '.tbi'
ver_path = os.path.join(OUTPUT_DIR, 'database_versions.txt')

check('TSV output exists',              os.path.exists(tsv_path))
check('VCF output exists',              os.path.exists(vcf_path))
check('VCF tabix index exists',         os.path.exists(tbi_path))
check('database_versions.txt exists',   os.path.exists(ver_path))

# =============================================================================
# 4. TSV content checks
# =============================================================================
print('\n=== TSV content ===')

if os.path.exists(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t', dtype=str, low_memory=False)
    check(f'TSV: {len(df)} rows',           len(df) > 0)
    check('TSV: wl_tier column present',    'wl_tier' in df.columns)
    check('TSV: chrom column present',      'chrom' in df.columns)
    check('TSV: no chr prefix in chrom',
          not df['chrom'].dropna().str.startswith('chr').any(),
          detail='chrom values still have chr prefix')

    # Check for tabs in any field
    tab_cols = [c for c in df.columns if df[c].dropna().str.contains('\t').any()]
    check('TSV: no embedded tabs in any field', len(tab_cols) == 0,
          detail=f'Tabs found in columns: {tab_cols}')

    # Tier distribution
    if 'wl_tier' in df.columns:
        tier_counts = df['wl_tier'].value_counts().to_dict()
        print(f'        Tier distribution: {tier_counts}')

    # Source distribution
    if 'sources' in df.columns:
        from collections import Counter
        source_counts = Counter()
        for val in df['sources'].dropna():
            for s in val.split('|'):
                source_counts[s.strip()] += 1
        print(f'        Source counts: {dict(source_counts.most_common())}')

# =============================================================================
# 5. VCF format checks
# =============================================================================
print('\n=== VCF format ===')

if os.path.exists(vcf_path):
    header_lines = []
    data_lines   = []
    with gzip.open(vcf_path, 'rt') as fh:
        for i, line in enumerate(fh):
            if i > 10000:
                break
            if line.startswith('#'):
                header_lines.append(line.rstrip())
            else:
                data_lines.append(line.rstrip())

    check('VCF: has ##fileformat line',  any('fileformat' in l for l in header_lines))
    check('VCF: has ##contig lines',     any('##contig' in l for l in header_lines))
    check('VCF: has #CHROM header',      any(l.startswith('#CHROM') for l in header_lines))
    check('VCF: has data rows',          len(data_lines) > 0)

    if data_lines:
        # Check field count
        fields = data_lines[0].split('\t')
        check('VCF: 8 fields per row', len(fields) == 8, detail=f'Got {len(fields)} fields')

        # Check for tabs in INFO field
        tab_in_info = sum(1 for l in data_lines if '\t' in l.split('\t', 7)[-1])
        check('VCF: no tabs in INFO field', tab_in_info == 0,
              detail=f'{tab_in_info} rows have tabs in INFO')

        # Check chrom prefix
        chroms = set(l.split('\t')[0] for l in data_lines)
        chr_prefixed = sum(1 for c in chroms if c.startswith('chr'))
        check('VCF: no chr prefix in CHROM',  chr_prefixed == 0,
              detail=f'chr-prefixed chroms: {[c for c in chroms if c.startswith("chr")][:5]}')

# =============================================================================
# Summary
# =============================================================================
print(f'\n{"="*50}')
if failures == 0:
    print(f'  {PASS}  All checks passed.')
else:
    print(f'  {FAIL}  {failures} check(s) failed.')
print('='*50 + '\n')

sys.exit(0 if failures == 0 else 1)
