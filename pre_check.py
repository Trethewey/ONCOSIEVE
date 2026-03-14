#!/usr/bin/env python3
"""
ONCOSIEVE — pre-run audit script
Author : Dr Christopher Trethewey
Email  : christopher.trethewey@nhs.net

Detects and reports the state of all input files before the pipeline runs.
Only fails on conditions that will definitely prevent the pipeline from running.
Column name detection is advisory — parsers handle their own adaptive mapping.

Run from the oncosieve directory:
    python3 pre_check.py
    python3 pre_check.py --config config.yaml

Exit code 0 = safe to run.
Exit code 1 = hard failure (missing file, wrong compression, missing required index).
"""

import os
import sys
import gzip
import argparse
import subprocess
import yaml

PASS  = '\033[92mPASS\033[0m'
FAIL  = '\033[91mFAIL\033[0m'
WARN  = '\033[93mWARN\033[0m'
INFO  = '\033[94mINFO\033[0m'
SKIP  = '\033[90mSKIP\033[0m'

failures = 0
warnings = 0


def fail(label, detail=''):
    global failures
    print(f'  {FAIL}  {label}')
    if detail:
        print(f'        {detail}')
    failures += 1


def warn(label, detail=''):
    global warnings
    print(f'  {WARN}  {label}')
    if detail:
        print(f'        {detail}')
    warnings += 1


def ok(label):
    print(f'  {PASS}  {label}')


def info(label, detail=''):
    print(f'  {INFO}  {label}')
    if detail:
        print(f'        {detail}')


def skip(label, reason=''):
    print(f'  {SKIP}  {label}' + (f'  ({reason})' if reason else ''))


def section(title):
    print(f'\n=== {title} ===')


def file_exists(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


def detect_columns(path):
    """Detect column names from first non-comment line. Auto-detects tab/comma."""
    opener = gzip.open if path.endswith('.gz') else open
    try:
        with opener(path, 'rt') as fh:
            for _ in range(500):
                line = fh.readline()
                if not line:
                    break
                s = line.strip()
                if not s or s.startswith('#'):
                    continue
                sep = '\t' if s.count('\t') >= s.count(',') else ','
                return s.split(sep)
    except Exception:
        return []
    return []


# Known chromosome column names across all supported sources
_CHROM_COLNAMES = (
    'CHROMOSOME',       # COSMIC
    'Chromosome',       # GENIE MAF
    'CHROM',            # VCF
    '#CHROM',           # VCF with hash
    'Chr',              # TP53
    'chrom',            # generic
    'chromosome',       # generic lowercase
)


def detect_chrom_style(path, n=2000):
    """
    Locate the chromosome column by name, then sample values to determine
    naming style. Falls back to column 0 if no known name is found.
    """
    opener = gzip.open if path.endswith('.gz') else open
    chr_count = plain_count = seen = 0
    sample = set()
    chrom_col_idx = None
    sep = '\t'

    try:
        with opener(path, 'rt') as fh:
            # Find header line
            for line in fh:
                s = line.strip()
                if not s or s.startswith('##'):
                    continue
                # Detect separator
                sep = '\t' if s.count('\t') >= s.count(',') else ','
                cols = s.split(sep)
                # Try to find a known chromosome column
                for name in _CHROM_COLNAMES:
                    if name in cols:
                        chrom_col_idx = cols.index(name)
                        break
                break  # header found

            if chrom_col_idx is None:
                return {'style': 'unknown', 'sample': [],
                        'note': 'No known chromosome column found in header'}

            # Sample data rows
            for line in fh:
                s = line.strip()
                if not s or s.startswith('#'):
                    continue
                parts = s.split(sep)
                if chrom_col_idx >= len(parts):
                    continue
                chrom = parts[chrom_col_idx].strip()
                if not chrom:
                    continue
                sample.add(chrom)
                if chrom.startswith('chr'):
                    chr_count += 1
                elif chrom[:1].isdigit() or chrom in ('X', 'Y', 'MT', 'M'):
                    plain_count += 1
                seen += 1
                if seen >= n:
                    break

    except Exception:
        return {'style': 'unknown', 'sample': []}

    if seen == 0:
        style = 'unknown'
    elif chr_count > 0 and plain_count == 0:
        style = 'chr-prefixed'
    elif plain_count > 0 and chr_count == 0:
        style = 'plain'
    else:
        style = 'mixed'
    return {'style': style, 'sample': sorted(sample)[:6]}


def is_bgzf(path):
    try:
        with open(path, 'rb') as fh:
            return fh.read(4) == b'\x1f\x8b\x08\x04'
    except Exception:
        return False


def has_index(path, extensions):
    return any(os.path.isfile(path + e) for e in extensions)


def main():
    parser = argparse.ArgumentParser(description='ONCOSIEVE pre-run audit')
    parser.add_argument('--config', default='config.yaml')
    parser.add_argument('--data-dir', default='',
                        help='Path to reference data directory (overrides relative paths in config)')
    parser.add_argument('--skip-sources', default='',
                        help='Comma-separated sources to skip, e.g. "oncokb,tcga"')
    args = parser.parse_args()
    skip_sources = {s.strip().lower() for s in args.skip_sources.split(',') if s.strip()}

    section('Configuration')
    if not file_exists(args.config):
        fail('config.yaml exists', args.config)
        sys.exit(1)
    ok('config.yaml exists')

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    # Apply data_dir prefix to all relative paths if provided
    data_dir = args.data_dir.rstrip('/') if args.data_dir else ''
    if data_dir:
        path_keys = {'tsv', 'vcf', 'maf', 'clinical_sample', 'somatic_tsv',
                     'germline_tsv', 'variants_file', 'dir', 'chain'}

        def _repath(v):
            if not v or os.path.isabs(v):
                return v
            parts = v.replace('\\', '/').split('/', 1)
            remainder = parts[1] if len(parts) > 1 else parts[0]
            return os.path.join(data_dir, remainder)

        for source in cfg.get('data_sources', {}).values():
            if not isinstance(source, dict):
                continue
            for k, v in source.items():
                if k in path_keys and isinstance(v, str):
                    source[k] = _repath(v)
        ok(f'Data directory: {data_dir}')

    settings_file = cfg.get('settings_file', 'settings.yaml')
    if file_exists(settings_file):
        ok(f'settings file: {settings_file}')
    else:
        fail(f'settings file missing: {settings_file}')

    ds = cfg.get('data_sources', {})

    section('Required tools')
    for tool in ['bgzip', 'tabix', 'bcftools', 'python3']:
        r = subprocess.run(['which', tool], capture_output=True)
        if r.returncode == 0:
            ok(f'{tool}')
        else:
            fail(f'{tool} not found on PATH')

    section('Python dependencies')
    pkgs = {'pandas': 'pandas', 'yaml': 'pyyaml', 'requests': 'requests',
            'pysam': 'pysam', 'polars': 'polars'}
    for pkg, install in pkgs.items():
        try:
            __import__(pkg)
            ok(pkg)
        except ImportError:
            if pkg == 'polars':
                warn('polars not installed — aggregation will use slow pandas fallback',
                     detail='pip install polars --break-system-packages')
            else:
                fail(f'{pkg} not installed',
                     detail=f'pip install {install} --break-system-packages')

    section('COSMIC')
    if not ds.get('cosmic', {}).get('enabled', False):
        skip('COSMIC', 'disabled')
    else:
        tsv = ds['cosmic'].get('tsv', '')
        vcf = ds['cosmic'].get('vcf', '')
        if file_exists(tsv):
            ok(f'TSV: {tsv}')
            cols = detect_columns(tsv)
            info(f'{len(cols)} columns detected', ', '.join(cols[:6]) + ' ...')
            if 'PRIMARY_SITE' not in cols and 'COSMIC_PHENOTYPE_ID' in cols:
                info('No PRIMARY_SITE column — cancer types will use numeric phenotype IDs',
                     'Download Cosmic_Classification file for readable names')
            chrom = detect_chrom_style(tsv)
            info(f'Chromosome style: {chrom["style"]}', str(chrom["sample"]))
        else:
            fail(f'TSV missing: {tsv}')
        if file_exists(vcf):
            ok(f'VCF: {vcf}')
            ok('VCF BGZF') if is_bgzf(vcf) else warn('VCF not BGZF')
            if not has_index(vcf, ['.tbi', '.csi']):
                info('VCF has no tabix index', 'Not required — pipeline uses TSV')
        else:
            warn(f'VCF missing: {vcf}', 'Not required — pipeline uses TSV')

    section('ClinVar')
    if not ds.get('clinvar', {}).get('enabled', False):
        skip('ClinVar', 'disabled')
    else:
        vcf = ds['clinvar'].get('vcf', '')
        if not file_exists(vcf):
            fail(f'VCF missing: {vcf}')
        else:
            ok(f'VCF: {vcf}')
            if is_bgzf(vcf):
                ok('BGZF compressed')
            else:
                fail('Not BGZF', 'bgzip -f data/clinvar/clinvar.vcf')
            if has_index(vcf, ['.csi']):
                ok('CSI index present')
            elif has_index(vcf, ['.tbi']):
                warn('TBI index only — CSI required for GRCh38',
                     'tabix -p vcf -C data/clinvar/clinvar.vcf.gz')
            else:
                fail('No index', 'tabix -p vcf -C data/clinvar/clinvar.vcf.gz')
            chrom = detect_chrom_style(vcf)
            info(f'Chromosome style: {chrom["style"]}', str(chrom["sample"]))

    section('GENIE')
    if not ds.get('genie', {}).get('enabled', False):
        skip('GENIE', 'disabled')
    else:
        maf = ds['genie'].get('maf', '')
        clinical = ds['genie'].get('clinical_sample', '')
        if file_exists(maf):
            ok(f'MAF: {maf}')
            cols = detect_columns(maf)
            info(f'{len(cols)} columns detected', ', '.join(cols[:6]) + ' ...')
            chrom = detect_chrom_style(maf)
            info(f'Chromosome style: {chrom["style"]}', str(chrom["sample"]))
            # Check NCBI_Build to confirm liftover has been run
            build_val = ''
            try:
                import gzip as _gz
                _opener = _gz.open if maf.endswith('.gz') else open
                with _opener(maf, 'rt') as _fh:
                    for _line in _fh:
                        if _line.startswith('#'):
                            continue
                        _hdrs = _line.rstrip('\n').split('\t')
                        _col_i = {c: i for i, c in enumerate(_hdrs)}
                        _nb_i = _col_i.get('NCBI_Build')
                        if _nb_i is None:
                            break
                        for _dline in _fh:
                            _parts = _dline.rstrip('\n').split('\t')
                            if _nb_i < len(_parts):
                                build_val = _parts[_nb_i].strip()
                            break
                        break
            except Exception:
                pass
            if build_val.upper() in ('GRCH38', 'HG38', '38'):
                ok(f'NCBI_Build: {build_val} (GRCh38 confirmed)')
            elif build_val:
                fail(
                    f'GENIE NCBI_Build: {build_val} — MAF appears to be GRCh37, not GRCh38',
                    'Run: python3 tools/db_fix.py --genie-maf data/genie/data_mutations_extended.txt\n'
                    '        Then update config.yaml genie.maf to point to the lifted file.'
                )
            else:
                warn('Could not detect NCBI_Build — verify MAF is GRCh38')
        else:
            fail(
                f'MAF missing: {maf}',
                'Download GENIE MAF from Synapse: https://www.synapse.org/#!Synapse:syn7222066\n'
                '        Then run: python3 tools/db_fix.py --genie-maf <path_to_downloaded_maf>'
            )
        if file_exists(clinical):
            ok(f'Clinical sample file: {clinical}')
        else:
            warn(f'Clinical sample file missing: {clinical}',
                 'Cancer type labels will be unavailable for GENIE variants')

    section('TCGA')
    if 'tcga' in skip_sources:
        skip('TCGA', '--skip-sources tcga specified')
    elif not ds.get('tcga', {}).get('enabled', False):
        skip('TCGA', 'disabled')
    else:
        maf = ds['tcga'].get('maf', '')
        if file_exists(maf):
            ok(f'MAF: {maf}')
            cols = detect_columns(maf)
            info(f'{len(cols)} columns detected', ', '.join(cols[:6]) + ' ...')
            chrom = detect_chrom_style(maf)
            info(f'Chromosome style: {chrom["style"]}', str(chrom["sample"]))
            # Check NCBI_Build to confirm liftover has been run
            build_val = ''
            try:
                import gzip as _gz
                opener = _gz.open if maf.endswith('.gz') else open
                with opener(maf, 'rt') as _fh:
                    for _line in _fh:
                        if _line.startswith('#'):
                            continue
                        _hdrs = _line.rstrip('\n').split('\t')
                        _col_i = {c: i for i, c in enumerate(_hdrs)}
                        _nb_i = _col_i.get('NCBI_Build')
                        if _nb_i is None:
                            break
                        # Read first data row
                        for _dline in _fh:
                            _parts = _dline.rstrip('\n').split('\t')
                            if _nb_i < len(_parts):
                                build_val = _parts[_nb_i].strip()
                            break
                        break
            except Exception:
                pass
            if build_val.upper() in ('GRCH38', 'HG38', '38'):
                ok(f'NCBI_Build: {build_val} (GRCh38 confirmed)')
            elif build_val:
                fail(
                    f'NCBI_Build: {build_val} — MAF appears to be GRCh37, not GRCh38',
                    'Run: python3 tools/db_fix.py --tcga-maf data/TCGA/mc3.v0.2.8.PUBLIC.maf.gz\n'
                    '        Then update config.yaml tcga.maf to point to the lifted file.'
                )
            else:
                warn('Could not detect NCBI_Build — verify MAF is GRCh38')
        else:
            fail(
                f'MAF missing: {maf}',
                'Download mc3.v0.2.8.PUBLIC.maf.gz from:\n'
                '        https://gdc.cancer.gov/about-data/publications/pancanatlas\n'
                '        Then run: python3 tools/db_fix.py --tcga-maf data/TCGA/mc3.v0.2.8.PUBLIC.maf.gz'
            )

    section('TP53 database')
    if not ds.get('tp53', {}).get('enabled', False):
        skip('TP53', 'disabled')
    else:
        somatic = ds['tp53'].get('somatic_tsv', '')
        if file_exists(somatic):
            ok(f'Somatic file: {somatic}')
            cols = detect_columns(somatic)
            info(f'{len(cols)} columns detected', ', '.join(cols[:6]) + ' ...')
            chrom = detect_chrom_style(somatic)
            info(f'Chromosome style: {chrom["style"]}', str(chrom["sample"]))
            info('Single-gene database (chr17 only) — chromosome style detection not applicable')
        else:
            fail(f'Somatic file missing: {somatic}')
        germline = ds['tp53'].get('germline_tsv', '')
        if ds['tp53'].get('include_germline', False) and not file_exists(germline):
            fail(f'Germline file missing: {germline}',
                 'Set include_germline: false in settings.yaml to skip')

    section('OncoKB')
    if 'oncokb' in skip_sources:
        skip('OncoKB', '--skip-sources oncokb specified')
    elif not ds.get('oncokb', {}).get('enabled', False):
        skip('OncoKB', 'disabled')
    else:
        variants_file = ds['oncokb'].get('variants_file', '')
        if file_exists(variants_file):
            ok(f'Variants file: {variants_file}')
            cols = detect_columns(variants_file)
            info(f'{len(cols)} columns detected', ', '.join(cols[:6]))
        else:
            info('Local variants file not found — will use API', variants_file)
        if file_exists(settings_file):
            with open(settings_file) as fh:
                settings = yaml.safe_load(fh)
            token = settings.get('oncokb', {}).get('api_token', '')
            if token and token not in ('', 'YOUR_TOKEN_HERE'):
                ok('API token set')
            elif file_exists(variants_file):
                ok('API token not set — using local variants file (no token required)')
            else:
                fail(
                    'OncoKB API token is missing or not set',
                    'No local variants file found and no API token provided.\n'
                    '        The pipeline will fire hundreds of API requests, all returning\n'
                    '        HTTP 401, and oncokb_oncogenicity will be empty in the output.\n'
                    '        Fix option 1: set oncokb.api_token in settings.yaml\n'
                    '                      Tokens are free: https://www.oncokb.org/account/register\n'
                    '                      Tokens expire every 6 months — renew at the same URL.\n'
                    '        Fix option 2: download allAnnotatedVariants.txt from OncoKB and\n'
                    '                      set oncokb.variants_file in config.yaml.'
                )

    section('Chromosome naming')
    print('  All sources are normalised to plain style (1, 2, X, Y, MT)')
    print('  at aggregation regardless of input format.')

    section('Output paths')
    for d in [cfg.get('output', {}).get('dir', 'output'),
              cfg.get('intermediate_dir', 'intermediate')]:
        os.makedirs(d, exist_ok=True)
        test = os.path.join(d, '.write_test')
        try:
            open(test, 'w').close()
            os.remove(test)
            ok(f'Writable: {d}')
        except Exception as e:
            fail(f'Not writable: {d}', str(e))

    print(f'\n{"="*55}')
    if failures == 0 and warnings == 0:
        print(f'  {PASS}  All checks passed. Safe to run pipeline.')
    elif failures == 0:
        print(f'  {WARN}  {warnings} warning(s). Pipeline can run — review above.')
    else:
        print(f'  {FAIL}  {failures} failure(s). Fix before running pipeline.')
    print('='*55 + '\n')

    sys.exit(0 if failures == 0 else 1)


if __name__ == '__main__':
    main()
