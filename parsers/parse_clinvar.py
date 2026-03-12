#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: ClinVar VCF (GRCh38)
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_clinvar.py
Parse ClinVar GRCh38 VCF for somatic/oncogenic pathogenic variants.

Required file:
  clinvar.vcf.gz
  Download: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
"""

import gzip
import os
import re

import pandas as pd

from parsers.common import (
    STANDARD_COLS,
    clean_allele,
    empty_standard_df,
    is_valid_allele,
    map_consequence,
    normalise_chrom,
    setup_logger,
)

log = setup_logger('ClinVar')

_CLINSIG_INCLUDE = {
    'Pathogenic',
    'Likely_pathogenic',
    'Pathogenic/Likely_pathogenic',
}

_CLNSIG_RE   = re.compile(r'CLNSIG=([^;]+)')
_CLNDN_RE    = re.compile(r'CLNDN=([^;]+)')
_GENEINFO_RE = re.compile(r'GENEINFO=([^;|:]+)')
_MC_RE       = re.compile(r'MC=([^;]+)')
_ORIGIN_RE   = re.compile(r'CLNORIGIN=([^;]+)')


def parse_clinvar(vcf_path: str,
                  include_clinsig: list | None = None,
                  somatic_only: bool = True) -> pd.DataFrame:
    if include_clinsig is None:
        include_clinsig = list(_CLINSIG_INCLUDE)

    if not os.path.exists(vcf_path):
        log.warning('ClinVar VCF not found: %s  — skipping', vcf_path)
        return empty_standard_df()

    log.info('Parsing ClinVar VCF: %s', vcf_path)
    rows = []
    opener = gzip.open if vcf_path.endswith('.gz') else open

    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue

            chrom, pos_str, _, ref, alt_str, _, _, info = parts[:8]
            chrom = normalise_chrom(chrom)
            ref   = clean_allele(ref)

            # Skip multi-allelic records (uncommon in ClinVar VCF)
            alts = alt_str.split(',')

            clnsig_m = _CLNSIG_RE.search(info)
            if not clnsig_m:
                continue
            clnsig = clnsig_m.group(1).replace('_', ' ').strip()
            # ClinVar may pipe-delimit multiple assertions
            clnsig_vals = set(v.strip() for v in re.split(r'[|/,]', clnsig_m.group(1)))

            if not clnsig_vals.intersection(include_clinsig):
                continue

            # Optionally restrict to somatic/oncogenic assertions
            if somatic_only:
                origin_m = _ORIGIN_RE.search(info)
                if origin_m:
                    # CLNORIGIN bitmask: 2 = somatic, 4 = inherited, 8 = maternal,
                    # 16 = paternal, 32 = de novo, 64 = biparental, 128 = uniparental,
                    # 256 = not-tested, 512 = tested-inconclusive, 1073741824 = other
                    origin_val = int(origin_m.group(1))
                    if not (origin_val & 2):     # bit 2 = somatic
                        continue

            gene_m    = _GENEINFO_RE.search(info)
            gene      = gene_m.group(1).strip() if gene_m else ''
            clndn_m   = _CLNDN_RE.search(info)
            cancer_type = clndn_m.group(1).replace('_', ' ').lower() \
                          if clndn_m else 'unspecified'

            mc_m      = _MC_RE.search(info)
            consequence = 'unknown'
            if mc_m:
                # MC field: e.g. SO:0001583|missense_variant
                mc_parts = mc_m.group(1).split('|')
                if len(mc_parts) > 1:
                    consequence = map_consequence(mc_parts[1])

            for alt in alts:
                alt = clean_allele(alt)
                if not is_valid_allele(ref) or not is_valid_allele(alt):
                    continue
                rows.append({
                    'chrom':        chrom,
                    'pos':          int(pos_str),
                    'ref':          ref,
                    'alt':          alt,
                    'gene':         gene,
                    'hgvsc':        '',
                    'hgvsp':        '',
                    'consequence':  consequence,
                    'cancer_type':  cancer_type,
                    'n_samples':    1,   # ClinVar doesn't provide sample counts
                    'source':       'ClinVar',
                })

    if not rows:
        log.warning('ClinVar: no rows produced')
        return empty_standard_df()

    df = pd.DataFrame(rows, columns=STANDARD_COLS)
    log.info('ClinVar: %d rows', len(df))
    return df
