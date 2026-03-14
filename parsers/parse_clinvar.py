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

Filtering strategy:
  1. CLNSIG: Pathogenic / Likely_pathogenic / Pathogenic/Likely_pathogenic only.
  2. CLNORIGIN bitmask:
       - If present and somatic bit (2) is set: include.
       - If present and somatic bit is NOT set: exclude.
       - If absent: fall through to cancer disease name filter.
  3. CLNDN cancer filter: when CLNORIGIN is absent, include only if the disease
     name contains a recognised cancer/oncology term.
  4. n_samples is set to 0 for all ClinVar entries; ClinVar does not provide
     cohort sample counts and tiering relies on evidence level alone.
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

# Cancer and oncology terms matched against CLNDN (lowercase).
# Used when CLNORIGIN is absent to restrict to cancer-relevant entries.
_CANCER_TERMS = {
    'cancer', 'carcinoma', 'adenocarcinoma', 'sarcoma', 'lymphoma',
    'leukemia', 'leukaemia', 'melanoma', 'glioma', 'glioblastoma',
    'mesothelioma', 'myeloma', 'neuroblastoma', 'hepatocellular',
    'cholangiocarcinoma', 'tumor', 'tumour', 'neoplasm', 'neoplasia',
    'malignant', 'malignancy', 'oncogenic', 'metastatic', 'metastasis',
    'blastoma', 'seminoma', 'teratoma', 'thymoma', 'meningioma',
    'myelodysplastic', 'myeloproliferative', 'polycythemia', 'polycythaemia',
    'thrombocythemia', 'thrombocythaemia', 'myelofibrosis',
    'paraganglioma', 'pheochromocytoma', 'phaeochromocytoma',
    'ependymoma', 'astrocytoma', 'medulloblastoma', 'retinoblastoma',
    'wilms', 'ewing', 'osteosarcoma', 'chondrosarcoma', 'rhabdomyosarcoma',
    'liposarcoma', 'leiomyosarcoma', 'angiosarcoma', 'fibrosarcoma',
    'pancreatic', 'colorectal', 'hereditary breast', 'hereditary ovarian',
    'lynch', 'cowden', 'li-fraumeni', 'von hippel', 'birt-hogg',
    'familial adenomatous', 'multiple endocrine', 'peutz',
    'juvenile polyposis', 'hereditary diffuse gastric',
}

_CLNSIG_RE   = re.compile(r'CLNSIG=([^;]+)')
_CLNDN_RE    = re.compile(r'CLNDN=([^;]+)')
_GENEINFO_RE = re.compile(r'GENEINFO=([^;|:]+)')
_MC_RE       = re.compile(r'MC=([^;]+)')
_ORIGIN_RE   = re.compile(r'CLNORIGIN=([^;]+)')


def _is_cancer_disease(clndn: str) -> bool:
    """Return True if the CLNDN field contains a recognised cancer term."""
    clndn_lower = clndn.lower().replace('_', ' ')
    return any(term in clndn_lower for term in _CANCER_TERMS)


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
    n_skipped_origin  = 0
    n_skipped_disease = 0
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

            alts = alt_str.split(',')

            # --- Clinical significance filter ---
            clnsig_m = _CLNSIG_RE.search(info)
            if not clnsig_m:
                continue
            clnsig_vals = set(
                v.strip() for v in re.split(r'[|/,]', clnsig_m.group(1))
            )
            if not clnsig_vals.intersection(include_clinsig):
                continue

            # --- Origin / cancer relevance filter ---
            if somatic_only:
                origin_m = _ORIGIN_RE.search(info)
                if origin_m:
                    # CLNORIGIN bitmask: 2 = somatic
                    origin_val = int(origin_m.group(1))
                    if not (origin_val & 2):
                        n_skipped_origin += 1
                        continue
                    # Somatic bit set — include regardless of disease name.
                else:
                    # CLNORIGIN absent — gate on cancer disease name.
                    clndn_m = _CLNDN_RE.search(info)
                    clndn   = clndn_m.group(1) if clndn_m else ''
                    if not _is_cancer_disease(clndn):
                        n_skipped_disease += 1
                        continue

            gene_m  = _GENEINFO_RE.search(info)
            gene    = gene_m.group(1).strip() if gene_m else ''

            clndn_m     = _CLNDN_RE.search(info)
            cancer_type = clndn_m.group(1).split('|')[0].replace('_', ' ').lower().strip() \
                          if clndn_m else 'unspecified'

            mc_m        = _MC_RE.search(info)
            consequence = 'unknown'
            if mc_m:
                mc_parts = mc_m.group(1).split('|')
                if len(mc_parts) > 1:
                    consequence = map_consequence(mc_parts[1])

            for alt in alts:
                alt = clean_allele(alt)
                if not is_valid_allele(ref) or not is_valid_allele(alt):
                    continue
                rows.append({
                    'chrom':       chrom,
                    'pos':         int(pos_str),
                    'ref':         ref,
                    'alt':         alt,
                    'gene':        gene,
                    'hgvsc':       '',
                    'hgvsp':       '',
                    'consequence': consequence,
                    'cancer_type': cancer_type,
                    'n_samples':   0,
                    'source':      'ClinVar',
                    'clinvar_clinical_significance': clnsig_m.group(1).strip(),
                })

    log.info(
        'ClinVar: %d rows kept | %d excluded by CLNORIGIN | %d excluded by disease name',
        len(rows), n_skipped_origin, n_skipped_disease
    )

    if not rows:
        log.warning('ClinVar: no rows produced — check filters and input file')
        return empty_standard_df()

    df = pd.DataFrame(rows)
    for col in STANDARD_COLS:
        if col not in df.columns:
            df[col] = ''
    df = df[STANDARD_COLS + ['clinvar_clinical_significance']]
    return df
