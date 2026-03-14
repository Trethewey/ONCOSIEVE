#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: TCGA PanCancer Atlas mc3 MAF
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_tcga.py
Parse the TCGA PanCancer Atlas mc3 somatic mutation MAF.

Required file:
  mc3.v0.2.8.PUBLIC.GRCh38.maf.gz
  Source: https://gdc.cancer.gov/about-data/publications/pancanatlas
  File:   mc3.v0.2.8.PUBLIC.maf.gz  (md5: 639ad8f8386e98dacc22e439188aa8fa)
  The GRCh37 MAF must be lifted to GRCh38 before use.
  Run: python3 tools/db_fix.py --config config.yaml

The mc3 MAF is the TCGA PanCancer Atlas harmonised somatic mutation call set,
produced by merging calls from 6 variant callers across 10,295 tumour-normal
pairs from 33 cancer types.

Key columns used:
  Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2
  Hugo_Symbol, HGVSc, HGVSp_Short, Transcript_ID
  Consequence      — VEP consequence term (used directly, no remapping needed)
  Tumor_Sample_Barcode — encodes TCGA project code for cancer type
  FILTER           — only PASS variants are included

Cancer type is derived from the TCGA project code in Tumor_Sample_Barcode
(e.g. TCGA-BRCA-... -> 'breast invasive carcinoma') using the TCGA project
code lookup table defined in this module.
"""

import gzip
import os

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

log = setup_logger('TCGA')

# Column names in the mc3 MAF
_COL_CHROM    = 'Chromosome'
_COL_POS      = 'Start_Position'
_COL_REF      = 'Reference_Allele'
_COL_ALT      = 'Tumor_Seq_Allele2'
_COL_GENE     = 'Hugo_Symbol'
_COL_HGVSC    = 'HGVSc'
_COL_HGVSP    = 'HGVSp_Short'
_COL_CSQ      = 'Consequence'       # VEP term — preferred over Variant_Classification
_COL_VARCLASS = 'Variant_Classification'  # fallback if Consequence absent
_COL_SAMPLE   = 'Tumor_Sample_Barcode'
_COL_FILTER   = 'FILTER'

# TCGA project code -> cancer type label
# Source: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
_TCGA_PROJECT_MAP: dict[str, str] = {
    'ACC':  'adrenocortical carcinoma',
    'BLCA': 'bladder urothelial carcinoma',
    'BRCA': 'breast invasive carcinoma',
    'CESC': 'cervical squamous cell carcinoma',
    'CHOL': 'cholangiocarcinoma',
    'COAD': 'colon adenocarcinoma',
    'DLBC': 'diffuse large b-cell lymphoma',
    'ESCA': 'esophageal carcinoma',
    'GBM':  'glioblastoma multiforme',
    'HNSC': 'head and neck squamous cell carcinoma',
    'KICH': 'kidney chromophobe',
    'KIRC': 'kidney renal clear cell carcinoma',
    'KIRP': 'kidney renal papillary cell carcinoma',
    'LAML': 'acute myeloid leukemia',
    'LGG':  'brain lower grade glioma',
    'LIHC': 'liver hepatocellular carcinoma',
    'LUAD': 'lung adenocarcinoma',
    'LUSC': 'lung squamous cell carcinoma',
    'MESO': 'mesothelioma',
    'OV':   'ovarian serous cystadenocarcinoma',
    'PAAD': 'pancreatic adenocarcinoma',
    'PCPG': 'pheochromocytoma and paraganglioma',
    'PRAD': 'prostate adenocarcinoma',
    'READ': 'rectum adenocarcinoma',
    'SARC': 'sarcoma',
    'SKCM': 'skin cutaneous melanoma',
    'STAD': 'stomach adenocarcinoma',
    'TGCT': 'testicular germ cell tumors',
    'THCA': 'thyroid carcinoma',
    'THYM': 'thymoma',
    'UCEC': 'uterine corpus endometrial carcinoma',
    'UCS':  'uterine carcinosarcoma',
    'UVM':  'uveal melanoma',
}


def _cancer_type_from_barcode(barcode: str) -> str:
    """
    Extract cancer type from TCGA tumour sample barcode.
    Format: TCGA-{project}-{participant}-{sample}-...
    e.g. TCGA-BRCA-A1B2-01A -> 'breast invasive carcinoma'
    """
    parts = str(barcode).split('-')
    if len(parts) >= 2:
        code = parts[1].upper()
        return _TCGA_PROJECT_MAP.get(code, code.lower())
    return 'unspecified'


def parse_tcga(maf_path: str) -> pd.DataFrame:
    """
    Parse TCGA mc3 PanCancer Atlas MAF (GRCh38 lifted).

    Parameters
    ----------
    maf_path : path to mc3.v0.2.8.PUBLIC.GRCh38.maf.gz

    Returns
    -------
    DataFrame with STANDARD_COLS schema.
    One row per mutation per sample. Aggregation to unique samples per variant
    is performed at the merge step.
    """
    if not os.path.exists(maf_path):
        log.warning('TCGA MAF not found: %s  — skipping', maf_path)
        return empty_standard_df()

    log.info('Parsing TCGA mc3 MAF: %s', maf_path)

    opener = gzip.open if maf_path.endswith('.gz') else open

    # Read header to get column indices
    col_idx: dict[str, int] = {}
    try:
        with opener(maf_path, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                col_idx = {c: i for i, c in enumerate(cols)}
                break
    except Exception as e:
        log.error('Failed to read TCGA MAF header: %s', e)
        return empty_standard_df()

    required = [_COL_CHROM, _COL_POS, _COL_REF, _COL_ALT,
                _COL_GENE, _COL_SAMPLE]
    missing = [c for c in required if c not in col_idx]
    if missing:
        log.error('TCGA MAF missing required columns: %s', missing)
        return empty_standard_df()

    chrom_i   = col_idx[_COL_CHROM]
    pos_i     = col_idx[_COL_POS]
    ref_i     = col_idx[_COL_REF]
    alt_i     = col_idx[_COL_ALT]
    gene_i    = col_idx[_COL_GENE]
    hgvsc_i   = col_idx.get(_COL_HGVSC)
    hgvsp_i   = col_idx.get(_COL_HGVSP)
    csq_i     = col_idx.get(_COL_CSQ)
    varclass_i = col_idx.get(_COL_VARCLASS)
    sample_i  = col_idx[_COL_SAMPLE]
    filter_i  = col_idx.get(_COL_FILTER)

    rows = []
    n_filtered = 0

    try:
        with opener(maf_path, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')

                # Skip header row
                if parts[chrom_i] == _COL_CHROM:
                    continue

                # FILTER column: keep PASS only
                if filter_i is not None and filter_i < len(parts):
                    if parts[filter_i].strip().upper() not in ('PASS', '.', ''):
                        n_filtered += 1
                        continue

                try:
                    chrom = normalise_chrom(parts[chrom_i])
                    pos   = int(float(parts[pos_i]))
                    ref   = clean_allele(parts[ref_i])
                    alt   = clean_allele(parts[alt_i])

                    if not is_valid_allele(ref) or not is_valid_allele(alt):
                        continue
                    if alt in ('', '-', '.'):
                        continue

                    gene   = str(parts[gene_i]).strip()
                    hgvsc  = str(parts[hgvsc_i]).strip()  if hgvsc_i  is not None and hgvsc_i  < len(parts) else ''
                    hgvsp  = str(parts[hgvsp_i]).strip()  if hgvsp_i  is not None and hgvsp_i  < len(parts) else ''
                    sample = str(parts[sample_i]).strip()

                    # Consequence: prefer VEP Consequence column, fall back to
                    # Variant_Classification mapped through map_consequence()
                    if csq_i is not None and csq_i < len(parts):
                        raw_csq = parts[csq_i].strip()
                        # VEP may return pipe-separated terms; take first
                        consequence = map_consequence(raw_csq.split('&')[0].split('|')[0])
                    elif varclass_i is not None and varclass_i < len(parts):
                        consequence = map_consequence(parts[varclass_i].strip())
                    else:
                        consequence = 'unknown'

                    cancer_type = _cancer_type_from_barcode(sample)

                    rows.append({
                        'chrom':       chrom,
                        'pos':         pos,
                        'ref':         ref,
                        'alt':         alt,
                        'gene':        gene,
                        'hgvsc':       hgvsc,
                        'hgvsp':       hgvsp,
                        'consequence': consequence,
                        'cancer_type': cancer_type,
                        'n_samples':   1,
                        'source':      'TCGA',
                    })

                except (ValueError, TypeError, IndexError):
                    continue

    except Exception as e:
        log.error('Failed to parse TCGA MAF: %s', e)
        return empty_standard_df()

    log.info(
        'TCGA: %d rows parsed, %d excluded by FILTER',
        len(rows), n_filtered
    )

    if not rows:
        log.warning('TCGA: no rows produced — check MAF path and FILTER column')
        return empty_standard_df()

    df = pd.DataFrame(rows, columns=STANDARD_COLS)
    log.info('TCGA: %d rows', len(df))
    return df
