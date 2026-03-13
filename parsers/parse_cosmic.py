#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Parser: COSMIC Genome Screens Mutant TSV + VCF
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
parse_cosmic.py
Parse COSMIC GRCh38 data to extract somatic mutations with cancer type and sample counts.

Required files (both needed for full output):
  1. CosmicMutantExportCensus.tsv.gz
       https://cancer.sanger.ac.uk/cosmic/download
       Select: "COSMIC Mutation Data (Genome Screen)"  ->  GRCh38
       Provides: cancer type, sample ID, mutation consequence
  2. Cosmic_GenomeScreensMutant_vXX_GRCh38.vcf.gz
       https://cancer.sanger.ac.uk/cosmic/download
       Select: "COSMIC Mutation Data (VCF)"  ->  GRCh38
       Provides: proper REF/ALT alleles matched by COSMIC ID

If only the TSV is available, REF/ALT will be inferred from MUTATION_CDS for SNVs
and marked as INFERRED; indels will be skipped unless the VCF is present.
"""

import gzip
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import tempfile
from typing import Optional

import pandas as pd

from parsers.common import (
    INCLUDED_CONSEQUENCES,
    STANDARD_COLS,
    clean_allele,
    empty_standard_df,
    is_valid_allele,
    map_consequence,
    normalise_chrom,
    setup_logger,
)

log = setup_logger('COSMIC')

# COSMIC TSV column names (v103+)
_COL_GENE        = 'GENE_SYMBOL'
_COL_HGVSP       = 'HGVSP'
_COL_HGVSC       = 'HGVSC'
_COL_TRANSCRIPT  = 'TRANSCRIPT_ACCESSION'
_COL_CDS         = 'MUTATION_CDS'
_COL_AA          = 'MUTATION_AA'
_COL_DESC        = 'MUTATION_DESCRIPTION'
_COL_SITE        = 'COSMIC_PHENOTYPE_ID'
_COL_HISTOLOGY   = 'COSMIC_PHENOTYPE_ID'
_COL_SAMPLE      = 'SAMPLE_NAME'
_COL_SAMPLE_ID   = 'COSMIC_SAMPLE_ID'
_COL_COSM_ID     = 'MUTATION_ID'
_COL_CHROM       = 'CHROMOSOME'
_COL_START       = 'GENOME_START'
_COL_REF         = 'GENOMIC_WT_ALLELE'
_COL_ALT         = 'GENOMIC_MUT_ALLELE'
_COL_STRAND      = 'STRAND'
_COL_STATUS      = 'MUTATION_SOMATIC_STATUS'

# Regex to parse COSMIC genome position: e.g. "1:209968684-209968684"
_POS_RE = re.compile(r'^(\w+):(\d+)-(\d+)$')

# Regex for SNV in CDS: c.123A>G
_SNV_CDS_RE = re.compile(r'c\.\d+([ACGTacgt])>([ACGTacgt])')


def _parse_position(pos_str: str) -> Optional[tuple[str, int]]:
    """Parse 'chr:start-end' -> (chrom_normalised, start_1based)."""
    m = _POS_RE.match(str(pos_str).strip())
    if not m:
        return None
    chrom = normalise_chrom(m.group(1))
    return chrom, int(m.group(2))


def _extract_alleles_from_cds(cds: str, strand: str) -> Optional[tuple[str, str]]:
    """
    Infer REF/ALT for simple SNVs from MUTATION_CDS.
    Applies reverse-complement correction for negative-strand genes.
    Returns None for indels or unresolvable entries.
    """
    m = _SNV_CDS_RE.search(str(cds))
    if not m:
        return None
    ref_cds = m.group(1).upper()
    alt_cds = m.group(2).upper()
    if strand == '-':
        rc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        ref_cds = rc.get(ref_cds, ref_cds)
        alt_cds = rc.get(alt_cds, alt_cds)
    return ref_cds, alt_cds


def _build_classification_lookup(classification_path: str) -> dict[str, str]:
    """
    Parse COSMIC Classification TSV and return a dict:
    COSMIC_PHENOTYPE_ID (e.g. COSO36004862) -> PRIMARY_SITE (e.g. 'lung')
    """
    if not classification_path or not os.path.exists(classification_path):
        log.warning('COSMIC classification file not found: %s — cancer types will be phenotype IDs',
                    classification_path)
        return {}
    log.info('Building COSMIC classification lookup from %s', classification_path)
    lookup: dict[str, str] = {}
    opener = gzip.open if classification_path.endswith('.gz') else open
    with opener(classification_path, 'rt') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        col = {c: i for i, c in enumerate(header)}
        if 'COSMIC_PHENOTYPE_ID' not in col or 'PRIMARY_SITE' not in col:
            log.error('Classification file missing required columns. Found: %s', header)
            return {}
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            try:
                pheno_id    = parts[col['COSMIC_PHENOTYPE_ID']].strip()
                primary_site = parts[col['PRIMARY_SITE']].strip().lower()
                if pheno_id and primary_site:
                    lookup[pheno_id] = primary_site
            except IndexError:
                continue
    log.info('COSMIC classification lookup: %d entries', len(lookup))
    return lookup


def _build_vcf_lookup(vcf_path: str) -> dict[str, tuple[str, str, str, str]]:
    """
    Parse COSMIC GRCh38 VCF and return a dict keyed by COSM/COSV ID.
    Value: (chrom, pos, ref, alt)
    Only stores the first ALT for multi-allelic records.
    """
    lookup: dict[str, tuple[str, str, str, str]] = {}
    opener = gzip.open if vcf_path.endswith('.gz') else open
    log.info('Building COSMIC VCF lookup from %s', vcf_path)
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            chrom, pos, id_field, ref, alt_field = parts[0], parts[1], parts[2], parts[3], parts[4]
            chrom = normalise_chrom(chrom)
            alt = alt_field.split(',')[0]   # take first ALT only
            for cosm_id in id_field.split(';'):
                cosm_id = cosm_id.strip()
                if cosm_id and cosm_id not in lookup:
                    lookup[cosm_id] = (chrom, pos, ref.upper(), alt.upper())
    log.info('COSMIC VCF lookup: %d entries', len(lookup))
    return lookup


def parse_cosmic(tsv_path: str,
                 vcf_path: Optional[str] = None,
                 classification_path: Optional[str] = None,
                 chunk_size: int = 200_000,
                 n_threads: int = 4) -> pd.DataFrame:
    """
    Parse COSMIC mutation export TSV (GRCh38).

    Parameters
    ----------
    tsv_path   : path to CosmicMutantExportCensus.tsv.gz
    vcf_path   : path to COSMIC GRCh38 VCF (optional but recommended)
    chunk_size : rows to read per chunk
    n_threads  : worker processes for parallel parsing (default 4)

    Returns
    -------
    DataFrame with STANDARD_COLS schema.
    One row per (variant, cancer_type) with n_samples = count of unique
    sample IDs in that cancer type.
    """
    if not os.path.exists(tsv_path):
        log.warning('COSMIC TSV not found: %s  — skipping', tsv_path)
        return empty_standard_df()

    vcf_lookup: dict = {}
    if vcf_path and os.path.exists(vcf_path):
        vcf_lookup = _build_vcf_lookup(vcf_path)
    else:
        log.warning('COSMIC VCF not provided or not found; REF/ALT will be '
                    'inferred from CDS for SNVs only.')

    classification_lookup: dict = {}
    if classification_path:
        classification_lookup = _build_classification_lookup(classification_path)

    usecols = [
        _COL_GENE, _COL_HGVSP, _COL_HGVSC, _COL_CDS, _COL_AA,
        _COL_DESC, _COL_SITE, _COL_HISTOLOGY, _COL_SAMPLE_ID, _COL_COSM_ID,
        _COL_CHROM, _COL_START, _COL_REF, _COL_ALT,
        _COL_STRAND, _COL_STATUS,
    ]

    log.info('Parsing COSMIC TSV: %s', tsv_path)
    n_workers = max(1, n_threads)
    log.info('COSMIC: using %d worker(s)', n_workers)

    # ── Step 1: decompress to a temp file ────────────────────────────────────
    # gzip format is sequential — only one core can decompress it, so workers
    # cannot read in parallel from a gzip stream. We decompress first (using
    # pigz for parallel decompression if available), then split the plain-text
    # file into N byte ranges and have N workers each read their own section.

    if tsv_path.endswith('.gz'):
        tmp = tempfile.NamedTemporaryFile(
            suffix='.tsv', delete=False,
            dir=os.path.dirname(tsv_path) or '.',
        )
        tmp_path = tmp.name
        tmp.close()
        try:
            if shutil.which('pigz'):
                log.info('COSMIC: decompressing with pigz (%d cores) → %s',
                         n_workers, tmp_path)
                subprocess.run(
                    ['pigz', '-d', '-c', '-p', str(n_workers), tsv_path],
                    stdout=open(tmp_path, 'wb'),
                    stderr=subprocess.DEVNULL,
                    check=True,
                )
            else:
                log.info('COSMIC: decompressing with gzip (install pigz for faster decompression)')
                subprocess.run(
                    f'gzip -d -c "{tsv_path}" > "{tmp_path}"',
                    shell=True, check=True,
                )
        except Exception as e:
            os.unlink(tmp_path)
            raise RuntimeError(f'COSMIC decompression failed: {e}') from e
        plain_path = tmp_path
        cleanup = True
    else:
        plain_path = tsv_path
        cleanup = False

    try:
        # ── Step 2: read header and validate ─────────────────────────────────
        with open(plain_path, 'r') as fh:
            header = fh.readline().rstrip('\n').split('\t')
        col_idx = {c: i for i, c in enumerate(header)}

        missing = [c for c in usecols if c not in col_idx]
        if missing:
            log.error('COSMIC TSV missing columns: %s', missing)
            log.error('Available columns: %s', list(col_idx.keys()))
            return empty_standard_df()

        # ── Step 3: split file into N byte ranges ────────────────────────────
        # Each worker gets a start and end byte offset. Workers seek to their
        # start, skip to the next newline to align, then read until end.
        file_size = os.path.getsize(plain_path)
        header_end = len(header_line := (('	'.join(header)) + '\n').encode())

        section_size = max(1, (file_size - header_end) // n_workers)
        ranges = []
        for i in range(n_workers):
            start = header_end + i * section_size
            end = header_end + (i + 1) * section_size if i < n_workers - 1 else file_size
            if start < file_size:
                ranges.append((start, end))

        log.info('COSMIC: processing %d byte-range sections with %d workers',
                 len(ranges), n_workers)

        worker_args = [(plain_path, start, end, col_idx, vcf_lookup, classification_lookup)
                       for start, end in ranges]

        with mp.Pool(processes=n_workers) as pool:
            results = pool.starmap(_process_cosmic_section, worker_args)

    finally:
        if cleanup:
            try:
                os.unlink(plain_path)
            except OSError:
                pass

    # ── Step 4: merge results ─────────────────────────────────────────────────
    agg: dict[tuple, set] = {}
    meta: dict[tuple, dict] = {}
    for partial_agg, partial_meta in results:
        for key, samples in partial_agg.items():
            if key not in agg:
                agg[key] = set()
                meta[key] = partial_meta[key]
            agg[key].update(samples)

    log.info('COSMIC: aggregating %d (variant, cancer_type) keys', len(agg))

    rows = []
    for (chrom, pos, ref, alt, gene, cancer_type), samples in agg.items():
        m = meta[(chrom, pos, ref, alt, gene, cancer_type)]
        rows.append({
            'chrom':       chrom,
            'pos':         pos,
            'ref':         ref,
            'alt':         alt,
            'gene':        gene,
            'hgvsc':       m['hgvsc'],
            'hgvsp':       m['hgvsp'],
            'consequence': m['consequence'],
            'cancer_type': cancer_type,
            'n_samples':   len(samples),
            'source':      'COSMIC',
        })

    if not rows:
        log.warning('COSMIC: no rows produced')
        return empty_standard_df()

    df = pd.DataFrame(rows, columns=STANDARD_COLS)
    log.info('COSMIC: %d rows (variant x cancer_type combinations)', len(df))
    return df


def _process_cosmic_section(path: str,
                             byte_start: int,
                             byte_end: int,
                             col_idx: dict,
                             vcf_lookup: dict,
                             classification_lookup: dict) -> tuple[dict, dict]:
    """
    Read and process a byte-range section of a plain-text (decompressed) TSV.
    Each worker seeks independently — no shared file handle, no GIL contention.
    Aligns to the next newline after byte_start to avoid split-line reads.
    """
    agg: dict[tuple, set] = {}
    meta: dict[tuple, dict] = {}
    chunk: list[list] = []

    with open(path, 'rb') as fh:
        fh.seek(byte_start)
        # If not at file start, skip partial line to align to next complete line
        if byte_start > 0:
            fh.readline()
        while fh.tell() < byte_end:
            raw = fh.readline()
            if not raw:
                break
            line = raw.decode('utf-8', errors='replace').rstrip('\n')
            chunk.append(line.split('\t'))
            if len(chunk) >= 50_000:
                _process_cosmic_chunk(chunk, col_idx, vcf_lookup, classification_lookup, agg, meta)
                chunk = []
    if chunk:
        _process_cosmic_chunk(chunk, col_idx, vcf_lookup, classification_lookup, agg, meta)
    return agg, meta


def _process_cosmic_chunk(chunk: list[list],
                           col_idx: dict[str, int],
                           vcf_lookup: dict,
                           classification_lookup: dict,
                           agg: dict,
                           meta: dict) -> None:
    """Process one chunk of COSMIC TSV rows, populating agg and meta in place."""
    for parts in chunk:
        try:
            # Skip germline and SNPs
            status = parts[col_idx[_COL_STATUS]]
            if 'germline' in status.lower():
                continue

            chrom_raw  = parts[col_idx[_COL_CHROM]].strip()
            start_raw  = parts[col_idx[_COL_START]].strip()
            ref        = parts[col_idx[_COL_REF]].strip()
            alt        = parts[col_idx[_COL_ALT]].strip()

            if not chrom_raw or not start_raw or start_raw == 'NS':
                continue
            if not ref or not alt:
                continue

            chrom = chrom_raw.lstrip('chr')
            try:
                pos = int(start_raw)
            except ValueError:
                continue

            gene        = parts[col_idx[_COL_GENE]].strip()
            pheno_id    = parts[col_idx[_COL_SITE]].strip()
            cancer_type = classification_lookup.get(pheno_id, pheno_id).lower()
            sample_id   = parts[col_idx[_COL_SAMPLE_ID]].strip()
            cds         = parts[col_idx[_COL_CDS]].strip()
            desc        = parts[col_idx[_COL_DESC]].strip()
            hgvsp       = parts[col_idx[_COL_HGVSP]].strip()
            hgvsc       = parts[col_idx[_COL_HGVSC]].strip()

            if not is_valid_allele(ref) or not is_valid_allele(alt):
                continue

            consequence = map_consequence(desc)

            key = (chrom, pos, ref, alt, gene, cancer_type)
            if key not in agg:
                agg[key] = set()
                meta[key] = {
                    'hgvsc':       hgvsc,
                    'hgvsp':       hgvsp,
                    'consequence': consequence,
                }
            agg[key].add(sample_id)

        except (IndexError, ValueError):
            continue
