#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Panel annotation: tag whitelist variants with diagnostic panel membership
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
annotate_panels.py — Annotate the pan-cancer whitelist with panel membership.

For each BED file in the panels/ directory, variants whose (chrom, pos) fall
within any interval are tagged with that panel name.

Usage:
    python annotate_panels.py \
        --whitelist output/pan_cancer_whitelist_GRCh38.tsv.gz \
        --panels-dir panels/ \
        --output output/pan_cancer_whitelist_GRCh38.annotated.tsv.gz

    # Filter to a specific panel after annotation:
    zcat output/pan_cancer_whitelist_GRCh38.annotated.tsv.gz \
        | awk -F'\\t' 'NR==1 || $0~/lymphoma/'

BED file naming:
    panels/lymphoma.bed   -> panel name = "lymphoma"
    panels/myeloid.bed    -> panel name = "myeloid"

BED format (0-based, half-open, as standard):
    chr1    100000    200000
    chr7    116411708 116411899
    ...

Chromosomes may use either 'chr' prefix or plain integers; the script
normalises both to match the whitelist format.
"""

import argparse
import bisect
import gzip
import logging
import os
import sys
from collections import defaultdict

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [annotate_panels] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# BED loading
# ---------------------------------------------------------------------------

def load_bed(path: str) -> dict:
    """
    Load a BED file into a dict of {chrom: [(start, end), ...]} intervals.
    Coordinates are 0-based, half-open (standard BED).
    Chromosomes are normalised to strip the 'chr' prefix.
    Intervals are sorted by start position for binary search.
    """
    intervals = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue
            chrom = parts[0].lstrip('chr')
            try:
                start = int(parts[1])
                end   = int(parts[2])
            except ValueError:
                continue
            intervals[chrom].append((start, end))
    # Sort intervals by start for binary search
    for chrom in intervals:
        intervals[chrom].sort()
    return dict(intervals)


def in_bed(chrom: str, pos: int, intervals: dict) -> bool:
    """
    Return True if pos (1-based VCF/TSV position) falls within any interval
    for chrom in the BED intervals dict.
    BED is 0-based half-open, so a 1-based pos P overlaps [start, end) when:
        start < P <= end  (equivalently: start <= P-1 < end)
    Uses binary search for O(log N) per query.
    """
    chrom = str(chrom).lstrip('chr')
    if chrom not in intervals:
        return False
    p0 = pos - 1  # convert to 0-based
    ivs = intervals[chrom]
    # Find rightmost interval whose start <= p0
    idx = bisect.bisect_right(ivs, (p0, float('inf'))) - 1
    # Check this interval and possibly the one before (overlapping intervals)
    for i in range(max(0, idx), min(idx + 2, len(ivs))):
        if ivs[i][0] <= p0 < ivs[i][1]:
            return True
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--whitelist',  required=True, help='Master whitelist TSV.gz')
    parser.add_argument('--panels-dir', default='panels', help='Directory containing .bed files (default: panels/)')
    parser.add_argument('--output',     required=True, help='Output annotated TSV.gz')
    args = parser.parse_args()

    # Load panel BED files
    if not os.path.isdir(args.panels_dir):
        log.error('Panels directory not found: %s', args.panels_dir)
        sys.exit(1)

    bed_files = [f for f in os.listdir(args.panels_dir) if f.endswith('.bed')]
    if not bed_files:
        log.error('No .bed files found in: %s', args.panels_dir)
        sys.exit(1)

    panels = {}
    for fname in sorted(bed_files):
        name = fname.replace('.bed', '')
        path = os.path.join(args.panels_dir, fname)
        panels[name] = load_bed(path)
        total = sum(len(v) for v in panels[name].values())
        log.info('Loaded panel: %-30s  %d intervals', name, total)

    # Load whitelist
    log.info('Loading whitelist: %s', args.whitelist)
    df = pd.read_csv(args.whitelist, sep='\t', dtype=str, low_memory=False)
    log.info('Whitelist rows: %d', len(df))

    if 'pos' not in df.columns or 'chrom' not in df.columns:
        log.error("Whitelist must contain 'chrom' and 'pos' columns.")
        sys.exit(1)

    df['pos'] = pd.to_numeric(df['pos'], errors='coerce')

    # Annotate
    log.info('Annotating variants against %d panel(s)...', len(panels))

    def get_panels(row):
        hits = []
        try:
            pos = int(row['pos'])
        except (ValueError, TypeError):
            return ''
        for panel_name, intervals in panels.items():
            if in_bed(row['chrom'], pos, intervals):
                hits.append(panel_name)
        return ';'.join(hits)

    df['panels'] = df.apply(get_panels, axis=1)

    n_tagged = (df['panels'] != '').sum()
    log.info('Variants overlapping at least one panel: %d / %d', n_tagged, len(df))

    for panel_name in panels:
        n = df['panels'].str.contains(panel_name, regex=False).sum()
        log.info('  %-30s  %d variants', panel_name, n)

    # Write output
    log.info('Writing annotated whitelist: %s', args.output)
    df.to_csv(args.output, sep='\t', index=False, compression='gzip')
    log.info('Done.')


if __name__ == '__main__':
    main()
