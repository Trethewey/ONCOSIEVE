#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — annotate_revel.py
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
#
# Standalone post-processing script. Run manually after the main pipeline.
#
# Adds a revel_score column to the final whitelist TSV by joining against
# the REVEL v1.3 precomputed score file (revel_with_transcript_ids).
#
# Join key: chrom + pos + ref + alt (GRCh38 coordinates).
# Only missense variants will have a score. All others will be empty.
#
# Usage:
#   python3 tools/annotate_revel.py \
#       --whitelist output/pan_cancer_whitelist_GRCh38_annotated.tsv.gz \
#       --revel     data/REVEL/revel_with_transcript_ids \
#       --out       output/pan_cancer_whitelist_GRCh38_revel.tsv.gz
# =============================================================================

import argparse
import os
import sys

import pandas as pd
import polars as pl


def main():
    ap = argparse.ArgumentParser(
        description='Annotate ONCOSIEVE whitelist with REVEL scores.'
    )
    ap.add_argument('--whitelist', required=True,
                    help='Whitelist TSV (output of post_process_whitelist.py)')
    ap.add_argument('--revel', required=True,
                    help='REVEL score file (revel_with_transcript_ids)')
    ap.add_argument('--out', required=True,
                    help='Output TSV path (.tsv.gz)')
    args = ap.parse_args()

    if not os.path.exists(args.whitelist):
        sys.exit(f'ERROR: whitelist not found: {args.whitelist}')
    if not os.path.exists(args.revel):
        sys.exit(f'ERROR: REVEL file not found: {args.revel}')

    print(f'Loading whitelist: {args.whitelist}')
    wl = pd.read_csv(args.whitelist, sep='\t', dtype=str, low_memory=False)
    print(f'  {len(wl)} variants, {len(wl.columns)} columns')

    # Normalise chrom: strip chr prefix to match REVEL format
    wl['_chrom_key'] = wl['chrom'].str.replace(r'^chr', '', regex=True)
    wl['_pos_key']   = pd.to_numeric(wl['pos'], errors='coerce')

    print(f'Loading REVEL scores: {args.revel}')
    print('  This may take 1-2 minutes...')
    revel = pl.read_csv(
        args.revel,
        separator=',',
        infer_schema_length=10000,
        null_values=['.'],
        schema_overrides={
            'chr':        pl.Utf8,
            'grch38_pos': pl.Int64,
            'ref':        pl.Utf8,
            'alt':        pl.Utf8,
            'REVEL':      pl.Float64,
        },
    ).select(['chr', 'grch38_pos', 'ref', 'alt', 'REVEL']).filter(
        pl.col('grch38_pos').is_not_null()
    )

    print(f'  {len(revel):,} REVEL entries loaded')

    # Deduplicate: keep max score per position
    revel = (
        revel
        .group_by(['chr', 'grch38_pos', 'ref', 'alt'])
        .agg(pl.col('REVEL').max().alias('revel_score'))
    ).to_pandas()

    revel['grch38_pos'] = revel['grch38_pos'].astype('Int64')

    print('Joining...')
    merged = wl.merge(
        revel,
        left_on  = ['_chrom_key', '_pos_key', 'ref', 'alt'],
        right_on = ['chr', 'grch38_pos', 'ref', 'alt'],
        how      = 'left',
    )

    # Drop helper columns
    merged = merged.drop(columns=['_chrom_key', '_pos_key', 'chr', 'grch38_pos'])

    # Round and fill
    merged['revel_score'] = merged['revel_score'].round(4).astype(str)
    merged.loc[merged['revel_score'] == 'nan', 'revel_score'] = ''

    n_scored = (merged['revel_score'] != '').sum()
    print(f'  Variants with REVEL score: {n_scored} / {len(merged)}')

    print(f'Writing output: {args.out}')
    merged.to_csv(args.out, sep='\t', index=False,
                  compression='gzip' if args.out.endswith('.gz') else None)
    print(f'  Columns: {list(merged.columns)}')
    print('Done.')


if __name__ == '__main__':
    main()
