#!/usr/bin/env python3
"""
mane_audit.py
Count how many variants in the whitelist are not on MANE Select transcripts.
Runs without any API calls -- just the MANE table download.

Matches on Ensembl transcript (ENST) IDs, which is what the whitelist hgvsc
column contains.

Usage:
  python3 mane_audit.py --whitelist output/pan_cancer_whitelist_GRCh38.tsv.gz
"""
import argparse
import gzip
import os
import urllib.request
import pandas as pd

_MANE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
    "MANE.GRCh38.v1.5.summary.txt.gz"
)

def load_mane_lookup(cache="mane_select.tsv.gz"):
    if not os.path.exists(cache):
        print("Downloading MANE Select table...")
        urllib.request.urlretrieve(_MANE_URL, cache)
    with gzip.open(cache, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", dtype=str, low_memory=False)

    print(f"MANE table columns: {list(df.columns)}")

    # Identify columns
    gene_col = next((c for c in df.columns if "symbol" in c.lower()), None)
    enst_col = next((c for c in df.columns if "ensembl_nuc" in c.lower()), None)
    nm_col   = next((c for c in df.columns if "refseq_nuc" in c.lower()), None)
    type_col = next((c for c in df.columns if "mane_status" in c.lower()), None)

    if not gene_col or not enst_col:
        raise RuntimeError(
            f"Cannot find required columns. Available: {list(df.columns)}"
        )

    # Keep only MANE Select rows
    if type_col:
        df = df[df[type_col].str.contains("Select", case=False, na=False)]

    # Build two lookups keyed by gene symbol:
    #   enst_lookup: gene -> ENST base ID (no version)
    #   nm_lookup:   gene -> NM base ID (no version)
    enst_lookup = {}
    nm_lookup   = {}
    for _, row in df.iterrows():
        gene = str(row[gene_col]).strip()
        enst = str(row[enst_col]).strip().split(".")[0]
        if gene and enst.startswith("ENST"):
            enst_lookup[gene] = enst
        if nm_col:
            nm = str(row[nm_col]).strip().split(".")[0]
            if gene and nm.startswith("NM_"):
                nm_lookup[gene] = nm

    print(f"MANE Select: {len(enst_lookup)} genes loaded (ENST lookup)")
    return enst_lookup, nm_lookup


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--whitelist", required=True)
    ap.add_argument("--mane-cache", default="mane_select.tsv.gz")
    args = ap.parse_args()

    enst_lookup, nm_lookup = load_mane_lookup(args.mane_cache)

    print(f"Loading whitelist: {args.whitelist}")
    df = pd.read_csv(args.whitelist, sep="\t", dtype=str, low_memory=False)
    print(f"Total variants: {len(df)}")

    no_hgvsc         = 0
    already_mane     = 0
    needs_remap      = 0
    gene_not_in_mane = 0
    unrecognised_tx  = 0

    for _, row in df.iterrows():
        hgvsc = str(row.get("hgvsc", "")).strip()
        gene  = str(row.get("gene",  "")).strip()

        if not hgvsc or hgvsc == "nan":
            no_hgvsc += 1
            continue

        if ":" not in hgvsc:
            unrecognised_tx += 1
            continue

        # Strip version number from transcript ID
        tx = hgvsc.split(":")[0].strip().split(".")[0]

        if gene not in enst_lookup and gene not in nm_lookup:
            gene_not_in_mane += 1
            continue

        mane_enst = enst_lookup.get(gene, "")
        mane_nm   = nm_lookup.get(gene, "")

        if tx == mane_enst or tx == mane_nm:
            already_mane += 1
        else:
            needs_remap += 1

    print(f"\nResults:")
    print(f"  No HGVSc annotation:          {no_hgvsc:>8}")
    print(f"  Unrecognised transcript fmt:  {unrecognised_tx:>8}")
    print(f"  Gene not in MANE table:       {gene_not_in_mane:>8}")
    print(f"  Already on MANE Select:       {already_mane:>8}")
    print(f"  Needs remapping:              {needs_remap:>8}")
    print(f"  Total:                        {len(df):>8}")


if __name__ == "__main__":
    main()
