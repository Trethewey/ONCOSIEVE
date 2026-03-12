#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Tool: ClinVar VEP transcript annotation
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
clinvar_vep_annotate.py

Adds transcript-level HGVSc annotation for ClinVar variants that lack it,
using the Ensembl VEP REST API. Prefers MANE Select transcripts.

Adds a new column 'hgvsc_vep' — does not overwrite the original 'hgvsc'.

Checkpoints progress after each batch so runs can be safely interrupted
and resumed.

Usage:
  python3 clinvar_vep_annotate.py \\
      --whitelist output/pan_cancer_whitelist_GRCh38.tsv.gz \\
      --output    output/pan_cancer_whitelist_GRCh38.annotated.tsv.gz

Optional:
  --batch-size   INT   Variants per VEP API call (default: 200, max: 200)
  --checkpoint   PATH  Checkpoint file path (default: vep_checkpoint.json)
  --delay        FLOAT Seconds between batches (default: 0.5)
  --all-sources        Annotate all sources, not just ClinVar
"""

import argparse
import gzip
import json
import os
import time
import urllib.error
import urllib.request

import pandas as pd

_VEP_URL  = "https://rest.ensembl.org/vep/human/region"
_HEADERS  = {
    "Content-Type": "application/json",
    "Accept":       "application/json",
}
_MAX_RETRIES   = 5
_BACKOFF_BASE  = 2   # seconds; doubles on each retry


# ---------------------------------------------------------------------------
# VEP API helpers
# ---------------------------------------------------------------------------

def _build_vep_input(row: pd.Series) -> str | None:
    """
    Convert a whitelist row to VEP region string format.
    Format: 'chrom start end allele strand'
    For SNVs and substitutions: start == end == pos
    For insertions: start = pos, end = pos - 1 (VEP convention)
    For deletions: start = pos + 1, end = pos + len(ref) - 1
    """
    chrom = str(row["chrom"]).replace("chr", "")
    pos   = int(row["pos"])
    ref   = str(row["ref"])
    alt   = str(row["alt"])

    if ref == "-" or alt == "-":
        return None  # skip symbolic alleles

    if len(ref) == 1 and len(alt) == 1:
        # SNV
        return f"{chrom} {pos} {pos} {ref}/{alt} 1"
    elif len(ref) == 1 and len(alt) > 1:
        # Insertion: VEP uses pos/pos-1 convention
        start = pos
        end   = pos - 1
        return f"{chrom} {start} {end} -/{alt[1:]} 1"
    elif len(ref) > 1 and len(alt) == 1:
        # Deletion
        start = pos + 1
        end   = pos + len(ref) - 1
        return f"{chrom} {start} {end} {ref[1:]}/- 1"
    else:
        # MNV / complex — pass as-is
        return f"{chrom} {pos} {pos + len(ref) - 1} {ref}/{alt} 1"


def _vep_batch(variants: list[str]) -> list[dict]:
    """
    POST a batch of VEP region strings. Returns list of VEP result dicts.
    Retries with exponential backoff on transient errors.
    """
    payload = json.dumps({"regions": variants}).encode("utf-8")
    params  = "?MANE=1&hgvs=1&transcript_version=1"
    url     = _VEP_URL + params

    for attempt in range(1, _MAX_RETRIES + 1):
        try:
            req  = urllib.request.Request(url, data=payload, headers=_HEADERS, method="POST")
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code in (429, 503) and attempt < _MAX_RETRIES:
                wait = _BACKOFF_BASE ** attempt
                print(f"  Rate limited (HTTP {e.code}), retrying in {wait}s "
                      f"(attempt {attempt}/{_MAX_RETRIES})")
                time.sleep(wait)
            else:
                print(f"  VEP HTTP error {e.code}: {e.reason}")
                return []
        except Exception as e:
            if attempt < _MAX_RETRIES:
                wait = _BACKOFF_BASE ** attempt
                print(f"  VEP error: {e}, retrying in {wait}s")
                time.sleep(wait)
            else:
                print(f"  VEP failed after {_MAX_RETRIES} attempts: {e}")
                return []
    return []


def _extract_hgvsc(vep_result: dict) -> str:
    """
    Extract HGVSc from a single VEP result dict.
    Preference order:
      1. MANE Select transcript
      2. MANE Plus Clinical transcript
      3. First transcript consequence with HGVSc
    """
    transcripts = vep_result.get("transcript_consequences", [])
    if not transcripts:
        return ""

    # Sort by preference
    mane_select  = [t for t in transcripts if t.get("mane_select")]
    mane_plus    = [t for t in transcripts if t.get("mane_plus_clinical")]
    any_with_hgvsc = [t for t in transcripts if t.get("hgvsc")]

    for candidate_list in (mane_select, mane_plus, any_with_hgvsc):
        for t in candidate_list:
            hgvsc = t.get("hgvsc", "")
            if hgvsc:
                return hgvsc

    return ""


# ---------------------------------------------------------------------------
# Checkpoint helpers
# ---------------------------------------------------------------------------

def _load_checkpoint(path: str) -> dict:
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {"completed_batches": [], "results": {}}


def _save_checkpoint(path: str, state: dict) -> None:
    with open(path, "w") as f:
        json.dump(state, f)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--whitelist",   required=True)
    ap.add_argument("--output",      required=True)
    ap.add_argument("--batch-size",  type=int,   default=200)
    ap.add_argument("--checkpoint",  default="vep_checkpoint.json")
    ap.add_argument("--delay",       type=float, default=0.5)
    ap.add_argument("--all-sources", action="store_true")
    args = ap.parse_args()

    batch_size = min(args.batch_size, 200)  # VEP API hard limit

    print(f"Loading whitelist: {args.whitelist}")
    df = pd.read_csv(args.whitelist, sep="\t", dtype=str, low_memory=False)
    print(f"Total variants: {len(df)}")

    # Initialise hgvsc_vep column if not present
    if "hgvsc_vep" not in df.columns:
        df["hgvsc_vep"] = ""

    # Select variants to annotate
    if args.all_sources:
        mask = df["hgvsc"].isna() | (df["hgvsc"] == "")
    else:
        mask = (
            df["sources"].str.contains("ClinVar", case=False, na=False) &
            (df["hgvsc"].isna() | (df["hgvsc"] == ""))
        )

    to_annotate = df[mask].copy()
    print(f"Variants needing VEP annotation: {len(to_annotate)}")

    if len(to_annotate) == 0:
        print("Nothing to annotate.")
        df.to_csv(args.output, sep="\t", index=False, compression="gzip")
        return

    # Load checkpoint
    state = _load_checkpoint(args.checkpoint)
    results_cache: dict = state.get("results", {})
    completed: set     = set(state.get("completed_batches", []))

    # Build batches
    indices = to_annotate.index.tolist()
    batches = [
        indices[i:i + batch_size]
        for i in range(0, len(indices), batch_size)
    ]
    n_batches = len(batches)
    print(f"Batches: {n_batches} x up to {batch_size} variants")

    # Process batches
    for batch_num, batch_idx in enumerate(batches):
        batch_key = str(batch_num)
        if batch_key in completed:
            print(f"  Batch {batch_num + 1}/{n_batches}: skipped (cached)")
            continue

        batch_rows = df.loc[batch_idx]
        vep_inputs = []
        idx_map    = {}  # vep input string -> dataframe index

        for idx, row in batch_rows.iterrows():
            vep_str = _build_vep_input(row)
            if vep_str:
                vep_inputs.append(vep_str)
                idx_map[vep_str] = idx

        if not vep_inputs:
            completed.add(batch_key)
            continue

        print(f"  Batch {batch_num + 1}/{n_batches}: {len(vep_inputs)} variants", end="", flush=True)
        vep_results = _vep_batch(vep_inputs)

        n_annotated = 0
        for result in vep_results:
            input_str = result.get("input", "")
            hgvsc     = _extract_hgvsc(result)
            if hgvsc:
                results_cache[input_str] = hgvsc
                n_annotated += 1

        # Apply results to dataframe
        for vep_str, idx in idx_map.items():
            hgvsc = results_cache.get(vep_str, "")
            if hgvsc:
                df.at[idx, "hgvsc_vep"] = hgvsc

        print(f" -> {n_annotated} annotated")
        completed.add(batch_key)

        # Save checkpoint after every batch
        _save_checkpoint(args.checkpoint, {
            "completed_batches": list(completed),
            "results": results_cache,
        })

        time.sleep(args.delay)

    # Write output
    print(f"\nWriting output: {args.output}")
    df.to_csv(args.output, sep="\t", index=False, compression="gzip")

    n_filled = (df["hgvsc_vep"] != "").sum()
    print(f"Variants with hgvsc_vep filled: {n_filled} / {len(to_annotate)}")
    print("Done.")


if __name__ == "__main__":
    main()
