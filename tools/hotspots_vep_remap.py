#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — pan-cancer variant curation and rescue tool
# Tool: CancerHotspots VEP remapping to GRCh38 / MANE Select
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
# =============================================================================

"""
hotspots_vep_remap.py

Reads CancerHotspots V1, V2, and V3 residue-level files, constructs
protein HGVS notation for every variant amino acid, and maps each to
GRCh38 genomic coordinates via the Ensembl VEP REST API.

MANE Select transcript matching is done against a local MANE table
(data/reference/mane_select.tsv.gz) rather than relying on the VEP
MANE flag, which is not reliably returned by all VEP versions.

Variants with no MANE Select result are excluded.
Checkpoints after each batch so runs can be safely interrupted and resumed.

Output: data/hotspots/hotspots_grch38_mane.tsv

Usage:
  python3 hotspots_vep_remap.py
  python3 hotspots_vep_remap.py --resume
"""

import argparse
import gzip
import json
import os
import re
import time
import urllib.error
import urllib.request

import pandas as pd

_VEP_URL      = "https://rest.ensembl.org/vep/human/hgvs"
_HEADERS      = {"Content-Type": "application/json", "Accept": "application/json"}
_MAX_RETRIES  = 5
_BACKOFF_BASE = 2
_BATCH_SIZE   = 50

_AA1TO3 = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Ter', 'X': 'Ter',
}


def _aa1to3(aa: str) -> str:
    return _AA1TO3.get(aa.upper(), aa)


# ---------------------------------------------------------------------------
# MANE lookup
# ---------------------------------------------------------------------------

def load_mane_lookup(mane_path: str) -> dict:
    """
    Returns dict: gene_symbol -> ENST base ID (no version)
    e.g. {'NRAS': 'ENST00000369535', 'BRAF': 'ENST00000288602', ...}
    """
    with gzip.open(mane_path, "rt") as f:
        df = pd.read_csv(f, sep="\t", dtype=str)

    gene_col  = next(c for c in df.columns if "symbol" in c.lower())
    enst_col  = next(c for c in df.columns if "ensembl_nuc" in c.lower())
    status_col = next(c for c in df.columns if "mane_status" in c.lower())

    df = df[df[status_col].str.contains("Select", case=False, na=False)]

    lookup = {}
    for _, row in df.iterrows():
        gene = str(row[gene_col]).strip()
        enst = str(row[enst_col]).strip().split(".")[0]
        if gene and enst.startswith("ENST"):
            lookup[gene] = enst

    print(f"MANE Select lookup: {len(lookup)} genes")
    return lookup


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def _parse_variant_aas(s: str) -> list[str]:
    """
    Parse variant amino acid strings. Handles formats:
      'E:520|K:33|R:4'   -> ['E', 'K', 'R']   (pipe-delimited with counts)
      'R,K,L'            -> ['R', 'K', 'L']   (comma-delimited)
    Single-letter AA codes only — counts are stripped.
    """
    s = str(s).strip()
    if not s or s == 'nan':
        return []
    aas = []
    if '|' in s:
        for part in s.split('|'):
            part = part.strip()
            if not part:
                continue
            # Strip count suffix: 'E:520' -> 'E'
            aa = part.split(':')[0].strip()
            if aa:
                aas.append(aa)
    else:
        for part in s.split(','):
            aa = part.strip().split(':')[0].strip()
            if aa:
                aas.append(aa)
    return aas


def read_v1(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', dtype=str)
    df.columns = df.columns.str.strip()
    rows = []
    for _, row in df.iterrows():
        gene    = str(row.get('Hugo Symbol', '')).strip()
        codon   = str(row.get('Codon', '')).strip()
        ref_aa  = re.sub(r'\d+', '', codon).strip()
        pos     = re.sub(r'\D', '', codon).strip()
        var_aas = _parse_variant_aas(row.get('Variant Amino Acid', ''))
        count   = str(row.get('Tumor Count', '0')).strip()
        qvalue  = str(row.get('Q-value', '')).strip()
        if not gene or not pos:
            continue
        for alt_aa in var_aas:
            rows.append({'gene': gene, 'ref_aa': ref_aa, 'pos': pos,
                         'alt_aa': alt_aa, 'n_samples': count,
                         'qvalue': qvalue, 'version': 'V1'})
    return pd.DataFrame(rows)


def read_v2(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', dtype=str)
    df.columns = df.columns.str.strip()
    rows = []
    for _, row in df.iterrows():
        gene   = str(row.get('Hugo_Symbol', '')).strip()
        pos    = str(row.get('Amino_Acid_Position', '')).strip()
        ref_aa = str(row.get('Reference_Amino_Acid', '')).strip()
        if ':' in ref_aa:
            ref_aa = ref_aa.split(':')[0].strip()
        var_aas = _parse_variant_aas(row.get('Variant_Amino_Acid', ''))
        count   = str(row.get('Mutation_Count', '0')).strip()
        qvalue  = str(row.get('qvalue', '')).strip()
        if not gene or not pos:
            continue
        for alt_aa in var_aas:
            rows.append({'gene': gene, 'ref_aa': ref_aa, 'pos': pos,
                         'alt_aa': alt_aa, 'n_samples': count,
                         'qvalue': qvalue, 'version': 'V2'})
    return pd.DataFrame(rows)


def read_v3(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', dtype=str)
    df.columns = df.columns.str.strip()
    rows = []
    for _, row in df.iterrows():
        gene   = str(row.get('Hugo_Symbol', '')).strip()
        codon  = str(row.get('Codon', '')).strip()
        ref_aa = re.sub(r'\D', '', codon) and re.sub(r'\d+', '', codon).strip()
        pos    = re.sub(r'\D', '', codon).strip()
        count  = str(row.get('# mut in MSK', '0')).strip()
        qvalue = str(row.get('Q value', '')).strip()
        if not gene or not pos:
            continue
        # V3 has no individual alt AAs listed — skip for now
        # These will be picked up from V2 overlap where possible
    return pd.DataFrame(rows) if rows else pd.DataFrame(
        columns=['gene','ref_aa','pos','alt_aa','n_samples','qvalue','version'])


# ---------------------------------------------------------------------------
# VEP helpers
# ---------------------------------------------------------------------------

def _build_hgvs(gene: str, ref_aa: str, pos: str, alt_aa: str) -> str | None:
    if not gene or not pos or not ref_aa or not alt_aa or alt_aa == '?':
        return None
    ref3 = _aa1to3(ref_aa)
    alt3 = _aa1to3(alt_aa)
    if not ref3 or not alt3:
        return None
    return f"{gene}:p.{ref3}{pos}{alt3}"


def _vep_single(hgvs: str) -> dict | None:
    """Query VEP for a single HGVS string. Returns result dict or None."""
    payload = json.dumps({"hgvs_notations": [hgvs]}).encode("utf-8")
    url     = _VEP_URL + "?hgvs=1&transcript_version=1"
    try:
        req = urllib.request.Request(
            url, data=payload, headers=_HEADERS, method="POST"
        )
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data[0] if data else None
    except Exception:
        return None


def _vep_batch(hgvs_list: list[str]) -> list[dict]:
    payload = json.dumps({"hgvs_notations": hgvs_list}).encode("utf-8")
    url     = _VEP_URL + "?hgvs=1&transcript_version=1"

    for attempt in range(1, _MAX_RETRIES + 1):
        try:
            req = urllib.request.Request(
                url, data=payload, headers=_HEADERS, method="POST"
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code == 400:
                # Bad input in batch -- fall back to one-by-one to isolate bad entries
                print(f"    HTTP 400 on batch, falling back to single queries")
                results = []
                for hgvs in hgvs_list:
                    res = _vep_single(hgvs)
                    if res:
                        results.append(res)
                    time.sleep(0.1)
                return results
            elif e.code in (429, 503) and attempt < _MAX_RETRIES:
                wait = _BACKOFF_BASE ** attempt
                print(f"    Rate limited ({e.code}), retrying in {wait}s")
                time.sleep(wait)
            else:
                print(f"    VEP HTTP error {e.code}")
                return []
        except Exception as e:
            if attempt < _MAX_RETRIES:
                wait = _BACKOFF_BASE ** attempt
                print(f"    VEP error: {e}, retrying in {wait}s")
                time.sleep(wait)
            else:
                print(f"    VEP failed: {e}")
                return []
    return []


def _extract_mane_coords(vep_result: dict, mane_lookup: dict) -> dict | None:
    """
    Match VEP transcript consequences against local MANE Select lookup.
    Returns coord dict if MANE Select transcript found, else None.
    """
    gene = None
    transcripts = vep_result.get("transcript_consequences", [])

    # Determine gene from input string e.g. "NRAS:p.Gln61Arg"
    inp = vep_result.get("input", "")
    if ":" in inp:
        gene = inp.split(":")[0].strip()

    if not gene or gene not in mane_lookup:
        return None

    mane_enst = mane_lookup[gene]  # base ID without version

    for t in transcripts:
        tx_id = str(t.get("transcript_id", "")).split(".")[0]
        if tx_id == mane_enst:
            return {
                "chrom":       str(vep_result.get("seq_region_name", "")),
                "pos":         str(vep_result.get("start", "")),
                "ref":         str(vep_result.get("allele_string", "/").split("/")[0]),
                "alt":         str(vep_result.get("allele_string", "/").split("/")[-1]),
                "hgvsc":       t.get("hgvsc", ""),
                "hgvsp":       t.get("hgvsp", ""),
                "consequence": t.get("consequence_terms", ["unknown"])[0],
            }
    return None


# ---------------------------------------------------------------------------
# Checkpoint helpers
# ---------------------------------------------------------------------------

def _load_checkpoint(path: str) -> dict:
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {"done": [], "results": {}}


def _save_checkpoint(path: str, state: dict) -> None:
    with open(path, "w") as f:
        json.dump(state, f)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hotspot-dir",  default="data/hotspots")
    ap.add_argument("--mane-table",   default="data/reference/mane_select.tsv.gz")
    ap.add_argument("--output",       default="data/hotspots/hotspots_grch38_mane.tsv")
    ap.add_argument("--checkpoint",   default="data/hotspots/vep_hotspot_checkpoint.json")
    ap.add_argument("--delay",        type=float, default=0.5)
    ap.add_argument("--resume",       action="store_true")
    args = ap.parse_args()

    # Load MANE lookup
    mane_lookup = load_mane_lookup(args.mane_table)

    # Read hotspot files
    print("Reading hotspot files...")
    d = args.hotspot_dir
    dfs = []
    for version, reader, fname in [
        ("V1", read_v1, "hotspots_v1.txt"),
        ("V2", read_v2, "hotspots_v2.txt"),
        ("V3", read_v3, "hotspots_v3.txt"),
    ]:
        path = os.path.join(d, fname)
        if os.path.exists(path):
            df = reader(path)
            print(f"  {version}: {len(df)} variant-residue entries")
            dfs.append(df)

    all_hotspots = pd.concat(dfs, ignore_index=True)

    # Build HGVS strings
    all_hotspots["hgvs_input"] = all_hotspots.apply(
        lambda r: _build_hgvs(r["gene"], r["ref_aa"], r["pos"], r["alt_aa"]),
        axis=1,
    )
    all_hotspots = all_hotspots[all_hotspots["hgvs_input"].notna()].copy()
    all_hotspots = all_hotspots.drop_duplicates(subset=["hgvs_input"], keep="last")
    print(f"  Total unique HGVS inputs: {len(all_hotspots)}")

    # Checkpoint
    if not args.resume and os.path.exists(args.checkpoint):
        os.remove(args.checkpoint)
    state    = _load_checkpoint(args.checkpoint)
    done_set = set(state.get("done", []))
    results  = state.get("results", {})

    # Batches
    hgvs_list = all_hotspots["hgvs_input"].tolist()
    batches   = [hgvs_list[i:i+_BATCH_SIZE]
                 for i in range(0, len(hgvs_list), _BATCH_SIZE)]
    n_batches = len(batches)
    print(f"Batches: {n_batches} x up to {_BATCH_SIZE} variants")

    for batch_num, batch in enumerate(batches):
        batch_key = str(batch_num)
        if batch_key in done_set:
            print(f"  Batch {batch_num+1}/{n_batches}: cached")
            continue

        to_query = [h for h in batch if h not in results]
        print(f"  Batch {batch_num+1}/{n_batches}: {len(to_query)} queries",
              end="", flush=True)

        if to_query:
            vep_results = _vep_batch(to_query)
            if not vep_results:
                # VEP returned nothing -- do not mark as done, will retry on resume
                print(f" -> VEP returned empty, will retry on resume")
                _save_checkpoint(args.checkpoint, {
                    "done":    list(done_set),
                    "results": results,
                })
                time.sleep(args.delay)
                continue
            n_mapped = 0
            for res in vep_results:
                inp    = res.get("input", "")
                coords = _extract_mane_coords(res, mane_lookup)
                if coords:
                    results[inp] = coords
                    n_mapped += 1
            print(f" -> {n_mapped} MANE Select mapped")
        else:
            print(" -> all cached")

        done_set.add(batch_key)
        _save_checkpoint(args.checkpoint, {
            "done":    list(done_set),
            "results": results,
        })
        time.sleep(args.delay)

    # Build output
    out_rows = []
    for _, row in all_hotspots.iterrows():
        coord = results.get(row["hgvs_input"])
        if not coord:
            continue
        out_rows.append({
            "chrom":       coord["chrom"],
            "pos":         coord["pos"],
            "ref":         coord["ref"],
            "alt":         coord["alt"],
            "gene":        row["gene"],
            "hgvsc":       coord["hgvsc"],
            "hgvsp":       coord["hgvsp"],
            "consequence": coord["consequence"],
            "n_samples":   row["n_samples"],
            "qvalue":      row["qvalue"],
            "version":     row["version"],
        })

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(args.output, sep='\t', index=False)

    print(f"\nOutput: {args.output}")
    print(f"  Mapped to MANE Select GRCh38: {len(out_df)}")
    print(f"  Not mapped (no MANE Select):  {len(all_hotspots) - len(out_df)}")
    print("Done.")


if __name__ == "__main__":
    main()
