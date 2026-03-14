#!/usr/bin/env python3
# =============================================================================
# ONCOSIEVE — per-database fix script
#
# Author : Dr Christopher Trethewey
# Email  : christopher.trethewey@nhs.net
#
# Phase 2 of 2. Actions the findings from db_audit.py.
# Run db_audit.py first to understand what will be fixed.
#
#   GENIE   — Lifts over coordinates from GRCh37 to GRCh38 using the UCSC
#             chain file.
#
#   TCGA    — Lifts over mc3.v0.2.8.PUBLIC.maf.gz from GRCh37 to GRCh38.
#             Writes mc3.v0.2.8.PUBLIC.GRCh38.maf.gz and a liftover log
#             alongside the input file. Updates config.yaml to point to
#             the lifted file. Rows that fail to map are discarded and logged.
#             HGVSc and HGVSp are transcript-relative and do not change at
#             this step. mane_remap.py updates them to MANE Select post-merge.
#             Writes a lifted MAF alongside the original and updates
#             config.yaml to point to the lifted file.
#
#   COSMIC        — No fix needed (GRCh38, ref 100%).
#   ClinVar       — No fix needed (GRCh38, ref 100%).
#   TP53db        — No fix needed (GRCh38, reader handles column detection).
#   CancerHotspots — No fix needed (GRCh38 MANE remapped).
#
# Usage:
#   python3 tools/db_fix.py --genie-maf data/genie/data_mutations_extended.txt
#   python3 tools/db_fix.py --tcga-maf  data/TCGA/mc3.v0.2.8.PUBLIC.maf.gz
#   python3 tools/db_fix.py --genie-maf <path> --tcga-maf <path>
#
# Optional flags:
#   --chain       Override chain file path (default: data/reference/hg19ToHg38.over.chain.gz)
#   --dry-run     Report what would be done without writing any files
#   --no-backup   Skip creating .orig backup of the original file
# =============================================================================

import argparse
import gzip
import os
import sys
import shutil
from datetime import datetime
from typing import Optional

import yaml

PASS = "\033[92mPASS\033[0m"
FAIL = "\033[91mFAIL\033[0m"
WARN = "\033[93mWARN\033[0m"
INFO = "\033[94mINFO\033[0m"
BOLD = "\033[1m"
RESET = "\033[0m"


def section(title):
    print(f"\n{BOLD}{'=' * 65}{RESET}\n{BOLD}  {title}{RESET}\n{BOLD}{'=' * 65}{RESET}")


def log(tag, msg, detail=""):
    print(f"  {tag}  {msg}")
    if detail:
        for line in detail.splitlines():
            print(f"       {line}")


# ── File utilities ────────────────────────────────────────────────────────────


def _open_r(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _open_w(path: str):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "wt")


def _backup(path: str, no_backup: bool) -> None:
    if no_backup:
        return
    bak = path + ".orig"
    if not os.path.exists(bak):
        shutil.copy2(path, bak)
        log(INFO, f"Backup: {bak}")
    else:
        log(INFO, f"Backup already exists, skipping: {bak}")


# ── Chromosome normalisation ──────────────────────────────────────────────────


def _to_ucsc(c: str) -> str:
    c = c.strip()
    if c.startswith("chr"):
        return c
    return "chrM" if c == "MT" else f"chr{c}"


def _from_ucsc(c: str) -> str:
    if c == "chrM":
        return "MT"
    return c[3:] if c.startswith("chr") else c


# ── GENIE liftover ────────────────────────────────────────────────────────────


def liftover_genie(
    maf_path: str, chain_path: str, dry_run: bool, no_backup: bool
) -> Optional[str]:
    section("Lifting over GENIE from GRCh37 to GRCh38")

    for path, label in [(maf_path, "GENIE MAF"), (chain_path, "Chain file")]:
        if not os.path.isfile(path):
            log(FAIL, f"{label} not found: {path}")
            return None

    try:
        from pyliftover import LiftOver
    except ImportError:
        log(
            FAIL,
            "pyliftover not installed.",
            "Install with: pip install pyliftover --break-system-packages",
        )
        return None

    log(INFO, f"Loading chain file: {chain_path}")
    try:
        lo = LiftOver(chain_path)
    except Exception as e:
        log(FAIL, f"Could not load chain file: {e}")
        return None

    base, ext = maf_path, ""
    for suffix in [".gz", ".txt", ".maf", ".tsv", ".csv"]:
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            ext = suffix + ext
    lifted_path = base + "_grch38_lifted" + ext

    log(INFO, f"Input:  {maf_path}")
    log(INFO, f"Output: {lifted_path}")

    sep = "\t"
    chrom_col = pos_col = end_col = ncbi_col = None

    with _open_r(maf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split(sep)
            col_i = {c: i for i, c in enumerate(cols)}
            chrom_col = col_i.get("Chromosome")
            pos_col   = col_i.get("Start_Position")
            end_col   = col_i.get("End_Position")
            ncbi_col  = col_i.get("NCBI_Build")
            break

    if chrom_col is None or pos_col is None:
        log(FAIL, "Could not find Chromosome or Start_Position in MAF header.")
        return None

    # ── GRCh38 detection: sample first data row ───────────────────────────────
    if ncbi_col is not None:
        with _open_r(maf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split(sep)
                if parts[ncbi_col] == "NCBI_Build":
                    continue  # header row
                if ncbi_col < len(parts):
                    build = parts[ncbi_col].strip().upper()
                    if build in ("GRCH38", "HG38", "38"):
                        log(WARN,
                            "GENIE MAF appears to already be GRCh38 "
                            f"(NCBI_Build={parts[ncbi_col]!r}).",
                            "Skipping liftover. If this is incorrect, check the "
                            "NCBI_Build column or point config.yaml to the original "
                            "pre-lift MAF.")
                        return None
                break

    log(INFO, "Scanning all rows (this may take a few minutes for 3.4M rows) ...")
    n_total = n_ok = n_fail_map = n_fail_chrom = n_parse_err = 0

    with _open_r(maf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(sep)
            if parts[chrom_col] == "Chromosome":
                continue
            if len(parts) <= max(chrom_col, pos_col):
                continue
            try:
                chrom = _to_ucsc(parts[chrom_col])
                pos = int(float(parts[pos_col])) - 1
            except (ValueError, TypeError):
                n_parse_err += 1
                continue

            n_total += 1
            result = lo.convert_coordinate(chrom, pos)
            if not result:
                n_fail_map += 1
                continue
            new_chrom = result[0][0]
            if _from_ucsc(new_chrom) != _from_ucsc(chrom):
                n_fail_chrom += 1
                continue
            n_ok += 1

    pct_ok = n_ok / n_total * 100 if n_total else 0
    pct_fail = (n_fail_map + n_fail_chrom) / n_total * 100 if n_total else 0

    log(INFO, f"Scan complete: {n_total:,} data rows")
    log(INFO, f"  Successfully lifted:          {n_ok:,}  ({pct_ok:.1f}%)")
    log(INFO, f"  Failed to map (discard):      {n_fail_map:,}")
    log(INFO, f"  Chromosome changed (discard): {n_fail_chrom:,}")
    log(INFO, f"  Parse errors (discard):       {n_parse_err:,}")
    log(
        INFO,
        f"  Total discarded:              "
        f"{n_fail_map + n_fail_chrom + n_parse_err:,}  ({pct_fail:.1f}%)",
    )

    if n_ok == 0:
        log(
            FAIL,
            "No rows lifted successfully. "
            "Check chain file is hg19ToHg38.over.chain.gz.",
        )
        return None

    if dry_run:
        log(INFO, f"DRY RUN — would write {n_ok:,} rows to: {lifted_path}")
        return lifted_path

    _backup(maf_path, no_backup)
    log(INFO, "Writing lifted MAF ...")

    n_written = n_skipped = 0
    with _open_r(maf_path) as fh_in, _open_w(lifted_path) as fh_out:
        for line in fh_in:
            if line.startswith("#"):
                fh_out.write(line)
                continue
            parts = line.rstrip("\n").split(sep)
            if parts[chrom_col] == "Chromosome":
                # Write liftover provenance comment before header row
                # NOTE: do NOT update parts[ncbi_col] here — that would rename
                # the column header from "NCBI_Build" to "GRCh38", breaking
                # any downstream tool that looks up the column by name.
                from datetime import datetime as _dt
                fh_out.write(
                    f"#liftover: GRCh37->GRCh38 via pyliftover  "
                    f"date={_dt.now().strftime('%Y-%m-%d')}  "
                    f"chain={chain_path}\n"
                )
                fh_out.write("\t".join(parts) + "\n")
                continue
            if len(parts) <= max(chrom_col, pos_col):
                fh_out.write(line)
                continue

            try:
                chrom = _to_ucsc(parts[chrom_col])
                orig_pos = int(float(parts[pos_col]))
                pos_0base = orig_pos - 1
            except (ValueError, TypeError):
                n_skipped += 1
                continue

            result = lo.convert_coordinate(chrom, pos_0base)
            if not result:
                n_skipped += 1
                continue

            new_chrom_ucsc, new_pos_0base, _, _ = result[0]
            if _from_ucsc(new_chrom_ucsc) != _from_ucsc(chrom):
                n_skipped += 1
                continue

            new_chrom = _from_ucsc(new_chrom_ucsc)
            new_pos = new_pos_0base + 1

            parts[chrom_col] = new_chrom
            parts[pos_col] = str(new_pos)

            # Update NCBI_Build value per data row
            if ncbi_col is not None and ncbi_col < len(parts):
                parts[ncbi_col] = "GRCh38"

            if end_col is not None and end_col < len(parts):
                try:
                    orig_end = int(float(parts[end_col]))
                    offset = orig_end - orig_pos
                    parts[end_col] = str(new_pos + offset)
                except (ValueError, TypeError):
                    pass

            fh_out.write("\t".join(parts) + "\n")
            n_written += 1

    log(
        PASS,
        f"GENIE lifted MAF written: {lifted_path}",
        f"{n_written:,} rows written, {n_skipped:,} rows discarded",
    )

    # Write liftover log
    log_path = lifted_path.replace(".txt", ".liftover.log").replace(".gz", ".liftover.log")
    if not log_path.endswith(".liftover.log"):
        log_path = lifted_path + ".liftover.log"
    _write_liftover_log(
        log_path     = log_path,
        source       = "AACR GENIE v19.0",
        chain_path   = chain_path,
        input_path   = maf_path,
        output_path  = lifted_path,
        n_total      = n_total,
        n_written    = n_written,
        n_skipped    = n_skipped,
        n_fail_map   = n_fail_map,
        n_fail_chrom = n_fail_chrom,
        n_parse_err  = n_parse_err,
    )
    log(PASS, f"Liftover log written: {log_path}")

    return lifted_path


# ── TCGA liftover ─────────────────────────────────────────────────────────────


def _write_liftover_log(log_path: str, source: str, chain_path: str,
                         input_path: str, output_path: str,
                         n_total: int, n_written: int, n_skipped: int,
                         n_fail_map: int, n_fail_chrom: int,
                         n_parse_err: int) -> None:
    """Write a liftover summary log alongside the lifted file."""
    from datetime import datetime
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    pct_ok   = n_written / n_total * 100 if n_total else 0
    pct_fail = (n_total - n_written) / n_total * 100 if n_total else 0
    lines = [
        f"ONCOSIEVE liftover log",
        f"======================",
        f"Source        : {source}",
        f"Date          : {ts}",
        f"Chain file    : {chain_path}",
        f"Input         : {input_path}",
        f"Output        : {output_path}",
        f"",
        f"Row counts",
        f"----------",
        f"Total input rows    : {n_total:,}",
        f"Successfully lifted : {n_written:,}  ({pct_ok:.2f}%)",
        f"Failed to map       : {n_fail_map:,}",
        f"Chromosome changed  : {n_fail_chrom:,}",
        f"Parse errors        : {n_parse_err:,}",
        f"Total discarded     : {n_total - n_written:,}  ({pct_fail:.2f}%)",
        f"",
        f"Method: pyliftover coordinate conversion (hg19->hg38)",
        f"Rows failing liftover are discarded.",
        f"HGVSc and HGVSp are transcript-relative and unchanged by liftover.",
    ]
    with open(log_path, 'w') as fh:
        fh.write("\n".join(lines) + "\n")


def liftover_tcga(
    maf_path: str, chain_path: str, dry_run: bool, no_backup: bool
) -> Optional[str]:
    section("Lifting over TCGA mc3 MAF from GRCh37 to GRCh38")

    for path, label in [(maf_path, "TCGA MAF"), (chain_path, "Chain file")]:
        if not os.path.isfile(path):
            log(FAIL, f"{label} not found: {path}")
            return None

    try:
        from pyliftover import LiftOver
    except ImportError:
        log(
            FAIL,
            "pyliftover not installed.",
            "Install with: pip install pyliftover --break-system-packages",
        )
        return None

    log(INFO, f"Loading chain file: {chain_path}")
    try:
        lo = LiftOver(chain_path)
    except Exception as e:
        log(FAIL, f"Could not load chain file: {e}")
        return None

    # Output path: same directory, GRCh38 suffix
    out_dir   = os.path.dirname(maf_path)
    lifted_path = os.path.join(out_dir, "mc3.v0.2.8.PUBLIC.GRCh38.maf.gz")
    log_path    = os.path.join(out_dir, "mc3.v0.2.8.PUBLIC.GRCh38.liftover.log")

    log(INFO, f"Input:  {maf_path}")
    log(INFO, f"Output: {lifted_path}")
    log(INFO, f"Log:    {log_path}")

    sep = "\t"
    chrom_col = pos_col = end_col = ncbi_col = None

    with _open_r(maf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split(sep)
            col_i = {c: i for i, c in enumerate(cols)}
            chrom_col = col_i.get("Chromosome")
            pos_col   = col_i.get("Start_Position")
            end_col   = col_i.get("End_Position")
            ncbi_col  = col_i.get("NCBI_Build")
            break

    if chrom_col is None or pos_col is None:
        log(FAIL, "Could not find Chromosome or Start_Position in MAF header.")
        return None

    # ── GRCh38 detection: sample first data row ───────────────────────────────
    if ncbi_col is not None:
        with _open_r(maf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split(sep)
                if parts[ncbi_col] == "NCBI_Build":
                    continue  # header row
                if ncbi_col < len(parts):
                    build = parts[ncbi_col].strip().upper()
                    if build in ("GRCH38", "HG38", "38"):
                        log(WARN,
                            "TCGA MAF appears to already be GRCh38 "
                            f"(NCBI_Build={parts[ncbi_col]!r}).",
                            "Skipping liftover. If this is incorrect, check the "
                            "NCBI_Build column or point config.yaml to the original "
                            "pre-lift MAF.")
                        return None
                break

    log(INFO, "Scanning all rows (this will take several minutes for 3.5M rows)...")
    n_total = n_ok = n_fail_map = n_fail_chrom = n_parse_err = 0

    with _open_r(maf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(sep)
            if parts[chrom_col] == "Chromosome":
                continue
            if len(parts) <= max(chrom_col, pos_col):
                continue
            try:
                chrom = _to_ucsc(parts[chrom_col])
                pos   = int(float(parts[pos_col])) - 1
            except (ValueError, TypeError):
                n_parse_err += 1
                continue

            n_total += 1
            result = lo.convert_coordinate(chrom, pos)
            if not result:
                n_fail_map += 1
                continue
            new_chrom = result[0][0]
            if _from_ucsc(new_chrom) != _from_ucsc(chrom):
                n_fail_chrom += 1
                continue
            n_ok += 1

    pct_ok   = n_ok / n_total * 100 if n_total else 0
    pct_fail = (n_fail_map + n_fail_chrom) / n_total * 100 if n_total else 0

    log(INFO, f"Scan complete: {n_total:,} data rows")
    log(INFO, f"  Successfully lifted:          {n_ok:,}  ({pct_ok:.1f}%)")
    log(INFO, f"  Failed to map (discard):      {n_fail_map:,}")
    log(INFO, f"  Chromosome changed (discard): {n_fail_chrom:,}")
    log(INFO, f"  Parse errors (discard):       {n_parse_err:,}")
    log(INFO, f"  Total discarded:              "
        f"{n_fail_map + n_fail_chrom + n_parse_err:,}  ({pct_fail:.1f}%)")

    if n_ok == 0:
        log(FAIL, "No rows lifted. Check chain file is hg19ToHg38.over.chain.gz.")
        return None

    if dry_run:
        log(INFO, f"DRY RUN — would write {n_ok:,} rows to: {lifted_path}")
        return lifted_path

    _backup(maf_path, no_backup)
    log(INFO, "Writing lifted MAF (gzip) ...")

    n_written = n_skipped = 0
    with _open_r(maf_path) as fh_in, _open_w(lifted_path) as fh_out:
        for line in fh_in:
            if line.startswith("#"):
                fh_out.write(line)
                continue
            parts = line.rstrip("\n").split(sep)
            if parts[chrom_col] == "Chromosome":
                # Write liftover provenance comment before header row
                # NOTE: do NOT update parts[ncbi_col] here — that would rename
                # the column header from "NCBI_Build" to "GRCh38", breaking
                # any downstream tool that looks up the column by name.
                from datetime import datetime as _dt
                fh_out.write(
                    f"#liftover: GRCh37->GRCh38 via pyliftover  "
                    f"date={_dt.now().strftime('%Y-%m-%d')}  "
                    f"chain={chain_path}\n"
                )
                fh_out.write("\t".join(parts) + "\n")
                continue
            if len(parts) <= max(chrom_col, pos_col):
                fh_out.write(line)
                continue

            try:
                chrom    = _to_ucsc(parts[chrom_col])
                orig_pos = int(float(parts[pos_col]))
                pos_0    = orig_pos - 1
            except (ValueError, TypeError):
                n_skipped += 1
                continue

            result = lo.convert_coordinate(chrom, pos_0)
            if not result:
                n_skipped += 1
                continue

            new_chrom_ucsc, new_pos_0, _, _ = result[0]
            if _from_ucsc(new_chrom_ucsc) != _from_ucsc(chrom):
                n_skipped += 1
                continue

            new_chrom = _from_ucsc(new_chrom_ucsc)
            new_pos   = new_pos_0 + 1

            parts[chrom_col] = new_chrom
            parts[pos_col]   = str(new_pos)

            # Update NCBI_Build per data row
            if ncbi_col is not None and ncbi_col < len(parts):
                parts[ncbi_col] = "GRCh38"

            if end_col is not None and end_col < len(parts):
                try:
                    orig_end       = int(float(parts[end_col]))
                    offset         = orig_end - orig_pos
                    parts[end_col] = str(new_pos + offset)
                except (ValueError, TypeError):
                    pass

            fh_out.write("\t".join(parts) + "\n")
            n_written += 1

    log(
        PASS,
        f"TCGA lifted MAF written: {lifted_path}",
        f"{n_written:,} rows written, {n_skipped:,} rows discarded",
    )

    # Write liftover log
    _write_liftover_log(
        log_path    = log_path,
        source      = "TCGA mc3 v0.2.8",
        chain_path  = chain_path,
        input_path  = maf_path,
        output_path = lifted_path,
        n_total     = n_total,
        n_written   = n_written,
        n_skipped   = n_skipped,
        n_fail_map  = n_fail_map,
        n_fail_chrom= n_fail_chrom,
        n_parse_err = n_parse_err,
    )
    log(PASS, f"Liftover log written: {log_path}")

    return lifted_path


# ── Config update ─────────────────────────────────────────────────────────────


def report_no_fix_needed(name: str, reason: str) -> None:
    section(f"{name} — no source-level fix required")
    log(INFO, reason)


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="ONCOSIEVE db_fix.py — Lift MAF files from GRCh37 to GRCh38.\n"
        "Completely independent of config.yaml. Takes MAF paths as direct arguments.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Lift both GENIE and TCGA:\n"
            "  python3 tools/db_fix.py \\\n"
            "      --genie-maf data/genie/data_mutations_extended.txt \\\n"
            "      --tcga-maf  data/TCGA/mc3.v0.2.8.PUBLIC.maf.gz\n\n"
            "  # Lift GENIE only:\n"
            "  python3 tools/db_fix.py --genie-maf data/genie/data_mutations_extended.txt\n\n"
            "  # Dry run to check what would happen:\n"
            "  python3 tools/db_fix.py --genie-maf data/genie/data_mutations_extended.txt --dry-run\n"
        )
    )
    parser.add_argument("--genie-maf",  default="",
                        help="Path to GENIE MAF (GRCh37). Output written alongside as *_grch38_lifted.txt")
    parser.add_argument("--tcga-maf",   default="",
                        help="Path to TCGA mc3 MAF (GRCh37). Output written alongside as *.GRCh38.maf.gz")
    parser.add_argument("--chain",      default="data/reference/hg19ToHg38.over.chain.gz",
                        help="hg19ToHg38 chain file (default: data/reference/hg19ToHg38.over.chain.gz)")
    parser.add_argument("--dry-run",    action="store_true",
                        help="Report what would be done without writing any files")
    parser.add_argument("--no-backup",  action="store_true",
                        help="Skip creating .orig backup of the original file")
    args = parser.parse_args()

    if not args.genie_maf and not args.tcga_maf:
        parser.error("Specify at least one of --genie-maf or --tcga-maf")

    if args.dry_run:
        print(f"\n{WARN}  DRY RUN MODE — no files will be written\n")

    chain_path = args.chain

    if args.genie_maf:
        liftover_genie(args.genie_maf, chain_path, args.dry_run, args.no_backup)

    if args.tcga_maf:
        liftover_tcga(args.tcga_maf, chain_path, args.dry_run, args.no_backup)

    section("Next steps")
    log(
        WARN,
        "Verify config.yaml points to the lifted files before running the pipeline:",
        "  genie.maf : data/genie/data_mutations_extended_grch38_lifted.txt\n"
        "  tcga.maf  : data/TCGA/mc3.v0.2.8.PUBLIC.GRCh38.maf.gz\n\n"
        "Then run: bash run_pipeline.sh",
    )

    section("Summary")
    ts = datetime.now().strftime("%Y-%m-%d %H:%M")
    if args.dry_run:
        log(INFO, f"Dry run completed at {ts}. No files written.")
    else:
        log(PASS, f"Liftover completed at {ts}")

    print(f"\n{'=' * 65}\n")


if __name__ == "__main__":
    main()
