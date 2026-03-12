#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh
#
# Orchestrates the full pan-cancer whitelist pipeline.
#
# Steps:
#   1. Check dependencies
#   2. Install Python requirements
#   3. Build the whitelist (all parsers -> merge -> filter -> tier -> output)
#   4. Normalise output VCF with bcftools norm
#   5. Liftover DoCM hg19 variants (if bcftools >= 1.17 and chain file present)
#   6. Index output VCF with tabix
#   7. Print summary statistics
#   8. (Optional) Run Mutect2 post-filter rescue
#
# Usage:
#   bash run_pipeline.sh [data_dir] [options]
#
# Arguments:
#   data_dir               Optional path to reference data directory.
#                          Overrides relative data paths in config.yaml.
#                          e.g. bash run_pipeline.sh /srv/data/reference/
#
# Options:
#   --skip-sources STR     Comma-separated sources to skip, e.g. "genie,cbioportal"
#   --intermediate-only    Run parsers and save intermediates; skip merge
#   --rescue               Run Mutect2 post-filter rescue after building whitelist
#   --mutect2-vcf PATH     Mutect2 FilterMutectCalls VCF (required if --rescue)
#   --tumor-sample NAME    Tumour sample name in the Mutect2 VCF
#   --rescue-output PATH   Output path for rescued VCF
#   --help                 Show this message and exit
# =============================================================================

set -euo pipefail

# Resolve DATA_DIR to absolute path BEFORE cd changes working directory
DATA_DIR=""
if [[ $# -gt 0 && ! "$1" =~ ^-- ]]; then
    DATA_DIR="$(realpath "$1")"
    shift
fi

# Always run relative to the script's own directory
cd "$(dirname "$(realpath "$0")")"

# =============================================================================
# DATA DIRECTORY — optional first positional argument
# =============================================================================

if [[ -n "$DATA_DIR" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Data directory: $DATA_DIR"
fi

# =============================================================================
# USER CONFIGURATION — set paths here before running
# =============================================================================

PYTHON="python3"
CONFIG="config.yaml"

OUTPUT_DIR="output"
INTERMEDIATE_DIR="intermediate"
PREFIX="pan_cancer_whitelist_GRCh38"

REFERENCE_FASTA="data/reference/GRCh38.fa"
CHAIN_HG19_HG38="data/reference/hg19ToHg38.over.chain.gz"
HG19_FASTA="data/reference/hg19.fa"

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

SKIP_SOURCES=""
INTERMEDIATE_ONLY=false
DO_RESCUE=false
MUTECT2_VCF=""
TUMOR_SAMPLE=""
RESCUE_OUTPUT=""

usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,3\}//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-sources)      SKIP_SOURCES="$2";       shift 2 ;;
        --intermediate-only) INTERMEDIATE_ONLY=true;  shift   ;;
        --rescue)            DO_RESCUE=true;           shift   ;;
        --mutect2-vcf)       MUTECT2_VCF="$2";        shift 2 ;;
        --tumor-sample)      TUMOR_SAMPLE="$2";        shift 2 ;;
        --rescue-output)     RESCUE_OUTPUT="$2";       shift 2 ;;
        --help|-h)           usage ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

# =============================================================================
# HELPERS
# =============================================================================

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die()  { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }

# =============================================================================
# STEP 1: PRE-RUN AUDIT
# =============================================================================

log "Running pre-run audit..."
"$PYTHON" pre_check.py --config "$CONFIG" ${DATA_DIR:+--data-dir "$DATA_DIR"} || die "Pre-run audit failed. Fix errors above before running pipeline."

# =============================================================================
# STEP 1b: CHECK DEPENDENCIES
# =============================================================================

log "Checking dependencies..."

"$PYTHON" --version  || die "python3 not found. Install Python >= 3.10."
tabix --version 2>&1 | head -1 || die "tabix not found. Install htslib."
bcftools --version | head -1   || die "bcftools not found."

BCFTOOLS_VERSION=$(bcftools --version | head -1 | grep -oP '\d+\.\d+' | head -1)
BCFTOOLS_MAJOR=$(echo "$BCFTOOLS_VERSION" | cut -d. -f1)
BCFTOOLS_MINOR=$(echo "$BCFTOOLS_VERSION" | cut -d. -f2)

HAVE_LIFTOVER=false
if [[ "$BCFTOOLS_MAJOR" -gt 1 ]] || \
   [[ "$BCFTOOLS_MAJOR" -eq 1 && "$BCFTOOLS_MINOR" -ge 17 ]]; then
    HAVE_LIFTOVER=true
fi

log "  Python   : $("$PYTHON" --version)"
log "  bcftools : $BCFTOOLS_VERSION  (liftover support: $HAVE_LIFTOVER)"
log "  tabix    : $(tabix --version 2>&1 | head -1)"

# =============================================================================
# STEP 2: INSTALL PYTHON REQUIREMENTS
# =============================================================================

log "Installing Python requirements..."
"$PYTHON" -m pip install --quiet -r requirements.txt \
    || die "pip install failed. Check requirements.txt."
log "  Done."

# =============================================================================
# STEP 3: BUILD WHITELIST
# =============================================================================

WL_RAW="${OUTPUT_DIR}/${PREFIX}.vcf.gz"

if [[ -f "$WL_RAW" ]]; then
    log "Whitelist VCF already exists at $WL_RAW — skipping build step."
    log "  Delete $WL_RAW to force a full rebuild."
else
    log "Building pan-cancer whitelist..."

    BUILD_ARGS="--config config.yaml${DATA_DIR:+ --data-dir $DATA_DIR}"
    [[ -n "$SKIP_SOURCES" ]] && BUILD_ARGS="$BUILD_ARGS --skip-sources $SKIP_SOURCES"
    $INTERMEDIATE_ONLY       && BUILD_ARGS="$BUILD_ARGS --intermediate-only"

    "$PYTHON" build_whitelist.py $BUILD_ARGS \
        || die "build_whitelist.py failed."
fi

if $INTERMEDIATE_ONLY; then
    log "Intermediate-only mode: stopping before merge."
    log "Intermediate files written to: $INTERMEDIATE_DIR/"
    exit 0
fi

WL_RAW="${OUTPUT_DIR}/${PREFIX}.vcf.gz"
[[ -f "$WL_RAW" ]] || die "Expected output VCF not found: $WL_RAW"

# =============================================================================
# STEP 4: INDEX VCF
# =============================================================================

WL_FINAL="$WL_RAW"
log "Indexing whitelist VCF..."
tabix -p vcf "$WL_FINAL" || die "tabix indexing failed."
log "  Indexed: ${WL_FINAL}.tbi"

# =============================================================================
# STEP 5: LIFTOVER DoCM (hg19 -> hg38)
# =============================================================================

DOCM_INTER="${INTERMEDIATE_DIR}/docm_hg19.tsv.gz"
DOCM_VCF_HG19="${INTERMEDIATE_DIR}/docm_hg19.vcf"
DOCM_VCF_HG38="${INTERMEDIATE_DIR}/docm_hg38.vcf"

if [[ -f "$DOCM_INTER" ]]; then
    if $HAVE_LIFTOVER && [[ -f "$CHAIN_HG19_HG38" && -f "$HG19_FASTA" ]]; then
        log "Lifting over DoCM variants (hg19 -> hg38)..."

        "$PYTHON" - <<'PYEOF'
import pandas as pd

df = pd.read_csv('intermediate/docm_hg19.tsv.gz', sep='\t', dtype=str)
df = df[
    df['ref'].notna() & df['ref'].ne('') &
    df['alt'].notna() & df['alt'].ne('') &
    df['pos'].notna()
]
with open('intermediate/docm_hg19.vcf', 'w') as fh:
    fh.write('##fileformat=VCFv4.2\n')
    fh.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for _, row in df.iterrows():
        chrom = str(row['chrom']).replace('chr', '')
        fh.write(f"{chrom}\t{row['pos']}\t.\t{row['ref']}\t{row['alt']}\t.\t.\t.\n")
print(f"DoCM: wrote {len(df)} rows to intermediate/docm_hg19.vcf")
PYEOF

        bcftools liftover \
            --src-fasta  "$HG19_FASTA" \
            --fasta-ref  "$REFERENCE_FASTA" \
            --chain      "$CHAIN_HG19_HG38" \
            --output     "$DOCM_VCF_HG38" \
            "$DOCM_VCF_HG19" 2>/dev/null \
            && log "  DoCM liftover complete: $DOCM_VCF_HG38" \
            || warn "DoCM liftover failed. DoCM variants excluded from hg38 whitelist."
    else
        [[ ! -f "$CHAIN_HG19_HG38" ]] && \
            warn "DoCM liftover skipped: chain file not found at $CHAIN_HG19_HG38."
        [[ ! -f "$HG19_FASTA" ]] && \
            warn "DoCM liftover skipped: hg19 FASTA not found at $HG19_FASTA."
        ! $HAVE_LIFTOVER && \
            warn "DoCM liftover skipped: bcftools >= 1.17 required (found $BCFTOOLS_VERSION)."
        warn "DoCM variants remain as hg19 coordinates. Exclude from analysis if liftover is not performed."
    fi
fi

# =============================================================================
# STEP 6: INDEX WITH tabix
# =============================================================================

log "Indexing whitelist VCF..."
tabix -f -p vcf "$WL_FINAL" \
    || die "tabix indexing failed for $WL_FINAL."
log "  Index written: ${WL_FINAL}.tbi"

# =============================================================================
# STEP 7: SUMMARY STATISTICS
# =============================================================================

log "Summary statistics:"
N_TOTAL=$(bcftools view -H "$WL_FINAL" | wc -l | tr -d ' ')
N_TIER1=$(bcftools view -H -i 'INFO/WL_TIER=1' "$WL_FINAL" | wc -l | tr -d ' ')
N_TIER2=$(bcftools view -H -i 'INFO/WL_TIER=2' "$WL_FINAL" | wc -l | tr -d ' ')
N_TIER3=$(bcftools view -H -i 'INFO/WL_TIER=3' "$WL_FINAL" | wc -l | tr -d ' ')

log "  Total whitelist variants : $N_TOTAL"
log "  Tier 1 (highest conf)    : $N_TIER1"
log "  Tier 2                   : $N_TIER2"
log "  Tier 3 (min threshold)   : $N_TIER3"
log ""
log "  TSV : ${OUTPUT_DIR}/${PREFIX}.tsv.gz"
log "  VCF : ${WL_FINAL}"

# =============================================================================
# STEP 8: MUTECT2 POST-FILTER RESCUE (optional)
# =============================================================================

if $DO_RESCUE; then
    [[ -n "$MUTECT2_VCF" ]] || die "--rescue requires --mutect2-vcf PATH"
    [[ -f "$MUTECT2_VCF"  ]] || die "Mutect2 VCF not found: $MUTECT2_VCF"

    if [[ -z "$RESCUE_OUTPUT" ]]; then
        RESCUE_OUTPUT="${MUTECT2_VCF%.vcf*}.rescued.vcf.gz"
    fi

    TUMOR_ARG=""
    [[ -n "$TUMOR_SAMPLE" ]] && TUMOR_ARG="--tumor-sample $TUMOR_SAMPLE"

    log "Running Mutect2 post-filter rescue..."
    log "  Input  : $MUTECT2_VCF"
    log "  Output : $RESCUE_OUTPUT"

    "$PYTHON" mutect2_rescue.py \
        --input     "$MUTECT2_VCF" \
        --whitelist "$WL_FINAL" \
        --output    "$RESCUE_OUTPUT" \
        --config    config.yaml \
        $TUMOR_ARG \
        || die "mutect2_rescue.py failed."

    tabix -f -p vcf "$RESCUE_OUTPUT" \
        || warn "tabix indexing of rescued VCF failed."

    log "Rescue complete: $RESCUE_OUTPUT"
fi

log "Pipeline complete."
