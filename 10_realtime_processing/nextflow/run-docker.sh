#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana Pipeline - Docker Run Helper
# ============================================================================
#
# Validates mount paths and builds the docker run command with correct
# volume mounts, user mapping, and Nextflow arguments.
#
# Usage:
#   ./run-docker.sh --input /path/to/data --outdir /path/to/output [options]
#
# Examples:
#   # Basic QC
#   ./run-docker.sh --input /data/run1 --outdir /data/output
#
#   # Full pipeline with Kraken2
#   ./run-docker.sh --input /data/run1 --outdir /data/output \
#       --run_kraken --kraken_db /path/to/krakendb \
#       --run_prokka --run_sketch --run_tetra
#
#   # With HMM profiling
#   ./run-docker.sh --input /data/run1 --outdir /data/output \
#       --run_prokka --hmm_databases /path/to/CANT-HYD.hmm
#
#   # Multiple HMM databases
#   ./run-docker.sh --input /data/run1 --outdir /data/output \
#       --run_prokka --hmm_databases /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm
#
#   # With DuckDB integration (load results into database)
#   ./run-docker.sh --input /data/run1 --outdir /data/output \
#       --run_kraken --kraken_db /path/to/krakendb \
#       --run_prokka --run_sketch --run_tetra \
#       --run_db_integration --danadir /path/to/danaSeq/10_realtime_processing
#
# ============================================================================

IMAGE="danaseq-realtime"
DOCKER_ARGS=()
NF_ARGS=()

# Paths to mount (host:container pairs)
MOUNTS=()

die() { echo "[ERROR] $1" >&2; exit 1; }
warn() { echo "[WARNING] $1" >&2; }

usage() {
    echo "Usage: $0 --input DIR --outdir DIR [pipeline options]"
    echo ""
    echo "Required:"
    echo "  --input DIR          Nanopore run directory (contains fastq_pass/)"
    echo "  --outdir DIR         Output directory (will be created if needed)"
    echo ""
    echo "Optional mounts:"
    echo "  --kraken_db DIR      Kraken2 database directory"
    echo "  --hmm_databases LIST Comma-separated HMM file paths"
    echo "  --danadir DIR        R scripts directory for DuckDB integration"
    echo "                       (10_realtime_processing/ in the repo)"
    echo ""
    echo "All other flags are passed to Nextflow (--run_prokka, --run_sketch, etc.)"
    echo "Run '$0 --help-pipeline' to see Nextflow pipeline options."
    exit 0
}

# ============================================================================
# Parse arguments — extract paths that need mounting, pass rest to Nextflow
# ============================================================================

INPUT_HOST=""
OUTDIR_HOST=""
KRAKEN_DB_HOST=""
DANADIR_HOST=""
HMM_HOST=""

while (( $# )); do
    case "$1" in
        -h|--help)
            usage ;;
        --help-pipeline)
            docker run --user "$(id -u):$(id -g)" "$IMAGE" \
                run /pipeline/main.nf --help
            exit 0 ;;
        --input)
            [[ -z "${2:-}" ]] && die "--input requires a directory path"
            INPUT_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            shift 2 ;;
        --outdir)
            [[ -z "${2:-}" ]] && die "--outdir requires a directory path"
            OUTDIR_HOST="$(realpath -m "$2")"
            shift 2 ;;
        --kraken_db)
            [[ -z "${2:-}" ]] && die "--kraken_db requires a directory path"
            KRAKEN_DB_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            shift 2 ;;
        --danadir)
            [[ -z "${2:-}" ]] && die "--danadir requires a directory path"
            DANADIR_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            shift 2 ;;
        --hmm_databases)
            [[ -z "${2:-}" ]] && die "--hmm_databases requires file path(s)"
            HMM_HOST="$2"
            shift 2 ;;
        *)
            NF_ARGS+=("$1")
            shift ;;
    esac
done

# ============================================================================
# Validate required paths
# ============================================================================

[[ -z "$INPUT_HOST" ]] && die "--input is required. Run with --help for usage."
[[ -z "$OUTDIR_HOST" ]] && die "--outdir is required. Run with --help for usage."

# Input directory must exist
[[ -d "$INPUT_HOST" ]] || die "Input directory does not exist: $INPUT_HOST"

# Output directory — create if needed
if [[ ! -d "$OUTDIR_HOST" ]]; then
    echo "[INFO] Creating output directory: $OUTDIR_HOST"
    mkdir -p "$OUTDIR_HOST" || die "Cannot create output directory: $OUTDIR_HOST"
fi
[[ -w "$OUTDIR_HOST" ]] || die "Output directory is not writable: $OUTDIR_HOST"

# Mount input (read-only) and output
MOUNTS+=("-v" "${INPUT_HOST}:/data/input:ro")
MOUNTS+=("-v" "${OUTDIR_HOST}:/data/output")
NF_ARGS=("--input" "/data/input" "--outdir" "/data/output" "${NF_ARGS[@]}")

# Persistent Nextflow directories for -resume support
# work/: task cache (inputs, outputs, scripts)
# .nextflow/: session history, run metadata
# Both must survive container restarts for -resume to work
NF_CACHE="${OUTDIR_HOST}/.nextflow-cache"
mkdir -p "${NF_CACHE}/work" "${NF_CACHE}/dotdir" 2>/dev/null || true
MOUNTS+=("-v" "${NF_CACHE}/work:/home/dana/work")
MOUNTS+=("-v" "${NF_CACHE}/dotdir:/home/dana/.nextflow")

# ============================================================================
# Optional: Kraken2 database
# ============================================================================

if [[ -n "$KRAKEN_DB_HOST" ]]; then
    [[ -d "$KRAKEN_DB_HOST" ]] || die "Kraken2 database directory does not exist: $KRAKEN_DB_HOST"
    # Check for .k2d files
    if ! ls "$KRAKEN_DB_HOST"/*.k2d &>/dev/null; then
        warn "No .k2d files found in $KRAKEN_DB_HOST — are you sure this is a Kraken2 database?"
    fi
    MOUNTS+=("-v" "${KRAKEN_DB_HOST}:/kraken_db:ro")
    NF_ARGS+=("--kraken_db" "/kraken_db")
fi

# ============================================================================
# Optional: DuckDB integration (R scripts directory)
# ============================================================================

if [[ -n "$DANADIR_HOST" ]]; then
    [[ -d "$DANADIR_HOST" ]] || die "R scripts directory does not exist: $DANADIR_HOST"
    # Sanity check: look for the R scripts
    if ! ls "$DANADIR_HOST"/4*_db.r &>/dev/null; then
        warn "No 4*_db.r scripts found in $DANADIR_HOST — expected 10_realtime_processing/ directory"
    fi
    MOUNTS+=("-v" "${DANADIR_HOST}:/danadir:ro")
    NF_ARGS+=("--danadir" "/danadir")
fi

# ============================================================================
# Optional: HMM databases
# ============================================================================

if [[ -n "$HMM_HOST" ]]; then
    container_hmm_paths=()
    hmm_counter=0
    IFS=',' read -ra HMM_FILES <<< "$HMM_HOST"
    for hmm_file in "${HMM_FILES[@]}"; do
        hmm_file="$(echo "$hmm_file" | xargs)"  # trim whitespace
        hmm_resolved="$(realpath "$hmm_file" 2>/dev/null || echo "$hmm_file")"
        [[ -f "$hmm_resolved" ]] || die "HMM file does not exist: $hmm_file"
        hmm_name="$(basename "$hmm_resolved")"
        container_path="/hmm/${hmm_name}"
        MOUNTS+=("-v" "${hmm_resolved}:${container_path}:ro")
        container_hmm_paths+=("$container_path")
        hmm_counter=$((hmm_counter + 1))
    done
    # Join with commas
    hmm_joined="$(IFS=','; echo "${container_hmm_paths[*]}")"
    NF_ARGS+=("--hmm_databases" "$hmm_joined")
    echo "[INFO] Mounting $hmm_counter HMM database(s)"
fi

# ============================================================================
# Build and run Docker command
# ============================================================================

DOCKER_CMD=(
    docker run
    --user "$(id -u):$(id -g)"
    "${MOUNTS[@]}"
    "$IMAGE"
    run /pipeline/main.nf
    "${NF_ARGS[@]}"
    -resume
)

echo "[INFO] Input:  $INPUT_HOST"
echo "[INFO] Output: $OUTDIR_HOST"
[[ -n "$KRAKEN_DB_HOST" ]] && echo "[INFO] Kraken DB: $KRAKEN_DB_HOST"
[[ -n "$DANADIR_HOST" ]] && echo "[INFO] R scripts: $DANADIR_HOST"
echo "[INFO] Running: ${DOCKER_CMD[*]}"
echo ""

exec "${DOCKER_CMD[@]}"
