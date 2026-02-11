#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana MAG Assembly Pipeline - Docker Run Helper
# ============================================================================
#
# Validates mount paths and builds the docker run command with correct
# volume mounts, user mapping, and Nextflow arguments.
#
# Usage:
#   ./run-docker.sh --input /path/to/reads --outdir /path/to/output [options]
#
# Examples:
#   # Basic assembly + binning
#   ./run-docker.sh --input /data/reads --outdir /data/output
#
#   # With Filtlong pre-filtering
#   ./run-docker.sh --input /data/reads --outdir /data/output \
#       --filtlong_size 40000000000
#
#   # Without MaxBin2
#   ./run-docker.sh --input /data/reads --outdir /data/output \
#       --run_maxbin false
#
# ============================================================================

IMAGE="danaseq-mag"
DOCKER_ARGS=()
NF_ARGS=()
MOUNTS=()

die() { echo "[ERROR] $1" >&2; exit 1; }

usage() {
    echo "Usage: $0 --input DIR --outdir DIR [pipeline options]"
    echo ""
    echo "Required:"
    echo "  --input DIR      Directory containing *.fastq.gz read files"
    echo "  --outdir DIR     Output directory (will be created if needed)"
    echo ""
    echo "All other flags are passed to Nextflow (--dedupe, --filtlong_size, etc.)"
    echo "Run '$0 --help-pipeline' to see Nextflow pipeline options."
    exit 0
}

# ============================================================================
# Parse arguments -- extract paths that need mounting, pass rest to Nextflow
# ============================================================================

INPUT_HOST=""
OUTDIR_HOST=""

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

# Check for FASTQ files
if ! ls "$INPUT_HOST"/*.fastq.gz &>/dev/null; then
    die "No *.fastq.gz files found in $INPUT_HOST"
fi

# Output directory -- create if needed
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
NF_CACHE="${OUTDIR_HOST}/.nextflow-cache"
mkdir -p "${NF_CACHE}/work" "${NF_CACHE}/dotdir" 2>/dev/null || true
MOUNTS+=("-v" "${NF_CACHE}/work:/home/dana/work")
MOUNTS+=("-v" "${NF_CACHE}/dotdir:/home/dana/.nextflow")

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
echo "[INFO] Running: ${DOCKER_CMD[*]}"
echo ""

exec "${DOCKER_CMD[@]}"
