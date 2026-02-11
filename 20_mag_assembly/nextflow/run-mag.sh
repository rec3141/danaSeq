#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana MAG Assembly Pipeline Launcher
# ============================================================================
#
# Runs the Nextflow pipeline locally (default) or via Docker (--docker).
#
# Usage:
#   ./run-mag.sh --input /path/to/reads --outdir /path/to/output [options]
#
# Examples:
#   # Basic assembly + binning (local conda)
#   ./run-mag.sh --input /data/reads --outdir /data/output
#
#   # With Filtlong pre-filtering
#   ./run-mag.sh --input /data/reads --outdir /data/output \
#       --filtlong_size 40000000000
#
#   # Without MaxBin2
#   ./run-mag.sh --input /data/reads --outdir /data/output \
#       --run_maxbin false
#
#   # Docker mode
#   ./run-mag.sh --docker --input /data/reads --outdir /data/output
#
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="danaseq-mag"
USE_DOCKER=false
NF_ARGS=()
MOUNTS=()

die() { echo "[ERROR] $1" >&2; exit 1; }

usage() {
    echo "Usage: $0 --input DIR --outdir DIR [pipeline options]"
    echo ""
    echo "Required:"
    echo "  --input DIR      Directory containing reads (*.fastq.gz or barcode structure)"
    echo "  --outdir DIR     Output directory (will be created if needed)"
    echo ""
    echo "Mode:"
    echo "  --docker         Run in Docker instead of local conda"
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
            if [[ "$USE_DOCKER" == true ]]; then
                docker run --user "$(id -u):$(id -g)" "$IMAGE" \
                    run /pipeline/main.nf --help
            else
                mamba run -p "${SCRIPT_DIR}/conda-envs/dana-mag-flye" \
                    nextflow run "${SCRIPT_DIR}/main.nf" --help
            fi
            exit 0 ;;
        --docker)
            USE_DOCKER=true
            shift ;;
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

# Output directory -- create if needed
if [[ ! -d "$OUTDIR_HOST" ]]; then
    echo "[INFO] Creating output directory: $OUTDIR_HOST"
    mkdir -p "$OUTDIR_HOST" || die "Cannot create output directory: $OUTDIR_HOST"
fi
[[ -w "$OUTDIR_HOST" ]] || die "Output directory is not writable: $OUTDIR_HOST"

# ============================================================================
# Docker mode
# ============================================================================

if [[ "$USE_DOCKER" == true ]]; then
    # Mount input (read-only) and output
    MOUNTS+=("-v" "${INPUT_HOST}:/data/input:ro")
    MOUNTS+=("-v" "${OUTDIR_HOST}:/data/output")
    NF_ARGS=("--input" "/data/input" "--outdir" "/data/output" "${NF_ARGS[@]}")

    # Persistent Nextflow directories for -resume support
    NF_CACHE="${OUTDIR_HOST}/.nextflow-cache"
    mkdir -p "${NF_CACHE}/work" "${NF_CACHE}/dotdir" 2>/dev/null || true
    MOUNTS+=("-v" "${NF_CACHE}/work:/home/dana/work")
    MOUNTS+=("-v" "${NF_CACHE}/dotdir:/home/dana/.nextflow")

    DOCKER_CMD=(
        docker run
        --user "$(id -u):$(id -g)"
        "${MOUNTS[@]}"
        "$IMAGE"
        run /pipeline/main.nf
        "${NF_ARGS[@]}"
        -resume
    )

    echo "[INFO] Mode:   Docker"
    echo "[INFO] Input:  $INPUT_HOST"
    echo "[INFO] Output: $OUTDIR_HOST"
    echo "[INFO] Running: ${DOCKER_CMD[*]}"
    echo ""

    exec "${DOCKER_CMD[@]}"
fi

# ============================================================================
# Local mode (default) — run via conda
# ============================================================================

LOCAL_CMD=(
    mamba run -p "${SCRIPT_DIR}/conda-envs/dana-mag-flye"
    nextflow run "${SCRIPT_DIR}/main.nf"
    --input "$INPUT_HOST"
    --outdir "$OUTDIR_HOST"
    "${NF_ARGS[@]}"
    -resume
)

echo "[INFO] Mode:   Local (conda)"
echo "[INFO] Input:  $INPUT_HOST"
echo "[INFO] Output: $OUTDIR_HOST"
echo "[INFO] Running: ${LOCAL_CMD[*]}"
echo ""

exec "${LOCAL_CMD[@]}"
