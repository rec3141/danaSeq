#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana Nanopore Assembly Pipeline Launcher
# ============================================================================
#
# Runs the Nextflow assembly pipeline locally (default) or via Docker.
# Produces assembly + depth table + BAMs for downstream mag_analysis.
#
# Usage:
#   ./run-nanopore-assembly.sh --input /path/to/reads --outdir /path/to/output [options]
#
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_IMAGE="ghcr.io/rec3141/danaseq-nanopore-assembly:latest"
USE_CONTAINER=false
CONTAINER_RUNTIME=""
SIF_PATH=""
PULL_SIF=false
NF_ARGS=()
MOUNTS=()
ORIGINAL_ARGS=("$@")

die() { echo "[ERROR] $1" >&2; exit 1; }

# Build a re-runnable self-invocation with the correct --session ID
save_run_command() {
    local outdir="$1" session="$2"
    local save_args=()
    local skip_next=false has_session=false
    for arg in "${ORIGINAL_ARGS[@]}"; do
        if $skip_next; then skip_next=false; continue; fi
        if [[ "$arg" == "--resume" ]]; then
            has_session=true
            continue
        fi
        if $has_session && [[ "$arg" =~ ^[0-9a-f-]{36}$ ]]; then
            has_session=false
            continue
        fi
        has_session=false
        save_args+=("$arg")
    done
    if [[ -n "$session" ]]; then
        save_args+=("--resume" "$session")
    fi
    printf '%s\n' "$(realpath "$0") ${save_args[*]}" >> "${outdir}/pipeline_info/run_command.txt"
}

usage() {
    echo "Usage: $0 --input DIR --outdir DIR [options]"
    echo ""
    echo "Required:"
    echo "  --input DIR      Directory containing reads (*.fastq.gz or barcode structure)"
    echo "  --outdir DIR     Output directory (will be created if needed)"
    echo ""
    echo "Mode:"
    echo "  --docker         Run in Docker container"
    echo "  --apptainer      Run in Apptainer/Singularity container"
    echo "  --container      Auto-detect container runtime"
    echo "  --image IMAGE    Override container image"
    echo "  --sif PATH       Use a specific .sif file"
    echo "  --pull           Pull/build SIF image if not found"
    echo ""
    echo "Caching & Resume:"
    echo "  --workdir DIR        Nextflow work directory (default: /tmp/nanopore_assembly_work)"
    echo "  --store_dir DIR      Persistent cache directory (storeDir)"
    echo "  --resume [ID]        Resume a previous run"
    echo ""
    echo "Pipeline flags (passed to Nextflow):"
    echo "  --assembler STR      flye, metamdbg, or myloasm [default: flye]"
    echo "  --dedupe             BBDuk deduplication before assembly"
    echo "  --filtlong_size N    Filtlong target bases; skip if not set"
    echo "  --min_overlap N      Flye --min-overlap [default: 1000]"
    echo "  --polish BOOL        Run Flye polisher [default: auto]"
    echo "  --assembly_cpus N    CPUs for assembly [default: 16]"
    echo "  --assembly_memory S  Memory for assembly [default: '60 GB']"
    echo ""
    echo "Output feeds into mag_analysis:"
    echo "  mag_analysis/nextflow/run-mag-analysis.sh \\"
    echo "      --assembly <outdir>/assembly/assembly.fasta \\"
    echo "      --depths <outdir>/mapping/depths.txt \\"
    echo "      --bam_dir <outdir>/mapping/"
    echo ""
    exit 0
}

# ============================================================================
# Parse arguments
# ============================================================================

INPUT_HOST=""
OUTDIR_HOST=""
WORKDIR_HOST="/tmp/nanopore_assembly_work"
RESUME_SESSION=""
DO_RESUME=false
STORE_DIR_HOST=""

while (( $# )); do
    case "$1" in
        -h|--help)
            usage ;;
        --docker)
            USE_CONTAINER=true; CONTAINER_RUNTIME=docker
            shift ;;
        --apptainer|--singularity)
            USE_CONTAINER=true; CONTAINER_RUNTIME="${1#--}"
            shift ;;
        --container)
            USE_CONTAINER=true; CONTAINER_RUNTIME=auto
            shift ;;
        --image)
            [[ -z "${2:-}" ]] && die "--image requires an image name"
            CONTAINER_IMAGE="$2"
            shift 2 ;;
        --sif)
            [[ -z "${2:-}" ]] && die "--sif requires a path to a .sif file"
            SIF_PATH="$2"
            shift 2 ;;
        --pull)
            PULL_SIF=true
            shift ;;
        --input)
            [[ -z "${2:-}" ]] && die "--input requires a directory path"
            INPUT_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            shift 2 ;;
        --outdir)
            [[ -z "${2:-}" ]] && die "--outdir requires a directory path"
            OUTDIR_HOST="$(realpath -m "$2")"
            shift 2 ;;
        --workdir|-w)
            [[ -z "${2:-}" ]] && die "--workdir requires a directory path"
            WORKDIR_HOST="$(realpath -m "$2")"
            shift 2 ;;
        --store_dir)
            [[ -z "${2:-}" ]] && die "--store_dir requires a directory path"
            STORE_DIR_HOST="$(realpath -m "$2")"
            NF_ARGS+=("--store_dir" "$STORE_DIR_HOST")
            shift 2 ;;
        --resume)
            DO_RESUME=true
            if [[ -n "${2:-}" && "${2}" != --* ]]; then
                RESUME_SESSION="$2"
                shift
            fi
            shift ;;
        --run_flye_polish)
            NF_ARGS+=("--polish" "$2")
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

if [[ -n "$OUTDIR_HOST" && -n "$STORE_DIR_HOST" ]]; then
    die "--outdir and --store_dir are mutually exclusive."
fi

if [[ -z "$OUTDIR_HOST" && -n "$STORE_DIR_HOST" ]]; then
    OUTDIR_HOST="$STORE_DIR_HOST"
    NF_ARGS+=("--outdir" "$OUTDIR_HOST")
fi
[[ -z "$OUTDIR_HOST" ]] && die "--outdir (or --store_dir) is required."

[[ -d "$INPUT_HOST" ]] || die "Input directory does not exist: $INPUT_HOST"

if [[ ! -d "$OUTDIR_HOST" ]]; then
    echo "[INFO] Creating output directory: $OUTDIR_HOST"
    mkdir -p "$OUTDIR_HOST" || die "Cannot create output directory: $OUTDIR_HOST"
fi
[[ -w "$OUTDIR_HOST" ]] || die "Output directory is not writable: $OUTDIR_HOST"

# ============================================================================
# Auto-detect session ID from previous runs
# ============================================================================

if [[ "$DO_RESUME" == true && -z "$RESUME_SESSION" ]]; then
    _base="${STORE_DIR_HOST:-$OUTDIR_HOST}"
    for _hist in "${_base}/.nextflow-cache/dotdir/history" \
                 "${_base}/.nextflow/history" \
                 "${OUTDIR_HOST}/.nextflow-cache/dotdir/history" \
                 "${OUTDIR_HOST}/.nextflow/history"; do
        if [[ -f "$_hist" ]]; then
            RESUME_SESSION=$(awk '{print $6}' "$_hist" | tail -1)
            break
        fi
    done
    if [[ -n "$RESUME_SESSION" ]]; then
        echo "[INFO] Auto-detected session from Nextflow history: $RESUME_SESSION"
    fi
fi

# ============================================================================
# Container runtime detection
# ============================================================================

if [[ "$USE_CONTAINER" == true ]]; then
    if [[ "$CONTAINER_RUNTIME" == "auto" ]]; then
        if command -v apptainer &>/dev/null; then
            CONTAINER_RUNTIME=apptainer
        elif command -v singularity &>/dev/null; then
            CONTAINER_RUNTIME=singularity
        elif command -v docker &>/dev/null; then
            CONTAINER_RUNTIME=docker
        else
            die "--container requires docker, apptainer, or singularity on PATH"
        fi
        echo "[INFO] Auto-detected container runtime: ${CONTAINER_RUNTIME}"
    fi

    if ! command -v "$CONTAINER_RUNTIME" &>/dev/null; then
        die "${CONTAINER_RUNTIME} not found on PATH"
    fi

    if [[ "$CONTAINER_RUNTIME" == "apptainer" || "$CONTAINER_RUNTIME" == "singularity" ]]; then
        if [[ -z "$SIF_PATH" ]]; then
            SIF_PATH="${SCRIPT_DIR}/.danaseq-nanopore-assembly.sif"
        fi
        if [[ ! -f "$SIF_PATH" ]]; then
            if [[ "$PULL_SIF" == true ]]; then
                echo "[INFO] Pulling container image..."
                "$CONTAINER_RUNTIME" pull "$SIF_PATH" "docker://${CONTAINER_IMAGE}"
            else
                die "SIF file not found: $SIF_PATH\n  Use --pull to download it, or --sif PATH."
            fi
        fi
    fi
fi

# ============================================================================
# Container mode
# ============================================================================

if [[ "$USE_CONTAINER" == true ]]; then
    BINDS=()
    BINDS+=("${INPUT_HOST}:/data/input:ro")
    BINDS+=("${OUTDIR_HOST}:/data/output")
    NF_ARGS=("--input" "/data/input" "--outdir" "/data/output" "${NF_ARGS[@]}")

    if [[ -n "${STORE_DIR_HOST:-}" ]]; then
        BINDS+=("${STORE_DIR_HOST}:/data/store")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${STORE_DIR_HOST}/\/data\/store}"
        done
    fi

    mkdir -p "${WORKDIR_HOST}" 2>/dev/null || true
    BINDS+=("${WORKDIR_HOST}:/data/work")

    NF_CACHE="${OUTDIR_HOST}/.nextflow-cache"
    mkdir -p "${NF_CACHE}/dotdir" 2>/dev/null || true

    if [[ "$DO_RESUME" == true && ! -f "${NF_CACHE}/dotdir/history" && -f "$HOME/.nextflow/history" ]]; then
        echo "[INFO] Migrating Nextflow cache from ~/.nextflow to per-run NXF_HOME"
        cp "$HOME/.nextflow/history" "${NF_CACHE}/dotdir/" 2>/dev/null || true
        if [[ -n "${RESUME_SESSION:-}" && -d "$HOME/.nextflow/cache/$RESUME_SESSION" ]]; then
            mkdir -p "${NF_CACHE}/dotdir/cache" 2>/dev/null || true
            ln -sfn "$HOME/.nextflow/cache/$RESUME_SESSION" "${NF_CACHE}/dotdir/cache/$RESUME_SESSION"
        elif [[ -d "$HOME/.nextflow/cache" ]]; then
            ln -sfn "$HOME/.nextflow/cache" "${NF_CACHE}/dotdir/cache"
        fi
    fi

    BINDS+=("${NF_CACHE}/dotdir:/home/dana/.nextflow")

    CONTAINER_CMD=()
    case "$CONTAINER_RUNTIME" in
        docker)
            CONTAINER_CMD+=(docker run --user "$(id -u):$(id -g)")
            CONTAINER_CMD+=("-e" "NXF_HOME=/home/dana/.nextflow")
            for bind in "${BINDS[@]}"; do
                CONTAINER_CMD+=("-v" "$bind")
            done
            CONTAINER_CMD+=("$CONTAINER_IMAGE" run /pipeline/main.nf)
            ;;
        apptainer|singularity)
            container_ca="/etc/ssl/certs/ca-certificates.crt"
            CONTAINER_CMD+=("$CONTAINER_RUNTIME" run)
            CONTAINER_CMD+=("--env" "NXF_HOME=/home/dana/.nextflow")
            CONTAINER_CMD+=("--env" "REQUESTS_CA_BUNDLE=${container_ca}")
            CONTAINER_CMD+=("--env" "SSL_CERT_FILE=${container_ca}")
            CONTAINER_CMD+=("--env" "CURL_CA_BUNDLE=${container_ca}")
            for bind in "${BINDS[@]}"; do
                CONTAINER_CMD+=("--bind" "$bind")
            done
            CONTAINER_CMD+=("$SIF_PATH" run /pipeline/main.nf)
            ;;
    esac
    CONTAINER_CMD+=(-w /data/work "${NF_ARGS[@]}")
    if [[ "$DO_RESUME" == true ]]; then
        CONTAINER_CMD+=(-resume ${RESUME_SESSION})
    fi

    echo "[INFO] Mode:   Container (${CONTAINER_RUNTIME})"
    echo "[INFO] Input:  $INPUT_HOST"
    echo "[INFO] Output: $OUTDIR_HOST"
    echo "[INFO] Running: ${CONTAINER_CMD[*]}"
    echo ""

    mkdir -p "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info" 2>/dev/null || true

    "${CONTAINER_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

    if [[ $NF_EXIT -ne 0 ]]; then
        _trace="${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info/trace.txt"
        if [[ -f "$_trace" ]] && ! grep -q 'FAILED' "$_trace" 2>/dev/null; then
            echo "[WARN] Container exited $NF_EXIT but Nextflow trace shows no failures"
            NF_EXIT=0
        fi
    fi

    NF_SESSION=$(awk '{print $6}' "${NF_CACHE}/dotdir/history" 2>/dev/null | tail -1)
    [[ -z "$NF_SESSION" ]] && NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${SCRIPT_DIR}/.nextflow.log" 2>/dev/null | tail -1)
    save_run_command "${STORE_DIR_HOST:-$OUTDIR_HOST}" "$NF_SESSION"

    exit $NF_EXIT
fi

# ============================================================================
# Local mode (default) — run via conda
# ============================================================================

CONDA_BASE_BIN="$(dirname "$(which mamba 2>/dev/null || which conda)")"
export PATH="${CONDA_BASE_BIN}:${PATH}"

WORKDIR_FLAG=()
if [[ -n "$WORKDIR_HOST" ]]; then
    mkdir -p "$WORKDIR_HOST" || die "Cannot create work directory: $WORKDIR_HOST"
    WORKDIR_FLAG=(-w "$WORKDIR_HOST")
fi

RESUME_FLAG=()
if [[ "$DO_RESUME" == true ]]; then
    RESUME_FLAG=(-resume ${RESUME_SESSION})
fi

LOCAL_CMD=(
    mamba run -p "${SCRIPT_DIR}/conda-envs/dana-mag-assembly"
    nextflow run "${SCRIPT_DIR}/main.nf"
    --input "$INPUT_HOST"
    --outdir "$OUTDIR_HOST"
    "${WORKDIR_FLAG[@]}"
    "${NF_ARGS[@]}"
    "${RESUME_FLAG[@]}"
)

echo "[INFO] Mode:   Local (conda)"
echo "[INFO] Input:  $INPUT_HOST"
echo "[INFO] Output: $OUTDIR_HOST"
echo "[INFO] Running: ${LOCAL_CMD[*]}"
echo ""

mkdir -p "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info" 2>/dev/null || true

"${LOCAL_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

# Capture session ID
NF_SESSION=$(awk '{print $6}' .nextflow/history 2>/dev/null | tail -1)
[[ -z "$NF_SESSION" ]] && NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' .nextflow.log 2>/dev/null | tail -1)
save_run_command "${STORE_DIR_HOST:-$OUTDIR_HOST}" "$NF_SESSION"

exit $NF_EXIT
