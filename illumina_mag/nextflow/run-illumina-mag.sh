#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Illumina MAG Pipeline Launcher
# ============================================================================
#
# Runs the Nextflow pipeline locally (default) or via Docker/Apptainer.
#
# Usage:
#   ./run-illumina-mag.sh --input /path/to/reads --outdir /path/to/output [options]
#
# Examples:
#   # Basic per-sample assembly (local conda)
#   ./run-illumina-mag.sh --input /data/reads --outdir /data/output
#
#   # Co-assembly mode
#   ./run-illumina-mag.sh --input /data/reads --outdir /data/output --coassembly
#
#   # Docker mode
#   ./run-illumina-mag.sh --docker --input /data/reads --outdir /data/output
#
#   # Apptainer (HPC)
#   ./run-illumina-mag.sh --apptainer --input /data/reads --outdir /data/output
#
#   # Skip human removal and some assemblers
#   ./run-illumina-mag.sh --input /data/reads --outdir /data/output \
#       --run_remove_human false --run_megahit false --run_spades false
#
#   # SLURM profile
#   ./run-illumina-mag.sh --input /data/reads --outdir /data/output \
#       -profile slurm --slurm_account def-rec3141 \
#       --conda_path ~/scratch/miniforge3/bin
#
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_IMAGE="ghcr.io/rec3141/danaseq-illumina-mag:latest"
USE_CONTAINER=false
CONTAINER_RUNTIME=""   # docker, apptainer, or singularity
SIF_PATH=""            # resolved path to .sif file (apptainer/singularity only)
PULL_SIF=false
NF_ARGS=()
MOUNTS=()
ORIGINAL_ARGS=("$@")

die() { echo "[ERROR] $1" >&2; exit 1; }

# Build a re-runnable self-invocation with the correct --session ID
save_run_command() {
    local outdir="$1" session="$2"
    local save_args=()
    local skip_next=false
    for arg in "${ORIGINAL_ARGS[@]}"; do
        if $skip_next; then skip_next=false; continue; fi
        if [[ "$arg" == "--session" ]]; then skip_next=true; continue; fi
        save_args+=("$arg")
    done
    if [[ -n "$session" ]]; then
        save_args+=("--session" "$session")
    fi
    printf '%s\n' "$(realpath "$0") ${save_args[*]}" >> "${outdir}/pipeline_info/run_command.sh"
}

usage() {
    echo "Usage: $0 --input DIR --outdir DIR [pipeline options]"
    echo ""
    echo "Required:"
    echo "  --input DIR      Directory containing *_R1_*.fastq.gz paired-end reads"
    echo "  --outdir DIR     Output directory (will be created if needed)"
    echo ""
    echo "Mode:"
    echo "  --docker         Run in Docker container"
    echo "  --apptainer      Run in Apptainer/Singularity container"
    echo "  --container      Auto-detect container runtime (apptainer > singularity > docker)"
    echo "  --image IMAGE    Override container image [default: $CONTAINER_IMAGE]"
    echo "  --sif PATH       Use a specific .sif file (apptainer/singularity only)"
    echo "  --pull           Pull/build SIF image if not found (apptainer/singularity)"
    echo ""
    echo "Pipeline flags (passed to Nextflow):"
    echo "  --coassembly              Co-assemble all samples (default: per-sample)"
    echo "  --human_ref PATH          Path to BBTools human reference index"
    echo "  --run_remove_human BOOL   Remove human reads [default: true]"
    echo "  --run_fastqc BOOL         Run FastQC on preprocessed reads [default: true]"
    echo "  --run_normalize BOOL      Enable bbnorm normalization [default: true]"
    echo "  --run_tadpole BOOL        Run Tadpole assembler [default: true]"
    echo "  --run_megahit BOOL        Run Megahit assembler [default: true]"
    echo "  --run_spades BOOL         Run SPAdes assembler [default: true]"
    echo "  --run_metaspades BOOL     Run metaSPAdes assembler [default: true]"
    echo "  --min_readlen N           Minimum read length [default: 70]"
    echo "  --dedupe_identity N       Deduplication identity threshold [default: 98]"
    echo "  --metabat_min_cls N       MetaBAT2 minimum cluster size [default: 2000]"
    echo "  --assembly_cpus N         CPUs for assembly [default: 24]"
    echo "  --assembly_memory S       Memory for assembly [default: '250 GB']"
    echo "  --store_dir DIR           Persistent cache directory (storeDir)"
    echo ""
    echo "SLURM flags:"
    echo "  --slurm_account STR  SLURM --account [default: def-rec3141]"
    echo "  --conda_path PATH    Path to conda/mamba bin/ dir for SLURM jobs"
    echo "                       (e.g. ~/scratch/miniforge3/bin)"
    echo ""
    echo "Run '$0 --help-pipeline' to see full Nextflow help."
    exit 0
}

# ============================================================================
# Parse arguments -- extract paths that need special handling
# ============================================================================

INPUT_HOST=""
OUTDIR_HOST=""
HUMAN_REF_HOST=""
STORE_DIR_HOST=""
RESUME_SESSION=""
AUTO_SESSION=true

while (( $# )); do
    case "$1" in
        -h|--help)
            usage ;;
        --help-pipeline)
            if [[ "$USE_CONTAINER" == true ]]; then
                docker run --user "$(id -u):$(id -g)" "$CONTAINER_IMAGE" \
                    run /pipeline/main.nf --help
            else
                mamba run -p "${SCRIPT_DIR}/conda-envs/dana-illumina-mag-bbmap" \
                    nextflow run "${SCRIPT_DIR}/main.nf" --help
            fi
            exit 0 ;;
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
        --human_ref)
            [[ -z "${2:-}" ]] && die "--human_ref requires a directory path"
            HUMAN_REF_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            NF_ARGS+=("--human_ref" "$HUMAN_REF_HOST")
            shift 2 ;;
        --store_dir)
            [[ -z "${2:-}" ]] && die "--store_dir requires a directory path"
            STORE_DIR_HOST="$(realpath -m "$2")"
            NF_ARGS+=("--store_dir" "$STORE_DIR_HOST")
            shift 2 ;;
        --session)
            [[ -z "${2:-}" ]] && die "--session requires a session ID"
            RESUME_SESSION="$2"
            AUTO_SESSION=false
            shift 2 ;;
        --assembly_memory)
            [[ -z "${2:-}" ]] && die "--assembly_memory requires a value (e.g. '50 GB')"
            mem_val="$2"; shift 2
            if [[ "${1:-}" =~ ^[KMGT]B?$ ]]; then
                mem_val="${mem_val} ${1}"
                shift
            fi
            NF_ARGS+=("--assembly_memory" "${mem_val}")
            ;;
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
# Auto-detect session ID from previous runs in this outdir
# ============================================================================

if [[ "$AUTO_SESSION" == true && -z "$RESUME_SESSION" ]]; then
    if [[ -f "${OUTDIR_HOST}/pipeline_info/run_command.sh" ]]; then
        RESUME_SESSION=$(grep -oP '(?<=-resume )[0-9a-f-]{36}' \
            "${OUTDIR_HOST}/pipeline_info/run_command.sh" | tail -1)
        if [[ -n "$RESUME_SESSION" ]]; then
            echo "[INFO] Auto-detected session from previous run: $RESUME_SESSION"
        fi
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

    # For apptainer/singularity, ensure we have a SIF image
    if [[ "$CONTAINER_RUNTIME" == "apptainer" || "$CONTAINER_RUNTIME" == "singularity" ]]; then
        if [[ -z "$SIF_PATH" ]]; then
            SIF_PATH="${SCRIPT_DIR}/.danaseq-illumina-mag.sif"
        fi
        if [[ ! -f "$SIF_PATH" ]]; then
            if [[ "$PULL_SIF" == true ]]; then
                echo "[INFO] Pulling container image..."
                echo "  Source: docker://${CONTAINER_IMAGE}"
                echo "  Destination: ${SIF_PATH}"
                "$CONTAINER_RUNTIME" pull "$SIF_PATH" "docker://${CONTAINER_IMAGE}"
            else
                die "SIF file not found: $SIF_PATH\n  Use --pull to download it, or --sif PATH to specify an existing file."
            fi
        fi
        [[ -f "$SIF_PATH" ]] || die "SIF file not found: $SIF_PATH"
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

    # Mount human reference if provided and rewrite path
    if [[ -n "${HUMAN_REF_HOST:-}" ]]; then
        BINDS+=("${HUMAN_REF_HOST}:/data/human_ref:ro")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${HUMAN_REF_HOST}/\/data\/human_ref}"
        done
    fi

    # Mount store directory (read-write) and rewrite paths
    if [[ -n "${STORE_DIR_HOST:-}" ]]; then
        BINDS+=("${STORE_DIR_HOST}:/data/store")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${STORE_DIR_HOST}/\/data\/store}"
        done
    fi

    # Persistent Nextflow .nextflow directory for -resume support
    NF_CACHE="${OUTDIR_HOST}/.nextflow-cache"
    mkdir -p "${NF_CACHE}/dotdir" 2>/dev/null || true
    BINDS+=("${NF_CACHE}/dotdir:/home/dana/.nextflow")

    # Build runtime-specific command
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
            CONTAINER_CMD+=("--env" "NXF_HOME=${NF_CACHE}/dotdir")
            CONTAINER_CMD+=("--env" "REQUESTS_CA_BUNDLE=${container_ca}")
            CONTAINER_CMD+=("--env" "SSL_CERT_FILE=${container_ca}")
            CONTAINER_CMD+=("--env" "CURL_CA_BUNDLE=${container_ca}")
            for bind in "${BINDS[@]}"; do
                CONTAINER_CMD+=("--bind" "$bind")
            done
            CONTAINER_CMD+=("$SIF_PATH" run /pipeline/main.nf)
            ;;
    esac
    CONTAINER_CMD+=("${NF_ARGS[@]}" -resume ${RESUME_SESSION})

    echo "[INFO] Mode:   Container (${CONTAINER_RUNTIME})"
    echo "[INFO] Input:  $INPUT_HOST"
    echo "[INFO] Output: $OUTDIR_HOST"
    echo "[INFO] Running: ${CONTAINER_CMD[*]}"
    echo ""

    mkdir -p "${OUTDIR_HOST}/pipeline_info" 2>/dev/null || true

    "${CONTAINER_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

    # Capture session ID from Nextflow history
    NF_SESSION=$(awk '{print $6}' "${NF_CACHE}/dotdir/history" 2>/dev/null | tail -1 || true)
    save_run_command "$OUTDIR_HOST" "$NF_SESSION"

    exit $NF_EXIT
fi

# ============================================================================
# Local mode (default) — run via conda
# ============================================================================

# Nextflow's conda activation calls `conda info --json` to find the base prefix.
# Ensure conda/mamba base bin/ is on PATH so this works inside `mamba run`.
CONDA_BASE_BIN="$(dirname "$(which mamba 2>/dev/null || which conda)")"
export PATH="${CONDA_BASE_BIN}:${PATH}"

# Point Nextflow at the conda env's Java so it never picks up system Java.
BBMAP_ENV="${SCRIPT_DIR}/conda-envs/dana-illumina-mag-bbmap"
if [[ -d "${BBMAP_ENV}/lib/jvm" ]]; then
    export NXF_JAVA_HOME="${BBMAP_ENV}/lib/jvm"
elif [[ -x "${BBMAP_ENV}/bin/java" ]]; then
    export NXF_JAVA_HOME="${BBMAP_ENV}"
fi

LOCAL_CMD=(
    mamba run -p "${SCRIPT_DIR}/conda-envs/dana-illumina-mag-bbmap"
    nextflow run "${SCRIPT_DIR}/main.nf"
    --input "$INPUT_HOST"
    --outdir "$OUTDIR_HOST"
    "${NF_ARGS[@]}"
    -resume ${RESUME_SESSION}
)

echo "[INFO] Mode:   Local (conda)"
echo "[INFO] Input:  $INPUT_HOST"
echo "[INFO] Output: $OUTDIR_HOST"
echo "[INFO] Running: ${LOCAL_CMD[*]}"
echo ""

mkdir -p "${OUTDIR_HOST}/pipeline_info" 2>/dev/null || true

"${LOCAL_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

# Capture session ID from Nextflow log and save self-invocation for reliable resume
NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${SCRIPT_DIR}/.nextflow.log" 2>/dev/null | tail -1 || true)
save_run_command "$OUTDIR_HOST" "$NF_SESSION"

exit $NF_EXIT
