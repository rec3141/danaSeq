#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Illumina RNA Pipeline Launcher
# ============================================================================
#
# Runs the Nextflow pipeline locally (default) or via Docker/Apptainer.
# Produces per-reference gene count matrices + viz JSONs.
#
# Usage:
#   ./run-illumina-rna.sh --input DIR --references DIR --outdir DIR [options]
#
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_IMAGE="ghcr.io/rec3141/danaseq-illumina-rna:latest"
USE_CONTAINER=false
CONTAINER_RUNTIME=""
SIF_PATH=""
PULL_SIF=false
NF_ARGS=()
MOUNTS=()
ORIGINAL_ARGS=("$@")

die() { echo "[ERROR] $1" >&2; exit 1; }

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
    cat <<EOF
Usage: $0 --input DIR --references DIR --outdir DIR [pipeline options]

Required:
  --input DIR       Directory containing *_R1_*.fastq.gz paired-end reads
  --references DIR  Directory of <name>.fasta (+ optional <name>.gff)
  --outdir DIR      Output directory (will be created if needed)

Mode:
  --docker          Run in Docker container
  --apptainer       Run in Apptainer/Singularity container
  --container       Auto-detect container runtime
  --image IMAGE     Override container image
  --sif PATH        Use a specific .sif file
  --pull            Pull/build SIF image if not found

Pipeline flags (passed to Nextflow):
  --human_ref PATH         Path to BBTools human reference index
  --sortmerna_refs PATH    Directory of SortMeRNA rRNA FASTAs
  --run_remove_human BOOL  Remove human reads [default: true]
  --run_remove_rrna BOOL   Strip rRNA via SortMeRNA [default: true]
  --run_fastqc BOOL        Run FastQC on preprocessed reads [default: true]
  --min_readlen N          Minimum read length [default: 50]
  --min_identity N         bbmap minid [default: 90]
  --strandedness STR       unstranded | forward | reverse [default: unstranded]
  --feature_type STR       GFF feature for featureCounts [default: CDS]
  --attr_type STR          GFF attribute for grouping [default: locus_tag]
  --rna_cpus N             CPUs per task [default: 16]
  --rna_memory S           Memory per task [default: '64 GB']
  --store_dir DIR          Persistent cache directory (storeDir)

SLURM flags:
  --slurm_account STR  SLURM --account [default: def-rec3141]
  --conda_path PATH    Path to conda/mamba bin/ dir for SLURM jobs

EOF
    exit 0
}

# ============================================================================
# Parse arguments
# ============================================================================

INPUT_HOST=""
REFS_HOST=""
OUTDIR_HOST=""
HUMAN_REF_HOST=""
SORTMERNA_HOST=""
STORE_DIR_HOST=""
RESUME_SESSION=""
AUTO_SESSION=true

while (( $# )); do
    case "$1" in
        -h|--help) usage ;;
        --help-pipeline)
            mamba run -p "${SCRIPT_DIR}/conda-envs/dana-illumina-rna-bbmap" \
                nextflow run "${SCRIPT_DIR}/main.nf" --help
            exit 0 ;;
        --docker)
            USE_CONTAINER=true; CONTAINER_RUNTIME=docker; shift ;;
        --apptainer|--singularity)
            USE_CONTAINER=true; CONTAINER_RUNTIME="${1#--}"; shift ;;
        --container)
            USE_CONTAINER=true; CONTAINER_RUNTIME=auto; shift ;;
        --image)
            [[ -z "${2:-}" ]] && die "--image requires an image name"
            CONTAINER_IMAGE="$2"; shift 2 ;;
        --sif)
            [[ -z "${2:-}" ]] && die "--sif requires a path"
            SIF_PATH="$2"; shift 2 ;;
        --pull)
            PULL_SIF=true; shift ;;
        --input)
            [[ -z "${2:-}" ]] && die "--input requires a directory path"
            INPUT_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"; shift 2 ;;
        --references)
            [[ -z "${2:-}" ]] && die "--references requires a directory path"
            REFS_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"; shift 2 ;;
        --outdir)
            [[ -z "${2:-}" ]] && die "--outdir requires a directory path"
            OUTDIR_HOST="$(realpath -m "$2")"; shift 2 ;;
        --human_ref)
            [[ -z "${2:-}" ]] && die "--human_ref requires a directory path"
            HUMAN_REF_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            NF_ARGS+=("--human_ref" "$HUMAN_REF_HOST"); shift 2 ;;
        --sortmerna_refs)
            [[ -z "${2:-}" ]] && die "--sortmerna_refs requires a directory path"
            SORTMERNA_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            NF_ARGS+=("--sortmerna_refs" "$SORTMERNA_HOST"); shift 2 ;;
        --store_dir)
            [[ -z "${2:-}" ]] && die "--store_dir requires a directory path"
            STORE_DIR_HOST="$(realpath -m "$2")"
            NF_ARGS+=("--store_dir" "$STORE_DIR_HOST"); shift 2 ;;
        --session)
            [[ -z "${2:-}" ]] && die "--session requires a session ID"
            RESUME_SESSION="$2"; AUTO_SESSION=false; shift 2 ;;
        --rna_memory)
            [[ -z "${2:-}" ]] && die "--rna_memory requires a value"
            mem_val="$2"; shift 2
            if [[ "${1:-}" =~ ^[KMGT]B?$ ]]; then mem_val="${mem_val} ${1}"; shift; fi
            NF_ARGS+=("--rna_memory" "${mem_val}") ;;
        *)
            NF_ARGS+=("$1"); shift ;;
    esac
done

# ============================================================================
# Validate
# ============================================================================

[[ -z "$INPUT_HOST" ]]  && die "--input is required."
[[ -z "$REFS_HOST" ]]   && die "--references is required."
[[ -z "$OUTDIR_HOST" ]] && die "--outdir is required."
[[ -d "$INPUT_HOST" ]]  || die "Input directory does not exist: $INPUT_HOST"
[[ -d "$REFS_HOST" ]]   || die "References directory does not exist: $REFS_HOST"

if [[ ! -d "$OUTDIR_HOST" ]]; then
    echo "[INFO] Creating output directory: $OUTDIR_HOST"
    mkdir -p "$OUTDIR_HOST" || die "Cannot create output directory: $OUTDIR_HOST"
fi
[[ -w "$OUTDIR_HOST" ]] || die "Output directory is not writable: $OUTDIR_HOST"

# Auto-detect session
if [[ "$AUTO_SESSION" == true && -z "$RESUME_SESSION" ]]; then
    if [[ -f "${OUTDIR_HOST}/pipeline_info/run_command.sh" ]]; then
        # save_run_command writes "--session <uuid>"; tolerate either prefix
        # and let an empty match through (|| true) so set -e doesn't trip.
        RESUME_SESSION=$(grep -oP '(?<=--session |--?resume )[0-9a-f-]{36}' \
            "${OUTDIR_HOST}/pipeline_info/run_command.sh" | tail -1 || true)
        if [[ -n "$RESUME_SESSION" ]]; then
            echo "[INFO] Auto-detected session from previous run: $RESUME_SESSION"
        fi
    fi
fi

# ============================================================================
# Container runtime
# ============================================================================

if [[ "$USE_CONTAINER" == true ]]; then
    if [[ "$CONTAINER_RUNTIME" == "auto" ]]; then
        if command -v apptainer &>/dev/null; then CONTAINER_RUNTIME=apptainer
        elif command -v singularity &>/dev/null; then CONTAINER_RUNTIME=singularity
        elif command -v docker &>/dev/null; then CONTAINER_RUNTIME=docker
        else die "--container requires docker, apptainer, or singularity on PATH"; fi
        echo "[INFO] Auto-detected container runtime: ${CONTAINER_RUNTIME}"
    fi

    command -v "$CONTAINER_RUNTIME" &>/dev/null \
        || die "${CONTAINER_RUNTIME} not found on PATH"

    if [[ "$CONTAINER_RUNTIME" == "apptainer" || "$CONTAINER_RUNTIME" == "singularity" ]]; then
        [[ -z "$SIF_PATH" ]] && SIF_PATH="${SCRIPT_DIR}/.danaseq-illumina-rna.sif"
        if [[ ! -f "$SIF_PATH" ]]; then
            if [[ "$PULL_SIF" == true ]]; then
                "$CONTAINER_RUNTIME" pull "$SIF_PATH" "docker://${CONTAINER_IMAGE}"
            else
                die "SIF file not found: $SIF_PATH. Use --pull to download."
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
    BINDS+=("${REFS_HOST}:/data/references:ro")
    BINDS+=("${OUTDIR_HOST}:/data/output")
    NF_ARGS=("--input" "/data/input" "--references" "/data/references" "--outdir" "/data/output" "${NF_ARGS[@]}")

    if [[ -n "${HUMAN_REF_HOST:-}" ]]; then
        BINDS+=("${HUMAN_REF_HOST}:/data/human_ref:ro")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${HUMAN_REF_HOST}/\/data\/human_ref}"
        done
    fi
    if [[ -n "${SORTMERNA_HOST:-}" ]]; then
        BINDS+=("${SORTMERNA_HOST}:/data/sortmerna_refs:ro")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${SORTMERNA_HOST}/\/data\/sortmerna_refs}"
        done
    fi
    if [[ -n "${STORE_DIR_HOST:-}" ]]; then
        BINDS+=("${STORE_DIR_HOST}:/data/store")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${STORE_DIR_HOST}/\/data\/store}"
        done
    fi

    WORKDIR_HOST="${OUTDIR_HOST}/work"
    mkdir -p "${WORKDIR_HOST}" 2>/dev/null || true
    BINDS+=("${WORKDIR_HOST}:/data/work")

    NF_CACHE="${OUTDIR_HOST}/.nextflow-cache"
    mkdir -p "${NF_CACHE}/dotdir" 2>/dev/null || true
    BINDS+=("${NF_CACHE}/dotdir:/home/dana/.nextflow")

    CONTAINER_CMD=()
    case "$CONTAINER_RUNTIME" in
        docker)
            CONTAINER_CMD+=(docker run --user "$(id -u):$(id -g)")
            CONTAINER_CMD+=("-e" "NXF_HOME=/home/dana/.nextflow")
            for bind in "${BINDS[@]}"; do CONTAINER_CMD+=("-v" "$bind"); done
            CONTAINER_CMD+=("$CONTAINER_IMAGE" run /pipeline/main.nf) ;;
        apptainer|singularity)
            container_ca="/etc/ssl/certs/ca-certificates.crt"
            CONTAINER_CMD+=("$CONTAINER_RUNTIME" run)
            CONTAINER_CMD+=("--env" "NXF_HOME=/home/dana/.nextflow")
            CONTAINER_CMD+=("--env" "REQUESTS_CA_BUNDLE=${container_ca}")
            CONTAINER_CMD+=("--env" "SSL_CERT_FILE=${container_ca}")
            CONTAINER_CMD+=("--env" "CURL_CA_BUNDLE=${container_ca}")
            for bind in "${BINDS[@]}"; do CONTAINER_CMD+=("--bind" "$bind"); done
            CONTAINER_CMD+=("$SIF_PATH" run /pipeline/main.nf) ;;
    esac
    CONTAINER_CMD+=(-w /data/work "${NF_ARGS[@]}" -resume ${RESUME_SESSION})

    echo "[INFO] Mode:       Container (${CONTAINER_RUNTIME})"
    echo "[INFO] Input:      $INPUT_HOST"
    echo "[INFO] References: $REFS_HOST"
    echo "[INFO] Output:     $OUTDIR_HOST"
    echo ""

    mkdir -p "${OUTDIR_HOST}/pipeline_info" 2>/dev/null || true
    "${CONTAINER_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

    NF_SESSION=$(awk '{print $6}' "${NF_CACHE}/dotdir/history" 2>/dev/null | tail -1 || true)
    save_run_command "$OUTDIR_HOST" "$NF_SESSION"
    exit $NF_EXIT
fi

# ============================================================================
# Local mode (default)
# ============================================================================

CONDA_BASE_BIN="$(dirname "$(which mamba 2>/dev/null || which conda)")"
export PATH="${CONDA_BASE_BIN}:${PATH}"

BBMAP_ENV="${SCRIPT_DIR}/conda-envs/dana-illumina-rna-bbmap"
if [[ -d "${BBMAP_ENV}/lib/jvm" ]]; then
    export NXF_JAVA_HOME="${BBMAP_ENV}/lib/jvm"
elif [[ -x "${BBMAP_ENV}/bin/java" ]]; then
    export NXF_JAVA_HOME="${BBMAP_ENV}"
fi

LOCAL_CMD=(
    mamba run -p "${BBMAP_ENV}"
    nextflow run "${SCRIPT_DIR}/main.nf"
    --input "$INPUT_HOST"
    --references "$REFS_HOST"
    --outdir "$OUTDIR_HOST"
    "${NF_ARGS[@]}"
    -resume ${RESUME_SESSION}
)

echo "[INFO] Mode:       Local (conda)"
echo "[INFO] Input:      $INPUT_HOST"
echo "[INFO] References: $REFS_HOST"
echo "[INFO] Output:     $OUTDIR_HOST"
echo ""

mkdir -p "${OUTDIR_HOST}/pipeline_info" 2>/dev/null || true
"${LOCAL_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${SCRIPT_DIR}/.nextflow.log" 2>/dev/null | tail -1 || true)
save_run_command "$OUTDIR_HOST" "$NF_SESSION"

exit $NF_EXIT
