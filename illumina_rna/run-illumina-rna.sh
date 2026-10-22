#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Illumina RNA Pipeline Launcher
# ============================================================================
#
# Runs the Nextflow pipeline locally (default) or via Docker/Apptainer.
# Produces per-reference expression matrices + viz JSONs.
#
# Usage:
#   ./run-illumina-rna.sh --input DIR --from_mag DIR --outdir DIR [options]
#
# In container mode every path argument is bound at the same path inside the
# container (plus its resolved target, for symlinked inputs), so Nextflow sees
# exactly the paths given here. Extra mounts: --bind PATH.
#
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_IMAGE="ghcr.io/rec3141/danaseq-illumina-rna:latest"
USE_CONTAINER=false
CONTAINER_RUNTIME=""
SIF_PATH=""
PULL_SIF=false
NF_ARGS=()
BIND_PATHS=()
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
Usage: $0 --input DIR (--from_mag DIR | --from_assembly DIR | --references DIR) --outdir DIR [options]

Reads (one of):
  --input DIR          Directory searched recursively for paired reads
  --samplesheet CSV    sample,r1,r2

References (one of):
  --from_mag DIR       mag_analysis output (assembly + Bakta GFF + bins)
  --from_assembly DIR  nanopore_assembly / illumina_assembly output
  --references DIR     <name>.fasta (+ <name>.gff, <name>.bins.tsv)

Required:
  --outdir DIR         Output directory (created if needed)

Mode:
  --docker             Run in Docker container
  --apptainer          Run in Apptainer/Singularity container
  --container          Auto-detect container runtime
  --image IMAGE        Override container image
  --sif PATH           Use a specific .sif file
  --pull               Pull/build SIF image if not found
  --bind PATH          Extra host path to bind into the container (repeatable)
  --workdir DIR        Nextflow work dir [default: <outdir>/work]

Pipeline flags (passed to Nextflow; see --help-pipeline for all):
  --human_ref PATH         Human reference FASTA (indexed once)
  --rrna_refs PATH         rRNA FASTA or dir (e.g. SILVA SSU+LSU NR99)
  --rrna_method STR        bbduk | sortmerna | none [default: bbduk]
  --run_remove_human BOOL  Remove human reads [default: true]
  --max_read_pairs N       Use only the first N pairs per sample
  --chunk I/N              Process one of N sample chunks (with --summarize false)
  --mapper STR             strobealign | bbmap [default: strobealign]
  --rna_cpus N / --rna_memory S   Per-task resources [default: 16 / '64 GB']
  --map_cpus N / --map_memory S   Mapping/index resources [default: rna_*]
  --store_dir DIR          Persistent cache directory (storeDir)

SLURM flags:
  --slurm_account STR  SLURM --account (required with -profile slurm)
  --conda_path PATH    Path to conda/mamba bin/ dir for SLURM jobs

EOF
    exit 0
}

# ============================================================================
# Parse arguments
# ============================================================================

OUTDIR_HOST=""
WORKDIR_HOST=""
RESUME_SESSION=""
AUTO_SESSION=true

# Flags whose value is a host path: resolved to absolute, and bound in
# container mode.
PATH_FLAGS=" --input --samplesheet --references --from_mag --from_assembly --human_ref --rrna_refs --store_dir --contig_explorer "

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
        --bind)
            [[ -z "${2:-}" ]] && die "--bind requires a path"
            BIND_PATHS+=("$(realpath -m "$2")"); shift 2 ;;
        --outdir)
            [[ -z "${2:-}" ]] && die "--outdir requires a directory path"
            OUTDIR_HOST="$(realpath -m "$2")"; shift 2 ;;
        --workdir)
            [[ -z "${2:-}" ]] && die "--workdir requires a directory path"
            WORKDIR_HOST="$(realpath -m "$2")"; shift 2 ;;
        --session)
            [[ -z "${2:-}" ]] && die "--session requires a session ID"
            RESUME_SESSION="$2"; AUTO_SESSION=false; shift 2 ;;
        --rna_memory|--map_memory)
            [[ -z "${2:-}" ]] && die "$1 requires a value"
            flag="$1"; mem_val="$2"; shift 2
            if [[ "${1:-}" =~ ^[KMGT]B?$ ]]; then mem_val="${mem_val} ${1}"; shift; fi
            NF_ARGS+=("$flag" "${mem_val}") ;;
        *)
            if [[ "$PATH_FLAGS" == *" $1 "* ]]; then
                [[ -z "${2:-}" ]] && die "$1 requires a path"
                p="$(realpath -m "$2")"
                [[ -e "$p" ]] || [[ "$1" == "--store_dir" ]] || die "$1: path does not exist: $2"
                NF_ARGS+=("$1" "$p")
                BIND_PATHS+=("$p")
                shift 2
            else
                NF_ARGS+=("$1"); shift
            fi ;;
    esac
done

# ============================================================================
# Validate
# ============================================================================

[[ -z "$OUTDIR_HOST" ]] && die "--outdir is required."
[[ -z "$WORKDIR_HOST" ]] && WORKDIR_HOST="${OUTDIR_HOST}/work"

mkdir -p "$OUTDIR_HOST" "$WORKDIR_HOST" "${OUTDIR_HOST}/pipeline_info" || die "Cannot create output directory: $OUTDIR_HOST"
[[ -w "$OUTDIR_HOST" ]] || die "Output directory is not writable: $OUTDIR_HOST"
for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
    [[ "${NF_ARGS[$i]}" == "--store_dir" ]] && mkdir -p "${NF_ARGS[$((i+1))]}"
done

# Auto-detect session
if [[ "$AUTO_SESSION" == true && -z "$RESUME_SESSION" ]]; then
    if [[ -f "${OUTDIR_HOST}/pipeline_info/run_command.sh" ]]; then
        # save_run_command writes "--session <uuid>"; `|| true` keeps set -e
        # happy when grep finds nothing.
        RESUME_SESSION=$(grep -oP -- '(--session|--?resume)\s+\K[0-9a-f-]{36}' \
            "${OUTDIR_HOST}/pipeline_info/run_command.sh" | tail -1 || true)
        if [[ -n "$RESUME_SESSION" ]]; then
            echo "[INFO] Auto-detected session from previous run: $RESUME_SESSION"
        fi
    fi
fi

# ============================================================================
# Container mode
# ============================================================================

if [[ "$USE_CONTAINER" == true ]]; then
    if [[ "$CONTAINER_RUNTIME" == "auto" ]]; then
        if command -v apptainer &>/dev/null; then CONTAINER_RUNTIME=apptainer
        elif command -v singularity &>/dev/null; then CONTAINER_RUNTIME=singularity
        elif command -v docker &>/dev/null; then CONTAINER_RUNTIME=docker
        else die "--container requires docker, apptainer, or singularity on PATH"; fi
        echo "[INFO] Auto-detected container runtime: ${CONTAINER_RUNTIME}"
    fi
    command -v "$CONTAINER_RUNTIME" &>/dev/null || die "${CONTAINER_RUNTIME} not found on PATH"

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

    # The pipeline code itself runs from this checkout (bound read-only), not
    # the copy baked into the image, so a branch can be run with a released image.
    BIND_PATHS+=("$SCRIPT_DIR" "$OUTDIR_HOST" "$WORKDIR_HOST")
    # The host TMPDIR is inherited by the container (Java dies without it).
    [[ -n "${TMPDIR:-}" && -d "${TMPDIR}" ]] && BIND_PATHS+=("$(realpath "$TMPDIR")")
    # Bind resolved targets too: inputs are often symlinks into other trees.
    ALL_BINDS=()
    for p in "${BIND_PATHS[@]}"; do
        ALL_BINDS+=("$p")
        r="$(realpath -m "$p")"; [[ "$r" != "$p" ]] && ALL_BINDS+=("$r")
        if [[ -L "$p" || -f "$p" ]]; then
            d="$(dirname "$(realpath -m "$p")")"; ALL_BINDS+=("$d")
        fi
    done
    # Drop duplicates and paths already covered by a bound parent.
    mapfile -t ALL_BINDS < <(printf '%s\n' "${ALL_BINDS[@]}" | awk 'NF' | sort -u)
    BINDS=()
    for p in "${ALL_BINDS[@]}"; do
        covered=false
        for q in "${ALL_BINDS[@]}"; do
            [[ "$q" != "$p" && "$p" == "$q"/* ]] && { covered=true; break; }
        done
        [[ -d "$p" || -f "$p" ]] && [[ "$covered" == false ]] && BINDS+=("$p")
    done

    NF_CACHE="${OUTDIR_HOST}/.nextflow-cache"
    mkdir -p "${NF_CACHE}/dotdir" 2>/dev/null || true

    CONTAINER_CMD=()
    case "$CONTAINER_RUNTIME" in
        docker)
            CONTAINER_CMD+=(docker run --rm --user "$(id -u):$(id -g)" -w "$OUTDIR_HOST")
            CONTAINER_CMD+=("-e" "NXF_HOME=/home/dana/.nextflow" "-v" "${NF_CACHE}/dotdir:/home/dana/.nextflow")
            for b in "${BINDS[@]}"; do CONTAINER_CMD+=("-v" "${b}:${b}"); done
            CONTAINER_CMD+=("$CONTAINER_IMAGE") ;;
        apptainer|singularity)
            container_ca="/etc/ssl/certs/ca-certificates.crt"
            CONTAINER_CMD+=("$CONTAINER_RUNTIME" run --pwd "$OUTDIR_HOST")
            CONTAINER_CMD+=("--env" "NXF_HOME=/home/dana/.nextflow")
            CONTAINER_CMD+=("--env" "REQUESTS_CA_BUNDLE=${container_ca}")
            CONTAINER_CMD+=("--env" "SSL_CERT_FILE=${container_ca}")
            CONTAINER_CMD+=("--env" "CURL_CA_BUNDLE=${container_ca}")
            CONTAINER_CMD+=("--bind" "${NF_CACHE}/dotdir:/home/dana/.nextflow")
            for b in "${BINDS[@]}"; do CONTAINER_CMD+=("--bind" "$b"); done
            CONTAINER_CMD+=("$SIF_PATH") ;;
    esac
    CONTAINER_CMD+=(run "${SCRIPT_DIR}/main.nf" -w "$WORKDIR_HOST" --outdir "$OUTDIR_HOST" "${NF_ARGS[@]}" -resume ${RESUME_SESSION})

    echo "[INFO] Mode:    Container (${CONTAINER_RUNTIME}, ${SIF_PATH:-$CONTAINER_IMAGE})"
    echo "[INFO] Output:  $OUTDIR_HOST"
    echo "[INFO] Binds:   ${BINDS[*]}"
    echo ""

    "${CONTAINER_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

    NF_SESSION=$(awk '{print $6}' "${NF_CACHE}/dotdir/history" 2>/dev/null | tail -1 || true)
    [[ "$NF_SESSION" =~ ^[0-9a-f-]{36}$ ]] || NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${OUTDIR_HOST}/.nextflow.log" 2>/dev/null | tail -1 || true)
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

cd "$OUTDIR_HOST"
LOCAL_CMD=(
    mamba run -p "${BBMAP_ENV}"
    nextflow run "${SCRIPT_DIR}/main.nf"
    -w "$WORKDIR_HOST"
    --outdir "$OUTDIR_HOST"
    "${NF_ARGS[@]}"
    -resume ${RESUME_SESSION}
)

echo "[INFO] Mode:    Local (conda)"
echo "[INFO] Output:  $OUTDIR_HOST"
echo ""

"${LOCAL_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${OUTDIR_HOST}/.nextflow.log" 2>/dev/null | tail -1 || true)
save_run_command "$OUTDIR_HOST" "$NF_SESSION"

exit $NF_EXIT
