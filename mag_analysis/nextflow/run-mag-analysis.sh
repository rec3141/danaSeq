#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana MAG Analysis Pipeline Launcher
# ============================================================================
#
# Runs the technology-agnostic downstream analysis pipeline on pre-computed
# assembly outputs. Accepts assembly FASTA + depth table from any assembler
# (nanopore_assembly, illumina_assembly, or external).
#
# Usage:
#   ./run-mag-analysis.sh --assembly /path/to/assembly.fasta \
#       --depths /path/to/depths.txt --outdir /path/to/output [options]
#
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_IMAGE="ghcr.io/rec3141/danaseq-mag-analysis:latest"
USE_CONTAINER=false
CONTAINER_RUNTIME=""
SIF_PATH=""
PULL_SIF=false
NF_ARGS=()
DB_ARGS=()
MOUNTS=()
ORIGINAL_ARGS=("$@")

die() { echo "[ERROR] $1" >&2; exit 1; }

# Resolve database paths from a standard base directory
resolve_db_dir() {
    local base="$1"
    local path

    local -A fixed_dbs=(
        [--bakta_db]="bakta/db"
        [--bakta_light_db]="bakta/db-light"
        [--genomad_db]="genomad_db"
        [--checkv_db]="checkv_db"
        [--checkm2_db]="checkm2"
        [--kofam_db]="kofam_db"
        [--eggnog_db]="eggnog_db"
        [--dbcan_db]="dbcan_db"
        [--macsyfinder_models]="macsyfinder_models"
        [--defensefinder_models]="defensefinder_models"
        [--marferret_db]="marferret_db"
        [--gtdbtk_db]="gtdbtk_db"
        [--metaeuk_db]="metaeuk_db/metaeuk_db"
        [--antismash_db]="antismash_db"
        [--magscot_hmm_dir]="magscot_hmm"
    )
    for flag in "${!fixed_dbs[@]}"; do
        path="${base}/${fixed_dbs[$flag]}"
        [[ -e "$path" ]] && DB_ARGS+=("$flag" "$path")
    done

    # Kaiju
    if [[ -d "${base}/kaiju/refseq_ref" ]]; then
        DB_ARGS+=(--kaiju_db "${base}/kaiju/refseq_ref")
    elif [[ -d "${base}/kaiju" ]]; then
        path=$(ls -d "${base}/kaiju"/*/2>/dev/null | head -1 || true)
        [[ -n "$path" ]] && DB_ARGS+=(--kaiju_db "$path")
    fi

    # Kraken2
    if [[ -d "${base}/krakendb/pluspfp_08gb" ]]; then
        DB_ARGS+=(--kraken2_db "${base}/krakendb/pluspfp_08gb")
    elif [[ -f "${base}/krakendb/hash.k2d" ]]; then
        DB_ARGS+=(--kraken2_db "${base}/krakendb")
    elif [[ -d "${base}/krakendb" ]]; then
        path=$(ls -d "${base}/krakendb"/*/2>/dev/null | head -1 || true)
        [[ -n "$path" ]] && DB_ARGS+=(--kraken2_db "$path")
    fi

    # SILVA
    path=$(ls "${base}/silva_db/SILVA_"*"_SSURef_NR99.fasta" 2>/dev/null | head -1 || true)
    [[ -n "$path" ]] && DB_ARGS+=(--silva_ssu_db "$path")
    path=$(ls "${base}/silva_db/SILVA_"*"_LSURef_NR99.fasta" 2>/dev/null | head -1 || true)
    [[ -n "$path" ]] && DB_ARGS+=(--silva_lsu_db "$path")
}

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
    echo "Usage: $0 --assembly FILE --depths FILE --outdir DIR [options]"
    echo ""
    echo "Required:"
    echo "  --assembly FILE    Path to assembly FASTA"
    echo "  --depths FILE      Path to MetaBAT2-format depth table"
    echo "  --outdir DIR       Output directory"
    echo ""
    echo "Optional inputs:"
    echo "  --bam_dir DIR      Directory with *.sorted.bam + .bai (for SemiBin2, LorBin, COMEBin)"
    echo "  --tnf FILE         Tetranucleotide frequency table (for viz)"
    echo ""
    echo "Mode:"
    echo "  --docker           Run in Docker container"
    echo "  --apptainer        Run in Apptainer/Singularity container"
    echo "  --container        Auto-detect container runtime"
    echo ""
    echo "Shortcuts:"
    echo "  --all              Enable all optional modules (still requires DB paths)"
    echo "  --db_dir DIR       Auto-resolve database paths from download-databases.sh layout"
    echo ""
    echo "Caching & Resume:"
    echo "  --workdir DIR      Nextflow work directory (default: /tmp/mag_analysis_work)"
    echo "  --store_dir DIR    Persistent cache directory (storeDir)"
    echo "  --resume [ID]      Resume a previous run"
    echo ""
    echo "Pipeline flags: --annotator, --run_metabolism, --run_genomad, --run_gtdbtk,"
    echo "  --run_kaiju, --run_kraken2, --run_sendsketch, --run_rrna, --run_eukaryotic,"
    echo "  --run_antismash, --run_viz, --run_semibin, --run_maxbin, --run_lorbin,"
    echo "  --run_comebin, --run_vamb, --run_binette, --run_magscot, etc."
    echo ""
    echo "Run '$0 --help-pipeline' for full Nextflow help."
    exit 0
}

# ============================================================================
# Parse arguments
# ============================================================================

ASSEMBLY_HOST=""
DEPTHS_HOST=""
BAM_DIR_HOST=""
OUTDIR_HOST=""
WORKDIR_HOST="/tmp/mag_analysis_work"
RESUME_SESSION=""
DO_RESUME=false
DB_DIR_HOST=""
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
            [[ -z "${2:-}" ]] && die "--sif requires a path"
            SIF_PATH="$2"
            shift 2 ;;
        --pull)
            PULL_SIF=true
            shift ;;
        --assembly)
            [[ -z "${2:-}" ]] && die "--assembly requires a file path"
            ASSEMBLY_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            NF_ARGS+=("--assembly" "$ASSEMBLY_HOST")
            shift 2 ;;
        --depths)
            [[ -z "${2:-}" ]] && die "--depths requires a file path"
            DEPTHS_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            NF_ARGS+=("--depths" "$DEPTHS_HOST")
            shift 2 ;;
        --bam_dir)
            [[ -z "${2:-}" ]] && die "--bam_dir requires a directory path"
            BAM_DIR_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            NF_ARGS+=("--bam_dir" "$BAM_DIR_HOST")
            shift 2 ;;
        --outdir)
            [[ -z "${2:-}" ]] && die "--outdir requires a directory path"
            OUTDIR_HOST="$(realpath -m "$2")"
            shift 2 ;;
        --workdir|-w)
            [[ -z "${2:-}" ]] && die "--workdir requires a directory path"
            WORKDIR_HOST="$(realpath -m "$2")"
            shift 2 ;;
        --db_dir)
            [[ -z "${2:-}" ]] && die "--db_dir requires a directory path"
            [[ -d "$2" ]] || die "--db_dir directory does not exist: $2"
            DB_DIR_HOST="$(realpath "$2")"
            resolve_db_dir "$DB_DIR_HOST"
            shift 2 ;;
        --store_dir)
            [[ -z "${2:-}" ]] && die "--store_dir requires a directory path"
            STORE_DIR_HOST="$(realpath -m "$2")"
            NF_ARGS+=("--store_dir" "$STORE_DIR_HOST")
            shift 2 ;;
        --all)
            NF_ARGS+=(
                --bakta_extra
                --run_semibin true
                --run_maxbin true
                --run_lorbin true
                --run_comebin true
                --run_vamb true
                --run_vamb_tax true
                --run_binette true
                --run_magscot true
                --run_kraken2 true
                --run_sendsketch true
                --run_rrna true
                --run_metabolism true
                --run_antismash true
                --run_eukaryotic true
                --run_metaeuk true
                --run_marferret true
                --run_gtdbtk true
                --run_kaiju true
                --run_genomad true
                --run_checkv true
                --run_defensefinder true
                --run_integronfinder true
                --run_islandpath true
                --run_macsyfinder true
                --run_viz true
            )
            shift ;;
        --resume)
            DO_RESUME=true
            if [[ -n "${2:-}" && "${2}" != --* ]]; then
                RESUME_SESSION="$2"
                shift
            fi
            shift ;;
        *)
            NF_ARGS+=("$1")
            shift ;;
    esac
done

# Prepend DB_ARGS so explicit user flags override auto-detected paths
NF_ARGS=("${DB_ARGS[@]}" "${NF_ARGS[@]}")

# ============================================================================
# Validate
# ============================================================================

[[ -z "$ASSEMBLY_HOST" ]] && die "--assembly is required."
[[ -z "$DEPTHS_HOST" ]] && die "--depths is required."

if [[ -n "$OUTDIR_HOST" && -n "$STORE_DIR_HOST" ]]; then
    die "--outdir and --store_dir are mutually exclusive."
fi

if [[ -z "$OUTDIR_HOST" && -n "$STORE_DIR_HOST" ]]; then
    OUTDIR_HOST="$STORE_DIR_HOST"
    NF_ARGS+=("--outdir" "$OUTDIR_HOST")
fi
[[ -z "$OUTDIR_HOST" ]] && die "--outdir (or --store_dir) is required."

[[ -f "$ASSEMBLY_HOST" ]] || die "Assembly file not found: $ASSEMBLY_HOST"
[[ -f "$DEPTHS_HOST" ]] || die "Depths file not found: $DEPTHS_HOST"

if [[ ! -d "$OUTDIR_HOST" ]]; then
    echo "[INFO] Creating output directory: $OUTDIR_HOST"
    mkdir -p "$OUTDIR_HOST" || die "Cannot create output directory: $OUTDIR_HOST"
fi

# ============================================================================
# Auto-detect session
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
        echo "[INFO] Auto-detected session: $RESUME_SESSION"
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
        else die "No container runtime found"; fi
    fi

    if [[ "$CONTAINER_RUNTIME" == "apptainer" || "$CONTAINER_RUNTIME" == "singularity" ]]; then
        [[ -z "$SIF_PATH" ]] && SIF_PATH="${SCRIPT_DIR}/.danaseq-mag-analysis.sif"
        if [[ ! -f "$SIF_PATH" ]]; then
            if [[ "$PULL_SIF" == true ]]; then
                "$CONTAINER_RUNTIME" pull "$SIF_PATH" "docker://${CONTAINER_IMAGE}"
            else
                die "SIF not found: $SIF_PATH"
            fi
        fi
    fi
fi

# ============================================================================
# Container mode
# ============================================================================

if [[ "$USE_CONTAINER" == true ]]; then
    BINDS=()
    BINDS+=("$(dirname "$ASSEMBLY_HOST"):/data/assembly:ro")
    BINDS+=("$(dirname "$DEPTHS_HOST"):/data/depths:ro")
    BINDS+=("${OUTDIR_HOST}:/data/output")

    # Rewrite assembly/depths paths
    for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
        NF_ARGS[$i]="${NF_ARGS[$i]//${ASSEMBLY_HOST}/\/data\/assembly\/$(basename "$ASSEMBLY_HOST")}"
        NF_ARGS[$i]="${NF_ARGS[$i]//${DEPTHS_HOST}/\/data\/depths\/$(basename "$DEPTHS_HOST")}"
    done
    NF_ARGS=("--outdir" "/data/output" "${NF_ARGS[@]}")

    if [[ -n "${BAM_DIR_HOST:-}" ]]; then
        BINDS+=("${BAM_DIR_HOST}:/data/bams:ro")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${BAM_DIR_HOST}/\/data\/bams}"
        done
    fi

    if [[ -n "${DB_DIR_HOST:-}" ]]; then
        BINDS+=("${DB_DIR_HOST}:/data/db:ro")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${DB_DIR_HOST}/\/data\/db}"
        done
    fi

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
    BINDS+=("${NF_CACHE}/dotdir:/home/dana/.nextflow")

    CONTAINER_CMD=()
    case "$CONTAINER_RUNTIME" in
        docker)
            CONTAINER_CMD+=(docker run --user "$(id -u):$(id -g)")
            CONTAINER_CMD+=("-e" "NXF_HOME=/home/dana/.nextflow")
            for bind in "${BINDS[@]}"; do CONTAINER_CMD+=("-v" "$bind"); done
            CONTAINER_CMD+=("$CONTAINER_IMAGE" run /pipeline/main.nf)
            ;;
        apptainer|singularity)
            container_ca="/etc/ssl/certs/ca-certificates.crt"
            CONTAINER_CMD+=("$CONTAINER_RUNTIME" run)
            CONTAINER_CMD+=("--env" "NXF_HOME=${NF_CACHE}/dotdir")
            CONTAINER_CMD+=("--env" "REQUESTS_CA_BUNDLE=${container_ca}")
            CONTAINER_CMD+=("--env" "SSL_CERT_FILE=${container_ca}")
            CONTAINER_CMD+=("--env" "CURL_CA_BUNDLE=${container_ca}")
            for bind in "${BINDS[@]}"; do CONTAINER_CMD+=("--bind" "$bind"); done
            CONTAINER_CMD+=("$SIF_PATH" run /pipeline/main.nf)
            ;;
    esac
    CONTAINER_CMD+=(-w /data/work "${NF_ARGS[@]}")
    if [[ "$DO_RESUME" == true ]]; then
        CONTAINER_CMD+=(-resume ${RESUME_SESSION})
    fi

    echo "[INFO] Mode:   Container (${CONTAINER_RUNTIME})"
    echo "[INFO] Assembly: $ASSEMBLY_HOST"
    echo "[INFO] Depths:   $DEPTHS_HOST"
    echo "[INFO] Output:   $OUTDIR_HOST"
    echo ""

    mkdir -p "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info" 2>/dev/null || true
    "${CONTAINER_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

    NF_SESSION=$(awk '{print $6}' "${NF_CACHE}/dotdir/history" 2>/dev/null | tail -1)
    save_run_command "${STORE_DIR_HOST:-$OUTDIR_HOST}" "$NF_SESSION"
    exit $NF_EXIT
fi

# ============================================================================
# Local mode (default)
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

# Use the flye env (has Nextflow + Java) or any available env with Nextflow
NF_ENV="${SCRIPT_DIR}/conda-envs/dana-mag-assembly"
if [[ ! -d "$NF_ENV" ]]; then
    # Fallback: look for any env with nextflow
    for env in "${SCRIPT_DIR}"/conda-envs/dana-mag-*; do
        if [[ -x "${env}/bin/nextflow" ]]; then
            NF_ENV="$env"
            break
        fi
    done
fi

LOCAL_CMD=(
    mamba run -p "$NF_ENV"
    nextflow run "${SCRIPT_DIR}/main.nf"
    --outdir "$OUTDIR_HOST"
    "${WORKDIR_FLAG[@]}"
    "${NF_ARGS[@]}"
    "${RESUME_FLAG[@]}"
)

echo "[INFO] Mode:     Local (conda)"
echo "[INFO] Assembly: $ASSEMBLY_HOST"
echo "[INFO] Depths:   $DEPTHS_HOST"
echo "[INFO] Output:   $OUTDIR_HOST"
echo ""

mkdir -p "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info" 2>/dev/null || true

"${LOCAL_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

NF_SESSION=$(awk '{print $6}' .nextflow/history 2>/dev/null | tail -1)
[[ -z "$NF_SESSION" ]] && NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' .nextflow.log 2>/dev/null | tail -1)
save_run_command "${STORE_DIR_HOST:-$OUTDIR_HOST}" "$NF_SESSION"

exit $NF_EXIT
