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
#   # Kitchen sink — all options with defaults
#   ./run-mag.sh --input /data/reads --outdir /data/output \
#       --dedupe \
#       --filtlong_size 40000000000 \
#       --min_overlap 1000 \
#       --run_maxbin true \
#       --run_lorbin true \
#       --run_comebin true \
#       --lorbin_min_length 80000 \
#       --metabat_min_cls 50000 \
#       --assembly_cpus 24 \
#       --assembly_memory '64 GB'
#
#   # Docker mode
#   ./run-mag.sh --docker --input /data/reads --outdir /data/output
#
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_IMAGE="ghcr.io/rec3141/danaseq-mag:latest"
USE_CONTAINER=false
CONTAINER_RUNTIME=""   # docker, apptainer, or singularity
SIF_PATH=""            # resolved path to .sif file (apptainer/singularity only)
NF_ARGS=()
DB_ARGS=()
MOUNTS=()
ORIGINAL_ARGS=("$@")

die() { echo "[ERROR] $1" >&2; exit 1; }

# Resolve database paths from a standard base directory (as created by download-databases.sh).
# Populates DB_ARGS with --flag path pairs for every subdirectory/file that exists.
# These are prepended to NF_ARGS so explicit user flags override them.
resolve_db_dir() {
    local base="$1"
    local path

    # Simple fixed-layout databases
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
        [--metaeuk_db]="metaeuk_db/metaeuk_db"
    )
    for flag in "${!fixed_dbs[@]}"; do
        path="${base}/${fixed_dbs[$flag]}"
        [[ -e "$path" ]] && DB_ARGS+=("$flag" "$path")
    done

    # Kaiju: refseq_ref preferred, then first available subdir
    if [[ -d "${base}/kaiju/refseq_ref" ]]; then
        DB_ARGS+=(--kaiju_db "${base}/kaiju/refseq_ref")
    elif [[ -d "${base}/kaiju" ]]; then
        path=$(ls -d "${base}/kaiju"/*/2>/dev/null | head -1 || true)
        [[ -n "$path" ]] && DB_ARGS+=(--kaiju_db "$path")
    fi

    # Kraken2: pluspfp_08gb preferred, then flat layout (hash.k2d at root), then first subdir
    if [[ -d "${base}/krakendb/pluspfp_08gb" ]]; then
        DB_ARGS+=(--kraken2_db "${base}/krakendb/pluspfp_08gb")
    elif [[ -f "${base}/krakendb/hash.k2d" ]]; then
        DB_ARGS+=(--kraken2_db "${base}/krakendb")
    elif [[ -d "${base}/krakendb" ]]; then
        path=$(ls -d "${base}/krakendb"/*/2>/dev/null | head -1 || true)
        [[ -n "$path" ]] && DB_ARGS+=(--kraken2_db "$path")
    fi

    # SILVA: find SSU and LSU by glob (version-independent)
    path=$(ls "${base}/silva_db/SILVA_"*"_SSURef_NR99.fasta" 2>/dev/null | head -1 || true)
    [[ -n "$path" ]] && DB_ARGS+=(--silva_ssu_db "$path")
    path=$(ls "${base}/silva_db/SILVA_"*"_LSURef_NR99.fasta" 2>/dev/null | head -1 || true)
    [[ -n "$path" ]] && DB_ARGS+=(--silva_lsu_db "$path")
}

# Build a re-runnable self-invocation with the correct --session ID
save_run_command() {
    local outdir="$1" session="$2"
    local save_args=()
    local skip_next=false has_session=false
    for arg in "${ORIGINAL_ARGS[@]}"; do
        if $skip_next; then skip_next=false; continue; fi
        if [[ "$arg" == "--session" ]]; then skip_next=true; has_session=true; continue; fi
        save_args+=("$arg")
    done
    if [[ -n "$session" ]]; then
        save_args+=("--session" "$session")
    fi
    printf '%s\n' "$(realpath "$0") ${save_args[*]}" >> "${outdir}/pipeline_info/run_command.txt"
}

usage() {
    echo "Usage: $0 --input DIR --outdir DIR [pipeline options]"
    echo ""
    echo "Required:"
    echo "  --input DIR      Directory containing reads (*.fastq.gz or barcode structure)"
    echo "  --outdir DIR     Output directory (will be created if needed)"
    echo ""
    echo "Mode:"
    echo "  --docker         Run in Docker container"
    echo "  --apptainer      Run in Apptainer/Singularity container"
    echo "  --container      Auto-detect container runtime (apptainer > singularity > docker)"
    echo "  --image IMAGE    Override container image [default: ghcr.io/rec3141/danaseq-mag:latest]"
    echo "  --sif PATH       Use a specific .sif file (apptainer/singularity only)"
    echo ""
    echo "Shortcuts:"
    echo "  --all            Enable all optional modules (dedupe, bakta_extra, kraken2,"
    echo "                   sendsketch, rrna, metabolism, eukaryotic, metaeuk, marferret)"
    echo "                   Still requires database paths to be provided separately."
    echo "  --db_dir DIR     Base directory of databases installed by download-databases.sh."
    echo "                   Auto-resolves all standard subdirectory paths (bakta, genomad,"
    echo "                   checkv, checkm2, kaiju, kraken2, silva, kofam, eggnog, dbcan,"
    echo "                   macsyfinder, defensefinder, metaeuk, marferret). Explicit"
    echo "                   --flag PATH overrides any auto-detected path."
    echo ""
    echo "Caching:"
    echo "  --workdir DIR        Nextflow work directory [-w] (default: /tmp/nanopore_mag_work)"
    echo "  --store_dir DIR      Persistent cache directory (storeDir); skips completed processes"
    echo "                       across runs even after work/ cleanup. Off by default."
    echo ""
    echo "Pipeline flags (passed to Nextflow):"
    echo "  --dedupe             BBDuk deduplication before assembly"
    echo "  --filtlong_size N    Filtlong target bases (e.g. 40000000000); skip if not set"
    echo "  --min_overlap N      Flye --min-overlap [default: 1000]"
    echo "  --run_maxbin BOOL    Include MaxBin2 in consensus [default: true]"
    echo "  --run_lorbin BOOL    Include LorBin in consensus [default: true]"
    echo "  --run_comebin BOOL   Include COMEBin in consensus [default: true]"
    echo "  --lorbin_min_length N LorBin minimum bin size in bp [default: 80000]"
    echo "  --metabat_min_cls N  MetaBAT2 minimum cluster size [default: 50000]"
    echo "  --assembly_cpus N    CPUs for assembly [default: 24]"
    echo "  --assembly_memory S  Memory for assembly [default: '64 GB']"
    echo ""
    echo "Kitchen sink — compact form (this system):"
    echo "  $0 --input /data/minknow/QEI2025 \\"
    echo "      --outdir /data/minknow/QEI2025/nanopore_mag/tmpdir \\"
    echo "      --store_dir /data/minknow/QEI2025/nanopore_mag/final \\"
    echo "      --workdir /data/scratch/work \\"
    echo "      --all \\"
    echo "      --db_dir /data/scratch/refdbs \\"
    echo "      --sendsketch_address http://10.151.50.41:3068/sketch \\"
    echo "      --filtlong_size 40000000000 \\"
    echo "      --annotator bakta \\"
    echo "      --assembly_cpus 24 \\"
    echo "      --assembly_memory '120 GB'"
    echo ""
    echo "Kitchen sink — explicit form (same command, all flags spelled out):"
    echo "  $0 --input /data/minknow/QEI2025 \\"
    echo "      --outdir /data/minknow/QEI2025/nanopore_mag/tmpdir \\"
    echo "      --store_dir /data/minknow/QEI2025/nanopore_mag/final \\"
    echo "      --workdir /data/scratch/work \\"
    echo "      --dedupe \\"
    echo "      --filtlong_size 40000000000 \\"
    echo "      --annotator bakta \\"
    echo "      --bakta_light_db /data/scratch/refdbs/bakta/db-light \\"
    echo "      --bakta_db /data/scratch/refdbs/bakta/db \\"
    echo "      --bakta_extra \\"
    echo "      --genomad_db /data/scratch/refdbs/genomad_db \\"
    echo "      --checkv_db /data/scratch/refdbs/checkv_db \\"
    echo "      --checkm2_db /data/scratch/refdbs/checkm2 \\"
    echo "      --kaiju_db /data/scratch/refdbs/kaiju/refseq_ref \\"
    echo "      --run_kraken2 true \\"
    echo "      --kraken2_db /data/scratch/refdbs/krakendb/pluspfp_08gb \\"
    echo "      --run_sendsketch true \\"
    echo "      --sendsketch_address http://10.151.50.41:3068/sketch \\"
    echo "      --run_rrna true \\"
    echo "      --silva_ssu_db /data/scratch/refdbs/silva_db/SILVA_138.2_SSURef_NR99.fasta \\"
    echo "      --silva_lsu_db /data/scratch/refdbs/silva_db/SILVA_138.2_LSURef_NR99.fasta \\"
    echo "      --run_metabolism true \\"
    echo "      --kofam_db /data/scratch/refdbs/kofam_db \\"
    echo "      --eggnog_db /data/scratch/refdbs/eggnog_db \\"
    echo "      --dbcan_db /data/scratch/refdbs/dbcan_db \\"
    echo "      --macsyfinder_models /data/scratch/refdbs/macsyfinder_models \\"
    echo "      --defensefinder_models /data/scratch/refdbs/defensefinder_models \\"
    echo "      --run_eukaryotic true \\"
    echo "      --run_metaeuk true \\"
    echo "      --metaeuk_db /data/scratch/refdbs/metaeuk_db/metaeuk_db \\"
    echo "      --run_marferret true \\"
    echo "      --marferret_db /data/scratch/refdbs/marferret_db \\"
    echo "      --assembly_cpus 24 \\"
    echo "      --assembly_memory '120 GB'"
    echo ""
    echo "Run '$0 --help-pipeline' to see full Nextflow help."
    exit 0
}

# ============================================================================
# Parse arguments -- extract paths that need mounting, pass rest to Nextflow
# ============================================================================

INPUT_HOST=""
OUTDIR_HOST=""
WORKDIR_HOST="/tmp/nanopore_mag_work"
RESUME_SESSION=""
AUTO_SESSION=true
DB_DIR_HOST=""
STORE_DIR_HOST=""

while (( $# )); do
    case "$1" in
        -h|--help)
            usage ;;
        --help-pipeline)
            if [[ "$USE_CONTAINER" == true ]]; then
                # help-pipeline before runtime detection — just use docker if available
                docker run --user "$(id -u):$(id -g)" "$CONTAINER_IMAGE" \
                    run /pipeline/main.nf --help
            else
                mamba run -p "${SCRIPT_DIR}/conda-envs/dana-mag-flye" \
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
                --dedupe
                --bakta_extra
                --run_kraken2 true
                --run_sendsketch true
                --run_rrna true
                --run_metabolism true
                --run_eukaryotic true
                --run_metaeuk true
                --run_marferret true
            )
            shift ;;
        --session)
            [[ -z "${2:-}" ]] && die "--session requires a session ID"
            RESUME_SESSION="$2"
            AUTO_SESSION=false
            shift 2 ;;
        *)
            NF_ARGS+=("$1")
            shift ;;
    esac
done

# Prepend DB_ARGS so explicit user flags (later in NF_ARGS) override auto-detected paths
NF_ARGS=("${DB_ARGS[@]}" "${NF_ARGS[@]}")

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
    # Extract session ID from last run_command.sh entry (UUID after -resume)
    if [[ -f "${OUTDIR_HOST}/pipeline_info/run_command.txt" ]]; then
        RESUME_SESSION=$(grep -oP '(?<=--session )[0-9a-f-]{36}' \
            "${OUTDIR_HOST}/pipeline_info/run_command.txt" | tail -1 || true)
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
            SIF_PATH="${SCRIPT_DIR}/.danaseq-mag.sif"
            if [[ ! -f "$SIF_PATH" ]]; then
                echo "[INFO] Pulling container image (one-time download)..."
                echo "  Source: docker://${CONTAINER_IMAGE}"
                echo "  Destination: ${SIF_PATH}"
                "$CONTAINER_RUNTIME" pull "$SIF_PATH" "docker://${CONTAINER_IMAGE}"
            fi
        fi
        [[ -f "$SIF_PATH" ]] || die "SIF file not found: $SIF_PATH"
    fi
fi

# ============================================================================
# Container mode
# ============================================================================

if [[ "$USE_CONTAINER" == true ]]; then
    # Collect bind mounts as "host:container[:opts]" pairs (runtime-agnostic)
    BINDS=()
    BINDS+=("${INPUT_HOST}:/data/input:ro")
    BINDS+=("${OUTDIR_HOST}:/data/output")
    NF_ARGS=("--input" "/data/input" "--outdir" "/data/output" "${NF_ARGS[@]}")

    # Mount database directory and rewrite all host paths in NF_ARGS to container paths
    if [[ -n "${DB_DIR_HOST:-}" ]]; then
        BINDS+=("${DB_DIR_HOST}:/data/db:ro")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${DB_DIR_HOST}/\/data\/db}"
        done
    fi

    # Mount store directory (read-write) and rewrite paths
    if [[ -n "${STORE_DIR_HOST:-}" ]]; then
        BINDS+=("${STORE_DIR_HOST}:/data/store")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${STORE_DIR_HOST}/\/data\/store}"
        done
    fi

    # Persistent Nextflow directories for -resume support
    NF_CACHE="${OUTDIR_HOST}/.nextflow-cache"
    mkdir -p "${NF_CACHE}/work" "${NF_CACHE}/dotdir" 2>/dev/null || true
    BINDS+=("${NF_CACHE}/work:/home/dana/work")
    BINDS+=("${NF_CACHE}/dotdir:/home/dana/.nextflow")

    # Build runtime-specific command
    CONTAINER_CMD=()
    case "$CONTAINER_RUNTIME" in
        docker)
            CONTAINER_CMD+=(docker run --user "$(id -u):$(id -g)")
            for bind in "${BINDS[@]}"; do
                CONTAINER_CMD+=("-v" "$bind")
            done
            CONTAINER_CMD+=("$CONTAINER_IMAGE" run /pipeline/main.nf)
            ;;
        apptainer|singularity)
            # Override host SSL env vars (RHEL CA paths don't exist in the Debian container)
            container_ca="/etc/ssl/certs/ca-certificates.crt"
            CONTAINER_CMD+=("$CONTAINER_RUNTIME" run)
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

    # Capture session ID from Nextflow log and save self-invocation for reliable resume
    NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${SCRIPT_DIR}/.nextflow.log" 2>/dev/null | tail -1)
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

WORKDIR_FLAG=()
if [[ -n "$WORKDIR_HOST" ]]; then
    mkdir -p "$WORKDIR_HOST" || die "Cannot create work directory: $WORKDIR_HOST"
    WORKDIR_FLAG=(-w "$WORKDIR_HOST")
fi

LOCAL_CMD=(
    mamba run -p "${SCRIPT_DIR}/conda-envs/dana-mag-flye"
    nextflow run "${SCRIPT_DIR}/main.nf"
    --input "$INPUT_HOST"
    --outdir "$OUTDIR_HOST"
    "${WORKDIR_FLAG[@]}"
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
NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${SCRIPT_DIR}/.nextflow.log" 2>/dev/null | tail -1)
save_run_command "$OUTDIR_HOST" "$NF_SESSION"

exit $NF_EXIT
