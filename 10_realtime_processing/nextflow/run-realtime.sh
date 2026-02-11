#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana Real-Time Pipeline Launcher
# ============================================================================
#
# Runs the Nextflow pipeline locally (default) or via Docker (--docker).
#
# Usage:
#   ./run-realtime.sh --input /path/to/data --outdir /path/to/output [options]
#
# Examples:
#   # Basic QC (local conda)
#   ./run-realtime.sh --input /data/run1 --outdir /data/output
#
#   # Full pipeline with Kraken2
#   ./run-realtime.sh --input /data/run1 --outdir /data/output \
#       --run_kraken --kraken_db /path/to/krakendb \
#       --run_prokka --run_sketch --run_tetra
#
#   # With HMM profiling
#   ./run-realtime.sh --input /data/run1 --outdir /data/output \
#       --run_prokka --hmm_databases /path/to/CANT-HYD.hmm
#
#   # Multiple HMM databases
#   ./run-realtime.sh --input /data/run1 --outdir /data/output \
#       --run_prokka --hmm_databases /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm
#
#   # With DuckDB integration
#   ./run-realtime.sh --input /data/run1 --outdir /data/output \
#       --run_kraken --kraken_db /path/to/krakendb \
#       --run_prokka --run_sketch --run_tetra \
#       --run_db_integration
#
#   # Kitchen sink — all modules, all options with defaults
#   ./run-realtime.sh --input /data/run1 --outdir /data/output \
#       --run_kraken --kraken_db /path/to/krakendb \
#       --run_prokka \
#       --hmm_databases /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm \
#       --run_sketch \
#       --run_tetra \
#       --run_db_integration \
#       --cleanup \
#       --min_readlen 1500 \
#       --keep_percent 80 \
#       --min_file_size 1000000
#
#   # Kitchen sink — watch mode for live sequencing
#   ./run-realtime.sh --input /data/runs --outdir /data/output \
#       --watch --db_sync_minutes 10 \
#       --run_kraken --kraken_db /path/to/krakendb \
#       --run_prokka \
#       --hmm_databases /path/to/CANT-HYD.hmm \
#       --run_sketch \
#       --run_tetra \
#       --run_db_integration
#
#   # Docker mode
#   ./run-realtime.sh --docker --input /data/run1 --outdir /data/output \
#       --run_kraken --kraken_db /path/to/krakendb \
#       --run_prokka
#
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="danaseq-realtime"
USE_DOCKER=false
DOCKER_ARGS=()
NF_ARGS=()
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
    echo "Mode:"
    echo "  --docker             Run in Docker instead of local conda"
    echo ""
    echo "Optional:"
    echo "  --kraken_db DIR      Kraken2 database directory"
    echo "  --hmm_databases LIST Comma-separated HMM file paths"
    echo ""
    echo "Pipeline flags (passed to Nextflow):"
    echo "  --run_kraken         Kraken2 taxonomic classification (requires --kraken_db)"
    echo "  --run_prokka         Prokka gene annotation"
    echo "  --run_sketch         Sendsketch taxonomic profiling"
    echo "  --run_tetra          Tetranucleotide frequency analysis"
    echo "  --run_db_integration Load results into DuckDB"
    echo "  --cleanup            Compress/delete source files after DuckDB import"
    echo "  --watch              Monitor for new FASTQ files (live sequencing)"
    echo "  --db_sync_minutes N  DB sync interval in watch mode [default: 10]"
    echo "  --min_readlen N      Minimum read length in bp [default: 1500]"
    echo "  --keep_percent N     Filtlong keep percent [default: 80]"
    echo "  --min_file_size N    Minimum FASTQ size in bytes [default: 1000000]"
    echo ""
    echo "Kitchen sink example (all modules, all options with defaults):"
    echo "  $0 --input /data/run1 --outdir /data/output \\"
    echo "      --run_kraken --kraken_db /path/to/krakendb \\"
    echo "      --run_prokka \\"
    echo "      --hmm_databases /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm \\"
    echo "      --run_sketch \\"
    echo "      --run_tetra \\"
    echo "      --run_db_integration \\"
    echo "      --cleanup \\"
    echo "      --min_readlen 1500 \\"
    echo "      --keep_percent 80 \\"
    echo "      --min_file_size 1000000"
    echo ""
    echo "Run '$0 --help-pipeline' to see full Nextflow help."
    exit 0
}

# ============================================================================
# Parse arguments — extract paths that need mounting, pass rest to Nextflow
# ============================================================================

INPUT_HOST=""
OUTDIR_HOST=""
KRAKEN_DB_HOST=""
HMM_HOST=""

while (( $# )); do
    case "$1" in
        -h|--help)
            usage ;;
        --help-pipeline)
            if [[ "$USE_DOCKER" == true ]]; then
                docker run --user "$(id -u):$(id -g)" "$IMAGE" \
                    run /pipeline/main.nf --help
            else
                mamba run -p "${SCRIPT_DIR}/conda-envs/dana-tools" \
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
        --kraken_db)
            [[ -z "${2:-}" ]] && die "--kraken_db requires a directory path"
            KRAKEN_DB_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            shift 2 ;;
        --danadir)
            # Accepted but ignored — auto-detected from script location
            warn "--danadir is deprecated and auto-detected; ignoring"
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

    # Kraken2 database
    if [[ -n "$KRAKEN_DB_HOST" ]]; then
        [[ -d "$KRAKEN_DB_HOST" ]] || die "Kraken2 database directory does not exist: $KRAKEN_DB_HOST"
        if ! ls "$KRAKEN_DB_HOST"/*.k2d &>/dev/null; then
            warn "No .k2d files found in $KRAKEN_DB_HOST — are you sure this is a Kraken2 database?"
        fi
        MOUNTS+=("-v" "${KRAKEN_DB_HOST}:/kraken_db:ro")
        NF_ARGS+=("--kraken_db" "/kraken_db")
    fi

    # DuckDB integration — auto-mount bin/ directory
    for arg in "${NF_ARGS[@]}"; do
        if [[ "$arg" == "--run_db_integration" ]]; then
            MOUNTS+=("-v" "${SCRIPT_DIR}/bin:/danadir:ro")
            NF_ARGS+=("--danadir" "/danadir")
            break
        fi
    done

    # HMM databases
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
        hmm_joined="$(IFS=','; echo "${container_hmm_paths[*]}")"
        NF_ARGS+=("--hmm_databases" "$hmm_joined")
        echo "[INFO] Mounting $hmm_counter HMM database(s)"
    fi

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
    [[ -n "$KRAKEN_DB_HOST" ]] && echo "[INFO] Kraken DB: $KRAKEN_DB_HOST"
    echo "[INFO] Running: ${DOCKER_CMD[*]}"
    echo ""

    exec "${DOCKER_CMD[@]}"
fi

# ============================================================================
# Local mode (default) — run via conda
# ============================================================================

# Auto-set --danadir for DB integration
for arg in "${NF_ARGS[@]}"; do
    if [[ "$arg" == "--run_db_integration" ]]; then
        NF_ARGS+=("--danadir" "${SCRIPT_DIR}/bin")
        break
    fi
done

# Kraken2 database — validate in local mode
if [[ -n "$KRAKEN_DB_HOST" ]]; then
    [[ -d "$KRAKEN_DB_HOST" ]] || die "Kraken2 database directory does not exist: $KRAKEN_DB_HOST"
    if ! ls "$KRAKEN_DB_HOST"/*.k2d &>/dev/null; then
        warn "No .k2d files found in $KRAKEN_DB_HOST — are you sure this is a Kraken2 database?"
    fi
    NF_ARGS+=("--kraken_db" "$KRAKEN_DB_HOST")
fi

# HMM databases — validate paths exist
if [[ -n "$HMM_HOST" ]]; then
    IFS=',' read -ra HMM_FILES <<< "$HMM_HOST"
    for hmm_file in "${HMM_FILES[@]}"; do
        hmm_file="$(echo "$hmm_file" | xargs)"
        hmm_resolved="$(realpath "$hmm_file" 2>/dev/null || echo "$hmm_file")"
        [[ -f "$hmm_resolved" ]] || die "HMM file does not exist: $hmm_file"
    done
    NF_ARGS+=("--hmm_databases" "$HMM_HOST")
fi

LOCAL_CMD=(
    mamba run -p "${SCRIPT_DIR}/conda-envs/dana-tools"
    nextflow run "${SCRIPT_DIR}/main.nf"
    --input "$INPUT_HOST"
    --outdir "$OUTDIR_HOST"
    "${NF_ARGS[@]}"
    -resume
)

echo "[INFO] Mode:   Local (conda)"
echo "[INFO] Input:  $INPUT_HOST"
echo "[INFO] Output: $OUTDIR_HOST"
[[ -n "$KRAKEN_DB_HOST" ]] && echo "[INFO] Kraken DB: $KRAKEN_DB_HOST"
echo "[INFO] Running: ${LOCAL_CMD[*]}"
echo ""

exec "${LOCAL_CMD[@]}"
