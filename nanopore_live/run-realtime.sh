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
#       --annotator bakta --run_sketch --run_tetra
#
#   # With HMM profiling
#   ./run-realtime.sh --input /data/run1 --outdir /data/output \
#       --annotator bakta --hmm_databases /path/to/CANT-HYD.hmm
#
#   # Multiple HMM databases
#   ./run-realtime.sh --input /data/run1 --outdir /data/output \
#       --annotator bakta --hmm_databases /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm
#
#   # With DuckDB integration
#   ./run-realtime.sh --input /data/run1 --outdir /data/output \
#       --run_kraken --kraken_db /path/to/krakendb \
#       --annotator bakta --run_sketch --run_tetra \
#       --run_db_integration
#
#   # Kitchen sink — all modules, all options with defaults
#   ./run-realtime.sh --input /data/run1 --outdir /data/output \
#       --run_kraken --kraken_db /path/to/krakendb \
#       --annotator bakta \
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
#       --annotator bakta \
#       --hmm_databases /path/to/CANT-HYD.hmm \
#       --run_sketch \
#       --run_tetra \
#       --run_db_integration
#
#   # Docker mode
#   ./run-realtime.sh --docker --input /data/run1 --outdir /data/output \
#       --run_kraken --kraken_db /path/to/krakendb \
#       --annotator bakta
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
    echo "  --session UUID       Resume a specific Nextflow session ID (otherwise"
    echo "                       auto-detected from <outdir>/pipeline_info/run_command.sh)"
    echo ""
    echo "Optional paths:"
    echo "  --kraken_db DIR      Kraken2 database directory"
    echo "  --bakta_db PATH      Bakta database [default: db-light from config]"
    echo "  --hmm_databases LIST Comma-separated HMM file paths"
    echo "  --metadata FILE      Sample metadata TSV (flowcell, barcode, lat, lon, etc.)"
    echo ""
    echo "Pipeline flags (passed to Nextflow):"
    echo "  --run_kraken         Kraken2 taxonomic classification (requires --kraken_db)"
    echo "  --annotator STR      Annotator: bakta (default), prokka, or none"
    echo "  --bakta_full         Also run full Bakta annotation (ncRNA/tRNA/CRISPR — slow)"
    echo "  --run_sketch         Sendsketch per-read GTDB classification (local server)"
    echo "  --sendsketch_address URL  Override sendsketch server URL"
    echo "  --run_tetra          Tetranucleotide frequency analysis"
    echo "  --run_db_integration Load results into DuckDB (Python loaders in bin/)"
    echo "  --cleanup            Compress/delete source files after DuckDB import"
    echo "  --watch              Monitor for new FASTQ files (live sequencing)"
    echo "  --db_sync_minutes N  DB sync interval in watch mode [default: 10]"
    echo "  --store_dir DIR      Persistent process output cache (skips completed tasks)"
    echo "  --min_readlen N      Minimum read length in bp [default: 1500]"
    echo "  --keep_percent N     Filtlong keep percent [default: 80]"
    echo "  --min_file_size N    Minimum FASTQ size in bytes [default: 1000000]"
    echo ""
    echo "microscape.app live deploy (watch mode + DB_SYNC):"
    echo "  --deploy_slug SLUG       URL slug for the run on microscape.app"
    echo "  --deploy_name NAME       Display name (quote if it contains spaces)"
    echo "  --deploy_visibility V    Run visibility: private (default) | shared | public"
    echo "                             private — lab-only"
    echo "                             shared  — any signed-in user (any lab)"
    echo "                             public  — also moves the run into the public"
    echo "                                       lab so anonymous web visitors can"
    echo "                                       read it (requires API key with"
    echo "                                       can_publish_public=1)"
    echo "  --deploy_public          DEPRECATED — alias of --deploy_visibility shared"
    echo "  Setting --deploy_slug + --deploy_name writes <outdir>/deploy.sh that"
    echo "  DB_SYNC fires every sync tick. Existing hook is overwritten. Requires"
    echo "  \$MICROSCAPE_API_KEY or ~/.config/microscape/api-key."
    echo ""
    echo "Notes:"
    echo "  - watch vs batch mode is locked per --outdir (recorded in"
    echo "    <outdir>/.pipeline_mode); switching modes on the same outdir is rejected."
    echo ""
    echo "Kitchen sink example (all modules, all options with defaults):"
    echo "  $0 --input /data/run1 --outdir /data/output \\"
    echo "      --run_kraken --kraken_db /path/to/krakendb \\"
    echo "      --annotator bakta \\"
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
RESUME_SESSION=""
AUTO_SESSION=true
DEPLOY_SLUG=""
DEPLOY_NAME=""
DEPLOY_VISIBILITY="private"

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
        --session)
            [[ -z "${2:-}" ]] && die "--session requires a session ID"
            RESUME_SESSION="$2"
            AUTO_SESSION=false
            shift 2 ;;
        --deploy_slug)
            [[ -z "${2:-}" ]] && die "--deploy_slug requires a value"
            DEPLOY_SLUG="$2"
            shift 2 ;;
        --deploy_name)
            [[ -z "${2:-}" ]] && die "--deploy_name requires a value"
            DEPLOY_NAME="$2"
            shift 2 ;;
        --deploy_visibility)
            [[ -z "${2:-}" ]] && die "--deploy_visibility requires a value (private|shared|public)"
            DEPLOY_VISIBILITY="$2"
            shift 2 ;;
        --deploy_public)
            warn "--deploy_public is deprecated; use --deploy_visibility shared"
            DEPLOY_VISIBILITY="shared"
            shift ;;
        *)
            NF_ARGS+=("$1")
            shift ;;
    esac
done

# --deploy_slug and --deploy_name are paired
if [[ -n "$DEPLOY_SLUG" && -z "$DEPLOY_NAME" ]] || [[ -z "$DEPLOY_SLUG" && -n "$DEPLOY_NAME" ]]; then
    die "--deploy_slug and --deploy_name must be provided together"
fi

case "$DEPLOY_VISIBILITY" in
    private|shared|public) ;;
    *) die "--deploy_visibility must be private|shared|public (got '$DEPLOY_VISIBILITY')" ;;
esac

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
# microscape.app deploy hook — DB_SYNC invokes <outdir>/deploy.sh each tick
# ============================================================================

if [[ -n "$DEPLOY_SLUG" ]]; then
    DEPLOY_HOOK="${OUTDIR_HOST}/deploy.sh"
    deploy_tail="    --slug $(printf '%q' "$DEPLOY_SLUG") --name $(printf '%q' "$DEPLOY_NAME") --visibility ${DEPLOY_VISIBILITY}"
    # Quote slug/name through printf %q so spaces and special chars survive
    printf '%s\n' \
        '#!/usr/bin/env bash' \
        "exec ${SCRIPT_DIR}/viz/deploy.sh \\" \
        '    --preprocess-dir "$1/viz" \' \
        "${deploy_tail}" \
        > "$DEPLOY_HOOK"
    chmod +x "$DEPLOY_HOOK"
    echo "[INFO] Wrote deploy hook: $DEPLOY_HOOK (slug=$DEPLOY_SLUG, visibility=$DEPLOY_VISIBILITY)"

    # Soft check for credentials — deploy will still run, just warn early
    if [[ -z "${MICROSCAPE_API_KEY:-}" && ! -f "${HOME}/.config/microscape/api-key" ]]; then
        warn "No \$MICROSCAPE_API_KEY and no ~/.config/microscape/api-key — deploy will fail until one is set."
    fi
fi

# ============================================================================
# Auto-detect session ID from previous runs in this outdir
#
# Three sources, in priority order:
#   1. pipeline_info/session.uuid  — written by the EXIT trap below, survives
#      SIGTERM/SIGINT of run-realtime.sh.
#   2. pipeline_info/run_command.sh — legacy, written only on clean exit.
#   3. <outdir>/.nextflow.log[.N] — last-ditch grep, works as long as the
#      previous run was launched with cwd=outdir (which we now enforce).
# ============================================================================

SESSION_FILE="${OUTDIR_HOST}/pipeline_info/session.uuid"

if [[ "$AUTO_SESSION" == true && -z "$RESUME_SESSION" ]]; then
    if [[ -f "$SESSION_FILE" ]]; then
        RESUME_SESSION=$(tr -d '[:space:]' < "$SESSION_FILE")
    elif [[ -f "${OUTDIR_HOST}/pipeline_info/run_command.sh" ]]; then
        RESUME_SESSION=$(grep -oP '(?<=-resume )[0-9a-f-]{36}' \
            "${OUTDIR_HOST}/pipeline_info/run_command.sh" | tail -1)
    else
        # Fallback: scrape the most recent Session UUID from nextflow's own
        # log in the outdir. Globs include rotated .nextflow.log.[1-9] so a
        # very recent prior session is still recoverable.
        RESUME_SESSION=$(grep -hoP 'Session UUID: \K[0-9a-f-]{36}' \
            "${OUTDIR_HOST}"/.nextflow.log "${OUTDIR_HOST}"/.nextflow.log.[0-9] 2>/dev/null | tail -1)
    fi
    if [[ -n "$RESUME_SESSION" ]]; then
        echo "[INFO] Auto-detected session from previous run: $RESUME_SESSION"
    fi
fi

# Persist the session UUID as soon as nextflow logs it, so a SIGTERM that
# kills the launcher mid-run still leaves a resume marker behind. Reads
# nextflow's own log in the outdir (we cd there below for predictability).
persist_session() {
    local uuid
    uuid=$(grep -hoP 'Session UUID: \K[0-9a-f-]{36}' \
        "${OUTDIR_HOST}"/.nextflow.log "${OUTDIR_HOST}"/.nextflow.log.[0-9] 2>/dev/null | tail -1)
    if [[ -n "$uuid" ]]; then
        mkdir -p "${OUTDIR_HOST}/pipeline_info" 2>/dev/null || true
        echo "$uuid" > "$SESSION_FILE"
    fi
}
trap persist_session EXIT

# nextflow drops .nextflow.log + .nextflow/cache/<uuid>/ in cwd. Pinning
# cwd to the outdir makes both predictable: the resume cache stays with
# the run rather than scattered across whatever shell launched us.
cd "$OUTDIR_HOST"

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
        -resume ${RESUME_SESSION}
    )

    echo "[INFO] Mode:   Docker"
    echo "[INFO] Input:  $INPUT_HOST"
    echo "[INFO] Output: $OUTDIR_HOST"
    [[ -n "$KRAKEN_DB_HOST" ]] && echo "[INFO] Kraken DB: $KRAKEN_DB_HOST"
    echo "[INFO] Running: ${DOCKER_CMD[*]}"
    echo ""

    mkdir -p "${OUTDIR_HOST}/pipeline_info" 2>/dev/null || true

    "${DOCKER_CMD[@]}"
    NF_EXIT=$?

    # Capture session ID from Nextflow log and save command with it for reliable resume
    NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${OUTDIR_HOST}/.nextflow.log" 2>/dev/null | tail -1)
    if [[ -n "$NF_SESSION" ]]; then
        DOCKER_CMD[-1]="-resume ${NF_SESSION}"
    fi
    printf '%s\n' "${DOCKER_CMD[*]}" >> "${OUTDIR_HOST}/pipeline_info/run_command.sh"

    exit $NF_EXIT
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
    -resume ${RESUME_SESSION}
)

echo "[INFO] Mode:   Local (conda)"
echo "[INFO] Input:  $INPUT_HOST"
echo "[INFO] Output: $OUTDIR_HOST"
[[ -n "$KRAKEN_DB_HOST" ]] && echo "[INFO] Kraken DB: $KRAKEN_DB_HOST"
echo "[INFO] Running: ${LOCAL_CMD[*]}"
echo ""

mkdir -p "${OUTDIR_HOST}/pipeline_info" 2>/dev/null || true

"${LOCAL_CMD[@]}"
NF_EXIT=$?

# Capture session ID from Nextflow log and save command with it for reliable resume
NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${OUTDIR_HOST}/.nextflow.log" 2>/dev/null | tail -1)
if [[ -n "$NF_SESSION" ]]; then
    LOCAL_CMD[-1]="-resume ${NF_SESSION}"
fi
printf '%s\n' "${LOCAL_CMD[*]}" >> "${OUTDIR_HOST}/pipeline_info/run_command.sh"

exit $NF_EXIT
