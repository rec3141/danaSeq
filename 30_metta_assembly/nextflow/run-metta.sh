#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# METTA Assembly Pipeline Launcher
# ============================================================================
#
# Runs the Nextflow pipeline locally via conda.
#
# Usage:
#   ./run-metta.sh --input /path/to/reads --outdir /path/to/output [options]
#
# Examples:
#   # Basic per-sample assembly
#   ./run-metta.sh --input /data/reads --outdir /data/output
#
#   # Co-assembly mode
#   ./run-metta.sh --input /data/reads --outdir /data/output --coassembly
#
#   # Skip some assemblers
#   ./run-metta.sh --input /data/reads --outdir /data/output \
#       --run_megahit false --run_spades false
#
#   # SLURM profile
#   ./run-metta.sh --input /data/reads --outdir /data/output -profile slurm
#
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NF_ARGS=()
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
    echo "Pipeline flags (passed to Nextflow):"
    echo "  --coassembly         Co-assemble all samples (default: per-sample)"
    echo "  --run_normalize BOOL Enable bbnorm normalization [default: true]"
    echo "  --run_tadpole BOOL   Run Tadpole assembler [default: true]"
    echo "  --run_megahit BOOL   Run Megahit assembler [default: true]"
    echo "  --run_spades BOOL    Run SPAdes assembler [default: true]"
    echo "  --run_metaspades BOOL Run metaSPAdes assembler [default: true]"
    echo "  --min_readlen N      Minimum read length [default: 70]"
    echo "  --dedupe_identity N  Deduplication identity threshold [default: 98]"
    echo "  --metabat_min_cls N  MetaBAT2 minimum cluster size [default: 2000]"
    echo "  --assembly_cpus N    CPUs for assembly [default: 24]"
    echo "  --assembly_memory S  Memory for assembly [default: '250 GB']"
    echo "  --store_dir DIR      Persistent cache directory (storeDir)"
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
RESUME_SESSION=""
AUTO_SESSION=true

while (( $# )); do
    case "$1" in
        -h|--help)
            usage ;;
        --help-pipeline)
            mamba run -p "${SCRIPT_DIR}/conda-envs/dana-metta-bbmap" \
                nextflow run "${SCRIPT_DIR}/main.nf" --help
            exit 0 ;;
        --input)
            [[ -z "${2:-}" ]] && die "--input requires a directory path"
            INPUT_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            shift 2 ;;
        --outdir)
            [[ -z "${2:-}" ]] && die "--outdir requires a directory path"
            OUTDIR_HOST="$(realpath -m "$2")"
            shift 2 ;;
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
# Local mode — run via conda
# ============================================================================

# Nextflow's conda activation calls `conda info --json` to find the base prefix.
# Ensure conda/mamba base bin/ is on PATH so this works inside `mamba run`.
CONDA_BASE_BIN="$(dirname "$(which mamba 2>/dev/null || which conda)")"
export PATH="${CONDA_BASE_BIN}:${PATH}"

LOCAL_CMD=(
    mamba run -p "${SCRIPT_DIR}/conda-envs/dana-metta-bbmap"
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
NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${SCRIPT_DIR}/.nextflow.log" 2>/dev/null | tail -1)
save_run_command "$OUTDIR_HOST" "$NF_SESSION"

exit $NF_EXIT
