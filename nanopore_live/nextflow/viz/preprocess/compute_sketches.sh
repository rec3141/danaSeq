#!/usr/bin/env bash
set -euo pipefail
# Compute all-vs-all minhash distances between samples using BBMap sketch tools.
# Phase 1: sketch.sh per sample (streaming, no temp FASTA copies)
# Phase 2: comparesketch.sh alltoall on small sketch files
#
# Usage: compute_sketches.sh --input <dir> --output <tsv> [--env <conda_env_path>] [--threads N]

usage() {
  cat <<EOF
Usage: $(basename "$0") --input <dir> --output <tsv> [--env <path>] [--threads N]

  --input DIR     Directory containing nanopore_live outputs (flowcell/barcode/fa/)
  --output TSV    Output pairwise distance file
  --env PATH      Path to conda environment containing BBMap (e.g. /path/to/envs/dana-bbmap)
  --threads N     Threads for comparesketch (default: 8)
  -h, --help      Show this help
EOF
  exit "${1:-0}"
}

INPUT_DIR=""
OUTPUT_TSV=""
CONDA_ENV=""
THREADS=8

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)   INPUT_DIR="$2"; shift 2 ;;
    --output)  OUTPUT_TSV="$2"; shift 2 ;;
    --env)     CONDA_ENV="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage 0 ;;
    *)         echo "[ERROR] Unknown argument: $1" >&2; usage 1 ;;
  esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_TSV" ]]; then
  echo "[ERROR] --input and --output are required" >&2
  usage 1
fi

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

# Persistent sketch dir next to output so work survives restarts
SKETCH_DIR="$(dirname "$OUTPUT_TSV")/sketches"
mkdir -p "$SKETCH_DIR"

# Locate BBMap tools
BBMAP_DIR=""

# 1. Explicit --env
if [[ -n "$CONDA_ENV" ]]; then
  if [[ -x "$CONDA_ENV/bin/comparesketch.sh" ]]; then
    BBMAP_DIR="$CONDA_ENV/bin"
    export PATH="$CONDA_ENV/bin:$PATH"
  else
    echo "[ERROR] comparesketch.sh not found in --env $CONDA_ENV/bin/" >&2
    exit 1
  fi
fi

# 2. Already on PATH
if [[ -z "$BBMAP_DIR" ]]; then
  local_cs="$(which comparesketch.sh 2>/dev/null || echo "")"
  if [[ -n "$local_cs" ]]; then
    BBMAP_DIR="$(dirname "$local_cs")"
  fi
fi

# 3. Try conda env by name (dana-bbmap)
if [[ -z "$BBMAP_DIR" ]]; then
  for base in "${CONDA_PREFIX:-}" "$(conda info --base 2>/dev/null || echo "")" "$HOME/miniconda3" "$HOME/miniforge3" "$HOME/anaconda3"; do
    [[ -z "$base" ]] && continue
    env_path="$base/envs/dana-bbmap"
    if [[ -x "$env_path/bin/comparesketch.sh" ]]; then
      BBMAP_DIR="$env_path/bin"
      export PATH="$env_path/bin:$PATH"
      break
    fi
  done
fi

if [[ -z "$BBMAP_DIR" ]]; then
  echo "[ERROR] BBMap tools (sketch.sh, comparesketch.sh) not found." >&2
  echo "[ERROR] Provide --env <conda_env_path>, activate dana-bbmap, or install BBMap." >&2
  exit 1
fi

SKETCH="$BBMAP_DIR/sketch.sh"
COMPARESKETCH="$BBMAP_DIR/comparesketch.sh"

echo "[INFO] Using BBMap tools from: $BBMAP_DIR"

# Phase 1: Sketch each sample (streaming, parallel)
echo "[INFO] Phase 1: Sketching individual samples ($THREADS parallel)..."

# Build job list, skipping samples that already have sketch files
JOBFILE="$TMPDIR/sketch_jobs.tsv"
> "$JOBFILE"
SKIPPED=0
for sample_dir in "$INPUT_DIR"/*/; do
  sample_name="$(basename "$sample_dir")"
  for barcode_dir in "$sample_dir"/barcode*/; do
    [[ -d "$barcode_dir" ]] || continue
    barcode="$(basename "$barcode_dir")"
    fa_dir="$barcode_dir/fa"
    [[ -d "$fa_dir" ]] || continue

    safe_name="${sample_name}_${barcode}"

    # Skip if sketch already exists
    if [[ -s "$SKETCH_DIR/${safe_name}.sketch" ]]; then
      SKIPPED=$((SKIPPED + 1))
      continue
    fi

    # Check for non-empty FASTA files
    if find "$fa_dir" -name "*.fa" -size +0c -maxdepth 1 -print -quit 2>/dev/null | grep -q .; then
      printf '%s\t%s\n' "$safe_name" "$fa_dir"
    fi
  done
done >> "$JOBFILE"

TOTAL_JOBS="$(wc -l < "$JOBFILE")"
echo "[INFO] Found $TOTAL_JOBS samples to sketch ($SKIPPED already done)"

# Run sketches in parallel using background jobs with a concurrency limiter
DONE_COUNT=0
RUNNING=0

while IFS=$'\t' read -r safe_name fa_dir; do
  (
    find "$fa_dir" -name "*.fa" -size +0c -print0 2>/dev/null \
      | xargs -0 cat \
      | "$SKETCH" in=stdin.fa out="$SKETCH_DIR/${safe_name}.sketch" \
          name="$safe_name" threads=1 -Xmx1g 2>/dev/null
  ) &
  RUNNING=$((RUNNING + 1))

  if (( RUNNING >= THREADS )); then
    wait -n 2>/dev/null || true
    RUNNING=$((RUNNING - 1))
    DONE_COUNT=$((DONE_COUNT + 1))
    if (( DONE_COUNT % 50 == 0 )); then
      echo "  [INFO] Sketched ~$DONE_COUNT / $TOTAL_JOBS samples..."
    fi
  fi
done < "$JOBFILE"

# Wait for remaining jobs
wait 2>/dev/null || true

SAMPLE_COUNT="$(find "$SKETCH_DIR" -name "*.sketch" | wc -l)"
echo "[INFO] Sketched $SAMPLE_COUNT total samples"

if [[ $SAMPLE_COUNT -lt 2 ]]; then
  echo "[WARNING] Need at least 2 samples for all-vs-all comparison. Skipping."
  exit 0
fi

# Check sketch dir size (should be tiny)
sketch_size="$(du -sh "$SKETCH_DIR" | cut -f1)"
echo "[INFO] Sketch files total: $sketch_size"

# Phase 2: Compare all sketches
echo "[INFO] Phase 2: Running comparesketch alltoall ($SAMPLE_COUNT samples, $THREADS threads)..."
cd "$SKETCH_DIR"
"$COMPARESKETCH" alltoall *.sketch \
  format=3 ow \
  out="$OUTPUT_TSV" \
  threads="$THREADS" \
  2>&1 | grep -v "^$" | head -20 || true

echo "[INFO] Sketch distances written to: $OUTPUT_TSV"
if [[ -f "$OUTPUT_TSV" ]]; then
  lines="$(wc -l < "$OUTPUT_TSV")"
  echo "[INFO] $lines pairwise comparisons"
fi
