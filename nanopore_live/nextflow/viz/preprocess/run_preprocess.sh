#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<EOF
Usage: $(basename "$0") --input <dir> --output <dir> [options]

Arguments:
  --input DIR         Directory containing nanopore_live run outputs (e.g. /matika/vistara/dana/out_dana_bc/)
  --output DIR        Output directory for JSON files (default: ../public/data)
  --metadata TSV      Optional sample metadata file (TSV with sample_id, barcode, lat, lon, etc.)
  --env PATH          Conda environment path for BBMap (passed to compute_sketches.sh)
  --skip-sketches     Skip comparesketch all-vs-all computation
  --skip-reads        Skip read-level t-SNE computation (saves time for large datasets)
  --max-reads N       Max reads for t-SNE (default: 200000)
  --threads N         Threads for comparesketch (default: 8)
  -h, --help          Show this help

Example:
  bash run_preprocess.sh \\
    --input /matika/vistara/dana/out_dana_bc/ \\
    --metadata /path/to/metadata.tsv \\
    --output /tmp/viz/live_run
EOF
  exit "${1:-0}"
}

INPUT_DIR=""
OUTPUT_DIR=""
METADATA=""
CONDA_ENV=""
SKIP_SKETCHES=false
SKIP_READS=false
MAX_READS=200000
THREADS=8

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input)         INPUT_DIR="$2"; shift 2 ;;
    --output)        OUTPUT_DIR="$2"; shift 2 ;;
    --metadata)      METADATA="$2"; shift 2 ;;
    --env)           CONDA_ENV="$2"; shift 2 ;;
    --skip-sketches) SKIP_SKETCHES=true; shift ;;
    --skip-reads)    SKIP_READS=true; shift ;;
    --max-reads)     MAX_READS="$2"; shift 2 ;;
    --threads)       THREADS="$2"; shift 2 ;;
    -h|--help)       usage 0 ;;
    *)               echo "[ERROR] Unknown argument: $1" >&2; usage 1 ;;
  esac
done

if [[ -z "$INPUT_DIR" ]]; then
  echo "[ERROR] --input is required" >&2
  usage 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
  OUTPUT_DIR="${SCRIPT_DIR}/../public/data"
fi

mkdir -p "$OUTPUT_DIR"

echo "[INFO] Input directory: $INPUT_DIR"
echo "[INFO] Output directory: $OUTPUT_DIR"
echo "[INFO] Metadata file: ${METADATA:-none}"

# Step 1: Comparesketch all-vs-all
if [[ "$SKIP_SKETCHES" != true ]]; then
  echo "[INFO] Computing all-vs-all sketches..."
  SKETCH_ARGS=(
    "$SCRIPT_DIR/compute_sketches.sh"
    --input "$INPUT_DIR"
    --output "$OUTPUT_DIR/sketch_distances.tsv"
    --threads "$THREADS"
  )
  if [[ -n "$CONDA_ENV" ]]; then
    SKETCH_ARGS+=(--env "$CONDA_ENV")
  fi
  bash "${SKETCH_ARGS[@]}"
else
  echo "[INFO] Skipping sketch computation"
fi

# Step 2: Main preprocessing (DuckDB -> JSON)
echo "[INFO] Running main preprocessing..."
PREPROCESS_ARGS=(
  "$SCRIPT_DIR/preprocess.py"
  --input "$INPUT_DIR"
  --output "$OUTPUT_DIR"
)
if [[ -n "$METADATA" ]]; then
  PREPROCESS_ARGS+=(--metadata "$METADATA")
fi
if [[ -f "$OUTPUT_DIR/sketch_distances.tsv" ]]; then
  PREPROCESS_ARGS+=(--sketch-distances "$OUTPUT_DIR/sketch_distances.tsv")
fi
python3 "${PREPROCESS_ARGS[@]}"

# Step 3: Read-level t-SNE
if [[ "$SKIP_READS" != true ]]; then
  echo "[INFO] Computing read-level t-SNE..."
  python3 "$SCRIPT_DIR/compute_read_tsne.py" \
    --input "$INPUT_DIR" \
    --output "$OUTPUT_DIR" \
    --max-reads "$MAX_READS"
else
  echo "[INFO] Skipping read t-SNE"
fi

# Compress large JSON files
echo "[INFO] Compressing JSON files..."
for f in "$OUTPUT_DIR"/*.json; do
  if [[ -f "$f" ]] && [[ $(stat -c%s "$f" 2>/dev/null || stat -f%z "$f") -gt 100000 ]]; then
    gzip -kf "$f"
  fi
done

echo "[SUCCESS] Preprocessing complete. Output in: $OUTPUT_DIR"
echo "[INFO] Start the dev server with: VIZ_DATA_DIR=$OUTPUT_DIR npm run dev"
