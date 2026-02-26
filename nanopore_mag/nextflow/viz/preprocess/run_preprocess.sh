#!/usr/bin/env bash
# Run preprocessing via the dana-mag-semibin conda env (has pandas, scipy, sklearn)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VIZ_DIR="$(dirname "$SCRIPT_DIR")"
NEXTFLOW_DIR="$(dirname "$VIZ_DIR")"

# Resolve results path: storeDir (persistent cache) preferred over publishDir
# Priority: RESULTS_DIR env var > --store_dir from run_command.txt > --outdir > glob fallback
if [[ -z "${RESULTS_DIR:-}" ]]; then
    # Find most recent run_command.txt under any results_*/pipeline_info/
    RUN_CMD=$(ls -t "${NEXTFLOW_DIR}"/results_*/pipeline_info/run_command.txt 2>/dev/null | head -1)
    if [[ -n "$RUN_CMD" ]]; then
        STORE_DIR=$(grep -oP '(?<=--store_dir )\S+' "$RUN_CMD" || true)
        OUTDIR=$(grep -oP '(?<=--outdir )\S+' "$RUN_CMD" || true)
    fi
    if [[ -n "${STORE_DIR:-}" && -d "$STORE_DIR" ]]; then
        RESULTS_DIR="$STORE_DIR"
    elif [[ -n "${OUTDIR:-}" && -d "$OUTDIR" ]]; then
        RESULTS_DIR="$OUTDIR"
    else
        RESULTS_DIR=$(ls -td "${NEXTFLOW_DIR}"/results_* 2>/dev/null | head -1)
    fi
fi
OUTPUT_DIR="${VIZ_DIR}/public/data"

CONDA_ENV="${NEXTFLOW_DIR}/conda-envs/dana-mag-semibin"

if [[ ! -d "$CONDA_ENV" ]]; then
    echo "[ERROR] Conda env not found: $CONDA_ENV"
    echo "Run: cd ${NEXTFLOW_DIR} && ./install.sh"
    exit 1
fi

if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "[ERROR] Results directory not found: $RESULTS_DIR"
    echo "Set RESULTS_DIR=/path/to/results or place results in ${NEXTFLOW_DIR}/results_*"
    exit 1
fi

echo "[INFO] Preprocessing pipeline results: $RESULTS_DIR"
echo "[INFO] Output directory: $OUTPUT_DIR"

EXTRA_ARGS=""
if [[ "${SKIP_TSNE:-}" == "1" ]]; then
    EXTRA_ARGS="--skip-tsne"
fi
if [[ -n "${STORE_DIR:-}" && -d "$STORE_DIR" ]]; then
    EXTRA_ARGS="$EXTRA_ARGS --store-dir $STORE_DIR"
fi

mamba run -p "$CONDA_ENV" python3 "$SCRIPT_DIR/preprocess.py" \
    --results "$RESULTS_DIR" \
    --output "$OUTPUT_DIR" \
    $EXTRA_ARGS

echo "[SUCCESS] Preprocessing complete. JSON files in $OUTPUT_DIR"
