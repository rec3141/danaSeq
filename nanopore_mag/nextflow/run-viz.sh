#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana MAG Viz Dashboard — build + serve
# ============================================================================
#
# Rebuilds the Svelte dashboard from preprocessed JSON and starts a preview
# server. Optionally re-runs the Python preprocessor first (--preprocess).
#
# Usage:
#   ./run-viz.sh --outdir /path/to/results [options]
#
# Examples:
#   # Rebuild frontend only (data already preprocessed)
#   ./run-viz.sh --outdir /data/minknow/QEI2025/nanopore_mag/tmpdir
#
#   # Re-preprocess data then rebuild
#   ./run-viz.sh --outdir /data/minknow/QEI2025/nanopore_mag/tmpdir --preprocess
#
#   # Re-preprocess, skip t-SNE/UMAP (reuse existing embeddings)
#   ./run-viz.sh --outdir ... --preprocess --skip-tsne --skip-umap
#
#   # Custom port
#   ./run-viz.sh --outdir /data/minknow/QEI2025/nanopore_mag/tmpdir --port 5175
#
#   # Build only (no server)
#   ./run-viz.sh --outdir /data/minknow/QEI2025/nanopore_mag/tmpdir --no-serve
#

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VIZ_DIR="${SCRIPT_DIR}/viz"
CONDA_ENV="${SCRIPT_DIR}/conda-envs/dana-mag-viz"

# Defaults
OUTDIR=""
STORE_DIR=""
PORT=5174
PREPROCESS=false
SERVE=true
SKIP_TSNE=false
SKIP_UMAP=false

usage() {
    sed -n '3,/^$/{ s/^# \{0,1\}//; p }' "$0"
    exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)       OUTDIR="$2"; shift 2 ;;
        --store-dir|--store_dir) STORE_DIR="$2"; shift 2 ;;
        --port)         PORT="$2"; shift 2 ;;
        --preprocess)   PREPROCESS=true; shift ;;
        --skip-tsne|--skip_tsne) SKIP_TSNE=true; shift ;;
        --skip-umap|--skip_umap) SKIP_UMAP=true; shift ;;
        --no-serve)     SERVE=false; shift ;;
        -h|--help)      usage 0 ;;
        *)              echo "Unknown option: $1" >&2; usage 1 ;;
    esac
done

if [[ -z "$OUTDIR" ]]; then
    echo "Error: --outdir is required" >&2
    usage 1
fi

# Resolve conda env binaries
if [[ -d "$CONDA_ENV/bin" ]]; then
    NPM="${CONDA_ENV}/bin/npm"
    PYTHON="${CONDA_ENV}/bin/python3"
else
    echo "Warning: conda env not found at ${CONDA_ENV}, using system binaries" >&2
    NPM="npm"
    PYTHON="python3"
fi

VIZ_DATA="${OUTDIR}/viz/data"
mkdir -p "$VIZ_DATA"

# --- Optional: re-run preprocessor -------------------------------------------
if $PREPROCESS; then
    echo "==> Preprocessing pipeline results → JSON..."
    store_flag=""
    [[ -n "$STORE_DIR" ]] && store_flag="--store-dir ${STORE_DIR}"
    tsne_flag=""
    $SKIP_TSNE && tsne_flag="--skip-tsne"
    umap_flag=""
    $SKIP_UMAP && umap_flag="--skip-umap"
    $PYTHON "${VIZ_DIR}/preprocess/preprocess.py" \
        --results "$OUTDIR" \
        --output "$VIZ_DATA/" \
        $store_flag $tsne_flag $umap_flag

    # Regenerate genes.json
    find_first() { for f in "$@"; do [ -f "$f" ] && echo "$f" && return; done; }
    ANNOT_TSV=$(find_first \
        ${STORE_DIR:+"${STORE_DIR}/annotation/bakta/extra/annotation.tsv"} \
        "${OUTDIR}/annotation/bakta/extra/annotation.tsv" \
        ${STORE_DIR:+"${STORE_DIR}/annotation/bakta/basic/annotation.tsv"} \
        "${OUTDIR}/annotation/bakta/basic/annotation.tsv" \
        ${STORE_DIR:+"${STORE_DIR}/annotation/prokka/annotation.tsv"} \
        "${OUTDIR}/annotation/prokka/annotation.tsv")
    RRNA_TSV=$(find_first \
        ${STORE_DIR:+"${STORE_DIR}/taxonomy/rrna/rrna_genes.tsv"} \
        "${OUTDIR}/taxonomy/rrna/rrna_genes.tsv") || true
    TRNA_TSV=$(find_first \
        ${STORE_DIR:+"${STORE_DIR}/taxonomy/rrna/trna_genes.tsv"} \
        "${OUTDIR}/taxonomy/rrna/trna_genes.tsv") || true
    GENE_DEPTHS=$(find_first \
        ${STORE_DIR:+"${STORE_DIR}/mapping/gene_depths.tsv"} \
        "${OUTDIR}/mapping/gene_depths.tsv") || true
    ASSEMBLY=$(find_first \
        ${STORE_DIR:+"${STORE_DIR}/assembly/assembly.fasta"} \
        "${OUTDIR}/assembly/assembly.fasta") || true

    if [[ -n "${ANNOT_TSV:-}" ]]; then
        echo "==> Building genes.json..."
        $PYTHON "${VIZ_DIR}/preprocess/genes_to_json.py" \
            "$ANNOT_TSV" "${VIZ_DATA}/genes.json" \
            "${RRNA_TSV:-}" "${TRNA_TSV:-}" "${GENE_DEPTHS:-}" "${ASSEMBLY:-}"
    else
        echo "Warning: no annotation TSV found, writing empty genes.json" >&2
        echo '{}' > "${VIZ_DATA}/genes.json"
        echo '{}' | gzip > "${VIZ_DATA}/genes.json.gz"
    fi
fi

# --- Copy data into public/ for build ----------------------------------------
echo "==> Copying data files..."
cp "${VIZ_DATA}"/* "${VIZ_DIR}/public/data/" 2>/dev/null || true

# --- Install deps if needed --------------------------------------------------
cd "$VIZ_DIR"
if [[ ! -d node_modules ]]; then
    echo "==> Installing npm dependencies..."
    $NPM ci --prefer-offline 2>/dev/null || $NPM install --no-audit --no-fund
fi

# --- Build --------------------------------------------------------------------
echo "==> Building Svelte app..."
$NPM run build

# Copy build output to results
cp -r dist "${OUTDIR}/viz/site" 2>/dev/null || true

# --- Serve --------------------------------------------------------------------
if $SERVE; then
    # Kill any existing preview server on the same port
    pkill -f "node.*vite preview.*--port ${PORT}" 2>/dev/null || true
    sleep 1

    echo "==> Starting preview server on port ${PORT}..."
    setsid nohup $NPM run serve -- --host 0.0.0.0 --port "$PORT" > /tmp/vite_preview.log 2>&1 &

    sleep 2
    SERVER_IP=$(hostname -I | awk '{print $1}')
    echo ""
    echo "  Viz dashboard: http://${SERVER_IP}:${PORT}/"
    echo ""
fi
