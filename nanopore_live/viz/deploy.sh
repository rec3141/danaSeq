#!/usr/bin/env bash
set -euo pipefail

# Deploy the nanopore_live viz SPA + a preprocessed-data directory to the
# microscape.app portal via the ingest API. Designed to be invoked:
#   (a) manually — after a one-off `run_preprocess.sh` run
#   (b) automatically from DB_SYNC's deploy-hook each sync tick
#
# The companion per-run hook at <outdir>/deploy.sh typically `exec`s this
# script with the appropriate flags; see README / examples/ for a template.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENDPOINT_DEFAULT="https://microscape.app/api/v1/deploy"
PIPELINE_DEFAULT="danaseq-nanopore-live"
API_KEY_FILE_DEFAULT="$HOME/.config/microscape/api-key"

usage() {
    cat <<EOF
Usage: $(basename "$0") --preprocess-dir DIR --slug NAME [options]

Required:
  --preprocess-dir DIR   Directory containing preprocess JSONs
                         (e.g. <outdir>/viz or a standalone preprocess output)
  --slug NAME            URL slug on microscape.app (also tarball name).
                         Lowercase alnum + dashes, <= 64 chars.

Optional:
  --name "Display"       Human-readable display name   [default: --slug value]
  --pipeline SLUG        Pipeline slug registered on the server
                         [default: ${PIPELINE_DEFAULT}]
  --endpoint URL         Ingest endpoint
                         [default: ${ENDPOINT_DEFAULT}]
  --api-key-file PATH    File containing MICROSCAPE_API_KEY
                         [default: ${API_KEY_FILE_DEFAULT}]
  --api-key KEY          Inline API key (overrides --api-key-file)
  --skip-build           Reuse existing viz/dist/ (useful for rapid iteration)
  --dry-run              Build + tarball but do not POST
  -h, --help             Show this help

Environment:
  MICROSCAPE_API_KEY     Overrides --api-key-file when set

Example (manual):
  $(basename "$0") --preprocess-dir /data/scratch/nanopore_live/genice_ci \\
                   --slug genice_ci --name "GenIce CI"

Example (DB_SYNC hook, <outdir>/deploy.sh):
  #!/usr/bin/env bash
  exec ${SCRIPT_DIR}/deploy.sh \\
      --preprocess-dir "\$1/viz" \\
      --slug genice_ci --name "GenIce CI"
EOF
}

PREPROCESS_DIR=""
SLUG=""
DISPLAY_NAME=""
PIPELINE="$PIPELINE_DEFAULT"
ENDPOINT="$ENDPOINT_DEFAULT"
API_KEY_FILE="$API_KEY_FILE_DEFAULT"
API_KEY_INLINE=""
SKIP_BUILD=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --preprocess-dir) PREPROCESS_DIR="$2"; shift 2 ;;
        --slug)           SLUG="$2"; shift 2 ;;
        --name)           DISPLAY_NAME="$2"; shift 2 ;;
        --pipeline)       PIPELINE="$2"; shift 2 ;;
        --endpoint)       ENDPOINT="$2"; shift 2 ;;
        --api-key-file)   API_KEY_FILE="$2"; shift 2 ;;
        --api-key)        API_KEY_INLINE="$2"; shift 2 ;;
        --skip-build)     SKIP_BUILD=true; shift ;;
        --dry-run)        DRY_RUN=true; shift ;;
        -h|--help)        usage; exit 0 ;;
        *)  echo "[ERROR] Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

die() { echo "[ERROR] $1" >&2; exit 1; }
info() { echo "[deploy] $*" >&2; }

[[ -n "$PREPROCESS_DIR" ]] || die "--preprocess-dir is required"
[[ -n "$SLUG" ]]            || die "--slug is required"
[[ -d "$PREPROCESS_DIR" ]]   || die "preprocess dir not found: $PREPROCESS_DIR"
[[ "$SLUG" =~ ^[a-z0-9][a-z0-9_-]{0,63}$ ]] \
    || die "slug must be lowercase alnum + '-'/'_', 1-64 chars: got '$SLUG'"
[[ -z "$DISPLAY_NAME" ]] && DISPLAY_NAME="$SLUG"

# Resolve API key — inline > env > file
if [[ -n "$API_KEY_INLINE" ]]; then
    API_KEY="$API_KEY_INLINE"
elif [[ -n "${MICROSCAPE_API_KEY:-}" ]]; then
    API_KEY="$MICROSCAPE_API_KEY"
else
    [[ -f "$API_KEY_FILE" ]] || die "api-key file not found: $API_KEY_FILE"
    API_KEY="$(tr -d '\n' < "$API_KEY_FILE")"
fi
[[ -n "$API_KEY" ]] || die "API key is empty"

# Count staged JSON files up-front so we fail fast on an empty preprocess dir
shopt -s nullglob
preproc_files=("$PREPROCESS_DIR"/*.json "$PREPROCESS_DIR"/*.json.gz)
shopt -u nullglob
[[ ${#preproc_files[@]} -gt 0 ]] || die "no *.json or *.json.gz under $PREPROCESS_DIR"

cd "$SCRIPT_DIR"

if [[ "$SKIP_BUILD" != true ]]; then
    if [[ ! -d node_modules ]]; then
        info "node_modules/ missing — running npm install"
        npm install --no-audit --no-fund
    fi
    info "building SPA (vite build)"
    npm run build >/dev/null
else
    [[ -f dist/index.html ]] || die "--skip-build set but no dist/index.html"
    info "skipping build — reusing existing dist/"
fi

info "staging preprocess output (${#preproc_files[@]} files) into dist/data/"
rm -rf dist/data
mkdir -p dist/data
cp "${preproc_files[@]}" dist/data/

TARBALL="$(mktemp -t "${SLUG}.XXXXXX.tgz")"
trap 'rm -f "$TARBALL"' EXIT
tar -czf "$TARBALL" -C dist .
tarball_bytes=$(stat -c%s "$TARBALL" 2>/dev/null || stat -f%z "$TARBALL")
info "tarball: $TARBALL ($((tarball_bytes / 1024 / 1024)) MB)"

if [[ "$DRY_RUN" == true ]]; then
    info "--dry-run: skipping POST"
    echo "tarball=$TARBALL"
    trap - EXIT   # keep tarball for inspection
    exit 0
fi

info "POST $ENDPOINT  slug=$SLUG  pipeline=$PIPELINE"
response=$(curl --fail -sS -X POST \
    -H "Authorization: Bearer $API_KEY" \
    -H "Content-Type: application/gzip" \
    -H "X-Microscape-Slug: $SLUG" \
    -H "X-Microscape-Pipeline: $PIPELINE" \
    -H "X-Microscape-Name: $DISPLAY_NAME" \
    --data-binary "@$TARBALL" \
    "$ENDPOINT")
echo "$response"
