#!/usr/bin/env bash
set -euo pipefail

# Generic microscape.app deployer for any dānaSeq viz SPA.
#
# Builds the Svelte SPA at --viz-dir, stages preprocess JSONs from
# --preprocess-dir into <viz-dir>/dist/data/, tars dist/, and POSTs to
# the microscape.app ingest endpoint.
#
# Per-pipeline wrappers (e.g. nanopore_live/viz/deploy.sh) pin --viz-dir
# and a default --pipeline; per-run hooks in each outdir then invoke
# the wrapper with --slug and --name. Lifted out of nanopore_live in
# 2026-06 so the same machinery works for mag_analysis, illumina_rna, etc.

ENDPOINT_DEFAULT="https://microscape.app/api/v1/deploy"
API_KEY_FILE_DEFAULT="$HOME/.config/microscape/api-key"

usage() {
    cat <<EOF
Usage: $(basename "$0") --viz-dir DIR --preprocess-dir DIR --slug NAME --pipeline SLUG [options]

Required:
  --viz-dir DIR          Vite SPA root (the dir with package.json + src/);
                         dist/ is built and staged here.
  --preprocess-dir DIR   Directory containing preprocess JSON / JSON.gz files.
  --slug NAME            URL slug on microscape.app (also tarball name).
                         Lowercase alnum + '-'/'_', 1-64 chars.
  --pipeline SLUG        Pipeline slug registered on the server
                         (e.g. danaseq-nanopore-live, danaseq-mag-analysis).

Optional:
  --name "Display"       Human-readable run name           [default: --slug value]
  --endpoint URL         Ingest endpoint                    [default: ${ENDPOINT_DEFAULT}]
  --api-key-file PATH    File with MICROSCAPE_API_KEY       [default: ${API_KEY_FILE_DEFAULT}]
  --api-key KEY          Inline API key (overrides file/env)
  --skip-build           Reuse existing <viz-dir>/dist/
  --visibility V         private (default) | shared | public
                           private — lab-only
                           shared  — any signed-in user
                           public  — also moves the run into the public lab
                                     (requires can_publish_public=1 on the key)
  --dry-run              Build + tarball but do not POST
  -h, --help             Show this help

Environment:
  MICROSCAPE_API_KEY     Overrides --api-key-file when set
EOF
}

VIZ_DIR=""
PREPROCESS_DIR=""
SLUG=""
PIPELINE=""
DISPLAY_NAME=""
ENDPOINT="$ENDPOINT_DEFAULT"
API_KEY_FILE="$API_KEY_FILE_DEFAULT"
API_KEY_INLINE=""
SKIP_BUILD=false
DRY_RUN=false
VISIBILITY="private"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --viz-dir)        VIZ_DIR="$2"; shift 2 ;;
        --preprocess-dir) PREPROCESS_DIR="$2"; shift 2 ;;
        --slug)           SLUG="$2"; shift 2 ;;
        --pipeline)       PIPELINE="$2"; shift 2 ;;
        --name)           DISPLAY_NAME="$2"; shift 2 ;;
        --endpoint)       ENDPOINT="$2"; shift 2 ;;
        --api-key-file)   API_KEY_FILE="$2"; shift 2 ;;
        --api-key)        API_KEY_INLINE="$2"; shift 2 ;;
        --skip-build)     SKIP_BUILD=true; shift ;;
        --visibility)     VISIBILITY="$2"; shift 2 ;;
        --dry-run)        DRY_RUN=true; shift ;;
        -h|--help)        usage; exit 0 ;;
        *)  echo "[ERROR] Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

case "$VISIBILITY" in
    private|shared|public) ;;
    *) echo "[ERROR] --visibility must be private|shared|public (got '$VISIBILITY')" >&2; exit 2 ;;
esac

die() { echo "[ERROR] $1" >&2; exit 1; }
info() { echo "[deploy] $*" >&2; }

[[ -n "$VIZ_DIR" ]]         || die "--viz-dir is required"
[[ -n "$PREPROCESS_DIR" ]]  || die "--preprocess-dir is required"
[[ -n "$SLUG" ]]            || die "--slug is required"
[[ -n "$PIPELINE" ]]        || die "--pipeline is required"
[[ -d "$VIZ_DIR" ]]         || die "viz dir not found: $VIZ_DIR"
[[ -f "$VIZ_DIR/package.json" ]] || die "no package.json under $VIZ_DIR (not a Vite SPA root)"
[[ -d "$PREPROCESS_DIR" ]]  || die "preprocess dir not found: $PREPROCESS_DIR"
[[ "$SLUG" =~ ^[a-z0-9][a-z0-9_-]{0,63}$ ]] \
    || die "slug must be lowercase alnum + '-'/'_', 1-64 chars: got '$SLUG'"
[[ -z "$DISPLAY_NAME" ]] && DISPLAY_NAME="$SLUG"

if [[ -n "$API_KEY_INLINE" ]]; then
    API_KEY="$API_KEY_INLINE"
elif [[ -n "${MICROSCAPE_API_KEY:-}" ]]; then
    API_KEY="$MICROSCAPE_API_KEY"
else
    [[ -f "$API_KEY_FILE" ]] || die "api-key file not found: $API_KEY_FILE"
    API_KEY="$(tr -d '\n' < "$API_KEY_FILE")"
fi
[[ -n "$API_KEY" ]] || die "API key is empty"

shopt -s nullglob
preproc_files=("$PREPROCESS_DIR"/*.json "$PREPROCESS_DIR"/*.json.gz)
shopt -u nullglob
[[ ${#preproc_files[@]} -gt 0 ]] || die "no *.json or *.json.gz under $PREPROCESS_DIR"

cd "$VIZ_DIR"

if [[ "$SKIP_BUILD" != true ]]; then
    if [[ ! -d node_modules ]]; then
        info "node_modules/ missing — running npm install"
        npm install --no-audit --no-fund
    fi
    info "building SPA at $VIZ_DIR (vite build)"
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
    trap - EXIT
    exit 0
fi

info "POST $ENDPOINT  slug=$SLUG  pipeline=$PIPELINE  visibility=$VISIBILITY"
response=$(curl --fail -sS -X POST \
    -H "Authorization: Bearer $API_KEY" \
    -H "Content-Type: application/gzip" \
    -H "X-Microscape-Slug: $SLUG" \
    -H "X-Microscape-Pipeline: $PIPELINE" \
    -H "X-Microscape-Name: $DISPLAY_NAME" \
    -H "X-Microscape-Visibility: $VISIBILITY" \
    --data-binary "@$TARBALL" \
    "$ENDPOINT")
echo "$response"
