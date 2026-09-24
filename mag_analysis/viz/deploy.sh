#!/usr/bin/env bash
# Publish a mag_analysis run's viz to microscape.app. Wraps
# tools/microscape_deploy.sh with this pipeline's SPA and slug.
#
#   deploy.sh --preprocess-dir <viz>/data --slug SLUG [--name N]
#             [--visibility private|shared|public] [--api-key-file F] [--dry-run]
#
# Runs where the pipeline runs, including inside the container, whose copy of
# this SPA is read-only. The SPA is therefore built once per pipeline build into
# <viz>/.deploy_build and reused by every later snapshot of the run. Where a
# file exists as both .json and .json.gz only the .gz is sent: the SPA fetches
# <file>.json.gz first and falls back to .json.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The image carries a copy of the shared deployer beside this script (CI stages
# it, since the image is built from mag_analysis/ alone); a checkout uses the
# repo's own.
DEPLOYER="$SCRIPT_DIR/microscape_deploy.sh"
[[ -f "$DEPLOYER" ]] || DEPLOYER="$SCRIPT_DIR/../../tools/microscape_deploy.sh"
[[ -f "$DEPLOYER" ]] || { echo "[ERROR] microscape_deploy.sh not found" >&2; exit 1; }

PRE=""
REST=()
while (( $# )); do
    case "$1" in
        --preprocess-dir) PRE="$2"; shift 2 ;;
        *) REST+=("$1"); shift ;;
    esac
done
[[ -d "$PRE" ]] || { echo "[ERROR] --preprocess-dir not found: $PRE" >&2; exit 1; }
PRE="$(cd "$PRE" && pwd)"

BUILD="$(dirname "$PRE")/.deploy_build"
STAMP="${DANASEQ_GIT_SHA:-$(git -C "$SCRIPT_DIR" rev-parse HEAD 2>/dev/null || echo dev)}"
if [[ ! -f "$BUILD/dist/index.html" || "$(cat "$BUILD/.stamp" 2>/dev/null)" != "$STAMP" ]]; then
    echo "[deploy] building SPA for $STAMP into $BUILD" >&2
    rm -rf "$BUILD"
    mkdir -p "$BUILD"
    cp -r "$SCRIPT_DIR"/. "$BUILD"/
    rm -rf "$BUILD/dist"
    (
        cd "$BUILD"
        [[ -d node_modules ]] || npm ci --no-audit --no-fund >/dev/null 2>&1 \
            || npm install --no-audit --no-fund >/dev/null
        npm run build >/dev/null
    )
    echo "$STAMP" > "$BUILD/.stamp"
fi

STAGE="$BUILD/.stage"
rm -rf "$STAGE"
mkdir -p "$STAGE"
shopt -s nullglob
for f in "$PRE"/*.json.gz; do ln -s "$f" "$STAGE/"; done
for f in "$PRE"/*.json; do [[ -e "$f.gz" ]] || ln -s "$f" "$STAGE/"; done
shopt -u nullglob

exec bash "$DEPLOYER" \
    --viz-dir "$BUILD" --skip-build \
    --preprocess-dir "$STAGE" \
    --pipeline danaseq-mag-analysis \
    "${REST[@]}"
