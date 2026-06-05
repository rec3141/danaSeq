#!/usr/bin/env bash
# Thin wrapper that pins the nanopore_live SPA dir + pipeline slug, then
# forwards everything else to tools/microscape_deploy.sh. Existing per-run
# hooks under <outdir>/deploy.sh keep working unchanged.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/../../tools/microscape_deploy.sh" \
    --viz-dir "$SCRIPT_DIR" \
    --pipeline danaseq-nanopore-live \
    "$@"
