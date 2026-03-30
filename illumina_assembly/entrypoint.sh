#!/bin/sh
# Entrypoint wrapper for Docker --user support.
# Ensures the Nextflow output directory exists and is writable before
# launching Nextflow. When Docker creates a bind-mount target that doesn't
# exist on the host, it creates it as root:root -- this script detects that
# and exits with a clear error instead of a cryptic publishDir failure.

# Parse --outdir from arguments (default: "results")
outdir="results"
for arg in "$@"; do
    case "$prev" in
        --outdir) outdir="$arg" ;;
    esac
    prev="$arg"
done

# If outdir looks like an absolute path (bind-mounted volume), verify it's writable
if [ "${outdir#/}" != "$outdir" ]; then
    mkdir -p "$outdir" 2>/dev/null
    if [ ! -w "$outdir" ]; then
        echo "[ERROR] Output directory is not writable: $outdir" >&2
        echo "        When using --user, create the output directory on the host first:" >&2
        echo "          mkdir -p /path/to/output" >&2
        echo "        Then mount it: -v /path/to/output:$outdir" >&2
        exit 1
    fi
fi

exec nextflow -c /pipeline/docker.config "$@"
