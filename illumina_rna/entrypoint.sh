#!/bin/sh
# Entrypoint wrapper for Docker --user support.
# Ensures the Nextflow output directory exists and is writable before launch.

outdir="results"
for arg in "$@"; do
    case "$prev" in
        --outdir) outdir="$arg" ;;
    esac
    prev="$arg"
done

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
