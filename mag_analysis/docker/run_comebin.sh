#!/usr/bin/env bash
# dana-tools wrapper for COMEBin in the container.
#
# run_comebin.sh resolves -a with realpath and writes <assembly>_lengths.txt and
# its marker-seed files beside the resolved file. A staged input resolves into
# a read-only bind mount, so when the assembly's directory is not writable the
# assembly is copied into the working directory first.
set -euo pipefail
ENV=/opt/conda/envs/dana-mag-comebin
args=()
while (( $# )); do
    case "$1" in
        -a) f=$(readlink -f "$2")
            if [ ! -w "$(dirname "$f")" ]; then
                cp "$f" "$PWD/comebin_input.fasta"
                f="$PWD/comebin_input.fasta"
            fi
            args+=(-a "$f"); shift 2 ;;
        *)  args+=("$1"); shift ;;
    esac
done
export PATH="$ENV/bin:$PATH"
exec "$ENV/bin/run_comebin.sh" "${args[@]}"
