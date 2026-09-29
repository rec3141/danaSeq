#!/bin/sh
# Pin the PyTorch build in conda env specs for one image variant.
#   torch_variant.sh cpu|cuda env.yml...
# Drops any pytorch / pytorch-cpu / pytorch-gpu line and appends the pin, so
# the env files themselves stay variant-neutral for conda installs.
#   cpu:  pytorch=*=cpu*
#   cuda: pytorch=*=cuda* with cuda-version=12.9. The packages carry their own
#         CUDA runtime, so hosts need only an NVIDIA driver >= 525; 12.9 is the
#         last CUDA that still builds for Volta (V100).
set -eu
variant=$1; shift
case "$variant" in
    cpu)  pins='  - pytorch=*=cpu*' ;;
    cuda) pins='  - pytorch=*=cuda*
  - cuda-version=12.9' ;;
    *) echo "torch_variant.sh: variant must be cpu or cuda, got '$variant'" >&2; exit 1 ;;
esac
for f in "$@"; do
    grep -vE '^[[:space:]]*-[[:space:]]*pytorch(-cpu|-gpu)?([[:space:]=<>]|$)' "$f" > "$f.tmp"
    # a trailing newline may be missing; the pins must start on their own line
    [ -z "$(tail -c1 "$f.tmp")" ] || echo >> "$f.tmp"
    printf '%s\n' "$pins" >> "$f.tmp"
    mv "$f.tmp" "$f"
done
