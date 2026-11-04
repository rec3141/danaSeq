#!/bin/bash
# Submit the CAMI III chain for one dataset on fir: read prep (array), then the
# long- and short-read co-assemblies, then a kitchen-sink MAG run on each.
#   ./submit.sh toy    (20 samples)
#   ./submit.sh full   (36 samples)
# Run from the directory that holds these scripts.
set -euo pipefail
DS=$1
case "$DS" in toy) LAST=19 ;; full) LAST=35 ;; *) echo "dataset must be toy or full" >&2; exit 1 ;; esac
mkdir -p /scratch/rec3141/cami3/logs
prep=$(sbatch --parsable --job-name="prep-$DS" --array=0-$LAST prep-reads.sh "$DS")
for t in long short; do
    asm=$(sbatch --parsable --job-name="asm-$DS-$t" --dependency=afterok:$prep asm.sh "$DS" "$t")
    mag=$(sbatch --parsable --job-name="mag-$DS-$t" --dependency=afterok:$asm mag.sh "$DS" "$t")
    echo "$DS $t: prep $prep -> asm $asm -> mag $mag"
done
