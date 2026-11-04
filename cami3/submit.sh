#!/bin/bash
# Submit the CAMI III chain for one dataset on fir or grex: read prep (array), then the
# long- and short-read co-assemblies, then a kitchen-sink MAG run on each.
#   ./submit.sh toy    (20 samples)
#   ./submit.sh full   (36 samples)
# Run on fir or grex (cluster.sh picks the paths); the scripts must sit in $BASE.
set -euo pipefail
DS=$1
. "$(dirname "$0")/cluster.sh"
cd "$BASE"
mkdir -p logs
BIG=(); [ -n "$BIG_PARTITION" ] && BIG=(--partition="$BIG_PARTITION")
case "$DS" in toy) LAST=19 ;; full) LAST=35 ;; *) echo "dataset must be toy or full" >&2; exit 1 ;; esac
prep=$(sbatch --parsable --job-name="prep-$DS" --array=0-$LAST prep-reads.sh "$DS")
for t in long short; do
    asm=$(sbatch --parsable --job-name="asm-$DS-$t" "${BIG[@]}" --dependency=afterok:$prep asm.sh "$DS" "$t")
    mag=$(sbatch --parsable --job-name="mag-$DS-$t" "${BIG[@]}" --dependency=afterok:$asm mag.sh "$DS" "$t")
    echo "$DS $t: prep $prep -> asm $asm -> mag $mag"
done
