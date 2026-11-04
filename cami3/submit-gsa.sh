#!/bin/bash
# Submit the scoreable baseline chain for one dataset on fir: for each read type,
# index the anonymous GSA, map every sample, build depths, then the kitchen-sink MAG
# run on the GSA. Needs the input_long/ and input_short/ staging from prep-reads.sh.
#   ./submit-gsa.sh toy    (20 samples)    ./submit-gsa.sh full   (36 samples)
set -euo pipefail
DS=$1
case "$DS" in toy) LAST=19 ;; full) LAST=35 ;; *) echo "dataset must be toy or full" >&2; exit 1 ;; esac
for t in long short; do
    idx=$(sbatch --parsable --job-name="gsaidx-$DS-$t" map-gsa.sh "$DS" "$t" index)
    map=$(sbatch --parsable --job-name="gsamap-$DS-$t" --array=0-$LAST --dependency=afterok:$idx map-gsa.sh "$DS" "$t" map)
    dep=$(sbatch --parsable --job-name="gsadep-$DS-$t" --dependency=afterok:$map map-gsa.sh "$DS" "$t" depths)
    mag=$(sbatch --parsable --job-name="mag-$DS-gsa-$t" --dependency=afterok:$dep mag.sh "$DS" "gsa_$t")
    echo "$DS gsa_$t: index $idx -> map $map -> depths $dep -> mag $mag"
done
