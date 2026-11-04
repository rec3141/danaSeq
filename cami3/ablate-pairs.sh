#!/bin/bash
# Matched two-sample depth control: do a subject's two samples separate strains
# better or worse than two samples from different subjects?
#   sbatch --job-name=pairs-toy-long --array=0-19 ablate-pairs.sh toy long
# Array tasks 0-9 bin with the two samples of subject S_(i+1) (within arm).
# Tasks 10-19 bin with the first sample of S_(j+1) and the second sample of
# S_(j+2 mod 10), j = i-10 (between arm), so every sample appears exactly once
# in each arm. "First" and "second" are the lower and higher sample numbers;
# the mapping file gives no collection order, and the design needs none.
# Everything else is held at the ablate-depth.sh settings (MetaBAT2, seed 1,
# -m 1500, the same contigs). Writes
# ablation_gsa_<type>/pairs/binning/<config>/<config>_bins.tsv for score_run.py.
#SBATCH --account=def-rec3141
#SBATCH --time=0-04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
DS=$1 TYPE=$2
. "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/cluster.sh"
C=$BASE/$DS
A=$C/assembly_gsa_$TYPE
O=$C/ablation_gsa_$TYPE/pairs
SIF=$GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif
container_setup
mkdir -p "$O/depths" "$O/binning"
DEP=$A/mapping/depths.txt
MAP=$C/sample_subject_mapping.tsv
[ -s "$MAP" ] || MAP=$C/short/sample_subject_mapping.tsv
[ -s "$DEP" ] && [ -s "$MAP" ] || { echo "[ERROR] need $DEP and the sample-subject mapping" >&2; exit 1; }
. "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/ablate-lib.sh"

# samples_of S<n>: that subject's samples, lower sample number first.
samples_of() { awk -F'\t' -v s="$1" '$2==s{sub(/^sample_/,"",$1); print $1}' "$MAP" | sort -n | sed 's/^/sample_/'; }
for s in $(seq 1 10); do
    [ "$(samples_of S$s | wc -l)" -eq 2 ] || { echo "[ERROR] S$s does not have exactly 2 samples" >&2; exit 1; }
done

i=${SLURM_ARRAY_TASK_ID:?run as an array job, --array=0-19}
if [ "$i" -lt 10 ]; then
    a=$((i + 1)); cfg=within_S$a
    set -- $(samples_of S$a)
elif [ "$i" -lt 20 ]; then
    j=$((i - 10)); a=$((j + 1)); b=$(( (j + 1) % 10 + 1 )); cfg=between_S${a}_S$b
    set -- "$(samples_of S$a | head -1)" "$(samples_of S$b | tail -1)"
else
    echo "[ERROR] array index $i out of range 0-19" >&2; exit 1
fi
echo "$cfg: $*"
subset_depth "$O/depths/$cfg.txt" "$@"
run_metabat "$cfg" "$O/depths/$cfg.txt"
