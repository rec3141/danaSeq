#!/bin/bash
# Sample-count curve: MetaBAT2 on depth tables restricted to fixed sample subsets.
#   sbatch --job-name=kcurve-toy-long --array=0-50 ablate-kcurve.sh toy long
# Array task i bins with row i of kcurve_manifest.tsv (written in advance by
# make-kcurve-manifest.py: 'spread' and 'paired' families at k = 2..10, five
# draws each, plus all 20 samples). Settings match ablate-depth.sh (seed 1,
# -m 1500, the same contigs). Writes
# ablation_gsa_<type>/kcurve/binning/<config>/<config>_bins.tsv for score_run.py.
#SBATCH --account=def-rec3141
#SBATCH --time=0-02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
DS=$1 TYPE=$2
. "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/cluster.sh"
C=$BASE/$DS
A=$C/assembly_gsa_$TYPE
O=$C/ablation_gsa_$TYPE/kcurve
SIF=$GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif
MANIFEST=${SLURM_SUBMIT_DIR:-$(dirname "$0")}/kcurve_manifest.tsv
container_setup
mkdir -p "$O/depths" "$O/binning"
DEP=$A/mapping/depths.txt
[ -s "$DEP" ] && [ -s "$MANIFEST" ] || { echo "[ERROR] need $DEP and $MANIFEST" >&2; exit 1; }
. "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/ablate-lib.sh"

i=${SLURM_ARRAY_TASK_ID:?run as an array job}
row=$(awk -F'\t' -v r=$((i + 2)) 'NR==r' "$MANIFEST")
[ -n "$row" ] || { echo "[ERROR] no manifest row for task $i" >&2; exit 1; }
cfg=$(cut -f1 <<<"$row"); samples=$(cut -f6 <<<"$row" | tr ',' ' ')
echo "$cfg: $samples"
subset_depth "$O/depths/$cfg.txt" $samples
run_metabat "$cfg" "$O/depths/$cfg.txt"
rm -f "$O/depths/$cfg.txt"
