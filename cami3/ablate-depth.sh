#!/bin/bash
# How much does multi-sample coverage help? MetaBAT2 bins one GSA with depth tables
# restricted to different sample subsets. Everything else is held fixed: the same
# contigs, the same seed (1) and the same minimum contig length (1500).
#   sbatch --job-name=ablate-toy-long ablate-depth.sh toy long
# Reads assembly_gsa_<type>/ from map-gsa.sh. Writes
# ablation_gsa_<type>/binning/<config>/<config>_bins.tsv (contig<TAB>bin, unbinned
# contigs left out), the layout cami3-eval/score_run.py scores.
# Configs: TNF only (metabat2 with no -a), the composition-only floor; one sample
# (sample_0); one subject's samples (S1, from sample_subject_mapping.tsv; the toy
# has 10 subjects with 2 samples each); every sample.
#SBATCH --account=def-rec3141
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail
DS=$1 TYPE=$2
. "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/cluster.sh"
C=$BASE/$DS
A=$C/assembly_gsa_$TYPE
O=$C/ablation_gsa_$TYPE
SIF=$GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif
container_setup
mkdir -p "$O/depths" "$O/binning"
DEP=$A/mapping/depths.txt
MAP=$C/sample_subject_mapping.tsv
[ -s "$MAP" ] || MAP=$C/short/sample_subject_mapping.tsv
[ -s "$DEP" ] && [ -s "$MAP" ] || { echo "[ERROR] need $DEP and the sample-subject mapping" >&2; exit 1; }

. "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/ablate-lib.sh"

S1=$(awk -F'\t' '$2=="S1"{print $1}' "$MAP" | tr '\n' ' ')
ALL=$(awk -F'\t' 'NR>1{print $1}' "$MAP" | tr '\n' ' ')
subset_depth "$O/depths/one_sample.txt" sample_0
subset_depth "$O/depths/one_subject.txt" $S1
# all_samples is the source table itself, copied exactly. Recomputing it would
# re-round the totals (the converter rounds totals and sample means separately to
# four decimals). subset_depth still runs on the full sample set first, as a check
# that every mapped sample has its column; awk exits 3 otherwise.
subset_depth "$O/depths/all_samples.check" $ALL
rm -f "$O/depths/all_samples.check"
cp "$DEP" "$O/depths/all_samples.txt"
echo "S1 samples: $S1"

run_metabat tnf_only
run_metabat one_sample "$O/depths/one_sample.txt"
run_metabat one_subject "$O/depths/one_subject.txt"
run_metabat all_samples "$O/depths/all_samples.txt"
