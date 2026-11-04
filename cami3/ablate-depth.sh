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

# subset_depth <out> <sample names...>: keep the given samples' mean/variance column
# pairs, in source order, and recompute totalAvgDepth the way the source table does.
# jgi_summarize (short reads) sums the per-sample means; nanopore's
# jgi_depth_table.awk (long reads) averages them.
case "$TYPE" in long) TOTAL=mean ;; short) TOTAL=sum ;; *) echo "type must be long or short" >&2; exit 1 ;; esac
subset_depth() {
    local out=$1; shift
    awk -F'\t' -v keep="$*" -v total="$TOTAL" 'BEGIN{OFS="\t"; n=split(keep,k," "); for(i=1;i<=n;i++) want[k[i]]=1}
        NR==1{ printf "%s\t%s\t%s", $1,$2,$3
               for(i=4;i<=NF;i+=2){ s=$i; sub(/\.sorted(\.bam)?$/,"",s); if(s in want){cols[++m]=i; printf "\t%s\t%s",$i,$(i+1)} }
               print ""; if(m!=n){print "[ERROR] matched "m" of "n" samples" > "/dev/stderr"; exit 3}; next }
        { t=0; for(j=1;j<=m;j++) t+=$(cols[j]); if(total=="mean") t/=m; printf "%s\t%s\t%s", $1,$2,t
          for(j=1;j<=m;j++) printf "\t%s\t%s", $(cols[j]), $(cols[j]+1); print "" }' "$DEP" > "$out"
}

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

run_metabat() {  # run_metabat <config> [depth file]
    local cfg=$1 dep=${2:-}
    local d=$O/binning/$cfg; mkdir -p "$d"
    rm -f "$d"/bin.*.fa
    cexec "$SIF" metabat2 -i "$A/assembly/assembly.fasta" ${dep:+-a "$dep"} -m 1500 --seed 1 \
        -t "$SLURM_CPUS_PER_TASK" -o "$d/bin"
    # The bin FASTAs are the authoritative membership: this MetaBAT2 build's
    # --saveCls matrix holds internal cluster IDs (including clusters below
    # minClsSize), not bin numbers. Keep only the contig<TAB>bin table.
    local f n
    for f in "$d"/bin.*.fa; do
        [ -e "$f" ] || continue
        n=${f##*/bin.}; n=${n%.fa}
        awk -v b="$(printf '%s_%04d' "$cfg" "$n")" '/^>/{sub(/^>/,""); sub(/[ \t].*/,""); print $0"\t"b}' "$f"
    done > "$d/${cfg}_bins.tsv"
    rm -f "$d"/bin.*.fa
    [ -s "$d/${cfg}_bins.tsv" ] || { echo "[ERROR] $cfg: no bins written" >&2; exit 1; }
    echo "$cfg: $(cut -f2 "$d/${cfg}_bins.tsv" | sort -u | wc -l) bins, $(wc -l < "$d/${cfg}_bins.tsv") contigs binned"
}
run_metabat tnf_only
run_metabat one_sample "$O/depths/one_sample.txt"
run_metabat one_subject "$O/depths/one_subject.txt"
run_metabat all_samples "$O/depths/all_samples.txt"
