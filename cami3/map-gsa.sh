#!/bin/bash
# Map one dataset's reads onto CAMI's provided anonymous pooled GSA. The result
# is the scoreable baseline input: bins use the official contig IDs, which AMBER needs.
#   sbatch --job-name=gsaidx-toy-long  map-gsa.sh toy long index
#   sbatch --job-name=gsamap-toy-long  --array=0-19 --dependency=afterok:<index> map-gsa.sh toy long map
#   sbatch --job-name=gsadep-toy-long  --dependency=afterok:<map> map-gsa.sh toy long depths
# (submit-gsa.sh chains these.) Output layout matches nanopore_assembly, so mag.sh
# runs on it unchanged with type gsa_long or gsa_short:
#   $C/assembly_gsa_<type>/assembly/assembly.fasta    anonymous GSA, uncompressed
#   $C/assembly_gsa_<type>/mapping/sample_N.sorted.bam (+ .bai)
#   $C/assembly_gsa_<type>/mapping/depths.txt          MetaBAT2 format
# Long reads map with minimap2 -x map-ont, short with -x sr (R1/R2 from prep-reads.sh).
# Each BAM is written to a .part name and renamed only after samtools quickcheck,
# so a rerun redoes any sample whose BAM is missing or bad.
#SBATCH --account=def-rec3141
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=/scratch/rec3141/cami3/logs/%x-%A_%a.out
#SBATCH --error=/scratch/rec3141/cami3/logs/%x-%A_%a.err
set -euo pipefail
DS=$1 TYPE=$2 STEP=$3
C=/scratch/rec3141/cami3/$DS
A=$C/assembly_gsa_$TYPE
SIF=/home/rec3141/GENICE/danaSeq/nanopore_assembly/.danaseq-nanopore-assembly.sif
module load apptainer
T=${SLURM_TMPDIR:-/tmp}
# Node-local TMP is bound explicitly; samtools sort writes its temp files there.
run() { apptainer exec --bind /scratch --bind "$T" "$SIF" "$@"; }
mkdir -p "$A/assembly" "$A/mapping"
case "$TYPE" in long) PRESET=map-ont ;; short) PRESET=sr ;; *) echo "type must be long or short" >&2; exit 1 ;; esac

case "$STEP" in
index)
    if [ ! -s "$A/assembly/assembly.fasta" ]; then
        gzip -dc "$C/$TYPE/anonymous_gsa_pooled.fasta.gz" > "$A/assembly/assembly.fasta.part"
        mv "$A/assembly/assembly.fasta.part" "$A/assembly/assembly.fasta"
    fi
    run minimap2 -t "$SLURM_CPUS_PER_TASK" -x "$PRESET" -d "$A/assembly/assembly.$PRESET.mmi.part" "$A/assembly/assembly.fasta"
    mv "$A/assembly/assembly.$PRESET.mmi.part" "$A/assembly/assembly.$PRESET.mmi" ;;
map)
    N=${SLURM_ARRAY_TASK_ID:?run as an array job}
    out=$A/mapping/sample_$N.sorted.bam
    if [ -s "$out" ] && run samtools quickcheck "$out"; then echo "sample $N already mapped"; exit 0; fi
    if [ "$TYPE" = long ]; then reads=("$C/input_long/sample_$N.fq.gz")
    else reads=("$C/input_short/sample_${N}_R1_001.fastq.gz" "$C/input_short/sample_${N}_R2_001.fastq.gz"); fi
    run bash -c "minimap2 -t $SLURM_CPUS_PER_TASK -ax $PRESET --secondary=no $A/assembly/assembly.$PRESET.mmi ${reads[*]} \
        | samtools sort -@ 8 -m 2G -T $T/sort -o $T/sample_$N.bam -"
    cp "$T/sample_$N.bam" "$out.part"
    run samtools quickcheck "$out.part"
    mv "$out.part" "$out"
    run samtools index "$out" ;;
depths)
    run jgi_summarize_bam_contig_depths --outputDepth "$A/mapping/depths.txt.part" "$A"/mapping/sample_*.sorted.bam
    mv "$A/mapping/depths.txt.part" "$A/mapping/depths.txt"
    echo "depth columns: $(head -1 "$A/mapping/depths.txt" | tr '\t' '\n' | grep -c 'sorted.bam$')" ;;
*) echo "step must be index, map or depths" >&2; exit 1 ;;
esac
