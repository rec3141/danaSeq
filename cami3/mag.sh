#!/bin/bash
# Kitchen-sink danaSeq MAG analysis (--all) of one CAMI III co-assembly on fir.
#   sbatch --job-name=mag-toy-long --dependency=afterok:<asm job> mag.sh toy long
# Reads assembly_<type>/ from asm.sh and writes mag_<type>/. The launcher keeps
# its work dir on node-local $SLURM_TMPDIR and stores finished task outputs, so
# a resubmission picks up where the last job stopped.
# IntegronFinder is off: the image runs it on the whole assembly at once, which
# holds two files per contig until it finishes (danaSeq #69, #70).
#SBATCH --account=def-rec3141
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=0
#SBATCH --output=/scratch/rec3141/cami3/logs/%x-%j.out
#SBATCH --error=/scratch/rec3141/cami3/logs/%x-%j.err
set -uo pipefail
DS=$1 TYPE=$2
C=/scratch/rec3141/cami3/$DS
ASM=$C/assembly_$TYPE
GENICE=/home/rec3141/GENICE

# The assembly pipeline links assembly/* to /data/store/..., its in-container
# path; point those at the same files relative to the assembly dir.
for l in "$ASM"/assembly/*; do
    t=$(readlink "$l" 2>/dev/null) || continue
    case "$t" in /data/store/*) ln -sfn "../${t#/data/store/}" "$l" ;; esac
done
for f in assembly/assembly.fasta mapping/depths.txt; do
    if [ ! -s "$ASM/$f" ]; then
        echo "[ERROR] $ASM/$f is missing or empty; the assembly did not finish" >&2
        exit 1
    fi
done
if [ -s "$ASM/mapping/unmapped_samples.txt" ]; then
    echo "[WARNING] $(wc -l < "$ASM/mapping/unmapped_samples.txt") samples are missing from depths.txt" >&2
fi
echo "[INFO] BAMs: $(ls "$ASM"/mapping/*.sorted.bam | wc -l)  depth columns: $(head -1 "$ASM/mapping/depths.txt" | tr "\t" "\n" | grep -c "sorted.bam$")"
TNF=(); [ -s "$ASM/assembly/tnf.tsv" ] && TNF=(--tnf "$ASM/assembly/tnf.tsv")

module load apptainer
export APPTAINER_CACHEDIR=/home/rec3141/scratch/apptainer_cache
export APPTAINER_TMPDIR=${SLURM_TMPDIR:-/tmp}
export NXF_SYNTAX_PARSER=v1
# OpenBLAS in the image segfaults above 128 threads; tasks inherit this.
export OPENBLAS_NUM_THREADS=64

"$GENICE/danaSeq/mag_analysis/run-mag-analysis.sh" --container \
    --assembly "$ASM/assembly/assembly.fasta" \
    --depths "$ASM/mapping/depths.txt" \
    --bam_dir "$ASM/mapping/" \
    "${TNF[@]}" \
    --outdir "$C/mag_$TYPE" \
    --db_dir /home/rec3141/scratch/databases \
    --all --run_integronfinder false \
    --binning_cpus 32 --binning_memory "160 GB" \
    --bakta_cpus 64 --bakta_memory "300 GB" \
    --medium_cpus 32 --medium_memory "96 GB"
