#!/bin/bash
# Standard danaSeq co-assembly of one CAMI III read set on fir: every sample of
# the dataset co-assembled, then each sample mapped back for per-contig depths.
#   sbatch --job-name=asm-toy-long  asm.sh toy  long     (nanopore_assembly, Flye)
#   sbatch --job-name=asm-toy-short asm.sh toy  short    (illumina_assembly, co-assembly)
# Inputs come from prep-reads.sh. To continue after a timeout or failure,
# resubmit with --resume as a third argument. A first run must not pass it: with
# no run history, the launcher adopts ~/.nextflow/history and resumes an
# unrelated session.
#SBATCH --account=def-rec3141
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=0
#SBATCH --output=/scratch/rec3141/cami3/logs/%x-%j.out
#SBATCH --error=/scratch/rec3141/cami3/logs/%x-%j.err
set -uo pipefail
DS=$1 TYPE=$2 RESUME=(${3:+--resume})
C=/scratch/rec3141/cami3/$DS
GENICE=/home/rec3141/GENICE
# Nanopore removes human reads with minimap2 against the FASTA; Illumina runs
# removehuman.sh, which takes the BBMap index directory (holding ref/).
HUMAN_FA=/home/rec3141/scratch/databases/human_ref/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz
HUMAN_BBMAP=/home/rec3141/scratch/databases/human_ref
module load apptainer
export APPTAINER_CACHEDIR=/home/rec3141/scratch/apptainer_cache
export APPTAINER_TMPDIR=${SLURM_TMPDIR:-/tmp}
# The image's nextflow.config needs Nextflow's legacy config parser.
export NXF_SYNTAX_PARSER=v1
# Each dataset/read type launches from its own directory, so its .nextflow
# history is separate: the illumina launcher always passes -resume, which with no
# session ID resumes the last run in the launch directory. Work dirs and history
# stay on /scratch so a resubmitted job can resume.
L=$C/launch_assembly_$TYPE
mkdir -p "$L" && cd "$L" || exit 1

case "$TYPE" in
    long)  "$GENICE/danaSeq/nanopore_assembly/run-nanopore-assembly.sh" --container \
               --input "$C/input_long" \
               --outdir "$C/assembly_long" \
               --workdir "$C/work_assembly_long" \
               --human_ref "$HUMAN_FA" \
               --read_type nano-hq \
               --assembly_cpus 128 --assembly_memory '700 GB' \
               "${RESUME[@]}" ;;
    short) "$GENICE/danaSeq/illumina_assembly/run-illumina-assembly.sh" --container \
               --input "$C/input_short" \
               --outdir "$C/assembly_short" \
               --workdir "$C/work_assembly_short" \
               --human_ref "$HUMAN_BBMAP" \
               --coassembly true \
               --assembly_cpus 128 --assembly_memory '700 GB' \
               "${RESUME[@]}" ;;
    *) echo "[ERROR] type must be long or short, got $TYPE" >&2; exit 1 ;;
esac
