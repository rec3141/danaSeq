#!/bin/bash
# Stage CAMI III human-gut reads into the input layout the danaSeq assembly
# pipelines expect, one sample per array task.
#   sbatch --array=0-19 prep-reads.sh toy      (toy: 20 samples)
#   sbatch --array=0-35 prep-reads.sh full     (full: 36 samples)
# Long reads -> $C/input_long/sample_N.fq.gz (nanopore_assembly flat mode).
# Short reads are interleaved pairs (@SxRy/1, @SxRy/2) -> split into
#   $C/input_short/sample_N_R1_001.fastq.gz and _R2_001.fastq.gz (illumina_assembly).
# The toy reads come in sample_N_reads.tar.gz (anonymous_reads.fq.gz inside);
# the full reads come as cami3_human_gut_{long,short}_read_sample_N.fq.gz.
# The split checks that every R1 and R2 name pairs up and fails the task otherwise.
# Under set -e, any failing step (extraction, pairing, gzip, copy) fails the task.
#SBATCH --account=def-rec3141
#SBATCH --time=0-06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
set -euo pipefail
DS=$1
N=${SLURM_ARRAY_TASK_ID:?run as an array job}
. "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/cluster.sh"
C=$BASE/$DS
mkdir -p "$C/input_long" "$C/input_short"
T=${SLURM_TMPDIR:-/tmp}

src() {  # print the gzipped FASTQ for read type $1 (long|short) on stdout
    if [ "$DS" = toy ]; then
        tar xzf "$C/$1/sample_${N}_reads.tar.gz" -O --wildcards '*anonymous_reads.fq.gz'
    else
        cat "$C/$1/cami3_human_gut_$1_read_sample_${N}.fq.gz"
    fi
}

# Each output is copied to a .part name on /scratch, checked with gzip -t, then
# renamed, so a killed copy never leaves a file under its final name. The
# sample's .done marker is written only after all three outputs are in place;
# a rerun skips the sample only when that marker exists.
done_mark=$C/input_short/.sample_$N.done
if [ -e "$done_mark" ]; then echo "sample $N already staged"; exit 0; fi
stage() {  # stage <tmp file> <final path>
    cp "$1" "$2.part" && gzip -t "$2.part" && mv -f "$2.part" "$2"
}

out=$C/input_long/sample_$N.fq.gz
src long > "$T/long.fq.gz"
stage "$T/long.fq.gz" "$out"

r1=$C/input_short/sample_${N}_R1_001.fastq.gz
r2=$C/input_short/sample_${N}_R2_001.fastq.gz
src short | gzip -dc | paste - - - - - - - - | awk -F'\t' -v r1="$T/r1" -v r2="$T/r2" '
    { a=$1; b=$5; sub(/\/1$/,"",a); sub(/\/2$/,"",b)
      if (a != b || $1 !~ /\/1$/ || $5 !~ /\/2$/) { print "unpaired record " NR ": " $1 " " $5 > "/dev/stderr"; exit 3 }
      print $1"\n"$2"\n"$3"\n"$4 > r1; print $5"\n"$6"\n"$7"\n"$8 > r2 }'
gzip -1 "$T/r1" & p1=$!
gzip -1 "$T/r2" & p2=$!
wait "$p1"; wait "$p2"
stage "$T/r1.gz" "$r1"
stage "$T/r2.gz" "$r2"
touch "$done_mark"
echo "sample $N ready: $(du -h "$out" "$r1" "$r2" | cut -f1 | tr '\n' ' ')"
