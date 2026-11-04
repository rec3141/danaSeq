#!/bin/bash
# Map one dataset's reads onto CAMI's provided anonymous pooled GSA. The result
# is the scoreable baseline input: bins use the official contig IDs, which AMBER needs.
#   sbatch --job-name=gsaidx-toy-long  map-gsa.sh toy long index    (stage the uncompressed GSA)
#   sbatch --job-name=gsamap-toy-long  --array=0-19 --dependency=afterok:<index> map-gsa.sh toy long map
#   sbatch --job-name=gsadep-toy-long  --dependency=afterok:<map> map-gsa.sh toy long depths
# (submit-gsa.sh chains these.) Output layout matches nanopore_assembly, so mag.sh
# runs on it unchanged with type gsa_long or gsa_short:
#   $C/assembly_gsa_<type>/assembly/assembly.fasta    anonymous GSA, uncompressed
#   $C/assembly_gsa_<type>/mapping/sample_N.sorted.bam (+ .bai)
#   $C/assembly_gsa_<type>/mapping/depths.txt          MetaBAT2 format
# Long reads map with minimap2 -x map-ont, short with -x sr (R1/R2 from prep-reads.sh).
# Depth policy: long reads use CoverM at 80% identity (the nanopore_assembly recipe),
# short reads use jgi_summarize_bam_contig_depths (the illumina_assembly recipe).
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
    # Mapping indexes the reference on the fly, with --split-prefix, as the pipeline does.
    echo "reference contigs: $(grep -c '^>' "$A/assembly/assembly.fasta")" ;;
map)
    # A sample is done only when its .done marker exists. The marker is written
    # after the BAM, its index, its @SQ count and its alignment count all check
    # out, so a job killed between renaming the BAM and indexing it reruns.
    N=${SLURM_ARRAY_TASK_ID:?run as an array job}
    out=$A/mapping/sample_$N.sorted.bam
    done_mark=$A/mapping/.sample_$N.done
    if [ -e "$done_mark" ]; then echo "sample $N already mapped"; exit 0; fi
    if [ "$TYPE" = long ]; then reads="$C/input_long/sample_$N.fq.gz"
    else reads="$C/input_short/sample_${N}_R1_001.fastq.gz $C/input_short/sample_${N}_R2_001.fastq.gz"; fi
    MM=$(( SLURM_CPUS_PER_TASK - 4 )); [ "$MM" -lt 1 ] && MM=1
    # The inner shell sets pipefail too, so a mapper that dies after writing some
    # records fails the task instead of leaving a valid-looking truncated BAM.
    # --split-prefix keeps @SQ lines when minimap2 has to use a multi-part index
    # (see nanopore_assembly/modules/mapping.nf). -F 0x104 drops unmapped and
    # secondary reads but keeps supplementary alignments, as the pipeline does.
    run bash -o pipefail -c "minimap2 -a -x $PRESET --secondary=no -t $MM --split-prefix $T/split_$N \
        $A/assembly/assembly.fasta $reads \
        | samtools view -u -F 0x104 - | samtools sort -@ 4 -m 2G -T $T/sort_$N -o $T/sample_$N.bam -"
    run samtools index "$T/sample_$N.bam"
    REF=$(grep -c '^>' "$A/assembly/assembly.fasta")
    SQ=$(run samtools view -H "$T/sample_$N.bam" | grep -c '^@SQ' || true)
    ALN=$(run samtools view -c "$T/sample_$N.bam")
    echo "sample $N: reference contigs $REF, BAM @SQ $SQ, alignments $ALN"
    if [ "$SQ" -ne "$REF" ] || [ "$ALN" -eq 0 ]; then echo "[ERROR] sample $N BAM does not match the reference or is empty" >&2; exit 1; fi
    cp "$T/sample_$N.bam" "$out.part"; cp "$T/sample_$N.bam.bai" "$out.bai.part"
    run samtools quickcheck "$out.part"
    mv "$out.part" "$out"; mv "$out.bai.part" "$out.bai"
    touch "$done_mark" ;;
depths)
    # Every expected sample must have its .done marker before depths are built.
    LAST=$( [ "$DS" = toy ] && echo 19 || echo 35 )
    for n in $(seq 0 "$LAST"); do
        [ -e "$A/mapping/.sample_$n.done" ] || { echo "[ERROR] sample $n is not mapped" >&2; exit 1; }
    done
    if [ "$TYPE" = long ]; then
        # The nanopore_assembly recipe: CoverM mean/variance at 80% identity,
        # reshaped into the MetaBAT table. jgi_summarize's 97% default drops most
        # raw-ONT coverage.
        run bash -o pipefail -c "coverm contig -b $A/mapping/sample_*.sorted.bam --methods length mean variance \
            --min-read-percent-identity 80 --min-read-aligned-percent 0 --min-covered-fraction 0 \
            --contig-end-exclusion 75 --threads $SLURM_CPUS_PER_TASK | /pipeline/bin/jgi_depth_table.awk" \
            > "$A/mapping/depths.txt.part"
    else
        # The illumina_assembly recipe: jgi_summarize_bam_contig_depths at its default identity.
        run jgi_summarize_bam_contig_depths --outputDepth "$A/mapping/depths.txt.part" "$A"/mapping/sample_*.sorted.bam
    fi
    # jgi_summarize names columns sample_N.sorted (no .bam); the CoverM converter keeps .bam.
    COLS=$(head -1 "$A/mapping/depths.txt.part" | tr '\t' '\n' | grep -cE 'sorted(\.bam)?$' || true)
    if [ "$COLS" -ne $(( LAST + 1 )) ]; then echo "[ERROR] depths has $COLS sample columns, expected $(( LAST + 1 ))" >&2; exit 1; fi
    mv "$A/mapping/depths.txt.part" "$A/mapping/depths.txt"
    echo "depths: $COLS sample columns, $(( $(wc -l < "$A/mapping/depths.txt") - 1 )) contigs" ;;
*) echo "step must be index, map or depths" >&2; exit 1 ;;
esac
