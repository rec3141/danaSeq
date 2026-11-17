#!/bin/bash
# Kitchen-sink danaSeq MAG analysis (--all) of one CAMI III co-assembly on fir or grex.
#   sbatch --job-name=mag-toy-long --dependency=afterok:<asm job> mag.sh toy long
# Types: long, short (danaSeq assemblies), gsa_long, gsa_short (CAMI anonymous GSA
# mapped by map-gsa.sh, the scoreable baseline).
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
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -uo pipefail
DS=$1 TYPE=$2
. "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/cluster.sh"
C=$BASE/$DS
ASM=$C/assembly_$TYPE

# The assembly pipeline links assembly/* to /data/store/..., its in-container
# path; point those at the same files relative to the assembly dir.
for l in "$ASM"/assembly/*; do
    t=$(readlink "$l" 2>/dev/null) || continue
    case "$t" in /data/store/*) ln -sfn "../${t#/data/store/}" "$l" ;; esac
done
# nanopore_assembly writes assembly/assembly.fasta, mapping/depths.txt and
# mapping/*.sorted.bam. illumina_assembly (co-assembly) writes
# assembly/coassembly/coassembly.dedupe.fasta, mapping/coassembly/coassembly.depths.txt
# and mapping/<sample>/<sample>.sorted.bam; its BAMs and indexes are hard-linked into
# one flat directory (same filesystem, so the links resolve inside the container).
case "$TYPE" in
    long|gsa_long|gsa_short)
           FASTA=$ASM/assembly/assembly.fasta; DEPTHS=$ASM/mapping/depths.txt; BAMS=$ASM/mapping ;;
    short) FASTA=$ASM/assembly/coassembly/coassembly.dedupe.fasta
           DEPTHS=$ASM/mapping/coassembly/coassembly.depths.txt
           BAMS=$ASM/mapping_flat; mkdir -p "$BAMS"
           for b in "$ASM"/mapping/*/*.sorted.bam; do
               ln -f "$b" "$BAMS/" && ln -f "$b.bai" "$BAMS/" || { echo "[ERROR] cannot link $b" >&2; exit 1; }
           done ;;
    *) echo "[ERROR] type must be long or short, got $TYPE" >&2; exit 1 ;;
esac
for f in "$FASTA" "$DEPTHS"; do
    if [ ! -s "$f" ]; then
        echo "[ERROR] $f is missing or empty; the assembly did not finish" >&2
        exit 1
    fi
done
if [ -s "$ASM/mapping/unmapped_samples.txt" ]; then
    echo "[WARNING] $(wc -l < "$ASM/mapping/unmapped_samples.txt") samples are missing from depths.txt" >&2
fi
echo "[INFO] BAMs: $(ls "$BAMS"/*.sorted.bam | wc -l)  depth columns: $(head -1 "$DEPTHS" | tr "\t" "\n" | grep -cE "sorted(\.bam)?$")"
TNF=(); [ -s "$ASM/assembly/tnf.tsv" ] && TNF=(--tnf "$ASM/assembly/tnf.tsv")

container_setup
export NXF_SYNTAX_PARSER=v1
# grex must keep its work dir on /project; on fir the launcher's default
# (node-local SLURM_TMPDIR plus a store dir) applies.
WORKDIR=(); [ "$CLUSTER" = grex ] && WORKDIR=(--workdir "$C/work_mag_$TYPE")
# OpenBLAS in the image segfaults above 128 threads; tasks inherit this.
export OPENBLAS_NUM_THREADS=64
L=$C/launch_mag_$TYPE
mkdir -p "$L" && cd "$L" || exit 1

"$GENICE/danaSeq/mag_analysis/run-mag-analysis.sh" --container \
    --assembly "$FASTA" \
    --depths "$DEPTHS" \
    --bam_dir "$BAMS/" \
    "${TNF[@]}" \
    --outdir "$C/mag_$TYPE" \
    --db_dir "$DB" \
    "${WORKDIR[@]}" \
    --all --run_integronfinder false \
    --binning_cpus 32 --binning_memory "160 GB" \
    --bakta_cpus 64 --bakta_memory "300 GB" \
    --medium_cpus 32 --medium_memory "96 GB"
