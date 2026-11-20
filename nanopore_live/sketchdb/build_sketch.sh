#!/usr/bin/env bash
# Build a GTDB sketch database and taxonomy tree for the sendsketch TaxServer.
#
#   build_sketch.sh <gtdb_genomes_reps_rNNN.tar.gz> <taxonomy dir> <workdir>
#
# <taxonomy dir> holds the release's bac120_taxonomy_rNNN.tsv.gz and
# ar53_taxonomy_rNNN.tsv.gz. Writes to <workdir>:
#   gtdb_rNNN.taxtree.gz      BBTools tree of the GTDB taxonomy (tree= on the server)
#   gtdb_rNNN_k2432.sketch    one k=32,24 sketch per representative genome, each
#                             carrying its GTDB species taxid and its accession
#   taxdump/                  the nodes/names dmp files and accession2taxid.tsv
#
# Finished sketch batches are kept, so a rerun resumes where it stopped.
# BBTools comes from sketch.sh on PATH, or from $BBTOOLS_CP (its current/ dir).
# Run the server against the result with:
#   java -Xmx<heap> -cp <bbtools>/current tax.TaxServer port=3068 k=32,24 \
#       index=auto tree=gtdb_rNNN.taxtree.gz gtdb_rNNN_k2432.sketch
set -euo pipefail
TARBALL=$(readlink -f "$1") TAXDIR=$(readlink -f "$2") WORK=$3
HERE=$(cd "$(dirname "$0")" && pwd)
REL=$(basename "$TARBALL" | grep -o 'r[0-9]\+' | head -1)
[ -n "$REL" ] || { echo "cannot tell the GTDB release from $TARBALL" >&2; exit 1; }
CP=${BBTOOLS_CP:-$(dirname "$(readlink -f "$(command -v sketch.sh)")")/current}
BATCH=${BATCH:-2000}   # genomes per SketchMaker run
JOBS=${JOBS:-6}        # concurrent SketchMaker runs, 5 threads and 12 GB heap each
mkdir -p "$WORK/genomes" "$WORK/batches"
cd "$WORK"

# Taxonomy tree and accession -> taxid map
python3 "$HERE/build_gtdb_taxdump.py" taxdump \
    "$TAXDIR/bac120_taxonomy_$REL.tsv.gz" "$TAXDIR/ar53_taxonomy_$REL.tsv.gz"
java -Xmx16g -cp "$CP" tax.TaxTree taxdump/names.dmp taxdump/nodes.dmp \
    taxdump/merged.dmp "gtdb_$REL.taxtree.gz" > taxtree.log 2>&1

# Genomes
if [ ! -f genomes/.unpacked ]; then
    echo "[$(date +%T)] unpacking $TARBALL"
    tar -I pigz -xf "$TARBALL" -C genomes
    touch genomes/.unpacked
fi
find "$PWD/genomes" -name '*.fna.gz' | sort > genome_list.txt
echo "[$(date +%T)] $(wc -l < genome_list.txt) genomes"

# Sketches. The genome list order is fixed, so batch N holds the same genomes
# on every run and a finished batch is never redone.
rm -f batches/list_*
split -l "$BATCH" -d -a 4 genome_list.txt batches/list_
sketch_one() {  # $1 = batch list file
    local list=$1 out=${1/list_/batch_}.sketch
    [ -s "$out" ] && return 0
    java -da -Xmx12g -cp "$CP" sketch.SketchMaker in="$(paste -sd, "$list")" \
        out="$out.tmp" k=32,24 autosize=t mode=perfile fixjunk=t threads=5 ow=t \
        2>"$out.log" && mv "$out.tmp" "$out"
}
export -f sketch_one; export CP
ls batches/list_* | xargs -P "$JOBS" -I{} bash -c \
    'sketch_one {} && echo "[$(date +%T)] done $(basename {})"'

cat batches/batch_*.sketch > raw.sketch
python3 "$HERE/set_sketch_taxids.py" taxdump/accession2taxid.tsv raw.sketch \
    "gtdb_${REL}_k2432.sketch"
n=$(grep -c '^#SZ:' "gtdb_${REL}_k2432.sketch")
want=$(wc -l < genome_list.txt)
echo "[$(date +%T)] $n sketches in $WORK/gtdb_${REL}_k2432.sketch ($want genomes)"
[ "$n" -eq "$want" ] || { echo "sketch count differs from genome count" >&2; exit 1; }
rm raw.sketch
