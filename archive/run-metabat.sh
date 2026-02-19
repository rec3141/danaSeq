#!/bin/bash
# this script runs metabat on co-assembly using mapped split reads
# outdir is the output directory from quick-nano.sh
#source ~/.bash_profile

outdir=$1
#assembly=$outdir/flye/assembly.fasta
assembly=$outdir/mapping.assembly.fasta
fadir=$outdir/fa

mkdir $outdir/metabat

PCTID=70

#run metabat
#~/apps/metabat/runMetaBat.sh -m 1500 -s 60000 --saveCls --unbinned $coassembly *.bam
# above only results in a few bins (9)

# using just the TNF (no depth coverage) results in many more bins (42)
echo "running metabat"
~/apps/metabat/bin/metabat2 -i <(gzip -c $assembly) -o $outdir/metabat/bin -m 1500 -s 60000 --saveCls --unbinned -v >> $outdir/log.txt 2>&1
echo "metabat complete"

