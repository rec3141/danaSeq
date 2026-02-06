#!/bin/bash
# runs assembly on QC'd reads
# outdir is quick-nano.sh output directory containing fq/*.fastq.gz

outdir=$1

# concatenate all fastqs
cat $outdir/fq/*.fastq.gz > /data/scratch/all.fastq.gz

source activate flye29
#ulimit -s 655360

echo "running flye on $outdir"
flye --nano-raw /data/scratch/all.fastq.gz --out-dir $outdir/flye --threads 8 --meta --min-overlap 2000 >> $outdir/log.txt 2>&1
#flye --nano-raw $outdir/all.fastq.gz --out-dir $outdir/flye --threads 8 --meta --min-overlap 5000 >> $outdir/log.txt 2>&1
#flye --nano-raw $outdir/fq/*.fastq.gz --out-dir $outdir/flye --threads 8 --meta --min-overlap 1000 >> $outdir/log.txt 2>&1

if [ -s $outdir/flye/assembly.fasta ]; then
 # sketch assembly
 ~/apps/bbmap/sendsketch.sh in=$outdir/flye/assembly.fasta address=nt

 # make assembly graph image
 ~/apps/Bandage image $outdir/flye/assembly_graph.gfa $outdir/flye/assembly_graph.png --depth --lengths --colour depth --fontsize 16

 # view assembly graph image
 eog $outdir/flye/assembly_graph.png &

fi;



