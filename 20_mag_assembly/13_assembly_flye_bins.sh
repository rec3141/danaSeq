#!/bin/bash
# runs assembly on mapped reads
# outdir is quick-nano.sh output directory containing minimap_bins/*.fq

source activate flye29

outdir=$1
bin=$2
mkdir $outdir/flye_bins

# concatenate all fastqs
cat $outdir/minimap_bins/${bin}*.fq > $outdir/minimap_bins/${bin}.all.fq
~/apps/bbmap/reformat.sh in=$outdir/minimap_bins/${bin}.all.fq out=$outdir/minimap_bins/${bin}.clean.fq minlength=1000 ow tossbrokenreads uniquenames
~/apps/bbmap/clumpify.sh in=$outdir/minimap_bins/${bin}.clean.fq out=$outdir/minimap_bins/${bin}.uniq.fq dedupe 
#~/apps/bbmap/bbduk.sh in=$outdir/minimap_bins/${bin}.uniq.fq out=$outdir/minimap_bins/${bin}.noribo.fq k=15 ref=adapters ref=/home/grid/apps/bbmap/resources/16S_consensus_sequence.fa ref=/home/grid/apps/bbmap/resources/18S_consensus_sequence.fa ref=/home/grid/apps/bbmap/resources/23S_consensus_sequence.fa kmask
#ulimit -s 655360

echo "running flye on $binfile"
flye --nano-raw $outdir/minimap_bins/${bin}.uniq.fq --out-dir $outdir/flye_bins/${bin}/ --threads 8 --meta --min-overlap 1000 ##>> $outdir/log.txt 2>&1
#flye --nano-raw $outdir/all.fastq.gz --out-dir $outdir/flye --threads 8 --meta --min-overlap 5000 >> $outdir/log.txt 2>&1
#flye --nano-raw $outdir/fq/*.fastq.gz --out-dir $outdir/flye --threads 8 --meta --min-overlap 1000 >> $outdir/log.txt 2>&1

rm $outdir/minimap_bins/${bin}.*.fq

if [ -s $outdir/flye_bins/${bin}/assembly.fasta ]; then
 # sketch assembly
 ~/apps/bbmap/sendsketch.sh in=$outdir/flye_bins/${bin}/assembly.fasta address=nt

 # make assembly graph image
 ~/apps/Bandage image $outdir/flye_bins/${bin}/assembly_graph.gfa $outdir/flye_bins/${bin}/assembly_graph.png --depth --lengths --colour depth --fontsize 16

 # view assembly graph image
 eog $outdir/flye_bins/${bin}/assembly_graph.png &

fi;


