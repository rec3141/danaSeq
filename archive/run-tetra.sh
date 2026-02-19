#!/bin/bash
# this script calculates tetranucleotide frequencies for binning
# outdir is the directory containing output from quick-nano.sh
outdir=$1
indir=$2

mkdir $outdir/tetra

for file in $outdir/$indir/*.fa; do 
 base=$(basename $file .fa)

 #make fake annotation file
 grep '>' $file | sed 's/>//' | paste - - -  > $outdir/tetra/annotation.$base.txt
 perl ~/apps/tetramer_freqs_esom.pl -f $file -a $outdir/tetra/annotation.$base.txt -min 1000 -max 10000000

 mv Tetra_$file* $outdir/tetra/

# for mapping the splits
# split=$outdir/tetra/${base}_1000_10000000_split.fasta
# ~/apps/minimap2/minimap2 -ax map-ont ref.mmi $split > /tmp/tmp.sam;
# samtools sort /tmp/tmp.sam -o $outdir/tetra/$(basename $split .fasta).bam;
# samtools index $outdir/tetra/$(basename $split .fasta).bam;
# rm /tmp/tmp.sam;
 echo "complete tetra for $base";
done

# tetraSOM
#grep '>' $assembly | sed 's/>//' | paste - - -  > annotation.txt
#perl ~/apps/tetramer_freqs_esom.pl -f $assembly -a annotation.txt -min 1000
