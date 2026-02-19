#!/bin/bash

mkdir tetra

for bin in bin*.fa; do 
grep '>' $bin | sed 's/>//' | paste - - -  > tetra/annotation.$bin.txt
perl ~/apps/tetramer_freqs_esom.pl -f $bin -a tetra/annotation.$bin.txt -min 1000 -max 1000

split=tetra/${bin}_1000_1000_split.fasta

~/apps/minimap2/minimap2 -ax map-ont ref.mmi $split > tmp.sam; 
samtools sort tmp.sam -o $(basename $split .fasta).bam; 
samtools index $(basename $split .fasta).bam; 
rm tmp.sam;
echo $bin;
done
