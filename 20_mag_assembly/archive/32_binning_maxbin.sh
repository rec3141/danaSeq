#!/bin/bash
outdir=$1

rm $outdir/maxbin_coverage.txt
rm $outdir/maxbin_reads_list.txt

for file in $outdir/minimap_fa/*.bam; do
echo $file
#bedtools bamtofastq -i $file -fq $file.fq
samtools depth -a $file >> $outdir/maxbin_coverage.txt
done

#ls $outdir/minimap_fa/*.fq > $outdir/maxbin_reads_list.txt

#~/apps/MaxBin-2.2.7/run_MaxBin.pl -contig $outdir/mapping.assembly.fasta -reads_list $outdir/maxbin_reads_list.txt -out $outdir/maxbin 
~/apps/MaxBin-2.2.7/run_MaxBin.pl -contig $outdir/mapping.assembly.fasta -abund $outdir/maxbin_coverage.txt -out $outdir/maxbin 

