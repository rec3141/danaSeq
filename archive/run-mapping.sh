#!/bin/bash
# this program maps long reads to the an assembly using minimap
# mapping is then used for binning with metabat2
# dir is directory containing long reads
# assembly is assembly; can also just use FASTA formatted read library
source ~/.bash_profile

assembly=$1
dir=$2

mkdir minimap

#assembly=~/Desktop/out_data/flye_assembly/assembly.fasta

# for each run
#for dir in July-25-2024 July-26-2024; do
# for each barcode
# for bc in barcode01 barcode02 barcode03 barcode04; do
# without barcode
# FAR94440_pass_1ac55b0d_c4a4702b_9.fastq.gz
# with barcode
# FAZ32473_pass_barcode02_54dec70f_7ecc2808_135.fastq.gz
fclist=$(find $dir -name "*.fastq.gz" -printf '%f\n' | cut -f1 -d"_" | sort -u)
# for each flowcell
for cell in fclist; do
# for each barcode
bclist=$(find $dir -name "*.fastq.gz" -printf '%f\n' | cut -f3 -d"_" | sort -u)
for bc in bclist; do

# concatenate fastq.gz files
 sqcat=/tmp/${cell}_${bc}.fastq.gz 
 sqcln=/tmp/${cell}_${bc}.clean.fastq.gz 
 alnout=minimap/${cell}_${bc}.aln.sam
 alnsort=minimap/${cell}_${bc}.aln.sort.bam

if [ ! -f $sqcat ]; then
cat $dir/${cell}*${bc}*.fastq.gz > $sqcat
#cat /data/$dir/*/*/fastq_pass/$bc/*.fastq.gz > $sqcat
fi

# filter fastq.gz

# ~/apps/bbmap/bbduk.sh in=$sqcat out=$sqcln ref=adapters,phix entropy=0.9 minlength=1000

# ~/apps/minimap2/minimap2 -ax map-ont $assembly $sqcln > $alnout
 ~/apps/minimap2/minimap2 -ax map-ont $assembly $sqcat > $alnout

 samtools sort $alnout -o $alnsort
 samtools index $alnsort

 echo complete ${cell} ${bc}
 done;
done;

#~/apps/metabat/bin/jgi_summarize_bam_contig_depths --outputDepth depths.txt --noIntraDepthVariance --referenceFasta $assembly *.bam
#~/apps/metabat/bin/metabat2 -i <(gzip -c $assembly) -o bin -m 1500 -s 60000 --saveCls --unbinned -v 
#cp ./../bin.*.fa ./; for file in bin.*.fa; do sed -i "s|>|>${file}.|g" $file; done

