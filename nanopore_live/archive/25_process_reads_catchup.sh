#!/bin/bash
# this script does quick qc and annotation of nanopore reads for real-time analysis
# run it from the project folder you want to analyze
# you can link folders of nanopore runs or barcodes into the project folder

# this file catches up the database without running any new analyses

dir=$1
outdir=$2

dir="/data/project_QEI2025"
outdir="out_dana"

if [ -z $dir ]; then echo "no input directory specified"; exit; fi
if [ -z $outdir ]; then outdir=$(mktemp -u out_$(basename ${dir})_XXXXXX); fi

# refresh database of links in /data/fastq_pass
#rm /data/fastq_pass/*.fastq.gz
#find /data/project_CMO -name "*.fastq.gz" -type f -regex ".*_pass_.*" -exec ln -s {} /data/fastq_pass/ \;
#filelist=$(find /data/project_CMO -name "*.fastq.gz" -type f -regex ".*_pass_.*")

runsketch=1
runkraken=1
runprokka=1
runtetra=1

apps="/work/apps"
danadir="/work/apps/dana"

filelist=$(find $dir -follow -type f -name '*fastq.gz' -size +10k -regex '.*fastq_pass.*'|grep -v "/data/minknow"|sort|grep barcode)
if [ -z "$filelist" ]; then echo "no barcodes found in $dir"; exit 0; fi
barcodelist=$(echo $filelist|xargs -n1 basename|cut -f3 -d'_'|sort -u)
fclist=$(echo $filelist|xargs -n1 basename|cut -f1 -d'_'|sort -u)

# start loop for each file off the nanopore in the requested directory

for fc in $fclist; do
for barcode in $barcodelist; do

bcdir=$outdir/$fc/$barcode

 # run sendsketch on fasta file
 if (( $runsketch )); then
  echo "$fc|$barcode running sendsketch update"
  Rscript $danadir/sketch-db.r $bcdir >> $bcdir/log.txt 2>&1
else echo "$fc|$barcode skipping sendsketch"
fi

 # run kraken on fasta file
if (( $runkraken )); then
 echo "$fc|$barcode running kraken db update"
 Rscript $danadir/kraken-db.r $bcdir >> $bcdir/log.txt 2>&1
 echo "$fc|$barcode running kraken report update"
 Rscript $danadir/krakenreport-db.r $bcdir >> $bcdir/log.txt 2>&1
else echo "$fc|$barcode skipping kraken"
fi

# run prokka on fasta file
if (( $runprokka )); then
echo "$fc|$barcode running prokka update"
Rscript $danadir/prokka-db.r $bcdir >> $bcdir/log.txt 2>&1
else echo "$fc|$barcode skipping prokka"
fi


if (( $runtetra )); then
  echo "$fc|$barcode running tetra update"
  Rscript $danadir/tetra-db.r $bcdir >> $bcdir/log.txt 2>&1
else echo "$fc|$barcode skipping tetra"
fi


done

done

exit 0


