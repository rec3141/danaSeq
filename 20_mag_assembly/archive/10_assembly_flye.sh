#!/bin/bash
OUTDIR=$1
PRJ="CMO2025"
#PRJDIR=$2 #/data/project_CMO2
#PRJDAT=$3 #CMO2

CONCATDIR="/data/project_"${PRJ}"/concat"

#need to edit this to follow symlinks
for d1 in `ls -d /data/minknow/${PRJ}/20*`; do 
 for d2 in `ls -d $d1/*/`; do 
  for d3 in `ls -d $d2/fastq_pass/barcode*/`; do 
   echo $d2 $d3; 
   cat $d3/*.fastq.gz | /work/apps/bbmap/dedupe.sh in=stdin.fastq.gz out=${CONCATDIR}/$(basename ${d2})_$(basename ${d3}).fastq.gz outd=${CONCATDIR}/dupe_$(basename ${d2})_$(basename ${d3}).fastq.gz ow
  done; 
 done; 
done;

cat ${CONCATDIR}/20*.fastq.gz > ${CONCATDIR}/all.fq.gz

date
filtlong -t 40000000000 ${CONCATDIR}/all.fq.gz | pigz > ${CONCATDIR}/all.filt.fq.gz

# conda activate flye 
#conda run -n flye --live-stream flye --out-dir tmp_all_flye1000_`date  +'%Y%m%d-%H%M%S'` --threads 12 --meta --min-overlap 1000 --nano-raw all.fq.gz
date
conda run -n flye --live-stream flye --out-dir ${OUTDIR}/flye1000 --threads 24 --meta --min-overlap 1000 --nano-hq ${CONCATDIR}/all.filt.fq.gz


