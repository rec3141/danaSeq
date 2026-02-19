#!/bin/bash

#date
for f1 in `ls 2025*.fastq.gz`; do
f2=$(basename $f1 .fastq.gz)
for d1 in `ls -d /data/minknow/QEI2025/$f2`; do for d2 in `ls -d $d1/*/`; do for d3 in `ls -d $d2/fastq_pass/barcode*/`; do echo $d2 $d3; cat $d3/*.fastq.gz > $(basename ${d2})_$(basename ${d3}).fastq.gz; done; done; done;

filtlong -p 80 $f1 | pigz > $f2.filt.fq.gz

date
conda run -n flye --live-stream flye --out-dir tmp_${f2}_flye1000_`date  +'%Y%m%d-%H%M%S'` --threads 24 --meta --min-overlap 1000 --nano-hq $f2.filt.fq.gz

done;

