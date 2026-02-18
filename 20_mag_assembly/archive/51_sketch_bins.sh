#!/bin/bash
# sketches taxonomy for bins
# ourdir is directory containing project files
outdir=$1

for file in $outdir/metabat/bin.*.fa; do ~/apps/bbmap/sendsketch.sh in=$file out=$outdir/sketch/refseq_$(basename $file .fa) address=refseq format=3 records=1 ow; done
for file in $outdir/metabat/bin.*.fa; do ~/apps/bbmap/sendsketch.sh in=$file out=$outdir/sketch/nt_$(basename $file .fa) address=nt format=3 records=1 ow; done
for file in $outdir/metabat/bin.*.fa; do ~/apps/bbmap/sendsketch.sh in=$file out=$outdir/sketch/prot_$(basename $file .fa) address=protein translate=t format=3 records=1 ow; done

rm $outdir/*.tax
cat $outdir/sketch/refseq_bin.* | grep '^bin' | sort > $outdir/sketch/refseq_bin.tax
cat $outdir/sketch/nt_bin.* | grep '^bin' | sort > $outdir/sketch/nt_bin.tax
cat $outdir/sketch/prot_bin.* | grep '^bin' | sort > $outdir/sketch/prot_bin.tax


