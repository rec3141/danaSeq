#!/bin/bash
# sketches taxonomy for bins
# dir is directory containing bins
dir=$1
cd $dir

mkdir sketch

for file in bin.*.fa; do ~/apps/bbmap/sendsketch.sh in=$file out=sketch/nt_$(basename $file .fa) address=nt format=3 records=1 ow; done
for file in bin.*.fa; do ~/apps/bbmap/sendsketch.sh in=$file out=sketch/prot_$(basename $file .fa) address=protein sixframes=t format=3; done

cat sketch/nt_bin.* | grep '^bin' | sort > sketch/nt_bin.tax
cat sketch/prot_bin.* | grep '^bin' | sort > sketch/prot_bin.tax


