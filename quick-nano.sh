#!/bin/bash
# this script does quick qc and annotation of nanopore reads for real-time analysis
export JAVA_TOOL_OPTIONS="-Dhttps.protocols=TLSv1.2"

dir=$1
outdir=$2

if [ -z $1 ]; then echo "no input directory specified"; exit; fi
if [ -z $2 ]; then outdir=$(mktemp -u out_$(basename ${dir})_XXXXXX); fi

echo "making output directory $outdir"
mkdir -p $outdir/fa $outdir/fq $outdir/sketch $outdir/prokka $outdir/tetra

# refresh database of links in /data/fastq_pass
rm /data/fastq_pass/*.fastq.gz
#find /data/ -name "*.fastq.gz" -type f -regex ".*_pass_.*" -exec ln -s {} /data/fastq_pass/ \;
filelist=$(find /data/ -name "*.fastq.gz" -type f -regex ".*_pass_.*")

#check files for errors
#for file in /data/fastq_pass/*.fastq.gz; do
for file in $filelist; do
 echo $file $(basename $file)
 if [ -s "/data/fastq_pass/$(basename $file)" ]; then 
  echo "already done $file";
  continue;
 else
 # Test the gzip file
 gzip -t $file

 # Check the exit status
 if [ $? -eq 0 ]; then
  echo "good $file" >> $outdir/log.txt
  ln -s $file /data/fastq_pass/
 else
  tmpfile=$(mktemp --suffix=.fastq.gz)
  ~/apps/bbmap/reformat.sh in=$file out=$tmpfile ow >> $outdir/log.txt 2>&1
  gzip -t $tmpfile
  if [ $? -eq 0 ]; then
   mv $tmpfile /data/fastq_pass/
   echo "bad file fixed: $file" >> $outdir/log.txt
  else
   rm $tmpfile
#   rm $file
   echo "bad file could not be fixed: $file" >> $outdir/log.txt
  fi
 fi
fi
done;


# start loop for each file off the nanopore in the requested directory
filelist=$(find $dir -follow -type f -name '*fastq.gz' -size +1k -regex '.*fastq_pass.*' | grep -v "/data/work" | sort)

for file in $filelist; do

 base=$(basename $file .fastq.gz)
 fafile=$outdir/fa/$base.fa
 fqfile=$outdir/fq/$base.fastq.gz

 # run bbduk on fastq file, trim and remove short reads
 if [ ! -s $fafile ]; then 
  echo "running bbduk on $base"
  ~/apps/bbmap/bbduk.sh in=$file out=$fafile ref=adapters,artifacts,phix,lambda qtrim=rl trimq=15 entropy=0.9 qin=33 minlength=1000>> $outdir/log.txt 2>&1
 else 
  echo "already completed bbduk on $base"
 fi;

 # soft link fastq file into fq folder
 if [ -s $fafile ]; then 
  ln -sf $file $fqfile;
 else
  continue
 fi

 # run sendsketch on fasta file
 if [ ! -s $outdir/sketch/${base}.txt ]; then 
  echo "running sketch on $base"
  ~/apps/bbmap/sendsketch.sh in=$fafile address=nt out=$outdir/sketch/${base}.txt format=3 >> $outdir/log.txt 2>&1
 else
  echo "already completed sketch on $base"
 fi

 # run prokka on fasta file
 if [ ! -s $outdir/prokka/$base/PROKKA_*.tsv ]; then
  echo "running prokka on $base"
  prokka --metagenome --fast --cpus 0 --evalue 1e-20 --outdir $outdir/prokka/$base --quiet $fafile >> $outdir/log.txt 2>&1
 else
  echo "already completed prokka on $base"
 fi

 # run tetra on fasta file
 if [ ! -s $outdir/tetra/annotation.$base.txt ]; then
  echo "running tetra for $base";

  #make fake annotation file
  grep '>' $fafile | sed 's/>//' | paste - - -  > $outdir/tetra/annotation.$base.txt
  perl ~/apps/tetramer_freqs_esom.pl -f $fafile -a $outdir/tetra/annotation.$base.txt -min 1000 -max 10000000 >> $outdir/log.txt 2>&1

  mv Tetra_${base}* $outdir/tetra/

 else
  echo "already completed tetra on $base"
 fi

done

# run mapping etc.

cat $outdir/tetra/Tetra_*.lrn > $outdir/Tetra_all.lrn
head -n4 $outdir/Tetra_all.lrn | tail -n1 | cut -f2- -d' ' > $outdir/tnfs.txt
cat $outdir/tetra/Tetra_*.names > $outdir/Tetra_all.names

Rscript --vanilla plotting.R $outdir

cp $outdir/tsne.pdf /tara/

#./run-flye.sh $outdir

./run-mapping.sh $outdir

./run-metabat.sh $outdir

./run-mapping-bins.sh $outdir

