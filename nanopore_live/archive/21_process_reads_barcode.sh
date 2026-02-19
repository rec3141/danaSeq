#!/bin/bash
# this script does quick qc and annotation of nanopore reads for real-time analysis
# run it from the project folder you want to analyze
# you can link folders of nanopore runs or barcodes into the project folder
export JAVA_TOOL_OPTIONS="-Dhttps.protocols=TLSv1.2 -Xmx2048M"

dir=$1
outdir=$2

if [ -z $1 ]; then echo "no input directory specified"; exit; fi
if [ -z $2 ]; then outdir=$(mktemp -u out_$(basename ${dir})_XXXXXX); fi


# refresh database of links in /data/fastq_pass
#rm /data/fastq_pass/*.fastq.gz
#find /data/project_CMO -name "*.fastq.gz" -type f -regex ".*_pass_.*" -exec ln -s {} /data/fastq_pass/ \;
#filelist=$(find /data/project_CMO -name "*.fastq.gz" -type f -regex ".*_pass_.*")

filelist=$(find $dir -follow -type f -name '*fastq.gz' -size +1k -regex '.*fastq_pass.*' | grep -v "/data/work" | sort | grep barcode)
if [ -z "$filelist" ]; then echo "no barcodes found in $dir"; exit 0; fi
barcodelist=$(echo $filelist | xargs -n1 basename | cut -f3 -d'_' | sort -u)
fclist=$(echo $filelist | xargs -n1 basename | cut -f1 -d'_' | sort -u)

echo "making output directory $outdir"
for fc in $fclist; do
for barcode in $barcodelist; do
mkdir -p $outdir/$fc/$barcode/fa $outdir/$fc/$barcode/fq $outdir/$fc/$barcode/sketch $outdir/$fc/$barcode/prokka $outdir/$fc/$barcode/tetra $outdir/$fc/$barcode/stats $outdir/$fc/$barcode/kraken
done
done

#check files for errors
#for file in /data/fastq_pass/*.fastq.gz; do
echo "checking fastq files"
for file in $filelist; do
# echo $file $(basename $file)
 if [ -s "/data/fastq_pass/$(basename $file)" ]; then 
echo -n '.'
#  echo "already done $file";
  continue;
 else
 # Test the gzip file
 gzip -t $file

 # Check the exit status
 if [ $? -eq 0 ]; then
  echo "good $file" >> $outdir/log.txt
  ln -sf $file /data/fastq_pass/
 else
  tmpfile=$(mktemp --suffix=.fastq.gz)
  ~/apps/bbmap/reformat.sh -Xmx2048M in=$file out=$tmpfile ow >> $outdir/log.txt 2>&1
  gzip -t $tmpfile
  if [ $? -eq 0 ]; then
   mv $tmpfile /data/fastq_pass/
   echo "bad file fixed: $file" >> $outdir/log.txt
  else
   rm -f $tmpfile
#   rm $file
   echo "bad file could not be fixed: $file" >> $outdir/log.txt
  fi
 fi
fi
done;
echo "done checking fastq files"

# start loop for each file off the nanopore in the requested directory
filelist=$(find $dir -follow -type f -name '*fastq.gz' -size +1k -regex '.*fastq_pass.*' | grep -v "/data/work" | sort)

for fc in $fclist; do
for barcode in $barcodelist; do
bcfiles=$(find $dir -follow -type f -name '*'${barcode}'*fastq.gz' -size +1k -regex '.*fastq_pass.*' | grep -v "/data/work" | sort)

for file in $bcfiles; do
 base=$(basename $file .fastq.gz)
 bcdir=$outdir/$fc/$barcode
 fafile=$bcdir/fa/$base.fa
 fqfile=$bcdir/fq/$base.fastq.gz
 ftfile=$bcdir/fq/$base.filt.fastq
#if [ ! -s "$fafile" ]; then
#    echo "File is empty: $fafile"
#    continue
#fi

 # run bbduk on fastq file, trim and remove short reads
 if [ ! -e "$fafile" ]; then 
  echo "running bbduk on $base"
  ~/apps/bbmap/bbduk.sh -Xmx2048M in=$file out=$fqfile ref=adapters,artifacts,phix,lambda qtrim=rl trimq=15 entropy=0.75 qin=33 minlength=1500 >> $bcdir/log.txt 2>&1
  # run filtlong to remove low quality reads
  ~/apps/Filtlong/bin/filtlong --min_length 1500 --keep_percent 80 $fqfile > $ftfile 2>> $bcdir/log.txt
  # reformat to fasta
  ~/apps/bbmap/reformat.sh -Xmx2048M in=$ftfile out=$fafile.tmp.fa fastawrap=0 >> $bcdir/log.txt 2>&1
  # remove excess headers
  cut -f1 -d' ' $fafile.tmp.fa > $fafile
  rm $ftfile $fafile.tmp.fa
 else 
  echo "already completed bbduk on $base"
 fi;

 # soft link fastq file into fq folder
# if [ -s $fafile ]; then 
#  ln -sf $file $fqfile;
# else
#  continue
# fi

if [ ! -s "$fafile" ]; then
    echo "File is empty, skipping analysis: $fafile"
    continue
fi


 # run sendsketch on fasta file
 if [ ! -s $bcdir/sketch/${base}.txt ]; then 
  echo "running sketch on $base"
  ~/apps/bbmap/sendsketch.sh -Xmx2048M in=$fafile address=nt out=$bcdir/sketch/${base}.txt format=3 >> $bcdir/log.txt 2>&1
  Rscript /data/work/dana/sketch-db.r $bcdir >> $bcdir/log.txt 2>&1
 else
  echo "already completed sketch on $base"
 fi

 # run stats on fasta file
# we're getting this from GFF instead
# if [ ! -s $outdir/stats/${base}.tsv ]; then 
#  echo "running stats on $base"
#  ~/apps/bbmap/stats.sh in=$fafile format=0 gc=$outdir/stats/${base}.tsv >> $outdir/log.txt 2>&1
#  Rscript /data/work/dana/stats-db.r $outdir
# else
#  echo "already completed stats on $base"
# fi


 # run kraken on fasta file
 if [ ! -e $bcdir/kraken/${base}.tsv ]; then
  DB=/data/scratch/refdbs/krakendb/pluspfp_08gb #runs in 3 seconds, 30% classification
#  DB=/data/scratch/refdbs/krakendb/standard #runs in 60 seconds, 60% classification

  echo "running kraken2 on $base using $DB"
  eval "$(micromamba shell hook -s bash)"
  micromamba activate
  kraken2 --db "$DB" --use-names --threads 1 --report $bcdir/kraken/${base}.report  $fafile 2>> $bcdir/log.txt | gawk -f /data/work/dana/kraken_parse.awk > $bcdir/kraken/${base}.tsv
  micromamba deactivate
 Rscript /data/work/dana/kraken-db.r $bcdir >> $bcdir/log.txt 2>&1
 Rscript /data/work/dana/krakenreport-db.r $bcdir >> $bcdir/log.txt 2>&1
 else
  echo "already completed kraken2 on $base"
 fi


 # run prokka on fasta file
shopt -s nullglob # Expand to empty if no match
prokfiles=(`pwd`/$bcdir/prokka/$base/PROKKA_*.tsv)

if [ "${#prokfiles[@]}" -eq 0 ]; then
  echo "running prokka on $base"
  prokdir=`pwd`/$bcdir/prokka/$base
  ~/apps/prokka/bin/prokka --metagenome --fast --cpus 0 --evalue 1e-20 --outdir $prokdir --force --quiet `pwd`/$fafile >> $bcdir/log.txt 2>&1
  if [ ! -s $prokdir/*.err ]; then rm $prokdir/*.err $prokdir/*.fna $prokdir/*.fsa $prokdir/*.gbk $prokdir/*.log $prokdir/*.sqn $prokdir/*.txt; fi

#  rm .err .fna .fsa .gbk .log .sqn .txt
# keep .faa .ffn .gff .tbl .tsv
Rscript /data/work/dana/prokka-db.r $bcdir >> $bcdir/log.txt 2>&1

else
  echo "already completed prokka on $base"
fi

 # run tetra on fasta file
 if [ ! -s $bcdir/tetra/${base}.lrn ]; then
  echo "running tetra for $base";

  #make fake annotation file
  grep '>' $fafile | sed 's/>//' | paste - - -  > $bcdir/tetra/annotation.$base.txt
  perl ~/apps/tetramer_freqs_esom.pl -f $fafile -a $bcdir/tetra/annotation.$base.txt -min 1500 -max 10000000 >> $bcdir/log.txt 2>&1

  if [ ! -s $bcdir/tnfs.txt ]; then
  paste <(echo "seqid") <(head -n4 Tetra_*.lrn | tail -n1 | cut -f2-) > $bcdir/tnfs.txt
  fi

  paste <(awk '$1 !~ /^%/' Tetra_${base}*.names) <(awk '$1 !~ /^%/' Tetra_${base}*.lrn) | cut -f3,5- > $bcdir/tetra/${base}.lrn
  rm Tetra_${base}.*
  rm $bcdir/tetra/annotation.$base.txt

  Rscript /data/work/dana/tetra-db.r $bcdir >> $bcdir/log.txt 2>&1

 else
  echo "already completed tetra on $base"
 fi

done

done

done


# run mapping etc.
#Rscript --vanilla /data/work/dana/edna_mapping.r $dir


# concatenate all sendsketches
#for file in $(find $outdir/sketch/ -name "*.txt"); do awk '$1 !~ /^#/' $file | cut -f1-3; done > $outdir/sendsketch_all.tsv
#Rscript /data/work/dana/sendsketch-db.r
#rm $outdir/sendsketch_all.tsv


# concatenate all TNFs
#for file in $outdir/tetra/*.lrn; do f1=${file/.lrn/.names}; paste <(awk '$1 !~ /^%/' $f1) <(awk '$1 !~ /^%/' $file) ; done | cut -f3,5- > $outdir/Tetra_all.lrn
#Rscript /data/work/dana/tetra-db.r
#rm $outdir/Tetra_all.lrn
#cat $outdir/tetra/*.lrn > $outdir/Tetra_all.lrn
#paste <(echo "seqid") <(head -n4 $outdir/Tetra_all.lrn | tail -n1 | cut -f2-) > $outdir/tnfs.txt

#Rscript /data/work/dana/prokka-db.r

#cat $outdir/tetra/Tetra_*.names > $outdir/Tetra_all.names

#allfafile=$outdir/all_clean_1500.fasta
#cat $outdir/fa/*.fa > $allfafile
#grep '>' $allfafile | sed 's/>//' | paste - - -  > $outdir/annotation.all_clean_1500.txt

#ls $outdir/fa/*.fa | parallel perl ~/apps/tetramer_freqs_esom.pl -f {} -a $outdir/annotation.all_clean_1500.txt -min 1500 -max 10000000 >> $outdir/log.txt 2>&1
#mv Tetra_all_clean_1500* $outdir/tetra/
#rm $outdir/fa/*.fa

#rm $outdir/metabat/tmp
#for file in $outdir/metabat/*.fa; do
#grep "^>" "$file" | sed 's/^>//' | awk -v filename="$file" '{print filename, $0}' >> $outdir/metabat/tmp; 
#done
#mv $outdir/metabat/tmp $outdir/metabat/bin.contigs.txt


exit 0

#./run-mapping.sh $outdir

#./run-metabat.sh $outdir

#./run-sketch-bins.sh $outdir

#Rscript --vanilla plotting.R $outdir

#./run-flye.sh $outdir

#./run-mapping-bins.sh $outdir

#./run-flye-bins.sh $outdir


