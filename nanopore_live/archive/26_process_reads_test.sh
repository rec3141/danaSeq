#!/bin/bash
# this script does quick qc and annotation of nanopore reads for real-time analysis
# run it from the project folder you want to analyze
# you can link folders of nanopore runs or barcodes into the project folder
export JAVA_TOOL_OPTIONS="-Dhttps.protocols=TLSv1.2"

dir=$1
outdir=$2

runsketch=$3
runkraken=$4
runtetra=$5
runprokka=$6
prokkathreads=0

dir="/data/project_CMO2"
outdir="out_dana"

if [ -z $dir ]; then echo "no input directory specified"; exit; fi
if [ -z $outdir ]; then outdir=$(mktemp -u out_$(basename ${dir})_XXXXXX); fi

# refresh database of links in /data/fastq_pass
#rm /data/fastq_pass/*.fastq.gz
#find /data/project_CMO -name "*.fastq.gz" -type f -regex ".*_pass_.*" -exec ln -s {} /data/fastq_pass/ \;
#filelist=$(find /data/project_CMO -name "*.fastq.gz" -type f -regex ".*_pass_.*")


apps="/work/apps"
bbmap="/work/apps/bbmap"
danadir="/work/apps/dana"
prokka="/work/apps/prokka/bin/prokka"
filtlong="/work/apps/Filtlong/bin/filtlong"
kraken="/usr/bin/kraken2"
KRAKEN_DB=/data/scratch/refdbs/krakendb/pluspfp_08gb #runs in 3 seconds, 30% classification

filelist=$(find $dir -follow -type f -name '*fastq.gz' -size +10k -regex '.*fastq_pass.*' | grep -v "/data/minknow" | grep barcode | grep -v '.tmp.' | sort)
if [ -z "$filelist" ]; then echo "no barcodes found in $dir"; exit 0; fi
barcodelist=$(echo $filelist | xargs -n1 basename | cut -f3 -d'_' | sort -u)
fclist=$(echo $filelist | xargs -n1 basename | cut -f1 -d'_' | sort -u)

numfiles=$(echo $filelist | wc -w)
echo "found $numfiles FASTQ files to process"

echo "making output directory $outdir"
for fc in $fclist; do
for barcode in $barcodelist; do
mkdir -p $outdir/$fc/$barcode/fa $outdir/$fc/$barcode/fq $outdir/$fc/$barcode/sketch $outdir/$fc/$barcode/prokka $outdir/$fc/$barcode/tetra $outdir/$fc/$barcode/stats $outdir/$fc/$barcode/kraken
done
done


#check files for errors
#for file in /data/fastq_pass/*.fastq.gz; do
echo "checking fastq files"
mkdir -p /data/.fastq_pass


i=0
for file in $filelist; do
    if [ -s "/data/.fastq_pass/$(basename "$file")" ]; then 
        ((i++))  # Correct increment syntax
        if (( i % 100 == 0 )); then  # Fix comparison operator
            echo -n '.'
	fi
    continue;
 else
 echo -n "+"
 # Test the gzip file
 gzip -t $file

 # Check the exit status
 if [ $? -eq 0 ]; then
#  echo "good $file" >> $outdir/log.txt
  ln -sf $file /data/.fastq_pass/
 else
  tmpfile=$(mktemp --suffix=.fastq.gz)
  $bbmap/reformat.sh in=$file out=$tmpfile ow >> $outdir/log.txt 2>&1
  gzip -t $tmpfile
  if [ $? -eq 0 ]; then
   mv $tmpfile /data/.fastq_pass/
   echo "bad file fixed: $file" >> $outdir/log.txt
  else
   rm -f $tmpfile
#   rm $file
   echo "bad file could not be fixed: $file" >> $outdir/log.txt
  fi
 fi
fi
done;
echo -n "done checking fastq files: $i / $numfiles were good"

#for fc in $fclist; do
#for barcode in $barcodelist; do

#bcfiles=$(find $dir -follow -type f -name '*fastq.gz' -size +10k -regex '.*fastq_pass.*' | grep -v "/data/minknow" | sort -R | grep barcode)
#sort by md5 so they are random but consistent across runs
bcfiles=$(find $dir -follow -type f -name '*fastq.gz' -size +10k -regex '.*fastq_pass.*' | grep -v "/data/minknow" | grep barcode | grep -v '.tmp.' | grep -v 'out_dana' | xargs -P 8 -I {} md5sum {} | sort | cut -f3- -d' ')

start_time=$(date +%s)
max_duration=3600  # 1 hour in seconds

for file in $bcfiles; do

    current_time=$(date +%s)
    elapsed=$((current_time - start_time))
    
    if [ $elapsed -gt $max_duration ]; then
        echo "Loop timeout after $elapsed seconds"
        break
    fi
    

base=""
barcode=""
fc=""
bcdir=""
fafile=""
fqfile=""
ftfile=""
 base=$(basename $file .fastq.gz)
 barcode=$(basename $file | cut -f3 -d'_')
 fc=$(basename $file | cut -f1 -d'_')
 bcdir=$outdir/$fc/$barcode
 # echo $base $fc $barcode $bcdir
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
  $bbmap/bbduk.sh in=$file out=$fqfile ref=adapters,artifacts,phix,lambda qtrim=rl trimq=15 entropy=0.75 qin=33 minlength=1500 >> $bcdir/log.txt 2>&1
  # run filtlong to remove low quality reads
  $filtlong --min_length 1500 --keep_percent 80 $fqfile > $ftfile 2>> $bcdir/log.txt
  # reformat to fasta
  $bbmap/reformat.sh in=$ftfile out=$fafile.tmp.fa fastawrap=0 >> $bcdir/log.txt 2>&1
  # remove excess headers
  cut -f1 -d' ' $fafile.tmp.fa > $fafile
  rm $ftfile $fafile.tmp.fa
# else   echo "already completed bbduk on $base"
 fi;

 # soft link fastq file into fq folder
# if [ -s $fafile ]; then 
#  ln -sf $file $fqfile;
# else
#  continue
# fi

if [ ! -s "$fafile" ]; then
#    echo "File is empty, skipping analysis: $fafile"
    continue
fi


 # run sendsketch on fasta file
 if (( $runsketch )); then
 if [ ! -s $bcdir/sketch/${base}.txt ]; then 
  echo "running sketch on $base"
  $bbmap/sendsketch.sh in=$fafile address=nt out=$bcdir/sketch/${base}.txt format=3 >> $bcdir/log.txt 2>&1
  Rscript $danadir/sketch-db.r $bcdir >> $bcdir/log.txt 2>&1
# else  echo "already completed sketch on $base"
 fi
#else  echo "skipping sendsketch"
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
if (( $runkraken )); then
 if [ ! -e $bcdir/kraken/${base}.tsv ]; then
#  DB=/data/scratch/refdbs/krakendb/standard #runs in 60 seconds, 60% classification

  echo "running kraken2 on $base using $KRAKEN_DB"
#  eval "$(micromamba shell hook -s bash)"
#  micromamba activate
  $kraken --db "$KRAKEN_DB" --use-names --threads 1 --report $bcdir/kraken/${base}.report  $fafile 2>> $bcdir/log.txt | gawk -f $danadir/kraken_parse.awk > $bcdir/kraken/${base}.tsv
#  micromamba deactivate
 Rscript $danadir/kraken-db.r $bcdir >> $bcdir/log.txt 2>&1
 Rscript $danadir/krakenreport-db.r $bcdir >> $bcdir/log.txt 2>&1
# else  echo "already completed kraken2 on $base"
 fi
#else  echo "skipping kraken"
fi

# run prokka on fasta file
if (( $runprokka )); then
shopt -s nullglob # Expand to empty if no match
prokfiles=(`pwd`/$bcdir/prokka/$base/PROKKA_*.tsv)

if [ "${#prokfiles[@]}" -eq 0 ]; then
  echo "running prokka on $base"
  prokdir=`pwd`/$bcdir/prokka/$base
  $prokka --metagenome --fast --cpus $prokkathreads --evalue 1e-20 --outdir $prokdir --force --quiet `pwd`/$fafile >> $bcdir/log.txt 2>&1
  rm $prokdir/*.err $prokdir/*.fna $prokdir/*.fsa $prokdir/*.gbk $prokdir/*.log $prokdir/*.sqn $prokdir/*.txt
#  rm .err .fna .fsa .gbk .log .sqn .txt
# keep .faa .ffn .gff .tbl .tsv
Rscript $danadir/prokka-db.r $bcdir >> $bcdir/log.txt 2>&1

#else  echo "already completed prokka on $base"
fi
#else  echo "skipping prokka"
fi


#  # run whokaryote on fasta file
# prokfiles=(`pwd`/$bcdir/prokka/$base/PROKKA_*.tsv)
# 
# if [ "${#prokfiles[@]}" -eq 0 ]; then
#   echo "running prokka on $base"
#   prokdir=`pwd`/$bcdir/prokka/$base
# /home/grid/apps/whokaryote-1.1.2/bin/whokaryote.py --contigs barcode23/fa/FAZ28569_pass_barcode23_9d7271ef_2e2b5cfb_2.fa --out barcode23/whokaryote_FAZ28569_pass_barcode23_9d7271ef_2e2b5cfb_20 --minsize 0 --threads 4
#   ~/apps/prokka/bin/prokka --metagenome --fast --cpus 0 --evalue 1e-20 --outdir $prokdir --force --quiet `pwd`/$fafile >> $bcdir/log.txt 2>&1
#   rm $prokdir/*.err $prokdir/*.fna $prokdir/*.fsa $prokdir/*.gbk $prokdir/*.log $prokdir/*.sqn $prokdir/*.txt
# #  rm .err .fna .fsa .gbk .log .sqn .txt
# # keep .faa .ffn .gff .tbl .tsv
# Rscript /data/work/dana/prokka-db.r $bcdir >> $bcdir/log.txt 2>&1
# 
# else
#   echo "already completed prokka on $base"
# fi

if (( $runtetra )); then
 # run tetra on fasta file
 if [ ! -s $bcdir/tetra/${base}.lrn ]; then
  echo "running tetra for $base";

  #make fake annotation file
  grep '>' $fafile | sed 's/>//' | paste - - -  > $bcdir/tetra/annotation.$base.txt
  perl $apps/tetramer_freqs_esom.pl -f $fafile -a $bcdir/tetra/annotation.$base.txt -min 1500 -max 10000000 >> $bcdir/log.txt 2>&1

  if [ ! -s $bcdir/tnfs.txt ]; then
  paste <(echo "seqid") <(head -n4 Tetra_*.lrn | tail -n1 | cut -f2-) > $bcdir/tnfs.txt
  fi

  paste <(awk '$1 !~ /^%/' Tetra_${base}*.names) <(awk '$1 !~ /^%/' Tetra_${base}*.lrn) | cut -f3,5- > $bcdir/tetra/${base}.lrn
  rm Tetra_${base}.*
  rm $bcdir/tetra/annotation.$base.txt

  Rscript $danadir/tetra-db.r $bcdir >> $bcdir/log.txt 2>&1

# else  echo "already completed tetra on $base"
 fi
#else  echo "skipping tetra"
fi


done

#done

#done


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

# exec "$0" "$@"

#/work/apps/dana/quick-nano-barcode.test.sh
