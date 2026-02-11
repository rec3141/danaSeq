#!/bin/bash
# this program maps long reads to the an assembly using minimap
# mapping is then used for binning with metabat2 or visualization
# outdir is directory containing long reads
source ~/.bash_profile

# edit need to combine all non-target barcodes into unclassified
# edit need to map non-barcoded reads

outdir=$1
fqdir=$outdir/fq
#assembly=$outdir/flye/assembly.fasta
#assembly=$outdir/all.coassembly.fasta
#mapdir=$outdir/minimap
mapdir=$outdir/minimap_fa
assembly=$outdir/all_clean_1500.fasta
fasout=$outdir/mapping.assembly.fasta
mkdir -p $mapdir

#cat $assembly $outdir/all_clean_1500.fasta > $fasout
ln -s $assembly $fasout

#make minimap index
echo "making index"
~/apps/minimap2/minimap2 -ax map-ont -d $mapdir/index.idx $fasout >> $outdir/log.txt 2>&1

# without barcode: FAR94440_pass_1ac55b0d_c4a4702b_9.fastq.gz
# with barcode:    FAZ32473_pass_barcode02_54dec70f_7ecc2808_135.fastq.gz

# for each flowcell
fclist=$(find $fqdir -name "*.fastq.gz" -printf '%f\n' | cut -f1 -d"_" | sort -u)
for cell in $fclist; do

 # for each barcode
 bclist=$(find $fqdir -name "*.fastq.gz" -printf '%f\n' | grep $cell | cut -f3 -d"_" | sort -u)
# bclist='barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08'
# bclist='barcode06'
 for bc in $bclist; do
  echo "mapping $cell $bc"

  # concatenate fastq.gz files
  sqcat=/tmp/${cell}_${bc}.fastq.gz 
  alnout=/tmp/${cell}_${bc}.aln.sam
  alnsort=$mapdir/${cell}_${bc}.aln.sort.bam

  if [ ! -f $sqcat ]; then
   cat $fqdir/${cell}*${bc}*.fastq.gz > $sqcat
  fi

  # map reads using minimap
  ~/apps/minimap2/minimap2 -ax map-ont $mapdir/index.idx $sqcat > $alnout 2>>$outdir/log.txt

  samtools sort $alnout -o $alnsort
  samtools index $alnsort

  rm $sqcat $alnout
  echo complete ${cell} ${bc}
 done;
done;

# run binning

