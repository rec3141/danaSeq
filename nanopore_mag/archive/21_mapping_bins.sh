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
mapdir=$outdir/minimap_bins
fadir=$outdir/metabat

mkdir -p $mapdir

# for each bin
binlist=$(find $fadir -type f -name "bin.*.fa" | grep -E '[0-9]' | sort -u)

for fasout in $binlist; do
bin=$(basename $fasout .fa)

#make minimap index
echo "indexing $bin"
~/apps/minimap2/minimap2 -ax map-ont -d $mapdir/${bin}_index.idx $fasout >> $outdir/log.txt 2>&1


# without barcode: FAR94440_pass_1ac55b0d_c4a4702b_9.fastq.gz
# with barcode:    FAZ32473_pass_barcode02_54dec70f_7ecc2808_135.fastq.gz

# for each flowcell
fclist=$(find $fqdir -name "*.fastq.gz" -printf '%f\n' | cut -f1 -d"_" | sort -u)
for cell in $fclist; do
echo "processing flowcell $cell"

 # for each barcode
 bclist=$(find $fqdir -name "${cell}*.fastq.gz" -printf '%f\n' | cut -f3 -d"_" | sort -u)
# bclist='barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08'
# bclist='barcode06'
 for bc in $bclist; do
  echo "mapping $fasout $cell $bc"

  # concatenate fastq.gz files
  sqcat=/tmp/${cell}_${bc}.fastq.gz 
  sqnorib=/tmp/${cell}_${bc}.noribo.fastq.gz 
  alnout=/tmp/${cell}_${bc}.aln.sam

  if [ ! -f $sqcat ]; then
   cat $fqdir/${cell}*${bc}*.fastq.gz > $sqcat
   ~/apps/bbmap/bbduk.sh in=$sqcat out=$sqnorib k=15 ref=adapters ref=/home/grid/apps/bbmap/resources/16S_consensus_sequence.fa ref=/home/grid/apps/bbmap/resources/18S_consensus_sequence.fa ref=/home/grid/apps/bbmap/resources/23S_consensus_sequence.fa kmask
  fi


  # map reads using minimap
  ~/apps/minimap2/minimap2 -ax map-ont $mapdir/${bin}_index.idx $sqnorib > $alnout 2>>$outdir/log.txt

  alnsort=$mapdir/${bin}_${cell}_${bc}.aln.sort.bam
  samtools sort $alnout -o $alnsort
  samtools index $alnsort

  alnfq=$mapdir/${bin}_${cell}_${bc}.aln.sort.fq
  bedtools bamtofastq -i $alnsort -fq $alnfq

  rm -f $alnout
  echo mapping complete for ${cell} ${bc} ${bin}

 done;
done;

# assembly bin
  echo beginning flye assembly for $cell $bc $bin

  ./run-flye-bins.sh $outdir $bin

done;

# rm -f $sqcat $sqnorib $alnfq $alnsort

