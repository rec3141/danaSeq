#!/bin/bash

# raw files live here:
RAW=/data/project_QEI2025
# collapsed per-barcode reads go here:
SAMP=/data/project_QEI2025/
mkdir -p "$SAMP"

# need an if mmin < 60 next
# for d1 in `ls -d /data/minknow/QEI2025/2025*`; do 
#   for d2 in `ls -d $d1/*/`; do 
#     for d3 in `ls -d $d2/fastq_pass/barcode*/`; do 
#       echo $d2 $d3; cat $d3/*.fastq.gz > $(basename ${d2})_$(basename ${d3}).fastq.gz; 
#       done; 
#     done; 
#   done;
# 
# cat 2025*.fastq.gz > all.fq.gz

# 
# # discover all barcode numbers in RAW (e.g., 01,02,03,12,23)
# BCLIST=$(ls "$RAW"/*barcode*.fastq.gz \
#   | sed -E 's/.*barcode([0-9]+)\.fastq\.gz/\1/' \
#   | sort -n | uniq)
# 
# echo "Found barcodes: $BCLIST"
# 
# # collapse each barcode across runs to BCXX.fastq.gz
# for b in $BCLIST; do
#   pad=$(printf "%02d" "$b")
#   out="$SAMP/BC${pad}.fastq.gz"
#   echo ">> Collapsing barcode $pad -> $out"
#   # concatenating .gz is valid gzip (concatenated members)
#   cat "$RAW"/*barcode${b}.fastq.gz > "$out"
# done


# concatenated contigs from your SemiBin plan (or co-assembly):
CONCAT=/data/project_QEI2025/polish_mag/medaka_consensus.fasta      # change to your path
BAMDIR=/data/project_QEI2025/semibin_bams
THREADS=24
mkdir -p "$BAMDIR"

for fq in "$SAMP"/2025*.fastq.gz; do
  base=$(basename "$fq" .fastq.gz)         # e.g., BC01
  bam="$BAMDIR/${base}.sorted.bam"
  echo ">> Mapping $base -> $bam"
  minimap2 -a -x map-ont --secondary=no -t "$THREADS" "$CONCAT" "$fq" \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index -@ "$THREADS" "$bam"
done

# for fq in "$SAMP"/2025*.fastq.gz; do
#   base=$(basename "$fq" .fastq.gz)         # e.g., BC01
#   bam="$BAMDIR/${base}.sorted.bam"
#   echo ">> Extracting from $base -> $bam"
#   gunzip -c $fq | awk 'NR%4==1' | cut -f1 -d' '| tr -d '@' > $base.txt
#   samtools view -N $base.txt -o $bam polish_mag/racon1.bam \
#     | samtools sort -@ "$THREADS" -o "$bam" -
# 
#   samtools index -@ "$THREADS" "$bam"
# done

echo "Done. BAMs in $BAMDIR:"
ls -1 "$BAMDIR"/*.sorted.bam | head
# 
# conda run -n SemiBin --live-stream SemiBin2 multi_easy_bin \
#   --sequencing-type long_read \
#   -i "$CONCAT" \
#   -b "$BAMDIR"/*.sorted.bam \
#   -o /data/project_QEI2025/semibin_multi_out
#   
#   
# conda run -n SemiBin --live-stream SemiBin2 single_easy_bin \
#   --environment marine \
#   -i polish_mag/medaka_consensus.fasta \
#   -b polish_mag/racon1.bam \
#   -o semibin_out
#   
#   
   
conda run -n SemiBin --live-stream SemiBin2 single_easy_bin \
    -i polish_mag/medaka_consensus.fasta \
    -b "$BAMDIR"/*.sorted.bam \
    -o coassembly_output
  
   for f in *.fa.gz; do stats.sh $f out=$(basename $f .fa.gz).txt; sendsketch.sh in=$f out=$(basename $f .fa.gz).tsv format=3; done
   

jgi_summarize_bam_contig_depths --outputDepth coassembly_output/depths_jgi.txt --percentIdentity 80 --minMapQual 5 --referenceFasta polish_mag/medaka_consensus.fasta  semibin_bams/*.bam  
metabat2 -i polish_mag/medaka_consensus.fasta -o metabat2bins/bin -s -0 --saveCls -a coassembly_output/depths_jgi.txt
rm metabat2bins/binlist.tsv
for file in metabat2bins/bin*.fa; do  awk -v v=$(basename "$file") -v OFS='\t' '{print $0, v}' <(grep '>' $file | tr -d '>' | cut -f1) >> metabat2bins/binlist.tsv; done

mkdir -p tetra
  grep '>' polish_mag/medaka_consensus.fasta | sed 's/>//' | paste - - -  > tetra/annotation.txt
  perl /work/apps/tetramer_freqs_esom.pl -f polish_mag/medaka_consensus.fasta -a tetra/annotation.txt -min 1500 -max 10000000

  paste <(echo "seqid") <(head -n4 Tetra_*.lrn | tail -n1 | cut -f2-) > tetra/tnfs.txt

paste <(awk '$1 !~ /^%/' Tetra_medaka_consensus*.names) <(awk '$1 !~ /^%/' Tetra_medaka_consensus*.lrn) | cut -f3,5- > tetra/medaka_consensus.lrn
mv Tetra_medaka_consensus* tetra/

# need to remove header line from coassembly_output/contig_bins.tsv

conda activate SemiBin
DAS_Tool -i coassembly_output/contig_bins.tsv,metabat2bins/binlist.tsv -l semibin2,metabat2 -c polish_mag/medaka_consensus.fasta -o dastool --threads=12 --write_bin_evals --write_bins

mv dastool_*.* dastool_DASTool_bins/

  
