#!/bin/bash

OUTDIR=$1
PRJ="CMO2025"
METAFILE="/data/project_CMO2025/CMO\ 2025\ MICB.xlsx"

# main project directory:
PRJDIR="/data/project_"${PRJ}

# write output files here
if [[ -z "${OUTDIR}" ]]; then
OUTDIR=${PRJDIR}"/map_and_bin_"`date  +'%Y%m%d-%H%M%S'`
fi

mkdir -p "$OUTDIR"

# collapsed per-barcode reads are here:
FASTQS="/data/project_"${PRJ}"/concat"
mkdir -p "${FASTQS}"

# concatenated contigs from your SemiBin plan (or co-assembly):
ASSEMBLY=${OUTDIR}"/flye1000/assembly.fasta" # change to your path

if [ -s "${ASSEMBLY}" ]; then echo "assembly completed"; 
else /data/dana/run-flye.sh ${OUTDIR} ${PRJ}
fi

if [ ! -s $ASSEMBLY ]; then echo "no assembly"; exit 1; fi

BAMDIR=$OUTDIR/bamdir
mkdir -p "$BAMDIR"

THREADS=24

for fq in "$FASTQS"/20*.fastq.gz; do
  base=$(basename "$fq" .fastq.gz)         # e.g., BC01
  bam="$BAMDIR/${base}.sorted.bam"
 if [ ! -s $bam ]; then 
  echo ">> Mapping $base -> $bam"
  minimap2 -a -x map-ont --secondary=no -t "$THREADS" "$ASSEMBLY" "$fq" \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index -@ "$THREADS" "$bam"
fi
done

echo "Done. BAMs in $BAMDIR:"
ls -l "$BAMDIR"/*.sorted.bam | head

if [ ! -s "$BAMDIR/depths_jgi.txt" ]; then
jgi_summarize_bam_contig_depths --outputDepth $BAMDIR/depths_jgi.txt --percentIdentity 80 --minMapQual 5 --referenceFasta "$ASSEMBLY"  "$BAMDIR"/*.bam
fi

if [ ! -s "$OUTDIR/semibin2/contig_bins.tsv" ]; then
conda run -n SemiBin --live-stream SemiBin2 single_easy_bin \
    -i "$ASSEMBLY" \
    -b "$BAMDIR"/*.sorted.bam \
    -o $OUTDIR/semibin2

# need to remove header line from contig_bins.tsv
awk 'NR > 1' $OUTDIR/semibin2/contig_bins.tsv > /tmp/tmpbins
mv /tmp/tmpbins $OUTDIR/semibin2/contig_bins.tsv

fi

if [ ! -s "$OUTDIR/metabat2/contig_bins.tsv" ]; then
metabat2 -i "$ASSEMBLY" -o $OUTDIR/metabat2/bin --saveCls --minClsSize 50000 -a $BAMDIR/depths_jgi.txt

rm $OUTDIR/metabat2/contig_bins.tsv
for file in $OUTDIR/metabat2/bin*.fa; do  awk -v v=$(basename "$file") -v OFS='\t' '{print $0, v}' <(grep '>' $file | tr -d '>' | cut -f1) >> $OUTDIR/metabat2/contig_bins.tsv; done
fi

if [ ! -s "$OUTDIR/dastool_DASTool_bins/dastool_DASTool_contig2bin.tsv" ]; then
conda run -n SemiBin --live-stream DAS_Tool -i $OUTDIR/semibin2/contig_bins.tsv,$OUTDIR/metabat2/contig_bins.tsv -l semibin2,metabat2 -c $ASSEMBLY -o $OUTDIR/dastool --threads=12 --write_bin_evals --write_bins
mv $OUTDIR/dastool*.* $OUTDIR/dastool_DASTool_bins/
cat $OUTDIR/dastool_DASTool_bins/*.fa > $OUTDIR/dastool_DASTool_bins/allbins.fa
fi

if [ ! -s "$OUTDIR/dastool_DASTool_bins/kaiju.allbins.summary.tsv" ]; then
echo "running kaiju"
conda run -n SemiBin --live-stream kaiju -t /data/scratch/refdbs/kaiju/progenomes/nodes.dmp -f /data/scratch/refdbs/kaiju/progenomes/kaiju_db_progenomes.fmi -i $OUTDIR/dastool_DASTool_bins/allbins.fa -z 12 -o $OUTDIR/dastool_DASTool_bins/kaiju.allbins.tsv
conda run -n SemiBin --live-stream kaiju-addTaxonNames -t /data/scratch/refdbs/kaiju/progenomes/nodes.dmp -n /data/scratch/refdbs/kaiju/progenomes/names.dmp -i $OUTDIR/dastool_DASTool_bins/kaiju.allbins.tsv -o $OUTDIR/dastool_DASTool_bins/kaiju.allbins-taxa.tsv -r superkingdom,phylum,class,order,family,genus,species
join -t $'\t' <(sort -k1,1 $OUTDIR/dastool_DASTool_bins/dastool.seqlength) <(paste <(sort -k1,1 $OUTDIR/dastool_DASTool_bins/dastool_DASTool_contig2bin.tsv) <(sort -k2,2 $OUTDIR/dastool_DASTool_bins/kaiju.allbins-taxa.tsv)) | sort -k3,3 -k2,2rn > $OUTDIR/dastool_DASTool_bins/kaiju.allbins.summary.tsv
fi

if [ ! -s "$OUTDIR/tetra/tnfs.txt" ]; then
mkdir -p $OUTDIR/tetra
paste <(grep '>' "$ASSEMBLY" | sed 's/>//') <(grep '>' "$ASSEMBLY" | sed 's/>//') <(grep '>' "$ASSEMBLY" | sed 's/>//')  > $OUTDIR/tetra/annotation.txt
perl /work/apps/tetramer_freqs_esom.pl -f $ASSEMBLY -a $OUTDIR/tetra/annotation.txt -min 1500 -max 10000000
paste <(echo "seqid") <(head -n4 Tetra_*.lrn | tail -n1 | cut -f2-) > $OUTDIR/tetra/tnfs.txt
paste <(awk '$1 !~ /^%/' Tetra_*.names) <(awk '$1 !~ /^%/' Tetra_*.lrn) | cut -f3,5- > $OUTDIR/tetra/assembly.lrn
mv Tetra_* $OUTDIR/tetra/
fi

for f in $OUTDIR/dastool_DASTool_bins/*.fa; do stats.sh $f out=$OUTDIR/dastool_DASTool_bins/$(basename $f .fa).txt; sendsketch.sh in=$f out=$OUTDIR/dastool_DASTool_bins/ss_refseq_$(basename $f .fa).tsv level=3 format=3 address=refseq; done
for f in $OUTDIR/dastool_DASTool_bins/*.fa; do stats.sh $f out=$OUTDIR/dastool_DASTool_bins/$(basename $f .fa).txt; sendsketch.sh in=$f out=$OUTDIR/dastool_DASTool_bins/ss_nt_$(basename $f .fa).tsv format=3 address=nt; done
for f in $OUTDIR/dastool_DASTool_bins/*.fa; do stats.sh $f out=$OUTDIR/dastool_DASTool_bins/$(basename $f .fa).txt; sendsketch.sh in=$f out=$OUTDIR/dastool_DASTool_bins/ss_prot_$(basename $f .fa).tsv level=3 format=3 address=protein; done


# need to download database before first run
conda run -n checkm2 --live-stream checkm2 predict --threads 24 --input $OUTDIR/dastool_DASTool_bins/ -x fa --output-directory $OUTDIR/dastool_DASTool_bins/checkm2

Rscript /data/dana/plot-bins.R $OUTDIR $METAFILE
