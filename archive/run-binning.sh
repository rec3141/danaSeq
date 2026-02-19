#!/bin/bash
source ~/.bash_profile

assembly=$1
#assembly=~/Desktop/out_data/flye_assembly/assembly.fasta

# jgi_summarize doesn't work well with long reads?
# but running it through runMetabat.sh seems to work so may just be options

#~/apps/metabat/bin/jgi_summarize_bam_contig_depths --outputDepth depths.txt --noIntraDepthVariance --referenceFasta $assembly *.bam
PCTID=70; ~/apps/metabat/runMetaBat.sh -m 1500 -s 60000 --saveCls --unbinned $assembly *.bam

# above only results in a few bins (9)

# using just the TNF (no depth coverage) results in many more bins (42)
~/apps/metabat/bin/metabat2 -i <(gzip -c $assembly) -o bin -m 1500 -s 60000 --saveCls --unbinned -v


# tetraSOM
grep '>' $assembly | sed 's/>//' | paste - - -  > annotation.txt
perl ~/apps/tetramer_freqs_esom.pl -f $assembly -a annotation.txt -min 1000

