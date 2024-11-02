cut -f1 -d' ' ./../out_data/all_clean.fasta > /tmp/blast
makeblastdb -in /tmp/blast -dbtype nucl -out all_clean_db -parse_seqids
tblastn -query McyA_C.fasta -db nucleotide_db -out results.out -outfmt 6 -num_threads 12

sort -k11,11g results.out | awk '!seen[$2]++ && $11 < 1e-20' > results_best.txt


makeblastdb -in ./../mapping/all_clean/Tetra_all_clean_1000_5000_split.fasta -dbtype nucl -out all_clean_split_db
tblastn -query myc_cluster.fasta -db all_clean_split_db -out myc_cluster.out -outfmt 6 -num_threads 12
sort -k11,11g myc_cluster.out | awk '!seen[$1,$2]++ && $11 < 1e-20' > myc_cluster_top.txt


#trying to maximize coverage over the target sequence
tblastn -query MycC.faa -db all_clean_db -word_size 2 -gapopen 12 -gapextend 2 -evalue 10 -max_target_seqs 500 -matrix BLOSUM45 -comp_based_stats 1 -out results.out -qcov_hsp_perc 40

awk -F$'\t' '{print ">"$1 "|" $2"|" $3 "\n" $4}' results.out > results.faa


tblastn -query guanitoxin.fasta -db ./../all_clean_db -out guanitoxin_tblastn.out -outfmt "6 qseqid sseqid evalue qcovs" -num_threads 14 -evalue 1e-20 &
sort -k3,3g guanitoxin_tblastn.out | awk '!seen[$1,$2]++ && $4 > 50' > guanitoxin.best.txt
