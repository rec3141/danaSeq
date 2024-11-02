

#bacteria
~/apps/prokka/bin/prokka --metagenome --fast --cpus 0 --evalue 1e-20 contigs.fa

# better to run on each fastq.gz that comes off the instrument
#paste <(grep product *.gff | cut -f1) <(grep product *.gff | grep -o 'ID=[^;]*' | sed 's/ID=//') <(grep -o 'product=[^;]*' *.gff | sed 's/product=//') | grep -v hypot | sort -k3,3 | sed 's/%2C/,/g'


#eukaryotes
~/apps/metaeuk/metaeuk

