#!/bin/bash

## BLAST a FASTA file and return all of the hits as a FASTA file

# Define variables
DB=all_clean_db
QUERY=$1
RESULTS=$(basename $QUERY .fasta).hits.fasta
OUTPUT=$(basename $QUERY .fasta).match.fasta

# run the blast search
#blastn -query $QUERY -db $DB -outfmt "6 qseqid sseqid sstart send" -out $RESULTS -evalue 1e-20
tblastn -query $QUERY -db $DB -outfmt "6 qseqid sseqid sstart send" -out $RESULTS -evalue 1e-20

# Initialize output file
> $OUTPUT

# Read the results file line by line
while read -r line; do
  # Extract query ID, subject ID, start and end positions
  qseqid=$(echo "$line" | awk '{print $1}')
  sseqid=$(echo "$line" | awk '{print $2}')
  sstart=$(echo "$line" | awk '{print $3}')
  send=$(echo "$line" | awk '{print $4}')

  # Determine the correct start and end (to handle reverse alignments)
  if [ "$sstart" -lt "$send" ]; then
    start=$sstart
    end=$send
  else
    start=$send
    end=$sstart
  fi

  # Fetch the sequence from the database
  blastdbcmd -db $DB -entry $sseqid -range $start-$end -outfmt "%f" >> $OUTPUT

done < $RESULTS

# Print completion message
echo "Matching sequences saved to $OUTPUT"
