#!/bin/bash

# Usage:
# ./kraken2_to_bandage.sh segments.fasta /path/to/kraken2db > taxonomy_colors.csv

FASTA="$1"
DB="$2"

if [[ -z "$FASTA" || -z "$DB" ]]; then
  echo "Usage: $0 <fastafile> </path/to/kraken2db>" >&2
  exit 1
fi

#!/bin/bash

# Usage: ./kraken2_to_bandage.sh contigs.fasta /path/to/kraken2db > taxonomy_colors.csv

FASTA="$1"
DB="$2"

if [[ -z "$FASTA" || -z "$DB" ]]; then
  echo "Usage: $0 contigs.fasta /path/to/kraken2db" >&2
  exit 1
fi

kraken2 --db "$DB" --use-names --report ${FASTA}.report "$FASTA" \
  | gawk -f /data/dana/kraken_parse.awk \
  | gawk -f /data/dana/taxonomy_colors.awk


