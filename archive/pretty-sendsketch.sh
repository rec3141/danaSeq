sendsketch.sh gt100k.fa address=protein format=3 mode=sequence rec
ords=1 | sort -t$'\t' -k2,2 -k4,4rn | column -t -s$'\t'
