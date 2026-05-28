#!/usr/bin/awk -f

BEGIN {OFS="\t"}
# Keep SAM headers
#$0 ~ /^@/ {print; next}
$0 ~ /^@/ {next}

{
    mapq = $5
    cigar = $6
    len = 0

    # Parse the CIGAR string to calculate aligned length (count M, =, X)
    while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
        op = substr(cigar, RSTART + RLENGTH - 1, 1)
        val = substr(cigar, RSTART, RLENGTH - 1) + 0
        if (op ~ /[M=X]/) len += val
        cigar = substr(cigar, RSTART + RLENGTH)
    }

    if (mapq >= 1 && len >= 10) print
}
