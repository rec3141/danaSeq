#!/usr/bin/awk -f
# Reshape `coverm contig --methods length mean variance` (dense, one
# "<sample> Length/Mean/Variance" column triple per sample) into the MetaBAT
# depth table that CoverM's own metabat method prints:
#   contigName  contigLen  totalAvgDepth  <sample>.bam  <sample>.bam-var ...
# totalAvgDepth is the mean depth across samples, as CoverM computes it.
# Values are rounded to 4 decimals without trailing zeros, also as CoverM does.
BEGIN { FS = OFS = "\t" }
function r4(x,   v) {
    v = sprintf("%.4f", x)
    if (v ~ /\./) { sub(/0+$/, "", v); sub(/\.$/, "", v) }
    if (v == "-0") v = "0"
    return v
}
NR == 1 {
    n = (NF - 1) / 3
    if (n != int(n) || n < 1 || $2 !~ / Length$/ || $3 !~ / Mean$/ || $4 !~ / Variance$/) {
        print "jgi_depth_table.awk: unexpected CoverM header: " substr($0, 1, 200) > "/dev/stderr"
        exit 2
    }
    line = "contigName" OFS "contigLen" OFS "totalAvgDepth"
    for (i = 0; i < n; i++) {
        s = $(2 + 3 * i); sub(/ Length$/, "", s)
        line = line OFS s ".bam" OFS s ".bam-var"
    }
    print line
    next
}
$1 == "unmapped" { next }
{
    total = 0; vals = ""
    for (i = 0; i < n; i++) {
        total += $(3 + 3 * i)
        vals = vals OFS r4($(3 + 3 * i)) OFS r4($(4 + 3 * i))
    }
    print $1, $2, r4(total / n) vals
}
