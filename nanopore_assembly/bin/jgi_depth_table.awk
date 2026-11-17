#!/usr/bin/awk -f
# Reshape `coverm contig --methods length mean variance` (dense, one
# "<sample> Length/Mean/Variance" column triple per sample) into the MetaBAT
# depth table that CoverM's own metabat method prints:
#   contigName  contigLen  totalAvgDepth  <sample>.bam  <sample>.bam-var ...
# totalAvgDepth is the mean depth across samples, as CoverM computes it.
# Values are rounded to 4 decimals without trailing zeros, also as CoverM does.
BEGIN { FS = OFS = "\t" }
function fail(message) { failed = 1; print "jgi_depth_table.awk: " message > "/dev/stderr"; exit 2 }
function depth(value) { return value ~ /^[+]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$/ }
function r4(x,   v) {
    v = sprintf("%.4f", x)
    if (v ~ /\./) { sub(/0+$/, "", v); sub(/\.$/, "", v) }
    if (v == "-0") v = "0"
    return v
}
NR == 1 {
    n = (NF - 1) / 3
    if ($1 != "Contig" || n != int(n) || n < 1)
        fail("unexpected CoverM header: " substr($0, 1, 200))
    line = "contigName" OFS "contigLen" OFS "totalAvgDepth"
    for (i = 0; i < n; i++) {
        s = $(2 + 3 * i)
        if (!sub(/ Length$/, "", s) || !length(s) || $(3 + 3 * i) != s " Mean" || $(4 + 3 * i) != s " Variance")
            fail("invalid sample column triple " (i + 1))
        if (s in seen) fail("duplicate sample " s)
        seen[s] = 1
        line = line OFS s ".bam" OFS s ".bam-var"
    }
    print line
    next
}
{
    if (NF != 1 + 3 * n) fail("incorrect field count at row " NR)
    if ($1 !~ /[^[:space:]]/) fail("missing contig name at row " NR)
    total = 0; vals = ""
    for (i = 0; i < n; i++) {
        sample_len = $(2 + 3 * i)
        if (sample_len !~ /^[0-9]+$/ || sample_len + 0 < 1 || sample_len + 0 > 2147483647)
            fail("invalid contig length at row " NR)
        if (sample_len != $2) fail("inconsistent contig lengths at row " NR)
        if (!depth($(3 + 3 * i)) || !depth($(4 + 3 * i))) fail("invalid depth at row " NR)
        total += $(3 + 3 * i)
        vals = vals OFS r4($(3 + 3 * i)) OFS r4($(4 + 3 * i))
    }
    print $1, $2, r4(total / n) vals
}
END { if (!failed && NR < 2) fail(NR ? "no contig rows in CoverM table" : "empty input") }
