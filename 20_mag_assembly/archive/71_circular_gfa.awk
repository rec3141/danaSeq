#!/usr/bin/awk -f
# Usage:
#   awk -f circular_gfa.awk file.gfa
#   awk -v method=selfloop -f circular_gfa.awk file.gfa
#   awk -v method=seteq   -f circular_gfa.awk file.gfa   # default
#
# Outputs FASTA of circular segments (skips '*' sequences).

BEGIN {
    FS = "\t"
    if (method == "") method = "seteq"   # your requested rule
}

# S <name> <sequence> ...
$1 == "S" {
    id = $2
    seq[id] = $3
    next
}

# L <from> <from_orient> <to> <to_orient> <CIGAR>
$1 == "L" {
    from = $2; to = $4
    outgoing[from, to] = 1
    incoming[to, from] = 1
    if (from == to) selfloop[from] = 1   # for optional selfloop mode
    next
}

END {
    if (method == "selfloop") {
        for (s in selfloop) {
            if (s in seq && seq[s] != "*") {
                print ">" s
                printwrap(seq[s], 60)
            }
        }
        exit
    }

    # method == "seteq": neighbor sets equal (ignoring orientation)
    for (s in seq) {
        if (seq[s] == "*" || seq[s] == "") continue

        same = 1
        nin = nout = 0

        # every outgoing neighbor has a matching incoming neighbor
        for (k in outgoing) {
            split(k, a, SUBSEP)
            if (a[1] == s) {
                nout++
                if (!((s, a[2]) in incoming)) same = 0
            }
        }
        # every incoming neighbor has a matching outgoing neighbor
        for (k in incoming) {
            split(k, a, SUBSEP)
            if (a[1] == s) {
                nin++
                if (!((s, a[2]) in outgoing)) same = 0
            }
        }

        if (same && nin == nout && nin > 0) {
            print ">" s
            printwrap(seq[s], 60)
        }
    }
}

function printwrap(str, width, i) {
    for (i = 1; i <= length(str); i += width)
        print substr(str, i, width)
}
