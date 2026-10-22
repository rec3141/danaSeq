#!/usr/bin/env python3
"""
Read name-grouped SAM (strobealign --interleaved, unmapped reads included) on
stdin and write the non-human pairs as interleaved FASTQ on stdout.

A mate is human when it aligns over at least --min-aligned of its length at
--min-identity or better (NM over the aligned bases), the stringency of BBTools
removehuman.sh. A pair goes if either mate is human. Reads are restored to
their original orientation. Prints "<kept> <removed>" read counts to --counts.
Stdlib only.
"""
import argparse, re, sys

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def is_human(f, min_aligned, min_identity):
    flag = int(f[1])
    if flag & 4:
        return False
    aligned = sum(int(n) for n, op in CIGAR.findall(f[5]) if op in "MI=X")
    if aligned < min_aligned * len(f[9]):
        return False
    nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), None)
    if nm is None:
        return True
    return nm <= (1 - min_identity) * aligned


def fastq(f):
    seq, qual = f[9], f[10]
    if int(f[1]) & 16:
        seq, qual = seq.translate(COMP)[::-1], qual[::-1]
    return f"@{f[0]}\n{seq}\n+\n{qual}\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-aligned", type=float, default=0.8)
    ap.add_argument("--min-identity", type=float, default=0.95)
    ap.add_argument("--counts", required=True)
    a = ap.parse_args()

    out = sys.stdout
    kept = removed = 0
    pending = None
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if int(f[1]) & 0x900:          # secondary / supplementary
            continue
        if pending is None:
            pending = f
            continue
        if pending[0] != f[0]:
            # unpaired leftover: judge it alone
            if is_human(pending, a.min_aligned, a.min_identity):
                removed += 1
            else:
                out.write(fastq(pending)); kept += 1
            pending = f
            continue
        r1, r2 = (pending, f) if int(pending[1]) & 64 else (f, pending)
        pending = None
        if is_human(r1, a.min_aligned, a.min_identity) or is_human(r2, a.min_aligned, a.min_identity):
            removed += 2
        else:
            out.write(fastq(r1)); out.write(fastq(r2)); kept += 2
    if pending is not None:
        if is_human(pending, a.min_aligned, a.min_identity):
            removed += 1
        else:
            out.write(fastq(pending)); kept += 1
    with open(a.counts, "w") as fh:
        fh.write(f"{kept}\t{removed}\n")


if __name__ == "__main__":
    main()
