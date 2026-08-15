#!/usr/bin/env python3
"""Split a bundled primer-pair FASTA into one file per end.

primers-<FWD>-<REV>.fa carries the forward record first and the reverse second.
Record order is the only thing that says which end a record belongs to — the
headers hold primer names, and a name does not say which end it primes.

A file with more than two records cannot be split this way and is copied to both
ends unchanged, which is what cutadapt was given before there were two files:
every adapter offered to both reads, the best match per read winning.

Usage:
    split_primer_pair.py primers-515F-806RB.fa out_fwd.fa out_rev.fa
"""
import sys


def records(text):
    out, cur = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                out.append("\n".join(cur))
            cur = [line]
        elif cur:
            cur.append(line)
    if cur:
        out.append("\n".join(cur))
    return out


def main(argv):
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    src, fwd_out, rev_out = argv[1], argv[2], argv[3]
    with open(src) as fh:
        text = fh.read()
    recs = records(text)

    if len(recs) == 2:
        halves = [recs[0], recs[1]]
    else:
        print(f"[WARN] {src} holds {len(recs)} records, not a pair — offering all "
              "of them to both ends", file=sys.stderr)
        halves = [text.strip(), text.strip()]

    for path, body in ((fwd_out, halves[0]), (rev_out, halves[1])):
        with open(path, "w") as fh:
            fh.write(body.rstrip("\n") + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
