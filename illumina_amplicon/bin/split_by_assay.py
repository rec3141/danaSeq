#!/usr/bin/env python3
"""Split a run that holds more than one amplicon into one stream per assay.

An SRA run is supposed to be one library. Some are not: PRJNA779070, PRJNA895866
and PRJNA942251 each deposit a 16S library and an 18S library concatenated into
one accession — flat for tens of thousands of reads, then a clean switch. The
pipeline cannot treat that as one sample. One assay per sample is what keys the
error models and the truncation lengths, and a 253bp V4 amplicon and a 123bp V9
amplicon need different ones.

It also breaks chimera detection, which is less obvious and worse. Consensus
bimera removal asks whether a sequence is flagged in nearly every sample it
appears in; a sequence from one amplicon is *present* in samples dominated by the
other, where its parents are not, so it cannot be flagged there and no sequence
ever reaches the threshold. PRJNA779070 reports 0 chimeras in 3,271 ASVs while
the pooled method on the same table finds 197.

Each read pair is assigned to the assay whose 5' block its R1 carries, and the
pairs are written to one FASTQ pair per assay. A run carrying a single assay is
passed through untouched and keeps its name, so nothing changes for the runs that
were deposited properly.

Usage:
    split_by_assay.py --primer-set primer_set.json --id SAMPLE \\
        --r1 in_R1.fastq.gz --r2 in_R2.fastq.gz [--min-share 0.05] [--min-reads 500]
"""
import argparse
import gzip
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from primer_db import IUPAC, _ssu_references  # noqa: E402

# How much of the read's 5' end decides which assay it belongs to. The assays
# differ within the first few bases — they start at different places on the gene
# — so this only has to be long enough to be specific, not to span a primer.
_MATCH_LEN = 18
# Absolute identity is not the test — the assays sit hundreds of bases apart on
# the gene, so the read resembles one far more than the other. A floor keeps
# junk out and a margin decides between them.
_MIN_IDENTITY = 0.60
_MIN_MARGIN = 0.15
# The read may open a base or two before the coordinate the assay is named for.
_MAX_START = 4

_SETS = {k: frozenset(v) for k, v in IUPAC.items()}


def _label(assay):
    """Filesystem-safe form of an assay key: ecoli_16S@534-786 -> ecoli_16S_534-786."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", assay)


def _matches(read, block):
    """Best IUPAC identity of an assay's block against the read's 5' window, or 0.

    Early exit per offset: a read from the other assay disagrees within a few
    bases, so the common case costs almost nothing.
    """
    n = min(len(block), _MATCH_LEN, len(read))
    if n < 8:
        return 0.0
    allowed = n  # score every offset; the margin decides, not a hard cutoff
    best = 0.0
    for off in range(_MAX_START + 1):
        if off + n > len(read):
            break
        bad = 0
        for i in range(n):
            if read[off + i] not in _SETS.get(block[i], _SETS["N"]):
                bad += 1
                if bad > allowed:
                    break
        else:
            score = (n - bad) / n
            if score > best:
                best = score
                if best == 1.0:
                    return best
    return best


def _fastq(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            yield h, fh.readline(), fh.readline(), fh.readline()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primer-set", required=True)
    ap.add_argument("--id", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--min-share", type=float, default=0.05,
                    help="assay must hold this fraction of reads to be split out")
    ap.add_argument("--min-reads", type=int, default=500)
    ap.add_argument("--positions", default=None,
                    help="per-stream assay coordinates, one row per emitted sample")
    args = ap.parse_args(argv)

    with open(args.primer_set) as fh:
        pset = json.load(fh)

    # The block for each assay comes off the reference at the coordinate the assay
    # is named for, not from the group's collapsed core. A core is what its
    # members happened to share and has drifted off the start — ecoli_16S@534's
    # longest core begins at 554, twenty bases into the read, and matches nothing
    # anchored there.
    refs = _ssu_references()
    blocks = {}
    for assay in sorted(set((pset.get("assay_of") or {}).values())):
        if "@" not in assay:
            continue
        ref, span = assay.split("@", 1)
        seq = refs.get(ref)
        if not seq:
            continue
        try:
            start = int(span.split("-")[0])
        except ValueError:
            continue
        block = seq[start - 1: start - 1 + _MATCH_LEN]
        if len(block) >= 8:
            blocks[assay] = block

    for k in sorted(blocks):
        print(f"[INFO] {args.id}: {k} -> {blocks[k]}", file=sys.stderr)

    if len(blocks) < 2:
        print(f"[INFO] {args.id}: one assay — passed through", file=sys.stderr)
        _passthrough(args)
        _write_positions(args.positions,
                         [(args.id, next(iter(blocks), None))])
        return 0

    keys = sorted(blocks)
    counts = {k: 0 for k in keys}
    counts[None] = 0
    pairs = []
    for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(_fastq(args.r1), _fastq(args.r2)):
        read = s1.strip().upper()
        scored = sorted(((_matches(read, blocks[k]), k) for k in keys), reverse=True)
        top, runner = scored[0], (scored[1] if len(scored) > 1 else (0.0, None))
        best = top[1] if (top[0] >= _MIN_IDENTITY
                          and top[0] - runner[0] >= _MIN_MARGIN) else None
        counts[best] += 1
        pairs.append((best, (h1, s1, p1, q1), (h2, s2, p2, q2)))

    total = sum(counts.values()) or 1
    keep = [k for k in keys
            if counts[k] >= args.min_reads and counts[k] / total >= args.min_share]

    for k in keys:
        print(f"[INFO] {args.id}: {counts[k]:>8,} reads ({100*counts[k]/total:5.1f}%) "
              f"{k}{'' if k in keep else '  — below threshold, dropped'}",
              file=sys.stderr)
    print(f"[INFO] {args.id}: {counts[None]:>8,} reads ({100*counts[None]/total:5.1f}%) "
          "matched no assay", file=sys.stderr)

    if len(keep) < 2:
        # One assay carries the run after all: keep it whole rather than renaming
        # a sample for a split that did not happen.
        print(f"[INFO] {args.id}: one assay above threshold — passed through",
              file=sys.stderr)
        _passthrough(args)
        _write_positions(args.positions, [(args.id, keep[0] if keep else None)])
        return 0

    handles = {}
    for k in keep:
        name = f"{args.id}__{_label(k)}"
        handles[k] = (gzip.open(f"{name}_R1.fastq.gz", "wt"),
                      gzip.open(f"{name}_R2.fastq.gz", "wt"))
    try:
        for k, r1rec, r2rec in pairs:
            if k not in handles:
                continue
            f1, f2 = handles[k]
            f1.write("".join(r1rec))
            f2.write("".join(r2rec))
    finally:
        for f1, f2 in handles.values():
            f1.close()
            f2.close()

    print(f"[INFO] {args.id}: split into {len(keep)} assays", file=sys.stderr)
    _write_positions(args.positions,
                     [(f"{args.id}__{_label(k)}", k) for k in keep])
    return 0


_POS_COLS = ["sample", "assay_set", "assay_reference", "assay_start", "assay_end"]


def _write_positions(path, rows):
    """One row per emitted stream, so a split sample still reaches samples.json.

    RESOLVE_PRIMER_SET writes these keyed on the run; after a split the run is
    two samples with new names, and the rows have to follow them.
    """
    if not path:
        return
    with open(path, "w") as fh:
        fh.write("\t".join(_POS_COLS) + "\n")
        for sample, assay in rows:
            ref = start = end = ""
            if assay and "@" in assay:
                ref, span = assay.split("@", 1)
                bits = span.split("-")
                start = bits[0]
                end = bits[1] if len(bits) > 1 else ""
            fh.write("\t".join([sample, assay or "", ref, start, end]) + "\n")


def _passthrough(args):
    """Link the reads through under their own name, so the run is unchanged."""
    for src, suffix in ((args.r1, "_R1.fastq.gz"), (args.r2, "_R2.fastq.gz")):
        dst = f"{args.id}{suffix}"
        if os.path.abspath(src) != os.path.abspath(dst):
            if os.path.lexists(dst):
                os.remove(dst)
            os.symlink(os.path.abspath(src), dst)


if __name__ == "__main__":
    sys.exit(main())
