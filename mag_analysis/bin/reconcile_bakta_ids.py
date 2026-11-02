#!/usr/bin/env python3
"""Give the Bakta extra run the same locus tags as the Bakta basic run.

BAKTA_BASIC and BAKTA_EXTRA call the same CDS, but each run numbers its locus
tags independently, and the extra run also tags sORFs. The tags therefore drift
apart after the first sORF. Every protein tool (KofamScan, eggNOG, dbCAN,
DefenseFinder, ...) runs on the basic proteins, so anything that joins those
results to extra-run genes by locus tag gets another gene's function.

This script rewrites the extra-run files so that each feature with the same
(contig, start, stop, strand) as a basic CDS takes the basic locus tag.
Features with no basic counterpart (sORFs, and any CDS only the extra run
called) get new tags numbered after the highest basic tag. Running it on
files that were already reconciled leaves them unchanged.

Usage:
    reconcile_bakta_ids.py --basic-tsv basic/annotation.tsv \\
        --extra-dir extra_raw/ --outdir . [--files annotation.tsv annotation.gff3 ...]

Writes the rewritten files plus locus_tag_map.tsv (extra_tag, new_tag,
matched) to --outdir.
"""

import argparse
import os
import re
import sys
from array import array

DEFAULT_FILES = ['annotation.tsv', 'annotation.gff3', 'annotation.faa',
                 'annotation.hypotheticals.faa', 'annotation.gbff']
TAG_RE = re.compile(r'^(.+)_(\d+)$')


def split_tag(tag):
    m = TAG_RE.match(tag)
    if not m:
        sys.exit(f"ERROR: unexpected locus tag format: {tag!r}")
    return m.group(1), m.group(2)


def read_features(path):
    """Yield (contig, start, stop, strand, type, locus_tag) for tagged features."""
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 6 or not f[5]:
                continue
            yield f[0], int(f[2]), int(f[3]), f[4], f[1], f[5]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--basic-tsv', required=True)
    ap.add_argument('--extra-dir', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--files', nargs='+', default=DEFAULT_FILES,
                    help='extra-run files to rewrite (missing ones are skipped)')
    args = ap.parse_args()

    # Pack (contig, strand, start, stop) into one int to keep 10M+ keys affordable.
    contig_idx = {}
    def key(contig, start, stop, strand):
        ci = contig_idx.setdefault(contig, len(contig_idx))
        return (((ci << 1 | (strand == '-')) << 40 | start) << 40) | stop

    basic = {}
    b_prefix = b_width = None
    b_max = 0
    for contig, start, stop, strand, ftype, tag in read_features(args.basic_tsv):
        prefix, num = split_tag(tag)
        if b_prefix is None:
            b_prefix, b_width = prefix, len(num)
        elif prefix != b_prefix:
            sys.exit(f"ERROR: mixed locus tag prefixes in basic TSV: {b_prefix}, {prefix}")
        n = int(num)
        b_max = max(b_max, n)
        basic[key(contig, start, stop, strand)] = n
    if b_prefix is None:
        sys.exit(f"ERROR: no locus tags in {args.basic_tsv}")

    # extra number -> new number (-1 = not a tag of this run)
    extra_tsv = os.path.join(args.extra_dir, 'annotation.tsv')
    os.makedirs(args.outdir, exist_ok=True)
    fmt = f"{b_prefix}_{{:0{b_width}d}}"
    remap = array('q')
    e_prefix = None
    next_new = b_max
    matched = unmatched_cds = 0
    unmatched_types = {}
    with open(os.path.join(args.outdir, 'locus_tag_map.tsv'), 'w') as map_out:
        map_out.write('extra_tag\tnew_tag\tmatched\n')
        for contig, start, stop, strand, ftype, tag in read_features(extra_tsv):
            prefix, num = split_tag(tag)
            if e_prefix is None:
                e_prefix, e_width = prefix, len(num)
            elif prefix != e_prefix:
                sys.exit(f"ERROR: mixed locus tag prefixes in extra TSV: {e_prefix}, {prefix}")
            n = int(num)
            if n >= len(remap):
                remap.extend([-1] * (n + 1 - len(remap)))
            # pop: each basic tag is handed out at most once
            new = basic.pop(key(contig, start, stop, strand), None) if ftype == 'cds' else None
            if new is not None:
                matched += 1
            else:
                next_new += 1
                new = next_new
                unmatched_types[ftype] = unmatched_types.get(ftype, 0) + 1
                unmatched_cds += ftype == 'cds'
            remap[n] = new
            map_out.write(f"{tag}\t{fmt.format(new)}\t{int(new <= b_max)}\n")
    del basic
    if e_prefix is None:
        sys.exit(f"ERROR: no locus tags in {extra_tsv}")

    # A tag is the prefix plus a fixed number of digits, not inside a longer word or number
    # (Bakta also writes IDs like <tag>_gene and gnl|Bakta|<tag>).
    tok = re.compile(r'(?<![A-Za-z0-9])' + re.escape(e_prefix) + rf'_(\d{{{e_width}}})(?!\d)')
    def sub(m):
        n = int(m.group(1))
        new = remap[n] if n < len(remap) else -1
        if new < 0:
            sys.exit(f"ERROR: locus tag {m.group(0)} is not in {extra_tsv}")
        return fmt.format(new)

    for name in args.files:
        src = os.path.join(args.extra_dir, name)
        if not os.path.exists(src):
            continue
        dst = os.path.join(args.outdir, name)
        if os.path.abspath(src) == os.path.abspath(dst):
            sys.exit(f"ERROR: --outdir must differ from --extra-dir ({name})")
        tmp = dst + '.tmp'
        with open(src) as fh, open(tmp, 'w') as out:
            for line in fh:
                if e_prefix in line:
                    line = tok.sub(sub, line)
                out.write(line)
        os.replace(tmp, dst)
        print(f"  rewrote {name}", file=sys.stderr)

    total = matched + sum(unmatched_types.values())
    print(f"[reconcile_bakta_ids] {total} extra features: {matched} took basic tags, "
          f"{total - matched} new tags after {fmt.format(b_max)} "
          f"({', '.join(f'{v} {k}' for k, v in sorted(unmatched_types.items())) or 'none'})",
          file=sys.stderr)
    if unmatched_cds > 0.01 * max(matched, 1):
        print(f"[WARNING] {unmatched_cds} extra-run CDS have no basic CDS at the same coordinates; "
              "protein-tool results (KofamScan, eggNOG, ...) will be missing for them.",
              file=sys.stderr)


if __name__ == '__main__':
    main()
