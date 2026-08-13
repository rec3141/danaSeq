#!/usr/bin/env python3
"""Detect the primer pair a sample carries, from the reads themselves.

Samples a few hundred reads once and matches their 5' ends against the curated
table (bin/primer_db.py), falling back to a de-novo degenerate consensus when
nothing scores well. Writes the primers it found as a FASTA for REMOVE_PRIMERS,
so trimming is not limited to the pairs that happen to have a file in primers/.

Because the table travels with the pipeline, the result names the gene and the
lineage it belongs to — "16S" alone is the prokaryotic SSU rRNA to one field and
the mitochondrial LSU rRNA to another, so an assay record that just says "16S" is
not interpretable on its own.

Usage:
    detect_primers.py R1.fastq.gz [R2.fastq.gz] [-o detected_primers.fa]
                      [--json detected_primers.json] [-n 500]
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import primer_db  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("r1")
    ap.add_argument("r2", nargs="?")
    ap.add_argument("-o", "--out", default="detected_primers.fa")
    ap.add_argument("--json", dest="json_out", default=None)
    ap.add_argument("-n", "--reads", type=int, default=500,
                    help="reads to sample per file [default: 500]")
    args = ap.parse_args(argv)

    hit = primer_db.detect_from_reads(args.r1, args.r2, n=args.reads)

    if not hit:
        # Write nothing rather than an empty FASTA: cutadapt with an empty
        # adapter file and --discard-untrimmed silently drops every read, which
        # looks like a clean run that produced no data.
        print("[WARN] no primer detected — too few reads, or no usable 5' "
              "consensus", file=sys.stderr)
        return 1

    described = primer_db.describe_pair(hit.get("fwd_name"), hit.get("rev_name")) or {}
    record = {**hit, **{f"assay_{k}": v for k, v in described.items()}}

    with open(args.out, "w") as fh:
        fh.write(f">{hit.get('fwd_name') or 'fwd'}\n{hit['fwd']}\n")
        if hit.get("rev"):
            fh.write(f">{hit.get('rev_name') or 'rev'}\n{hit['rev']}\n")

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(record, fh, indent=2)

    gene = described.get("gene", "unknown gene")
    lineage = described.get("lineage")
    print(f"[INFO] detected {hit.get('fwd_name')}/{hit.get('rev_name')} "
          f"({gene}{', ' + lineage if lineage else ''}) "
          f"via {hit.get('source')} confidence={hit.get('confidence')}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
