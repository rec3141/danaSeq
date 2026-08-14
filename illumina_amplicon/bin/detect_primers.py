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
        # Report nothing and succeed. A sample is trimmed with the set resolved
        # across the whole run, not with its own detection, so one that yields no
        # consensus still gets the primer its assay-mates supplied; exiting non-
        # zero here would drop it from the run instead.
        print("[WARN] no primer detected — too few reads, or no usable 5' "
              "consensus", file=sys.stderr)
        return 0

    described = primer_db.describe_pair(hit.get("fwd_name"), hit.get("rev_name")) or {}
    record = {**hit, **{f"assay_{k}": v for k, v in described.items()}}

    # The FASTA header is the primer's own sequence, because it is the only field
    # that survives the trip downstream — cutadapt reports the header it matched
    # and nothing else — and because plateFor() keys the AUTO_TRIM and
    # LEARN_ERRORS groups on it. Keying on a sequence groups samples by what was
    # actually trimmed; keying on a name splits them whenever two synonyms in the
    # catalogue describe the same primer.
    with open(args.out, "w") as fh:
        fh.write(f">{hit['fwd']}\n{hit['fwd']}\n")
        if hit.get("rev"):
            fh.write(f">{hit['rev']}\n{hit['rev']}\n")

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump(record, fh, indent=2)

    def _where(locs):
        return ", ".join(f"{l['reference']}@{l['start']}" for l in locs) or "no SSU match"

    print(f"[INFO] {hit['fwd']} ({hit.get('fwd_source')}, "
          f"{_where(hit.get('fwd_location') or [])}) / "
          f"{hit.get('rev') or '-'} ({hit.get('rev_source')}, "
          f"{_where(hit.get('rev_location') or [])})", file=sys.stderr)
    if not hit.get("ribosomal", True):
        print("[WARN] primers match neither the 16S nor the 18S reference — this "
              "assay is not ribosomal, so an rRNA taxonomy database would return "
              "confident nonsense for it", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
