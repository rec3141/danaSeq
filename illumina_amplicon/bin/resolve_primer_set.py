#!/usr/bin/env python3
"""Reduce per-sample primer detections to the one set the whole run is trimmed with.

DETECT_PRIMERS runs per sample, so a run arrives here with one detection each.
Handing all of them to cutadapt would let it pick a different adapter per sample,
and plateFor() keys the AUTO_TRIM and LEARN_ERRORS groups on which one it picked
— so one entry per sample splits a single assay into as many error-model groups
as there are samples. What varies between samples is usually the barcode or
spacer in front of the primer rather than the primer, so compatible detections
are merged and each group keeps only the stretch its members share.

A catalogue sequence is preferred as a group's representative: it has the right
ends, and it already passed the read check in DETECT_PRIMERS, whereas a merged
core is bounded by whatever its members happened to agree on.

Usage:
    resolve_primer_set.py -o primer_set.fa --json primer_set.json *_detected.json
"""
import argparse
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import primer_db  # noqa: E402


def _representatives(records, end):
    """Collapse one end's detections into (sequence, samples, source) entries."""
    seqs = [r.get(end) for r in records if r.get(end)]
    if not seqs:
        return []
    catalogue = {r[end] for r in records
                 if r.get(end) and r.get("source") == "inferred-catalogue"}
    # Two consensus groups can resolve to the same catalogue primer, so the
    # sequences are merged after substitution rather than before it. Leaving them
    # separate puts the same adapter in the FASTA twice, and plateFor() keys the
    # error-model groups on the adapter cutadapt matched — so one assay would be
    # split in two by an entry duplicated against itself.
    merged: dict[str, dict] = {}
    for core, n in primer_db.collapse_primers(seqs):
        members = [s for s in seqs if primer_db.shared_region(s, core)]
        hits = [s for s in members if s in catalogue]
        if hits:
            seq, source = Counter(hits).most_common(1)[0][0], "catalogue"
        else:
            seq, source = core, "de-novo"
        if seq in merged:
            merged[seq]["samples"] += n
        else:
            merged[seq] = {"sequence": seq, "samples": n, "source": source,
                           "location": primer_db.locate_on_ssu(seq)}
    return sorted(merged.values(), key=lambda e: -e["samples"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("detections", nargs="+")
    ap.add_argument("-o", "--out", default="primer_set.fa")
    ap.add_argument("--json", dest="json_out", default="primer_set.json")
    args = ap.parse_args(argv)

    records = []
    for path in args.detections:
        try:
            with open(path) as fh:
                records.append(json.load(fh))
        except (OSError, ValueError) as exc:
            print(f"[WARN] unreadable detection {path}: {exc}", file=sys.stderr)

    if not records:
        print("[ERROR] no usable primer detections", file=sys.stderr)
        return 1

    fwd = _representatives(records, "fwd")
    rev = _representatives(records, "rev")
    entries = fwd + rev
    if not entries:
        print("[ERROR] detections carried no primer sequences", file=sys.stderr)
        return 1

    # A run whose primers sit on neither SSU reference is not a ribosomal assay,
    # so the taxonomy step has nothing meaningful to say about it.
    ribosomal = any(e["location"] for e in entries)

    # The header is the sequence: it is what cutadapt reports back in its log,
    # and therefore what the assay is keyed on downstream.
    with open(args.out, "w") as fh:
        for e in entries:
            fh.write(f">{e['sequence']}\n{e['sequence']}\n")

    summary = {"samples": len(records), "ribosomal": ribosomal,
               "fwd": fwd, "rev": rev}
    with open(args.json_out, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[INFO] {len(records)} samples -> {len(entries)} primer(s)", file=sys.stderr)
    for e in entries:
        where = ", ".join(f"{l['reference']}@{l['start']}" for l in e["location"])
        print(f"[INFO]   {e['sequence']} ({e['source']}, {e['samples']} samples) "
              f"{where or 'no SSU match'}", file=sys.stderr)
    if not ribosomal:
        print("[WARN] no primer matches the 16S or 18S reference — this run is not "
              "a ribosomal assay, and an rRNA taxonomy database would label it "
              "confidently and wrongly", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
