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

The two ends are written separately, and an end whose reads already begin past
where a primer would sit gets an empty file — a study can reach SRA trimmed on
one end and not the other, and one set for both ends can only be wrong for one
of them.

Usage:
    resolve_primer_set.py --fwd-out primer_set_fwd.fa --rev-out primer_set_rev.fa \
        --json primer_set.json *_detected.json
"""
import argparse
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import primer_db  # noqa: E402


def _representatives(records, end):
    """Collapse one end's detections into entries, and say which sample is in which.

    Group membership is the run's assay key. It is decided here, on the primer
    sequences themselves, which works whatever the amplicon is — cutadapt's own
    per-adapter counts answer the same question, but only for a run that is
    trimmed, and only for one that is trimmed with something real.
    """
    pairs = [(r.get("sample"), r[end], r.get(f"{end}_state", "unplaced"),
              r.get(f"{end}_abuts") or [])
             for r in records if r.get(end)]
    if not pairs:
        return []
    seqs = [s for _, s, _, _ in pairs]
    # Both "catalogue" and "catalogue-scan" hand back a sequence copied from the
    # table, so both earn the same protection from being shortened by collapsing.
    catalogue = {r[end] for r in records
                 if r.get(end)
                 and str(r.get(f"{end}_source", "")).startswith("inferred-catalogue")}

    # A sequence can share a core with more than one group, so it goes to the
    # first that takes it — the same first-match rule collapse_primers built the
    # groups under. A sample in two groups is not a key.
    cores = [core for core, _ in primer_db.collapse_primers(seqs)]
    membership: dict[str, list] = {core: [] for core in cores}
    for sample, seq, state, abuts in pairs:
        for core in cores:
            if primer_db.shared_region(seq, core):
                membership[core].append((sample, seq, state, abuts))
                break

    # Two consensus groups can resolve to the same catalogue primer, so the
    # sequences are merged after substitution rather than before it. Leaving them
    # separate puts the same adapter in the FASTA twice, and plateFor() keys the
    # error-model groups on the adapter cutadapt matched — so one assay would be
    # split in two by an entry duplicated against itself.
    merged: dict[str, dict] = {}
    for core in cores:
        members = membership[core]
        hits = [s for _, s, _, _ in members if s in catalogue]
        if hits:
            seq, source = Counter(hits).most_common(1)[0][0], "catalogue"
        else:
            seq, source = core, "de-novo"
        samples = [smp for smp, _, _, _ in members if smp]
        # The group's state is its members' verdict, taken per sample before
        # collapsing moved the edge off the boundary it was measured against.
        state = Counter(st for _, _, st, _ in members).most_common(1)[0][0]
        abuts = sorted({n for _, _, st, ab in members if st == "pre-trimmed"
                        for n in ab})
        if seq in merged:
            merged[seq]["samples"] += len(members)
            merged[seq]["members"].extend(samples)
        else:
            merged[seq] = {"sequence": seq, "samples": len(members), "source": source,
                           "state": state,
                           "location": primer_db.locate_on_ssu(seq),
                           "abuts": abuts,
                           "members": samples}
    return sorted(merged.values(), key=lambda e: -e["samples"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("detections", nargs="+")
    ap.add_argument("--fwd-out", default="primer_set_fwd.fa")
    ap.add_argument("--rev-out", default="primer_set_rev.fa")
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

    # Trimming is a separate question from grouping, and it is answered here.
    # An end whose reads begin at a cut rather than at a primer has nothing left
    # to remove, and handing it to cutadapt takes the front of the amplicon
    # instead. An end that sits on no reference cannot be judged this way and does
    # not block trimming: that is every non-ribosomal assay, whose primer is real
    # and still attached.
    #
    # Per end, weighted by samples rather than by entries. A run splinters into
    # singleton groups wherever one sample's consensus drifted, and counting
    # groups lets three stray samples outvote sixty.
    def _cut(end_entries):
        placed = [e for e in end_entries if e["state"] != "unplaced"]
        if not placed:
            return None
        cut = sum(e["samples"] for e in placed if e["state"] == "pre-trimmed")
        return cut > sum(e["samples"] for e in placed) / 2

    # Per end, so a run trimmed on one end only keeps the primer it still has
    # and keeps the amplicon it no longer has a primer in front of.
    fwd_cut, rev_cut = _cut(fwd), _cut(rev)
    trim_fwd = fwd_cut is not None and not fwd_cut
    trim_rev = rev_cut is not None and not rev_cut
    pre_trimmed = bool(fwd_cut or rev_cut)
    cut_at = [e for e in entries if e["state"] == "pre-trimmed"]

    # Which group each sample fell into — the assay key, derived from the primer
    # sequences rather than from cutadapt's per-adapter counts, so it survives a
    # run that is never handed to cutadapt at all.
    assay_of = {}
    for e in fwd:
        for smp in e["members"]:
            assay_of[smp] = e["sequence"]

    # One file per end, because the ends are trimmed independently. An end with
    # nothing to remove gets an empty file, and REMOVE_PRIMERS passes cutadapt no
    # adapter option for it — an empty adapter file with --discard-untrimmed
    # matches nothing and would drop every read.
    #
    # The header is the sequence: it is what cutadapt reports back in its log,
    # and therefore what the assay is keyed on downstream.
    for path, end_entries, wanted in ((args.fwd_out, fwd, trim_fwd),
                                      (args.rev_out, rev, trim_rev)):
        with open(path, "w") as fh:
            if wanted:
                for e in end_entries:
                    fh.write(f">{e['sequence']}\n{e['sequence']}\n")

    # Reads that arrive in a random orientation carry the forward primer on half
    # of R1 and the reverse on the other half, so neither end can be trimmed as
    # an end. The consensus still resolves — to whichever orientation won the
    # majority — and quietly describes half the run.
    def _median(key):
        vals = sorted(v for v in (r.get(key) for r in records) if v is not None)
        return vals[len(vals) // 2] if vals else None

    # A primer belongs to one end. Finding the forward primer in R2, and the
    # reverse in R1, at the same time and in a good share of reads, is not a
    # degenerate match — it is both orientations in both files.
    fwd_in_r2, rev_in_r1 = _median("fwd_in_r2"), _median("rev_in_r1")
    scrambled = (fwd_in_r2 is not None and rev_in_r1 is not None
                 and fwd_in_r2 >= 0.25 and rev_in_r1 >= 0.25)

    summary = {"samples": len(records), "ribosomal": ribosomal,
               "pre_trimmed": pre_trimmed,
               "trim_fwd": trim_fwd, "trim_rev": trim_rev,
               "cross_orientation": {"fwd_in_r2": fwd_in_r2, "rev_in_r1": rev_in_r1},
               "mixed_orientation": scrambled,
               "assay_of": assay_of, "fwd": fwd, "rev": rev}
    with open(args.json_out, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[INFO] {len(records)} samples -> {len(entries)} primer(s)", file=sys.stderr)
    for e in entries:
        where = ", ".join(f"{l['reference']}@{l['start']}-{l['end']}"
                          for l in e["location"])
        print(f"[INFO]   {e['sequence']} ({e['source']}, {e['state']}, "
              f"{e['samples']} samples) {where or 'no SSU match'}", file=sys.stderr)
    if pre_trimmed:
        # Only the ends that carry a cut, and only the reference each one is
        # placed best on — a V4 primer lands on both genes, and quoting both
        # coordinates for one boundary reads as two boundaries.
        spans = sorted({f"{e['location'][0]['reference']} "
                        f"{e['location'][0]['start']}-{e['location'][0]['end']}"
                        for e in cut_at if e["location"]})
        names = sorted({n for e in cut_at for n in e["abuts"]})
        print("[WARN] these reads begin where a primer ends, not where one starts: "
              "they were trimmed before submission, so nothing is removed here. "
              "Amplicon spans " + "; ".join(spans), file=sys.stderr)
        print(f"[WARN] which primer was used is not recoverable — {len(names)} "
              "catalogued primers share that boundary, including "
              + ", ".join(names[:6]), file=sys.stderr)
    if scrambled:
        print(f"[WARN] the forward primer is on {fwd_in_r2:.0%} of R2 and the reverse "
              f"on {rev_in_r1:.0%} of R1 — these reads are in mixed orientation, so "
              "neither file is an end. -g is given the forward primer and -G the "
              "reverse, so the half that is the other way round matches no adapter "
              "and --discard-untrimmed drops it", file=sys.stderr)
    if not ribosomal:
        print("[WARN] no primer matches the 16S or 18S reference — this run is not "
              "a ribosomal assay, and an rRNA taxonomy database would label it "
              "confidently and wrongly", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
