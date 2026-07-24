#!/usr/bin/env python3
"""Turn cutadapt's own logs into a per-sample record of which primer amplified what.

REMOVE_PRIMERS runs cutadapt with `-g file:` against every candidate primer at
once, and cutadapt already reports, per adapter, how many reads it trimmed. That
is the ground truth for which assay a sample actually carries — better than the
submitter's label, which is wrong often enough to matter (PRJNA1473294 labels all
84 runs "16S" when 40 are eukaryotic 18S). Nothing read those logs, so the answer
was computed and thrown away on every run.

This emits one row per sample: the primer names cutadapt actually matched, the
read counts, and — resolved through the curated table in primer_db.py — which
gene those primers target and whose lineage it belongs to. The gene matters more
than the region, and naming it precisely matters more still: "16S" alone is the
prokaryotic SSU rRNA to a microbial ecologist and the mitochondrial LSU rRNA to a
zoologist, so a record that just says "16S" cannot be read safely.

Names that resolve to no table entry keep their counts and simply carry no gene,
rather than being dropped or guessed at.

Usage:
    parse_primer_assignment.py OUT.tsv LOG [LOG ...]
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from primer_db import describe_pair
except ImportError:  # table unavailable — still report what was observed
    describe_pair = lambda *_a, **_k: None

# "=== First read: Adapter 341Fv3 ===" / "=== Second read: Adapter Bakt_805R ==="
# and the single-end form "=== Adapter 341Fv3 ===".
_ADAPTER_RE = re.compile(r"^===\s*(?:(First|Second) read:\s*)?Adapter\s+(.+?)\s*===\s*$")
_TRIMMED_RE = re.compile(r"Trimmed:\s*([\d,]+)\s*times")
_TOTAL_RE = re.compile(r"^Total read(?: pair)?s processed:\s*([\d,]+)", re.M)

COLUMNS = [
    "sample",
    "assay_gene",
    "assay_gene_lineage",
    "assay_region",
    "assay_primer_fwd",
    "assay_primer_rev",
    "assay_reads_in",
    "assay_reads_matched",
    "assay_match_fraction",
    "assay_source",
]


def _int(text):
    return int(text.replace(",", ""))


def parse_log(text):
    """Return the winning adapter per read direction, plus read counts.

    The winner is the adapter cutadapt trimmed most often. Ties and zero-match
    samples resolve to None rather than an arbitrary pick — a sample nothing
    matched has no assay to report, which is itself worth seeing.
    """
    total_match = _TOTAL_RE.search(text)
    reads_in = _int(total_match.group(1)) if total_match else None

    # direction -> {adapter_name: trimmed_count}
    counts = {"First": {}, "Second": {}}
    current = None
    for line in text.splitlines():
        header = _ADAPTER_RE.match(line.strip())
        if header:
            # Single-end logs omit the direction; treat them as the forward read.
            current = (header.group(1) or "First", header.group(2))
            continue
        if current:
            trimmed = _TRIMMED_RE.search(line)
            if trimmed:
                direction, name = current
                counts[direction][name] = counts[direction].get(name, 0) + _int(trimmed.group(1))
                current = None

    def winner(direction):
        hits = {k: v for k, v in counts[direction].items() if v > 0}
        if not hits:
            return None, 0
        top = max(hits.values())
        best = [k for k, v in hits.items() if v == top]
        return (best[0] if len(best) == 1 else None), top

    fwd, fwd_n = winner("First")
    rev, _rev_n = winner("Second")
    return {
        "assay_primer_fwd": fwd,
        "assay_primer_rev": rev,
        "assay_reads_in": reads_in,
        "assay_reads_matched": fwd_n or None,
        "assay_match_fraction": (round(fwd_n / reads_in, 4)
                                 if reads_in and fwd_n else None),
        "assay_source": "cutadapt",
    }


def sample_id(path):
    """`SRR38958117_cutadapt.log` -> `SRR38958117`, matching the pipeline's ids."""
    name = os.path.basename(path)
    for suffix in ("_cutadapt.log", ".cutadapt.log", ".log"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def main(argv):
    if len(argv) < 3:
        sys.exit(__doc__.strip())
    out_path, logs = argv[1], argv[2:]

    rows = []
    for path in logs:
        try:
            with open(path) as fh:
                rec = parse_log(fh.read())
        except OSError as e:
            print(f"[WARN] unreadable log {path}: {e}", file=sys.stderr)
            continue
        described = describe_pair(rec["assay_primer_fwd"], rec["assay_primer_rev"]) or {}
        rec["assay_gene"] = described.get("gene")
        rec["assay_gene_lineage"] = described.get("lineage")
        rec["assay_region"] = described.get("region")
        rec["sample"] = sample_id(path)
        rows.append(rec)

    rows.sort(key=lambda r: r["sample"])
    with open(out_path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c))
                               for c in COLUMNS) + "\n")

    resolved = sum(1 for r in rows if r["assay_primer_fwd"])
    named = sum(1 for r in rows if r.get("assay_gene"))
    print(f"[INFO] primer assignment: {resolved}/{len(rows)} samples matched a "
          f"primer, {named} resolved to a gene -> {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
