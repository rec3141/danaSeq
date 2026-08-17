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

Reads that reach us with their primers already removed produce no cutadapt log
at all, and then the coordinates from RESOLVE_PRIMER_SET are the whole record:
where on the SSU reference each sample's amplicon begins and ends. Those columns
are merged in whether or not there are logs, so a run that is never trimmed still
reports what each sample carries.

Usage:
    parse_primer_assignment.py OUT.tsv [--positions POS.tsv] [LOG ...]
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

# Coordinates first when they exist, because for a pre-trimmed run they are the
# only description of the assay there is.
POSITION_COLUMNS = [
    "assay_reference",
    "assay_start",
    "assay_end",
]

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


def read_positions(path):
    """Per-sample coordinates from RESOLVE_PRIMER_SET, keyed by sample."""
    out = {}
    try:
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                if not line.strip():
                    continue
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                sample = row.get("sample")
                if sample:
                    out[sample] = row
    except OSError as e:
        print(f"[WARN] unreadable positions {path}: {e}", file=sys.stderr)
    return out


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__.strip())
    out_path = argv[1]
    args = argv[2:]
    positions = {}
    if args and args[0] == "--positions":
        if len(args) < 2:
            sys.exit("--positions needs a file")
        positions = read_positions(args[1])
        args = args[2:]
    logs = args

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

    # Samples the logs never mentioned — every sample, on a run that was not
    # trimmed — still get a row from their coordinates.
    seen = {r["sample"] for r in rows}
    for sample, pos in positions.items():
        if sample not in seen:
            rows.append({"sample": sample,
                         "assay_primer_fwd": pos.get("assay_primer_fwd") or None,
                         "assay_primer_rev": pos.get("assay_primer_rev") or None})
    for r in rows:
        pos = positions.get(r["sample"])
        if pos:
            for c in POSITION_COLUMNS:
                r[c] = pos.get(c) or None

    columns = COLUMNS + [c for c in POSITION_COLUMNS if any(r.get(c) for r in rows)]
    rows.sort(key=lambda r: r["sample"])
    with open(out_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c))
                               for c in columns) + "\n")

    resolved = sum(1 for r in rows if r.get("assay_primer_fwd"))
    named = sum(1 for r in rows if r.get("assay_gene"))
    placed = sum(1 for r in rows if r.get("assay_reference"))
    print(f"[INFO] primer assignment: {resolved}/{len(rows)} samples matched a "
          f"primer, {named} resolved to a gene, {placed} placed on a reference "
          f"-> {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
