#!/usr/bin/env python3
"""Per-sample read accounting across every stage of the pipeline.

Builds one table tracing each sample's reads from the raw FASTQ all the way to
the final ASV table, so where (and why) a sample loses reads — or drops out
entirely — is visible at a glance. Counting starts at the raw FASTQ rather than
at the stats files, which begin post-cutadapt: without that, a sample cutadapt
discarded is indistinguishable from a shallow one. Stages tracked:

    raw      Raw read pairs                 (cutadapt logs)
    primer   After primer removal           (cutadapt logs)
    filter   After DADA2 quality filter      (*_filt_stats.tsv)
    denoise  After denoise + merge pairs     (seqtabs/*.seqtab.tsv, per ROW)
    merge    After combining seqtabs         (merge_sample_reads.tsv)
    chimera  After chimera removal           (chimera_sample_reads.tsv)
    final    Final filtered table            (final_sample_reads.tsv)

Comparing `denoise` (every published seqtab) against `merge` exposes any sample
dropped between the two. A seqtab file is NOT one sample: error models are keyed
on (run, assay), so one file holds every sample in its group and samples are read
off the row labels.

Usage: read_tracking.py <out.tsv> [cutadapt_logs...] [filt_stats...]

The seqtabs/ dir and the seqtab_final/*_sample_reads.tsv files are discovered
relative to <out.tsv> (expected at <outdir>/seqtab_final/read_tracking.tsv), so
only the cutadapt logs and filter stats need to be passed explicitly.
"""
import glob
import json
import os
import re
import sys

RAW_RE = re.compile(r"Total read pairs processed:\s*([\d,]+)")
RAW_SE_RE = re.compile(r"Total reads processed:\s*([\d,]+)")
KEPT_RE = re.compile(r"Pairs written \(passing filters\):\s*([\d,]+)")
KEPT_SE_RE = re.compile(r"Reads written \(passing filters\):\s*([\d,]+)")

STAGES = [
    ("raw",     "Raw FASTQ",             "reads_raw"),
    ("primer",  "After primer removal",  "reads_after_primer"),
    ("filter",  "After quality filter",  "reads_after_filter"),
    ("denoise", "After denoising",       "reads_after_denoise"),
    ("merge",   "After merge",           "reads_after_merge"),
    ("chimera", "After chimera removal",  "reads_after_chimera"),
    ("final",   "After table filter",     "reads_final"),
]


def _num(m):
    return int(m.group(1).replace(",", "")) if m else 0


def _seqtab_row_totals(path):
    """Per-sample read totals from a wide seqtab TSV (row = sample, cols = ASVs).

    Keyed on the ROW LABEL, never the filename: DENOISE publishes one seqtab per
    (run, assay) GROUP, whose name is not a sample at all.
    """
    totals = {}
    try:
        with open(path) as fh:
            next(fh, None)  # header: sample, then one column per ASV
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                n = 0
                for cell in parts[1:]:
                    try:
                        n += int(float(cell))
                    except ValueError:
                        pass
                totals[parts[0]] = totals.get(parts[0], 0) + n
    except OSError:
        pass
    return totals


def _read_sample_reads_tsv(path, samples, col):
    """Load a `sample<TAB>reads` TSV into samples[acc][col]."""
    try:
        with open(path) as fh:
            header = fh.readline()  # sample\treads
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                acc, val = parts[0], parts[1]
                try:
                    samples.setdefault(acc, {})[col] = int(float(val))
                except ValueError:
                    pass
    except OSError:
        pass


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: read_tracking.py <out.tsv> [inputs...]")
    out_path = sys.argv[1]
    inputs = sys.argv[2:] or (sorted(glob.glob("*_cutadapt.log"))
                              + sorted(glob.glob("*_filt_stats.tsv")))

    seqtab_final_dir = os.path.dirname(out_path) or "."
    outdir = os.path.dirname(seqtab_final_dir) or "."
    seqtabs_dir = os.path.join(outdir, "seqtabs")
    viz_dir = os.path.join(outdir, "viz")

    samples = {}

    # 1-3: cutadapt logs (raw, primer) + filter stats (filter)
    for f in inputs:
        base = os.path.basename(f)
        if base.endswith("_cutadapt.log"):
            acc = base[: -len("_cutadapt.log")]
            try:
                txt = open(f, errors="ignore").read()
            except OSError:
                continue
            d = samples.setdefault(acc, {})
            d["reads_raw"] = _num(RAW_RE.search(txt)) or _num(RAW_SE_RE.search(txt))
            d["reads_after_primer"] = _num(KEPT_RE.search(txt)) or _num(KEPT_SE_RE.search(txt))
        elif base.endswith("_filt_stats.tsv"):
            acc = base[: -len("_filt_stats.tsv")]
            try:
                lines = [l.rstrip("\n") for l in open(f) if l.strip()]
            except OSError:
                continue
            if len(lines) < 2:
                continue
            row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
            try:
                samples.setdefault(acc, {})["reads_after_filter"] = int(float(row.get("reads_out", 0)))
            except ValueError:
                pass

    # 4: denoise — every published seqtab, attributed row by row (all of them, so
    # a sample present here but missing from `merge` reveals a DENOISE->MERGE loss).
    # Accumulate rather than assign: a sample appearing in two groups must not
    # overwrite itself, which would under-report that same loss.
    for f in sorted(glob.glob(os.path.join(seqtabs_dir, "*.seqtab.tsv"))):
        for acc, total in _seqtab_row_totals(f).items():
            d = samples.setdefault(acc, {})
            d["reads_after_denoise"] = d.get("reads_after_denoise", 0) + total

    # 5-7: merge / chimera / final per-sample read TSVs
    _read_sample_reads_tsv(os.path.join(seqtab_final_dir, "merge_sample_reads.tsv"),
                           samples, "reads_after_merge")
    _read_sample_reads_tsv(os.path.join(seqtab_final_dir, "chimera_sample_reads.tsv"),
                           samples, "reads_after_chimera")
    _read_sample_reads_tsv(os.path.join(seqtab_final_dir, "final_sample_reads.tsv"),
                           samples, "reads_final")

    cols = [col for _, _, col in STAGES]
    totals = dict.fromkeys(cols, 0)
    with open(out_path, "w") as out:
        out.write("sample\t" + "\t".join(cols) + "\tpct_primer\tpct_overall\n")
        for acc in sorted(samples):
            d = samples[acc]
            vals = [d.get(c, 0) for c in cols]
            for c, v in zip(cols, vals):
                totals[c] += v
            raw = vals[0] or 1
            out.write(
                f"{acc}\t" + "\t".join(str(v) for v in vals)
                + f"\t{100 * vals[1] / raw:.1f}\t{100 * vals[-1] / raw:.1f}\n"
            )
        raw = totals["reads_raw"] or 1
        out.write(
            "TOTAL\t" + "\t".join(str(totals[c]) for c in cols)
            + f"\t{100 * totals['reads_after_primer'] / raw:.1f}"
            + f"\t{100 * totals['reads_final'] / raw:.1f}\n"
        )

    # Attrition guard: flag any stage that loses more than half the samples that
    # entered it. Expected attrition (filter reads, prevalence pruning) stays
    # well under this; a silent 75% DENOISE->MERGE loss would trip it. We warn
    # loudly and record it in the provenance — never fail the run, so a sparse
    # but real dataset still produces (and deploys) a result.
    LOSS_FRACTION = 0.5
    stage_samples = {sid: sum(1 for acc in samples if samples[acc].get(col, 0) > 0)
                     for sid, _, col in STAGES}
    warnings = []
    prev = None
    for sid, lab, _ in STAGES:
        n = stage_samples[sid]
        if prev and prev[1] > 0 and n < prev[1] * (1 - LOSS_FRACTION):
            lost = prev[1] - n
            msg = (f"{lost}/{prev[1]} samples ({100 * lost / prev[1]:.0f}%) lost "
                   f"between '{prev[0]}' and '{sid}'")
            warnings.append({"stage": sid, "message": msg})
            print(f"[WARN] read_tracking: {msg} — inspect the provenance funnel",
                  file=sys.stderr)
        prev = (sid, n)

    # Provenance JSON for the viz Provenance tab. It must land in viz/ (deployed
    # as data/provenance.json); fall back to the tsv's dir if viz/ isn't there.
    payload = {
        "stages": [{"id": sid, "label": lab} for sid, lab, _ in STAGES],
        "total": {sid: totals[col] for sid, _, col in STAGES},
        "sample_counts": stage_samples,
        "warnings": warnings,
        "samples": {
            acc: {sid: samples[acc].get(col, 0) for sid, _, col in STAGES}
            for acc in sorted(samples)
        },
    }
    targets = [seqtab_final_dir]
    if os.path.isdir(viz_dir):
        targets.append(viz_dir)
    for d in targets:
        try:
            with open(os.path.join(d, "provenance.json"), "w") as jf:
                json.dump(payload, jf)
        except OSError as e:
            print(f"[WARN] could not write provenance.json to {d}: {e}")

    print(f"[INFO] read_tracking: {len(samples)} samples, "
          f"{totals['reads_raw']:,} raw -> {totals['reads_after_primer']:,} primer "
          f"-> {totals['reads_after_filter']:,} filter -> "
          f"{totals['reads_after_denoise']:,} denoise -> "
          f"{totals['reads_after_merge']:,} merge -> {totals['reads_final']:,} final")


if __name__ == "__main__":
    main()
