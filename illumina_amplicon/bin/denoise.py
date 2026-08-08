#!/usr/bin/env python3
#
# dada2_denoise.py — Denoise, merge pairs, and build a sequence table (per plate)
#
# Python mirror of dada2_denoise.R. For each sample on a plate:
#   1. Dereplicate forward and reverse filtered reads via dada2py
#   2. Apply the DADA2 denoising algorithm using pre-learned error models
#   3. Merge denoised forward/reverse pairs via C-level NW alignment
#   4. Combine all merged samples into a sequence table (pandas DataFrame)
#   5. Per-plate chimera removal via C-level consensus bimera detection
#
# Usage:
#   dada2_denoise.py <plate_id> <errF.pkl> <errR.pkl> <min_overlap> <cpus>
#
# Inputs (discovered automatically):
#   *_R1.filt.fastq.gz   Forward filtered reads
#   *_R2.filt.fastq.gz   Reverse filtered reads
#
# Outputs:
#   <plate_id>.seqtab.pkl   Sequence table (pickle, long-format DataFrame)
#   <plate_id>.seqtab.tsv   Sequence table (tab-delimited, for inspection)

import sys
import os
import re
import glob
import pickle
import numpy as np
import pandas as pd

import papa2 as dada2py


# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
if len(sys.argv) < 6:
    print("Usage: dada2_denoise.py <plate_id> <errF.pkl> <errR.pkl> "
          "<min_overlap> <cpus>", file=sys.stderr)
    sys.exit(1)

plate_id    = sys.argv[1]
errF_path   = sys.argv[2]
errR_path   = sys.argv[3]
min_overlap = int(sys.argv[4])
cpus        = int(sys.argv[5])


def write_empty_seqtab(reason):
    """Emit a valid, empty seqtab and exit 0 instead of crashing.

    A hard exit here is poison: DENOISE runs under errorStrategy 'ignore', and a
    single ignored task prematurely closes the DENOISE output channel, so the
    downstream `.collect()` fires MERGE with only the seqtabs finished so far —
    silently dropping every later sample. Degenerate inputs (too few reads to
    merge, dada() failures) must therefore succeed with an empty table, which is
    then dropped cleanly at the filter step while still being counted in the
    provenance. The sample row (no ASV columns) is kept so it shows up with 0.
    """
    print(f"[WARN] Plate {plate_id}: {reason} — emitting empty seqtab", file=sys.stderr)
    # Explicit dtypes matter: an all-object empty frame (the pandas default)
    # upcasts `count` to object when concatenated in MERGE, which then poisons
    # the numeric matrices downstream (t-SNE fails with "Unsupported dtype
    # object"). Keep count int64 so the merge stays numeric.
    empty = pd.DataFrame({
        "sample": pd.Series([], dtype="object"),
        "sequence": pd.Series([], dtype="object"),
        "count": pd.Series([], dtype="int64"),
    })
    with open(f"{plate_id}.seqtab.pkl", "wb") as f:
        pickle.dump(empty, f)
    with open(f"{plate_id}.seqtab.tsv", "w") as f:
        f.write("sample\n" + f"{plate_id}\n")
    sys.exit(0)

# ---------------------------------------------------------------------------
# Load pre-learned error models
# ---------------------------------------------------------------------------
with open(errF_path, "rb") as f:
    errF = pickle.load(f)
with open(errR_path, "rb") as f:
    errR = pickle.load(f)

# LEARN_ERRORS emits a sentinel when it couldn't build a model for the plate.
# Without a usable error model the samples can't be denoised — drop them via an
# empty seqtab (visible in the provenance) rather than crashing.
if (isinstance(errF, dict) and errF.get("__learn_errors_failed__")) or \
   (isinstance(errR, dict) and errR.get("__learn_errors_failed__")):
    write_empty_seqtab("no usable error model (learn_errors sentinel)")

# ---------------------------------------------------------------------------
# Discover filtered FASTQ files
# ---------------------------------------------------------------------------
fwd_files = sorted(glob.glob("*_R1.filt.fastq.gz"))
rev_files = sorted(glob.glob("*_R2.filt.fastq.gz"))

# Remove empty/missing files (keep pairs aligned)
valid_pairs = []
for f, r in zip(fwd_files, rev_files):
    f_ok = os.path.exists(f) and os.path.getsize(f) > 100
    r_ok = os.path.exists(r) and os.path.getsize(r) > 100
    if f_ok and r_ok:
        valid_pairs.append((f, r))

fwd_files = [p[0] for p in valid_pairs]
rev_files = [p[1] for p in valid_pairs]

# Derive sample names by stripping the _R1.filt.fastq.gz suffix
sample_names = [re.sub(r"_R1\.filt\.fastq\.gz$", "", os.path.basename(f))
                for f in fwd_files]

if len(fwd_files) == 0:
    write_empty_seqtab("no valid filtered files found")

print(f"[INFO] Plate {plate_id} : denoising {len(fwd_files)} samples")


# ---------------------------------------------------------------------------
# Merge pairs and chimera removal use dada2's C-level implementations
# via the dada2py library (py/paired.py and py/chimera.py).
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Per-sample: dereplicate -> denoise -> merge pairs
# ---------------------------------------------------------------------------
all_merged = {}  # sample_name -> {merged_seq: abundance}

for i, sample_name in enumerate(sample_names):
    print(f"[INFO] Processing: {sample_name}")

    # Dereplicate + denoise. Any failure here (degenerate/tiny sample) must not
    # crash the task — see write_empty_seqtab. Skip the sample instead.
    try:
        derepF = dada2py.derep_fastq(fwd_files[i])
        derepR = dada2py.derep_fastq(rev_files[i])
        ddF = dada2py.dada(derepF, err=errF)
        ddR = dada2py.dada(derepR, err=errR)
    except Exception as e:
        print(f"[WARNING] {sample_name}: derep/dada() failed: {e}")
        continue

    # Check for denoised output
    if not ddF.get("cluster_seqs") or not ddR.get("cluster_seqs"):
        print(f"[WARNING] {sample_name}: no denoised sequences, skipping")
        continue

    # Merge pairs using C-level NW alignment (matches R's mergePairs)
    merged_list = dada2py.merge_pairs(ddF, derepF, ddR, derepR,
                                       min_overlap=min_overlap,
                                       trim_overhang=True, verbose=False)

    # Collect accepted merged reads into {sequence: abundance}
    accepted = {}
    for m in merged_list:
        if m["accept"]:
            accepted[m["sequence"]] = accepted.get(m["sequence"], 0) + m["abundance"]

    if accepted:
        all_merged[sample_name] = accepted

# Drop samples that produced no merged reads — empty table, never a crash.
if len(all_merged) == 0:
    write_empty_seqtab("no merged reads produced")

# ---------------------------------------------------------------------------
# Build sequence table
# ---------------------------------------------------------------------------
# Collect all data in long format
rows = []
for sample_name, seq_dict in all_merged.items():
    for seq, count in seq_dict.items():
        rows.append({"sample": sample_name, "sequence": seq, "count": int(count)})

dt = pd.DataFrame(rows)
# Aggregate any duplicates
dt = dt.groupby(["sample", "sequence"], as_index=False)["count"].sum()

n_samples = dt["sample"].nunique()
n_asvs = dt["sequence"].nunique()
total_reads = dt["count"].sum()

print(f"[INFO] Plate {plate_id} : raw table has {n_samples} samples, "
      f"{n_asvs} unique ASVs, {total_reads} total reads")

# ---------------------------------------------------------------------------
# Per-plate chimera removal using dada2's C-level consensus method
# ---------------------------------------------------------------------------
all_seqs = sorted(set(s for d in all_merged.values() for s in d))
sample_list = list(all_merged.keys())
seq_to_idx = {s: i for i, s in enumerate(all_seqs)}

mat = np.zeros((len(sample_list), len(all_seqs)), dtype=np.int32)
for si, sam in enumerate(sample_list):
    for seq, ab in all_merged[sam].items():
        mat[si, seq_to_idx[seq]] = ab

seqtab = {"table": mat, "seqs": all_seqs}
chim_result = dada2py.remove_bimera_denovo(seqtab, verbose=True)

is_chimera = np.array(chim_result["is_chimera"])
chimeric_seqs = set(s for s, c in zip(all_seqs, is_chimera) if c)
dt_nochim = dt[~dt["sequence"].isin(chimeric_seqs)].copy()

n_chimeras = int(is_chimera.sum())
n_asvs_after = dt_nochim["sequence"].nunique()
reads_after = dt_nochim["count"].sum()
pct = round(reads_after / max(total_reads, 1) * 100, 1)

print(f"[INFO] Plate {plate_id} : removed {n_chimeras} chimeras, "
      f"retained {pct} % of reads ( {n_asvs_after} ASVs)")

# ---------------------------------------------------------------------------
# Save in both formats: pickle for downstream Python steps, TSV for inspection
# ---------------------------------------------------------------------------
with open(f"{plate_id}.seqtab.pkl", "wb") as f:
    pickle.dump(dt_nochim, f)

# Wide TSV for inspection
wide = dt_nochim.pivot_table(index="sample", columns="sequence",
                              values="count", fill_value=0, aggfunc="sum")
wide.to_csv(f"{plate_id}.seqtab.tsv", sep="\t")
