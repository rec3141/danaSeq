#!/usr/bin/env python3
#
# merge_seqtabs.py — Merge per-plate sequence tables using long-format approach
#
# Python mirror of merge_seqtabs.R. Discovers per-plate pickle files, melts
# each to long format (sample, sequence, count), concatenates, and aggregates
# any duplicate (sample, sequence) pairs by summing counts.
#
# Usage:
#   merge_seqtabs.py
#
# Inputs (discovered automatically):
#   *.seqtab.pkl   One or more per-plate sequence tables in the working directory
#
# Outputs:
#   seqtab_merged.pkl   Long-format DataFrame (sample, sequence, count)
#   merge_stats.tsv     Summary statistics

import sys
import os
import glob
import pickle
import pandas as pd


# ---------------------------------------------------------------------------
# Discover per-plate sequence tables
# ---------------------------------------------------------------------------
pkl_files = sorted(glob.glob("*.seqtab.pkl"))

if len(pkl_files) == 0:
    print("[ERROR] No .seqtab.pkl files found in working directory",
          file=sys.stderr)
    sys.exit(1)

print(f"[INFO] Found {len(pkl_files)} sequence tables to merge")

# ---------------------------------------------------------------------------
# Load each plate table and stack
#
# Each seqtab.pkl is a DataFrame in long format (sample, sequence, count)
# as produced by dada2_denoise.py. We concatenate all of them.
# ---------------------------------------------------------------------------
long_list = []

for pkl_file in pkl_files:
    with open(pkl_file, "rb") as f:
        st = pickle.load(f)

    # Handle both wide (samples x sequences matrix) and long format
    if isinstance(st, pd.DataFrame):
        if "sample" in st.columns and "sequence" in st.columns:
            # Already long format
            plate_long = st[["sample", "sequence", "count"]].copy()
        else:
            # Wide format — melt to long
            st_reset = st.reset_index()
            id_col = st_reset.columns[0]
            plate_long = st_reset.melt(id_vars=[id_col],
                                        var_name="sequence",
                                        value_name="count")
            plate_long = plate_long.rename(columns={id_col: "sample"})
            plate_long = plate_long[plate_long["count"] > 0]
    else:
        print(f"[ERROR] Unexpected type in {pkl_file}: {type(st)}",
              file=sys.stderr)
        sys.exit(1)

    n_samples = plate_long["sample"].nunique()
    n_asvs = plate_long["sequence"].nunique()
    n_reads = plate_long["count"].sum()
    print(f"  {os.path.basename(pkl_file)} -> {n_samples} samples, "
          f"{n_asvs} ASVs, {n_reads} reads")

    long_list.append(plate_long)

print("[INFO] Concatenating long tables...")
dt = pd.concat(long_list, ignore_index=True)
del long_list

# Defensive: keep counts numeric. A stray object-dtype input (e.g. an empty
# per-plate table) can upcast the whole column to object under concat, which
# breaks every numeric consumer downstream (t-SNE, renormalize, viz).
dt["count"] = pd.to_numeric(dt["count"], errors="coerce").fillna(0).astype("int64")

# ---------------------------------------------------------------------------
# Aggregate: sum counts for (sample, sequence) pairs that appear on multiple
# plates (e.g., same sample re-sequenced). In most cases there are no
# duplicates and this is a no-op, but it's necessary for correctness.
# ---------------------------------------------------------------------------
n_before = len(dt)
dt = dt.groupby(["sample", "sequence"], as_index=False)["count"].sum()
n_dupes = n_before - len(dt)
if n_dupes > 0:
    print(f"[INFO] Aggregated {n_dupes} duplicate (sample, sequence) entries")

n_samples = dt["sample"].nunique()
n_asvs    = dt["sequence"].nunique()
n_reads   = int(dt["count"].sum())

print(f"[INFO] Merged total: {n_samples} samples, {n_asvs} unique ASVs, "
      f"{n_reads} reads, {len(dt)} non-zero entries")

# ---------------------------------------------------------------------------
# Save as long-format DataFrame (no dense matrix)
# ---------------------------------------------------------------------------
with open("seqtab_merged.pkl", "wb") as f:
    pickle.dump(dt, f)

stats = pd.DataFrame([{
    "step": "merged",
    "samples": n_samples,
    "sequences": n_asvs,
    "total_reads": n_reads,
    "nonzero_cells": len(dt),
}])
stats.to_csv("merge_stats.tsv", sep="\t", index=False)

# Per-sample read counts after merge — feeds the provenance tracker. Comparing
# this against the denoise stage (published per-sample seqtabs) is what exposes
# any sample lost between DENOISE and MERGE.
(dt.groupby("sample")["count"].sum().rename("reads")
   .reset_index().to_csv("merge_sample_reads.tsv", sep="\t", index=False))
