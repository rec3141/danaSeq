#!/usr/bin/env python3
#
# remove_chimeras.py — Consensus chimera removal using dada2's C library
#
# Python mirror of remove_chimeras.R. Converts long-format count data to a
# dense matrix and calls dada2's C-level table_bimera() for consensus
# chimera detection (identical algorithm to R's removeBimeraDenovo).
#
# Usage:
#   remove_chimeras.py <seqtab.pkl> <cpus>
#
# Outputs:
#   seqtab_nochim.pkl    Chimera-free sequence table (long-format DataFrame)
#   chimera_stats.tsv    Summary statistics

import sys
import os
import pickle
import pandas as pd
import numpy as np

import papa2 as dada2py


# ===========================================================================
# Main script
# ===========================================================================

if len(sys.argv) < 3:
    print("Usage: remove_chimeras.py <seqtab_long.pkl> <cpus>",
          file=sys.stderr)
    sys.exit(1)

input_path = sys.argv[1]
cpus       = int(sys.argv[2])

# ---------------------------------------------------------------------------
# Load input (long-format DataFrame from merge_seqtabs.py)
# ---------------------------------------------------------------------------
with open(input_path, "rb") as f:
    dt = pickle.load(f)

if not isinstance(dt, pd.DataFrame):
    print("[ERROR] Expected long-format DataFrame with columns: "
          "sample, sequence, count", file=sys.stderr)
    sys.exit(1)

n_input_asvs  = dt["sequence"].nunique()
n_input_reads = int(dt["count"].sum())
n_samples     = dt["sample"].nunique()

print(f"[INFO] Input: {n_samples} samples, {n_input_asvs} ASVs, "
      f"{n_input_reads} reads, {len(dt)} non-zero entries")

# ---------------------------------------------------------------------------
# Pre-filter: remove obvious noise before the expensive chimera check
# ---------------------------------------------------------------------------
seq_stats = dt.groupby("sequence").agg(
    total=("count", "sum"),
    seq_len=("sequence", lambda x: len(x.iloc[0]))
).reset_index()

pre_singleton = set(seq_stats.loc[seq_stats["total"] <= 1, "sequence"])
pre_short     = set(seq_stats.loc[seq_stats["seq_len"] < 20, "sequence"])
pre_remove    = pre_singleton | pre_short
n_pre_removed = len(pre_remove)

if n_pre_removed > 0:
    n_short_only = len(pre_short - pre_singleton)
    print(f"[INFO] Pre-filter: removing {len(pre_singleton)} singletons and "
          f"{n_short_only} ultra-short ASVs")
    dt = dt[~dt["sequence"].isin(pre_remove)].copy()

print(f"[INFO] Chimera input: {dt['sequence'].nunique()} ASVs, "
      f"{len(dt)} non-zero entries")

# ---------------------------------------------------------------------------
# Run chimera detection using dada2's C-level consensus method
# ---------------------------------------------------------------------------
# Build dense matrix from long-format DataFrame
all_seqs = sorted(dt["sequence"].unique())
all_samples = sorted(dt["sample"].unique())
seq_to_idx = {s: i for i, s in enumerate(all_seqs)}
sam_to_idx = {s: i for i, s in enumerate(all_samples)}

mat = np.zeros((len(all_samples), len(all_seqs)), dtype=np.int32)
row_idx = dt["sample"].map(sam_to_idx).values
col_idx = dt["sequence"].map(seq_to_idx).values
mat[row_idx, col_idx] = dt["count"].values.astype(np.int32)

print(f"[INFO] Checking {len(all_seqs)} ASVs for chimeras "
      f"({len(all_samples)} samples)")

seqtab = {"table": mat, "seqs": all_seqs}
chim_result = dada2py.remove_bimera_denovo(seqtab, verbose=True)

is_chimera = np.array(chim_result["is_chimera"])
chimeric_seqs = set(s for s, c in zip(all_seqs, is_chimera) if c)

# Remove chimeric sequences
dt_clean = dt[~dt["sequence"].isin(chimeric_seqs)].copy()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
chimeras_removed = len(chimeric_seqs)
total_removed    = n_pre_removed + chimeras_removed
n_output_asvs    = dt_clean["sequence"].nunique()
n_output_reads   = int(dt_clean["count"].sum())
pct_retained     = round(n_output_reads / max(n_input_reads, 1) * 100, 1)

print(f"[INFO] Removed {chimeras_removed} chimeras + "
      f"{n_pre_removed} pre-filtered = {total_removed} total")
print(f"[INFO] Retained {pct_retained} % of reads ( "
      f"{n_output_asvs} ASVs across {dt_clean['sample'].nunique()} samples)")

# ---------------------------------------------------------------------------
# Save output — stays in long format
# ---------------------------------------------------------------------------
with open("seqtab_nochim.pkl", "wb") as f:
    pickle.dump(dt_clean, f)

stats = pd.DataFrame([{
    "step": "chimera_removal",
    "input_asvs": n_input_asvs,
    "pre_filtered": n_pre_removed,
    "chimeras": chimeras_removed,
    "output_asvs": n_output_asvs,
    "samples": dt_clean["sample"].nunique(),
    "total_reads": n_output_reads,
    "pct_retained": pct_retained,
}])
stats.to_csv("chimera_stats.tsv", sep="\t", index=False)

# Per-sample read counts after chimera removal — feeds the provenance tracker.
(dt_clean.groupby("sample")["count"].sum().rename("reads")
         .reset_index().to_csv("chimera_sample_reads.tsv", sep="\t", index=False))
