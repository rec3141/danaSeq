#!/usr/bin/env python3
#
# filter_seqtab.py — Quality-control filtering on long-format sequence data
#
# Python mirror of filter_seqtab.R. Applies a cascade of filters to remove
# noise and low-confidence observations. Operates entirely in long format
# (pandas DataFrame) — never builds the full dense matrix, so memory use
# is proportional to non-zero entries.
#
# Filters (in order):
#   1. Length      — ASVs shorter than a threshold (can't assign taxonomy)
#   2. Prevalence  — ASVs in fewer than N samples (likely artefacts)
#   3. Abundance   — ASVs with fewer than N total reads (singletons/doubletons)
#   4. Depth       — samples with too few total reads (insufficient coverage)
#
# Usage:
#   filter_seqtab.py <seqtab_long.pkl> <min_length> <min_samples> <min_seqs> <min_reads>
#
# Outputs:
#   seqtab_final.pkl          Filtered long-format DataFrame
#   seqtab_final_wide.csv     Filtered dense matrix (CSV for compatibility)
#   seqtab_orphans.pkl        ASVs removed by prevalence filter
#   seqtab_small.pkl          Samples removed by depth filter
#   filter_stats.tsv          Per-step filtering statistics
#   sequence_summaries.pdf    Diagnostic rank-abundance plots

import sys
import os
import pickle
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
if len(sys.argv) < 6:
    print("Usage: filter_seqtab.py <seqtab_long.pkl> "
          "<min_length> <min_samples> <min_seqs> <min_reads>",
          file=sys.stderr)
    sys.exit(1)

input_path  = sys.argv[1]
min_length  = int(sys.argv[2])
min_samples = int(sys.argv[3])
min_seqs    = int(sys.argv[4])
min_reads   = int(sys.argv[5])

# ---------------------------------------------------------------------------
# Load input (long-format DataFrame)
# ---------------------------------------------------------------------------
with open(input_path, "rb") as f:
    dt = pickle.load(f)

if not isinstance(dt, pd.DataFrame):
    print("[ERROR] Expected long-format DataFrame with columns: "
          "sample, sequence, count", file=sys.stderr)
    sys.exit(1)

total_input  = int(dt["count"].sum())
n_input_asvs = dt["sequence"].nunique()

print(f"[INFO] Input: {dt['sample'].nunique()} samples, "
      f"{n_input_asvs} ASVs, {total_input} reads")

# ---------------------------------------------------------------------------
# 1. Remove short sequences
#    ASVs shorter than ~50 bp cannot be reliably assigned taxonomy.
# ---------------------------------------------------------------------------
seq_lengths = dt.groupby("sequence").agg(
    seq_len=("sequence", lambda x: len(x.iloc[0]))
).reset_index()
short_seqs = set(seq_lengths.loc[seq_lengths["seq_len"] < min_length, "sequence"])
n_short = len(short_seqs)

if n_short > 0:
    print(f"[INFO] Removing {n_short} ASVs shorter than {min_length} bp")
    dt = dt[~dt["sequence"].isin(short_seqs)].copy()

# ---------------------------------------------------------------------------
# 2. Remove low-prevalence ASVs (orphans)
#    ASVs appearing in very few samples are likely artefacts.
# ---------------------------------------------------------------------------
seq_prevalence = dt.groupby("sequence")["sample"].nunique().reset_index()
seq_prevalence.columns = ["sequence", "n_samples"]
orphan_seqs = set(seq_prevalence.loc[seq_prevalence["n_samples"] < min_samples,
                                      "sequence"])
n_orphans = len(orphan_seqs)

# Save orphans for review before removing
dt_orphans = dt[dt["sequence"].isin(orphan_seqs)].copy()
dt = dt[~dt["sequence"].isin(orphan_seqs)].copy()

print(f"[INFO] Removed {n_orphans} orphan ASVs (present in < "
      f"{min_samples} samples)")

# ---------------------------------------------------------------------------
# 3. Remove low-abundance ASVs
#    ASVs with very few total reads are unreliable.
# ---------------------------------------------------------------------------
seq_abundance = dt.groupby("sequence")["count"].sum().reset_index()
seq_abundance.columns = ["sequence", "total"]
rare_seqs = set(seq_abundance.loc[seq_abundance["total"] < min_seqs, "sequence"])
n_rare = len(rare_seqs)

dt = dt[~dt["sequence"].isin(rare_seqs)].copy()

print(f"[INFO] Removed {n_rare} rare ASVs (< {min_seqs} total reads)")

# ---------------------------------------------------------------------------
# 4. Remove shallow samples
#    Samples with too few reads lack sufficient sequencing depth.
# ---------------------------------------------------------------------------
sample_depth = dt.groupby("sample")["count"].sum().reset_index()
sample_depth.columns = ["sample", "total"]
small_samples = set(sample_depth.loc[sample_depth["total"] < min_reads, "sample"])
n_small = len(small_samples)

# Save small samples for review before removing
dt_small = dt[dt["sample"].isin(small_samples)].copy()
dt = dt[~dt["sample"].isin(small_samples)].copy()

print(f"[INFO] Removed {n_small} samples (< {min_reads} reads)")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
n_final_samples = dt["sample"].nunique()
n_final_asvs    = dt["sequence"].nunique()
n_final_reads   = int(dt["count"].sum())
pct_retained    = round(n_final_reads / max(total_input, 1) * 100, 1)

print(f"[INFO] Final: {n_final_samples} samples, {n_final_asvs} ASVs, "
      f"{n_final_reads} reads ( {pct_retained} % of input)")

# ---------------------------------------------------------------------------
# Diagnostic plots — four-panel rank-abundance overview
# ---------------------------------------------------------------------------
with PdfPages("sequence_summaries.pdf") as pdf:
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    if len(dt) > 0:
        reads_per_asv = dt.groupby("sequence")["count"].sum().sort_values()
        reads_per_sample = dt.groupby("sample")["count"].sum().sort_values()
        samps_per_asv = dt.groupby("sequence")["sample"].nunique().sort_values()
        asvs_per_sample = dt.groupby("sample")["sequence"].nunique().sort_values()

        ax = axes[0, 0]
        ax.plot(range(len(reads_per_asv)),
                np.log10(reads_per_asv.values.astype(float)),
                'o', markersize=1, color='steelblue')
        ax.set_title("Reads per ASV")
        ax.set_ylabel("log10(reads)")
        ax.set_xlabel("ASV rank")

        ax = axes[0, 1]
        ax.plot(range(len(reads_per_sample)),
                np.log10(reads_per_sample.values.astype(float)),
                'o', markersize=1, color='steelblue')
        ax.set_title("Reads per Sample")
        ax.set_ylabel("log10(reads)")
        ax.set_xlabel("Sample rank")

        ax = axes[1, 0]
        ax.plot(range(len(samps_per_asv)),
                np.log10(samps_per_asv.values.astype(float)),
                'o', markersize=1, color='darkgreen')
        ax.set_title("Samples per ASV")
        ax.set_ylabel("log10(samples)")
        ax.set_xlabel("ASV rank")

        ax = axes[1, 1]
        ax.plot(range(len(asvs_per_sample)),
                np.log10(asvs_per_sample.values.astype(float)),
                'o', markersize=1, color='darkgreen')
        ax.set_title("ASVs per Sample")
        ax.set_ylabel("log10(ASVs)")
        ax.set_xlabel("Sample rank")
    else:
        axes[0, 0].text(0.5, 0.5, "No data remaining after filtering",
                         ha='center', va='center', fontsize=14,
                         transform=axes[0, 0].transAxes)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close('all')

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
with open("seqtab_final.pkl", "wb") as f:
    pickle.dump(dt, f)
with open("seqtab_orphans.pkl", "wb") as f:
    pickle.dump(dt_orphans, f)
with open("seqtab_small.pkl", "wb") as f:
    pickle.dump(dt_small, f)

# Per-sample read counts in the final table — the last provenance stage.
if len(dt) > 0:
    (dt.groupby("sample")["count"].sum().rename("reads")
       .reset_index().to_csv("final_sample_reads.tsv", sep="\t", index=False))
else:
    with open("final_sample_reads.tsv", "w") as f:
        f.write("sample\treads\n")

# Wide matrix for compatibility
print("[INFO] Casting final table to wide matrix...")
if len(dt) > 0:
    wide = dt.pivot_table(index="sample", columns="sequence",
                           values="count", fill_value=0, aggfunc="sum")
    wide.to_csv("seqtab_final_wide.csv")
    with open("seqtab_final_wide.pkl", "wb") as f:
        pickle.dump(wide, f)
else:
    pd.DataFrame().to_csv("seqtab_final_wide.csv")
    with open("seqtab_final_wide.pkl", "wb") as f:
        pickle.dump(pd.DataFrame(), f)

# Stats
stats = pd.DataFrame([
    {"step": "length",     "asvs_removed": n_short,   "samples_removed": None,
     "remaining_samples": None, "remaining_asvs": None, "pct_reads_retained": None},
    {"step": "prevalence", "asvs_removed": n_orphans,  "samples_removed": None,
     "remaining_samples": None, "remaining_asvs": None, "pct_reads_retained": None},
    {"step": "abundance",  "asvs_removed": n_rare,     "samples_removed": None,
     "remaining_samples": None, "remaining_asvs": None, "pct_reads_retained": None},
    {"step": "depth",      "asvs_removed": None,       "samples_removed": n_small,
     "remaining_samples": None, "remaining_asvs": None, "pct_reads_retained": None},
    {"step": "final",      "asvs_removed": None,       "samples_removed": None,
     "remaining_samples": n_final_samples, "remaining_asvs": n_final_asvs,
     "pct_reads_retained": pct_retained},
])
stats.to_csv("filter_stats.tsv", sep="\t", index=False)
