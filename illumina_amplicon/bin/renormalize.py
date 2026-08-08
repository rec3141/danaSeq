#!/usr/bin/env python3
#
# renormalize.py — Group ASVs by taxonomy and normalize within groups
#
# Assigns each ASV to a biological group (prokaryote, eukaryote,
# chloroplast, mitochondria, unknown) based on its taxonomy, then
# computes within-group proportions per sample.
#
# Usage:
#   renormalize.py <seqtab.pkl> <taxonomy.pkl> <tax_db_name>
#
# Outputs:
#   renorm_by_group.pkl     dict of DataFrames per group
#   renorm_merged.pkl       single DataFrame with group + proportion columns
#   renorm_table_list.pkl   sequence-to-group mapping
#   renorm_stats.tsv        per-group summary statistics

import sys
import os
import pickle
import time
import pandas as pd

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
if len(sys.argv) < 4:
    print("Usage: renormalize.py <seqtab.pkl> <taxonomy.pkl> <tax_db_name>",
          file=sys.stderr)
    sys.exit(1)

seqtab_path = sys.argv[1]
tax_path    = sys.argv[2]
db_name     = sys.argv[3]

t0 = time.time()

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
print(f"[INFO] Loading count data from {seqtab_path}")
with open(seqtab_path, "rb") as fh:
    seqtab = pickle.load(fh)

if not isinstance(seqtab, pd.DataFrame):
    print("[ERROR] seqtab.pkl must contain a pandas DataFrame", file=sys.stderr)
    sys.exit(1)

for col in ("sample", "sequence", "count"):
    if col not in seqtab.columns:
        print(f"[ERROR] Missing required column: {col}", file=sys.stderr)
        sys.exit(1)

print(f"[INFO] Count table: {len(seqtab)} rows, "
      f"{seqtab['sample'].nunique()} samples, "
      f"{seqtab['sequence'].nunique()} ASVs")

print(f"[INFO] Loading taxonomy from {tax_path}")
with open(tax_path, "rb") as fh:
    tax_data = pickle.load(fh)

if isinstance(tax_data, dict) and "tax" in tax_data:
    tax_df = tax_data["tax"]
elif isinstance(tax_data, pd.DataFrame):
    tax_df = tax_data
else:
    print("[ERROR] taxonomy.pkl must be a dict with 'tax' key or a DataFrame",
          file=sys.stderr)
    sys.exit(1)

print(f"[INFO] Taxonomy table: {len(tax_df)} sequences, "
      f"columns: {list(tax_df.columns)}")

# ---------------------------------------------------------------------------
# Map taxonomy columns (flexible naming)
# ---------------------------------------------------------------------------

def _find_col(df, candidates):
    """Find first matching column name (case-insensitive)."""
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None

col_domain = _find_col(tax_df, ["Kingdom", "Domain", "Level1"])
col_order  = _find_col(tax_df, ["Order", "Level4"])
col_family = _find_col(tax_df, ["Family", "Level5"])

print(f"[INFO] Taxonomy columns — Domain: {col_domain}, "
      f"Order: {col_order}, Family: {col_family}")

# ---------------------------------------------------------------------------
# Assign groups
# ---------------------------------------------------------------------------

def assign_group(row):
    """Classify a sequence into a biological group based on taxonomy."""
    domain = str(row.get(col_domain, "")).strip() if col_domain else ""
    order  = str(row.get(col_order, "")).strip() if col_order else ""
    family = str(row.get(col_family, "")).strip() if col_family else ""

    # Check organellar origins first (these override domain)
    if order == "Chloroplast":
        return "chloroplast"
    if family == "Mitochondria":
        return "mitochondria"
    if domain == "Eukaryota":
        return "eukaryote"
    if domain in ("Bacteria", "Archaea"):
        return "prokaryote"
    return "unknown"


# Build sequence -> group mapping
seq_groups = {}
for seq in tax_df.index:
    row = tax_df.loc[seq]
    seq_groups[seq] = assign_group(row)

# Sequences in count table but missing from taxonomy
all_seqs = seqtab["sequence"].unique()
n_missing = 0
for seq in all_seqs:
    if seq not in seq_groups:
        seq_groups[seq] = "unknown"
        n_missing += 1

if n_missing > 0:
    print(f"[INFO] {n_missing} sequences not in taxonomy — assigned to 'unknown'")

group_counts = pd.Series(seq_groups).value_counts()
print("[INFO] Group assignments:")
for grp, cnt in group_counts.items():
    print(f"[INFO]   {grp}: {cnt} ASVs")

# ---------------------------------------------------------------------------
# Add group column to count data
# ---------------------------------------------------------------------------
seqtab = seqtab.copy()
seqtab["group"] = seqtab["sequence"].map(seq_groups)

# ---------------------------------------------------------------------------
# Normalize: proportion = count / sum(count for group in sample)
# ---------------------------------------------------------------------------
print("[INFO] Computing within-group proportions ...")

group_sample_totals = seqtab.groupby(["group", "sample"])["count"].transform("sum")
seqtab["proportion"] = seqtab["count"] / group_sample_totals
# Handle division by zero (empty groups)
seqtab["proportion"] = seqtab["proportion"].fillna(0.0)

# ---------------------------------------------------------------------------
# Split by group
# ---------------------------------------------------------------------------
by_group = {}
for grp, sub in seqtab.groupby("group"):
    by_group[grp] = sub.reset_index(drop=True)
    print(f"[INFO] Group '{grp}': {len(sub)} rows, "
          f"{sub['sample'].nunique()} samples, "
          f"{sub['sequence'].nunique()} ASVs")

# ---------------------------------------------------------------------------
# Build summary statistics
# ---------------------------------------------------------------------------
stats_rows = []
for grp in sorted(by_group.keys()):
    sub = by_group[grp]
    stats_rows.append({
        "group": grp,
        "n_asvs": sub["sequence"].nunique(),
        "n_samples": sub["sample"].nunique(),
        "total_reads": int(sub["count"].sum()),
        "mean_reads_per_sample": round(
            sub.groupby("sample")["count"].sum().mean(), 1),
        "median_reads_per_sample": round(
            sub.groupby("sample")["count"].sum().median(), 1),
    })

stats_df = pd.DataFrame(stats_rows)

# ---------------------------------------------------------------------------
# Sequence-to-group mapping table
# ---------------------------------------------------------------------------
table_list = pd.DataFrame({
    "sequence": list(seq_groups.keys()),
    "group": list(seq_groups.values()),
})

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
with open("renorm_by_group.pkl", "wb") as fh:
    pickle.dump(by_group, fh, protocol=4)
print("[INFO] Saved renorm_by_group.pkl")

with open("renorm_merged.pkl", "wb") as fh:
    pickle.dump(seqtab, fh, protocol=4)
print("[INFO] Saved renorm_merged.pkl")

with open("renorm_table_list.pkl", "wb") as fh:
    pickle.dump(table_list, fh, protocol=4)
print("[INFO] Saved renorm_table_list.pkl")

stats_df.to_csv("renorm_stats.tsv", sep="\t", index=False)
print("[INFO] Saved renorm_stats.tsv")

elapsed = time.time() - t0
print(f"[INFO] Renormalization complete in {elapsed:.1f}s")
print(f"[INFO] Total reads: {int(seqtab['count'].sum()):,}")
print(f"[INFO] Groups: {', '.join(sorted(by_group.keys()))}")
