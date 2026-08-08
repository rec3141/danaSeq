#!/usr/bin/env python3
#
# load_metadata.py — Load sample metadata from MIMARKS-compliant TSV/CSV
#
# Usage:
#   load_metadata.py <seqtab.pkl> <metadata_file> [sample_id_column]

import sys
import os
import pickle
import pandas as pd

if len(sys.argv) < 3:
    print("Usage: load_metadata.py <seqtab.pkl> <metadata_file> [sample_id_column]",
          file=sys.stderr)
    sys.exit(1)

seqtab_path = sys.argv[1]
meta_path   = sys.argv[2]
id_col      = sys.argv[3] if len(sys.argv) > 3 else "sample_name"

# Load seqtab to get sample IDs
print(f"[INFO] Loading sequence table: {seqtab_path}")
dt = pickle.load(open(seqtab_path, "rb"))
pipeline_samples = sorted(dt["sample"].unique()) if isinstance(dt, pd.DataFrame) else []
print(f"[INFO] {len(pipeline_samples)} samples in pipeline")

# Load metadata
print(f"[INFO] Loading metadata: {meta_path}")
ext = os.path.splitext(meta_path)[1].lower()
if ext in (".tsv", ".txt"):
    meta = pd.read_csv(meta_path, sep="\t")
else:
    meta = pd.read_csv(meta_path)

print(f"[INFO] Metadata: {len(meta)} rows, {len(meta.columns)} columns")

# MIMARKS fields
mimarks = ["sample_name", "collection_date", "geo_loc_name", "lat_lon",
           "depth", "env_broad_scale", "env_local_scale", "env_medium"]
found = [f for f in mimarks if f in meta.columns]
if found:
    print(f"[INFO] MIMARKS fields found: {', '.join(found)}")

# Match samples
if id_col in meta.columns:
    meta = meta.set_index(id_col)
else:
    print(f"[WARNING] Column '{id_col}' not found, using first column")
    meta = meta.set_index(meta.columns[0])

matched = set(pipeline_samples) & set(meta.index)
in_pipe_only = set(pipeline_samples) - set(meta.index)
in_meta_only = set(meta.index) - set(pipeline_samples)

print(f"[INFO] Matched: {len(matched)}")
print(f"[INFO] In pipeline only: {len(in_pipe_only)}")
print(f"[INFO] In metadata only: {len(in_meta_only)}")

# Save
with open("metadata.pkl", "wb") as f:
    pickle.dump(meta, f)

stats = pd.DataFrame([{
    "matched": len(matched),
    "pipeline_only": len(in_pipe_only),
    "metadata_only": len(in_meta_only),
    "total_pipeline": len(pipeline_samples),
    "total_metadata": len(meta)
}])
stats.to_csv("match_stats.tsv", sep="\t", index=False)
