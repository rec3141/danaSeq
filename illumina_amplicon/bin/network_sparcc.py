#!/usr/bin/env python3
"""
network_sparcc.py - Run SparCC-style correlation analysis using CLR approach.

Usage: network_sparcc.py <renorm_merged.pkl> <min_prevalence> <cpus>

Outputs:
  sparcc_correlations.pkl  - Melted edge list DataFrame (node1, node2, correlation, weight, color)
  sparcc_stats.tsv         - Summary statistics
"""

import sys
import pickle
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


def log_info(msg):
    print(f"[INFO] {msg}", flush=True)


def log_error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr, flush=True)


def clr_transform(count_matrix, pseudocount=0.5):
    """Centered log-ratio transformation.

    Args:
        count_matrix: numpy array (samples x ASVs)
        pseudocount: value to add before log transformation

    Returns:
        CLR-transformed matrix (same shape)
    """
    # Add pseudocount
    mat = count_matrix + pseudocount

    # Geometric mean per sample (row)
    log_mat = np.log(mat)
    geo_mean_log = log_mat.mean(axis=1, keepdims=True)

    # CLR = log(x) - mean(log(x))
    clr = log_mat - geo_mean_log

    return clr


def compute_clr_correlation(clr_matrix):
    """Compute Pearson correlation on CLR-transformed data (between ASVs = columns)."""
    n_asvs = clr_matrix.shape[1]
    cor_matrix = np.corrcoef(clr_matrix.T)
    return cor_matrix


def main():
    if len(sys.argv) < 4:
        log_error("Usage: network_sparcc.py <renorm_merged.pkl> <min_prevalence> <cpus>")
        sys.exit(1)

    renorm_pkl = sys.argv[1]
    min_prevalence = float(sys.argv[2])
    n_cpus = int(sys.argv[3])

    log_info(f"Loading renormalized data from {renorm_pkl}")
    with open(renorm_pkl, "rb") as f:
        renorm = pickle.load(f)

    log_info(f"Loaded DataFrame with {len(renorm)} rows")

    # Expect long-format: sample, sequence, count
    required_cols = {"sample", "sequence", "count"}
    if not required_cols.issubset(set(renorm.columns)):
        col_map = {}
        for col in renorm.columns:
            cl = col.lower()
            if cl in ("sample", "sample_id", "sampleid"):
                col_map[col] = "sample"
            elif cl in ("sequence", "seq", "asv"):
                col_map[col] = "sequence"
            elif cl in ("count", "abundance", "reads"):
                col_map[col] = "count"
        if set(col_map.values()) >= required_cols:
            renorm = renorm.rename(columns=col_map)
            log_info(f"Renamed columns: {col_map}")
        else:
            log_error(f"Expected columns {required_cols}, got {set(renorm.columns)}")
            sys.exit(1)

    # Pivot to wide count matrix (samples x sequences)
    log_info("Pivoting to wide format")
    wide = renorm.pivot_table(index="sample", columns="sequence", values="count", fill_value=0)

    n_samples = wide.shape[0]
    n_asvs_total = wide.shape[1]
    log_info(f"Count matrix: {n_samples} samples x {n_asvs_total} ASVs")

    # Filter ASVs by prevalence (fraction of samples where ASV is present)
    prevalence = (wide > 0).sum(axis=0) / n_samples
    keep_mask = prevalence >= min_prevalence
    n_kept = keep_mask.sum()
    log_info(f"Prevalence filter (>= {min_prevalence}): keeping {n_kept}/{n_asvs_total} ASVs")

    if n_kept < 2:
        log_error("Fewer than 2 ASVs pass the prevalence filter. Cannot compute correlations.")
        # Save empty outputs
        empty_df = pd.DataFrame(columns=["node1", "node2", "correlation", "weight", "color"])
        with open("sparcc_correlations.pkl", "wb") as f:
            pickle.dump(empty_df, f)
        with open("sparcc_stats.tsv", "w") as f:
            f.write("metric\tvalue\n")
            f.write(f"n_samples\t{n_samples}\n")
            f.write(f"n_asvs_total\t{n_asvs_total}\n")
            f.write(f"n_asvs_filtered\t{n_kept}\n")
            f.write("n_edges\t0\n")
        log_info("Saved empty outputs (insufficient ASVs)")
        sys.exit(0)

    wide_filtered = wide.loc[:, keep_mask]
    asv_labels = list(wide_filtered.columns)

    # CLR transform
    log_info("Computing CLR transformation")
    clr_matrix = clr_transform(wide_filtered.values)

    # Compute Pearson correlation on CLR-transformed data
    log_info("Computing CLR-based correlations (SparCC approximation)")
    cor_matrix = compute_clr_correlation(clr_matrix)

    # Create correlation DataFrame
    cor_df = pd.DataFrame(cor_matrix, index=asv_labels, columns=asv_labels)

    # Melt to edge list (upper triangle only, no self-loops)
    log_info("Melting correlation matrix to edge list")
    edges = []
    n = len(asv_labels)
    for i in range(n):
        for j in range(i + 1, n):
            corr_val = cor_matrix[i, j]
            if np.isnan(corr_val):
                continue
            if abs(corr_val) > 0.1:
                edges.append({
                    "node1": asv_labels[i],
                    "node2": asv_labels[j],
                    "correlation": corr_val,
                })

    edge_df = pd.DataFrame(edges)
    log_info(f"Edges with |correlation| > 0.1: {len(edge_df)}")

    if len(edge_df) > 0:
        # Compute weight (absolute correlation) and color (positive=blue, negative=red)
        edge_df["weight"] = edge_df["correlation"].abs()
        edge_df["color"] = edge_df["correlation"].apply(
            lambda x: "positive" if x > 0 else "negative"
        )

        # Sort by absolute correlation descending
        edge_df = edge_df.sort_values("weight", ascending=False).reset_index(drop=True)
    else:
        edge_df = pd.DataFrame(columns=["node1", "node2", "correlation", "weight", "color"])

    # Save outputs
    log_info("Saving sparcc_correlations.pkl")
    with open("sparcc_correlations.pkl", "wb") as f:
        pickle.dump(edge_df, f)

    # Summary statistics
    n_positive = (edge_df["correlation"] > 0).sum() if len(edge_df) > 0 else 0
    n_negative = (edge_df["correlation"] < 0).sum() if len(edge_df) > 0 else 0
    mean_abs_corr = edge_df["weight"].mean() if len(edge_df) > 0 else 0.0
    max_corr = edge_df["correlation"].max() if len(edge_df) > 0 else 0.0
    min_corr = edge_df["correlation"].min() if len(edge_df) > 0 else 0.0

    log_info("Saving sparcc_stats.tsv")
    with open("sparcc_stats.tsv", "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"n_samples\t{n_samples}\n")
        f.write(f"n_asvs_total\t{n_asvs_total}\n")
        f.write(f"n_asvs_filtered\t{n_kept}\n")
        f.write(f"min_prevalence\t{min_prevalence}\n")
        f.write(f"n_edges\t{len(edge_df)}\n")
        f.write(f"n_positive\t{n_positive}\n")
        f.write(f"n_negative\t{n_negative}\n")
        f.write(f"mean_abs_correlation\t{mean_abs_corr:.6f}\n")
        f.write(f"max_correlation\t{max_corr:.6f}\n")
        f.write(f"min_correlation\t{min_corr:.6f}\n")

    log_info("network_sparcc.py complete")


if __name__ == "__main__":
    main()
