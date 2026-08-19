#!/usr/bin/env python3
"""
cluster_tsne.py - Compute t-SNE ordinations for samples and ASVs from a sequence table.

Usage: cluster_tsne.py <seqtab.pkl> <cpus> [primer_assignment.tsv]

Outputs:
  sample_bray_tsne.pkl  - DataFrame (label, tSNE1, tSNE2)
  seq_bray_tsne.pkl     - DataFrame (label, tSNE1, tSNE2)
  sample_bray_dist.pkl  - Bray-Curtis distance matrix (samples)
  seq_bray_dist.pkl     - Bray-Curtis distance matrix (ASVs)
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis, squareform, pdist
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


def log_info(msg):
    print(f"[INFO] {msg}", flush=True)


def log_error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr, flush=True)


def compute_bray_curtis_matrix(prop_matrix):
    """Compute pairwise Bray-Curtis distance matrix from a proportional matrix."""
    n = prop_matrix.shape[0]
    dist_vec = pdist(prop_matrix, metric="braycurtis")
    dist_mat = squareform(dist_vec)
    return dist_mat


def run_ordination(dist_matrix, labels, entity_name, n_cpus):
    """Run PCA then t-SNE on a distance matrix. Returns a DataFrame with label, tSNE1, tSNE2."""
    n = dist_matrix.shape[0]
    log_info(f"{entity_name}: {n} entities for ordination")

    # PCA on distance matrix
    n_pca_components = min(50, n - 1) if n > 1 else 1
    log_info(f"{entity_name}: Running PCA with {n_pca_components} components")
    pca = PCA(n_components=n_pca_components)
    pca_coords = pca.fit_transform(dist_matrix)
    log_info(f"{entity_name}: PCA explained variance ratio sum = {pca.explained_variance_ratio_.sum():.4f}")

    if n < 3:
        # Too few entities for t-SNE; use first two PCA components
        log_info(f"{entity_name}: Fewer than 3 entities, skipping t-SNE, using PCA coordinates")
        dim1 = pca_coords[:, 0] if pca_coords.shape[1] >= 1 else np.zeros(n)
        dim2 = pca_coords[:, 1] if pca_coords.shape[1] >= 2 else np.zeros(n)
        result = pd.DataFrame({
            "label": labels,
            "tSNE1": dim1,
            "tSNE2": dim2,
        })
        return result

    # Determine perplexity
    perplexity = 30.0
    if n <= 30:
        # Perplexity must be less than number of samples
        perplexity = max(1.0, float(n - 1) / 2.0)
        log_info(f"{entity_name}: Reducing perplexity to {perplexity} (n={n} < 30)")

    n_iter = 5000
    log_info(f"{entity_name}: Running t-SNE (perplexity={perplexity}, max_iter={n_iter})")
    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        max_iter=n_iter,
        random_state=42,
    )
    tsne_coords = tsne.fit_transform(pca_coords)
    log_info(f"{entity_name}: t-SNE complete")

    result = pd.DataFrame({
        "label": labels,
        "tSNE1": tsne_coords[:, 0],
        "tSNE2": tsne_coords[:, 1],
    })
    return result


def _blocks(dist_matrix, labels, groups):
    """Sets of samples that share nothing with the rest, largest first.

    Two samples amplified with primers for different genes cannot hold a
    sequence in common, so every distance between them is the maximum Bray-Curtis
    can return. An ordination given a matrix that is saturated everywhere except
    within blocks has no gradient to descend between them and arranges the whole
    run on a circle — which reads as a result, and is the absence of one.

    Assay assignment names the blocks when it is available. It is checked rather
    than trusted: two assays of the same gene do share sequences, and splitting
    those would discard the comparison the run was for.
    """
    if not groups:
        return None
    by_assay = {}
    for i, lab in enumerate(labels):
        by_assay.setdefault(groups.get(lab, ""), []).append(i)
    by_assay.pop("", None)
    if len(by_assay) < 2:
        return None

    keys = sorted(by_assay, key=lambda k: -len(by_assay[k]))
    for a, b in ((x, y) for i, x in enumerate(keys) for y in keys[i + 1:]):
        pairs = [dist_matrix[i][j] for i in by_assay[a] for j in by_assay[b]]
        # Bray-Curtis reaches 1 only when two samples share no sequence at all.
        if not pairs or min(pairs) < 0.999:
            return None
    return [(k, by_assay[k]) for k in keys]


def run_ordination_by_block(dist_matrix, labels, groups, entity_name, n_cpus):
    """Ordinate each block on its own, then set the blocks side by side.

    Within a block the layout means what it always meant. Between blocks it means
    nothing, and cannot: there is no shared sequence to measure a distance from.
    Laying them out separately at least stops the emptiness between them from
    destroying the structure inside them.
    """
    blocks = _blocks(dist_matrix, labels, groups)
    if not blocks:
        return run_ordination(dist_matrix, labels, entity_name, n_cpus)

    log_info(f"{entity_name}: {len(blocks)} assay(s) share no sequences with each "
             f"other; ordinating each on its own")
    frames, offset = [], 0.0
    for name, idx in blocks:
        sub = dist_matrix[np.ix_(idx, idx)]
        sub_labels = [labels[i] for i in idx]
        log_info(f"{entity_name}: {name} ({len(idx)} entities)")
        part = run_ordination(sub, sub_labels, f"{entity_name}/{name}", n_cpus)
        span = float(part["tSNE1"].max() - part["tSNE1"].min()) or 1.0
        part["tSNE1"] = part["tSNE1"] - part["tSNE1"].min() + offset
        offset += span * 1.35   # a gap wider than either block, so it reads as one
        frames.append(part)
    return pd.concat(frames, ignore_index=True)


def read_assays(path):
    """{sample: assay} from primer_assignment.tsv, or {} when there is none."""
    if not path or not os.path.exists(path):
        return {}
    try:
        table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    except Exception as exc:
        log_info(f"assay assignment unreadable ({exc}); ordinating the run as one")
        return {}
    if "sample" not in table.columns:
        return {}
    keys = [c for c in ("assay_set", "assay_gene", "assay_region",
                        "assay_primer_fwd", "assay_primer_rev") if c in table.columns]
    if not keys:
        return {}
    out = {}
    for _, row in table.iterrows():
        # Placed assays are named by where they sit; the rest by their primers.
        name = row.get("assay_set") or "|".join(row.get(k, "") for k in keys)
        if name.strip("|"):
            out[row["sample"]] = name
    return out

def main():
    if len(sys.argv) < 3:
        log_error("Usage: cluster_tsne.py <seqtab.pkl> <cpus>")
        sys.exit(1)

    seqtab_pkl = sys.argv[1]
    n_cpus = int(sys.argv[2])
    assay_path = sys.argv[3] if len(sys.argv) > 3 else None
    assays = read_assays(assay_path)
    if assays:
        log_info(f"Assay assignment for {len(assays)} samples")

    log_info(f"Loading sequence table from {seqtab_pkl}")
    with open(seqtab_pkl, "rb") as f:
        seqtab = pickle.load(f)

    log_info(f"Loaded DataFrame with {len(seqtab)} rows")

    # Expect long-format: sample, sequence, count
    required_cols = {"sample", "sequence", "count"}
    if not required_cols.issubset(set(seqtab.columns)):
        # Try common alternative column names
        col_map = {}
        for col in seqtab.columns:
            cl = col.lower()
            if cl in ("sample", "sample_id", "sampleid"):
                col_map[col] = "sample"
            elif cl in ("sequence", "seq", "asv"):
                col_map[col] = "sequence"
            elif cl in ("count", "abundance", "reads"):
                col_map[col] = "count"
        if set(col_map.values()) >= required_cols:
            seqtab = seqtab.rename(columns=col_map)
            log_info(f"Renamed columns: {col_map}")
        else:
            log_error(f"Expected columns {required_cols}, got {set(seqtab.columns)}")
            sys.exit(1)

    # Pivot to wide matrix (samples x sequences)
    log_info("Pivoting to wide format")
    wide = seqtab.pivot_table(index="sample", columns="sequence", values="count", fill_value=0)

    # Normalize rows to proportions, then 4th-root transform
    row_sums = wide.sum(axis=1)
    row_sums = row_sums.replace(0, 1)  # avoid division by zero
    prop_matrix = wide.div(row_sums, axis=0)
    transformed = prop_matrix.apply(lambda x: np.power(x, 0.25))

    sample_labels = list(prop_matrix.index)
    seq_labels = list(prop_matrix.columns)

    log_info(f"Proportional matrix: {prop_matrix.shape[0]} samples x {prop_matrix.shape[1]} sequences")
    log_info("Applied 4th-root transformation for ordination")

    # --- Sample ordination (on 4th-root transformed proportions) ---
    log_info("Computing Bray-Curtis distances for samples")
    sample_dist = compute_bray_curtis_matrix(transformed.values)
    sample_dist_df = pd.DataFrame(sample_dist, index=sample_labels, columns=sample_labels)

    sample_tsne_df = run_ordination_by_block(sample_dist, sample_labels, assays,
                                             "Samples", n_cpus)

    # --- ASV ordination (transpose of 4th-root transformed proportions) ---
    log_info("Computing Bray-Curtis distances for ASVs")
    transformed_t = transformed.T
    seq_dist = compute_bray_curtis_matrix(transformed_t.values)
    seq_dist_df = pd.DataFrame(seq_dist, index=seq_labels, columns=seq_labels)

    seq_tsne_df = run_ordination(seq_dist, seq_labels, "ASVs", n_cpus)

    # --- Save outputs ---
    log_info("Saving sample_bray_tsne.pkl")
    with open("sample_bray_tsne.pkl", "wb") as f:
        pickle.dump(sample_tsne_df, f)

    log_info("Saving seq_bray_tsne.pkl")
    with open("seq_bray_tsne.pkl", "wb") as f:
        pickle.dump(seq_tsne_df, f)

    log_info("Saving sample_bray_dist.pkl")
    with open("sample_bray_dist.pkl", "wb") as f:
        pickle.dump(sample_dist_df, f)

    log_info("Saving seq_bray_dist.pkl")
    with open("seq_bray_dist.pkl", "wb") as f:
        pickle.dump(seq_dist_df, f)

    log_info("cluster_tsne.py complete")


if __name__ == "__main__":
    main()
