#!/usr/bin/env python3
"""
build_viz.py - Convert pipeline outputs (pickle files) into JSON files for the
Svelte frontend visualization.

Replaces the R build_shiny.R script.

Usage:
  build_viz.py <seqtab.pkl> <renorm.pkl> <taxonomy_dir> <metadata.pkl_or_NONE> \
               <sample_tsne.pkl> <seq_tsne.pkl> <network.pkl>

Outputs (written to current working directory):
  samples.json        Sample metadata with t-SNE coordinates
  asvs.json.gz        ASV info with t-SNE coordinates and taxonomy
  counts.json.gz      Sparse count matrix (gzipped)
  network.json        Correlation edge list
  taxonomy.json       Per-database taxonomy assignments
  renorm_stats.json   Group-level summary statistics
"""

import sys
import os
import json
import gzip
import glob
import pickle
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def log_info(msg):
    print(f"[INFO] {msg}", flush=True)


def log_warn(msg):
    print(f"[WARN] {msg}", flush=True)


def log_error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr, flush=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_pickle(path):
    """Load a pickle file and return its contents."""
    log_info(f"Loading {path}")
    with open(path, "rb") as fh:
        return pickle.load(fh)


def write_json(obj, path, compress=False):
    """Write an object as JSON, optionally gzip-compressed.

    When compress=True, writes both .json.gz and .json versions.
    When compress=False, writes only .json.
    """
    json_str = json.dumps(obj, separators=(",", ":"), allow_nan=False)

    if compress:
        gz_path = path if path.endswith(".gz") else path + ".gz"
        with gzip.open(gz_path, "wt", encoding="utf-8") as fh:
            fh.write(json_str)
        log_info(f"Wrote {gz_path} ({os.path.getsize(gz_path):,} bytes gzipped)")
    else:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(json_str)
        log_info(f"Wrote {path} ({os.path.getsize(path):,} bytes)")


def round_float(val, decimals=4):
    """Round a float value, handling NaN/Inf gracefully."""
    if val is None or (isinstance(val, float) and (np.isnan(val) or np.isinf(val))):
        return 0.0
    return round(float(val), decimals)


# ---------------------------------------------------------------------------
# Build samples.json
# ---------------------------------------------------------------------------

def build_samples(seqtab, sample_tsne, metadata, primer_assignment=None):
    """Build the samples array with t-SNE coordinates and metadata.

    `primer_assignment` is the per-sample record of what cutadapt actually
    matched (see parse_primer_assignment.py). It merges exactly like metadata
    does, and is optional in both directions: absent for pre-trimmed input, and
    absent entirely for pipelines with no amplification step, so consumers must
    treat missing assay_* fields as normal.

    Returns:
        list of dicts, one per sample
    """
    log_info("Building samples.json")

    # Aggregate per-sample stats from count table
    sample_stats = seqtab.groupby("sample").agg(
        total_reads=("count", "sum"),
        n_asvs=("sequence", "nunique"),
    ).reset_index()

    # Build t-SNE lookup: label -> (x, y)
    tsne_lookup = {}
    for _, row in sample_tsne.iterrows():
        tsne_lookup[row["label"]] = (
            round_float(row["tSNE1"]),
            round_float(row["tSNE2"]),
        )

    # Build metadata lookup if available
    meta_lookup = {}
    meta_fields = []
    if metadata is not None:
        if isinstance(metadata, pd.DataFrame):
            # Determine the sample ID column
            id_col = None
            for candidate in ("sample", "Sample", "sample_id", "SampleID", "id"):
                if candidate in metadata.columns:
                    id_col = candidate
                    break
            if id_col is None and metadata.index.name:
                # Use index as sample ID
                metadata = metadata.reset_index()
                id_col = metadata.columns[0]
            elif id_col is None:
                id_col = metadata.columns[0]

            meta_fields = [c for c in metadata.columns if c != id_col]
            for _, row in metadata.iterrows():
                sid = str(row[id_col])
                meta_lookup[sid] = {
                    col: _safe_json_value(row[col]) for col in meta_fields
                }
            log_info(f"Metadata: {len(meta_lookup)} samples, fields: {meta_fields}")

    # Observed primer assignment, keyed by the same sample ids
    assay_lookup = {}
    if primer_assignment is not None and not primer_assignment.empty:
        id_col = "sample" if "sample" in primer_assignment.columns else primer_assignment.columns[0]
        cols = [c for c in primer_assignment.columns if c != id_col]
        for _, row in primer_assignment.iterrows():
            assay_lookup[str(row[id_col])] = {
                c: _safe_json_value(row[c]) for c in cols
                if row[c] is not None and str(row[c]) != "" and not pd.isna(row[c])
            }
        log_info(f"Primer assignment: {len(assay_lookup)} samples, fields: {cols}")

    # Assemble sample records
    records = []
    for _, row in sample_stats.iterrows():
        sid = row["sample"]
        x, y = tsne_lookup.get(sid, (0.0, 0.0))
        rec = {
            "id": sid,
            "x": x,
            "y": y,
            "total_reads": int(row["total_reads"]),
            "n_asvs": int(row["n_asvs"]),
        }
        # Merge metadata fields
        if sid in meta_lookup:
            rec.update(meta_lookup[sid])
        # …then the observed assay, same mechanism, different source
        if sid in assay_lookup:
            rec.update(assay_lookup[sid])
        records.append(rec)

    log_info(f"samples.json: {len(records)} samples")
    return records


def _safe_json_value(val):
    """Convert a value to a JSON-safe type."""
    if val is None:
        return None
    if isinstance(val, (np.integer,)):
        return int(val)
    if isinstance(val, (np.floating,)):
        if np.isnan(val) or np.isinf(val):
            return None
        return round(float(val), 4)
    if isinstance(val, float):
        if np.isnan(val) or np.isinf(val):
            return None
        return round(val, 4)
    if isinstance(val, (np.bool_,)):
        return bool(val)
    return str(val)


# ---------------------------------------------------------------------------
# Build asvs.json
# ---------------------------------------------------------------------------

def build_asvs(seqtab, seq_tsne, renorm_table_list, taxonomy_dict):
    """Build the ASVs array with t-SNE coordinates, group, and taxonomy.

    Returns:
        tuple of (list of dicts, list of ASV IDs in order)
    """
    log_info("Building asvs.json")

    # Sort sequences by total abundance (descending) so ASV_000001 is most abundant
    seq_abundance = seqtab.groupby("sequence")["count"].sum().sort_values(ascending=False)
    sequences = list(seq_abundance.index)
    seq_to_id = {seq: f"ASV_{i+1:06d}" for i, seq in enumerate(sequences)}

    # Aggregate per-ASV stats
    asv_stats = seqtab.groupby("sequence").agg(
        total_reads=("count", "sum"),
        n_samples=("sample", "nunique"),
    ).reset_index()

    # t-SNE lookup
    tsne_lookup = {}
    for _, row in seq_tsne.iterrows():
        tsne_lookup[row["label"]] = (
            round_float(row["tSNE1"]),
            round_float(row["tSNE2"]),
        )

    # Group lookup from renorm
    group_lookup = {}
    if renorm_table_list is not None:
        for _, row in renorm_table_list.iterrows():
            group_lookup[row["sequence"]] = row["group"]

    # Taxonomy lookup: pick first available database for the concatenated string
    tax_string_lookup = {}
    if taxonomy_dict:
        # Use the first database
        first_db = next(iter(taxonomy_dict))
        tax_data = taxonomy_dict[first_db]
        if "tax" in tax_data and isinstance(tax_data["tax"], pd.DataFrame):
            tax_df = tax_data["tax"]
            for seq in tax_df.index:
                vals = [str(v) if pd.notna(v) else "" for v in tax_df.loc[seq]]
                tax_string_lookup[seq] = ";".join(vals)

    # Assemble ASV records
    records = []
    asv_id_list = []
    for _, row in asv_stats.iterrows():
        seq = row["sequence"]
        asv_id = seq_to_id[seq]
        x, y = tsne_lookup.get(seq, (0.0, 0.0))

        rec = {
            "id": asv_id,
            "sequence": seq,
            "x": x,
            "y": y,
            "total_reads": int(row["total_reads"]),
            "n_samples": int(row["n_samples"]),
            "group": group_lookup.get(seq, "unknown"),
            "taxonomy": tax_string_lookup.get(seq, ""),
        }
        records.append(rec)
        asv_id_list.append(asv_id)

    # Sort by ASV ID so ASV_000001 (most abundant) is first
    records.sort(key=lambda r: r["id"])
    log_info(f"asvs.json: {len(records)} ASVs")
    return records, seq_to_id


# ---------------------------------------------------------------------------
# Build counts.json.gz
# ---------------------------------------------------------------------------

def build_counts(seqtab, seq_to_id):
    """Build the sparse count matrix.

    Returns:
        dict with 'data', 'samples', 'asvs' keys
    """
    log_info("Building counts.json.gz")

    # Ordered sample and ASV lists
    samples = sorted(seqtab["sample"].unique())
    asvs = sorted(seq_to_id.keys(), key=lambda s: seq_to_id[s])
    asv_ids = [seq_to_id[s] for s in asvs]

    sample_idx = {s: i for i, s in enumerate(samples)}
    asv_idx = {s: i for i, s in enumerate(asvs)}

    # Compute per-sample totals for proportions
    sample_totals = seqtab.groupby("sample")["count"].sum()

    # Build sparse data array: [sample_idx, asv_idx, count, proportion]
    data = []
    for _, row in seqtab.iterrows():
        cnt = int(row["count"])
        if cnt == 0:
            continue
        si = sample_idx[row["sample"]]
        ai = asv_idx[row["sequence"]]
        total = sample_totals[row["sample"]]
        prop = round_float(cnt / total if total > 0 else 0.0, 6)
        data.append([si, ai, cnt, prop])

    result = {
        "data": data,
        "samples": samples,
        "asvs": asv_ids,
    }

    log_info(f"counts.json.gz: {len(data)} non-zero entries, "
             f"{len(samples)} samples, {len(asv_ids)} ASVs")
    return result


# ---------------------------------------------------------------------------
# Build aggregated counts per taxonomy level
# ---------------------------------------------------------------------------

def build_aggregated_counts(seqtab, seq_to_id, taxonomy_dict):
    """Pre-aggregate counts by each taxonomy level + group.

    Returns dict keyed by level name, each containing:
        {data: [[si, taxon_idx, count, proportion], ...], samples: [...], taxa: [...]}
    """
    log_info("Building aggregated counts")

    db = next(iter(taxonomy_dict), None)
    if not db or not taxonomy_dict[db]:
        return {}

    tax_data = taxonomy_dict[db]
    if "tax" in tax_data and isinstance(tax_data["tax"], pd.DataFrame):
        tax_df = tax_data["tax"]
    else:
        return {}

    levels = list(tax_df.columns)
    samples = sorted(seqtab["sample"].unique())
    sample_idx = {s: i for i, s in enumerate(samples)}
    sample_totals = seqtab.groupby("sample")["count"].sum()

    result = {}

    for level_name in levels:
        # Map each sequence to its taxon at this level
        seq_to_taxon = {}
        for seq in seqtab["sequence"].unique():
            if seq in tax_df.index:
                val = tax_df.loc[seq, level_name]
                seq_to_taxon[seq] = str(val) if pd.notna(val) and val else "unclassified"
            else:
                seq_to_taxon[seq] = "unclassified"

        # Aggregate counts: (sample, taxon) -> total count
        agg = seqtab.copy()
        agg["taxon"] = agg["sequence"].map(seq_to_taxon)
        grouped = agg.groupby(["sample", "taxon"])["count"].sum().reset_index()

        # Build sparse array
        taxa = sorted(grouped["taxon"].unique())
        taxon_idx = {t: i for i, t in enumerate(taxa)}

        data = []
        for _, row in grouped.iterrows():
            cnt = int(row["count"])
            if cnt == 0:
                continue
            si = sample_idx[row["sample"]]
            ti = taxon_idx[row["taxon"]]
            total = sample_totals[row["sample"]]
            prop = round_float(cnt / total if total > 0 else 0.0, 6)
            data.append([si, ti, cnt, prop])

        result[level_name] = {
            "data": data,
            "samples": samples,
            "taxa": taxa,
        }
        log_info(f"  {level_name}: {len(taxa)} taxa, {len(data)} entries")

    # Also add 'group' level using the ASV group field (prokaryote/eukaryote/etc)
    # We need the group info from the renorm data or asvs
    log_info(f"  Aggregated {len(levels)} levels")
    return result


# ---------------------------------------------------------------------------
# Build network.json
# ---------------------------------------------------------------------------

def build_network(network_df, seq_to_id):
    """Build the network edge list.

    Returns:
        dict with 'edges' key
    """
    log_info("Building network.json")

    # Build ASV index lookup
    asv_ids_sorted = sorted(seq_to_id.values())
    asv_id_to_idx = {aid: i for i, aid in enumerate(asv_ids_sorted)}

    # Also need sequence-to-ASV-index
    seq_to_idx = {}
    for seq, aid in seq_to_id.items():
        if aid in asv_id_to_idx:
            seq_to_idx[seq] = asv_id_to_idx[aid]

    edges = []
    if network_df is not None and len(network_df) > 0:
        for _, row in network_df.iterrows():
            n1 = row["node1"]
            n2 = row["node2"]
            corr = round_float(row["correlation"])

            idx1 = seq_to_idx.get(n1)
            idx2 = seq_to_idx.get(n2)
            if idx1 is not None and idx2 is not None:
                edges.append([idx1, idx2, corr])

    log_info(f"network.json: {len(edges)} edges")
    return {"edges": edges}


# ---------------------------------------------------------------------------
# Build taxonomy.json
# ---------------------------------------------------------------------------

def build_taxonomy(taxonomy_dict, seq_to_id):
    """Build per-database taxonomy assignments.

    Returns:
        dict keyed by database name
    """
    log_info("Building taxonomy.json")

    result = {}
    for db_name, tax_data in taxonomy_dict.items():
        if "tax" in tax_data and isinstance(tax_data["tax"], pd.DataFrame):
            tax_df = tax_data["tax"]
        elif isinstance(tax_data, pd.DataFrame):
            tax_df = tax_data
        else:
            log_warn(f"Skipping taxonomy database '{db_name}': unexpected format")
            continue

        levels = list(tax_df.columns)
        assignments = {}
        for seq in tax_df.index:
            if seq in seq_to_id:
                asv_id = seq_to_id[seq]
                vals = [str(v) if pd.notna(v) else "" for v in tax_df.loc[seq]]
                assignments[asv_id] = vals

        # Include bootstraps if available
        bootstraps = {}
        boot_df = tax_data.get("boot") if isinstance(tax_data, dict) else None
        if boot_df is not None and isinstance(boot_df, pd.DataFrame):
            for seq in boot_df.index:
                if seq in seq_to_id:
                    asv_id = seq_to_id[seq]
                    bootstraps[asv_id] = [int(v) if pd.notna(v) else 0 for v in boot_df.loc[seq]]

        entry = {
            "levels": levels,
            "assignments": assignments,
        }
        if bootstraps:
            entry["bootstraps"] = bootstraps

        result[db_name] = entry
        log_info(f"  {db_name}: {len(levels)} levels, "
                 f"{len(assignments)} ASVs assigned")

    return result


# ---------------------------------------------------------------------------
# Build renorm_stats.json
# ---------------------------------------------------------------------------

def build_renorm_stats(renorm_data):
    """Build group-level summary statistics.

    Accepts either:
      - renorm_merged DataFrame (with group, count columns)
      - renorm_by_group dict of DataFrames

    Returns:
        dict keyed by group name
    """
    log_info("Building renorm_stats.json")

    result = {}

    if isinstance(renorm_data, dict):
        # renorm_by_group format
        for grp, df in renorm_data.items():
            result[grp] = {
                "n_asvs": int(df["sequence"].nunique()) if "sequence" in df.columns else 0,
                "n_samples": int(df["sample"].nunique()) if "sample" in df.columns else 0,
                "n_reads": int(df["count"].sum()) if "count" in df.columns else 0,
            }
    elif isinstance(renorm_data, pd.DataFrame):
        # renorm_merged format
        if "group" in renorm_data.columns:
            for grp, sub in renorm_data.groupby("group"):
                result[grp] = {
                    "n_asvs": int(sub["sequence"].nunique()),
                    "n_samples": int(sub["sample"].nunique()),
                    "n_reads": int(sub["count"].sum()),
                }
        else:
            log_warn("renorm data has no 'group' column; stats will be empty")
    else:
        log_warn(f"Unexpected renorm data type: {type(renorm_data)}")

    for grp, stats in result.items():
        log_info(f"  {grp}: {stats['n_asvs']} ASVs, "
                 f"{stats['n_reads']:,} reads")

    return result


# ---------------------------------------------------------------------------
# Load taxonomy files from a directory
# ---------------------------------------------------------------------------

def load_taxonomy_dir(taxonomy_dir):
    """Load all *_taxonomy.pkl files from a directory.

    Returns:
        dict of {db_name: tax_data}
    """
    taxonomy_dict = {}

    # Search provided directory
    patterns = [
        os.path.join(taxonomy_dir, "*_taxonomy.pkl"),
        os.path.join(taxonomy_dir, "*_taxonomy.pickle"),
    ]
    tax_files = []
    for pattern in patterns:
        tax_files.extend(glob.glob(pattern))

    # Fallback: search current working directory
    if not tax_files:
        log_warn(f"No taxonomy files in {taxonomy_dir}, searching current directory")
        tax_files = glob.glob("*_taxonomy.pkl")

    if not tax_files:
        log_warn("No taxonomy files found anywhere")
        return taxonomy_dict

    for tax_file in sorted(tax_files):
        # Extract database name from filename: <db_name>_taxonomy.pkl
        basename = os.path.basename(tax_file)
        db_name = basename.replace("_taxonomy.pkl", "").replace("_taxonomy.pickle", "")
        log_info(f"Loading taxonomy: {tax_file} (db={db_name})")
        try:
            tax_data = load_pickle(tax_file)
            taxonomy_dict[db_name] = tax_data
        except Exception as e:
            log_warn(f"Failed to load {tax_file}: {e}")

    return taxonomy_dict


# ---------------------------------------------------------------------------
# Build heatmap.json.gz — pre-clustered matrix with dendrogram coordinates
# ---------------------------------------------------------------------------

def build_heatmap(seqtab, seq_to_id, taxonomy_dict):
    """Build clustered heatmap data with scipy dendrograms.

    Returns dict with: z, sampleIds, asvIds, rowDendro, colDendro, colColors
    """
    from scipy.spatial.distance import braycurtis, pdist, squareform
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, fcluster, to_tree

    log_info("Building heatmap.json.gz")

    # Build wide matrix: samples x ASVs (4th-root relative abundance)
    samples = sorted(seqtab["sample"].unique())
    seqs = sorted(seq_to_id.keys(), key=lambda s: seq_to_id[s])
    asv_ids = [seq_to_id[s] for s in seqs]

    sample_idx = {s: i for i, s in enumerate(samples)}
    seq_idx = {s: i for i, s in enumerate(seqs)}

    mat = np.zeros((len(samples), len(seqs)))
    for _, row in seqtab.iterrows():
        si = sample_idx.get(row["sample"])
        ji = seq_idx.get(row["sequence"])
        if si is not None and ji is not None:
            mat[si, ji] = row["count"]

    # Convert to relative abundance, then 4th root
    row_sums = mat.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    rel = mat / row_sums
    transformed = np.power(rel, 0.25)

    # Cluster rows (samples) by Bray-Curtis
    if len(samples) > 1:
        row_dist = pdist(transformed, metric='braycurtis')
        row_dist = np.nan_to_num(row_dist, nan=1.0)
        row_link = linkage(row_dist, method='ward')
        row_dendro = dendrogram(row_link, no_plot=True)
        row_order = row_dendro['leaves']
        # Extract dendrogram coordinates
        row_dendro_data = {
            'icoord': [list(map(float, x)) for x in row_dendro['icoord']],
            'dcoord': [list(map(float, x)) for x in row_dendro['dcoord']],
        }
    else:
        row_order = [0]
        row_dendro_data = {'icoord': [], 'dcoord': []}

    # Cluster columns (ASVs) by Bray-Curtis on transposed matrix
    if len(seqs) > 1:
        col_dist = pdist(transformed.T, metric='braycurtis')
        col_dist = np.nan_to_num(col_dist, nan=1.0)
        col_link = linkage(col_dist, method='ward')
        col_dendro = dendrogram(col_link, no_plot=True)
        col_order = col_dendro['leaves']
        col_dendro_data = {
            'icoord': [list(map(float, x)) for x in col_dendro['icoord']],
            'dcoord': [list(map(float, x)) for x in col_dendro['dcoord']],
        }
    else:
        col_order = [0]
        col_dendro_data = {'icoord': [], 'dcoord': []}

    # Reorder matrix
    ordered_z = transformed[np.ix_(row_order, col_order)].tolist()
    ordered_samples = [samples[i] for i in row_order]
    ordered_asvs = [asv_ids[i] for i in col_order]

    # ASV labels (deepest taxonomy)
    asv_labels = []
    db = next(iter(taxonomy_dict), None)
    tax_data = taxonomy_dict.get(db, {}) if db else {}
    tax_df = tax_data.get("tax") if isinstance(tax_data, dict) else None

    for asv_id in ordered_asvs:
        label = asv_id
        if tax_df is not None and isinstance(tax_df, pd.DataFrame):
            # Find the sequence for this ASV ID
            seq = None
            for s, aid in seq_to_id.items():
                if aid == asv_id:
                    seq = s
                    break
            if seq and seq in tax_df.index:
                vals = tax_df.loc[seq]
                for v in reversed(list(vals)):
                    if pd.notna(v) and v:
                        label = str(v)
                        break
        asv_labels.append(label)

    # Pre-compute cluster assignments for k=2..max_k
    max_k = min(20, len(samples))
    sample_clusters = {}  # k -> {sampleId: clusterIdx}
    if len(samples) > 1:
        for k in range(2, max_k + 1):
            labels = fcluster(row_link, k, criterion='maxclust')
            sample_clusters[str(k)] = {
                samples[i]: int(labels[i]) for i in range(len(samples))
            }

    max_k_asv = min(20, len(seqs))
    asv_clusters = {}  # k -> {asvId: clusterIdx}
    if len(seqs) > 1:
        for k in range(2, max_k_asv + 1):
            labels = fcluster(col_link, k, criterion='maxclust')
            asv_clusters[str(k)] = {
                asv_ids[i]: int(labels[i]) for i in range(len(seqs))
            }

    # Convert ASV linkage to Newick for the Ward tree view
    def linkage_to_newick(Z, labels):
        """Convert scipy linkage matrix to Newick string."""
        root = to_tree(Z)
        def _to_nwk(node):
            if node.is_leaf():
                return f"{labels[node.id]}:{node.dist:.6f}"
            left = _to_nwk(node.get_left())
            right = _to_nwk(node.get_right())
            return f"({left},{right}):{node.dist:.6f}"
        return _to_nwk(root) + ";"

    asv_ward_newick = ''
    if len(seqs) > 1:
        asv_ward_newick = linkage_to_newick(col_link, asv_ids)

    result = {
        'z': ordered_z,
        'sampleIds': ordered_samples,
        'asvIds': ordered_asvs,
        'asvLabels': asv_labels,
        'rowDendro': row_dendro_data,
        'colDendro': col_dendro_data,
        'sampleClusters': sample_clusters,
        'asvClusters': asv_clusters,
        'asvWardNewick': asv_ward_newick,
        'nSamples': len(samples),
        'nAsvs': len(seqs),
    }

    log_info(f"  heatmap: {len(ordered_samples)} samples x {len(ordered_asvs)} ASVs, "
             f"cluster cuts k=2..{max_k}")
    return result


# ---------------------------------------------------------------------------
# Copy Svelte app
# ---------------------------------------------------------------------------

def copy_svelte_app():
    """Copy the Svelte app's index.html to the output directory if it exists."""
    # Look for the built app relative to this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, "..", "..", "viz", "dist", "index.html"),
        os.path.join(script_dir, "..", "..", "..", "viz", "dist", "index.html"),
        os.path.join(script_dir, "..", "..", "viz", "public", "index.html"),
    ]

    for candidate in candidates:
        candidate = os.path.normpath(candidate)
        if os.path.isfile(candidate):
            import shutil
            dest = os.path.join(os.getcwd(), "index.html")
            shutil.copy2(candidate, dest)
            log_info(f"Copied Svelte app from {candidate}")
            return

    log_warn("Svelte app dist/index.html not found; skipping copy")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 8:
        log_error(
            "Usage: build_viz.py <seqtab.pkl> <renorm.pkl> <taxonomy_dir> "
            "<metadata.pkl_or_NONE> <sample_tsne.pkl> <seq_tsne.pkl> <network.pkl> "
            "[primer_assignment.tsv_or_NONE]"
        )
        sys.exit(1)

    seqtab_path = sys.argv[1]
    renorm_path = sys.argv[2]
    taxonomy_dir = sys.argv[3]
    metadata_path = sys.argv[4]
    sample_tsne_path = sys.argv[5]
    seq_tsne_path = sys.argv[6]
    network_path = sys.argv[7]
    # Optional 8th arg so older callers keep working.
    assignment_path = sys.argv[8] if len(sys.argv) > 8 else None

    # --- Load data ---
    seqtab = load_pickle(seqtab_path)
    if not isinstance(seqtab, pd.DataFrame):
        log_error("seqtab.pkl must contain a pandas DataFrame")
        sys.exit(1)

    # Normalize column names
    for col in ("sample", "sequence", "count"):
        if col not in seqtab.columns:
            log_error(f"seqtab missing required column: {col}")
            sys.exit(1)

    renorm_data = load_pickle(renorm_path)

    # Extract the renorm_table_list (sequence-to-group mapping)
    renorm_table_list = None
    if isinstance(renorm_data, pd.DataFrame) and "group" in renorm_data.columns:
        # renorm_merged format: extract unique sequence->group mapping
        renorm_table_list = renorm_data[["sequence", "group"]].drop_duplicates()
    elif isinstance(renorm_data, dict):
        # Could be renorm_by_group; build the mapping from it
        rows = []
        for grp, df in renorm_data.items():
            if "sequence" in df.columns:
                for seq in df["sequence"].unique():
                    rows.append({"sequence": seq, "group": grp})
        if rows:
            renorm_table_list = pd.DataFrame(rows).drop_duplicates()

    # Load taxonomy
    taxonomy_dict = load_taxonomy_dir(taxonomy_dir)

    # Load metadata (optional)
    metadata = None
    if metadata_path and metadata_path.upper() != "NONE":
        if os.path.isfile(metadata_path):
            try:
                metadata = load_pickle(metadata_path)
                log_info(f"Metadata loaded: {type(metadata)}")
            except Exception as e:
                log_warn(f"Failed to load metadata: {e}")
        else:
            log_warn(f"Metadata file not found: {metadata_path}")

    # t-SNE coords are an optional enhancement: CLUSTER_TSNE runs under
    # errorStrategy 'ignore', so its output may be a `NO_*` sentinel (the
    # channel was empty). An empty frame is fine — build_samples/build_asvs
    # already default un-located points to the origin.
    def load_optional(path, label, cols):
        if os.path.basename(path).startswith("NO_") or not os.path.isfile(path):
            log_warn(f"No {label} ({os.path.basename(path)}) — this view will be degraded")
            return pd.DataFrame(columns=cols)
        try:
            return load_pickle(path)
        except Exception as e:
            log_warn(f"Failed to load {label} ({path}): {e}")
            return pd.DataFrame(columns=cols)

    sample_tsne = load_optional(sample_tsne_path, "sample t-SNE", ["label", "x", "y"])
    seq_tsne = load_optional(seq_tsne_path, "ASV t-SNE", ["label", "x", "y"])

    # Load network (optional — empty/sentinel yields a network with no edges)
    network_df = None
    if os.path.isfile(network_path) and not os.path.basename(network_path).startswith("NO_"):
        try:
            network_df = load_pickle(network_path)
        except Exception as e:
            log_warn(f"Failed to load network ({network_path}): {e}")
            network_df = None
        if isinstance(network_df, pd.DataFrame):
            log_info(f"Network: {len(network_df)} edges")
        elif network_df is not None:
            log_warn(f"Network file has unexpected type: {type(network_df)}")
            network_df = None
    else:
        log_warn(f"No network ({os.path.basename(network_path)}) — ASV network will be empty")

    # Observed primer assignment (optional — absent when input was pre-trimmed,
    # and absent by nature for pipelines with no amplification step).
    primer_assignment = None
    if (assignment_path and os.path.isfile(assignment_path)
            and not os.path.basename(assignment_path).startswith("NO_")):
        try:
            primer_assignment = pd.read_csv(assignment_path, sep="\t")
            log_info(f"Primer assignment: {len(primer_assignment)} rows "
                     f"from {os.path.basename(assignment_path)}")
        except (OSError, ValueError) as e:
            log_warn(f"Could not read primer assignment: {e}")

    # --- Build outputs ---

    # samples.json
    samples = build_samples(seqtab, sample_tsne, metadata, primer_assignment)
    write_json(samples, "samples.json")

    # asvs.json.gz
    asvs, seq_to_id = build_asvs(seqtab, seq_tsne, renorm_table_list, taxonomy_dict)
    write_json(asvs, "asvs.json.gz", compress=True)

    # counts.json.gz
    counts = build_counts(seqtab, seq_to_id)
    write_json(counts, "counts.json.gz", compress=True)

    # network.json
    network = build_network(network_df, seq_to_id)
    write_json(network, "network.json")

    # taxonomy.json
    taxonomy = build_taxonomy(taxonomy_dict, seq_to_id)
    write_json(taxonomy, "taxonomy.json")

    # renorm_stats.json
    renorm_stats = build_renorm_stats(renorm_data)
    write_json(renorm_stats, "renorm_stats.json")

    # heatmap.json.gz
    heatmap = build_heatmap(seqtab, seq_to_id, taxonomy_dict)
    write_json(heatmap, "heatmap.json.gz", compress=True)

    # aggregated_counts.json.gz — pre-aggregated per taxonomy level
    agg_counts = build_aggregated_counts(seqtab, seq_to_id, taxonomy_dict)
    write_json(agg_counts, "aggregated_counts.json.gz", compress=True)

    # --- Copy Svelte app ---
    copy_svelte_app()

    log_info("build_viz.py complete")


if __name__ == "__main__":
    main()
