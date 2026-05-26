#!/usr/bin/env python3
"""
Augment illumina_rna viz JSONs with mag_analysis-derived functional info:
  - CheckM2 quality per MAG
  - Bakta gene annotations (locus_tag → product/gene/dbxrefs)
  - Computed t-SNE embedding per MAG on the log2(count+1) gene × sample matrix

Writes:
  - functions.json.gz   (gene × sample expression with annotations)
  - mag_quality.json    (per-MAG checkm2 stats merged with the existing references)
  - tsne.json.gz        (per-MAG 2D embedding)

Also rewrites references.json to add quality fields.
"""
import argparse, csv, gzip, json, os, re, sys
from collections import defaultdict

try:
    import numpy as np
    from sklearn.manifold import TSNE
except ImportError:
    print("sklearn/numpy missing — run inside the dana-illumina-rna-rnaseq env", file=sys.stderr)
    sys.exit(1)


def parse_checkm2(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        header = next(csv.reader(f, delimiter="\t"))
        idx = {h: i for i, h in enumerate(header)}
        for row in csv.reader(f, delimiter="\t"):
            if not row or len(row) < len(header):
                continue
            name = row[idx["Name"]]
            out[name] = {
                "completeness":   float(row[idx["Completeness"]]),
                "contamination":  float(row[idx["Contamination"]]),
                "coding_density": float(row[idx["Coding_Density"]]),
                "contig_n50":     int(row[idx["Contig_N50"]]),
                "avg_gene_len":   float(row[idx["Average_Gene_Length"]]),
                "genome_size":    int(row[idx["Genome_Size"]]),
                "gc":             float(row[idx["GC_Content"]]),
                "n_cds":          int(row[idx["Total_Coding_Sequences"]]),
                "n_contigs":      int(row[idx["Total_Contigs"]]),
                "max_contig":     int(row[idx["Max_Contig_Length"]]),
            }
    return out


def parse_bakta_tsv(path):
    """Return {locus_tag: {seqid, type, start, stop, strand, gene, product, dbxrefs}}."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            seqid, ftype, start, stop, strand, locus, gene, product = parts[:8]
            dbxrefs = parts[8] if len(parts) > 8 else ""
            if not locus:
                continue
            out[locus] = {
                "seqid":   seqid.split()[0],   # strip Megahit trailing metadata
                "type":    ftype,
                "start":   int(start),
                "stop":    int(stop),
                "strand":  strand,
                "gene":    gene or None,
                "product": product or None,
                "dbxrefs": dbxrefs or None,
            }
    return out


def parse_contig2bin(path):
    """Return {contig: bin_name}."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
    return out


def load_gene_counts(path):
    """Return (samples, gene_ids, matrix_int)."""
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        gene_ids, rows = [], []
        for r in reader:
            if not r:
                continue
            gene_ids.append(r[0])
            try:
                rows.append([int(x) for x in r[1:]])
            except ValueError:
                rows.append([0] * len(samples))
    return samples, gene_ids, np.array(rows, dtype=np.int32)


def compute_tsne(matrix, n_genes_max=500, perplexity=None):
    """t-SNE on log2(count+1)-transformed top-variance genes."""
    if matrix.size == 0 or matrix.shape[0] < 5:
        return None, []
    lg = np.log2(matrix.astype(np.float32) + 1.0)
    # rank by row variance, keep top n
    var = lg.var(axis=1)
    top = np.argsort(-var)[: min(n_genes_max, lg.shape[0])]
    sub = lg[top]
    # If all rows are zero-variance (no expression), skip
    if sub.var() < 1e-10:
        return None, top.tolist()
    perp = perplexity or max(5, min(30, sub.shape[0] // 4))
    tsne = TSNE(
        n_components=2,
        perplexity=perp,
        learning_rate="auto",
        init="pca",
        max_iter=500,
        random_state=42,
    )
    coords = tsne.fit_transform(sub)
    return coords.tolist(), top.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--viz-dir",     required=True, help="illumina_rna_out/viz dir to enrich")
    ap.add_argument("--mag-dir",     required=True, help="mag_analysis_out root")
    ap.add_argument("--expression-dir", required=True, help="illumina_rna_out/expression root")
    args = ap.parse_args()

    checkm2  = parse_checkm2(os.path.join(args.mag_dir, "binning/checkm2/quality_report.tsv"))
    bakta    = parse_bakta_tsv(os.path.join(args.mag_dir, "annotation/bakta/basic/annotation.tsv"))
    c2b      = parse_contig2bin(os.path.join(args.mag_dir, "binning/dastool/contig2bin.tsv"))

    print(f"[augment] checkm2:    {len(checkm2)} bins", file=sys.stderr)
    print(f"[augment] bakta:      {len(bakta)} CDS annotations", file=sys.stderr)
    print(f"[augment] contig2bin: {len(c2b)} contig→bin links", file=sys.stderr)

    # ── enrich references.json with quality ──────────────────────────────────
    refs_path = os.path.join(args.viz_dir, "references.json")
    with open(refs_path) as f:
        refs = json.load(f)
    for r in refs:
        r["quality"] = checkm2.get(r["name"])
    with open(refs_path, "w") as f:
        json.dump(refs, f)
    print(f"[augment] wrote references.json with quality on {sum(1 for r in refs if r['quality'])}/{len(refs)} MAGs", file=sys.stderr)

    # ── per-MAG functions + t-SNE ────────────────────────────────────────────
    functions = {}   # ref -> {samples, genes, matrix, annotations[]}
    tsne_out  = {}   # ref -> {coords, gene_ids}

    for ref in refs:
        name = ref["name"]
        counts_path = os.path.join(args.expression_dir, name, f"{name}.gene_counts.tsv")
        if not os.path.exists(counts_path):
            continue
        samples, gene_ids, mat = load_gene_counts(counts_path)
        if not gene_ids:
            continue

        annotations = []
        for g in gene_ids:
            a = bakta.get(g)
            if a:
                annotations.append({
                    "product": a["product"],
                    "gene":    a["gene"],
                    "seqid":   a["seqid"],
                })
            else:
                annotations.append({"product": None, "gene": None, "seqid": None})

        functions[name] = {
            "samples":     samples,
            "genes":       gene_ids,
            "matrix":      mat.tolist(),
            "annotations": annotations,
        }

        coords, idx = compute_tsne(mat)
        if coords is not None:
            tsne_out[name] = {
                "samples":   samples,
                "gene_ids":  [gene_ids[i] for i in idx],
                "products":  [annotations[i]["product"] for i in idx],
                "totals":    [int(mat[i].sum()) for i in idx],
                "coords":    coords,
            }

    # gzip these — they can get large
    with gzip.open(os.path.join(args.viz_dir, "functions.json.gz"), "wt") as f:
        json.dump(functions, f)
    with gzip.open(os.path.join(args.viz_dir, "tsne.json.gz"), "wt") as f:
        json.dump(tsne_out, f)

    print(f"[augment] wrote functions.json.gz for {len(functions)} MAGs", file=sys.stderr)
    print(f"[augment] wrote tsne.json.gz for {len(tsne_out)} MAGs", file=sys.stderr)


if __name__ == "__main__":
    main()
