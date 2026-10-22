#!/usr/bin/env python3
"""
Merge per-sample featureCounts outputs for one reference into matrices and
aggregate them to MAGs and functions, plus the viz JSONs.

Per-sample inputs share one row order (one SAF per reference), so they are
read as bare count columns and stacked.

Outputs (in --out):
  <ref>.gene_counts.mtx.gz      genes × samples, MatrixMarket (sparse; rows = genes.tsv.gz order)
  <ref>.gene_counts.tsv.gz      dense copy, only when genes × samples is small
  <ref>.contig_counts.mtx.gz    contigs × samples
  <ref>.samples.tsv             column order of the .mtx files
  <ref>.sample_summary.tsv      read flow, mapping and assignment per sample
  <ref>.bin_counts.tsv          read pairs on each bin's contigs × samples
  <ref>.bin_tpm.tsv             sum of gene TPM per bin × samples
  <ref>.<ko|cog|cazy>_tpm.tsv   function TPM × samples (when annotated)
  viz/<ref>.*.json.gz           small tables for the viz
"""
import argparse, csv, glob, gzip, json, os, re, sys
from collections import defaultdict

import numpy as np
import pandas as pd

DENSE_LIMIT = 2e8   # cells; above this only the sparse matrix is written


def log(msg):
    print(f"[summarize] {msg}", file=sys.stderr, flush=True)


def sample_of(path, ref):
    base = os.path.basename(path)
    return re.sub(rf"_vs_{re.escape(ref)}\..*$", "", base)


def read_count_column(path, n_expected=None):
    col = pd.read_csv(path, sep="\t", usecols=[1], dtype=np.int32, engine="c").iloc[:, 0].to_numpy()
    if n_expected is not None and len(col) != n_expected:
        raise SystemExit(f"[summarize] ERROR: {path} has {len(col)} rows, expected {n_expected}")
    return col


def stack_counts(files, ref, n_rows):
    samples = [sample_of(f, ref) for f in files]
    order = np.argsort(samples)
    samples = [samples[i] for i in order]
    files = [files[i] for i in order]
    mat = np.zeros((n_rows, len(files)), dtype=np.int32)
    for j, f in enumerate(files):
        mat[:, j] = read_count_column(f, n_rows)
        if (j + 1) % 20 == 0:
            log(f"  read {j + 1}/{len(files)}")
    return samples, mat


def write_mtx(path, mat):
    nnz = int(np.count_nonzero(mat))
    with gzip.open(path, "wt", compresslevel=3) as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{mat.shape[0]} {mat.shape[1]} {nnz}\n")
        for j in range(mat.shape[1]):
            rows = np.flatnonzero(mat[:, j])
            if len(rows):
                block = np.column_stack([rows + 1, np.full(len(rows), j + 1), mat[rows, j]])
                np.savetxt(f, block, fmt="%d")
    return nnz


def tpm(counts, lengths_bp):
    # float32: a co-assembly matrix is ~10M genes x 100+ samples
    rpk = counts.astype(np.float32) / (np.maximum(lengths_bp, 1)[:, None] / 1e3).astype(np.float32)
    tot = rpk.sum(axis=0)
    tot[tot == 0] = 1
    return rpk / tot * 1e6


def aggregate(idx, n_groups, mat):
    """Sum rows of mat into groups; idx < 0 means 'no group'."""
    keep = idx >= 0
    out = np.zeros((n_groups, mat.shape[1]), dtype=np.float64)
    for j in range(mat.shape[1]):
        out[:, j] = np.bincount(idx[keep], weights=mat[keep, j], minlength=n_groups)
    return out


def write_table(path, row_name, rows, samples, mat, fmt="%.4g"):
    df = pd.DataFrame(mat, index=pd.Index(rows, name=row_name), columns=samples)
    df.to_csv(path, sep="\t", float_format=fmt)


def parse_flagstat(path):
    out = {}
    with open(path) as f:
        for line in f:
            m = re.match(r"(\d+) \+ \d+ (.*)", line.strip())
            if not m:
                continue
            n, what = int(m.group(1)), m.group(2)
            if what.startswith("primary mapped"):
                out["mapped_primary"] = n
            elif what.startswith("properly paired"):
                out["properly_paired"] = n
    return out


def parse_fc_summary(path):
    out = {}
    with open(path) as f:
        next(f, None)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1].isdigit():
                out[p[0]] = int(p[1])
    return out


def load_bins(path, contig_index):
    if not path or not os.path.getsize(path):
        return None, None
    bins, idx = [], np.full(len(contig_index), -1, dtype=np.int64)
    bin_id = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or p[0] in ("contig", "contig_id", "scaffold"):
                continue
            ci = contig_index.get(p[0])
            if ci is None or idx[ci] >= 0:
                continue
            if p[1] not in bin_id:
                bin_id[p[1]] = len(bins)
                bins.append(p[1])
            idx[ci] = bin_id[p[1]]
    return bins, idx


def load_tsv_dict(path, key, cols):
    if not path or not os.path.exists(path) or not os.path.getsize(path):
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    if key not in df.columns:
        return {}
    keep = [c for c in cols if c in df.columns]
    return df.set_index(key)[keep].to_dict("index")


def functions_from_annotations(path, gene_ids):
    """KO / COG category / CAZy per gene from mag_analysis merged_annotations.tsv,
    falling back to nothing. Returns {kind: list-of-label-per-gene}."""
    if not path or not os.path.getsize(path):
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, usecols=lambda c: c in
                     ("protein_id", "KO", "COG_category", "CAZy"))
    df = df.drop_duplicates("protein_id").set_index("protein_id")
    out = {}
    for col, kind in (("KO", "ko"), ("COG_category", "cog"), ("CAZy", "cazy")):
        if col not in df.columns:
            continue
        s = df[col].reindex(gene_ids)
        # one label per gene: the first of a comma list; COG keeps the first letter
        s = s.fillna("").str.split(r"[,;]", n=1).str[0].str.strip()
        s = s.replace({"-": "", "nan": ""}).str.replace(r"^ko:", "", regex=True)
        if kind == "cog":
            s = s.str[:1]
        out[kind] = s.to_numpy()
    return out


def functions_from_dbxref(dbxref):
    ko = np.array([next((x[5:] for x in d.split(",") if x.startswith("KEGG:")), "") for d in dbxref])
    cog = np.array([next((x[4:] for x in d.split(",") if x.startswith("COG:")), "") for d in dbxref])
    out = {}
    if (ko != "").any():
        out["ko"] = ko
    if (cog != "").any():
        out["cog_id"] = cog
    return out


def load_contig_embedding(viz_dir, contig_index):
    """t-SNE coordinates per contig from mag_analysis's viz data
    (contig_explorer.part-NNN.json[.gz]: delta-encoded contig numbers plus
    columns). Returns (x, y) float arrays aligned to contig_index, NaN where absent."""
    if not viz_dir or not os.path.isdir(viz_dir):
        return None
    parts = sorted(glob.glob(os.path.join(viz_dir, "contig_explorer.part-*.json.gz"))) \
        or sorted(glob.glob(os.path.join(viz_dir, "contig_explorer.part-*.json")))
    man = next((m for m in (os.path.join(viz_dir, "contig_explorer_manifest.json.gz"),
                            os.path.join(viz_dir, "contig_explorer_manifest.json")) if os.path.exists(m)), None)
    if not parts or not man:
        return None
    op = lambda f: gzip.open(f, "rt") if f.endswith(".gz") else open(f)
    with op(man) as fh:
        manifest = json.load(fh)
    if not manifest.get("has_tsne"):
        return None
    enc = manifest.get("id_encoding", {"prefix": "contig_", "width": 7})
    x = np.full(len(contig_index), np.nan, dtype=np.float32)
    y = np.full(len(contig_index), np.nan, dtype=np.float32)
    hit = 0
    for f in parts:
        with op(f) as fh:
            c = json.load(fh)
        nums = np.cumsum(np.asarray(c["ids"], dtype=np.int64))
        tx, ty = c["columns"]["tsne_x"], c["columns"]["tsne_y"]
        for n, a, b in zip(nums, tx, ty):
            i = contig_index.get(f"{enc['prefix']}{int(n):0{enc['width']}d}")
            if i is not None and a is not None and b is not None:
                x[i], y[i] = a, b
                hit += 1
    log(f"contig t-SNE: {hit} of {len(contig_index)} contigs placed from {len(parts)} parts")
    return (x, y) if hit else None


def jdump(path, obj):
    with gzip.open(path, "wt", compresslevel=6) as f:
        json.dump(obj, f, separators=(",", ":"), allow_nan=False)


def r(x, nd=3):
    return float(round(float(x), nd))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--genes-dir", required=True)
    ap.add_argument("--contigs-dir", required=True)
    ap.add_argument("--summary-dir", required=True)
    ap.add_argument("--flagstat-dir", required=True)
    ap.add_argument("--readcounts-dir", required=True)
    ap.add_argument("--gene-annot", required=True)
    ap.add_argument("--contig-lengths", required=True)
    ap.add_argument("--bins")
    ap.add_argument("--checkm2")
    ap.add_argument("--gtdbtk")
    ap.add_argument("--annotations")
    ap.add_argument("--contig-explorer", help="mag_analysis viz/data dir with contig_explorer parts")
    ap.add_argument("--max-map-contigs", type=int, default=400000)
    ap.add_argument("--top-genes", type=int, default=2000)
    ap.add_argument("--top-functions", type=int, default=500)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    ref = a.ref
    os.makedirs(os.path.join(a.out, "viz"), exist_ok=True)

    # ---- reference tables ----
    genes = pd.read_csv(a.gene_annot, sep="\t", dtype={"gene_id": str, "contig": str, "product": str,
                                                         "gene": str, "dbxref": str, "type": str},
                        keep_default_na=False)
    clen = pd.read_csv(a.contig_lengths, sep="\t", dtype={"contig": str})
    contig_index = {c: i for i, c in enumerate(clen["contig"])}
    log(f"{ref}: {len(genes)} genes, {len(clen)} contigs")

    # ---- counts ----
    gene_files = sorted(glob.glob(os.path.join(a.genes_dir, "*.genes.tsv.gz")))
    contig_files = sorted(glob.glob(os.path.join(a.contigs_dir, "*.contigs.tsv.gz")))
    samples, gmat = stack_counts(gene_files, ref, len(genes))
    csamples, cmat = stack_counts(contig_files, ref, len(clen))
    assert samples == csamples, "gene and contig sample sets differ"
    n_s = len(samples)
    log(f"{n_s} samples")

    pd.Series(samples, name="sample").to_csv(os.path.join(a.out, f"{ref}.samples.tsv"), index=False)
    nnz = write_mtx(os.path.join(a.out, f"{ref}.gene_counts.mtx.gz"), gmat)
    log(f"gene matrix: {nnz} non-zero of {gmat.size}")
    write_mtx(os.path.join(a.out, f"{ref}.contig_counts.mtx.gz"), cmat)
    if gmat.size <= DENSE_LIMIT and len(genes):
        pd.DataFrame(gmat, index=pd.Index(genes["gene_id"], name="gene_id"), columns=samples) \
          .to_csv(os.path.join(a.out, f"{ref}.gene_counts.tsv.gz"), sep="\t", compression={"method": "gzip", "compresslevel": 3})

    # ---- per-sample summary ----
    stages = defaultdict(dict)
    for f in glob.glob(os.path.join(a.readcounts_dir, "*")):
        with open(f) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) == 3 and p[2].isdigit():
                    stages[p[0]][p[1]] = int(p[2])
    flag = {sample_of(f, ref): parse_flagstat(f) for f in glob.glob(os.path.join(a.flagstat_dir, "*"))}
    fcs = {sample_of(f, ref): parse_fc_summary(f) for f in glob.glob(os.path.join(a.summary_dir, "*"))}
    stage_order = ["raw", "trimmed", "filtered", "norrna", "nohuman"]
    rows = []
    for s in samples:
        st, fl, fc = stages.get(s, {}), flag.get(s, {}), fcs.get(s, {})
        final = next((st[k] for k in reversed(stage_order) if k in st), None)
        row = {"sample": s}
        for k in stage_order:
            row[f"reads_{k}"] = st.get(k)
        row["reads_mapped"] = fl.get("mapped_primary")
        row["pairs_assigned"] = fc.get("Assigned")
        row["pairs_no_feature"] = fc.get("Unassigned_NoFeatures")
        row["pairs_ambiguous"] = fc.get("Unassigned_Ambiguity")
        row["pct_rrna"] = (100 * (1 - st["norrna"] / st["filtered"])) if st.get("filtered") and "norrna" in st else None
        row["pct_mapped"] = (100 * fl["mapped_primary"] / final) if final and "mapped_primary" in fl else None
        row["pct_assigned"] = (100 * 2 * fc["Assigned"] / final) if final and "Assigned" in fc else None
        rows.append(row)
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(a.out, f"{ref}.sample_summary.tsv"), sep="\t", index=False, float_format="%.3f")

    # ---- TPM + MAG aggregation ----
    g_len = genes["length"].to_numpy(dtype=np.float64)
    g_tpm = tpm(gmat, g_len) if len(genes) else np.zeros((0, n_s), dtype=np.float32)
    g_contig = genes["contig"].map(contig_index).fillna(-1).to_numpy(dtype=np.int64)

    bins, contig_bin = load_bins(a.bins, contig_index)
    bin_rows = []
    if bins:
        n_b = len(bins)
        bin_counts = aggregate(contig_bin, n_b, cmat)
        gene_bin = np.where(g_contig >= 0, contig_bin[np.maximum(g_contig, 0)], -1)
        bin_tpm = aggregate(gene_bin, n_b, g_tpm)
        # unbinned share, so the MAG table accounts for every assigned read
        unb_counts = cmat[contig_bin < 0].sum(axis=0)
        unb_tpm = g_tpm[gene_bin < 0].sum(axis=0)
        write_table(os.path.join(a.out, f"{ref}.bin_counts.tsv"), "bin", bins + ["unbinned"], samples,
                    np.vstack([bin_counts, unb_counts]), fmt="%.0f")
        write_table(os.path.join(a.out, f"{ref}.bin_tpm.tsv"), "bin", bins + ["unbinned"], samples,
                    np.vstack([bin_tpm, unb_tpm]))

        checkm = load_tsv_dict(a.checkm2, "Name", ["Completeness", "Contamination", "Genome_Size", "GC_Content"])
        gtdb = load_tsv_dict(a.gtdbtk, "user_genome", ["classification"])
        bin_len = np.bincount(contig_bin[contig_bin >= 0], weights=clen["length"].to_numpy()[contig_bin >= 0], minlength=n_b)
        bin_ngenes = np.bincount(gene_bin[gene_bin >= 0], minlength=n_b)
        for i, b in enumerate(bins):
            q, t = checkm.get(b, {}), gtdb.get(b, {})
            bin_rows.append({
                "bin": b,
                "length": int(bin_len[i]),
                "n_genes": int(bin_ngenes[i]),
                "completeness": r(q["Completeness"], 2) if q.get("Completeness") else None,
                "contamination": r(q["Contamination"], 2) if q.get("Contamination") else None,
                "classification": t.get("classification"),
                "reads": [int(x) for x in bin_counts[i]],
                "tpm": [r(x, 1) for x in bin_tpm[i]],
            })
        bin_rows.append({"bin": "unbinned", "reads": [int(x) for x in unb_counts],
                         "tpm": [r(x, 1) for x in unb_tpm]})
        log(f"{n_b} bins; unbinned TPM share (mean) {unb_tpm.mean() / 1e4:.1f}%")
    else:
        gene_bin = np.full(len(genes), -1)

    # ---- functions ----
    funcs = functions_from_annotations(a.annotations, genes["gene_id"]) if len(genes) else {}
    if not funcs and len(genes):
        funcs = functions_from_dbxref(genes["dbxref"].to_numpy())
    func_json = {}
    for kind, labels in funcs.items():
        cats, idx = np.unique(labels, return_inverse=True)
        idx = np.where(cats[idx] == "", -1, idx)
        ftpm = aggregate(idx, len(cats), g_tpm)
        keep = cats != ""
        cats, ftpm = cats[keep], ftpm[keep]
        write_table(os.path.join(a.out, f"{ref}.{kind}_tpm.tsv"), kind, list(cats), samples, ftpm)
        top = np.argsort(-ftpm.mean(axis=1))[: a.top_functions]
        func_json[kind] = [{"id": str(cats[i]), "tpm": [r(x, 2) for x in ftpm[i]]} for i in top]
        log(f"{kind}: {len(cats)} categories")

    # ---- top genes for the viz ----
    top_rows = []
    if len(genes):
        mean_tpm = g_tpm.mean(axis=1)
        top = np.argsort(-mean_tpm)[: a.top_genes]
        bnames = (bins or [])
        for i in top:
            gb = gene_bin[i]
            top_rows.append({
                "id": genes.at[i, "gene_id"], "contig": genes.at[i, "contig"],
                "type": genes.at[i, "type"], "length": int(g_len[i]),
                "product": genes.at[i, "product"], "gene": genes.at[i, "gene"],
                "bin": bnames[gb] if gb >= 0 else None,
                **({k: (str(v[i]) or None) for k, v in funcs.items()}),
                "reads": [int(x) for x in gmat[i]],
                "tpm": [r(x, 1) for x in g_tpm[i]],
            })

    # ---- viz JSONs ----
    vz = os.path.join(a.out, "viz")
    jdump(os.path.join(vz, f"{ref}.samples.json.gz"), {
        "ref": ref, "samples": samples,
        "rows": json.loads(summ.to_json(orient="records")),
        "n_genes": int(len(genes)), "n_contigs": int(len(clen)),
        "genes_detected": [int(x) for x in np.count_nonzero(gmat, axis=0)],
    })
    jdump(os.path.join(vz, f"{ref}.bins.json.gz"), {"ref": ref, "samples": samples, "bins": bin_rows})
    jdump(os.path.join(vz, f"{ref}.functions.json.gz"), {"ref": ref, "samples": samples, "kinds": func_json})
    jdump(os.path.join(vz, f"{ref}.top_genes.json.gz"), {"ref": ref, "samples": samples, "genes": top_rows})
    # ---- contig t-SNE with per-sample expression ----
    emb = load_contig_embedding(a.contig_explorer, contig_index)
    if emb is not None:
        ex, ey = emb
        tot = cmat.sum(axis=1)
        cand = np.flatnonzero((tot > 0) & ~np.isnan(ex))
        if len(cand) > a.max_map_contigs:
            cand = cand[np.argsort(-tot[cand], kind="stable")[: a.max_map_contigs]]
        cand.sort()
        clen_arr = clen["length"].to_numpy()
        bnames = bins or []
        jdump(os.path.join(vz, f"{ref}.contig_map.json.gz"), {
            "ref": ref, "samples": samples,
            "n_contigs_total": int(len(clen)), "n": int(len(cand)),
            "x": [round(float(v), 2) for v in ex[cand]],
            "y": [round(float(v), 2) for v in ey[cand]],
            "length": [int(v) for v in clen_arr[cand]],
            "bins": bnames,
            "bin": [int(contig_bin[i]) if bins else -1 for i in cand],
            # reads on these contigs per sample, for per-million scaling
            "sample_totals": [int(v) for v in cmat[cand].sum(axis=0)],
            "files": [f"{ref}.contig_expr.{j:04d}.json.gz" for j in range(n_s)],
        })
        sub = cmat[cand]
        for j in range(n_s):
            nz = np.flatnonzero(sub[:, j])
            jdump(os.path.join(vz, f"{ref}.contig_expr.{j:04d}.json.gz"), {
                "sample": samples[j],
                "idx": np.diff(nz, prepend=0).astype(int).tolist(),   # delta-encoded
                "reads": sub[nz, j].astype(int).tolist(),
            })
        log(f"contig map: {len(cand)} contigs with RNA reads, {n_s} sample layers")

    log("done")


if __name__ == "__main__":
    main()
