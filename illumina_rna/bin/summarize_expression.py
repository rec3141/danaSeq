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


def tpm(counts, lengths_bp, norm_rows=None):
    """TPM with the per-million denominator taken over norm_rows (protein-coding
    genes): rRNA and other RNA genes get values on the same scale but do not
    dilute the CDS."""
    # float32: a co-assembly matrix is ~10M genes x 100+ samples
    rpk = counts.astype(np.float32) / (np.maximum(lengths_bp, 1)[:, None] / 1e3).astype(np.float32)
    # masked dot product: no copy of the selected rows
    tot = (norm_rows.astype(np.float32) @ rpk) if norm_rows is not None else rpk.sum(axis=0)
    rpk_sum = tot.astype(np.float64)
    tot[tot == 0] = 1
    rpk /= tot
    rpk *= 1e6
    return rpk, rpk_sum


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
        if kind == "cazy":
            # dbCAN strings mix families, domain ranges, subfamilies and
            # GenBank hits ("AAZ20843.1,CE11,CE11(6-276),CE11_e24(7-277)");
            # keep the first real family. Family 0 means unclassified.
            s = s.fillna("").str.extract(r"\b((?:GH|GT|PL|CE|AA|CBM)(?!0(?!\d))\d+)(?!\d)", expand=False).fillna("")
            out[kind] = s.to_numpy()
            continue
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
    # Taxon labels ride along for search and colour: coded columns index the
    # manifest vocabulary (or hold strings in older runs).
    vocab = manifest.get("vocab") or {}
    if isinstance(vocab, str):
        vocab = {}
    tax_cols = tuple(f"kaiju_{r}" for r in ("domain", "phylum", "class", "order", "family", "genus")) + \
               ("rrna_domain", "rrna_phylum", "rrna_genus")   # SILVA class of the contig's rRNA genes
    tax = {k: np.full(len(contig_index), "", dtype=object) for k in tax_cols}
    hit = 0
    for f in parts:
        with op(f) as fh:
            c = json.load(fh)
        nums = np.cumsum(np.asarray(c["ids"], dtype=np.int64))
        cols = c["columns"]
        tx, ty = cols["tsne_x"], cols["tsne_y"]
        for k, n in enumerate(nums):
            i = contig_index.get(f"{enc['prefix']}{int(n):0{enc['width']}d}")
            if i is None or tx[k] is None or ty[k] is None:
                continue
            x[i], y[i] = tx[k], ty[k]
            hit += 1
            for col in tax_cols:
                v = cols.get(col, [None] * len(nums))[k]
                if isinstance(v, int) and col in vocab:
                    v = vocab[col][v] if 0 <= v < len(vocab[col]) else ""
                tax[col][i] = v or ""
    log(f"contig t-SNE: {hit} of {len(contig_index)} contigs placed from {len(parts)} parts")
    return (x, y, tax) if hit else None


RANKS = ["domain", "phylum", "class", "order", "family", "genus"]
# taxon columns carried into the contig map (colour, search, contig panel)
MAP_TAX = ("kaiju_phylum", "kaiju_genus", "rrna_domain", "rrna_phylum", "rrna_genus")
GTDB_PREFIX = {"d": "domain", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus"}


def rank_labels(rank_cols):
    """Per item and rank, the lineage down to that rank; an empty rank keeps
    its parent's lineage plus ';unclassified'. rank_cols: 6 arrays of names."""
    n = len(rank_cols[0])
    out = []
    path = np.array([""] * n, dtype=object)   # classified ranks only
    for col in rank_cols:
        col = np.asarray(col, dtype=object)
        has = col != ""
        path = np.where(has, np.where(path == "", col, path + ";" + col), path)
        out.append(np.where(has, path, np.where(path == "", "unclassified", path + ";unclassified")))
    return out


def taxon_matrix(labels, mat, top):
    """Sum rows of mat by label; return (all_labels, sums) and the top rows by mean share."""
    cats, idx = np.unique(labels.astype(str), return_inverse=True)
    sums = aggregate(idx, len(cats), mat)
    tot = sums.sum(axis=0)
    tot[tot == 0] = 1
    share = (sums / tot).mean(axis=1)
    order = np.argsort(-share)[:top]
    return cats, sums, order


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
    ap.add_argument("--rrna-tax-dir", help="per-sample <sample>.rrna_taxonomy.tsv from classify_rrna.py")
    ap.add_argument("--top-taxa", type=int, default=300)
    ap.add_argument("--gene-shard-contigs", type=int, default=2000)
    ap.add_argument("--gene-annot", required=True)
    ap.add_argument("--contig-lengths", required=True)
    ap.add_argument("--bins")
    ap.add_argument("--checkm2")
    ap.add_argument("--gtdbtk")
    ap.add_argument("--annotations")
    ap.add_argument("--contig-explorer", help="mag_analysis viz/data dir with contig_explorer parts")
    ap.add_argument("--max-map-contigs", type=int, default=400000)
    ap.add_argument("--ko-names", help="KO<TAB>definition")
    ap.add_argument("--top-genes", type=int, default=2000)
    ap.add_argument("--top-functions", type=int, default=2000)
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

    gtype = genes["type"].to_numpy() if len(genes) else np.array([], dtype=object)
    is_cds = gtype == "CDS"
    is_rrna = gtype == "rRNA"
    cds_pairs = gmat[is_cds].sum(axis=0) if len(genes) else np.zeros(n_s)
    rrna_feat_pairs = gmat[is_rrna].sum(axis=0) if len(genes) else np.zeros(n_s)

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
        j = samples.index(s)
        # pairs on protein-coding genes vs on rRNA loci the rRNA filter missed
        row["pairs_cds"] = int(cds_pairs[j])
        row["pairs_rrna_features"] = int(rrna_feat_pairs[j])
        row["pct_cds"] = (100 * 2 * cds_pairs[j] / final) if final else None
        if st.get("filtered") and "norrna" in st:
            row["pct_rrna_total"] = 100 * (st["filtered"] - st["norrna"] + 2 * rrna_feat_pairs[j]) / st["filtered"]
        rows.append(row)
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(a.out, f"{ref}.sample_summary.tsv"), sep="\t", index=False, float_format="%.3f")

    # ---- TPM + MAG aggregation ----
    g_len = genes["length"].to_numpy(dtype=np.float64)
    if len(genes):
        g_tpm, cds_rpk_sum = tpm(gmat, g_len, is_cds)
    else:
        g_tpm, cds_rpk_sum = np.zeros((0, n_s), dtype=np.float32), np.zeros(n_s)
    g_contig = genes["contig"].map(contig_index).fillna(-1).to_numpy(dtype=np.int64)
    # Read pairs on protein-coding genes per contig: contig, MAG and taxon
    # expression use these, so rRNA (filtered or not) never counts as activity
    ccds = aggregate(np.where(is_cds, g_contig, -1), len(clen), gmat) if len(genes) else cmat.astype(np.float64)
    rrna_contig = np.zeros(len(clen), dtype=bool)
    rrna_contig[g_contig[is_rrna & (g_contig >= 0)]] = True

    bins, contig_bin = load_bins(a.bins, contig_index)
    bin_rows = []
    if bins:
        n_b = len(bins)
        bin_counts = aggregate(contig_bin, n_b, ccds)
        gene_bin = np.where(g_contig >= 0, contig_bin[np.maximum(g_contig, 0)], -1)
        bin_tpm = aggregate(np.where(is_cds, gene_bin, -1), n_b, g_tpm)
        # unbinned share, so the MAG table accounts for every CDS read
        unb_counts = ccds[contig_bin < 0].sum(axis=0)
        unb_tpm = g_tpm[is_cds & (gene_bin < 0)].sum(axis=0)
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
    func_json, func_totals = {}, {}
    ko_names = {}
    if a.ko_names and os.path.exists(a.ko_names):
        with open(a.ko_names) as fh:
            for line in fh:
                k, _, d = line.rstrip("\n").partition("\t")
                ko_names[k] = d
        log(f"{len(ko_names)} KO definitions")
    for kind, labels in funcs.items():
        labels = np.where(is_cds, labels, "")
        cats, idx = np.unique(labels, return_inverse=True)
        idx = np.where(cats[idx] == "", -1, idx)
        ftpm = aggregate(idx, len(cats), g_tpm)
        keep = cats != ""
        cats, ftpm = cats[keep], ftpm[keep]
        write_table(os.path.join(a.out, f"{ref}.{kind}_tpm.tsv"), kind, list(cats), samples, ftpm)
        top = np.argsort(-ftpm.mean(axis=1))[: a.top_functions]
        named = kind == "ko" and ko_names
        func_json[kind] = [{"id": str(cats[i]), **({"name": ko_names.get(str(cats[i]))} if named else {}),
                            "tpm": [r(x, 2) for x in ftpm[i]]} for i in top]
        func_totals[kind] = int(len(cats))
        log(f"{kind}: {len(cats)} categories")

    # ---- top genes for the viz ----
    top_rows = []
    if len(genes):
        mean_tpm = np.where(is_cds, g_tpm.mean(axis=1), -1)
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
    jdump(os.path.join(vz, f"{ref}.functions.json.gz"), {"ref": ref, "samples": samples, "kinds": func_json,
                                                           "n_total": func_totals})
    jdump(os.path.join(vz, f"{ref}.top_genes.json.gz"), {"ref": ref, "samples": samples, "genes": top_rows})
    # ---- contig t-SNE with per-sample expression ----
    emb = load_contig_embedding(a.contig_explorer, contig_index)
    vz = os.path.join(a.out, "viz")
    if emb is not None:
        ex, ey, etax = emb
        tot = ccds.sum(axis=1)
        cand = np.flatnonzero((tot > 0) & ~np.isnan(ex))
        if len(cand) > a.max_map_contigs:
            cand = cand[np.argsort(-tot[cand], kind="stable")[: a.max_map_contigs]]
        cand.sort()
        clen_arr = clen["length"].to_numpy()
        bnames = bins or []
        sub = np.rint(ccds[cand]).astype(np.int64)
        kb = clen_arr[cand] / 1e3
        rpkm_p99 = []
        for j in range(n_s):
            nz = sub[:, j] > 0
            tot_j = max(int(cds_pairs[j]), 1)
            rpkm_p99.append(round(float(np.percentile(sub[nz, j] / kb[nz] / (tot_j / 1e6), 99)), 2) if nz.any() else 0.0)
        jdump(os.path.join(vz, f"{ref}.contig_map.json.gz"), {
            "ref": ref, "samples": samples,
            "n_contigs_total": int(len(clen)), "n": int(len(cand)),
            "x": [round(float(v), 2) for v in ex[cand]],
            "y": [round(float(v), 2) for v in ey[cand]],
            "length": [int(v) for v in clen_arr[cand]],
            "bins": bnames,
            "bin": [int(contig_bin[i]) if bins else -1 for i in cand],
            "id": [str(c) for c in clen["contig"].to_numpy()[cand]],
            "total_reads": [int(v) for v in tot[cand]],
            "n_detected": [int(v) for v in np.count_nonzero(sub, axis=1)],
            # contigs carrying rRNA loci (Bakta or barrnap); their CDS reads are
            # often spurious ORFs inside rRNA
            "rrna": [int(v) for v in rrna_contig[cand]],
            # taxon codes into tax_vocab (empty string = unclassified)
            **{f"{col}_vocab": sorted(set(etax[col][cand])) for col in MAP_TAX},
            **{col: (lambda voc: [voc[t] for t in etax[col][cand]])({t: j for j, t in enumerate(sorted(set(etax[col][cand])))})
               for col in MAP_TAX},
            # reads on these contigs per sample, for per-million scaling
            # per-million basis: all CDS pairs in the sample
            "sample_totals": [int(v) for v in cds_pairs],
            # per-sample sum of CDS reads per kb, for gene TPM in the browser
            "rpk_sum": [round(float(v), 3) for v in cds_rpk_sum],
            # 99th percentile RPKM of each sample's detected contigs: one
            # colour ceiling that holds across samples
            "rpkm_p99": rpkm_p99,
            "files": [f"{ref}.contig_expr.{j:04d}.json.gz" for j in range(n_s)],
            "gene_shard_size": a.gene_shard_contigs,
            "gene_shards": [f"{ref}.contig_genes.part-{k:03d}.json.gz"
                            for k in range(-(-len(cand) // a.gene_shard_contigs))] if len(genes) else [],
            "gene_search": f"{ref}.gene_search.json.gz" if len(genes) else None,
        })
        for j in range(n_s):
            nz = np.flatnonzero(sub[:, j])
            jdump(os.path.join(vz, f"{ref}.contig_expr.{j:04d}.json.gz"), {
                "sample": samples[j],
                "idx": np.diff(nz, prepend=0).astype(int).tolist(),   # delta-encoded
                "reads": sub[nz, j].astype(int).tolist(),
            })
        log(f"contig map: {len(cand)} contigs with RNA reads, {n_s} sample layers")

        # ---- genes of the map contigs: shards for the contig panel, search index ----
        if len(genes):
            S = a.gene_shard_contigs
            mappos = np.full(len(clen), -1, dtype=np.int64)
            mappos[cand] = np.arange(len(cand))
            gp = np.where(g_contig >= 0, mappos[np.maximum(g_contig, 0)], -1)
            sel = np.flatnonzero(gp >= 0)
            gstart = genes["start"].to_numpy()
            order = sel[np.lexsort((gstart[sel], gp[sel]))]
            gp_o = gp[order]
            nsh = -(-len(cand) // S)
            bounds = np.searchsorted(gp_o, np.arange(0, nsh * S + 1, S))
            gid = genes["gene_id"].to_numpy(dtype=object)
            gend = genes["end"].to_numpy()
            gstrand = genes["strand"].to_numpy(dtype=object)
            gname = genes["gene"].to_numpy(dtype=object)
            gprod = genes["product"].to_numpy(dtype=object)
            gko = funcs.get("ko")
            for k in range(nsh):
                R = order[bounds[k]:bounds[k + 1]]
                first = k * S
                ncont = min(S, len(cand) - first)
                offs = np.searchsorted(gp_o[bounds[k]:bounds[k + 1]], np.arange(first, first + ncont + 1))
                pv, pc = np.unique(gprod[R].astype(str), return_inverse=True)
                sub_g = gmat[R]
                rr, cc = np.nonzero(sub_g)
                co = np.concatenate([[0], np.cumsum(np.bincount(rr, minlength=len(R)))])
                jdump(os.path.join(vz, f"{ref}.contig_genes.part-{k:03d}.json.gz"), {
                    "first": first, "n": int(ncont), "offsets": offs.tolist(),
                    "id": gid[R].tolist(), "s": gstart[R].tolist(), "e": gend[R].tolist(),
                    "d": [1 if x == "+" else -1 for x in gstrand[R]],
                    "t": [str(x).lower() for x in gtype[R]],
                    "g": gname[R].tolist(), "pv": pv.tolist(), "p": pc.tolist(),
                    **({"ko": [str(x) for x in gko[R]]} if gko is not None else {}),
                    # sparse read pairs per gene: samples ci[co[i]:co[i+1]], counts cv[...]
                    "ci": cc.tolist(), "cv": sub_g[rr, cc].tolist(), "co": co.tolist(),
                })
            # search terms per map contig: products (not hypothetical), gene names, KOs
            terms = []
            for arr in (gprod, gname, gko if gko is not None else None):
                if arr is None:
                    continue
                t = pd.Series(arr[order].astype(str)).str.strip()
                terms.append(pd.DataFrame({"pos": gp_o, "term": t}))
            if gko is not None and ko_names:
                kn = pd.Series(gko[order].astype(str)).map(lambda k: f"{k} {ko_names[k]}" if k in ko_names else "")
                terms.append(pd.DataFrame({"pos": gp_o, "term": kn}))
            tdf = pd.concat(terms, ignore_index=True)
            tdf = tdf[(tdf["term"] != "") & (tdf["term"] != "hypothetical protein") & (tdf["term"] != "nan")]
            tdf = tdf.drop_duplicates()
            codes, vocab = pd.factorize(tdf["term"])
            tdf = tdf.assign(code=codes).sort_values(["pos", "code"])
            toff = np.concatenate([[0], np.cumsum(np.bincount(tdf["pos"].to_numpy(), minlength=len(cand)))])
            jdump(os.path.join(vz, f"{ref}.gene_search.json.gz"), {
                "vocab": [str(v) for v in vocab], "offsets": toff.tolist(), "codes": tdf["code"].tolist(),
            })
            log(f"gene shards: {len(order)} genes in {nsh} shards; search index {len(vocab)} terms")

    # ---- taxonomy: rRNA reads (SILVA) and mRNA activity (contig Kaiju, MAG GTDB) ----
    sources = {}
    if a.rrna_tax_dir and os.path.isdir(a.rrna_tax_dir):
        per = {}
        totals = {}
        for f in glob.glob(os.path.join(a.rrna_tax_dir, "*.rrna_taxonomy.tsv")):
            smp = os.path.basename(f)[: -len(".rrna_taxonomy.tsv")]
            with open(f) as fh:
                for line in fh:
                    if line.startswith("# reads"):
                        totals[smp] = int(line.split("\t")[1])
                    if line.startswith("#") or line.startswith("rank\t"):
                        continue
                    rk, lin, c = line.rstrip("\n").split("\t")
                    per.setdefault(rk, {}).setdefault(lin, {})[smp] = int(c)
        if per:
            ranks_json = {}
            long_rows = []
            for rk in RANKS:
                labs = sorted(per.get(rk, {}))
                m = np.array([[per[rk][l].get(smp, 0) for smp in samples] for l in labs], dtype=np.float64) \
                    if labs else np.zeros((0, n_s))
                tot = np.array([max(totals.get(smp, 0), 1) for smp in samples], dtype=np.float64)
                order = np.argsort(-(m / tot).mean(axis=1))[: a.top_taxa] if len(labs) else []
                ranks_json[rk] = [{"taxon": labs[i], "counts": [int(x) for x in m[i]]} for i in order]
                long_rows += [[rk, labs[i], *[int(x) for x in m[i]]] for i in range(len(labs))]
            sources["rrna"] = {"label": "rRNA reads (SILVA)", "unit": "reads",
                               "totals": [totals.get(smp, 0) for smp in samples], "ranks": ranks_json}
            pd.DataFrame(long_rows, columns=["rank", "lineage", *samples]).to_csv(
                os.path.join(a.out, f"{ref}.taxonomy_rrna.tsv"), sep="\t", index=False)
            log(f"rRNA taxonomy: {len(totals)} samples")

    if emb is not None:
        labs_by_rank = rank_labels([etax[f"kaiju_{r}"] for r in RANKS])
        ranks_json, long_rows = {}, []
        for rk, labs in zip(RANKS, labs_by_rank):
            cats, sums, order = taxon_matrix(labs, ccds, a.top_taxa)
            ranks_json[rk] = [{"taxon": str(cats[i]), "counts": [int(round(x)) for x in sums[i]]} for i in order]
            long_rows += [[rk, str(cats[i]), *[int(round(x)) for x in sums[i]]] for i in range(len(cats))]
        sources["mrna_kaiju"] = {"label": "mRNA by contig (Kaiju)", "unit": "CDS read pairs",
                                 "totals": [int(x) for x in ccds.sum(axis=0)], "ranks": ranks_json}
        pd.DataFrame(long_rows, columns=["rank", "lineage", *samples]).to_csv(
            os.path.join(a.out, f"{ref}.taxonomy_mrna_kaiju.tsv"), sep="\t", index=False)

    if bins and any(br.get("classification") for br in bin_rows):
        cls = [br.get("classification") or "" for br in bin_rows if br["bin"] != "unbinned"]
        cols = [[] for _ in RANKS]
        for c in cls:
            d = {GTDB_PREFIX.get(p.split("__")[0]): p.split("__")[1] for p in c.split(";") if "__" in p}
            for i, rk in enumerate(RANKS):
                cols[i].append(d.get(rk, "") or "")
        labs_by_rank = rank_labels(cols)
        bmat = np.vstack([bin_counts, unb_counts])
        ranks_json, long_rows = {}, []
        for rk, labs in zip(RANKS, labs_by_rank):
            labs = np.append(labs, "unbinned")
            cats, sums, order = taxon_matrix(labs, bmat, a.top_taxa)
            ranks_json[rk] = [{"taxon": str(cats[i]), "counts": [int(round(x)) for x in sums[i]]} for i in order]
            long_rows += [[rk, str(cats[i]), *[int(round(x)) for x in sums[i]]] for i in range(len(cats))]
        sources["mrna_gtdb"] = {"label": "mRNA by MAG (GTDB-Tk)", "unit": "CDS read pairs",
                                "totals": [int(x) for x in bmat.sum(axis=0)], "ranks": ranks_json}
        pd.DataFrame(long_rows, columns=["rank", "lineage", *samples]).to_csv(
            os.path.join(a.out, f"{ref}.taxonomy_mrna_gtdb.tsv"), sep="\t", index=False)

    if sources:
        jdump(os.path.join(vz, f"{ref}.taxonomy.json.gz"), {"ref": ref, "samples": samples, "sources": sources})
        log(f"taxonomy sources: {', '.join(sources)}")

    log("done")


if __name__ == "__main__":
    main()
