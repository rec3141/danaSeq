#!/usr/bin/env python3
"""
Assemble the viz data dir from the per-reference tables written by
summarize_expression.py:

  <ref>.samples.json.gz    read flow / mapping / assignment per sample
  <ref>.bins.json.gz       per-MAG reads + TPM, CheckM2, GTDB-Tk
  <ref>.functions.json.gz  top KO / COG / CAZy TPM
  <ref>.top_genes.json.gz  most expressed genes with annotation
  <ref>.contig_map.json.gz contig t-SNE (from mag_analysis) + per-sample
                           <ref>.contig_expr.NNNN.json.gz layers

Copies them into --out and writes index.json, which the SPA loads first.
Stdlib only.
"""
import argparse, glob, gzip, json, os, shutil


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--mag-url", default="")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()

    refs = []
    for f in sorted(glob.glob(os.path.join(a.tables, "*.samples.json.gz"))):
        name = os.path.basename(f)[: -len(".samples.json.gz")]
        with gzip.open(f, "rt") as fh:
            s = json.load(fh)
        files = {}
        for src in sorted(glob.glob(os.path.join(a.tables, f"{name}.*.json.gz"))):
            dst = os.path.join(a.out, os.path.basename(src))
            if os.path.abspath(src) != os.path.abspath(dst):
                shutil.copyfile(src, dst)
        for kind in ("samples", "bins", "functions", "top_genes", "contig_map"):
            if os.path.exists(os.path.join(a.tables, f"{name}.{kind}.json.gz")):
                files[kind] = f"{name}.{kind}.json.gz"
        refs.append({
            "name": name,
            "n_samples": len(s["samples"]),
            "n_genes": s.get("n_genes"),
            "n_contigs": s.get("n_contigs"),
            "files": files,
        })

    index = {"title": a.title or (refs[0]["name"] if len(refs) == 1 else "RNA-seq"),
             "mag_url": a.mag_url or None,
             "references": refs}
    with open(os.path.join(a.out, "index.json"), "w") as fh:
        json.dump(index, fh, indent=1)
    # microscape.app serves the .gz copy of each data file
    with gzip.open(os.path.join(a.out, "index.json.gz"), "wt") as fh:
        json.dump(index, fh)
    print(f"[viz] {len(refs)} reference(s): {', '.join(r['name'] for r in refs)}")


if __name__ == "__main__":
    main()
