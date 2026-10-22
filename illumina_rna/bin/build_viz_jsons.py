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
It also writes samples.json ([{id, name}] per sample), which microscape.app
matches metadata uploaded in its run editor against. With --metadata, a
sample metadata TSV/CSV is shipped as metadata.json in the shape that editor
serves: {"<sample>": {"<column>": "<value>"}} (the SPA types the columns).
Stdlib only.
"""
import argparse, csv, glob, gzip, json, os, re, shutil


ID_COLUMNS = {"sample", "sample_id", "sample_name", "samp_name", "sampleid", "name", "id", "library", "library_id"}


def loose(s):
    s = s.strip().lower()
    s = re.sub(r"\.(fastq|fq)(\.gz)?$", "", s)
    return re.sub(r"[_.-](r?[12])$", "", s)


def read_metadata(path, names):
    """Rows of a TSV/CSV keyed by the sample they name (exact, then loose)."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        text = fh.read()
    delim = "\t" if "\t" in text.split("\n", 1)[0] else ","
    rows = [r for r in csv.reader(text.splitlines(), delimiter=delim) if any(v.strip() for v in r)]
    if len(rows) < 2:
        raise SystemExit(f"[viz] metadata {path}: needs a header and at least one row")
    header = [h.strip() for h in rows[0]]
    idc = next((j for j, h in enumerate(header) if h.lower().replace(" ", "_") in ID_COLUMNS), 0)
    exact = set(names)
    by_loose = {loose(n): n for n in names}
    out, unmatched = {}, []
    for r in rows[1:]:
        sid = r[idc].strip() if idc < len(r) else ""
        if not sid:
            continue
        name = sid if sid in exact else by_loose.get(loose(sid))
        if not name:
            unmatched.append(sid)
            continue
        out.setdefault(name, {h: r[j].strip() for j, h in enumerate(header)
                              if j != idc and h and j < len(r) and r[j].strip()})
    if unmatched:
        print(f"[viz] metadata: {len(unmatched)} rows name no sample: {', '.join(unmatched[:5])}")
    if not out:
        raise SystemExit(f"[viz] metadata {path}: no row names a sample (e.g. {names[0] if names else '?'})")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--mag-url", default="")
    ap.add_argument("--metadata", default="")
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
        for kind in ("samples", "bins", "functions", "top_genes", "contig_map", "taxonomy"):
            if os.path.exists(os.path.join(a.tables, f"{name}.{kind}.json.gz")):
                files[kind] = f"{name}.{kind}.json.gz"
        refs.append({
            "name": name,
            "n_samples": len(s["samples"]),
            "n_genes": s.get("n_genes"),
            "n_contigs": s.get("n_contigs"),
            "files": files,
            "_samples": s["samples"],
        })

    index = {"title": a.title or (refs[0]["name"] if len(refs) == 1 else "RNA-seq"),
             "mag_url": a.mag_url or None,
             "references": refs}
    # Sample list for microscape.app's run-editor metadata upload
    names = list(dict.fromkeys(n for r in refs for n in r.pop("_samples")))
    with gzip.open(os.path.join(a.out, "samples.json.gz"), "wt") as fh:
        json.dump([{"id": n, "name": n} for n in names], fh)
    if a.metadata:
        meta = read_metadata(a.metadata, names)
        with gzip.open(os.path.join(a.out, "metadata.json.gz"), "wt") as fh:
            json.dump(meta, fh)
        print(f"[viz] metadata: {len(meta)} of {len(names)} samples")
    with open(os.path.join(a.out, "index.json"), "w") as fh:
        json.dump(index, fh, indent=1)
    # microscape.app serves the .gz copy of each data file
    with gzip.open(os.path.join(a.out, "index.json.gz"), "wt") as fh:
        json.dump(index, fh)
    print(f"[viz] {len(refs)} reference(s): {', '.join(r['name'] for r in refs)}")


if __name__ == "__main__":
    main()
