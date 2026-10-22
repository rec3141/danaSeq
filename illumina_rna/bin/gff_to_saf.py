#!/usr/bin/env python3
"""
Turn a reference FASTA index + GFF3 into featureCounts SAF files and a slim
per-gene annotation table, in one streaming pass (Bakta GFFs of whole
co-assemblies run to several GB).

Outputs (prefix = --out):
  <out>.genes.saf          GeneID Chr Start End Strand   (one row per counted feature)
  <out>.contigs.saf        one row per contig, whole length
  <out>.genes.tsv.gz       gene_id contig start end strand type length product gene dbxref
  <out>.contig_lengths.tsv contig length

Genes are named by --attr (locus_tag) and fall back to ID. Features whose
seqid is not in the FASTA are dropped with a warning (BAM @SQ names come from
the FASTA). Stdlib only, so it runs under any python3.
"""
import argparse, gzip, sys
from urllib.parse import unquote


def read_fai(path):
    lengths = {}
    order = []
    with open(path) as f:
        for line in f:
            name, length = line.split("\t", 2)[:2]
            lengths[name] = int(length)
            order.append(name)
    return order, lengths


def parse_attrs(col):
    out = {}
    for kv in col.strip().split(";"):
        if not kv:
            continue
        k, _, v = kv.partition("=")
        out[k] = unquote(v)
    return out


def add_rrna_mask(path, lengths, saf, ann):
    by_contig = {}
    with open(path) as f:
        header = next(f).rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header) or p[ix["contig_id"]] not in lengths:
                continue
            c = p[ix["contig_id"]]
            by_contig.setdefault(c, []).append((int(p[ix["start"]]), int(p[ix["end"]]), p[ix["rrna_type"]], p[ix["kingdom"]]))
    n = 0
    for c, ivs in by_contig.items():
        ivs.sort()
        merged = []
        for s, e, typ, king in ivs:
            if merged and s <= merged[-1][1]:
                m = merged[-1]
                m[1] = max(m[1], e); m[2].add(typ); m[3].add(king)
            else:
                merged.append([s, e, {typ}, {king}])
        for s, e, typs, kings in merged:
            gid = f"rRNA:{c}:{s}-{e}"
            saf.write(f"{gid}\t{c}\t{s}\t{e}\t+\n")
            product = "/".join(sorted(t.replace("_rRNA", "") for t in typs)) + " rRNA (barrnap)"
            ann.write(f"{gid}\t{c}\t{s}\t{e}\t+\trRNA\t{e - s + 1}\t{product}\t\t\n")
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fai", required=True)
    ap.add_argument("--gff", help="GFF3 (optional; contig SAF only without it)")
    ap.add_argument("--types", default="CDS,rRNA,tRNA,tmRNA,ncRNA")
    ap.add_argument("--attr", default="locus_tag")
    ap.add_argument("--rrna-tsv", help="mag_analysis taxonomy/rrna/rrna_genes.tsv (barrnap calls)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    order, lengths = read_fai(a.fai)
    with open(f"{a.out}.contigs.saf", "w") as saf:
        saf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for c in order:
            saf.write(f"{c}\t{c}\t1\t{lengths[c]}\t+\n")
    with open(f"{a.out}.contig_lengths.tsv", "w") as f:
        f.write("contig\tlength\n")
        for c in order:
            f.write(f"{c}\t{lengths[c]}\n")

    types = set(t.strip() for t in a.types.split(",") if t.strip())
    n_genes = n_missing = 0
    seen = set()
    with open(f"{a.out}.genes.saf", "w") as saf, \
         gzip.open(f"{a.out}.genes.tsv.gz", "wt", compresslevel=3) as ann:
        saf.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        ann.write("gene_id\tcontig\tstart\tend\tstrand\ttype\tlength\tproduct\tgene\tdbxref\n")
        if a.gff:
            opener = gzip.open if a.gff.endswith(".gz") else open
            with opener(a.gff, "rt") as f:
                for line in f:
                    if line.startswith("##FASTA"):
                        break
                    if line.startswith("#") or not line.strip():
                        continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 9 or p[2] not in types:
                        continue
                    seqid = p[0].split()[0]
                    if seqid not in lengths:
                        n_missing += 1
                        continue
                    at = parse_attrs(p[8])
                    gid = at.get(a.attr) or at.get("ID")
                    if not gid or gid in seen:
                        continue
                    seen.add(gid)
                    start, end = int(p[3]), int(p[4])
                    strand = p[6] if p[6] in "+-" else "+"
                    saf.write(f"{gid}\t{seqid}\t{start}\t{end}\t{strand}\n")
                    product = at.get("product", "").replace("\t", " ")
                    ann.write(f"{gid}\t{seqid}\t{start}\t{end}\t{strand}\t{p[2]}\t{end - start + 1}\t"
                              f"{product}\t{at.get('gene', '')}\t{at.get('Dbxref', '')}\n")
                    n_genes += 1

        # rRNA genes called by barrnap in the MAG analysis (all kingdoms). Bakta
        # annotates bacterial/archaeal rRNA but not eukaryotic or divergent
        # rRNA, whose reads would otherwise land on spurious ORFs. Overlapping
        # calls are merged into one feature per locus; featureCounts'
        # --largestOverlap then gives those reads to the rRNA, not the ORF.
        if a.rrna_tsv:
            n_mask = add_rrna_mask(a.rrna_tsv, lengths, saf, ann)
            print(f"[gff_to_saf] {n_mask} rRNA loci from {a.rrna_tsv}", file=sys.stderr)

    print(f"[gff_to_saf] {len(order)} contigs, {n_genes} features ({','.join(sorted(types))})",
          file=sys.stderr)
    if n_missing:
        print(f"[gff_to_saf] WARNING: {n_missing} features on seqids absent from the FASTA were dropped",
              file=sys.stderr)


if __name__ == "__main__":
    main()
