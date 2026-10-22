#!/usr/bin/env python3
"""
Classify rRNA reads from their best SILVA hit (SAM on stdin) and count reads
per taxon at each rank.

A read is assigned down to the deepest rank whose identity cutoff its
alignment meets (Yarza et al. 2014, Nat Rev Microbiol: genus 94.5%, family
86.5%, order 82.0%, class 78.5%, phylum 75.0%); below that it stays at the
domain. Reads that do not align over --min-aligned of their length count as
unclassified.

Output TSV: rank  lineage  reads   (lineage is ';'-joined down to that rank)
Stdlib only.
"""
import argparse, re, sys
from collections import Counter

RANKS = ["domain", "phylum", "class", "order", "family", "genus"]
CUTOFF = {"phylum": 75.0, "class": 78.5, "order": 82.0, "family": 86.5, "genus": 94.5}
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def load_taxonomy(path):
    tax = {}
    with open(path) as f:
        for line in f:
            acc, _, lin = line.rstrip("\n").partition("\t")
            tax[acc] = lin.split(";")
    return tax


def ranks_of(lineage):
    """SILVA lineages are 6-7 levels for prokaryotes; eukaryote lineages are
    deeper and uneven, so map them to domain, then the next four levels, then
    the genus (second last, as SILVA ends in the organism name)."""
    if not lineage:
        return []
    if lineage[0] in ("Bacteria", "Archaea"):
        return lineage[:6]
    if len(lineage) <= 6:
        return lineage[:6]
    return lineage[:5] + [lineage[-2]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxonomy", required=True, help="accession<TAB>lineage")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--min-aligned", type=float, default=0.8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    tax = load_taxonomy(a.taxonomy)
    counts = {r: Counter() for r in RANKS}
    n = unclassified = 0
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 0x900:
            continue
        n += 1
        if flag & 4:
            unclassified += 1
            continue
        aligned = sum(int(k) for k, op in CIGAR.findall(f[5]) if op in "MI=X")
        if aligned < a.min_aligned * len(f[9]):
            unclassified += 1
            continue
        nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), 0)
        ident = 100.0 * (1 - nm / aligned)
        lin = ranks_of(tax.get(f[2].split()[0], []))
        if not lin:
            unclassified += 1
            continue
        for i, rank in enumerate(RANKS):
            if i >= len(lin) or (rank in CUTOFF and ident < CUTOFF[rank]):
                # unassigned at this rank and below: keep the parent's name
                label = ";".join(lin[:i]) + ";unclassified" if i else "unclassified"
                for r in RANKS[i:]:
                    counts[r][label] += 1
                break
            counts[rank][";".join(lin[: i + 1])] += 1

    with open(a.out, "w") as out:
        out.write(f"# sample\t{a.sample}\n# reads\t{n}\n# unclassified\t{unclassified}\n")
        out.write("rank\tlineage\treads\n")
        for rank in RANKS:
            for lin, c in counts[rank].most_common():
                out.write(f"{rank}\t{lin}\t{c}\n")
    print(f"[classify_rrna] {a.sample}: {n} reads, {unclassified} unclassified", file=sys.stderr)


if __name__ == "__main__":
    main()
