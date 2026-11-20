#!/usr/bin/env python3
"""Build an NCBI-style taxdump (nodes.dmp, names.dmp, merged.dmp) from GTDB
taxonomy tables, plus an accession -> species taxid map, for a BBTools tree.

Usage: build_gtdb_taxdump.py OUTDIR bac120_taxonomy.tsv.gz ar53_taxonomy.tsv.gz

Ranks are written as domain/phylum/class/order/family/genus/species. BBTools
treats "superkingdom" as a separate level from "domain", and danaSeq's
sendsketch parsers keep only the d:/p:/c:/o:/f:/g:/s: lineage tokens, so the
top rank must be "domain".

Taxid 1 is the root, d__Bacteria is 2 and d__Archaea is 2157 (BBTools looks
them up by those NCBI ids). Every other taxon is numbered from 3,000,000 in
sorted order of its prefixed GTDB name, so a rerun on the same release gives
the same ids. Ids are not comparable across GTDB releases.
Names are written without the rank prefix; BBTools prints its own (d:, p:, ...).

accession2taxid.tsv maps both forms of each accession (GCA_x.y and GCF_x.y,
taken from the RS_/GB_ entry) to its species taxid.
"""
import gzip
import os
import sys

RANKS = {'d': 'domain', 'p': 'phylum', 'c': 'class', 'o': 'order',
         'f': 'family', 'g': 'genus', 's': 'species'}


def main():
    outdir, *tables = sys.argv[1:]
    os.makedirs(outdir, exist_ok=True)
    parent = {}       # prefixed name -> prefixed parent name ('' for root)
    genome_species = {}
    for path in tables:
        with gzip.open(path, 'rt') as f:
            for line in f:
                acc, lineage = line.rstrip('\n').split('\t')
                tokens = lineage.split(';')
                assert len(tokens) == 7 and all(t[1:3] == '__' for t in tokens), line
                up = ''
                for t in tokens:
                    prev = parent.setdefault(t, up)
                    assert prev == up, f'{t} has parents {prev} and {up}'
                    up = t
                genome_species[acc[3:]] = tokens[-1]   # drop RS_/GB_

    # BBTools hard-codes NCBI taxids: Bacteria must be 2 and Archaea 2157, and
    # ids such as 2759 (Eukaryota) or 10239 (Viruses) drive its filters. The
    # other GTDB taxa are numbered from 3,000,000 so they hit none of those.
    fixed = {'d__Bacteria': 2, 'd__Archaea': 2157}
    others = sorted(n for n in parent if n not in fixed)
    taxid = dict(fixed, **{n: i for i, n in enumerate(others, start=3_000_000)})
    with open(os.path.join(outdir, 'nodes.dmp'), 'w') as nodes, \
         open(os.path.join(outdir, 'names.dmp'), 'w') as names:
        # NCBI nodes.dmp has 13 columns; BBTools' parser needs at least 4.
        pad = '\t|\t'.join([''] * 10) + '\t|\n'
        nodes.write('1\t|\t1\t|\tno rank\t|\t' + pad)
        names.write('1\t|\troot\t|\t\t|\tscientific name\t|\n')
        for name, tid in taxid.items():
            ptid = taxid[parent[name]] if parent[name] else 1
            nodes.write(f'{tid}\t|\t{ptid}\t|\t{RANKS[name[0]]}\t|\t' + pad)
            names.write(f'{tid}\t|\t{name[3:]}\t|\t\t|\tscientific name\t|\n')
    open(os.path.join(outdir, 'merged.dmp'), 'w').close()

    with open(os.path.join(outdir, 'accession2taxid.tsv'), 'w') as f:
        f.write('accession\ttaxid\tgtdb_species\n')
        for acc, sp in sorted(genome_species.items()):
            base = acc.split('_', 1)[1]          # 000005845.2
            for prefix in ('GCA_', 'GCF_'):
                f.write(f'{prefix}{base}\t{taxid[sp]}\t{sp}\n')
    print(f'{len(taxid)} taxa, {len(genome_species)} genomes', file=sys.stderr)


if __name__ == '__main__':
    main()
