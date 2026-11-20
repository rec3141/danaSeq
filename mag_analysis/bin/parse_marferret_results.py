#!/usr/bin/env python3
"""Join DIAMOND blastp hits against MarFERReT with taxonomy and Pfam annotations.

Reads:
  - DIAMOND tabular output (blast format 6: qseqid sseqid pident length
    evalue bitscore stitle)
  - MarFERReT taxonomies file (*.taxonomies.tab.gz): NCBI accession2taxid
    layout, columns accession, accession.version, taxid, gi. The MarFERReT
    protein ID (mft...) is in accession.version; accession is always NA.
  - MarFERReT Pfam annotations file (*.best_pfam_annotations.csv.gz): columns
    aa_id, pfam_name, pfam_id, pfam_eval, pfam_score
  - NCBI taxdump directory (nodes.dmp, names.dmp, optional merged.dmp), which
    turns taxids into lineages. Without it only taxon_id is filled.

Writes:
  - marferret_proteins.tsv  Per-protein: protein_id, contig_id, best_hit,
                            pident, ..., taxon_id, taxonomy, pfam
  - marferret_contigs.tsv   Per-contig: contig_id, n_proteins, n_classified,
                            top_taxonomy, pfam_domains

The taxonomy column is "domain;phylum;class;order;family;genus;species", with
an empty field for each rank the NCBI lineage lacks (the viz reads field 2 as
the phylum).

The reference tables hold ~34M rows; only rows for subjects that DIAMOND hit
are kept in memory.
"""

import sys
import os
import csv
import gzip
import argparse
from collections import defaultdict, Counter

RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
# NCBI names the top rank "superkingdom" in older taxdumps, "domain" in newer.
RANK_ALIASES = {'superkingdom': 'domain'}


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def load_taxids(path, wanted):
    """MarFERReT protein ID -> NCBI taxid, for the IDs in `wanted`."""
    taxids = {}
    if not path or not os.path.isfile(path):
        return taxids
    with _open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        id_col = header.index('accession.version') if 'accession.version' in header else 1
        tax_col = header.index('taxid') if 'taxid' in header else 2
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) > max(id_col, tax_col) and parts[id_col] in wanted:
                taxids[parts[id_col]] = parts[tax_col]
    print(f'[INFO] MarFERReT: taxids found for {len(taxids)}/{len(wanted)} hit subjects',
          file=sys.stderr)
    return taxids


def load_pfam(path, wanted):
    """MarFERReT protein ID -> Pfam accession, for the IDs in `wanted`."""
    pfam = {}
    if not path or not os.path.isfile(path):
        return pfam
    with _open(path) as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        id_col = next((c for c in ['aa_id', 'ref_name', 'seq_id'] if c in headers),
                      headers[0] if headers else None)
        pfam_col = next((c for c in ['pfam_id', 'pfam_accession', 'pfam'] if c in headers), None)
        if not id_col or not pfam_col:
            print('[WARNING] MarFERReT: could not identify columns in Pfam file', file=sys.stderr)
            return pfam
        for row in reader:
            seq_id = row.get(id_col, '')
            val = row.get(pfam_col, '')
            if seq_id in wanted and val and val not in ('nan', 'NA'):
                pfam[seq_id] = val
    print(f'[INFO] MarFERReT: Pfam annotations found for {len(pfam)} hit subjects',
          file=sys.stderr)
    return pfam


class Taxdump:
    """Taxid -> ranked lineage from NCBI nodes.dmp / names.dmp / merged.dmp."""

    def __init__(self, taxdump_dir):
        self.parent, self.rank, self.name, self.merged = {}, {}, {}, {}
        with open(os.path.join(taxdump_dir, 'nodes.dmp')) as f:
            for line in f:
                p = line.split('\t|\t')
                self.parent[p[0]] = p[1]
                self.rank[p[0]] = RANK_ALIASES.get(p[2], p[2])
        with open(os.path.join(taxdump_dir, 'names.dmp')) as f:
            for line in f:
                p = line.split('\t|\t')
                if p[3].startswith('scientific name'):
                    self.name[p[0]] = p[1]
        merged_path = os.path.join(taxdump_dir, 'merged.dmp')
        if os.path.isfile(merged_path):
            with open(merged_path) as f:
                for line in f:
                    p = line.split('\t|\t')
                    self.merged[p[0]] = p[1].rstrip('\t|\n')
        self._cache = {}
        print(f'[INFO] MarFERReT: loaded taxdump with {len(self.parent)} nodes',
              file=sys.stderr)

    def lineage(self, taxid):
        if taxid in self._cache:
            return self._cache[taxid]
        tid = self.merged.get(taxid, taxid)
        by_rank = {}
        seen = set()
        while tid in self.parent and tid not in seen:
            seen.add(tid)
            r = self.rank.get(tid)
            if r in RANKS and r not in by_rank:
                by_rank[r] = self.name.get(tid, '')
            nxt = self.parent[tid]
            if nxt == tid:
                break
            tid = nxt
        lin = ';'.join(by_rank.get(r, '') for r in RANKS) if by_rank else ''
        self._cache[taxid] = lin
        return lin


def extract_contig_id(protein_id):
    """Source contig of a MetaEuk protein.

    MetaEuk headers are TargetID|ContigID|strand|..., so the contig is the
    second field; the first is the MetaEuk reference protein that was matched.
    """
    parts = protein_id.split('|')
    return parts[1] if len(parts) > 1 else parts[0]


def read_diamond(path):
    hits = []
    if not (os.path.isfile(path) and os.path.getsize(path) > 0):
        return hits
    with open(path) as f:
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 6:
                continue
            hits.append(cols + [''] * (7 - len(cols)))
    return hits


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--diamond', required=True, help='DIAMOND blast6 output file')
    parser.add_argument('--taxonomy', default=None, help='MarFERReT taxonomies.tab.gz file')
    parser.add_argument('--pfam', default=None, help='MarFERReT best_pfam_annotations.csv.gz file')
    parser.add_argument('--taxdump', default=None,
                        help='NCBI taxdump directory (nodes.dmp, names.dmp, merged.dmp)')
    parser.add_argument('--out-proteins', default='marferret_proteins.tsv')
    parser.add_argument('--out-contigs', default='marferret_contigs.tsv')
    args = parser.parse_args()

    hits = read_diamond(args.diamond)
    wanted = {h[1] for h in hits}

    taxids = load_taxids(args.taxonomy, wanted)
    pfam_map = load_pfam(args.pfam, wanted)
    taxdump = None
    if args.taxdump and os.path.isfile(os.path.join(args.taxdump, 'nodes.dmp')):
        taxdump = Taxdump(args.taxdump)
    elif taxids:
        print('[WARNING] MarFERReT: no NCBI taxdump given; taxonomy column left empty',
              file=sys.stderr)

    proteins = []
    for qseqid, sseqid, pident, length, evalue, bitscore, stitle, *_ in hits:
        taxon_id = taxids.get(sseqid, '')
        proteins.append({
            'protein_id': qseqid,
            'contig_id': extract_contig_id(qseqid),
            'best_hit': sseqid,
            'pident': float(pident),
            'aln_length': int(length),
            'evalue': evalue,
            'bitscore': bitscore,
            'stitle': stitle,
            'taxon_id': taxon_id,
            'taxonomy': taxdump.lineage(taxon_id) if (taxdump and taxon_id) else '',
            'pfam': pfam_map.get(sseqid, ''),
        })

    with open(args.out_proteins, 'w') as f:
        f.write('protein_id\tcontig_id\tbest_hit\tpident\taln_length\tevalue\t'
                'bitscore\tstitle\ttaxon_id\ttaxonomy\tpfam\n')
        for p in sorted(proteins, key=lambda x: (x['contig_id'], x['protein_id'])):
            f.write(f"{p['protein_id']}\t{p['contig_id']}\t{p['best_hit']}\t"
                    f"{p['pident']:.1f}\t{p['aln_length']}\t{p['evalue']}\t"
                    f"{p['bitscore']}\t{p['stitle']}\t{p['taxon_id']}\t"
                    f"{p['taxonomy']}\t{p['pfam']}\n")

    contig_data = defaultdict(lambda: {'n_proteins': 0, 'n_classified': 0,
                                       'taxonomies': [], 'pfam_domains': set()})
    for p in proteins:
        cd = contig_data[p['contig_id']]
        cd['n_proteins'] += 1
        if p['taxonomy']:
            cd['n_classified'] += 1
            cd['taxonomies'].append(p['taxonomy'])
        for domain in p['pfam'].split(';'):
            if domain.strip():
                cd['pfam_domains'].add(domain.strip())

    with open(args.out_contigs, 'w') as f:
        f.write('contig_id\tn_proteins\tn_classified\ttop_taxonomy\tpfam_domains\n')
        for cid in sorted(contig_data):
            cd = contig_data[cid]
            # Majority vote over the full lineage strings of the contig's proteins
            top_tax = Counter(cd['taxonomies']).most_common(1)[0][0] if cd['taxonomies'] else ''
            f.write(f"{cid}\t{cd['n_proteins']}\t{cd['n_classified']}\t"
                    f"{top_tax}\t{';'.join(sorted(cd['pfam_domains']))}\n")

    n_classified = sum(1 for p in proteins if p['taxonomy'])
    print(f'[INFO] MarFERReT: {len(proteins)} proteins with hits, '
          f'{n_classified}/{len(proteins)} with taxonomy, '
          f'{len(contig_data)} contigs represented', file=sys.stderr)


if __name__ == '__main__':
    main()
