#!/usr/bin/env python3
"""Parse barrnap + vsearch results into per-gene and per-contig TSVs."""

import sys
import os
import re
import argparse
from collections import defaultdict

# Expected full-length rRNA gene sizes (bp) for completeness calculation
EXPECTED_LENGTHS = {
    '16S_rRNA': 1542,
    '23S_rRNA': 2904,
    '5S_rRNA': 120,
    '5_8S_rRNA': 157,
    '18S_rRNA': 1800,
    '28S_rRNA': 3400,
}

KINGDOM_MAP = {'bac': 'Bacteria', 'arc': 'Archaea', 'euk': 'Eukaryota'}
SSU_TYPES = {'16S_rRNA', '18S_rRNA'}
LSU_TYPES = {'23S_rRNA', '28S_rRNA'}

# ── Parse arguments ──────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument('--taxonomy-map', help='TSV mapping accession to taxonomy (from SILVA headers)')
args = parser.parse_args()

# ── Load accession → taxonomy map ────────────────────────────────────
acc_to_tax = {}
if args.taxonomy_map and os.path.isfile(args.taxonomy_map):
    with open(args.taxonomy_map) as f:
        for line in f:
            parts = line.strip().split('\t', 1)
            if len(parts) == 2:
                acc_to_tax[parts[0]] = parts[1]
    print(f'[INFO] Loaded {len(acc_to_tax)} accession→taxonomy mappings', file=sys.stderr)

# ── Parse barrnap GFF files ──────────────────────────────────────────
# barrnap --outseq headers: >TYPE::CONTIG:START-END(STRAND)
# barrnap GFF: contig \t barrnap:VERSION \t rRNA \t start \t end \t evalue \t strand \t . \t Name=TYPE;product=...
genes = []  # list of dicts

for kingdom_code in ['bac', 'arc', 'euk']:
    gff_path = f'barrnap_{kingdom_code}.gff'
    if not os.path.isfile(gff_path):
        continue
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#') or '\t' not in line:
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9 or cols[2] != 'rRNA':
                continue
            contig = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            score = cols[5]  # barrnap e-value
            strand = cols[6]
            attrs = cols[8]

            # Extract Name= attribute for rRNA type
            m = re.search(r'Name=([^;]+)', attrs)
            rrna_type = m.group(1) if m else 'unknown_rRNA'

            gene_length = abs(end - start) + 1
            expected = EXPECTED_LENGTHS.get(rrna_type, gene_length)
            completeness = min(gene_length / expected, 1.0) if expected > 0 else 0.0

            # Build gene_id matching barrnap --outseq header format
            # barrnap outseq uses 0-based start; GFF uses 1-based
            gene_id = f'{rrna_type}::{contig}:{start - 1}-{end}({strand})'

            genes.append({
                'gene_id': gene_id,
                'contig_id': contig,
                'start': start,
                'end': end,
                'strand': strand,
                'rrna_type': rrna_type,
                'kingdom': KINGDOM_MAP[kingdom_code],
                'kingdom_code': kingdom_code,
                'gene_length': gene_length,
                'gene_completeness': completeness,
                'barrnap_score': score,
                'best_match': '',
                'vsearch_identity': 0.0,
                'vsearch_coverage': 0.0,
                'taxonomy': '',
            })

# ── Parse vsearch blast6 output ──────────────────────────────────────
# blast6 columns: query, target, identity, alnlen, mismatch, gaps, qstart, qend, tstart, tend, evalue, bitscore
def parse_vsearch(tsv_path, db_type):
    """Parse vsearch blast6out into dict keyed by query ID."""
    hits = {}
    if not os.path.isfile(tsv_path) or os.path.getsize(tsv_path) == 0:
        return hits
    with open(tsv_path) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue
            query = cols[0]
            if query.startswith('*'):
                continue  # vsearch no-hit marker
            target = cols[1]
            identity = float(cols[2])
            alnlen = int(cols[3])

            # Look up taxonomy from accession→taxonomy map (built from SILVA headers)
            # blast6out only stores the accession (first word of the header)
            tax = acc_to_tax.get(target, '')

            hits[query] = {
                'target': target,
                'identity': identity,
                'alnlen': alnlen,
                'taxonomy': tax,
                'db_type': db_type,
            }
    return hits

# Collect all vsearch hits (SSU + LSU, per kingdom)
vsearch_hits = {}  # gene_id -> best hit dict
for kingdom_code in ['bac', 'arc', 'euk']:
    for db_type, prefix in [('SSU', 'vsearch_ssu'), ('LSU', 'vsearch_lsu')]:
        tsv_path = f'{prefix}_{kingdom_code}.tsv'
        hits = parse_vsearch(tsv_path, db_type)
        for qid, hit in hits.items():
            # vsearch query IDs match barrnap --outseq headers
            # Keep the best hit by identity if duplicated
            if qid not in vsearch_hits or hit['identity'] > vsearch_hits[qid]['identity']:
                vsearch_hits[qid] = hit

# ── Match vsearch hits to barrnap genes ──────────────────────────────
# barrnap --outseq uses slightly different ID format from GFF
# outseq: >TYPE::CONTIG:START-END(STRAND) — but may truncate long contig names
# We match by building both formats and trying both

# Build a lookup from outseq-style IDs
for gene in genes:
    gid = gene['gene_id']
    hit = vsearch_hits.get(gid)
    if hit:
        gene['best_match'] = hit['target']
        gene['vsearch_identity'] = hit['identity']
        gene['vsearch_coverage'] = hit['alnlen'] / gene['gene_length'] if gene['gene_length'] > 0 else 0.0
        gene['taxonomy'] = hit['taxonomy']

# ── Write rrna_genes.tsv ─────────────────────────────────────────────
with open('rrna_genes.tsv', 'w') as f:
    f.write('gene_id\tcontig_id\tstart\tend\tstrand\trrna_type\tkingdom\t'
            'gene_length\tgene_completeness\tbarrnap_score\t'
            'best_match\tvsearch_identity\tvsearch_coverage\ttaxonomy\n')
    for g in sorted(genes, key=lambda x: (x['contig_id'], x['start'])):
        f.write(f"{g['gene_id']}\t{g['contig_id']}\t{g['start']}\t{g['end']}\t"
                f"{g['strand']}\t{g['rrna_type']}\t{g['kingdom']}\t"
                f"{g['gene_length']}\t{g['gene_completeness']:.3f}\t{g['barrnap_score']}\t"
                f"{g['best_match']}\t{g['vsearch_identity']:.1f}\t"
                f"{g['vsearch_coverage']:.3f}\t{g['taxonomy']}\n")

# ── Write rrna_contigs.tsv ───────────────────────────────────────────
contig_data = defaultdict(lambda: {
    'genes': [],
    'best_ssu': None,
    'best_lsu': None,
})

for g in genes:
    cid = g['contig_id']
    contig_data[cid]['genes'].append(g)

    # Track best SSU and LSU hit per contig
    if g['rrna_type'] in SSU_TYPES and g['vsearch_identity'] > 0:
        cur = contig_data[cid]['best_ssu']
        if cur is None or g['vsearch_identity'] > cur['vsearch_identity']:
            contig_data[cid]['best_ssu'] = g

    if g['rrna_type'] in LSU_TYPES and g['vsearch_identity'] > 0:
        cur = contig_data[cid]['best_lsu']
        if cur is None or g['vsearch_identity'] > cur['vsearch_identity']:
            contig_data[cid]['best_lsu'] = g

with open('rrna_contigs.tsv', 'w') as f:
    f.write('contig_id\tn_rrna_genes\trrna_types\tkingdoms\t'
            'best_ssu_taxonomy\tbest_ssu_identity\t'
            'best_lsu_taxonomy\tbest_lsu_identity\n')
    for cid in sorted(contig_data.keys()):
        cd = contig_data[cid]
        n = len(cd['genes'])
        types = sorted(set(g['rrna_type'] for g in cd['genes']))
        kingdoms = sorted(set(g['kingdom'] for g in cd['genes']))

        ssu_tax = cd['best_ssu']['taxonomy'] if cd['best_ssu'] else ''
        ssu_id = f"{cd['best_ssu']['vsearch_identity']:.1f}" if cd['best_ssu'] else '0.0'
        lsu_tax = cd['best_lsu']['taxonomy'] if cd['best_lsu'] else ''
        lsu_id = f"{cd['best_lsu']['vsearch_identity']:.1f}" if cd['best_lsu'] else '0.0'

        f.write(f"{cid}\t{n}\t{','.join(types)}\t{','.join(kingdoms)}\t"
                f"{ssu_tax}\t{ssu_id}\t{lsu_tax}\t{lsu_id}\n")

# Summary stats
n_genes = len(genes)
n_classified = sum(1 for g in genes if g['vsearch_identity'] > 0)
n_contigs = len(contig_data)
print(f'[INFO] rRNA: {n_genes} genes detected on {n_contigs} contigs, '
      f'{n_classified}/{n_genes} classified by vsearch', file=sys.stderr)
