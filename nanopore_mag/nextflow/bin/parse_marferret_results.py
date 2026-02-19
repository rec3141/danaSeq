#!/usr/bin/env python3
"""Join DIAMOND blastp hits against MarFERReT with taxonomy and Pfam annotations.

Reads:
  - DIAMOND tabular output (blast format 6)
  - MarFERReT taxonomies file (*.taxonomies.tab.gz)
  - MarFERReT Pfam annotations file (*.best_pfam_annotations.csv.gz)

Writes:
  - marferret_proteins.tsv  Per-protein: protein_id, best_hit, pident, evalue, taxonomy, pfam
  - marferret_contigs.tsv   Per-contig: contig_id, n_proteins, n_classified, top_taxonomy, pfam_domains
"""

import sys
import os
import gzip
import argparse
from collections import defaultdict, Counter

import pandas as pd


def load_taxonomy(path):
    """Load MarFERReT taxonomy mapping: seq_id -> (taxon_id, lineage).

    The taxonomies.tab.gz has columns: ref_name, taxon_id, lineage
    (tab-separated, with header).
    """
    tax = {}
    if not path or not os.path.isfile(path):
        return tax

    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                seq_id = parts[0]
                taxon_id = parts[1]
                lineage = parts[2]
                tax[seq_id] = (taxon_id, lineage)
            elif len(parts) == 2:
                tax[parts[0]] = (parts[1], '')
    print(f'[INFO] MarFERReT: loaded {len(tax)} taxonomy entries', file=sys.stderr)
    return tax


def load_pfam(path):
    """Load MarFERReT Pfam annotations: seq_id -> pfam_str.

    The best_pfam_annotations.csv.gz is CSV with columns including
    ref_name and pfam_accession (or similar).
    """
    pfam = {}
    if not path or not os.path.isfile(path):
        return pfam

    try:
        df = pd.read_csv(path, compression='infer', low_memory=False)
    except Exception as e:
        print(f'[WARNING] MarFERReT: could not parse Pfam file: {e}', file=sys.stderr)
        return pfam

    # Identify the sequence ID column (first column or 'ref_name')
    id_col = None
    for candidate in ['ref_name', 'sequence_id', 'seq_id', 'protein_id']:
        if candidate in df.columns:
            id_col = candidate
            break
    if id_col is None:
        id_col = df.columns[0]

    # Identify Pfam column
    pfam_col = None
    for candidate in ['best_pfam', 'pfam_accession', 'pfam', 'pfam_id', 'Pfam']:
        if candidate in df.columns:
            pfam_col = candidate
            break
    if pfam_col is None:
        # Use the last column as a fallback
        pfam_col = df.columns[-1]

    for _, row in df.iterrows():
        seq_id = str(row[id_col])
        pfam_val = str(row[pfam_col]) if pd.notna(row[pfam_col]) else ''
        if pfam_val and pfam_val != 'nan':
            pfam[seq_id] = pfam_val

    print(f'[INFO] MarFERReT: loaded {len(pfam)} Pfam annotations', file=sys.stderr)
    return pfam


def extract_contig_id(protein_id):
    """Extract contig ID from MetaEuk protein header.

    MetaEuk headers follow the pattern: contig_name|... (pipe-delimited),
    where the first field is the source contig.
    """
    return protein_id.split('|')[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--diamond', required=True,
                        help='DIAMOND blast6 output file')
    parser.add_argument('--taxonomy', default=None,
                        help='MarFERReT taxonomies.tab.gz file')
    parser.add_argument('--pfam', default=None,
                        help='MarFERReT best_pfam_annotations.csv.gz file')
    parser.add_argument('--out-proteins', default='marferret_proteins.tsv',
                        help='Output per-protein TSV')
    parser.add_argument('--out-contigs', default='marferret_contigs.tsv',
                        help='Output per-contig TSV')
    args = parser.parse_args()

    # Load reference data
    tax_map = load_taxonomy(args.taxonomy)
    pfam_map = load_pfam(args.pfam)

    # Parse DIAMOND output (blast format 6)
    # Columns: qseqid sseqid pident length evalue bitscore stitle
    proteins = []
    if os.path.isfile(args.diamond) and os.path.getsize(args.diamond) > 0:
        with open(args.diamond) as f:
            for line in f:
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 6:
                    continue
                qseqid = cols[0]
                sseqid = cols[1]
                pident = float(cols[2])
                length = int(cols[3])
                evalue = cols[4]
                bitscore = cols[5]
                stitle = cols[6] if len(cols) > 6 else ''

                # Look up taxonomy and Pfam for the subject sequence
                taxon_id, lineage = tax_map.get(sseqid, ('', ''))
                pfam_annot = pfam_map.get(sseqid, '')

                contig_id = extract_contig_id(qseqid)

                proteins.append({
                    'protein_id': qseqid,
                    'contig_id': contig_id,
                    'best_hit': sseqid,
                    'pident': pident,
                    'aln_length': length,
                    'evalue': evalue,
                    'bitscore': bitscore,
                    'stitle': stitle,
                    'taxon_id': taxon_id,
                    'taxonomy': lineage,
                    'pfam': pfam_annot,
                })

    # Write per-protein TSV
    with open(args.out_proteins, 'w') as f:
        f.write('protein_id\tcontig_id\tbest_hit\tpident\taln_length\tevalue\t'
                'bitscore\tstitle\ttaxon_id\ttaxonomy\tpfam\n')
        for p in sorted(proteins, key=lambda x: (x['contig_id'], x['protein_id'])):
            f.write(f"{p['protein_id']}\t{p['contig_id']}\t{p['best_hit']}\t"
                    f"{p['pident']:.1f}\t{p['aln_length']}\t{p['evalue']}\t"
                    f"{p['bitscore']}\t{p['stitle']}\t{p['taxon_id']}\t"
                    f"{p['taxonomy']}\t{p['pfam']}\n")

    # Aggregate per-contig
    contig_data = defaultdict(lambda: {
        'n_proteins': 0,
        'n_classified': 0,
        'taxonomies': [],
        'pfam_domains': set(),
    })

    for p in proteins:
        cid = p['contig_id']
        contig_data[cid]['n_proteins'] += 1
        if p['taxonomy']:
            contig_data[cid]['n_classified'] += 1
            contig_data[cid]['taxonomies'].append(p['taxonomy'])
        if p['pfam']:
            # Pfam field may contain multiple domains separated by ;
            for domain in p['pfam'].split(';'):
                domain = domain.strip()
                if domain:
                    contig_data[cid]['pfam_domains'].add(domain)

    # Write per-contig TSV
    with open(args.out_contigs, 'w') as f:
        f.write('contig_id\tn_proteins\tn_classified\ttop_taxonomy\tpfam_domains\n')
        for cid in sorted(contig_data.keys()):
            cd = contig_data[cid]
            # Majority-vote taxonomy (most common lineage)
            if cd['taxonomies']:
                top_tax = Counter(cd['taxonomies']).most_common(1)[0][0]
            else:
                top_tax = ''
            pfam_str = ';'.join(sorted(cd['pfam_domains']))
            f.write(f"{cid}\t{cd['n_proteins']}\t{cd['n_classified']}\t"
                    f"{top_tax}\t{pfam_str}\n")

    n_proteins = len(proteins)
    n_classified = sum(1 for p in proteins if p['taxonomy'])
    n_contigs = len(contig_data)
    print(f'[INFO] MarFERReT: {n_proteins} proteins classified, '
          f'{n_classified}/{n_proteins} with taxonomy, '
          f'{n_contigs} contigs represented', file=sys.stderr)


if __name__ == '__main__':
    main()
