#!/usr/bin/env python3
"""Predict protein-coding genes with Pyrodigal (metagenomic mode).

Replacement for BAKTA_BASIC when only protein FASTA and GFF3 are needed
downstream. Uses the same Pyrodigal engine that Bakta uses internally,
but skips all annotation (DIAMOND PSC/PSCC, AMRFinderPlus, etc.).

Pyrodigal's find_genes() releases the GIL, so a ThreadPool gives true
parallel speedup across contigs.

Usage:
    pyrodigal_predict.py assembly.fasta --faa proteins.faa --gff genes.gff3 --threads 8
"""
import argparse
import io
import sys
from multiprocessing.pool import ThreadPool

import pyrodigal


def parse_fasta(path):
    """Return list of (seq_id, sequence_bytes) from a FASTA file."""
    sequences = []
    seq_id = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if seq_id is not None:
                    sequences.append((seq_id, ''.join(parts).encode()))
                seq_id = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if seq_id is not None:
        sequences.append((seq_id, ''.join(parts).encode()))
    return sequences


def main():
    parser = argparse.ArgumentParser(
        description='Predict genes with Pyrodigal (metagenomic mode)')
    parser.add_argument('assembly', help='Input assembly FASTA')
    parser.add_argument('--faa', required=True, help='Output protein FASTA')
    parser.add_argument('--gff', required=True, help='Output GFF3')
    parser.add_argument('--threads', type=int, default=1,
                        help='Number of parallel threads (default: 1)')
    parser.add_argument('--min-gene-length', type=int, default=90,
                        help='Minimum gene length in bp (default: 90)')
    args = parser.parse_args()

    gene_finder = pyrodigal.GeneFinder(meta=True, min_gene=args.min_gene_length)

    print(f'[INFO] Loading assembly: {args.assembly}', file=sys.stderr)
    sequences = parse_fasta(args.assembly)
    print(f'[INFO] Loaded {len(sequences)} sequences, '
          f'predicting genes with {args.threads} threads...', file=sys.stderr)

    # Pyrodigal is thread-safe and releases the GIL — ThreadPool gives true parallelism
    # https://pyrodigal.readthedocs.io/en/stable/guide/parallel.html
    seq_bytes_list = [seq_bytes for _, seq_bytes in sequences]

    with ThreadPool(processes=args.threads) as pool:
        all_genes = pool.map(gene_finder.find_genes, seq_bytes_list)

    n_genes = 0
    with open(args.faa, 'w') as faa_fh, open(args.gff, 'w') as gff_fh:
        gff_fh.write('##gff-version 3\n')
        for (seq_id, _), genes in zip(sequences, all_genes):
            n_genes += len(genes)
            genes.write_translations(faa_fh, sequence_id=seq_id)
            genes.write_gff(gff_fh, sequence_id=seq_id, header=False)

    print(f'[INFO] Pyrodigal: {n_genes} genes from {len(sequences)} sequences',
          file=sys.stderr)


if __name__ == '__main__':
    main()
