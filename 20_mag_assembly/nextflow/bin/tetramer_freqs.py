#!/usr/bin/env python3
"""Calculate tetranucleotide frequencies for metagenomic contigs.

Drop-in replacement for tetramer_freqs_esom.pl (Dick et al., 2009).
Outputs a headerless TSV with contig ID and 136 reverse-complement-
collapsed tetranucleotide frequencies, compatible with the DuckDB
integration pipeline (44_tetra_db.r).

Usage:
    tetramer_freqs.py -f assembly.fasta [-min 1500] [-max 10000000] -o output.lrn
    tetramer_freqs.py -f assembly.fasta -min 1500 > output.lrn

Citation:
    Dick, G.J. et al. (2009) Community-wide analysis of microbial genome
    sequence signatures. Genome Biology 10: R85.
"""

import argparse
import sys
from itertools import product

BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def build_canonical_kmers(k=4):
    """Return sorted list of 136 canonical (RC-collapsed) tetranucleotides."""
    seen = set()
    canonical = []
    for kmer in product(BASES, repeat=k):
        kmer = "".join(kmer)
        rc = revcomp(kmer)
        if rc not in seen:
            canonical.append(kmer)
            seen.add(kmer)
    return canonical


def parse_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def mask_ambiguous(seq, min_region=50):
    """Replace short regions between ambiguous bases with N (skip in counting)."""
    seq = list(seq.upper())
    ambig_positions = [0]
    for i in range(1, len(seq)):
        if seq[i] not in "ACGT":
            ambig_positions.append(i)
    ambig_positions.append(len(seq))

    for i in range(1, len(ambig_positions)):
        length = ambig_positions[i] - ambig_positions[i - 1] + 1
        if length < min_region:
            for j in range(ambig_positions[i - 1], ambig_positions[i]):
                seq[j] = "N"
    return "".join(seq)


def split_sequence(seq, window_size):
    """Split sequence into window-sized fragments, merging short tail."""
    if len(seq) < 2 * window_size:
        return [seq]
    fragments = []
    for i in range(0, len(seq), window_size):
        fragments.append(seq[i : i + window_size])
    if len(fragments) > 1 and len(fragments[-1]) < window_size:
        fragments[-2] = fragments[-2] + fragments[-1]
        fragments.pop()
    return fragments


def count_tetras(seq, canonical, kmer_to_canonical):
    """Count tetranucleotide frequencies in a sequence, return normalized freqs."""
    counts = {k: 0 for k in canonical}
    total = 0
    for i in range(len(seq) - 3):
        kmer = seq[i : i + 4]
        canon = kmer_to_canonical.get(kmer)
        if canon is not None:
            counts[canon] += 1
            total += 1
    if total == 0:
        return [0.0] * len(canonical)
    return [counts[k] / total for k in canonical]


def main():
    parser = argparse.ArgumentParser(
        description="Calculate tetranucleotide frequencies for metagenomic contigs."
    )
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument(
        "-min", type=int, default=2500, help="Minimum contig length [2500]"
    )
    parser.add_argument(
        "-max", type=int, default=5000, help="Window size for splitting [5000]"
    )
    parser.add_argument(
        "-o", "--output", default=None, help="Output file [stdout]"
    )
    args = parser.parse_args()

    canonical = build_canonical_kmers(4)

    # Build lookup: every 4-mer (and its RC) → canonical representative
    kmer_to_canonical = {}
    for kmer in canonical:
        kmer_to_canonical[kmer] = kmer
        kmer_to_canonical[revcomp(kmer)] = kmer

    out = open(args.output, "w") if args.output else sys.stdout
    n_contigs = 0
    n_written = 0

    for header, seq in parse_fasta(args.fasta):
        n_contigs += 1
        if len(seq) < args.min:
            continue

        seq = mask_ambiguous(seq)
        fragments = split_sequence(seq, args.max)

        for idx, frag in enumerate(fragments, 1):
            contig_id = header if len(fragments) == 1 else f"{header}_{idx}"
            freqs = count_tetras(frag, canonical, kmer_to_canonical)
            freq_str = "\t".join(f"{f:.10g}" for f in freqs)
            out.write(f"{contig_id}\t{freq_str}\n")
            n_written += 1

    if out is not sys.stdout:
        out.close()

    print(
        f"[INFO] tetramer_freqs.py: {n_contigs} contigs read, {n_written} written",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
