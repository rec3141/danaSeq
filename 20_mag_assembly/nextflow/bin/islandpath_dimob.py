#!/usr/bin/env python3
"""
IslandPath-DIMOB — Genomic island detection via dinucleotide bias + mobility genes.

Python reimplementation of the Perl IslandPath-DIMOB tool (Bertelli et al., 2018;
Hsiao et al., 2003). Works directly with GFF3 + FASTA + FAA outputs from Prokka,
removing the dependency on GenBank (.gbk) input and broken Perl BioPerl modules.

Algorithm:
  1. Per-contig processing (each contig = different organism in metagenome)
  2. Extract CDS nucleotide sequences from assembly using GFF coordinates
  3. Dinucleotide relative abundance (rho*) via 4-frame non-overlapping counting
  4. 6-gene sliding window delta* = mean |window_rho - contig_rho| / 16
  5. Threshold: strong (median + 2*SD → 6 genes), moderate (median + SD → 3 genes)
  6. Cluster consecutive marked genes, merge gaps < 6, keep clusters >= 8 genes
  7. HMM scan (hmmscan, E-value <= 1e-7) + keyword matching for mobility genes
  8. Keep islands with >= 1 mobility gene, discard islands < 2000 bp

Reference:
  Bertelli C, et al. (2018) IslandViewer 4. Nucleic Acids Res 46(W1):W59-W64.
  Hsiao W, et al. (2003) IslandPath. Bioinformatics 19(3):418-420.
"""

import argparse
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

# ---------------------------------------------------------------------------
# FASTA / GFF parsing
# ---------------------------------------------------------------------------

def parse_fasta(path):
    """Read a FASTA file and return {name: sequence} dict."""
    seqs = {}
    name = None
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = ''.join(parts)
    return seqs


def parse_gff_cds(path):
    """Extract CDS records per contig from a GFF3 file.

    Returns dict: contig_id -> list of CDS dicts sorted by start position.
    Each CDS dict has keys: contig, start, end, strand, gene_id, product.
    Stops at ##FASTA section (Prokka embeds sequence there).
    """
    cds_by_contig = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue
            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]
            gene_id = ''
            product = ''
            for attr in attrs.split(';'):
                if attr.startswith('ID='):
                    gene_id = attr[3:]
                elif attr.startswith('product='):
                    product = attr[8:]
            cds_by_contig[contig].append({
                'contig': contig,
                'start': start,
                'end': end,
                'strand': strand,
                'gene_id': gene_id,
                'product': product,
            })
    # Sort CDS by start position within each contig
    for contig in cds_by_contig:
        cds_by_contig[contig].sort(key=lambda c: c['start'])
    return dict(cds_by_contig)


def revcomp(seq):
    """Reverse complement a DNA sequence."""
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]


def extract_cds_nuc(contig_seq, cds):
    """Extract CDS nucleotide sequence from contig using GFF coordinates.

    GFF is 1-based inclusive. Returns upper-case sequence on the coding strand.
    """
    # GFF coordinates are 1-based, Python slicing is 0-based
    seq = contig_seq[cds['start'] - 1 : cds['end']]
    if cds['strand'] == '-':
        seq = revcomp(seq)
    return seq.upper()

# ---------------------------------------------------------------------------
# Dinucleotide bias calculation
# ---------------------------------------------------------------------------

DINUC_KEYS = [
    'AA', 'AC', 'AG', 'AT',
    'CC', 'CA', 'CG', 'CT',
    'GG', 'GA', 'GC', 'GT',
    'TT', 'TA', 'TC', 'TG',
]

VALID_BASES = set('ACGT')


def count_dimers_nonoverlap(seq):
    """Count non-overlapping dimers in seq (step=2).

    Returns (mono_counts, di_counts) as dicts.
    Matches original Perl SeqWords2: step-2 counting, NOT overlapping.
    """
    mono = defaultdict(int)
    di = defaultdict(int)
    # Ensure even length (truncate last base if odd)
    end = len(seq) if len(seq) % 2 == 0 else len(seq) - 1
    for i in range(0, end, 2):
        b1 = seq[i]
        b2 = seq[i + 1] if i + 1 < len(seq) else None
        if b1 in VALID_BASES:
            mono[b1] += 1
        if b2 is not None and b2 in VALID_BASES:
            mono[b2] += 1
        if b2 is not None and b1 in VALID_BASES and b2 in VALID_BASES:
            di[b1 + b2] += 1
    return mono, di


def count_dimers_4frames(seq):
    """Count dimers across 4 frames: original, shifted+1, revcomp, revcomp+shifted.

    This matches the original Perl implementation which combines counts from all
    4 reading frames to get a robust estimate.
    Returns combined (mono_counts, di_counts).
    """
    total_mono = defaultdict(int)
    total_di = defaultdict(int)

    rc = revcomp(seq)

    for s in [seq, seq[1:], rc, rc[1:]]:
        mono, di = count_dimers_nonoverlap(s)
        for k, v in mono.items():
            total_mono[k] += v
        for k, v in di.items():
            total_di[k] += v

    # Ensure all 16 dinucleotide keys are present
    for dk in DINUC_KEYS:
        if dk not in total_di:
            total_di[dk] = 0
    return total_mono, total_di


def calc_relative_abundance(mono, di):
    """Calculate rho*(XY) = f(XY) / (f(X) * f(Y)) for all 16 dinucleotides.

    f(XY) = count(XY) / total_dinucs
    f(X)  = count(X)  / total_mononucs
    """
    total_mono = sum(mono.values())
    total_di = sum(di.values())
    if total_mono == 0 or total_di == 0:
        return {dk: 1.0 for dk in DINUC_KEYS}

    rho = {}
    for dk in DINUC_KEYS:
        b1, b2 = dk[0], dk[1]
        f_di = di.get(dk, 0) / total_di
        f_b1 = mono.get(b1, 0) / total_mono
        f_b2 = mono.get(b2, 0) / total_mono
        denom = f_b1 * f_b2
        rho[dk] = f_di / denom if denom > 0 else 1.0
    return rho


def calc_dinuc_bias(cds_list, contig_seq):
    """Calculate per-CDS dinucleotide bias using 6-gene sliding windows.

    For each contig:
    1. Count mono/di across all CDS → genome-level rho*
    2. Sliding window of 6 CDS → window-level rho*
    3. delta* = mean(|window_rho - genome_rho|) across 16 dinucs * 1000

    Returns list of bias values (one per sliding window position).
    Length = max(0, n_cds - 5).
    """
    n = len(cds_list)
    if n < 6:
        return []

    # Extract nucleotide sequences and per-CDS counts
    cds_monos = []
    cds_dis = []
    for cds in cds_list:
        nuc_seq = extract_cds_nuc(contig_seq, cds)
        mono, di = count_dimers_4frames(nuc_seq)
        cds_monos.append(mono)
        cds_dis.append(di)

    # Genome-level (whole contig) counts: sum all CDS
    genome_mono = defaultdict(int)
    genome_di = defaultdict(int)
    for mono in cds_monos:
        for k, v in mono.items():
            genome_mono[k] += v
    for di in cds_dis:
        for k, v in di.items():
            genome_di[k] += v
    genome_rho = calc_relative_abundance(genome_mono, genome_di)

    # 6-gene sliding window bias
    biases = []
    for i in range(n - 5):
        win_mono = defaultdict(int)
        win_di = defaultdict(int)
        for j in range(6):
            for k, v in cds_monos[i + j].items():
                win_mono[k] += v
            for k, v in cds_dis[i + j].items():
                win_di[k] += v
        win_rho = calc_relative_abundance(win_mono, win_di)

        # delta* = mean |window - genome| across 16 dinucleotides, * 1000
        delta = sum(abs(win_rho[dk] - genome_rho[dk]) for dk in DINUC_KEYS) / 16.0
        biases.append(delta * 1000)

    return biases

# ---------------------------------------------------------------------------
# Island detection (thresholding + clustering)
# ---------------------------------------------------------------------------

def find_dinuc_islands(biases, min_cluster=8):
    """Identify clusters of biased genes from sliding window bias values.

    Thresholds (from original):
      - strong: median + 2*SD → mark all 6 genes in the window
      - moderate: median + SD → mark first 3 genes in the window

    Clustering:
      - Collect marked gene indices
      - Find consecutive runs
      - Keep runs >= min_cluster genes
      - Merge runs with gap < 6 genes apart
      - Exclude last gene of each island (per original precision tuning)

    Returns list of (start_idx, end_idx) tuples into the CDS list (inclusive).
    """
    if len(biases) < 1:
        return []

    # Median and SD
    sorted_vals = sorted(biases)
    n = len(sorted_vals)
    if n % 2 == 0:
        median = (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2.0
    else:
        median = sorted_vals[n // 2]

    mean_val = sum(biases) / n
    variance = sum((x - mean_val) ** 2 for x in biases) / n
    sd = math.sqrt(variance)

    if sd == 0:
        return []

    strong_thresh = median + 2 * sd
    moderate_thresh = median + sd

    # Mark gene indices that show bias
    marked = set()
    for i, bias in enumerate(biases):
        if bias > strong_thresh:
            # Mark all 6 genes in this window
            for j in range(6):
                marked.add(i + j)
        elif bias > moderate_thresh:
            # Mark first 3 genes
            for j in range(3):
                marked.add(i + j)

    if not marked:
        return []

    # Sort and find consecutive clusters
    marked_sorted = sorted(marked)

    clusters = []
    current = [marked_sorted[0]]
    for idx in marked_sorted[1:]:
        if idx == current[-1] + 1:
            current.append(idx)
        else:
            clusters.append(current)
            current = [idx]
    clusters.append(current)

    # Keep clusters >= min_cluster genes
    clusters = [c for c in clusters if len(c) >= min_cluster]

    if not clusters:
        return []

    # Merge clusters that are < 6 genes apart
    merged = [clusters[0]]
    for cluster in clusters[1:]:
        gap = cluster[0] - merged[-1][-1] - 1
        if gap < 6:
            # Fill the gap
            bridge = list(range(merged[-1][-1] + 1, cluster[0]))
            merged[-1] = merged[-1] + bridge + cluster
        else:
            merged.append(cluster)

    # Exclude last gene of each island (per original: improves precision)
    # Return (start_idx, end_idx) into CDS list, inclusive
    islands = []
    for cluster in merged:
        if len(cluster) > 1:
            islands.append((cluster[0], cluster[-2]))
        # Single-gene cluster after exclusion: skip (can't happen with min_cluster=8)

    return islands

# ---------------------------------------------------------------------------
# Mobility gene detection
# ---------------------------------------------------------------------------

MOBILITY_KEYWORDS = [
    'transposase',
    'istb',
    'insertion element',
    'recombinase',
    'insertion sequence',
    'resolvase',
    'integrase',
    'phage',
    'transposon',
    'transposable element',
    'excisionase',
]

MOBILITY_PATTERNS = [re.compile(kw, re.IGNORECASE) for kw in MOBILITY_KEYWORDS]


def find_keyword_mobility(cds_list):
    """Return set of gene_ids whose product annotation matches mobility keywords."""
    mob_ids = set()
    for cds in cds_list:
        product = cds.get('product', '')
        if any(pat.search(product) for pat in MOBILITY_PATTERNS):
            mob_ids.add(cds['gene_id'])
    return mob_ids


def run_hmmscan(faa_path, hmm_db, cpus=1, evalue=1e-7):
    """Run hmmscan and return set of protein IDs with significant HMM hits.

    Uses --domtblout for easy parsing. Returns set of query names (Prokka locus tags).
    """
    mob_ids = set()

    if not os.path.isfile(faa_path) or os.path.getsize(faa_path) == 0:
        return mob_ids

    with tempfile.NamedTemporaryFile(mode='w', suffix='.domtbl', delete=False) as tmp:
        domtbl_path = tmp.name

    try:
        cmd = [
            'hmmscan',
            '--domtblout', domtbl_path,
            '--noali',
            '-E', str(evalue),
            '--cpu', str(cpus),
            hmm_db,
            faa_path,
        ]
        subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=False,
        )

        if os.path.isfile(domtbl_path):
            with open(domtbl_path) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    fields = line.split()
                    if len(fields) >= 13:
                        query_name = fields[3]
                        domain_evalue = float(fields[12])
                        if domain_evalue <= evalue:
                            mob_ids.add(query_name)
    finally:
        if os.path.isfile(domtbl_path):
            os.unlink(domtbl_path)

    return mob_ids

# ---------------------------------------------------------------------------
# Per-contig pipeline
# ---------------------------------------------------------------------------

def process_contig(contig_id, contig_seq, cds_list, hmm_mob_ids, min_island_bp=2000):
    """Run full IslandPath-DIMOB on one contig.

    Returns list of dicts: {contig, start, end} for each island.
    """
    n_cds = len(cds_list)
    if n_cds < 8:
        return []

    # 1. Dinucleotide bias
    biases = calc_dinuc_bias(cds_list, contig_seq)
    if not biases:
        return []

    # 2. Dinucleotide islands
    dinuc_islands = find_dinuc_islands(biases, min_cluster=8)
    if not dinuc_islands:
        return []

    # 3. Keyword-based mobility genes
    keyword_mob = find_keyword_mobility(cds_list)

    # 4. Combined mobility gene set
    all_mob_ids = hmm_mob_ids | keyword_mob

    # 5. Filter: keep islands with >= 1 mobility gene
    islands = []
    for start_idx, end_idx in dinuc_islands:
        island_cds = cds_list[start_idx:end_idx + 1]
        island_gene_ids = {cds['gene_id'] for cds in island_cds}

        has_mobility = bool(island_gene_ids & all_mob_ids)
        if not has_mobility:
            continue

        # Size filter
        island_start = island_cds[0]['start']
        island_end = island_cds[-1]['end']
        island_size = island_end - island_start + 1
        if island_size < min_island_bp:
            continue

        islands.append({
            'contig': contig_id,
            'start': island_start,
            'end': island_end,
        })

    return islands

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='IslandPath-DIMOB: genomic island detection via dinucleotide bias + mobility genes',
    )
    parser.add_argument('--gff', required=True, help='Prokka GFF3 file')
    parser.add_argument('--fasta', required=True, help='Assembly FASTA')
    parser.add_argument('--faa', required=True, help='Prokka protein FASTA (.faa)')
    parser.add_argument('--hmm_db', required=True, help='Path to HMM database (Pfam mobility gene profiles)')
    parser.add_argument('--cpus', type=int, default=1, help='CPUs for hmmscan')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    # Validate inputs
    for path, label in [(args.gff, 'GFF'), (args.fasta, 'FASTA'), (args.faa, 'FAA')]:
        if not os.path.isfile(path):
            print(f'[ERROR] {label} file not found: {path}', file=sys.stderr)
            sys.exit(1)

    # Parse inputs
    print(f'[INFO] Parsing assembly FASTA: {args.fasta}', file=sys.stderr)
    assembly = parse_fasta(args.fasta)

    print(f'[INFO] Parsing GFF3: {args.gff}', file=sys.stderr)
    cds_by_contig = parse_gff_cds(args.gff)

    # Run batch HMM scan on all proteins
    print(f'[INFO] Running hmmscan on {args.faa} ({args.cpus} CPUs)', file=sys.stderr)
    hmm_mob_ids = run_hmmscan(args.faa, args.hmm_db, cpus=args.cpus, evalue=1e-7)
    print(f'[INFO] HMM mobility genes found: {len(hmm_mob_ids)}', file=sys.stderr)

    # Process each contig
    all_islands = []
    contigs_processed = 0
    contigs_with_islands = 0

    for contig_id in sorted(cds_by_contig.keys()):
        if contig_id not in assembly:
            print(f'[WARNING] Contig {contig_id} in GFF but not in FASTA — skipping', file=sys.stderr)
            continue

        cds_list = cds_by_contig[contig_id]
        contig_seq = assembly[contig_id]
        contigs_processed += 1

        islands = process_contig(contig_id, contig_seq, cds_list, hmm_mob_ids)
        if islands:
            contigs_with_islands += 1
            all_islands.extend(islands)

    # Write output
    print(f'[INFO] Contigs processed: {contigs_processed}', file=sys.stderr)
    print(f'[INFO] Contigs with islands: {contigs_with_islands}', file=sys.stderr)
    print(f'[INFO] Total islands found: {len(all_islands)}', file=sys.stderr)

    with open(args.output, 'w') as fh:
        fh.write('island_id\tcontig\tstart\tend\n')
        for i, island in enumerate(all_islands, 1):
            fh.write(f"GI_{i}\t{island['contig']}\t{island['start']}\t{island['end']}\n")

    print(f'[INFO] Results written to {args.output}', file=sys.stderr)


if __name__ == '__main__':
    main()
