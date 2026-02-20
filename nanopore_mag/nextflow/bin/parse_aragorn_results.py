#!/usr/bin/env python3
"""Parse Aragorn batch (-w) output into a TSV of tRNA/tmRNA genes."""
import re
import sys

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <aragorn_output> <output_tsv>", file=sys.stderr)
        sys.exit(1)

    aragorn_file = sys.argv[1]
    output_file = sys.argv[2]

    # Aragorn -w format:
    # >contig_name
    # N genes found
    # 1  tRNA-Ala  [start,end]  loop_len  (anticodon)
    # 1  tRNA-Ala  c[start,end]  loop_len  (anticodon)   <- complement
    # 1  tmRNA     c[start,end]  tag_offsets  TAG_PEPTIDE*

    # Pattern for tRNA lines
    trna_pat = re.compile(
        r'^\s*\d+\s+'
        r'(tRNA-\w+)\s+'
        r'(c?)\[(\d+),(\d+)\]\s+'
        r'(\d+)\s+'
        r'\(([a-z]+)\)'
    )
    # Pattern for tmRNA lines
    tmrna_pat = re.compile(
        r'^\s*\d+\s+'
        r'(tmRNA)\s+'
        r'(c?)\[(\d+),(\d+)\]\s+'
        r'(\S+)\s+'
        r'(\S+)'
    )

    genes = []
    contig = None
    gene_num = 0

    with open(aragorn_file) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                contig = line[1:].split()[0]
                gene_num = 0
                continue

            m = trna_pat.match(line)
            if m:
                gene_num += 1
                gene_type = m.group(1)  # tRNA-Ala
                strand = '-' if m.group(2) == 'c' else '+'
                start = int(m.group(3))
                end = int(m.group(4))
                anticodon = m.group(6)
                amino_acid = gene_type.split('-')[1] if '-' in gene_type else ''
                gene_id = f"{contig}_trna_{gene_num}"
                genes.append((gene_id, contig, start, end, strand,
                              'tRNA', amino_acid, anticodon, end - start + 1))
                continue

            m = tmrna_pat.match(line)
            if m:
                gene_num += 1
                strand = '-' if m.group(2) == 'c' else '+'
                start = int(m.group(3))
                end = int(m.group(4))
                tag_peptide = m.group(6).rstrip('*')
                gene_id = f"{contig}_tmrna_{gene_num}"
                genes.append((gene_id, contig, start, end, strand,
                              'tmRNA', '', tag_peptide, end - start + 1))

    with open(output_file, 'w') as out:
        out.write('gene_id\tcontig_id\tstart\tend\tstrand\tgene_type\tamino_acid\tanticodon_or_tag\tgene_length\n')
        for g in genes:
            out.write('\t'.join(str(x) for x in g) + '\n')

    print(f"  Parsed {len(genes)} genes from Aragorn output", file=sys.stderr)


if __name__ == '__main__':
    main()
