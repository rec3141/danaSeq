#!/usr/bin/env python3
"""
Aggregate per-task RNA-seq outputs into the JSONs consumed by the Svelte viz.

Inputs (directories of files emitted by the Nextflow pipeline):
  --idxstats-dir    *.idxstats.tsv   (per sample × reference)
  --flagstat-dir    *.flagstat.txt
  --covstats-dir    *.covstats.txt
  --readcounts-dir  *.count.tsv      (per stage of preprocessing)
  --genecounts-dir  <ref>.gene_counts.tsv

Outputs (written into --out):
  overview.json     — per-stage read totals, per-reference mapping rates
  samples.json      — per-sample summary rows
  references.json   — per-reference metadata + contig list
  expression.json.gz — gene × sample count matrix per reference
  read_flow.json    — read counts by (sample, stage) for the funnel chart
"""
import argparse, csv, glob, gzip, json, os, re, sys
from collections import defaultdict

FILE_RE = re.compile(r'(?P<sample>[^/]+)_vs_(?P<ref>[^/]+?)\.(?:idxstats\.tsv|flagstat\.txt|covstats\.txt|covhist\.txt)$')


def parse_flagstat(path: str) -> dict:
    out = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split(' ', 2)
            if not parts: continue
            try:
                n = int(parts[0])
            except ValueError:
                continue
            if 'in total' in line:                 out['total']     = n
            elif 'primary' == parts[-1]:           out['primary']   = n
            elif 'mapped (' in line and 'primary' not in line:
                out['mapped']    = n
            elif 'properly paired' in line:        out['proper_pair'] = n
            elif 'duplicates' in line:             out['duplicates'] = n
    return out


def parse_idxstats(path: str) -> list:
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4: continue
            contig, length, mapped, unmapped = parts[0], parts[1], parts[2], parts[3]
            if contig == '*': continue
            try:
                rows.append({
                    'contig': contig,
                    'length': int(length),
                    'mapped': int(mapped),
                    'unmapped': int(unmapped),
                })
            except ValueError:
                continue
    return rows


def parse_covstats(path: str) -> list:
    """BBmap covstats: header line then per-contig rows."""
    rows = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return rows
    with open(path) as f:
        header = None
        for line in f:
            line = line.rstrip('\n')
            if not line: continue
            if header is None:
                header = line.lstrip('#').split('\t')
                continue
            parts = line.split('\t')
            if len(parts) != len(header): continue
            row = dict(zip(header, parts))
            try:
                rows.append({
                    'contig': row.get('ID') or row.get('#ID') or parts[0],
                    'avg_cov': float(row.get('Avg_fold', 0) or 0),
                    'covered_pct': float(row.get('Covered_percent', 0) or 0),
                    'length': int(row.get('Length', 0) or 0),
                    'reads': int(row.get('Plus_reads', 0) or 0) + int(row.get('Minus_reads', 0) or 0),
                })
            except ValueError:
                continue
    return rows


def parse_readcount(path: str) -> tuple:
    with open(path) as f:
        row = next(csv.reader(f, delimiter='\t'), None)
        if not row or len(row) < 3: return None
        try:
            return row[0], row[1], int(row[2])
        except ValueError:
            return None


def parse_genecounts(path: str) -> dict:
    """Read merged gene_counts.tsv into a dict with genes, samples, matrix."""
    ref = os.path.basename(path).replace('.gene_counts.tsv', '')
    genes, samples, matrix = [], None, []
    with open(path) as f:
        header = next(csv.reader(f, delimiter='\t'), None)
        if not header or len(header) < 2:
            return {'ref': ref, 'genes': [], 'samples': [], 'matrix': []}
        samples = header[1:]
        for row in csv.reader(f, delimiter='\t'):
            if len(row) < 2: continue
            genes.append(row[0])
            try:
                matrix.append([int(x) for x in row[1:]])
            except ValueError:
                matrix.append([0] * (len(row) - 1))
    return {'ref': ref, 'genes': genes, 'samples': samples or [], 'matrix': matrix}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--idxstats-dir',   required=True)
    ap.add_argument('--flagstat-dir',   required=True)
    ap.add_argument('--covstats-dir',   required=True)
    ap.add_argument('--readcounts-dir', required=True)
    ap.add_argument('--genecounts-dir', required=True)
    ap.add_argument('--out',            required=True)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # --- read counts by (sample, stage) ---
    flow = defaultdict(dict)
    for path in glob.glob(os.path.join(args.readcounts_dir, '*.count.tsv')):
        parsed = parse_readcount(path)
        if parsed:
            sample, stage, n = parsed
            flow[sample][stage] = n
    samples_seen = sorted(flow.keys())

    # --- flagstat + idxstats + covstats per (sample, ref) ---
    sample_ref_stats = defaultdict(dict)  # (sample, ref) -> dict
    refs_seen = set()
    for path in glob.glob(os.path.join(args.flagstat_dir, '*.flagstat.txt')):
        m = FILE_RE.search(path)
        if not m: continue
        s, r = m.group('sample'), m.group('ref')
        refs_seen.add(r)
        sample_ref_stats[(s, r)]['flagstat'] = parse_flagstat(path)
    for path in glob.glob(os.path.join(args.idxstats_dir, '*.idxstats.tsv')):
        m = FILE_RE.search(path)
        if not m: continue
        s, r = m.group('sample'), m.group('ref')
        refs_seen.add(r)
        sample_ref_stats[(s, r)]['idxstats'] = parse_idxstats(path)
    for path in glob.glob(os.path.join(args.covstats_dir, '*.covstats.txt')):
        m = FILE_RE.search(path)
        if not m: continue
        s, r = m.group('sample'), m.group('ref')
        refs_seen.add(r)
        sample_ref_stats[(s, r)]['covstats'] = parse_covstats(path)

    refs_seen = sorted(refs_seen)

    # --- samples.json ---
    samples_out = []
    for s in samples_seen:
        per_ref = {}
        for r in refs_seen:
            fs = sample_ref_stats.get((s, r), {}).get('flagstat') or {}
            total  = fs.get('primary') or fs.get('total') or 0
            mapped = fs.get('mapped') or 0
            per_ref[r] = {
                'total_reads': total,
                'mapped_reads': mapped,
                'mapping_rate': (mapped / total) if total else 0.0,
            }
        samples_out.append({
            'id': s,
            'read_flow': flow.get(s, {}),
            'per_reference': per_ref,
        })

    # --- references.json ---
    refs_out = []
    for r in refs_seen:
        contigs = []
        # Aggregate per-contig totals across samples (use first sample's idxstats for contig metadata)
        contig_meta = {}
        contig_reads = defaultdict(int)
        for s in samples_seen:
            stats = sample_ref_stats.get((s, r), {})
            for row in stats.get('idxstats', []) or []:
                contig_meta[row['contig']] = row['length']
                contig_reads[row['contig']] += row['mapped']
        for contig, length in sorted(contig_meta.items()):
            contigs.append({
                'name': contig,
                'length': length,
                'reads': contig_reads.get(contig, 0),
            })
        refs_out.append({
            'name': r,
            'n_contigs': len(contigs),
            'total_length': sum(c['length'] for c in contigs),
            'contigs': contigs,
        })

    # --- overview.json ---
    overview = {
        'n_samples': len(samples_seen),
        'n_references': len(refs_seen),
        'samples': samples_seen,
        'references': refs_seen,
        'mapping_rate_matrix': [
            [(samples_out[i]['per_reference'].get(r, {}).get('mapping_rate') or 0.0)
             for r in refs_seen]
            for i in range(len(samples_out))
        ],
    }

    with open(os.path.join(args.out, 'overview.json'), 'w') as f:
        json.dump(overview, f)
    with open(os.path.join(args.out, 'samples.json'), 'w') as f:
        json.dump(samples_out, f)
    with open(os.path.join(args.out, 'references.json'), 'w') as f:
        json.dump(refs_out, f)
    with open(os.path.join(args.out, 'read_flow.json'), 'w') as f:
        json.dump({'samples': samples_seen, 'stages': flow}, f)

    # --- expression.json.gz ---
    expression = {}
    for path in sorted(glob.glob(os.path.join(args.genecounts_dir, '*.gene_counts.tsv'))):
        gc = parse_genecounts(path)
        if gc['genes']:
            expression[gc['ref']] = gc
    with gzip.open(os.path.join(args.out, 'expression.json.gz'), 'wt') as f:
        json.dump(expression, f)

    print(f"[viz] wrote overview.json ({len(samples_seen)} samples, {len(refs_seen)} refs)", file=sys.stderr)
    print(f"[viz] wrote expression for {len(expression)} reference(s)", file=sys.stderr)


if __name__ == '__main__':
    main()
