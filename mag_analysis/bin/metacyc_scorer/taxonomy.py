"""GTDB lineage -> NCBI taxids, for MetaCyc TAXONOMIC-RANGE rules.

Reference build (open data): GTDB bac120/ar53 metadata (gtdb_taxonomy, ncbi_taxid per genome)
plus NCBI taxdump. For every GTDB taxon at every rank, each NCBI taxid in the lineage of its
genomes is counted; taxids carried by > 50% of the taxon's genomes form the NCBI lineage of
that GTDB taxon (the taxid version of GTDB-Tk's gtdb_to_ncbi_majority_vote).
Output: gtdb2ncbi.tsv.gz  gtdb_taxon  n_genomes  taxid:fraction,...

MAG lineages come from GTDB-Tk (read_gtdbtk). read_sendsketch_bins() is a fallback for runs
without GTDB-Tk (--metacyc_taxonomy sendsketch): a length-weighted majority over each bin's
SendSketch contig calls, which are less reliable than placement in the GTDB tree.
"""
import gzip
import os
import sys
from collections import Counter, defaultdict

RANKS = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')


RANK_LEVEL = {'superkingdom': 0, 'domain': 0, 'acellular root': 0, 'phylum': 1, 'class': 2, 'order': 3,
              'family': 4, 'genus': 5, 'species': 6}


def _taxdump(d):
    parent, rank = {}, {}
    with open(os.path.join(d, 'nodes.dmp')) as fh:
        for line in fh:
            f = line.split('\t|\t')
            parent[f[0]] = f[1]
            rank[f[0]] = f[2]
    merged = {}
    p = os.path.join(d, 'merged.dmp')
    if os.path.exists(p):
        with open(p) as fh:
            for line in fh:
                f = line.rstrip('\t|\n').split('\t|\t')
                merged[f[0]] = f[1]
    return parent, merged, rank


def build(metadata_files, taxdump_dir, out_path, min_frac=0.5, log=sys.stderr):
    parent, merged, rank = _taxdump(taxdump_dir)
    per_taxon = defaultdict(Counter)            # gtdb taxon -> Counter(ncbi_taxid of genomes)
    for path in metadata_files:
        with gzip.open(path, 'rt') as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            ig, it = hdr.index('gtdb_taxonomy'), hdr.index('ncbi_taxid')
            for line in fh:
                f = line.rstrip('\n').split('\t')
                tid = merged.get(f[it], f[it])
                for t in f[ig].split(';'):
                    if len(t) > 3:
                        per_taxon[t][tid] += 1
    lineage_cache = {}

    def lineage(tid):
        if tid not in lineage_cache:
            out, t = [], tid
            while t in parent and t not in out:
                out.append(t)
                if parent[t] == t:
                    break
                t = parent[t]
            lineage_cache[tid] = out
        return lineage_cache[tid]
    n = 0
    with gzip.open(out_path, 'wt') as out:
        out.write('gtdb_taxon\tn_genomes\tncbi_taxids\n')
        for taxon, cnt in sorted(per_taxon.items()):
            total = sum(cnt.values())
            votes = Counter()
            for tid, k in cnt.items():
                for a in lineage(tid):
                    votes[a] += k
            kept = {a: k / total for a, k in votes.items() if k / total > min_frac}
            # do not go deeper than the GTDB rank itself (a GTDB phylum must not inherit an
            # NCBI class just because most of its genomes sit in one class)
            level = RANKS.index(taxon[:3])
            ok = [a for a in kept if RANK_LEVEL.get(rank.get(a), -1) <= level and RANK_LEVEL.get(rank.get(a), -1) >= 0]
            deepest = max(ok, key=lambda a: len(lineage(a)), default=None)
            chain = set(lineage(deepest)) if deepest else set()
            keep = sorted(((a, fr) for a, fr in kept.items() if a in chain), key=lambda x: -x[1])
            out.write(f'{taxon}\t{total}\t' + ','.join(f'{a}:{fr:.3f}' for a, fr in keep) + '\n')
            n += 1
    print(f'[gtdb2ncbi] {n:,} GTDB taxa written to {out_path}', file=log)


def load(path):
    out = {}
    with gzip.open(path, 'rt') as fh:
        next(fh)
        for line in fh:
            t, _n, ids = line.rstrip('\n').split('\t')
            out[t] = {'TAX-' + x.split(':')[0] for x in ids.split(',') if x}
    return out


def lowest_taxon(lineage_str):
    """'d__Bacteria;p__X;...;g__;s__' -> lowest named rank, e.g. 'p__X'."""
    best = None
    for t in lineage_str.split(';'):
        t = t.strip()
        if len(t) > 3 and t[:3] in RANKS:
            best = t
    return best


def ncbi_lineage(table, lineage_str):
    """Walk up from the lowest named GTDB rank to the first one the table knows."""
    ranks = [t.strip() for t in lineage_str.split(';') if len(t.strip()) > 3]
    for t in reversed(ranks):
        if t in table:
            return t, table[t]
    return None, set()


def read_gtdbtk(path):
    """GTDB-Tk summary (gtdbtk_taxonomy.tsv, or a *.summary.tsv): {bin: classification}.
    Bin names drop .fa/.fasta/.fna so they match DAS Tool's contig2bin names."""
    out = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ig = hdr.index('user_genome') if 'user_genome' in hdr else 0
        ic = hdr.index('classification') if 'classification' in hdr else 1
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) <= max(ig, ic):
                continue
            b = f[ig]
            for suf in ('.fasta', '.fa', '.fna'):
                if b.endswith(suf):
                    b = b[:-len(suf)]
            if f[ic] and not f[ic].startswith('Unclassified'):
                out[b] = f[ic]
    return out


def contig_lengths_gff(gff):
    out = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith('##sequence-region'):
                _, c, _s, e = line.split()
                out[c] = int(e)
            elif not line.startswith('#'):
                break
    return out


def read_sendsketch_bins(tsv, contig2bin, lengths, min_frac=0.5):
    """Fallback without GTDB-Tk. Per-bin GTDB-style lineage from SendSketch contig calls: descend
    rank by rank while one name holds > min_frac of the bin's total contig length."""
    calls = defaultdict(list)
    with open(tsv) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            b = contig2bin.get(f[0])
            if b and len(f) >= 5 and f[4]:
                calls[b].append((lengths.get(f[0], 1), f[4].split(';')))
    out = {}
    for b, cs in calls.items():
        total_bin = sum(lengths.get(c, 1) for c, bb in contig2bin.items() if bb == b)
        chosen = []
        for i in range(len(RANKS)):
            w = Counter()
            for L, lin in cs:
                if i < len(lin) and len(lin[i]) > 3 and lin[:i] == chosen:
                    w[lin[i]] += L
            if not w:
                break
            name, wt = w.most_common(1)[0]
            if wt / total_bin <= min_frac:
                break
            chosen.append(name)
        if chosen:
            out[b] = ';'.join(chosen)
    return out


