"""Protein -> Rhea -> MetaCyc reaction assignments, each labelled with its evidence class.

Evidence classes, strongest first:
  SEQ        DIAMOND hit to a Rhea-annotated UniProtKB protein (refs.py builds the reference:
             Swiss-Prot plus one Rhea-annotated member per UniRef90 cluster). Carries identity,
             coverage, the hit accession, whether it is Swiss-Prot, and the share of the hit's
             UniRef90 cluster members annotated with the same reaction.
  EC         complete EC number from Bakta or eggNOG -> Rhea (rhea2ec)
  KO         KO -> KEGG reaction -> Rhea, and the KO definition's [EC:] -> Rhea
  EC_DIRECT  EC -> MetaCyc reactions.dat EC-NUMBER, for MetaCyc reactions without a Rhea link
All routes are kept; downstream confidence uses the strongest class per reaction.

Memory: inputs are hash-partitioned by protein ID into chunks of about --chunk proteins, and
each chunk is resolved and written before the next is read. What stays in memory across
chunks is the reference maps, the DIAMOND subjects' Rhea annotations, and per-bin best
evidence per reaction.
"""
import gzip
import os
import re
import zlib
from collections import defaultdict

CLASSES = ('SEQ', 'EC', 'KO', 'EC_DIRECT')
EC_FULL = re.compile(r'^\d+\.\d+\.\d+\.\d+$')
DEFAULTS = dict(min_id=35.0, min_qcov=0.6, min_scov=0.6, min_support=0.5)
DIAMOND_FIELDS = 'qseqid sseqid pident length qlen slen evalue bitscore'
ASSIGN_HEADER = ['protein', 'contig', 'bin', 'metacyc_rxn', 'rhea_masters', 'evidence', 'all_evidence',
                 'seq_identity', 'seq_qcov', 'seq_hit', 'seq_hit_reviewed', 'seq_support',
                 'seq_cluster_agreement', 'ec', 'ko']


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def _rhea_tsv(path, key, val):
    out = defaultdict(set)
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            out[f[key]].add(f[val])
    return out


class References:
    """Open reference maps (Rhea, KEGG) plus MetaCyc's own Rhea and EC links."""

    def __init__(self, refs_dir, mc, ko_list=None):
        rd = os.path.join(refs_dir, 'rhea')
        self.ec2rhea = _rhea_tsv(f'{rd}/rhea2ec.tsv', 3, 2)
        self.rhea2ec = _rhea_tsv(f'{rd}/rhea2ec.tsv', 2, 3)
        self.kr2rhea = _rhea_tsv(f'{rd}/rhea2kegg_reaction.tsv', 3, 2)
        self.m2mc = _rhea_tsv(f'{rd}/rhea2metacyc.tsv', 2, 3)
        dirs = {}
        with open(f'{rd}/rhea-directions.tsv') as fh:
            next(fh)
            for line in fh:
                f = line.rstrip('\n').split('\t')
                for x in f:
                    dirs[x] = f[0]
        for rid, rxns in mc.rhea2rxn.items():         # MetaCyc's own Rhea links, folded to masters
            self.m2mc[dirs.get(rid, rid)] |= rxns
        self.ko2kr = defaultdict(set)
        p = os.path.join(refs_dir, 'kegg', 'kegg_ko2rn.tsv')
        if os.path.exists(p):
            with open(p) as fh:
                for line in fh:
                    k, r = line.rstrip('\n').split('\t')
                    self.ko2kr[k.replace('ko:', '')].add(r.replace('rn:', ''))
        self.ko_ec = defaultdict(set)
        if ko_list and os.path.exists(ko_list):     # KofamScan ko_list: definition is the last column
            with open(ko_list) as fh:
                for line in fh:
                    f = line.rstrip('\n').split('\t')
                    m = re.search(r'\[EC:([^\]]+)\]', f[-1])
                    if m:
                        self.ko_ec[f[0]].update(m.group(1).split())
        rhea_linked = set().union(*self.m2mc.values())
        self.ec2mc_direct = defaultdict(set)
        for r, es in mc.ec.items():
            if r not in rhea_linked:
                for e in es:
                    self.ec2mc_direct[e].add(r)


def load_seqref(path, wanted):
    """seqref_rhea.tsv.gz restricted to DIAMOND subjects -> {acc: (reviewed, {master: (n_with, n_rhea)})}"""
    out = {}
    with gzip.open(path, 'rt') as fh:
        next(fh)
        for line in fh:
            a, rest = line.split('\t', 1)
            if a not in wanted:
                continue
            _ur, rev, m, nr, nw = rest.rstrip('\n').split('\t')
            if a not in out:
                out[a] = (rev == '1', {})
            out[a][1][m] = (int(nw) if nw else None, int(nr) if nr else None)
    return out


def seq_calls(hits, seqref, o):
    """hits: iterable of DIAMOND tabular rows (lists). Hits passing identity and both coverages are
    pooled per protein; a Rhea reaction is accepted when the hits carrying it hold >= min_support
    of the pooled bitscore (DIAMOND --top keeps hits within 10% of each protein's best score)."""
    per = defaultdict(list)
    for q, s, pid, ln, ql, sl, _ev, bs in hits:
        pid, ln, ql, sl, bs = float(pid), int(ln), int(ql), int(sl), float(bs)
        if pid < o['min_id'] or ln / ql < o['min_qcov'] or ln / sl < o['min_scov'] or s not in seqref:
            continue
        per[q].append((s, pid, ln / ql, bs))
    out = {}
    for q, hs in per.items():
        tot = sum(h[3] for h in hs)
        by_m = defaultdict(list)
        for h in hs:
            for m in seqref[h[0]][1]:
                by_m[m].append(h)
        calls = {}
        for m, mh in by_m.items():
            support = sum(h[3] for h in mh) / tot
            if support < o['min_support']:
                continue
            best = max(mh, key=lambda h: h[3])
            rev, ann = seqref[best[0]]
            nw, nr = ann.get(m, (None, None))
            calls[m] = dict(hit=best[0], identity=round(best[1], 1), qcov=round(best[2], 2),
                            reviewed=int(rev), support=round(support, 3),
                            cluster_agreement='' if not nr else round((nw or 0) / nr, 3))
        if calls:
            out[q] = calls
    return out


def parse_bakta_rows(lines):
    """Bakta annotation.tsv rows -> {locus: dict(contig, ec, ko)}"""
    prot = {}
    for line in lines:
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[1] != 'cds':
            continue
        d = dict(contig=f[0], ec=set(), ko=set())
        for x in f[8].split(', '):
            x = x.strip()
            if x.startswith('EC:'):
                d['ec'].add(x[3:])
            elif x.startswith('KEGG:'):
                d['ko'].add(x[5:])
        prot[f[5]] = d
    return prot


def add_merged_rows(prot, lines, iko, iec):
    """merged_annotations.tsv rows: KO (KofamScan, eggNOG fallback) and EC (eggNOG)."""
    for line in lines:
        f = line.rstrip('\n').split('\t')
        d = prot.get(f[0])
        if d is None:
            continue
        for k in re.split('[,;]', f[iko]):
            k = k.strip().replace('ko:', '')
            if k.startswith('K'):
                d['ko'].add(k)
        for e in re.split('[,;]', f[iec]):
            e = e.strip()
            if e and e != '-':
                d['ec'].add(e)


def assign(prot, seq, refs):
    """One chunk's proteins -> merged assignment dicts (protein, metacyc_rxn, rhea, classes, best)."""
    ev = defaultdict(dict)                            # (protein, master) -> {class: detail}
    for p, d in prot.items():
        for m, info in seq.get(p, {}).items():
            ev[(p, m)]['SEQ'] = info
        for e in sorted(e for e in d['ec'] if EC_FULL.match(e)):
            ms = refs.ec2rhea.get(e, ())
            for m in ms:
                ev[(p, m)].setdefault('EC', dict(ec=[]))['ec'].append(e)
        for k in sorted(d['ko']):
            ms = set()
            for r in refs.ko2kr.get(k, ()):
                ms |= refs.kr2rhea.get(r, set())
            for e in refs.ko_ec.get(k, ()):
                if EC_FULL.match(e):
                    ms |= refs.ec2rhea.get(e, set())
            for m in ms:
                ev[(p, m)].setdefault('KO', dict(ko=[]))['ko'].append(k)
    merged = {}

    def add(p, rxn, master, classes):
        key = (p, rxn or 'RHEA:' + master)
        t = merged.setdefault(key, dict(protein=p, metacyc_rxn=rxn, rhea=set(), classes={}))
        if master:
            t['rhea'].add(master)
        for c, det in classes.items():
            have = t['classes'].get(c)
            if have is None:
                t['classes'][c] = dict(det)
            elif c == 'SEQ':
                if (det['identity'], have['hit']) > (have['identity'], det['hit']):   # best identity, then smaller accession
                    t['classes'][c] = dict(det)
            else:
                for k in ('ec', 'ko'):
                    if k in det:
                        have[k] = sorted(set(have[k]) | set(det[k]))

    for (p, m), cls in ev.items():
        rxns = refs.m2mc.get(m)
        for r in sorted(rxns or ()):
            add(p, r, m, cls)
        if not rxns:
            add(p, '', m, cls)
    for p, d in prot.items():
        ecs = {e for e in d['ec'] if EC_FULL.match(e)}
        for m in seq.get(p, {}):                       # ECs implied by sequence-supported Rhea reactions
            ecs |= refs.rhea2ec.get(m, set())
        for e in sorted(ecs):
            for r in sorted(refs.ec2mc_direct.get(e, ())):
                add(p, r, '', {'EC_DIRECT': dict(ec=[e])})
    for t in merged.values():
        t['best'] = next(c for c in CLASSES if c in t['classes'])
    return list(merged.values())


def assign_row(t, prot, c2b):
    p = t['protein']
    s = t['classes'].get('SEQ', {})
    ecs = sorted(set(t['classes'].get('EC', {}).get('ec', []) + t['classes'].get('EC_DIRECT', {}).get('ec', [])))
    return [p, prot[p]['contig'], c2b.get(prot[p]['contig'], ''), t['metacyc_rxn'], ';'.join(sorted(t['rhea'])),
            t['best'], ','.join(c for c in CLASSES if c in t['classes']),
            s.get('identity', ''), s.get('qcov', ''), s.get('hit', ''), s.get('reviewed', ''),
            s.get('support', ''), s.get('cluster_agreement', ''), ','.join(ecs),
            ','.join(sorted(t['classes'].get('KO', {}).get('ko', [])))]


def _part(key, k):
    return zlib.crc32(key.encode()) % k


def partition(bakta_tsv, merged_tsv, hits_tsv, k, tmp):
    """Split the three per-protein inputs into k hash partitions; returns (paths, DIAMOND subjects, header info)."""
    os.makedirs(tmp, exist_ok=True)
    names = ('bakta', 'merged', 'hits')
    fhs = {(n, i): open(os.path.join(tmp, f'{n}.{i}.tsv'), 'w') for n in names for i in range(k)}
    subjects = set()
    with _open(bakta_tsv) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 6)
            if len(f) > 5:
                fhs[('bakta', _part(f[5], k))].write(line)
    merged_cols = None
    if merged_tsv:
        with _open(merged_tsv) as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            merged_cols = (hdr.index('KO'), hdr.index('EC'))
            for line in fh:
                fhs[('merged', _part(line.split('\t', 1)[0], k))].write(line)
    with _open(hits_tsv) as fh:
        for line in fh:
            f = line.split('\t', 2)
            if len(f) > 2:
                subjects.add(f[1])
                fhs[('hits', _part(f[0], k))].write(line)
    for h in fhs.values():
        h.close()
    return [{n: os.path.join(tmp, f'{n}.{i}.tsv') for n in names} for i in range(k)], subjects, merged_cols


def iter_chunks(parts, merged_cols, seqref, refs, opts=None):
    """Yield (prot, assignments) per partition, deleting each partition's files once read."""
    o = dict(DEFAULTS, **(opts or {}))
    for paths in parts:
        with open(paths['bakta']) as fh:
            prot = parse_bakta_rows(fh)
        if merged_cols:
            with open(paths['merged']) as fh:
                add_merged_rows(prot, fh, *merged_cols)
        with open(paths['hits']) as fh:
            seq = seq_calls((line.rstrip('\n').split('\t') for line in fh), seqref, o)
        for p in paths.values():
            os.remove(p)
        yield prot, assign(prot, seq, refs)
