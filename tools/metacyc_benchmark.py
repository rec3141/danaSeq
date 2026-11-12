#!/usr/bin/env python3
"""Benchmark the MetaCyc pathway scorer (mag_analysis/bin/metacyc_scorer) against curated calls.

Ground truth is a curated organism PGDB exported to flat files, e.g. EcoCyc for E. coli K-12.
Nothing licensed is stored here: the PGDB export, MetaCyc and every pathway-level output stay
wherever the user keeps them (owner-only), and the tool prints only counts unless --detail is
given.

Subcommands
  export   Export a PGDB from an installed Pathway Tools to flat files, the way
           `metacyc_scorer.py setup` exports MetaCyc (batch-mode Lisp, xvfb).
             metacyc_benchmark.py export --ptools PT_INSTALL_DIR --org ecoli --out DIR
           PT_INSTALL_DIR comes from `metacyc_scorer.py setup --installer X --keep-install
           --work-dir W` (it is W/pt). The tier-1 bundle holds EcoCyc (org ID ecoli) and MetaCyc.
  compare  Precision and recall of one bin's calls against the curated pathways.
             metacyc_benchmark.py compare --truth DIR/pathways.dat --metacyc METACYC_DIR \\
                 --refs REFS_DIR --pathways metacyc_pathways.tsv --bin BIN [--detail out.tsv]
  degrade  Copy scorer inputs with random blocks of genes removed (held-out incompleteness test).
             metacyc_benchmark.py degrade --bakta-tsv T --merged M --hits H --drop 0.4 --seed 1 --out DIR
  species  Recall-only truth from MetaCyc itself: the pathways whose SPECIES slot (organisms with
           experimental evidence for the pathway) names one of the given taxa. Writes an ID-only
           pathways.dat for compare; precision against it is not meaningful.
             metacyc_benchmark.py species --metacyc METACYC_DIR --taxa TAX-1423,TAX-224308 --out truth.dat
  agree    Agreement of scorer calls with PathoLogic's (*_pwy.tsv from Pathway Tools / MetaPathways),
           per tier, summed over bins (the community is reported on its own).
             metacyc_benchmark.py agree --pathways metacyc_pathways.tsv --pathologic 'GLOB'

Definitions used by compare
  truth          curated pathways that also exist in MetaCyc (pathways only in the curated PGDB
                 cannot be predicted and are listed as 'not_in_metacyc')
  unreachable    truth pathways the scorer cannot reach from any protein: no scored reactions,
                 or none of them (or of their composite parts) has a Rhea link or an EC number
  base / super   super-pathways (MetaCyc SUB-PATHWAYS) are reported separately; they are called
                 from their components, so they mostly double-count base-pathway results
  precision      calls in the curated set / calls (per tier and cumulative high, high+medium, all)
  recall         calls in the curated set / reachable truth pathways
"""
import argparse
import csv
import glob
import gzip
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'mag_analysis', 'bin'))
from metacyc_scorer import ptdat  # noqa: E402

TIERS = ('high', 'medium', 'low')


def _ids(path):
    return {r['UNIQUE-ID']: r for r in ptdat.iter_records(path) if 'UNIQUE-ID' in r}


def linkable_reactions(mc, refs_dir=None):
    """MetaCyc reactions that some protein evidence can reach: a Rhea link (MetaCyc's own DBLINKS
    or Rhea's rhea2metacyc) or an EC number; composites through any linkable part."""
    linked = set().union(*mc.rhea2rxn.values()) if mc.rhea2rxn else set()
    linked |= {r for r, ecs in mc.ec.items() if ecs}
    if refs_dir:
        p = os.path.join(refs_dir, 'rhea', 'rhea2metacyc.tsv')
        if os.path.exists(p):
            with open(p) as fh:
                next(fh)
                linked |= {line.rstrip('\n').split('\t')[3] for line in fh}
    for r, parts in mc.composite.items():
        if parts & linked:
            linked.add(r)
    return linked


def truth_sets(truth_pathways, mc, linked):
    """-> dict of pathway-ID sets: truth (in MetaCyc), base, super, reachable, and the
    unreachable reasons."""
    cur = _ids(truth_pathways)
    out = dict(curated=set(cur), not_in_metacyc=set(cur) - set(mc.pw))
    t = set(cur) & set(mc.pw)
    out['super'] = {p for p in t if p in mc.is_super}
    out['base'] = t - out['super']
    no_rxn, no_link = set(), set()
    for p in t:
        rx = mc.base.get(p, set())
        scored = (rx - mc.hypothetical.get(p, set()) or rx) - mc.spontaneous
        if not scored:
            no_rxn.add(p)
        elif not scored & linked:
            no_link.add(p)
    out['no_reactions'] = no_rxn
    out['no_reaction_links'] = no_link
    out['reachable'] = t - no_rxn - no_link
    out['truth'] = t
    return out


def read_calls(path, bin_name):
    """{pathway: (tier, row)} for rows of bin_name with call == 1."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['bin'] == bin_name and r['call'] == '1':
                out[r['pathway']] = (r['confidence_tier'], r)
    return out


def metrics(calls, ts, subset):
    """calls {pw: (tier, row)}; subset 'base' or 'super'. Per tier and cumulative."""
    truth = ts[subset]
    reach = truth & ts['reachable']
    c = {p: t for p, (t, row) in calls.items() if (row.get('super') == '1') == (subset == 'super')}
    res = dict(subset=subset, n_truth=len(truth), n_reachable=len(reach))
    cum = set()
    for t in TIERS:
        sel = {p for p, x in c.items() if x == t}
        tp = len(sel & truth)
        cum |= sel
        ctp = len(cum & truth)
        res[t] = dict(n=len(sel), tp=tp, precision=_div(tp, len(sel)), recall=_div(tp, len(reach)),
                      cum_n=len(cum), cum_tp=ctp, cum_precision=_div(ctp, len(cum)),
                      cum_recall=_div(ctp, len(reach)))
    return res


def _div(a, b):
    return round(a / b, 4) if b else None


def cmd_compare(a):
    mc = ptdat.MetaCyc(ptdat.load(os.path.join(a.metacyc, 'pathways.dat')),
                       ptdat.load(os.path.join(a.metacyc, 'reactions.dat')),
                       ptdat.load(os.path.join(a.metacyc, 'classes.dat')))
    ts = truth_sets(a.truth, mc, linkable_reactions(mc, a.refs))
    calls = read_calls(a.pathways, a.bin)
    res = dict(bin=a.bin, n_curated=len(ts['curated']), n_not_in_metacyc=len(ts['not_in_metacyc']),
               n_no_reactions=len(ts['no_reactions']), n_no_reaction_links=len(ts['no_reaction_links']),
               base=metrics(calls, ts, 'base'), super=metrics(calls, ts, 'super'))
    if a.json:
        print(json.dumps(res, indent=1))
    else:
        print(format_result(res))
    if a.detail:
        write_detail(a.detail, a.pathways, a.bin, ts)
    return res


def format_result(res):
    lines = [f"bin {res['bin']}: {res['n_curated']} curated pathways; {res['n_not_in_metacyc']} not in "
             f"MetaCyc, {res['n_no_reactions']} without scored reactions, {res['n_no_reaction_links']} "
             f"without Rhea/EC-linked reactions"]
    for s in ('base', 'super'):
        m = res[s]
        lines.append(f"  {s}: truth {m['n_truth']}, reachable {m['n_reachable']}")
        lines.append('    tier     calls   TP  precision  recall | cumulative: calls   TP  precision  recall')
        for t in TIERS:
            x = m[t]
            lines.append(f"    {t:7s} {x['n']:6d} {x['tp']:4d}  {_f(x['precision'])}  {_f(x['recall'])}  |"
                         f"             {x['cum_n']:5d} {x['cum_tp']:4d}  {_f(x['cum_precision'])}  "
                         f"{_f(x['cum_recall'])}")
    return '\n'.join(lines)


def _f(x):
    return '   -    ' if x is None else f'{x:8.3f}'


def write_detail(out, pathways_tsv, bin_name, ts):
    """Per-pathway table (licensed content: keep it private)."""
    seen = set()
    with open(pathways_tsv) as fh, open(out, 'w') as o:
        rd = csv.DictReader(fh, delimiter='\t')
        cols = ['pathway', 'super', 'in_truth', 'reachable', 'call', 'confidence_tier', 'confidence',
                'reason', 'n_rxn', 'n_present', 'frac', 'n_seq', 'n_enzymatic', 'key_seq', 'seq_id_median',
                'n_ec_only', 'n_ko_only', 'n_ec_direct_only', 'n_linked_only']
        o.write('\t'.join(cols) + '\n')
        for r in rd:
            if r['bin'] != bin_name:
                continue
            seen.add(r['pathway'])
            r['in_truth'] = int(r['pathway'] in ts['truth'])
            r['reachable'] = int(r['pathway'] in ts['reachable'])
            o.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')
        for p in sorted(ts['truth'] - seen):
            why = ('no_reactions' if p in ts['no_reactions'] else
                   'no_reaction_links' if p in ts['no_reaction_links'] else 'no_evidence')
            o.write('\t'.join([p, str(int(p in ts['super'])), '1', str(int(p in ts['reachable'])), '0', '', '',
                               why] + [''] * (len(cols) - 8)) + '\n')
        for p in sorted(ts['not_in_metacyc']):
            o.write('\t'.join([p, '', '1', '0', '0', '', '', 'not_in_metacyc'] + [''] * (len(cols) - 8)) + '\n')
    os.chmod(out, 0o600)


# ---------------------------------------------------------------------------------------- degrade

def _open(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def cmd_degrade(a):
    """Remove a fraction of CDS in contiguous blocks (genome order), so whole operons go missing
    together as they do in an incomplete MAG; the same loci are removed from all three inputs."""
    with open(a.bakta_tsv) as fh:
        cds = [line.rstrip('\n').split('\t')[5] for line in fh
               if not line.startswith('#') and line.split('\t')[1:2] == ['cds']]
    n = len(cds)
    n_drop = round(n * a.drop)
    rng = random.Random(a.seed)
    k = max(1, round(n_drop / a.block))
    sizes = [n_drop // k + (i < n_drop % k) for i in range(k)]
    drop = [False] * n
    for s in sizes:                       # blocks placed at random, no overlap, circular genome
        for _ in range(1000):
            st = rng.randrange(n)
            idx = [(st + i) % n for i in range(s)]
            if not any(drop[i] for i in idx):
                for i in idx:
                    drop[i] = True
                break
    gone = {cds[i] for i in range(n) if drop[i]}
    os.makedirs(a.out, exist_ok=True)
    with open(a.bakta_tsv) as fh, open(os.path.join(a.out, 'annotation.tsv'), 'w') as o:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(f) < 6 or f[5] not in gone:
                o.write(line)
    with open(a.merged) as fh, open(os.path.join(a.out, 'merged_annotations.tsv'), 'w') as o:
        o.write(fh.readline())
        for line in fh:
            if line.split('\t', 1)[0] not in gone:
                o.write(line)
    with _open(a.hits) as fh, gzip.open(os.path.join(a.out, 'metacyc_hits.tsv.gz'), 'wt') as o:
        for line in fh:
            if line.split('\t', 1)[0] not in gone:
                o.write(line)
    print(f'{len(gone)} of {n} CDS removed in {k} blocks (seed {a.seed}) -> {a.out}', file=sys.stderr)


# ------------------------------------------------------------------------------------------ agree

def cmd_agree(a):
    pl = {}
    for f in sorted(glob.glob(a.pathologic, recursive=True)):
        ent = os.path.basename(f)[:-len('_pwy.tsv')]
        ent = '_community' if ent == a.community_name else ent
        with open(f) as fh:
            pl[ent] = {r['PWY_NAME'] for r in csv.DictReader(fh, delimiter='\t')}
    sc = {}
    with open(a.pathways) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['call'] == '1':
                sc.setdefault(r['bin'], {})[r['pathway']] = r['confidence_tier']
    tot = {t: [0, 0] for t in TIERS}
    n_pl = n_both = n_sc = 0
    for ent in sorted(pl):
        if ent == '_community':
            continue
        s = sc.get(ent, {})
        n_pl += len(pl[ent])
        n_sc += len(s)
        n_both += len(pl[ent] & set(s))
        for p, t in s.items():
            tot[t][0] += p in pl[ent]
            tot[t][1] += 1
    res = dict(bins=len([e for e in pl if e != '_community']), pathologic=n_pl, scorer=n_sc, both=n_both,
               recovered=_div(n_both, n_pl),
               tiers={t: dict(agree=v[0], n=v[1], frac=_div(v[0], v[1])) for t, v in tot.items()})
    print(json.dumps(res, indent=1))
    return res


def cmd_species(a):
    want = set(a.taxa.split(','))
    n = 0
    with open(a.out, 'w') as o:
        for r in ptdat.iter_records(os.path.join(a.metacyc, 'pathways.dat')):
            if 'UNIQUE-ID' in r and want & set(r.get('SPECIES', [])):
                o.write(f"UNIQUE-ID - {r['UNIQUE-ID']}\n//\n")
                n += 1
    os.chmod(a.out, 0o600)
    print(f'{n} pathways with experimental evidence in {a.taxa} -> {a.out}', file=sys.stderr)


def cmd_export(a):
    from metacyc_scorer import install
    pt = os.path.realpath(a.ptools)
    work = os.path.dirname(pt)
    os.makedirs(os.path.join(work, 'home'), exist_ok=True)
    install._export(install.PT_BASE if a.runtime is None else a.runtime, [work], f'{pt}/pathway-tools',
                    None, work, os.path.join(work, f'export_{a.org}.log'), org=a.org)
    data = sorted(glob.glob(f'{pt}/aic-export/pgdbs/biocyc/*/*/data/pathways.dat'), key=os.path.getmtime)
    if not data:
        sys.exit('no exported pathways.dat found')
    src = os.path.dirname(data[-1])
    os.makedirs(a.out, exist_ok=True)
    os.chmod(a.out, 0o700)
    for name in ('pathways.dat', 'reactions.dat', 'classes.dat', 'genes.dat', 'proteins.dat', 'enzrxns.dat'):
        p = os.path.join(src, name)
        if os.path.exists(p):
            dst = os.path.join(a.out, name)
            with open(p, 'rb') as i, open(dst, 'wb') as o:
                o.write(i.read())
            os.chmod(dst, 0o600)
    print(f'{src} -> {a.out}', file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd')
    sub.required = True
    e = sub.add_parser('export')
    e.add_argument('--ptools', required=True, help='Pathway Tools install dir (setup --keep-install: WORK/pt)')
    e.add_argument('--org', default='ecoli', help='PGDB org ID [ecoli = EcoCyc]')
    e.add_argument('--runtime', help='Apptainer sandbox with the Pathway Tools runtime [MetaPathways ptools-base]')
    e.add_argument('--out', required=True)
    c = sub.add_parser('compare')
    c.add_argument('--truth', required=True, help='curated PGDB pathways.dat (e.g. EcoCyc export)')
    c.add_argument('--metacyc', required=True, help='MetaCyc dir the scorer used')
    c.add_argument('--refs', help='open references dir (rhea2metacyc for reachability)')
    c.add_argument('--pathways', required=True, help='scorer metacyc_pathways.tsv')
    c.add_argument('--bin', required=True)
    c.add_argument('--detail', help='per-pathway TSV (licensed content; written owner-only)')
    c.add_argument('--json', action='store_true')
    d = sub.add_parser('degrade')
    d.add_argument('--bakta-tsv', required=True)
    d.add_argument('--merged', required=True)
    d.add_argument('--hits', required=True)
    d.add_argument('--drop', type=float, required=True, help='fraction of CDS to remove')
    d.add_argument('--block', type=int, default=50, help='genes per removed block [50]')
    d.add_argument('--seed', type=int, default=1)
    d.add_argument('--out', required=True)
    sp = sub.add_parser('species')
    sp.add_argument('--metacyc', required=True)
    sp.add_argument('--taxa', required=True, help='comma-separated MetaCyc TAX- frames')
    sp.add_argument('--out', required=True)
    g = sub.add_parser('agree')
    g.add_argument('--pathways', required=True)
    g.add_argument('--pathologic', required=True, help="glob of *_pwy.tsv files, quoted ('**' allowed)")
    g.add_argument('--community-name', default='assembly', help='entity name of the community PGDB')
    a = ap.parse_args(argv)
    {'compare': cmd_compare, 'degrade': cmd_degrade, 'agree': cmd_agree, 'export': cmd_export,
     'species': cmd_species}[a.cmd](a)


if __name__ == '__main__':
    main()
