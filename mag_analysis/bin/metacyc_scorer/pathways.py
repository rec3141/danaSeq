"""Stage 2: MetaCyc reactions -> pathways.

mode 'pathologic' (default) follows the published PathoLogic pathway-inference rules
(Karp, Latendresse & Caspi 2011, Stand Genomic Sci 5:424, doi:10.4056/sigs.1794338):
a base pathway is inferred unless a rule rejects it.
  hard rejections
    taxonomic range   outside the range -> inferred only if every reaction is present
    electron transfer inferred only if complete
    key reactions     any key reaction missing -> rejected; all present -> inferred
    mostly absent     > 1 reaction missing and <= 1 present -> rejected
    energy            energy-metabolism pathway with no unique reaction and < half present -> rejected
  soft rejections (a present reaction unique to the pathway among candidate pathways overrides)
    subset            biosynthetic pathway that is a proper subset of another pathway and
                      misses its final reactions; degradative one that misses its initial step
    superset          reactions are a superset of an inferred pathway's and more are missing
  super-pathways      inferred iff every component pathway is inferred (recursively)
mode 'frac' is a simple fraction threshold kept for comparison and for ML features
(Dale, Popescu & Karp 2010, BMC Bioinformatics 11:15 use similar features).

This re-implements the published rules; SRI's code has undocumented details (e.g. how ties
between unique-reaction pathways are broken), so calls will not match PathoLogic exactly.
"""
from collections import defaultdict
from itertools import combinations

DEFAULTS = dict(mode='pathologic', min_frac=0.5, composite_min=2 / 3)


def infer_present(mc, enzymatic, composite_min=2 / 3):
    """Reactions present = enzymatic evidence + spontaneous + composite reactions whose
    components are (mostly) present. Returns (present, enzymatic_or_composite)."""
    enz = set(enzymatic)
    for r, parts in mc.composite.items():
        if r in enz or not parts:
            continue
        have = len(parts & (enz | mc.spontaneous))
        if have == len(parts) or (len(parts) >= 3 and have / len(parts) >= composite_min):
            enz.add(r)
    return enz | mc.spontaneous, enz


def _tax_compat(mc, tr, lineage):
    """1 = organism lineage inside the pathway's range; 'U' = range lies below the deepest
    resolved rank of the organism (e.g. MAG resolved only to Bacteria, range is
    Pseudomonadota), so it cannot be excluded; 0 = disjoint."""
    if tr & lineage:
        return 1
    known = {t for t in lineage if t in mc.cls}
    for t in tr:
        if known and known <= mc.class_ancestors(t):
            return 'U'
    return 0


def score_pathways(mc, present_counts, lineage=None, opts=None):
    """mc: ptdat.MetaCyc. present_counts: {rxn: n_proteins} with enzyme evidence.
    lineage: set of 'TAX-<id>' for the organism (None disables taxonomic pruning).
    Returns {pathway: row}; row['call'] 1/0 and row['reason'] the deciding rule."""
    o = dict(DEFAULTS, **(opts or {}))
    present, enz = infer_present(mc, {r for r, n in present_counts.items() if n}, o['composite_min'])
    base = [p for p in mc.base if p not in mc.is_super and mc.base[p]]
    cand = {p for p in base if mc.base[p] & enz}
    # uniqueness among candidate pathways ("used in only one pathway in the generated PGDB")
    used_in = defaultdict(set)
    for p in cand:
        for r in mc.base[p] & present:
            used_in[r].add(p)
    rows = {}
    for p in base:
        rxns = mc.base[p]
        scored = rxns - mc.hypothetical.get(p, set()) or rxns
        pres = scored & present
        pe = scored & enz
        missing = scored - present
        uniq_global = {r for r in scored if len(mc.rxn2pw.get(r, ())) <= 1}
        uniq_present = {r for r in pres if used_in.get(r) == {p}}
        key = mc.key.get(p, set())
        tr = mc.taxrange.get(p, set())
        in_range = _tax_compat(mc, tr, lineage) if (lineage and tr) else None
        classes = mc.pathway_classes(p)
        rows[p] = dict(
            pathway=p, name=mc.name.get(p, p), super=0, n_rxn=len(scored), n_present=len(pres),
            n_enzymatic=len(pe), n_spontaneous=len(pres - pe), n_missing=len(missing),
            frac=round(len(pres) / len(scored), 4), n_unique=len(uniq_global),
            n_unique_present=len(uniq_present), n_key=len(key), n_key_present=len(key & present),
            tax_range=';'.join(sorted(tr)), in_tax_range='' if in_range is None else in_range,
            klass=_top_class(classes),
            present_rxns=';'.join(sorted(pres)), missing_rxns=';'.join(sorted(missing)),
            _sets=dict(present=pres, enzymatic=pe, key=key, unique_present=uniq_present, scored=scored))
        if o['mode'] == 'frac':
            rows[p]['call'], rows[p]['reason'] = _decide_frac(rows[p], o, in_range)
        else:
            rows[p]['call'], rows[p]['reason'] = _decide_pathologic(mc, p, rows[p], pres, missing,
                                                                     in_range, classes, cand)
    if o['mode'] != 'frac':
        _superset_pass(mc, rows, base)
    _super_pass(mc, rows, present, enz)
    return rows


def _top_class(classes):
    for c in ('Electron-Transfer', 'Energy-Metabolism', 'Biosynthesis', 'Degradation'):
        if c in classes:
            return c
    return ''


def _decide_frac(r, o, in_range):
    if not r['n_enzymatic']:
        return 0, 'no_enzymatic_evidence'
    if r['n_key'] and r['n_key_present'] < r['n_key']:
        return 0, 'key_reaction_missing'
    if in_range == 0 and r['n_missing']:
        return 0, 'outside_taxonomic_range'
    if r['frac'] >= o['min_frac']:
        return 1, 'frac>=%.2f' % o['min_frac']
    return 0, 'frac<%.2f' % o['min_frac']


def _decide_pathologic(mc, p, r, pres, missing, in_range, classes, cand):
    if p not in cand:
        return 0, 'no_enzymatic_evidence'
    if in_range == 0 and missing:
        return 0, 'outside_taxonomic_range'
    if 'Electron-Transfer' in classes and missing:
        return 0, 'electron_transfer_incomplete'
    if r['n_key']:
        if r['n_key_present'] < r['n_key']:
            return 0, 'key_reaction_missing'
        return 1, 'all_key_reactions_present'
    if len(missing) > 1 and len(pres) <= 1:
        return 0, 'mostly_absent'
    if 'Energy-Metabolism' in classes and not r['n_unique_present'] and len(pres) < len(missing):
        return 0, 'energy_no_unique_lt_half'
    soft = None
    if mc.is_proper_subset(p):
        if 'Biosynthesis' in classes and mc.final_rxns(p) and mc.final_rxns(p) <= missing:
            soft = 'biosynthetic_subset_missing_final'
        elif 'Degradation' in classes and mc.initial_rxns(p) and mc.initial_rxns(p) <= missing:
            soft = 'degradative_subset_missing_initial'
    if r['n_unique_present']:
        return 1, 'unique_reaction_present'
    if soft:
        return 0, soft
    return 1, 'no_rejection_rule'


def _superset_pass(mc, rows, base):
    """P1 is rejected if its reactions are a superset of an inferred P2 and P1 misses more.
    Judged against the first-pass calls so pathway order does not matter."""
    inferred = [p for p in base if rows[p]['call']]
    drop = {}
    for p1 in inferred:
        r1 = rows[p1]
        if r1['reason'] in ('all_key_reactions_present', 'unique_reaction_present'):
            continue
        s1 = mc.base[p1]
        for p2 in inferred:
            if p2 != p1 and mc.base[p2] < s1 and rows[p2]['n_missing'] < r1['n_missing']:
                drop[p1] = 'superset_of_%s_more_missing' % p2
                break
    for p, why in drop.items():
        rows[p]['call'], rows[p]['reason'] = 0, why


def _super_pass(mc, rows, present, enz):
    """Super-pathway inferred iff all components are inferred; resolved bottom-up."""
    memo = {}

    def inferred(p, stack=()):
        if p in memo:
            return memo[p]
        if p not in mc.is_super:
            memo[p] = bool(rows.get(p, {}).get('call'))
            return memo[p]
        if p in stack:
            return False
        subs = [s for s in mc.pw[p].get('SUB-PATHWAYS', []) if s in mc.pw]
        memo[p] = bool(subs) and all(inferred(s, stack + (p,)) for s in subs)
        return memo[p]

    for p in mc.is_super:
        rx = mc.base[p]
        if not rx & enz:
            continue
        subs = [s for s in mc.pw[p].get('SUB-PATHWAYS', []) if s in mc.pw]
        ok = inferred(p)
        n_ok = sum(inferred(s) for s in subs)
        pres = rx & present
        rows[p] = dict(pathway=p, name=mc.name.get(p, p), super=1, n_rxn=len(rx), n_present=len(pres),
                       n_enzymatic=len(rx & enz), n_spontaneous=len(pres - enz), n_missing=len(rx - present),
                       frac=round(len(pres) / len(rx), 4), n_unique='', n_unique_present='',
                       n_key='', n_key_present='', tax_range='', in_tax_range='', klass='Super-Pathways',
                       present_rxns=';'.join(sorted(pres)), missing_rxns=';'.join(sorted(rx - present)),
                       call=int(ok), reason='all_components_inferred' if ok else
                       'components_inferred_%d_of_%d' % (n_ok, len(subs)),
                       _sets=dict(present=pres, enzymatic=rx & enz, key=set(), unique_present=set(), scored=rx))


def minpath(mc, present_counts, max_exact=18, composite_min=2 / 3):
    """MinPath-style parsimony: smallest set of base pathways that explains every present
    reaction belonging to any pathway. Exact search per connected component when the
    component has <= max_exact candidate pathways, greedy set cover otherwise."""
    _, enz = infer_present(mc, {r for r, n in present_counts.items() if n}, composite_min)
    cand = defaultdict(set)
    for r in enz:
        for p in mc.rxn2pw.get(r, ()):
            cand[p].add(r)
    parent = {p: p for p in cand}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    owner = {}
    for p, rs in sorted(cand.items()):
        for r in sorted(rs):
            if r in owner:
                parent[find(p)] = find(owner[r])
            else:
                owner[r] = p
    comps = defaultdict(list)
    for p in sorted(cand):                   # sorted: the chosen cover must not depend on set order
        comps[find(p)].append(p)
    selected, methods = set(), defaultdict(int)
    for ps in comps.values():
        universe = set().union(*(cand[p] for p in ps))
        if len(ps) <= max_exact:
            selected |= _exact_cover(ps, cand, universe)
            methods['exact'] += 1
        else:
            selected |= _greedy_cover(ps, cand, universe)
            methods['greedy'] += 1
    return selected, dict(methods)


def _exact_cover(ps, cand, universe):
    forced = {p for p in ps if any(sum(r in cand[q] for q in ps) == 1 for r in cand[p])}
    left = universe - set().union(*(cand[p] for p in forced)) if forced else set(universe)
    if not left:
        return forced
    rest = sorted(p for p in ps if p not in forced and cand[p] & left)
    for k in range(1, len(rest) + 1):
        for combo in combinations(rest, k):
            if left <= set().union(*(cand[p] for p in combo)):
                return forced | set(combo)
    return forced | set(rest)


def _greedy_cover(ps, cand, universe):
    left, chosen = set(universe), []
    while left:
        best = max(ps, key=lambda p: (len(cand[p] & left), -len(cand[p]), p))
        chosen.append(best)
        left -= cand[best]
    for p in list(chosen):
        others = set().union(*(cand[q] for q in chosen if q != p)) if len(chosen) > 1 else set()
        if universe <= others:
            chosen.remove(p)
    return set(chosen)
