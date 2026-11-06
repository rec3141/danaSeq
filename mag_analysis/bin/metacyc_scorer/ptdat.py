"""Parser for Pathway Tools attribute-value flat files (pathways.dat, reactions.dat, classes.dat).

Format (http://bioinformatics.ai.sri.com/ptools/flatfile-format.html):
  '#' lines are comments; records end with a line '//'.
  'ATTRIBUTE - value' lines; an attribute may repeat (multi-valued slot).
  Lines starting with '/' continue the previous value (joined with a newline).
  '^ANNOT - value' lines annotate the preceding attribute value.
Each record becomes {'UNIQUE-ID': str, attr: [values], '_annot': {(attr, idx): {annot: [values]}}}.
"""
import re

_LINE = re.compile(r'^(\^?[A-Za-z0-9?_+.-]+) - ?(.*)$')


def iter_records(path):
    rec = None
    last = None              # (attr, index) of the most recent plain value
    with open(path, encoding='latin-1') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            if line == '//':
                if rec:
                    yield rec
                rec, last = None, None
                continue
            if rec is None:
                rec = {'_annot': {}}
            if line.startswith('/') and last is not None:
                attr, i = last
                if attr.startswith('^'):
                    base, annot = attr[1:].split('|', 1)
                    vals = rec['_annot'][(base, i)][annot]
                    vals[-1] += '\n' + line[1:]
                else:
                    rec[attr][i] += '\n' + line[1:]
                continue
            m = _LINE.match(line)
            if not m:
                continue
            attr, val = m.group(1), m.group(2)
            if attr.startswith('^'):
                if last is None:
                    continue
                base, i = last
                if base.startswith('^'):
                    base = base[1:].split('|', 1)[0]
                rec['_annot'].setdefault((base, i), {}).setdefault(attr[1:], []).append(val)
                last = ('^' + base + '|' + attr[1:], i)
                continue
            if attr == 'UNIQUE-ID':
                rec['UNIQUE-ID'] = val
                last = None
                continue
            vals = rec.setdefault(attr, [])
            vals.append(val)
            last = (attr, len(vals) - 1)
    if rec:
        yield rec


def load(path):
    return {r['UNIQUE-ID']: r for r in iter_records(path) if 'UNIQUE-ID' in r}


_TOK = re.compile(r'\(|\)|"[^"]*"|[^\s()]+')


def parse_sexpr(s):
    """Minimal Lisp s-expression reader for REACTION-LAYOUT / PREDECESSORS values."""
    toks = _TOK.findall(s)
    pos = 0

    def read():
        nonlocal pos
        t = toks[pos]
        pos += 1
        if t == '(':
            out = []
            while toks[pos] != ')':
                out.append(read())
            pos += 1
            return out
        return t.strip('"')
    return read() if toks else None


class MetaCyc:
    """Pathway/reaction model needed for scoring."""

    def __init__(self, pathways, reactions, classes=None):
        self.pw = pathways
        self.rxn = reactions
        self.cls = classes or {}
        self.spontaneous = {r for r, rec in reactions.items() if rec.get('SPONTANEOUS?', ['NIL'])[0] == 'T'}
        self.base = {}             # pathway -> set(reactions) (super-pathways expanded)
        for p in pathways:
            self.base[p] = self._expand(p, set())
        self.is_super = {p for p, rec in pathways.items() if rec.get('SUB-PATHWAYS')}
        self.key = {p: set(rec.get('KEY-REACTIONS', [])) for p, rec in pathways.items()}
        self.hypothetical = {p: set(rec.get('HYPOTHETICAL-REACTIONS', [])) for p, rec in pathways.items()}
        self.taxrange = {p: set(rec.get('TAXONOMIC-RANGE', [])) for p, rec in pathways.items()}
        self.name = {p: _clean(rec.get('COMMON-NAME', [p])[0]) for p, rec in pathways.items()}
        self.types = {p: rec.get('TYPES', []) for p, rec in pathways.items()}
        # reaction -> base (non-super) pathways containing it
        self.rxn2pw = {}
        for p, rs in self.base.items():
            if p in self.is_super:
                continue
            for r in rs:
                self.rxn2pw.setdefault(r, set()).add(p)
        self.rhea2rxn = {}
        for r, rec in reactions.items():
            for link in rec.get('DBLINKS', []):
                m = re.match(r'\(RHEA "(\d+)"', link)
                if m:
                    self.rhea2rxn.setdefault(m.group(1), set()).add(r)
        self.ec = {r: [e.replace('EC-', '') for e in rec.get('EC-NUMBER', [])] for r, rec in reactions.items()}
        # composite reactions (e.g. 2-oxoglutarate dehydrogenase = E1 + E2 + E3 sub-reactions)
        self.composite = {r: set(rec['REACTION-LIST']) for r, rec in reactions.items() if rec.get('REACTION-LIST')}
        self._subset = {}
        self._order = {}

    def is_proper_subset(self, p):
        """True if p's reactions are a proper subset of another base pathway's."""
        if p not in self._subset:
            rs = self.base[p]
            others = None
            for r in rs:
                s = self.rxn2pw.get(r, set()) - {p}
                others = s if others is None else others & s
                if not others:
                    break
            self._subset[p] = bool(others) and any(rs < self.base[q] for q in others)
        return self._subset[p]

    def _edges(self, p):
        """Predecessor edges within p from PREDECESSORS ("(RXN-B RXN-A)" = A precedes B)."""
        if p not in self._order:
            rs = self.base[p]
            preds = {r: set() for r in rs}
            for v in self.pw[p].get('PREDECESSORS', []):
                e = parse_sexpr(v)
                if isinstance(e, list) and e and e[0] in rs:
                    preds[e[0]].update(x for x in e[1:] if isinstance(x, str) and x in rs)
            self._order[p] = preds
        return self._order[p]

    def final_rxns(self, p):
        """Sink reactions plus their direct predecessors ("final two" for a linear pathway)."""
        preds = self._edges(p)
        if not any(preds.values()):
            return set()
        used = set().union(*preds.values())
        sinks = {r for r in preds if r not in used}
        return sinks | set().union(*(preds[s] for s in sinks))

    def initial_rxns(self, p):
        preds = self._edges(p)
        if not any(preds.values()):
            return set()
        return {r for r, ps in preds.items() if not ps}

    def _expand(self, p, seen):
        if p in seen:
            return set()
        seen.add(p)
        out = set()
        for x in self.pw[p].get('REACTION-LIST', []):
            if x in self.pw:
                out |= self._expand(x, seen)
            else:
                out.add(x)
        return out

    def class_ancestors(self, cid, _seen=None):
        """All ancestor classes of a class or TAX frame via TYPES (classes.dat)."""
        _seen = set() if _seen is None else _seen
        stack = [cid]
        while stack:
            c = stack.pop()
            for parent in self.cls.get(c, {}).get('TYPES', []):
                if parent not in _seen:
                    _seen.add(parent)
                    stack.append(parent)
        return _seen

    def pathway_classes(self, p):
        out = set()
        for t in self.types.get(p, []):
            out.add(t)
            out |= self.class_ancestors(t)
        return out


def _clean(s):
    return re.sub(r'<[^>]+>', '', s)
