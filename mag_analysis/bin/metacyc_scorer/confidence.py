"""Evidence-based confidence for pathway calls.

Calls follow PathoLogic's permissive published rules (pathways.py). These columns carry the
uncertainty. For each scored (non-hypothetical) reaction of a pathway in one bin, the strongest
evidence among the bin's proteins sets a weight w:

  SEQ         0.5 + 0.5 * clip((identity - 35) / 45, 0, 1)    35% -> 0.50, >= 80% -> 1.00
  EC          0.45   (Bakta/eggNOG EC -> Rhea)
  KO          0.45   (KofamScan/eggNOG KO -> Rhea)
  EC_DIRECT   0.30   (EC -> MetaCyc reaction without a Rhea link; broad ECs reach many reactions)
  spontaneous 1.0; composite reaction: weakest of its parts (weight and class); missing 0.

  confidence = mean of w over the pathway's scored reactions (0..1; equals frac when every
               present reaction has >= 80% identity sequence support)
  tier high   : confidence >= 0.6, at least half of the present enzymatic reactions have SEQ
                evidence, and every key reaction (if the pathway names any) is SEQ-supported
  tier medium : confidence >= 0.35
  tier low    : otherwise
The weights and cut-offs are first estimates, not fitted to a reference set.
"""
import math
import statistics

PARAMS = dict(
    w_ec=0.45, w_ko=0.45, w_ec_direct=0.30,     # fixed weights of the non-sequence classes
    seq_id_lo=35.0, seq_id_hi=80.0,             # SEQ weight rises 0.5 -> 1.0 between these identities
    high_conf=0.6, high_seq_share=0.5,          # tier high: confidence and SEQ share of present reactions
    high_key_seq=True,                          # tier high: every key reaction SEQ-supported
    medium_conf=0.35)
CLASS_ORDER = ('SEQ', 'EC', 'KO', 'EC_DIRECT')


def w_seq(identity, p=PARAMS):
    span = p['seq_id_hi'] - p['seq_id_lo']
    return 0.5 + 0.5 * min(1.0, max(0.0, (identity - p['seq_id_lo']) / span))


def class_weight(c, identity, p=PARAMS):
    if c == 'SEQ':
        return w_seq(identity, p)
    return p[{'EC': 'w_ec', 'KO': 'w_ko', 'EC_DIRECT': 'w_ec_direct'}[c]]


class BinEvidence:
    """Evidence per (bin, MetaCyc reaction), accumulated chunk by chunk:
    rec = [best SEQ identity or None, 'protein:hit:identity' of that hit, set of evidence classes].
    Weights are applied in annotate(), so one evidence table can be scored under any PARAMS."""

    def __init__(self):
        self.bins = {}

    def add(self, b, t):
        r = t['metacyc_rxn']
        if not r:
            return
        rec = self.bins.setdefault(b, {}).get(r)
        if rec is None:
            rec = self.bins[b][r] = [None, '', set()]
        for c, det in t['classes'].items():
            rec[2].add(c)
            if c == 'SEQ':
                # ties go to the smaller hit label, so results do not depend on protein order
                hit = f"{t['protein']}:{det['hit']}:{det['identity']}"
                if rec[0] is None or det['identity'] > rec[0] or (det['identity'] == rec[0] and hit < rec[1]):
                    rec[0], rec[1] = det['identity'], hit


def best(rec, p=PARAMS):
    """(weight, class) of the strongest evidence in rec; ties go to the stronger class."""
    out = (0.0, None)
    for c in CLASS_ORDER:
        if c in rec[2]:
            w = class_weight(c, rec[0], p)
            if w > out[0]:
                out = (w, c)
    return out


def _weight(mc, r, ev, present, p):
    if r in ev:
        return best(ev[r], p)
    if r in mc.spontaneous:
        return 1.0, 'SPONT'
    parts = mc.composite.get(r)
    if parts and r in present:
        bs = [best(ev[x], p) for x in parts if x in ev]
        ws = [w for w, _ in bs] + [1.0 for x in parts if x in mc.spontaneous]
        cl = [c for _, c in bs]
        return (min(ws) if ws else 0.0), (max(cl, key=CLASS_ORDER.index) if cl else 'COMPOSITE')
    return 0.0, None


def annotate(mc, row, ev, max_hits=12, params=None):
    """Add the confidence columns to one pathway row; ev = BinEvidence.bins[bin]."""
    p = dict(PARAMS, **(params or {}))
    s = row['_sets']
    scored, present = s['scored'], s['present']
    ws, cls = {}, {}
    for r in scored:
        ws[r], cls[r] = _weight(mc, r, ev, present, p)
    conf = math.fsum(ws.values()) / (len(scored) or 1)     # fsum: independent of set order
    enz = [r for r in scored if cls[r] not in (None, 'SPONT')]
    seq = [r for r in enz if cls[r] == 'SEQ' or (r in ev and ev[r][0] is not None)]
    ids = [ev[r][0] for r in seq if r in ev and ev[r][0] is not None]
    only = {c: sum(1 for r in enz if cls[r] == c) for c in ('EC', 'KO', 'EC_DIRECT')}
    key = s['key']
    key_seq = sum(1 for r in key if r in seq)
    seq_share = len(seq) / len(enz) if enz else 0.0
    if conf >= p['high_conf'] and seq_share >= p['high_seq_share'] and \
            (key_seq == len(key) or not p['high_key_seq']):
        tier = 'high'
    elif conf >= p['medium_conf']:
        tier = 'medium'
    else:
        tier = 'low'
    ranked = sorted((r for r in seq if r in ev and ev[r][1]), key=lambda r: (-ev[r][0], r))
    row.update(confidence=round(conf, 3), confidence_tier=tier, n_seq=len(seq),
               n_seq_id60=sum(1 for i in ids if i >= 60), n_ec_only=only['EC'], n_ko_only=only['KO'],
               n_ec_direct_only=only['EC_DIRECT'],
               seq_id_median=round(statistics.median(ids), 1) if ids else '',
               seq_id_min=min(ids) if ids else '',
               key_seq=f'{key_seq}/{len(key)}' if key else '',
               n_unique_seq=sum(1 for r in s['unique_present'] if r in seq),
               supporting_hits=';'.join(f'{r}={ev[r][1]}' for r in ranked[:max_hits]))
    return row
