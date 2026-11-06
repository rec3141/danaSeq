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

W_FIXED = {'EC': 0.45, 'KO': 0.45, 'EC_DIRECT': 0.3}
CLASS_ORDER = ('SEQ', 'EC', 'KO', 'EC_DIRECT')


def w_seq(identity):
    return 0.5 + 0.5 * min(1.0, max(0.0, (identity - 35.0) / 45.0))


class BinEvidence:
    """Best evidence per (bin, MetaCyc reaction), accumulated chunk by chunk:
    rec = [weight, class, best SEQ identity, 'protein:hit:identity' of that hit]."""

    def __init__(self):
        self.bins = {}

    def add(self, b, t):
        r = t['metacyc_rxn']
        if not r:
            return
        rec = self.bins.setdefault(b, {}).get(r)
        if rec is None:
            rec = self.bins[b][r] = [0.0, None, None, '']
        for c, det in t['classes'].items():
            w = w_seq(det['identity']) if c == 'SEQ' else W_FIXED[c]
            # ties go to the stronger class and the smaller hit label, so results do not
            # depend on the order proteins arrive in
            if w > rec[0] or (w == rec[0] and rec[1] and CLASS_ORDER.index(c) < CLASS_ORDER.index(rec[1])):
                rec[0], rec[1] = w, c
            if c == 'SEQ':
                hit = f"{t['protein']}:{det['hit']}:{det['identity']}"
                if rec[2] is None or det['identity'] > rec[2] or (det['identity'] == rec[2] and hit < rec[3]):
                    rec[2], rec[3] = det['identity'], hit


def _weight(mc, r, ev, present):
    if r in ev:
        return ev[r][0], ev[r][1]
    if r in mc.spontaneous:
        return 1.0, 'SPONT'
    parts = mc.composite.get(r)
    if parts and r in present:
        ws = [ev[x][0] for x in parts if x in ev] + [1.0 for x in parts if x in mc.spontaneous]
        best = [ev[x][1] for x in parts if x in ev]
        return (min(ws) if ws else 0.0), (max(best, key=CLASS_ORDER.index) if best else 'COMPOSITE')
    return 0.0, None


def annotate(mc, row, ev, max_hits=12):
    """Add the confidence columns to one pathway row; ev = BinEvidence.bins[bin]."""
    s = row['_sets']
    scored, present = s['scored'], s['present']
    ws, cls = {}, {}
    for r in scored:
        ws[r], cls[r] = _weight(mc, r, ev, present)
    conf = math.fsum(ws.values()) / (len(scored) or 1)     # fsum: independent of set order
    enz = [r for r in scored if cls[r] not in (None, 'SPONT')]
    seq = [r for r in enz if cls[r] == 'SEQ' or (r in ev and ev[r][2] is not None)]
    ids = [ev[r][2] for r in seq if r in ev and ev[r][2] is not None]
    only = {c: sum(1 for r in enz if cls[r] == c) for c in ('EC', 'KO', 'EC_DIRECT')}
    key = s['key']
    key_seq = sum(1 for r in key if r in seq)
    seq_share = len(seq) / len(enz) if enz else 0.0
    if conf >= 0.6 and seq_share >= 0.5 and key_seq == len(key):
        tier = 'high'
    elif conf >= 0.35:
        tier = 'medium'
    else:
        tier = 'low'
    ranked = sorted((r for r in seq if r in ev and ev[r][3]), key=lambda r: (-ev[r][2], r))
    row.update(confidence=round(conf, 3), confidence_tier=tier, n_seq=len(seq),
               n_seq_id60=sum(1 for i in ids if i >= 60), n_ec_only=only['EC'], n_ko_only=only['KO'],
               n_ec_direct_only=only['EC_DIRECT'],
               seq_id_median=round(statistics.median(ids), 1) if ids else '',
               seq_id_min=min(ids) if ids else '',
               key_seq=f'{key_seq}/{len(key)}' if key else '',
               n_unique_seq=sum(1 for r in s['unique_present'] if r in seq),
               supporting_hits=';'.join(f'{r}={ev[r][3]}' for r in ranked[:max_hits]))
    return row
