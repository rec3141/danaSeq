#!/usr/bin/env python3
"""Write the fixed sample subsets for the sample-count curve (ablate-kcurve.sh).

    make-kcurve-manifest.py sample_subject_mapping.tsv > kcurve_manifest.tsv

Two families at each k in 2, 4, 6, 8, 10, five draws each, seeded from k and
the draw number so the manifest is reproducible:
  spread  k samples from k different subjects (one sample each, picked at random)
  paired  both samples of k/2 subjects
Plus one 'all' row with all 20 samples. Comparing the families at equal k asks
whether more subjects or more time points per subject add more.
"""
import csv, random, sys
from collections import defaultdict

subj = defaultdict(list)
for r in csv.DictReader(open(sys.argv[1]), delimiter='\t'):
    subj[r['subject_id']].append(r['sample_id'])
subjects = sorted(subj, key=lambda s: int(s[1:]))
for s in subjects:
    subj[s].sort(key=lambda x: int(x.split('_')[1]))
    assert len(subj[s]) == 2, s
print('config\tfamily\tk\tdraw\tseed\tsamples')
for k in (2, 4, 6, 8, 10):
    for draw in range(5):
        seed = 1000 * k + draw
        rng = random.Random(seed)
        picked = rng.sample(subjects, k)
        spread = [rng.choice(subj[s]) for s in picked]
        print(f'spread_k{k}_d{draw}\tspread\t{k}\t{draw}\t{seed}\t{",".join(spread)}')
        rng = random.Random(seed + 500)
        pair_subj = rng.sample(subjects, k // 2)
        paired = [x for s in pair_subj for x in subj[s]]
        print(f'paired_k{k}_d{draw}\tpaired\t{k}\t{draw}\t{seed + 500}\t{",".join(paired)}')
all_s = [x for s in subjects for x in subj[s]]
print(f'all_k20\tall\t20\t0\t0\t{",".join(all_s)}')
