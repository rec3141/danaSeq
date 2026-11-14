#!/usr/bin/env python3
"""Regression test for the refinement winner streak in mag_analysis/bin/magscot.py.

Bins of the binner that holds the current best bin are accepted in a streak
without rescoring, but every accepted bin must reach the score cutoff
(--threshold, default 0.5), as in upstream MAGScoT refine_bins() (R/refine.R).
Needs pandas.

Toy input, ten single-copy markers, merging skipped to isolate refinement:

    set1  A = a1 (M1-M10)                   score 1.0
    set1  D = d1 (M1-M7)                    score 0.7   control: accepted
    set1  B = b1 + b2 (M1-M6 twice)         score 0.1   below cutoff
    set2  E = e1 (M1-M2)                    completeness 0.2, excluded

After E is excluded, A, D and B are consecutive bins of set1; B must not be
accepted.
"""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'mag_analysis/bin/magscot.py'

C2B = 'A\ta1\tset1\nD\td1\tset1\nB\tb1\tset1\nB\tb2\tset1\nE\te1\tset2\n'
HMM = ''.join(
    [f'a1_{m}\tM{m}\t1e-10\n' for m in range(1, 11)] +
    [f'd1_{m}\tM{m}\t1e-10\n' for m in range(1, 8)] +
    [f'{c}_{m}\tM{m}\t1e-10\n' for m in range(1, 7) for c in ('b1', 'b2')] +
    [f'e1_{m}\tM{m}\t1e-10\n' for m in (1, 2)])


class MagscotRefineStreak(unittest.TestCase):
    def test_streak_stops_at_score_cutoff(self):
        with tempfile.TemporaryDirectory(prefix='magscot-test-') as tmp:
            tmp = Path(tmp)
            (tmp / 'profiles').mkdir()
            (tmp / 'profiles/gtdb_rel207_default_markers.tsv').write_text(
                'marker\tset\n' + ''.join(f'M{i}\tA\n' for i in range(1, 11)))
            (tmp / 'c2b.tsv').write_text(C2B)
            (tmp / 'hmm.tsv').write_text(HMM)
            out = tmp / 'out'
            subprocess.run([sys.executable, str(SCRIPT), '-i', str(tmp / 'c2b.tsv'),
                            '--hmm', str(tmp / 'hmm.tsv'),
                            '--profile_dir', str(tmp / 'profiles'),
                            '--skip_merge_bins', '-o', str(out)],
                           check=True, capture_output=True)
            with open(f'{out}.scores.out') as f:
                scores = {r['bin']: float(r['score'])
                          for r in csv.DictReader(f, delimiter='\t')}
            with open(f'{out}.refined.out') as f:
                refined = [r['bin'] for r in csv.DictReader(f, delimiter='\t')]
        self.assertLess(scores['B'], 0.5)
        self.assertEqual(refined, ['A', 'D'])


if __name__ == '__main__':
    unittest.main()
