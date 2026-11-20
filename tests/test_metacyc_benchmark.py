#!/usr/bin/env python3
"""tools/metacyc_benchmark.py on the hand-made MetaCyc fixture (no licensed content)."""
import contextlib
import gzip
import io
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'tools'))
import metacyc_benchmark as B  # noqa: E402

FX = os.path.join(HERE, 'fixtures', 'metacyc')
COLS = ['bin', 'pathway', 'super', 'call', 'confidence_tier', 'confidence', 'reason']


class TestCompare(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.truth = os.path.join(self.d, 'truth.dat')
        with open(self.truth, 'w') as fh:
            for p in ('PWY-A', 'PWY-T', 'PWY-S', 'PWY-ORG-ONLY'):
                fh.write(f'UNIQUE-ID - {p}\n//\n')
        self.calls = os.path.join(self.d, 'metacyc_pathways.tsv')
        rows = [('b', 'PWY-A', '0', '1', 'high'), ('b', 'PWY-B', '0', '1', 'high'),
                ('b', 'PWY-S', '1', '1', 'medium'), ('b', 'PWY-D', '0', '0', 'low'),
                ('other', 'PWY-T', '0', '1', 'high')]
        with open(self.calls, 'w') as fh:
            fh.write('\t'.join(COLS) + '\n')
            for r in rows:
                fh.write('\t'.join(r + ('0.5', 'x')) + '\n')

    def test_truth_and_metrics(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            B.main(['compare', '--truth', self.truth, '--metacyc', FX, '--pathways', self.calls,
                    '--bin', 'b', '--detail', os.path.join(self.d, 'detail.tsv')])
        self.assertIn('base: truth 2, reachable 1', out.getvalue())
        from metacyc_scorer import ptdat
        mc = ptdat.MetaCyc(ptdat.load(os.path.join(FX, 'pathways.dat')),
                           ptdat.load(os.path.join(FX, 'reactions.dat')))
        ts = B.truth_sets(self.truth, mc, B.linkable_reactions(mc))
        self.assertEqual(ts['not_in_metacyc'], {'PWY-ORG-ONLY'})
        self.assertEqual(ts['no_reaction_links'], {'PWY-T'})            # no Rhea link, no EC
        self.assertEqual((ts['base'], ts['super']), ({'PWY-A', 'PWY-T'}, {'PWY-S'}))
        calls = B.read_calls(self.calls, 'b')
        m = B.metrics(calls, ts, 'base')
        self.assertEqual((m['n_truth'], m['n_reachable']), (2, 1))
        self.assertEqual((m['high']['n'], m['high']['tp'], m['high']['precision'], m['high']['recall']),
                         (2, 1, 0.5, 1.0))
        self.assertEqual(m['low']['cum_n'], 2)
        s = B.metrics(calls, ts, 'super')
        self.assertEqual((s['medium']['n'], s['medium']['tp'], s['high']['n']), (1, 1, 0))
        with open(os.path.join(self.d, 'detail.tsv')) as fh:
            body = fh.read()
        self.assertIn('PWY-T\t0\t1\t0\t0\t\t\tno_reaction_links', body)
        self.assertIn('not_in_metacyc', body)


class TestSpecies(unittest.TestCase):
    def test_species_truth(self):
        d = tempfile.mkdtemp()
        with open(os.path.join(d, 'pathways.dat'), 'w') as fh:
            fh.write('UNIQUE-ID - P1\nSPECIES - TAX-9\nSPECIES - TAX-7\n//\nUNIQUE-ID - P2\nSPECIES - TAX-8\n//\n'
                     'UNIQUE-ID - P3\n//\n')
        out = os.path.join(d, 't.dat')
        B.main(['species', '--metacyc', d, '--taxa', 'TAX-7,TAX-5', '--out', out])
        self.assertEqual(set(B._ids(out)), {'P1'})


class TestDegrade(unittest.TestCase):
    def test_blocks_removed_everywhere(self):
        d = tempfile.mkdtemp()
        tsv, merged, hits = (os.path.join(d, x) for x in ('a.tsv', 'm.tsv', 'h.tsv.gz'))
        with open(tsv, 'w') as fh:
            fh.write('#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\n')
            for i in range(100):
                fh.write(f'c1\tcds\t{i}\t{i + 1}\t+\tL{i:03d}\n')
            fh.write('c1\ttRNA\t1\t2\t+\tT1\n')
        with open(merged, 'w') as fh:
            fh.write('protein_id\tKO\tEC\n' + ''.join(f'L{i:03d}\tK1\t-\n' for i in range(100)))
        with gzip.open(hits, 'wt') as fh:
            fh.write(''.join(f'L{i:03d}\tP1\t90\t1\t1\t1\t0\t1\n' for i in range(100)))
        out = os.path.join(d, 'out')
        B.main(['degrade', '--bakta-tsv', tsv, '--merged', merged, '--hits', hits, '--drop', '0.4',
                '--block', '10', '--seed', '3', '--out', out])
        with open(os.path.join(out, 'annotation.tsv')) as fh:
            kept = [line.split('\t')[5].strip() for line in fh if not line.startswith('#')]
        self.assertEqual(len(kept), 61)                                  # 60 CDS + the tRNA
        with open(os.path.join(out, 'merged_annotations.tsv')) as fh:
            self.assertEqual(len(fh.readlines()) - 1, 60)
        with gzip.open(os.path.join(out, 'metacyc_hits.tsv.gz'), 'rt') as fh:
            self.assertEqual({line.split('\t')[0] for line in fh}, set(kept) - {'T1'})


if __name__ == '__main__':
    unittest.main()
