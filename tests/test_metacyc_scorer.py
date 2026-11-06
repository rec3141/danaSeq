#!/usr/bin/env python3
"""MetaCyc scorer (mag_analysis/bin/metacyc_scorer) on a hand-made fixture.

tests/fixtures/metacyc/*.dat are invented records in the Pathway Tools attribute-value
format; they hold no MetaCyc content, which is licensed and must never enter the repo.
"""
import gzip
import os
import subprocess
import sys
import tarfile
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(os.path.dirname(HERE), 'mag_analysis', 'bin')
sys.path.insert(0, BIN)
from metacyc_scorer import confidence, evidence, install, pathways, ptdat, refs, taxonomy  # noqa: E402

FX = os.path.join(HERE, 'fixtures', 'metacyc')


def model():
    return ptdat.MetaCyc(ptdat.load(os.path.join(FX, 'pathways.dat')),
                         ptdat.load(os.path.join(FX, 'reactions.dat')),
                         ptdat.load(os.path.join(FX, 'classes.dat')))


class TestParser(unittest.TestCase):
    def test_records(self):
        pw = ptdat.load(os.path.join(FX, 'pathways.dat'))
        self.assertEqual(sorted(pw), ['PWY-A', 'PWY-B', 'PWY-C', 'PWY-D', 'PWY-E', 'PWY-K', 'PWY-S', 'PWY-T'])
        self.assertEqual(pw['PWY-A']['REACTION-LIST'], ['RXN-A1', 'RXN-A2', 'RXN-A3', 'RXN-SHARED'])
        self.assertEqual(pw['PWY-A']['COMMENT'], ['first line\nsecond line'])

    def test_annotation(self):
        rx = ptdat.load(os.path.join(FX, 'reactions.dat'))
        self.assertEqual(rx['RXN-A1']['_annot'][('EC-NUMBER', 0)], {'OFFICIAL?': ['T']})

    def test_sexpr(self):
        self.assertEqual(ptdat.parse_sexpr('("RXN-A2" "RXN-A1")'), ['RXN-A2', 'RXN-A1'])
        self.assertEqual(ptdat.parse_sexpr('(R (:L X) (:D :L2R))')[1], [':L', 'X'])

    def test_model(self):
        mc = model()
        self.assertEqual(mc.base['PWY-S'], {'RXN-A1', 'RXN-A2', 'RXN-A3', 'RXN-SHARED', 'RXN-B1', 'RXN-SPONT'})
        self.assertEqual(mc.rxn2pw['RXN-SHARED'], {'PWY-A', 'PWY-B', 'PWY-C', 'PWY-K'})
        self.assertEqual(mc.composite['RXN-COMP'], {'RXN-P1', 'RXN-P2', 'RXN-P3'})
        self.assertNotIn('PWY-S', mc.rxn2pw['RXN-A1'])          # super-pathways excluded
        self.assertEqual(mc.rhea2rxn['200'], {'RXN-SHARED'})
        self.assertEqual(mc.spontaneous, {'RXN-SPONT'})
        self.assertEqual(mc.name['PWY-A'], 'fixture A biosynthesis')
        self.assertIn('TAX-131567', mc.class_ancestors('TAX-1224'))


class TestScoring(unittest.TestCase):
    BACT = {'TAX-1224', 'TAX-2', 'TAX-131567', 'TAX-1'}

    def setUp(self):
        self.mc = model()

    def call(self, present, pwy, lineage=None):
        r = pathways.score_pathways(self.mc, {x: 1 for x in present}, lineage)[pwy]
        return r['call'], r['reason']

    def test_key_reactions(self):
        self.assertEqual(self.call(['RXN-A2', 'RXN-A3', 'RXN-SHARED'], 'PWY-A'), (0, 'key_reaction_missing'))
        # all key reactions present -> inferred even though 3 of 4 reactions are missing
        self.assertEqual(self.call(['RXN-A1'], 'PWY-A'), (1, 'all_key_reactions_present'))

    def test_spontaneous_counts_but_is_not_evidence(self):
        rows = pathways.score_pathways(self.mc, {})
        self.assertEqual(rows['PWY-B']['frac'], round(1 / 3, 4))
        self.assertEqual(rows['PWY-B']['reason'], 'no_enzymatic_evidence')

    def test_mostly_absent_and_unique(self):
        self.assertEqual(self.call(['RXN-C1'], 'PWY-C'), (0, 'mostly_absent'))   # 1 present, 2 missing
        rows = pathways.score_pathways(self.mc, {'RXN-C1': 1, 'RXN-C2': 1})
        self.assertEqual(rows['PWY-C']['n_rxn'], 3)                                # hypothetical C3 excluded
        self.assertEqual((rows['PWY-C']['call'], rows['PWY-C']['reason']), (1, 'unique_reaction_present'))

    def test_taxonomic_range(self):
        self.assertEqual(self.call(['RXN-B1'], 'PWY-B')[0], 1)                      # no lineage given
        self.assertEqual(self.call(['RXN-B1'], 'PWY-B', self.BACT), (0, 'outside_taxonomic_range'))
        self.assertEqual(self.call(['RXN-B1', 'RXN-SHARED'], 'PWY-B', self.BACT)[0], 1)   # complete
        # lineage resolved only to cellular organisms: a Bacteria range cannot be excluded
        rows = pathways.score_pathways(self.mc, {'RXN-A1': 1, 'RXN-A2': 1}, {'TAX-131567', 'TAX-1'})
        self.assertEqual((rows['PWY-A']['in_tax_range'], rows['PWY-A']['call']), ('U', 1))

    def test_subset_missing_final(self):
        self.assertTrue(self.mc.is_proper_subset('PWY-D'))
        self.assertFalse(self.mc.is_proper_subset('PWY-E'))
        self.assertEqual(self.mc.final_rxns('PWY-D'), {'RXN-D3', 'RXN-D2'})
        self.assertEqual(self.mc.initial_rxns('PWY-E'), {'RXN-D0'})
        self.assertEqual(self.call(['RXN-D0', 'RXN-D1'], 'PWY-D'), (0, 'biosynthetic_subset_missing_final'))
        self.assertEqual(self.call(['RXN-D0', 'RXN-D1'], 'PWY-E'), (1, 'no_rejection_rule'))

    def test_superset_more_missing(self):
        rows = pathways.score_pathways(self.mc, {x: 1 for x in ['RXN-D0', 'RXN-D1', 'RXN-D2', 'RXN-D3']})
        self.assertEqual(rows['PWY-D']['call'], 1)
        self.assertEqual((rows['PWY-E']['call'], rows['PWY-E']['reason']), (0, 'superset_of_PWY-D_more_missing'))

    def test_electron_transfer_needs_all(self):
        self.assertEqual(self.call(['RXN-T1'], 'PWY-T'), (0, 'electron_transfer_incomplete'))
        self.assertEqual(self.call(['RXN-T1', 'RXN-T2'], 'PWY-T')[0], 1)

    def test_energy_and_composite(self):
        self.assertEqual(self.call(['RXN-B1', 'RXN-SHARED'], 'PWY-K'), (0, 'energy_no_unique_lt_half'))
        present, enz = pathways.infer_present(self.mc, {'RXN-P1', 'RXN-P2'})
        self.assertIn('RXN-COMP', enz)                                             # 2 of 3 parts
        self.assertNotIn('RXN-COMP', pathways.infer_present(self.mc, {'RXN-P1'})[1])
        self.assertEqual(self.call(['RXN-P1', 'RXN-P2', 'RXN-K1'], 'PWY-K'), (1, 'unique_reaction_present'))

    def test_super_pathway(self):
        self.assertEqual(self.call(['RXN-A1', 'RXN-B1', 'RXN-SHARED'], 'PWY-S'), (1, 'all_components_inferred'))
        self.assertEqual(self.call(['RXN-A1', 'RXN-B1'], 'PWY-S', self.BACT), (0, 'components_inferred_1_of_2'))

    def test_frac_mode(self):
        rows = pathways.score_pathways(self.mc, {'RXN-A1': 1}, None, {'mode': 'frac'})
        self.assertEqual(rows['PWY-A']['reason'], 'frac<0.50')

    def test_minpath(self):
        sel, _ = pathways.minpath(self.mc, {'RXN-A1': 1, 'RXN-SHARED': 1})
        self.assertEqual(sel, {'PWY-A'})
        sel, _ = pathways.minpath(self.mc, {'RXN-SHARED': 1})
        self.assertEqual(len(sel), 1)
        sel, _ = pathways.minpath(self.mc, {'RXN-A1': 1, 'RXN-B1': 1, 'RXN-C1': 1})
        self.assertEqual(len(sel), 3)
        self.assertTrue({'PWY-A', 'PWY-C'} <= sel)
        self.assertEqual(pathways._greedy_cover(['PWY-A', 'PWY-B', 'PWY-C'],
                                              {'PWY-A': {'x', 'y'}, 'PWY-B': {'y'}, 'PWY-C': {'z'}},
                                              {'x', 'y', 'z'}), {'PWY-A', 'PWY-C'})


class TestEvidence(unittest.TestCase):
    def test_seq_calls(self):
        seqref = {'P1': (True, {'100': (9, 10)}), 'Q2': (False, {'100': (5, 5), '200': (1, 5)}),
                  'Q3': (False, {'300': (1, 1)})}
        rows = [['prot1', 'P1', '82.0', '300', '310', '305', '1e-90', '500'],
                ['prot1', 'Q2', '70.0', '300', '310', '300', '1e-80', '450'],
                ['prot1', 'Q3', '90.0', '100', '310', '110', '1e-30', '200'],     # query coverage 0.32: dropped
                ['prot2', 'Q3', '30.0', '300', '310', '300', '1e-20', '100']]     # identity below 35: dropped
        calls = evidence.seq_calls(rows, seqref, dict(evidence.DEFAULTS))
        self.assertEqual(set(calls), {'prot1'})
        self.assertEqual(set(calls['prot1']), {'100'})              # '200' holds 450/950 < 0.5
        c = calls['prot1']['100']
        self.assertEqual((c['hit'], c['identity'], c['reviewed'], c['cluster_agreement']), ('P1', 82.0, 1, 0.9))


def ev_for(mc, assignments):
    be = confidence.BinEvidence()
    for t in assignments:
        be.add('b', t)
    return be.bins['b']


class TestConfidence(unittest.TestCase):
    def test_weights_and_tiers(self):
        mc = model()
        self.assertAlmostEqual(confidence.w_seq(35), 0.5)
        self.assertAlmostEqual(confidence.w_seq(90), 1.0)

        def seq(r, i):
            return dict(protein='p', metacyc_rxn=r, classes={'SEQ': dict(identity=i, hit='H')})

        def ec(r):
            return dict(protein='q', metacyc_rxn=r, classes={'EC': dict(ec=['1.1.1.1'])})
        ev = ev_for(mc, [seq('RXN-A1', 90), seq('RXN-A2', 80), seq('RXN-A3', 85), ec('RXN-SHARED')])
        rows = pathways.score_pathways(mc, {r: 1 for r in ev})
        r = confidence.annotate(mc, rows['PWY-A'], ev)
        self.assertEqual((r['n_seq'], r['n_ec_only'], r['key_seq']), (3, 1, '1/1'))
        self.assertAlmostEqual(r['confidence'], (3 + 0.45) / 4, places=3)
        self.assertEqual(r['confidence_tier'], 'high')
        ev = ev_for(mc, [ec('RXN-A1'), ec('RXN-A2')])
        rows = pathways.score_pathways(mc, {r: 1 for r in ev})
        r = confidence.annotate(mc, rows['PWY-A'], ev)
        self.assertEqual((r['call'], r['confidence_tier'], r['key_seq']), (1, 'low', '0/1'))

    def test_order_independent(self):
        a = dict(protein='p1', metacyc_rxn='R', classes={'KO': dict(ko=['K1'])})
        b = dict(protein='p2', metacyc_rxn='R', classes={'EC': dict(ec=['1.1.1.1'])})
        self.assertEqual(ev_for(None, [a, b]), ev_for(None, [b, a]))


class TestTaxonomy(unittest.TestCase):
    def test_lineage_lookup(self):
        table = {'p__Actinomycetota': {'TAX-2', 'TAX-201174'}}
        self.assertEqual(taxonomy.lowest_taxon('d__Bacteria;p__Actinomycetota;c__;o__'), 'p__Actinomycetota')
        self.assertEqual(taxonomy.ncbi_lineage(table, 'd__Bacteria;p__Actinomycetota;c__Foo')[0], 'p__Actinomycetota')
        c2b = {'c1': 'b1', 'c2': 'b1', 'c3': 'b1'}
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 's.tsv')
            with open(p, 'w') as fh:
                fh.write('contig_id\tstatus\tANI\tref\tlineage\n'
                         'c1\tC\t95\tx\td__Bacteria;p__A;c__X\nc2\tC\t95\tx\td__Bacteria;p__A;c__Y\n'
                         'c3\tC\t95\tx\td__Bacteria;p__B;c__Z\n')
            out = taxonomy.read_sendsketch_bins(p, c2b, {'c1': 60, 'c2': 30, 'c3': 10})
            g = os.path.join(d, 'gtdbtk.tsv')
            with open(g, 'w') as fh:
                fh.write('user_genome\tclassification\nbin1.fa\td__Bacteria;p__A\nbin2\tUnclassified Bacteria\n')
            self.assertEqual(taxonomy.read_gtdbtk(g), {'bin1': 'd__Bacteria;p__A'})
        self.assertEqual(out, {'b1': 'd__Bacteria;p__A;c__X'})     # c__X holds 60% of bin length

    def test_build_caps_at_gtdb_rank(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'nodes.dmp'), 'w') as fh:
                for tid, par, rank in (('1', '1', 'no rank'), ('2', '1', 'superkingdom'), ('10', '2', 'phylum'),
                                       ('100', '10', 'class'), ('101', '10', 'class'), ('1000', '100', 'species')):
                    fh.write(f'{tid}\t|\t{par}\t|\t{rank}\t|\t\t|\n')   # nodes.dmp has more columns after rank
            meta = os.path.join(d, 'meta.tsv.gz')
            with gzip.open(meta, 'wt') as fh:
                fh.write('accession\tgtdb_taxonomy\tncbi_taxid\n')
                for i in range(3):
                    fh.write(f'g{i}\td__Bacteria;p__P;c__C\t1000\n')
                fh.write('g9\td__Bacteria;p__P;c__D\t101\n')
            out = os.path.join(d, 't.tsv.gz')
            taxonomy.build([meta], d, out, log=open(os.devnull, 'w'))
            t = taxonomy.load(out)
        self.assertEqual(t['p__P'], {'TAX-1', 'TAX-2', 'TAX-10'})       # class 100 has 3/4 but is below phylum
        self.assertEqual(t['c__C'], {'TAX-1', 'TAX-2', 'TAX-10', 'TAX-100'})


def write_inputs(d):
    """A two-contig, four-protein toy run with its own tiny open reference set."""
    rd = os.path.join(d, 'refs', 'rhea')
    os.makedirs(rd)
    os.makedirs(os.path.join(d, 'refs', 'seqref'))
    hdr = 'RHEA_ID\tDIRECTION\tMASTER_ID\tID\n'
    with open(os.path.join(rd, 'rhea2metacyc.tsv'), 'w') as fh:
        fh.write(hdr + '101\tLR\t100\tRXN-A1\n301\tLR\t300\tRXN-A2\n401\tLR\t400\tRXN-B1\n')
    with open(os.path.join(rd, 'rhea2ec.tsv'), 'w') as fh:
        fh.write(hdr + '300\tUN\t300\t1.1.1.2\n')
    with open(os.path.join(rd, 'rhea2kegg_reaction.tsv'), 'w') as fh:
        fh.write(hdr + '400\tUN\t400\tR00001\n')
    with open(os.path.join(rd, 'rhea-directions.tsv'), 'w') as fh:
        fh.write('RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n100\t101\t102\t103\n200\t201\t202\t203\n')
    os.makedirs(os.path.join(d, 'refs', 'kegg'))
    with open(os.path.join(d, 'refs', 'kegg', 'kegg_ko2rn.tsv'), 'w') as fh:
        fh.write('ko:K00001\trn:R00001\n')
    with gzip.open(os.path.join(d, 'refs', 'seqref', 'seqref_rhea.tsv.gz'), 'wt') as fh:
        fh.write('accession\tuniref90\treviewed\trhea_master\tcluster_n_rhea_members\tcluster_n_with_rhea\n'
                 'P1\tUniRef90_P1\t1\t100\t4\t4\nP2\tUniRef90_P2\t0\t200\t2\t1\n')
    with open(os.path.join(d, 'bakta.tsv'), 'w') as fh:
        fh.write('#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n'
                 'c1\tcds\t1\t9\t+\tL1\t\tx\t\n'
                 'c1\tcds\t1\t9\t+\tL2\t\tx\tEC:1.1.1.2\n'
                 'c2\tcds\t1\t9\t+\tL3\t\tx\tKEGG:K00001\n'
                 'c2\tcds\t1\t9\t+\tL4\t\tx\t\n')
    with open(os.path.join(d, 'merged.tsv'), 'w') as fh:
        fh.write('protein_id\tcontig_id\tKO\tEC\nL4\tc2\tko:K00001\t-\n')
    with open(os.path.join(d, 'hits.tsv'), 'w') as fh:
        fh.write('L1\tP1\t88.0\t300\t310\t305\t1e-90\t500\nL4\tP2\t50.0\t300\t310\t305\t1e-50\t300\n')
    with open(os.path.join(d, 'c2b.tsv'), 'w') as fh:
        fh.write('c1\tbinA\nc2\tbinB\n')


class TestRun(unittest.TestCase):
    def run_cli(self, d, out, chunk):
        cmd = [sys.executable, os.path.join(BIN, 'metacyc_scorer.py'), 'run',
               '--bakta-tsv', os.path.join(d, 'bakta.tsv'), '--merged', os.path.join(d, 'merged.tsv'),
               '--hits', os.path.join(d, 'hits.tsv'), '--contig2bin', os.path.join(d, 'c2b.tsv'),
               '--metacyc', FX, '--refs', os.path.join(d, 'refs'), '--chunk', str(chunk), '--outdir', out]
        return subprocess.run(cmd, capture_output=True, text=True)

    def test_end_to_end_and_chunking(self):
        with tempfile.TemporaryDirectory() as d:
            write_inputs(d)
            r1 = self.run_cli(d, os.path.join(d, 'o1'), 1000)
            r2 = self.run_cli(d, os.path.join(d, 'o2'), 1)
            self.assertEqual((r1.returncode, r2.returncode), (0, 0), r1.stderr + r2.stderr)
            with gzip.open(os.path.join(d, 'o1', 'metacyc_reactions.tsv.gz'), 'rt') as fh:
                rows = [line.rstrip('\n').split('\t') for line in fh]
            with gzip.open(os.path.join(d, 'o2', 'metacyc_reactions.tsv.gz'), 'rt') as fh:
                self.assertEqual(sorted(rows), sorted(line.rstrip('\n').split('\t') for line in fh))
            got = {(r[0], r[3]): r[5] for r in rows[1:]}
            self.assertEqual(got[('L1', 'RXN-A1')], 'SEQ')
            self.assertEqual(got[('L2', 'RXN-A2')], 'EC')
            self.assertEqual(got[('L3', 'RXN-B1')], 'KO')
            self.assertNotIn(('L4', 'RXN-A1'), got)                       # P2's reaction has no MetaCyc link
            with open(os.path.join(d, 'o1', 'metacyc_pathways.tsv')) as fh:
                p1 = fh.read()
            with open(os.path.join(d, 'o2', 'metacyc_pathways.tsv')) as fh:
                self.assertEqual(p1, fh.read())
            calls = {(x[0], x[1]): x for x in (line.split('\t') for line in p1.splitlines()[1:])}
            self.assertIn(('binA', 'PWY-A'), calls)
            self.assertIn(('_community', 'PWY-A'), calls)

    def test_missing_metacyc_exit_3(self):
        with tempfile.TemporaryDirectory() as d:
            write_inputs(d)
            cmd = [sys.executable, os.path.join(BIN, 'metacyc_scorer.py'), 'run',
                   '--bakta-tsv', os.path.join(d, 'bakta.tsv'), '--hits', os.path.join(d, 'hits.tsv'),
                   '--contig2bin', os.path.join(d, 'c2b.tsv'), '--metacyc', os.path.join(d, 'none'),
                   '--refs', os.path.join(d, 'refs'), '--outdir', os.path.join(d, 'o')]
            r = subprocess.run(cmd, capture_output=True, text=True)
            self.assertEqual(r.returncode, 3)
            self.assertIn('download-databases.sh --metacyc', r.stderr)


class TestRefs(unittest.TestCase):
    def test_select_reference(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, 'rhea2uniprot_sprot.tsv'), 'w') as fh:
                fh.write('RHEA_ID\tDIRECTION\tMASTER_ID\tID\n11\tLR\t10\tS1\n')
            with gzip.open(os.path.join(d, 'rhea2uniprot_trembl.tsv.gz'), 'wt') as fh:
                fh.write('RHEA_ID\tDIRECTION\tMASTER_ID\tID\n10\tUN\t10\tT1\n20\tUN\t20\tT2\n'
                         '20\tUN\t20\tT3\n30\tUN\t30\tT4\n')
            idm = os.path.join(d, 'idm.tab.gz')
            cols = lambda a, ur: '\t'.join([a] + [''] * 7 + [ur] + [''] * 13)   # col 9 = UniRef90
            with gzip.open(idm, 'wt') as fh:
                fh.write('\n'.join([cols('S1', 'UniRef90_S1'), cols('T1', 'UniRef90_S1'),
                                     cols('T2', 'UniRef90_T3'), cols('T3', 'UniRef90_T3'),
                                     cols('T4', 'UniRef90_X'), cols('N1', 'UniRef90_X')]) + '\n')
            chosen = refs.select_reference(d, idm, d, d)
            with gzip.open(os.path.join(d, 'seqref_rhea.tsv.gz'), 'rt') as fh:
                rows = [line.rstrip('\n').split('\t') for line in fh][1:]
        # Swiss-Prot member stands for its cluster; else the representative (T3), else the first
        self.assertEqual(chosen, {'S1', 'T3', 'T4'})
        self.assertIn(['S1', 'UniRef90_S1', '1', '10', '2', '2'], rows)
        self.assertIn(['T3', 'UniRef90_T3', '0', '20', '2', '2'], rows)


class TestInstall(unittest.TestCase):
    def test_flatfiles_tarball(self):
        with tempfile.TemporaryDirectory() as d:
            tb = os.path.join(d, 'mc.tar.gz')
            with tarfile.open(tb, 'w:gz') as tf:
                for f in ('pathways.dat', 'reactions.dat', 'classes.dat'):
                    tf.add(os.path.join(FX, f), arcname=f'metacyc/test/data/{f}')
            install.main(['--db-dir', os.path.join(d, 'db'), '--flatfiles', tb])
            out = os.path.join(d, 'db', 'metacyc')
            self.assertTrue(os.path.exists(os.path.join(out, 'pathways.dat')))
            self.assertEqual(os.stat(out).st_mode & 0o777, 0o700)
            with open(os.path.join(out, 'VERSION.json')) as fh:
                self.assertIn('"release": "test"', fh.read())


if __name__ == '__main__':
    unittest.main()
