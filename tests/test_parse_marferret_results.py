#!/usr/bin/env python3
"""Regression test for mag_analysis/bin/parse_marferret_results.py contig IDs.

DIAMOND's qseqid is the MetaEuk protein header, T_acc|C_acc|strand|..., so the
contig_id columns must carry C_acc, never the reference target T_acc. The first
header below is verbatim MetaEuk easy-predict output (MetaEuk commit 1f7d7023)
for a contig named my_contig_42 and a reference protein TARGET_ACC_1.
"""
import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'mag_analysis/bin/parse_marferret_results.py'

HITS = [
    # (qseqid, expected contig_id)
    ('TARGET_ACC_1|my_contig_42|+|440|1.409e-130|1|300|899|300[300]:899[899]:600[600]',
     'my_contig_42'),
    ('TARGET_ACC_1|my_contig_42|-|120|2e-30|2|10|500|500[500]:400[400]:101[101]|200[200]:10[10]:191[191]',
     'my_contig_42'),
    ('REF9|odd|name|+|50|1e-5|1|1|90|1[1]:90[90]:90[90]', 'odd|name'),
    # Control: a plain protein ID (no strand token) keeps its first field.
    ('contigB_7', 'contigB_7'),
]


class MarferretContigIds(unittest.TestCase):
    def test_contig_id_is_metaeuk_contig_not_target(self):
        with tempfile.TemporaryDirectory(prefix='marferret-test-') as tmp:
            tmp = Path(tmp)
            diamond = tmp / 'hits.tsv'
            diamond.write_text(''.join(
                f'{q}\tSUBJ{i}\t90.0\t100\t1e-20\t200\tsubject {i}\n'
                for i, (q, _) in enumerate(HITS)))
            prot, contig = tmp / 'proteins.tsv', tmp / 'contigs.tsv'
            subprocess.run([sys.executable, str(SCRIPT), '--diamond', str(diamond),
                            '--out-proteins', str(prot), '--out-contigs', str(contig)],
                           check=True, capture_output=True)
            with open(prot) as f:
                rows = {r['protein_id']: r['contig_id']
                        for r in csv.DictReader(f, delimiter='\t')}
            self.assertEqual(rows, {q: c for q, c in HITS})
            with open(contig) as f:
                contigs = {r['contig_id']: int(r['n_proteins'])
                           for r in csv.DictReader(f, delimiter='\t')}
            self.assertEqual(contigs, {'my_contig_42': 2, 'odd|name': 1, 'contigB_7': 1})


class MarferretTaxonomy(unittest.TestCase):
    """MarFERReT.v1.1.1.taxonomies.tab.gz is NCBI accession2taxid layout:
    accession is NA on every row and the MarFERReT ID is accession.version.
    Lineages come from an NCBI taxdump, as domain;phylum;...;species."""

    def test_taxid_join_and_lineage(self):
        with tempfile.TemporaryDirectory(prefix='marferret-test-') as tmp:
            tmp = Path(tmp)
            with gzip.open(tmp / 'tax.tab.gz', 'wt') as f:
                f.write('accession\taccession.version\ttaxid\tgi\n'
                        'NA\tmft0000000001\t2850\tNA\n'
                        'NA\tmft0000000002\t9999\tNA\n'   # merged into 2850
                        'NA\tmft0000000003\t555\tNA\n')   # never hit
            dump = tmp / 'taxdump'
            dump.mkdir()
            nodes = [('1', '1', 'no rank'), ('2759', '1', 'superkingdom'),
                     ('2836', '2759', 'phylum'), ('2850', '2836', 'class'),
                     ('555', '2759', 'phylum')]
            (dump / 'nodes.dmp').write_text(''.join(
                f'{t}\t|\t{p}\t|\t{r}\t|\n' for t, p, r in nodes))
            names = [('1', 'root'), ('2759', 'Eukaryota'), ('2836', 'Bacillariophyta'),
                     ('2850', 'Bacillariophyceae'), ('555', 'Other')]
            (dump / 'names.dmp').write_text(''.join(
                f'{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n' for t, n in names))
            (dump / 'merged.dmp').write_text('9999\t|\t2850\t|\n')
            pfam = tmp / 'pfam.csv.gz'
            with gzip.open(pfam, 'wt') as f:
                f.write('aa_id,pfam_name,pfam_id,pfam_eval,pfam_score\n'
                        'mft0000000001,Cu2_monooxygen,PF01082.22,8e-18,75.2\n')
            diamond = tmp / 'hits.tsv'
            diamond.write_text(
                'T1|c1|+|1|1e-5|1|1|90|x\tmft0000000001\t90.0\t100\t1e-20\t200\tmft0000000001\n'
                'T2|c1|-|1|1e-5|1|1|90|x\tmft0000000002\t80.0\t100\t1e-20\t150\tmft0000000002\n'
                'T3|c2|+|1|1e-5|1|1|90|x\tmft0000009999\t70.0\t100\t1e-20\t100\tmft0000009999\n')
            prot, contig = tmp / 'proteins.tsv', tmp / 'contigs.tsv'
            subprocess.run([sys.executable, str(SCRIPT), '--diamond', str(diamond),
                            '--taxonomy', str(tmp / 'tax.tab.gz'), '--pfam', str(pfam),
                            '--taxdump', str(dump),
                            '--out-proteins', str(prot), '--out-contigs', str(contig)],
                           check=True, capture_output=True)
            lineage = 'Eukaryota;Bacillariophyta;Bacillariophyceae;;;;'
            with open(prot) as f:
                rows = {r['best_hit']: r for r in csv.DictReader(f, delimiter='\t')}
            self.assertEqual(rows['mft0000000001']['taxon_id'], '2850')
            self.assertEqual(rows['mft0000000001']['taxonomy'], lineage)
            self.assertEqual(rows['mft0000000001']['pfam'], 'PF01082.22')
            self.assertEqual(rows['mft0000000002']['taxonomy'], lineage)
            self.assertEqual(rows['mft0000009999']['taxonomy'], '')
            with open(contig) as f:
                contigs = {r['contig_id']: r for r in csv.DictReader(f, delimiter='\t')}
            self.assertEqual(contigs['c1']['n_classified'], '2')
            self.assertEqual(contigs['c1']['top_taxonomy'], lineage)
            self.assertEqual(contigs['c2']['n_classified'], '0')


if __name__ == '__main__':
    unittest.main()
