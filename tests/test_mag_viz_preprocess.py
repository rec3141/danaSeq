#!/usr/bin/env python3
"""Regression and control tests for mag_analysis taxonomy and viz preprocessing.

Covers:
- F2: Kraken2 fixed-rank lineage formatting and rank-prefix alignment in parse_lineage
- F3: DefenseFinder contig resolution for multi-underscore and zero-underscore contigs
- F4: geNomad provirus '|provirus_<start>_<end>' host-contig resolution across joins
- F5: antiSMASH region coordinate overlap filtering for CDS gene_kind counts
- F6: Prokka 7-column TSV + sibling GFF loading and reconciled BAKTA_EXTRA selection
"""

from __future__ import annotations

import gzip
import importlib.util
import io
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
MAG_DIR = REPO_ROOT / "mag_analysis"
BIN_DIR = MAG_DIR / "bin"
MODULES_DIR = MAG_DIR / "modules"
PREPROCESS_DIR = MAG_DIR / "viz" / "preprocess"


def _load_module(module_name: str, file_path: Path):
    parent = str(file_path.parent)
    if parent not in sys.path:
        sys.path.insert(0, parent)
    spec = importlib.util.spec_from_file_location(module_name, str(file_path))
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PREPROCESS_MOD = _load_module("danaseq_preprocess_test", PREPROCESS_DIR / "preprocess.py")
GENES_MOD = _load_module("danaseq_genes_to_json_test", PREPROCESS_DIR / "genes_to_json.py")
PARALLEL_DF_MOD = _load_module("danaseq_parallel_df_test", BIN_DIR / "parallel_defensefinder.py")


class MagVizPreprocessTests(unittest.TestCase):
    """Desired-behavior and control tests for code-03 F2, F3, F4, F5, and F6."""

    def test_f2_kraken2_sparse_and_dual_domain_lineages_align_ranks(self):
        """Sparse viral and dual superkingdom/domain bacterial lineages map to exact ranks."""
        tax_nf = (MODULES_DIR / "taxonomy.nf").read_text(encoding="utf-8")
        m = re.search(
            r'python3 - kraken2_raw\.tsv "\\\$taxfmt" "\\\$nodes" "\\\$names" <<\'PYEOF\'\n(.*?)\nPYEOF',
            tax_nf,
            flags=re.DOTALL,
        )
        self.assertIsNotNone(m, "Kraken2 Python post-processor block missing in taxonomy.nf")
        py_script = m.group(1).replace("\\\\", "\\")

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            nodes_dmp = tmp_path / "nodes.dmp"
            nodes_dmp.write_text(
                "1\t|\t1\t|\tno rank\t|\n"
                "10\t|\t1\t|\tsuperkingdom\t|\n"
                "11\t|\t10\t|\tclade\t|\n"
                "12\t|\t11\t|\tfamily\t|\n"
                "13\t|\t12\t|\tgenus\t|\n"
                "20\t|\t1\t|\tsuperkingdom\t|\n"
                "21\t|\t20\t|\tdomain\t|\n"
                "22\t|\t21\t|\tphylum\t|\n"
                "23\t|\t22\t|\tclass\t|\n"
                "24\t|\t23\t|\torder\t|\n"
                "25\t|\t24\t|\tfamily\t|\n"
                "26\t|\t25\t|\tgenus\t|\n",
                encoding="utf-8",
            )
            names_dmp = tmp_path / "names.dmp"
            names_dmp.write_text(
                "10\t|\tViruses\t|\t\t|\tscientific name\t|\n"
                "11\t|\tCaudoviricetes\t|\t\t|\tscientific name\t|\n"
                "12\t|\tSiphoviridae\t|\t\t|\tscientific name\t|\n"
                "13\t|\tLambdavirus\t|\t\t|\tscientific name\t|\n"
                "20\t|\tcellular organisms\t|\t\t|\tscientific name\t|\n"
                "21\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
                "22\t|\tPseudomonadota\t|\t\t|\tscientific name\t|\n"
                "23\t|\tGammaproteobacteria\t|\t\t|\tscientific name\t|\n"
                "24\t|\tEnterobacterales\t|\t\t|\tscientific name\t|\n"
                "25\t|\tEnterobacteriaceae\t|\t\t|\tscientific name\t|\n"
                "26\t|\tEscherichia\t|\t\t|\tscientific name\t|\n",
                encoding="utf-8",
            )
            raw_tsv = tmp_path / "kraken2_raw.tsv"
            raw_tsv.write_text(
                "C\tctg_virus\t13\t5000\t13:40\n"
                "C\tctg_dual\t26\t8000\t26:80\n",
                encoding="utf-8",
            )

            k2_dir = tmp_path / "results" / "taxonomy" / "kraken2"
            k2_dir.mkdir(parents=True)
            subprocess.run(
                [sys.executable, "-B", "-", str(raw_tsv), "dmp", str(nodes_dmp), str(names_dmp)],
                input=py_script,
                cwd=str(k2_dir),
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
                capture_output=True,
                text=True,
                check=True,
            )

            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                maps = PREPROCESS_MOD.build_taxonomy_maps(str(tmp_path / "results"), kaiju_df=None)

            self.assertEqual(
                maps["kraken2"]["ctg_virus"],
                {
                    "domain": "Viruses",
                    "phylum": "",
                    "class": "",
                    "order": "",
                    "family": "Siphoviridae",
                    "genus": "Lambdavirus",
                },
            )
            self.assertEqual(
                maps["kraken2"]["ctg_dual"],
                {
                    "domain": "Bacteria",
                    "phylum": "Pseudomonadota",
                    "class": "Gammaproteobacteria",
                    "order": "Enterobacterales",
                    "family": "Enterobacteriaceae",
                    "genus": "Escherichia",
                },
            )

        # Direct parse_lineage check on sparse rank-prefixed strings without empty placeholders
        self.assertEqual(
            PREPROCESS_MOD.parse_lineage("d__Viruses;f__Siphoviridae;g__Lambdavirus"),
            {
                "domain": "Viruses",
                "phylum": "",
                "class": "",
                "order": "",
                "family": "Siphoviridae",
                "genus": "Lambdavirus",
            },
        )

    def test_f2_control_dense_gtdb_and_unprefixed_kaiju_lineages_unchanged(self):
        """Complete 6-rank GTDB Kraken2 and unprefixed Kaiju strings map identically."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            k2_dir = tmp_path / "taxonomy" / "kraken2"
            k2_dir.mkdir(parents=True)
            (k2_dir / "kraken2_contigs.tsv").write_text(
                "contig_id\tstatus\ttaxid\ttaxon_name\tlineage\n"
                "ctg_dense\tC\t1578\tLactobacillus\t"
                "d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus\n",
                encoding="utf-8",
            )
            kaiju_df = pd.DataFrame(
                [
                    {
                        "contig_id": "ctg_kaiju",
                        "lineage": "Bacteria;Pseudomonadota;Gammaproteobacteria;"
                        "Enterobacterales;Enterobacteriaceae;Escherichia",
                    }
                ]
            )
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                maps = PREPROCESS_MOD.build_taxonomy_maps(str(tmp_path), kaiju_df=kaiju_df)

        self.assertEqual(
            maps["kraken2"]["ctg_dense"],
            {
                "domain": "Bacteria",
                "phylum": "Firmicutes",
                "class": "Bacilli",
                "order": "Lactobacillales",
                "family": "Lactobacillaceae",
                "genus": "Lactobacillus",
            },
        )
        self.assertEqual(
            maps["kaiju"]["ctg_kaiju"],
            {
                "domain": "Bacteria",
                "phylum": "Pseudomonadota",
                "class": "Gammaproteobacteria",
                "order": "Enterobacterales",
                "family": "Enterobacteriaceae",
                "genus": "Escherichia",
            },
        )

    def test_f3_defensefinder_multi_and_zero_underscore_contigs_join_bins(self):
        """Multi-underscore SPAdes contigs and zero-underscore contigs resolve accurately."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            spades_ctg = "NODE_1_length_5000_cov_12.3"
            plain_ctg = "ctg0001"
            unbinned_ctg = "NODE_9_length_3000_cov_4.1"

            gff = tmp_path / "annot.gff"
            gff.write_text(
                f"{spades_ctg}\tBakta\tCDS\t10\t900\t.\t+\t0\tID=BAKTA_00001\n"
                f"{plain_ctg}\tBakta\tCDS\t20\t800\t.\t+\t0\tID=BAKTA_00002\n"
                f"{unbinned_ctg}\tBakta\tCDS\t30\t700\t.\t+\t0\tID=BAKTA_00003\n",
                encoding="utf-8",
            )
            faa = tmp_path / "annot.faa"
            faa.write_text(
                ">BAKTA_00001\nMKTAYIAKQRQISFVK\n"
                ">BAKTA_00002\nMKTAYIAKQRQISFVK\n"
                ">BAKTA_00003\nMKTAYIAKQRQISFVK\n",
                encoding="utf-8",
            )

            prot2contig = PARALLEL_DF_MOD.parse_gff_protein_map(str(gff))
            _, rev_map = PARALLEL_DF_MOD.rewrite_and_group_fasta(str(faa), prot2contig)

            chunk_dir = tmp_path / "chunk_0_out"
            chunk_dir.mkdir()
            (chunk_dir / "chunk_0_defense_finder_systems.tsv").write_text(
                PARALLEL_DF_MOD.SYSTEMS_HEADER
                + f"{spades_ctg}_RM_Type_I_1\tRM\tRM_Type_I\tDefense\t"
                f"{spades_ctg}_BAKTA00001\t{spades_ctg}_BAKTA00001\t{spades_ctg}_BAKTA00001\t1\tRM_M\n"
                + f"{plain_ctg}_AbiD_1\tAbiD\tAbiD\tDefense\t"
                f"{plain_ctg}_BAKTA00002\t{plain_ctg}_BAKTA00002\t{plain_ctg}_BAKTA00002\t1\tAbiD\n"
                + f"{unbinned_ctg}_Cas_1\tCas\tCas\tDefense\t"
                f"{unbinned_ctg}_BAKTA00003\t{unbinned_ctg}_BAKTA00003\t{unbinned_ctg}_BAKTA00003\t1\tCas1\n",
                encoding="utf-8",
            )
            (chunk_dir / "chunk_0_defense_finder_genes.tsv").write_text(
                PARALLEL_DF_MOD.GENES_HEADER, encoding="utf-8"
            )
            (chunk_dir / "chunk_0_defense_finder_hmmer.tsv").write_text(
                PARALLEL_DF_MOD.HMMER_HEADER, encoding="utf-8"
            )

            merged_dir = tmp_path / "merged"
            PARALLEL_DF_MOD.merge_all([str(chunk_dir)], str(merged_dir), rev_map)
            defense_df = PREPROCESS_MOD.load_tsv(str(merged_dir / "systems.tsv"))

            contig2bin = pd.DataFrame(
                {"contig": [spades_ctg, plain_ctg], "bin": ["dastool-semibin_001", "dastool-semibin_001"]}
            )
            dastool_summary = pd.DataFrame(
                [
                    {
                        "bin": "dastool-semibin_001",
                        "SCG_completeness": 95.0,
                        "SCG_redundancy": 1.0,
                        "size": 10000,
                        "contigs": 2,
                        "N50": 5000,
                        "bin_score": 0.9,
                    }
                ]
            )
            checkm2_df = pd.DataFrame(
                [
                    {
                        "Name": "dastool-semibin_001",
                        "Completeness": 95.0,
                        "Contamination": 1.0,
                        "Genome_Size": 10000,
                        "Contig_N50": 5000,
                        "GC_Content": 0.51,
                    }
                ]
            )

            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                mags = PREPROCESS_MOD.build_mags(
                    dastool_summary, checkm2_df, contig2bin, None, None, None, None, defense_df, None
                )
                checkm2_all = PREPROCESS_MOD.build_checkm2_all(
                    str(tmp_path), checkm2_df, dastool_summary, contig2bin, None, None, None, defense_df, None
                )
                per_bin = PREPROCESS_MOD.build_mge_per_bin(contig2bin, None, None, defense_df, None)

            self.assertEqual(mags[0]["n_defense"], 2)
            self.assertEqual(checkm2_all[0]["n_defense"], 2)
            self.assertEqual(
                [d["contig"] for d in per_bin["dastool-semibin_001"]["defense"]],
                [spades_ctg, plain_ctg],
            )
            self.assertEqual(
                [d["contig"] for d in per_bin["_unbinned"]["defense"]],
                [unbinned_ctg],
            )

    def test_f3_control_two_token_contig_ids_unchanged(self):
        """Standard two-token contig_42 IDs join bins identically."""
        defense_df = pd.DataFrame(
            [{"sys_id": "contig_42_RM_Type_I_1", "type": "RM", "subtype": "RM_Type_I"}]
        )
        contig2bin = pd.DataFrame({"contig": ["contig_42"], "bin": ["dastool-semibin_001"]})
        dastool_summary = pd.DataFrame(
            [
                {
                    "bin": "dastool-semibin_001",
                    "SCG_completeness": 95.0,
                    "SCG_redundancy": 1.0,
                    "size": 5000,
                    "contigs": 1,
                    "N50": 5000,
                    "bin_score": 0.9,
                }
            ]
        )
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            mags = PREPROCESS_MOD.build_mags(
                dastool_summary, None, contig2bin, None, None, None, None, defense_df, None
            )
            per_bin = PREPROCESS_MOD.build_mge_per_bin(contig2bin, None, None, defense_df, None)
        self.assertEqual(mags[0]["n_defense"], 1)
        self.assertEqual(per_bin["dastool-semibin_001"]["defense"][0]["contig"], "contig_42")

    def test_f4_genomad_provirus_seq_name_joins_host_contig_and_bin(self):
        """geNomad '<contig>|provirus_<start>_<end>' entries join their host contig and bin."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            asm_dir = tmp_path / "assembly"
            asm_dir.mkdir(parents=True)
            (asm_dir / "tnf.tsv").write_text(
                "contig_42\t0.25\t0.25\t0.25\t0.25\n"
                "contig_43\t0.10\t0.40\t0.40\t0.10\n",
                encoding="utf-8",
            )
            genomad_dir = tmp_path / "mge" / "genomad"
            genomad_dir.mkdir(parents=True)
            (genomad_dir / "virus_summary.tsv").write_text(
                "seq_name\tlength\ttopology\tvirus_score\ttaxonomy\tn_hallmarks\n"
                "contig_42|provirus_101_2100\t2000\tProvirus\t0.96\t"
                "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota\t4\n",
                encoding="utf-8",
            )
            (genomad_dir / "taxonomy.tsv").write_text(
                "seq_name\tn_genes_with_taxonomy\tagreement\ttaxid\tlineage\n"
                "contig_42|provirus_101_2100\t5\t0.9\t12345\t"
                "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes\n",
                encoding="utf-8",
            )

            assembly_info = pd.DataFrame(
                {
                    "#seq_name": ["contig_42", "contig_43"],
                    "length": [10000, 4000],
                    "cov.": [25.0, 10.0],
                    "circ.": ["N", "N"],
                }
            )
            contig2bin = pd.DataFrame({"contig": ["contig_42"], "bin": ["dastool-semibin_001"]})
            dastool_summary = pd.DataFrame(
                [
                    {
                        "bin": "dastool-semibin_001",
                        "SCG_completeness": 92.0,
                        "SCG_redundancy": 1.5,
                        "size": 10000,
                        "contigs": 1,
                        "N50": 10000,
                        "bin_score": 0.88,
                    }
                ]
            )
            checkm2_df = pd.DataFrame(
                [
                    {
                        "Name": "dastool-semibin_001",
                        "Completeness": 92.0,
                        "Contamination": 1.5,
                        "Genome_Size": 10000,
                        "Contig_N50": 10000,
                        "GC_Content": 0.48,
                    }
                ]
            )
            virus_df = PREPROCESS_MOD.load_tsv(str(genomad_dir / "virus_summary.tsv"))

            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                mags = PREPROCESS_MOD.build_mags(
                    dastool_summary, checkm2_df, contig2bin, None, None, virus_df, None, None, None
                )
                checkm2_all = PREPROCESS_MOD.build_checkm2_all(
                    str(tmp_path), checkm2_df, dastool_summary, contig2bin, None, virus_df, None, None, None
                )
                per_bin = PREPROCESS_MOD.build_mge_per_bin(contig2bin, virus_df, None, None, None)
                explorer, _, _, _ = PREPROCESS_MOD.build_contig_explorer(
                    str(tmp_path), assembly_info, None, contig2bin, None, skip_tsne=True, skip_umap=True
                )

            self.assertEqual(mags[0]["n_virus"], 1)
            self.assertEqual(checkm2_all[0]["n_virus"], 1)
            self.assertIn("dastool-semibin_001", per_bin)
            self.assertEqual(per_bin["dastool-semibin_001"]["viruses"][0]["contig"], "contig_42")
            by_ctg = {c["id"]: c for c in explorer["contigs"]}
            self.assertEqual(by_ctg["contig_42"]["replicon"], "provirus")
            self.assertEqual(by_ctg["contig_42"]["genomad_domain"], "Viruses")

    def test_f4_control_standalone_viral_contig_unchanged(self):
        """Standalone viral contigs without '|provirus_' join bins and mark replicon='virus'."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            asm_dir = tmp_path / "assembly"
            asm_dir.mkdir(parents=True)
            (asm_dir / "tnf.tsv").write_text(
                "contig_99\t0.30\t0.20\t0.20\t0.30\n"
                "contig_100\t0.15\t0.35\t0.35\t0.15\n",
                encoding="utf-8",
            )
            genomad_dir = tmp_path / "mge" / "genomad"
            genomad_dir.mkdir(parents=True)
            (genomad_dir / "virus_summary.tsv").write_text(
                "seq_name\tlength\ttopology\tvirus_score\ttaxonomy\tn_hallmarks\n"
                "contig_99\t8000\tDirect terminal repeats\t0.99\tViruses;Duplodnaviria\t6\n",
                encoding="utf-8",
            )
            assembly_info = pd.DataFrame(
                {
                    "#seq_name": ["contig_99", "contig_100"],
                    "length": [8000, 3000],
                    "cov.": [30.0, 12.0],
                    "circ.": ["Y", "N"],
                }
            )
            contig2bin = pd.DataFrame({"contig": ["contig_99"], "bin": ["dastool-semibin_001"]})
            dastool_summary = pd.DataFrame(
                [
                    {
                        "bin": "dastool-semibin_001",
                        "SCG_completeness": 80.0,
                        "SCG_redundancy": 0.0,
                        "size": 8000,
                        "contigs": 1,
                        "N50": 8000,
                        "bin_score": 0.8,
                    }
                ]
            )
            virus_df = PREPROCESS_MOD.load_tsv(str(genomad_dir / "virus_summary.tsv"))
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                mags = PREPROCESS_MOD.build_mags(
                    dastool_summary, None, contig2bin, None, None, virus_df, None, None, None
                )
                explorer, _, _, _ = PREPROCESS_MOD.build_contig_explorer(
                    str(tmp_path), assembly_info, None, contig2bin, None, skip_tsne=True, skip_umap=True
                )
            self.assertEqual(mags[0]["n_virus"], 1)
            by_ctg = {c["id"]: c for c in explorer["contigs"]}
            self.assertEqual(by_ctg["contig_99"]["replicon"], "virus")

    def test_f5_antismash_gene_kinds_scoped_to_region_coordinates(self):
        """build_biosynthetic counts only CDS features overlapping each BGC region."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            as_dir = tmp_path / "metabolism" / "antismash" / "antismash_json"
            as_dir.mkdir(parents=True)

            antismash_payload = {
                "records": [
                    {
                        "id": "contig_1",
                        "areas": [
                            {"start": 1000, "end": 5000, "products": ["NRPS"]},
                            {"start": 20000, "end": 25000, "products": ["terpene"]},
                        ],
                        "features": [
                            {
                                "type": "CDS",
                                "location": "[2000:3500](+)",
                                "qualifiers": {"gene_kind": ["biosynthetic"]},
                            },
                            {
                                "type": "CDS",
                                "location": "[5000:6000](+)",
                                "qualifiers": {"gene_kind": ["other"]},
                            },
                            {
                                "type": "CDS",
                                "location": "[21000:22000](+)",
                                "qualifiers": {"gene_kind": ["transport"]},
                            },
                            {
                                "type": "CDS",
                                "location": "join{[23000:23500](-), [23600:24000](-)}",
                                "qualifiers": {"gene_kind": ["transport"]},
                            },
                            {
                                "type": "CDS",
                                "location": "[45000:46000](+)",
                                "qualifiers": {"gene_kind": ["regulatory"]},
                            },
                        ],
                    }
                ]
            }
            (as_dir / "antismash.json").write_text(json.dumps(antismash_payload), encoding="utf-8")
            contig2bin = pd.DataFrame({"contig": ["contig_1"], "bin": ["dastool-semibin_001"]})

            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                bio = PREPROCESS_MOD.build_biosynthetic(str(tmp_path), contig2bin)

            self.assertEqual(bio["n_regions"], 2)
            self.assertEqual(bio["regions"][0]["gene_kinds"], {"biosynthetic": 1})
            self.assertEqual(bio["regions"][1]["gene_kinds"], {"transport": 2})

    def test_f5_compound_origin_wrapping_cds_and_region_checked_per_part(self):
        """Compound join{...} CDS parts are tested individually rather than by bounding box; origin-crossing regions work."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            as_dir = tmp_path / "metabolism" / "antismash" / "antismash_json"
            as_dir.mkdir(parents=True)

            payload = {
                "records": [
                    {
                        "id": "circular_contig_1",
                        "areas": [
                            {"start": 1000, "end": 5000, "products": ["NRPS"]},
                            {"start": 9000, "end": 800, "products": ["bacteriocin"]},
                        ],
                        "features": [
                            {
                                "type": "CDS",
                                "location": "join{[9900:10000](+), [0:100](+)}",
                                "qualifiers": {"gene_kind": ["biosynthetic"]},
                            },
                            {
                                "type": "CDS",
                                "location": "[9500:9800](+)",
                                "qualifiers": {"gene_kind": ["transport"]},
                            },
                            {
                                "type": "CDS",
                                "location": "[200:600](-)",
                                "qualifiers": {"gene_kind": ["regulatory"]},
                            },
                        ],
                    }
                ]
            }
            (as_dir / "antismash.json").write_text(json.dumps(payload), encoding="utf-8")
            contig2bin = pd.DataFrame({"contig": ["circular_contig_1"], "bin": ["dastool-semibin_001"]})

            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                bio = PREPROCESS_MOD.build_biosynthetic(str(tmp_path), contig2bin)

            self.assertEqual(bio["n_regions"], 2)
            # Linear region [1000, 5000) does not overlap [9900, 10000) or [0, 100)
            self.assertEqual(bio["regions"][0]["gene_kinds"], {})
            self.assertEqual(bio["regions"][0]["length"], 4000)
            # Origin-crossing region [9000, 10000) U [0, 800) overlaps all 3 features once each
            self.assertEqual(
                bio["regions"][1]["gene_kinds"],
                {"biosynthetic": 1, "transport": 1, "regulatory": 1},
            )
            self.assertEqual(bio["regions"][1]["length"], 1800)

    def test_f5_control_single_region_spanning_all_cds_unchanged(self):
        """A single region spanning all CDS features counts every gene_kind on that contig."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            as_dir = tmp_path / "metabolism" / "antismash" / "antismash_json"
            as_dir.mkdir(parents=True)
            payload = {
                "records": [
                    {
                        "id": "contig_2",
                        "areas": [{"start": 0, "end": 10000, "products": ["T1PKS"]}],
                        "features": [
                            {
                                "type": "CDS",
                                "location": "[500:2500](+)",
                                "qualifiers": {"gene_kind": ["biosynthetic"]},
                            },
                            {
                                "type": "CDS",
                                "location": "[3000:4200](-)",
                                "qualifiers": {"gene_kind": ["transport"]},
                            },
                        ],
                    }
                ]
            }
            (as_dir / "antismash.json").write_text(json.dumps(payload), encoding="utf-8")
            contig2bin = pd.DataFrame({"contig": ["contig_2"], "bin": ["dastool-semibin_001"]})
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                bio = PREPROCESS_MOD.build_biosynthetic(str(tmp_path), contig2bin)
            self.assertEqual(bio["regions"][0]["gene_kinds"], {"biosynthetic": 1, "transport": 1})

    def test_f6_prokka_tsv_with_gff_and_unreconciled_bakta_extra_guard(self):
        """Prokka 7-col TSV joins sibling GFF coords; main() skips unreconciled bakta/extra."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)

            # Part A: Prokka annotation.tsv + sibling annotation.gff
            prokka_dir = tmp_path / "prokka_results" / "annotation" / "prokka"
            prokka_dir.mkdir(parents=True)
            prokka_tsv = prokka_dir / "annotation.tsv"
            prokka_tsv.write_text(
                "locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG\tproduct\n"
                "PROKKA_00001\tCDS\t990\tdnaA\t\tCOG0593\tChromosomal replication initiator protein DnaA\n"
                "PROKKA_00002\tCDS\t600\t\t\t\thypothetical protein\n"
                "PROKKA_00003\ttRNA\t76\t\t\t\ttRNA-Ala(GGC)\n"
                "PROKKA_00004\tmisc_RNA\t110\t\t\t\t6S RNA\n",
                encoding="utf-8",
            )
            prokka_gff = prokka_dir / "annotation.gff"
            prokka_gff.write_text(
                "##gff-version 3\n"
                "##sequence-region contig_1 1 5000\n"
                "contig_1\tProdigal\tCDS\t100\t1089\t.\t+\t0\tID=PROKKA_00001;locus_tag=PROKKA_00001;gene=dnaA\n"
                "contig_1\tProdigal\tCDS\t1200\t1799\t.\t-\t0\tID=PROKKA_00002;locus_tag=PROKKA_00002\n"
                "contig_1\tAragorn\ttRNA\t2000\t2075\t.\t+\t.\tID=PROKKA_00003;locus_tag=PROKKA_00003\n"
                "contig_1\tInfernal\tmisc_RNA\t2200\t2309\t.\t-\t.\tID=PROKKA_00004;locus_tag=PROKKA_00004\n"
                "##FASTA\n"
                ">contig_1\n"
                "ACGTACGT\n",
                encoding="utf-8",
            )

            genes, n_loaded = GENES_MOD.load_bakta(str(prokka_tsv))
            self.assertEqual(n_loaded, 4)
            self.assertEqual(len(genes["contig_1"]), 4)
            self.assertEqual(genes["contig_1"][0]["id"], "PROKKA_00001")
            self.assertEqual(genes["contig_1"][0]["s"], 100)
            self.assertEqual(genes["contig_1"][0]["e"], 1089)
            self.assertEqual(genes["contig_1"][0]["d"], 1)
            self.assertEqual(genes["contig_1"][3]["t"], "ncrna")

            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                annot_df = PREPROCESS_MOD.load_annotation_tsv(str(tmp_path / "prokka_results"))
            self.assertEqual(len(annot_df), 4)
            self.assertEqual(list(annot_df["contig"]), ["contig_1"] * 4)

            contig2bin = pd.DataFrame({"contig": ["contig_1"], "bin": ["dastool-semibin_001"]})
            dastool_summary = pd.DataFrame(
                [
                    {
                        "bin": "dastool-semibin_001",
                        "SCG_completeness": 90.0,
                        "SCG_redundancy": 1.0,
                        "size": 5000,
                        "contigs": 1,
                        "N50": 5000,
                        "bin_score": 0.85,
                    }
                ]
            )
            with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
                mags = PREPROCESS_MOD.build_mags(
                    dastool_summary, None, contig2bin, None, None, None, None, None, None, annot_df
                )
            self.assertEqual(mags[0]["n_genes"], 2)
            self.assertEqual(mags[0]["hypo_pct"], 50.0)
            self.assertEqual(mags[0]["n_trna"], 1)
            self.assertEqual(mags[0]["n_ncrna"], 1)

            # Part B: preprocess.py main() ignores unreconciled bakta/extra without locus_tag_map.tsv
            res_dir = tmp_path / "bakta_results"
            out_dir = tmp_path / "viz_out"
            (res_dir / "assembly").mkdir(parents=True)
            (res_dir / "annotation" / "bakta" / "basic").mkdir(parents=True)
            (res_dir / "annotation" / "bakta" / "extra").mkdir(parents=True)
            (res_dir / "mapping").mkdir(parents=True)

            (res_dir / "assembly" / "assembly_info.txt").write_text(
                "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
                "contig_1\t5000\t20.0\tN\tN\t1\t*\t1\n",
                encoding="utf-8",
            )
            (res_dir / "annotation" / "bakta" / "basic" / "annotation.tsv").write_text(
                "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                "contig_1\tcds\t100\t999\t+\tBASIC_00001\tdnaA\treplication initiator\t\n",
                encoding="utf-8",
            )
            (res_dir / "annotation" / "bakta" / "extra" / "annotation.tsv").write_text(
                "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                "contig_1\tcds\t100\t999\t+\tEXTRA_00001\tdnaA\treplication initiator\t\n",
                encoding="utf-8",
            )
            (res_dir / "mapping" / "gene_depths.tsv").write_text(
                "locus_tag\tcontig\tstart\tend\tmean_depth\n"
                "BASIC_00001\tcontig_1\t100\t999\t42.5\n",
                encoding="utf-8",
            )

            subprocess.run(
                [
                    sys.executable, "-B", str(PREPROCESS_DIR / "preprocess.py"),
                    "--results", str(res_dir),
                    "--output", str(out_dir),
                    "--skip-tsne", "--skip-umap", "--no-fastsne-viewer",
                ],
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
                capture_output=True,
                text=True,
                check=True,
            )
            with gzip.open(out_dir / "genes.part-000.json.gz", "rt", encoding="utf-8") as gz:
                shard = json.load(gz)
            feat = shard["contig_1"][0]
            self.assertEqual(feat["id"], "BASIC_00001")
            self.assertEqual(feat["dp"], 42.5)

    def test_f6_control_reconciled_bakta_extra_with_locus_tag_map_selected(self):
        """When locus_tag_map.tsv exists in bakta/extra, main() selects extra and merges depths."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            res_dir = tmp_path / "bakta_reconciled"
            out_dir = tmp_path / "viz_out"
            (res_dir / "assembly").mkdir(parents=True)
            (res_dir / "annotation" / "bakta" / "basic").mkdir(parents=True)
            (res_dir / "annotation" / "bakta" / "extra").mkdir(parents=True)
            (res_dir / "mapping").mkdir(parents=True)

            (res_dir / "assembly" / "assembly_info.txt").write_text(
                "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
                "contig_1\t5000\t20.0\tN\tN\t1\t*\t1\n",
                encoding="utf-8",
            )
            (res_dir / "annotation" / "bakta" / "basic" / "annotation.tsv").write_text(
                "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                "contig_1\tcds\t100\t999\t+\tBASIC_00001\t\thypothetical protein\t\n",
                encoding="utf-8",
            )
            (res_dir / "annotation" / "bakta" / "extra" / "annotation.tsv").write_text(
                "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                "contig_1\tcds\t100\t999\t+\tBASIC_00001\tdnaA\tChromosomal replication initiator\t\n",
                encoding="utf-8",
            )
            (res_dir / "annotation" / "bakta" / "extra" / "locus_tag_map.tsv").write_text(
                "extra_locus_tag\tbasic_locus_tag\nEXTRA_00001\tBASIC_00001\n",
                encoding="utf-8",
            )
            (res_dir / "mapping" / "gene_depths.tsv").write_text(
                "locus_tag\tcontig\tstart\tend\tmean_depth\n"
                "BASIC_00001\tcontig_1\t100\t999\t42.5\n",
                encoding="utf-8",
            )

            subprocess.run(
                [
                    sys.executable, "-B", str(PREPROCESS_DIR / "preprocess.py"),
                    "--results", str(res_dir),
                    "--output", str(out_dir),
                    "--skip-tsne", "--skip-umap", "--no-fastsne-viewer",
                ],
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
                capture_output=True,
                text=True,
                check=True,
            )
            with gzip.open(out_dir / "genes.part-000.json.gz", "rt", encoding="utf-8") as gz:
                shard = json.load(gz)
            feat = shard["contig_1"][0]
            self.assertEqual(feat["id"], "BASIC_00001")
            self.assertEqual(feat["p"], "Chromosomal replication initiator")
            self.assertEqual(feat["dp"], 42.5)


if __name__ == "__main__":
    unittest.main()
