#!/usr/bin/env python3
"""Regression and control tests for mag_analysis storeDir failure propagation (code-03 F7).

Verifies that storeDir-backed Nextflow processes propagate non-zero tool exit codes and
missing-database errors without writing placeholder output files, while legitimate
empty-input and zero-hit cases continue to exit 0 with valid header tables.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

REPO_ROOT = Path(__file__).resolve().parents[1]
MAG_DIR = REPO_ROOT / "mag_analysis"
MODULES_DIR = MAG_DIR / "modules"
MGE_TEST_DIR = MAG_DIR / "tests" / "mge_module"

NXF_PKG = Path("/opt/miniforge3/pkgs/nextflow-25.10.4-h2a3209d_0")
NXF_TEMPLATE = NXF_PKG / "bin" / "nextflow"
NXF_DIST = NXF_PKG / "share" / "nextflow" / "dist"


def extract_nf_process_block(nf_path: Path, process_name: str) -> str:
    text = nf_path.read_text(encoding="utf-8")
    m = re.search(rf"\bprocess\s+{re.escape(process_name)}\s*\{{", text)
    if not m:
        raise ValueError(f"Process {process_name} not present in {nf_path}")
    depth = 0
    start = m.start()
    for idx in range(m.end() - 1, len(text)):
        ch = text[idx]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return text[start : idx + 1]
    raise ValueError(f"Unbalanced braces for {process_name} in {nf_path}")


def extract_nf_bash_body(process_block: str) -> str:
    m = re.search(r'"""\n(.*?)\n\s*"""', process_block, flags=re.DOTALL)
    if not m:
        raise ValueError("Triple-quoted script block not present")
    body = m.group(1)
    return body.replace("\\$", "$").replace("\\\\", "\\")


def resolve_nextflow(dest_dir: Path) -> tuple[Path | None, str]:
    jdk_bins = sorted(Path("/opt/miniforge3/pkgs").glob("openjdk-*/bin"))
    extra_path = f"{jdk_bins[-1]}:/opt/miniforge3/bin" if jdk_bins else "/opt/miniforge3/bin"
    which_nxf = shutil.which("nextflow")
    if which_nxf:
        return Path(which_nxf), extra_path
    if NXF_TEMPLATE.is_file() and NXF_DIST.is_dir():
        dest_bin = dest_dir / "nextflow"
        lines = [
            f"NXF_DIST={NXF_DIST}" if ln.startswith("NXF_DIST=") else ln
            for ln in NXF_TEMPLATE.read_text(encoding="utf-8").splitlines()
        ]
        dest_bin.write_text("\n".join(lines) + "\n", encoding="utf-8")
        dest_bin.chmod(0o755)
        return dest_bin, extra_path
    return None, extra_path


class MagStoreDirFailurePropagationTests(unittest.TestCase):
    """Tests that failed tool executions never emit placeholder files into storeDir."""

    def _run_process_script(
        self,
        nf_path: Path,
        proc_name: str,
        subs: dict[str, str],
        files: dict[str, str],
        stubs: dict[str, str],
    ) -> tuple[subprocess.CompletedProcess[str], Path]:
        block = extract_nf_process_block(nf_path, proc_name)
        self.assertIn("storeDir params.store_dir", block)
        script_body = extract_nf_bash_body(block)
        for key, val in subs.items():
            script_body = script_body.replace(key, val)

        tmp_dir = Path(tempfile.mkdtemp(prefix=f"danaseq_{proc_name.lower()}_"))
        self.addCleanup(shutil.rmtree, tmp_dir, True)
        stub_dir = tmp_dir / "stubs"
        work_dir = tmp_dir / "work"
        stub_dir.mkdir()
        work_dir.mkdir()

        for tool_name, stub_body in stubs.items():
            stub_bin = stub_dir / tool_name
            stub_bin.write_text(stub_body, encoding="utf-8")
            stub_bin.chmod(0o755)

        for rel, content in files.items():
            dest = work_dir / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(content, encoding="utf-8")

        res = subprocess.run(
            ["bash", "-e", "-c", script_body],
            cwd=str(work_dir),
            env={**os.environ, "PATH": f"{stub_dir}:{os.environ['PATH']}"},
            capture_output=True,
            text=True,
        )
        return res, work_dir

    def _magscot(self, inputs, marker="M1", scorer=None, hmmsearch=None):
        files = {
            "asm.fa": ">c1\nACGT\n>c2\nTGCA\n",
            "db/gtdbtk_rel207_tigrfam.hmm": "stub\n",
            "db/gtdbtk_rel207_Pfam-A.hmm": "stub\n",
            "db/profiles/gtdb_rel207_default_markers.tsv": "marker\tset\nM1\tA\n",
            **inputs,
        }
        project = MAG_DIR
        if scorer is not None:
            files["project/bin/magscot.py"] = scorer
            project = "project"
        return self._run_process_script(
            MODULES_DIR / "binning.nf", "MAGSCOT_CONSENSUS",
            {"${task.cpus}": "1", "${params.magscot_hmm_dir}": "db",
             "${assembly}": "asm.fa", "${files_csv}": ",".join(inputs),
             "${labels_csv}": ",".join(inputs), "${projectDir}": str(project)},
            files,
            {"prodigal": '#!/bin/bash\nwhile [[ $1 != -a ]]; do shift; done\nprintf ">c1_1\\nMA\\n" > "$2"\n',
             "hmmsearch": hmmsearch or f'#!/bin/bash\nprintf "c1_1 - {marker} {marker} 1e-10\\n" > "$2"\n'},
        )

    def test_magscot_single_nonempty_binner_cannot_supply_consensus(self):
        res, work = self._magscot({"one.tsv": "c1\tb1\n", "two.tsv": ""})
        self.assertNotEqual(res.returncode, 0, res.stderr)
        self.assertIn("requires at least two nonempty binners", res.stderr)
        self.assertFalse((work / "magscot_bins.tsv").exists())

    def test_magscot_missing_mapping_after_zero_exit_is_a_failure(self):
        res, work = self._magscot({"one.tsv": "c1\tb1\n", "two.tsv": "c2\tb2\n"}, scorer="pass\n")
        self.assertNotEqual(res.returncode, 0, res.stderr)
        self.assertIn("missing refined mapping", res.stderr)
        self.assertFalse((work / "magscot_bins.tsv").exists())

    def test_magscot_rejects_invalid_mapping(self):
        inputs = {"one.tsv": "c1\tb1\n", "two.tsv": "c2\tb2\n"}
        for content in ["unexpected header\n", "binnew\tcontig\nb1\n"]:
            with self.subTest(content=content):
                scorer = f'open("magscot_out.refined.contig_to_bin.out", "w").write({content!r})'
                res, work = self._magscot(inputs, scorer=scorer)
                self.assertNotEqual(res.returncode, 0, res.stderr)
                self.assertIn("invalid refined mapping", res.stderr)
                self.assertFalse((work / "status.tsv").exists())

    def test_magscot_requires_hmmsearch_outputs(self):
        inputs = {"one.tsv": "c1\tb1\n", "two.tsv": "c2\tb2\n"}
        res, work = self._magscot(inputs, hmmsearch="#!/bin/bash\nexit 0\n")
        self.assertNotEqual(res.returncode, 0, res.stderr)
        self.assertIn("tigr.tblout", res.stderr)
        self.assertFalse((work / "status.tsv").exists())

    def test_magscot_accepts_completed_searches_without_hits(self):
        inputs = {"one.tsv": "c1\tb1\n", "two.tsv": "c2\tb2\n"}
        res, work = self._magscot(inputs, hmmsearch='#!/bin/bash\necho "# no hits" > "$2"\n')
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual((work / "status.tsv").read_text(), "status\nno_qualifying_bins\n")

    def test_magscot_propagates_scorer_exit(self):
        inputs = {"one.tsv": "c1\tb1\n", "two.tsv": "c2\tb2\n"}
        res, work = self._magscot(inputs, scorer="raise SystemExit(17)\n")
        self.assertEqual(res.returncode, 17, res.stderr)
        self.assertFalse((work / "status.tsv").exists())
        self.assertFalse((work / "magscot_bins.tsv").exists())

    def test_magscot_reports_qualified_and_empty_results(self):
        for inputs, marker, status in [
            ({"one.tsv": "", "two.tsv": ""}, "M1", "no_input_bins"),
            ({"one.tsv": "c1\tb1\n", "two.tsv": "c2\tb2\n"}, "M9", "no_qualifying_bins"),
            ({"one.tsv": "c1\tb1\n", "two.tsv": "c2\tb2\n"}, "M1", "ok"),
        ]:
            with self.subTest(status=status):
                res, work = self._magscot(inputs, marker)
                self.assertEqual(res.returncode, 0, res.stderr)
                self.assertEqual((work / "status.tsv").read_text(), f"status\n{status}\n")
                self.assertEqual((work / "magscot_bins.tsv").read_text(),
                                 "c1\tmagscot_001\n" if status == "ok" else "")
                if status == "ok":
                    self.assertEqual((work / "bins/magscot_001.fa").read_text(), ">c1\nACGT\n")

    def test_storedir_processes_propagate_nonzero_tool_exit_without_placeholders(self):
        """Wrapped tool exit 137 propagates non-zero and leaves zero output placeholders."""
        cases = [
            (
                MODULES_DIR / "mge.nf",
                "GENOMAD_CLASSIFY",
                "genomad",
                {
                    "${task.cpus}": "1",
                    "${assembly}": "asm.fa",
                    "${db_path}": "db",
                    "${emit_failed_status}": "false",
                },
                {"asm.fa": ">c1\nACGTACGT\n"},
                ["virus_summary.tsv", "plasmid_summary.tsv", "provirus.tsv", "genomad_summary.tsv", "status.tsv"],
            ),
            (
                MODULES_DIR / "mge.nf",
                "CHECKV_QUALITY",
                "checkv",
                {
                    "${task.cpus}": "1",
                    "${virus_fasta}": "virus.fna",
                    "${db_path}": "db",
                    "${emit_failed_status}": "false",
                },
                {"virus.fna": ">v1\nACGTACGT\n"},
                ["quality_summary.tsv", "viruses.fna", "proviruses.fna", "status.tsv"],
            ),
            (
                MODULES_DIR / "mge.nf",
                "ISLANDPATH_DIMOB",
                "islandpath_dimob.py",
                {
                    "${task.cpus}": "1",
                    "${gff}": "a.gff",
                    "${assembly}": "asm.fa",
                    "${faa}": "a.faa",
                    "${hmm_db}": "missing.hmm",
                },
                {"a.gff": "##gff-version 3\n", "asm.fa": ">c1\nACGT\n", "a.faa": ">p1\nMKTA\n"},
                ["genomic_islands.tsv"],
            ),
            (
                MODULES_DIR / "mge.nf",
                "MACSYFINDER",
                "macsyfinder",
                {"${task.cpus}": "1", "${proteins}": "a.faa", "${models_dir}": "models"},
                {"a.faa": ">p1\nMKTA\n"},
                ["all_systems.tsv", "all_systems.txt"],
            ),
            (
                MODULES_DIR / "mge.nf",
                "DEFENSEFINDER",
                "parallel_defensefinder.py",
                {"${task.cpus}": "1", "${proteins}": "a.faa", "${gff}": "a.gff", "${models_opt}": "--models-dir m"},
                {"a.faa": ">p1\nMKTA\n", "a.gff": "##gff-version 3\n"},
                ["systems.tsv", "genes.tsv", "hmmer.tsv"],
            ),
            (
                MODULES_DIR / "eukaryotic.nf",
                "TIARA_CLASSIFY",
                "tiara",
                {"${task.cpus}": "1", "${contigs}": "asm.fa", "${min_len}": "100"},
                {"asm.fa": ">c1\n" + "ACGT" * 50 + "\n"},
                ["tiara_output.tsv"],
            ),
            (
                MODULES_DIR / "eukaryotic.nf",
                "WHOKARYOTE_CLASSIFY",
                "whokaryote.py",
                {"${task.cpus}": "1", "${contigs}": "asm.fa", "${gff_arg}": "", "${min_len}": "100"},
                {"asm.fa": ">c1\n" + "ACGT" * 50 + "\n"},
                ["whokaryote_classifications.tsv"],
            ),
            (
                MODULES_DIR / "eukaryotic.nf",
                "METAEUK_PREDICT",
                "metaeuk",
                {
                    "${task.cpus}": "1",
                    "${contigs}": "asm.fa",
                    "${tiara_tsv}": "tiara.tsv",
                    "${whokaryote_tsv}": "whok.tsv",
                    "${params.metaeuk_db}": "db",
                    "${mem_limit}": "1G",
                    "${min_len}": "20",
                    "${max_intron}": "1000",
                },
                {
                    "asm.fa": ">c1\nACGTACGT\n",
                    "tiara.tsv": "sequence_id\tclass_fst_stage\tclass_snd_stage\nc1\teukarya\tunknown\n",
                    "whok.tsv": "contig\tprediction\nc1\teukaryote\n",
                },
                ["metaeuk_proteins.fas", "metaeuk_codon.fas", "metaeuk.gff", "metaeuk_headers.tsv"],
            ),
            (
                MODULES_DIR / "eukaryotic.nf",
                "MARFERRET_CLASSIFY",
                "diamond",
                {"${task.cpus}": "1", "${proteins}": "a.faa", "${marferret_db}": "mdb"},
                {"a.faa": ">p1\nMKTA\n", "mdb/ref.dmnd": "x"},
                ["marferret_proteins.tsv", "marferret_contigs.tsv"],
            ),
            (
                MODULES_DIR / "metabolism.nf",
                "KOFAMSCAN",
                "exec_annotation",
                {"${task.cpus}": "1", "${proteins}": "a.faa", "${kofam_profiles}": "p", "${ko_list}": "k"},
                {"a.faa": ">p1\nMKTA\n"},
                ["kofamscan_results.tsv"],
            ),
            (
                MODULES_DIR / "metabolism.nf",
                "EMAPPER",
                "emapper.py",
                {"${cpus}": "1", "${batch_size}": "50000", "${proteins}": "a.faa", "${eggnog_db}": "db"},
                {"a.faa": ">p1\nMKTA\n"},
                ["emapper_results.emapper.annotations"],
            ),
            (
                MODULES_DIR / "metabolism.nf",
                "DBCAN",
                "run_dbcan",
                {"${task.cpus}": "1", "${proteins}": "a.faa", "${dbcan_db}": "db"},
                {"a.faa": ">p1\nMKTA\n"},
                ["overview.tsv"],
            ),
            (
                MODULES_DIR / "taxonomy.nf",
                "KAIJU_CONTIG_CLASSIFY",
                "kaiju",
                {"${task.cpus}": "1", "${contigs}": "asm.fa", "${db_dir}": "kdb"},
                {"asm.fa": ">c1\nACGT\n", "kdb/db.fmi": "x", "kdb/nodes.dmp": "x", "kdb/names.dmp": "x"},
                ["kaiju_contigs.tsv"],
            ),
            (
                MODULES_DIR / "taxonomy.nf",
                "KAIJU_CLASSIFY",
                "kaiju",
                {"${task.cpus}": "1", "${proteins}": "a.faa", "${gff}": "a.gff", "${db_dir}": "kdb"},
                {
                    "a.faa": ">p1\nMKTA\n",
                    "a.gff": "##gff-version 3\n",
                    "kdb/db.fmi": "x",
                    "kdb/nodes.dmp": "x",
                    "kdb/names.dmp": "x",
                },
                ["kaiju_genes.tsv", "kaiju_protein_contigs.tsv"],
            ),
            (
                MODULES_DIR / "taxonomy.nf",
                "KRAKEN2_CLASSIFY",
                "kraken2",
                {"${task.cpus}": "1", "${contigs}": "asm.fa", "${db_dir}": "krdb", "${confidence}": "0.0"},
                {"asm.fa": ">c1\nACGT\n", "krdb/hash.k2d": "x", "krdb/nodes.dmp": "x", "krdb/names.dmp": "x"},
                ["kraken2_contigs.tsv", "kraken2_report.txt"],
            ),
            (
                MODULES_DIR / "taxonomy.nf",
                "SENDSKETCH_CLASSIFY",
                "sendsketch.sh",
                {"${contigs}": "asm.fa", "${address}": "https://example.invalid"},
                {"asm.fa": ">c1\nACGT\n"},
                ["sendsketch_contigs.tsv"],
            ),
        ]

        for nf_path, proc_name, stub_tool, subs, files, expected_outputs in cases:
            with self.subTest(process=proc_name):
                res, work_dir = self._run_process_script(
                    nf_path,
                    proc_name,
                    subs,
                    files,
                    {stub_tool: "#!/bin/sh\nexit 137\n"},
                )
                self.assertEqual(
                    res.returncode,
                    137,
                    f"{proc_name} returned {res.returncode} instead of 137: {res.stderr}",
                )
                for out_file in expected_outputs:
                    self.assertFalse(
                        (work_dir / out_file).exists(),
                        f"{proc_name} created placeholder {out_file} on failure",
                    )

    def test_missing_db_guards_fail_without_poisoning_storedir(self):
        """Missing database files return exit 1 without writing placeholder TSVs."""
        cases = [
            (
                MODULES_DIR / "eukaryotic.nf",
                "MARFERRET_CLASSIFY",
                {"${task.cpus}": "1", "${proteins}": "a.faa", "${marferret_db}": "empty_db"},
                {"a.faa": ">p1\nMKTA\n", "empty_db/.keep": ""},
                ["marferret_proteins.tsv", "marferret_contigs.tsv"],
            ),
            (
                MODULES_DIR / "taxonomy.nf",
                "KAIJU_CONTIG_CLASSIFY",
                {"${task.cpus}": "1", "${contigs}": "asm.fa", "${db_dir}": "empty_db"},
                {"asm.fa": ">c1\nACGT\n", "empty_db/.keep": ""},
                ["kaiju_contigs.tsv"],
            ),
            (
                MODULES_DIR / "taxonomy.nf",
                "KRAKEN2_CLASSIFY",
                {"${task.cpus}": "1", "${contigs}": "asm.fa", "${db_dir}": "empty_db", "${confidence}": "0.0"},
                {"asm.fa": ">c1\nACGT\n", "empty_db/.keep": ""},
                ["kraken2_contigs.tsv", "kraken2_report.txt"],
            ),
        ]
        for nf_path, proc_name, subs, files, expected_outputs in cases:
            with self.subTest(process=proc_name):
                res, work_dir = self._run_process_script(nf_path, proc_name, subs, files, {})
                self.assertEqual(res.returncode, 1)
                for out_file in expected_outputs:
                    self.assertFalse((work_dir / out_file).exists())

    def test_genomad_storedir_mode_under_real_nextflow_leaves_store_clean(self):
        """Under Nextflow with --store_dir set, failed geNomad leaves store_dir/mge/genomad empty."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            nxf_bin, extra_path = resolve_nextflow(tmp_path)
            if nxf_bin is None:
                self.skipTest("Nextflow binary not available")

            with gzip.open(tmp_path / "asm.fa.gz", "wt", encoding="utf-8") as gz:
                gz.write(">host1\n" + "A" * 4000 + "\n")
            (tmp_path / "c2b.tsv").write_text("host1\tbinA\n", encoding="utf-8")
            store_dir = tmp_path / "store"
            out_dir = tmp_path / "out"

            env = {
                **os.environ,
                "FAKE_GENOMAD_MODE": "fail",
                "NXF_HOME": str(tmp_path / ".nxf"),
                "NXF_OFFLINE": "true",
                "PATH": f"{MGE_TEST_DIR / 'fakebin'}:{MAG_DIR / 'bin'}:{extra_path}:{os.environ['PATH']}",
            }
            res = subprocess.run(
                [
                    str(nxf_bin),
                    "-q",
                    "run",
                    str(MGE_TEST_DIR / "test.nf"),
                    "-work-dir",
                    str(tmp_path / "work"),
                    "--assembly",
                    str(tmp_path / "asm.fa.gz"),
                    "--contig2bin",
                    str(tmp_path / "c2b.tsv"),
                    "--outdir",
                    str(out_dir),
                    "--store_dir",
                    str(store_dir),
                    "--run_checkv",
                    "false",
                ],
                cwd=str(tmp_path),
                env=env,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(res.returncode, 0)
            genomad_store = store_dir / "mge" / "genomad"
            stored_files = list(genomad_store.iterdir()) if genomad_store.exists() else []
            self.assertEqual(stored_files, [])

    def test_control_legitimate_empty_or_zero_hit_inputs_exit_zero_with_headers(self):
        """Empty inputs and valid zero-hit runs continue to exit 0 and emit valid headers."""
        # 1. CHECKV_QUALITY with 0 viral contigs
        res_cv, wd_cv = self._run_process_script(
            MODULES_DIR / "mge.nf",
            "CHECKV_QUALITY",
            {
                "${task.cpus}": "1",
                "${virus_fasta}": "virus.fna",
                "${db_path}": "db",
                "${emit_failed_status}": "false",
            },
            {"virus.fna": ""},
            {"checkv": "#!/bin/sh\necho 'CheckV v1.0.0'\nexit 0\n"},
        )
        self.assertEqual(res_cv.returncode, 0)
        self.assertTrue((wd_cv / "quality_summary.tsv").is_file())
        self.assertIn("no_hits", (wd_cv / "status.tsv").read_text(encoding="utf-8"))

        # 2. METAEUK_PREDICT when all contigs are prokaryotic (0 filtered eukaryotic contigs)
        res_me, wd_me = self._run_process_script(
            MODULES_DIR / "eukaryotic.nf",
            "METAEUK_PREDICT",
            {
                "${task.cpus}": "1",
                "${contigs}": "asm.fa",
                "${tiara_tsv}": "tiara.tsv",
                "${whokaryote_tsv}": "whok.tsv",
                "${params.metaeuk_db}": "db",
                "${mem_limit}": "1G",
                "${min_len}": "20",
                "${max_intron}": "1000",
            },
            {
                "asm.fa": ">c1\nACGTACGT\n",
                "tiara.tsv": "sequence_id\tclass_fst_stage\tclass_snd_stage\nc1\tbacteria\tNA\n",
                "whok.tsv": "contig\tprediction\nc1\tprokaryote\n",
            },
            {"metaeuk": "#!/bin/sh\nexit 137\n"},
        )
        self.assertEqual(res_me.returncode, 0)
        self.assertTrue((wd_me / "metaeuk_proteins.fas").exists())

        # 3. KOFAMSCAN with empty protein FASTA
        res_ko, wd_ko = self._run_process_script(
            MODULES_DIR / "metabolism.nf",
            "KOFAMSCAN",
            {"${task.cpus}": "1", "${proteins}": "empty.faa", "${kofam_profiles}": "p", "${ko_list}": "k"},
            {"empty.faa": ""},
            {"exec_annotation": "#!/bin/sh\nexit 137\n"},
        )
        self.assertEqual(res_ko.returncode, 0)
        self.assertTrue((wd_ko / "kofamscan_results.tsv").is_file())

        # 4. EMAPPER with non-empty proteins where Stage 1 succeeds with 0 hits
        res_em, wd_em = self._run_process_script(
            MODULES_DIR / "metabolism.nf",
            "EMAPPER",
            {"${cpus}": "1", "${batch_size}": "50000", "${proteins}": "a.faa", "${eggnog_db}": "db"},
            {"a.faa": ">p1\nMKTAYIAK\n"},
            {"emapper.py": "#!/bin/sh\nexit 0\n"},
        )
        self.assertEqual(res_em.returncode, 0)
        self.assertIn("#query\tseed_ortholog", (wd_em / "emapper_results.emapper.annotations").read_text(encoding="utf-8"))

    def test_control_antismash_and_gtdbtk_propagate_failure(self):
        """ANTISMASH and GTDBTK_CLASSIFY continue to propagate exit 137 with no outputs."""
        for nf_path, proc_name, tool, subs, files, outs in [
            (
                MODULES_DIR / "metabolism.nf",
                "ANTISMASH",
                "antismash",
                {
                    "${task.cpus}": "1",
                    "${db_arg}": "",
                    "${genefinding}": "--genefinding-tool prodigal-m",
                    "${assembly}": "asm.fa",
                },
                {"asm.fa": ">c1\n" + "ACGT" * 300 + "\n"},
                ["antismash_summary.tsv"],
            ),
            (
                MODULES_DIR / "phylogeny.nf",
                "GTDBTK_CLASSIFY",
                "gtdbtk",
                {
                    "${task.cpus}": "1",
                    "${params.gtdbtk_max_bin_bp}": "50000000",
                    "${params.gtdbtk_db}": "gdb",
                    "${params.gtdbtk_pplacer_cpus}": "1",
                },
                {"bins_1/bin1.fa": ">c1\nACGTACGT\n"},
                ["gtdbtk_taxonomy.tsv"],
            ),
        ]:
            with self.subTest(process=proc_name):
                res, wd = self._run_process_script(
                    nf_path, proc_name, subs, files, {tool: "#!/bin/sh\nexit 137\n"}
                )
                self.assertEqual(res.returncode, 137)
                for out_f in outs:
                    self.assertFalse((wd / out_f).exists())


if __name__ == "__main__":
    unittest.main()
