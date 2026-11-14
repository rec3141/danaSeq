#!/usr/bin/env python3
"""Regression and control tests for nanopore_assembly polish, read_type, GFA remap, and TNF masking."""

import contextlib
import gzip
import importlib.util
import io
import os
from pathlib import Path
import random
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
NANOPORE_DIR = REPO / "nanopore_assembly"
MAIN_NF = NANOPORE_DIR / "main.nf"
ASSEMBLY_NF = NANOPORE_DIR / "modules" / "assembly.nf"
NORMALIZE_ASSEMBLY_PY = NANOPORE_DIR / "bin" / "normalize_assembly.py"
TETRAMER_FREQS_C = NANOPORE_DIR / "bin" / "tetramer_freqs.c"

NXF_TEMPLATE = Path("/opt/miniforge3/pkgs/nextflow-25.10.4-h2a3209d_0/bin/nextflow")
NXF_DIST = Path("/opt/miniforge3/pkgs/nextflow-25.10.4-h2a3209d_0/share/nextflow/dist")
FLYE_BIN = Path("/opt/miniforge3/envs/flye/bin/flye")


def extract_nf_process_script(nf_path: Path, process_name: str, subs: dict[str, str]) -> str:
    """Render a Nextflow process script block with GString unescaping and variable substitutions."""
    text = nf_path.read_text(encoding="utf-8")
    proc_match = re.search(rf"^process\s+{re.escape(process_name)}\s*\{{", text, flags=re.MULTILINE)
    if not proc_match:
        raise ValueError(f"Process {process_name} not found in {nf_path}")
    sub = text[proc_match.end():]
    script_match = re.search(r'\bscript:\s*(?:.*?\n)*?\s*"""\n(.*?)\n\s*"""', sub, re.DOTALL)
    if not script_match:
        raise ValueError(f"Script block for {process_name} not found in {nf_path}")
    raw = script_match.group(1)
    rendered = re.sub(r"(?<!\\)\$\{([^}]+)\}", lambda m: subs.get(m.group(1), m.group(0)), raw)
    return rendered.replace(r"\$", "$").replace(r"\\", "\\").replace(r'\"', '"')


def write_executable(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def load_normalize_assembly():
    spec = importlib.util.spec_from_file_location("normalize_assembly", str(NORMALIZE_ASSEMBLY_PY))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def resolve_nextflow(dest_dir: Path) -> Path | None:
    which_nxf = shutil.which("nextflow")
    if which_nxf:
        return Path(which_nxf)
    if NXF_TEMPLATE.is_file() and NXF_DIST.is_dir():
        dest_bin = dest_dir / "nextflow"
        lines = [
            f"NXF_DIST={NXF_DIST}" if ln.startswith("NXF_DIST=") else ln
            for ln in NXF_TEMPLATE.read_text(encoding="utf-8").splitlines()
        ]
        write_executable(dest_bin, "\n".join(lines) + "\n")
        return dest_bin
    return None


class FlyePolishMetadataAndFallbackTests(unittest.TestCase):
    """FLYE_POLISH keeps info/GFA lengths consistent with polished sequences and handles staged symlinks."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="danaseq_f3_")
        self.tdir = Path(self.tmp.name)
        self.staging_src = self.tdir / "upstream"
        self.staging_src.mkdir()
        (self.staging_src / "assembly.fasta").write_text(
            ">contig_1\n" + "A" * 1500 + "\n", encoding="utf-8"
        )
        (self.staging_src / "assembly_info.txt").write_text(
            "#seq_name\tlength\tcov.\tcirc.\ncontig_1\t1500\t25\tN\n",
            encoding="utf-8",
        )
        (self.staging_src / "assembly_graph.gfa").write_text(
            "H\tVN:Z:1.0\n"
            f"S\tcontig_1\t{'A' * 1500}\tLN:i:1500\n"
            "S\tedge_1\tACGTACGT\n",
            encoding="utf-8",
        )
        (self.staging_src / "all_reads.fastq").write_text(
            "@r1\nACGT\n+\nIIII\n", encoding="utf-8"
        )

    def tearDown(self):
        self.tmp.cleanup()

    def _stage_symlinks(self, name: str) -> Path:
        wdir = self.tdir / name
        wdir.mkdir()
        for fn in ("assembly.fasta", "assembly_info.txt", "assembly_graph.gfa", "all_reads.fastq"):
            (wdir / fn).symlink_to(self.staging_src / fn)
        return wdir

    def test_polished_length_updates_info_and_gfa_without_mutating_upstream(self):
        rendered = extract_nf_process_script(
            ASSEMBLY_NF,
            "FLYE_POLISH",
            {
                "params.read_type": "nano-hq",
                "task.cpus": "2",
                "assembly": "assembly.fasta",
                "info": "assembly_info.txt",
                "graph": "assembly_graph.gfa",
                "reads": "all_reads.fastq",
            },
        )
        bin_dir = self.tdir / "bin_ok"
        bin_dir.mkdir()
        write_executable(
            bin_dir / "flye",
            "#!/usr/bin/env bash\n"
            "mkdir -p polish_out\n"
            "python3 -c 'open(\"polish_out/polished_1.fasta\",\"w\").write(\">contig_1\\n\"+\"G\"*1508+\"\\n\")'\n",
        )
        wdir = self._stage_symlinks("work_ok")
        env = dict(os.environ, PATH=f"{bin_dir}:{os.environ.get('PATH', '')}")
        res = subprocess.run(
            ["bash", "-ue", "-c", rendered], cwd=wdir, env=env, capture_output=True, text=True
        )
        self.assertEqual(res.returncode, 0, res.stderr)

        polished_seq = "".join(
            ln.strip()
            for ln in (wdir / "assembly.fasta").read_text(encoding="utf-8").splitlines()
            if not ln.startswith(">")
        )
        self.assertEqual(len(polished_seq), 1508)
        self.assertIn(
            "contig_1\t1508\t25\tN",
            (wdir / "assembly_info.txt").read_text(encoding="utf-8"),
        )
        gfa_lines = (wdir / "assembly_graph.gfa").read_text(encoding="utf-8").splitlines()
        s_by_id = {ln.split("\t")[1]: ln.split("\t") for ln in gfa_lines if ln.startswith("S\t")}
        self.assertEqual(len(s_by_id["contig_1"][2]), 1508)
        self.assertIn("LN:i:1508", s_by_id["contig_1"][3:])
        self.assertEqual(s_by_id["edge_1"][2], "ACGTACGT")

        # Upstream staged targets must remain untouched at 1500 bp.
        self.assertIn(
            "contig_1\t1500\t25\tN",
            (self.staging_src / "assembly_info.txt").read_text(encoding="utf-8"),
        )
        self.assertIn(
            "LN:i:1500",
            (self.staging_src / "assembly_graph.gfa").read_text(encoding="utf-8"),
        )

    def test_empty_polish_fallback_with_same_filename_succeeds(self):
        rendered = extract_nf_process_script(
            ASSEMBLY_NF,
            "FLYE_POLISH",
            {
                "params.read_type": "nano-hq",
                "task.cpus": "2",
                "assembly": "assembly.fasta",
                "info": "assembly_info.txt",
                "graph": "assembly_graph.gfa",
                "reads": "all_reads.fastq",
            },
        )
        bin_dir = self.tdir / "bin_empty"
        bin_dir.mkdir()
        write_executable(bin_dir / "flye", "#!/usr/bin/env bash\nmkdir -p polish_out\n")
        wdir = self._stage_symlinks("work_empty")
        env = dict(os.environ, PATH=f"{bin_dir}:{os.environ.get('PATH', '')}")
        res = subprocess.run(
            ["bash", "-ue", "-c", rendered], cwd=wdir, env=env, capture_output=True, text=True
        )
        self.assertEqual(res.returncode, 0, res.stderr)
        unpolished_seq = "".join(
            ln.strip()
            for ln in (wdir / "assembly.fasta").read_text(encoding="utf-8").splitlines()
            if not ln.startswith(">")
        )
        self.assertEqual(len(unpolished_seq), 1500)

    def test_control_unchanged_length_polish_from_distinct_input_names(self):
        rendered = extract_nf_process_script(
            ASSEMBLY_NF,
            "FLYE_POLISH",
            {
                "params.read_type": "nano-hq",
                "task.cpus": "2",
                "assembly": "draft.fasta",
                "info": "draft_info.txt",
                "graph": "draft_graph.gfa",
                "reads": "all_reads.fastq",
            },
        )
        wdir = self.tdir / "work_ctrl"
        wdir.mkdir()
        shutil.copy(self.staging_src / "assembly.fasta", wdir / "draft.fasta")
        shutil.copy(self.staging_src / "assembly_info.txt", wdir / "draft_info.txt")
        shutil.copy(self.staging_src / "assembly_graph.gfa", wdir / "draft_graph.gfa")
        shutil.copy(self.staging_src / "all_reads.fastq", wdir / "all_reads.fastq")
        bin_dir = self.tdir / "bin_ctrl"
        bin_dir.mkdir()
        write_executable(
            bin_dir / "flye",
            "#!/usr/bin/env bash\n"
            "mkdir -p polish_out\n"
            "cp draft.fasta polish_out/polished_1.fasta\n",
        )
        env = dict(os.environ, PATH=f"{bin_dir}:{os.environ.get('PATH', '')}")
        res = subprocess.run(
            ["bash", "-ue", "-c", rendered], cwd=wdir, env=env, capture_output=True, text=True
        )
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertIn(
            "contig_1\t1500\t25\tN",
            (wdir / "assembly_info.txt").read_text(encoding="utf-8"),
        )

    def test_flye_polish_module_under_real_nextflow_updates_lengths(self):
        nxf = resolve_nextflow(self.tdir)
        if nxf is None:
            self.skipTest("Nextflow is not available in this environment")
        bin_dir = self.tdir / "bin_nf"
        bin_dir.mkdir()
        write_executable(
            bin_dir / "flye",
            "#!/usr/bin/env bash\n"
            "mkdir -p polish_out\n"
            "python3 -c 'open(\"polish_out/polished_1.fasta\",\"w\").write(\">contig_1\\n\"+\"G\"*1508+\"\\n\")'\n",
        )
        nf_script = self.tdir / "run_polish.nf"
        nf_script.write_text(
            "nextflow.enable.dsl = 2\n"
            f"include {{ FLYE_POLISH }} from '{ASSEMBLY_NF}'\n"
            "workflow {\n"
            f"    FLYE_POLISH(file('{self.staging_src / 'assembly.fasta'}'),\n"
            f"                file('{self.staging_src / 'assembly_info.txt'}'),\n"
            f"                file('{self.staging_src / 'assembly_graph.gfa'}'),\n"
            f"                file('{self.staging_src / 'all_reads.fastq'}'))\n"
            "}\n",
            encoding="utf-8",
        )
        cfg = self.tdir / "nf.config"
        cfg.write_text(
            "conda.enabled = false\n"
            "params { outdir = 'out'; store_dir = null; read_type = 'nano-hq' }\n"
            "process { cpus = 1; memory = '256 MB' }\n",
            encoding="utf-8",
        )
        nxf_home = self.tdir / "nxf_home"
        nxf_home.mkdir()
        env = dict(
            os.environ,
            PATH=f"{bin_dir}:{os.environ.get('PATH', '')}",
            NXF_HOME=str(nxf_home),
            NXF_OFFLINE="true",
            NXF_ANSI_LOG="false",
            NXF_JVM_ARGS="-Xms64m -Xmx512m",
        )
        res = subprocess.run(
            [
                str(nxf),
                "-log", str(self.tdir / "nf.log"),
                "run", str(nf_script),
                "-c", str(cfg),
                "-w", str(self.tdir / "work"),
            ],
            cwd=self.tdir,
            env=env,
            capture_output=True,
            text=True,
            timeout=90,
        )
        self.assertEqual(res.returncode, 0, res.stdout + "\n" + res.stderr)
        pub_info = (self.tdir / "out" / "assembly" / "assembly_info.txt").read_text(encoding="utf-8")
        pub_gfa = (self.tdir / "out" / "assembly" / "assembly_graph.gfa").read_text(encoding="utf-8")
        self.assertIn("contig_1\t1508\t25\tN", pub_info)
        self.assertIn("LN:i:1508", pub_gfa)


class ReadTypePolishValidationTests(unittest.TestCase):
    """Polishing non-Flye assemblies without --read_type must fail parameter validation."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix="danaseq_f4_bins_")
        cls.bdir = Path(cls._tmp.name)
        cls.nextflow_bin = resolve_nextflow(cls.bdir)
        cls.nxf_home = cls.bdir / "nxf_home"
        cls.nxf_home.mkdir()
        if TETRAMER_FREQS_C.is_file():
            subprocess.run(
                ["gcc", "-O2", "-o", str(cls.bdir / "tetramer_freqs"), str(TETRAMER_FREQS_C), "-lm"],
                check=True,
            )

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_flye_polish_script_rejects_null_read_type(self):
        with tempfile.TemporaryDirectory(prefix="danaseq_f4_script_") as tmp:
            tdir = Path(tmp)
            rendered = extract_nf_process_script(
                ASSEMBLY_NF,
                "FLYE_POLISH",
                {
                    "params.read_type": "null",
                    "task.cpus": "1",
                    "assembly": "assembly.fasta",
                    "info": "assembly_info.txt",
                    "graph": "assembly_graph.gfa",
                    "reads": "reads.fq",
                },
            )
            res = subprocess.run(
                ["bash", "-ue", "-c", rendered], cwd=tdir, capture_output=True, text=True
            )
            self.assertNotEqual(res.returncode, 0)
            self.assertIn("read_type", res.stderr)

    def _run_nextflow(self, tdir: Path, extra_args: list[str], extra_bins: list[Path] = ()):
        if self.nextflow_bin is None:
            self.skipTest("Nextflow is not available in this environment")
        cfg = tdir / "override.config"
        cfg.write_text(
            "conda.enabled = false\n"
            "timeline.enabled = false\n"
            "report.enabled = false\n"
            "dag.enabled = false\n"
            "process {\n"
            "    withLabel: 'process_low'    { cpus = 1; memory = '256 MB' }\n"
            "    withLabel: 'process_medium' { cpus = 1; memory = '256 MB' }\n"
            "    withLabel: 'process_high'   { cpus = 1; memory = '256 MB' }\n"
            "}\n",
            encoding="utf-8",
        )
        path_parts = [str(p) for p in extra_bins] + [str(self.bdir)]
        if FLYE_BIN.parent.is_dir():
            path_parts.append(str(FLYE_BIN.parent))
        path_parts.append(os.environ.get("PATH", ""))
        env = dict(
            os.environ,
            PATH=":".join(p for p in path_parts if p),
            NXF_HOME=str(self.nxf_home),
            NXF_OFFLINE="true",
            NXF_ANSI_LOG="false",
            NXF_JVM_ARGS="-Xms64m -Xmx512m",
        )
        return subprocess.run(
            [
                str(self.nextflow_bin),
                "-log",
                str(tdir / "nextflow.log"),
                "run",
                str(MAIN_NF),
                "-c",
                str(cfg),
                "-w",
                str(tdir / "work"),
                *extra_args,
            ],
            cwd=tdir,
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )

    def test_metamdbg_with_polish_true_and_no_read_type_exits_nonzero(self):
        with tempfile.TemporaryDirectory(prefix="danaseq_f4_nf_") as tmp:
            tdir = Path(tmp)
            reads_dir = tdir / "reads_in"
            reads_dir.mkdir()
            stub_bin = tdir / "stub_bin"
            stub_bin.mkdir()

            rng = random.Random(42)
            contig_seq = "".join(rng.choice("ACGT") for _ in range(1500))
            with gzip.open(reads_dir / "sample1.fastq.gz", "wb") as gz:
                for i in range(1, 10):
                    rseq = contig_seq[(i * 40):(i * 40) + 300]
                    gz.write(f"@r{i}\n{rseq}\n+\n{'I' * len(rseq)}\n".encode("ascii"))

            contigs_gz = tdir / "stub_contigs.fasta.gz"
            with gzip.open(contigs_gz, "wb") as gz:
                gz.write(f">ctg1 length=1500 circular=no coverage=20\n{contig_seq}\n".encode("ascii"))
            write_executable(
                stub_bin / "metaMDBG",
                "#!/usr/bin/env bash\n"
                'if [[ "${1:-}" == "asm" ]]; then\n'
                "  mkdir -p metamdbg_out\n"
                f'  cp "{contigs_gz}" metamdbg_out/contigs.fasta.gz\n'
                "fi\n",
            )

            res = self._run_nextflow(
                tdir,
                [
                    "--input", str(reads_dir),
                    "--outdir", str(tdir / "results"),
                    "--assembler", "metamdbg",
                    "--polish", "true",
                    "--run_remove_human", "false",
                    "--metamdbg_gfa", "false",
                ],
                extra_bins=[stub_bin],
            )
            self.assertNotEqual(res.returncode, 0)
            combined = res.stdout + "\n" + res.stderr
            self.assertIn("--read_type is required", combined)

    def test_control_metamdbg_unpolished_does_not_require_read_type(self):
        with tempfile.TemporaryDirectory(prefix="danaseq_f4_ctrl_") as tmp:
            tdir = Path(tmp)
            empty_reads = tdir / "empty_in"
            empty_reads.mkdir()
            res = self._run_nextflow(
                tdir,
                [
                    "--input", str(empty_reads),
                    "--outdir", str(tdir / "results"),
                    "--assembler", "metamdbg",
                    "--run_remove_human", "false",
                ],
            )
            combined = res.stdout + "\n" + res.stderr
            self.assertNotIn("--read_type is required", combined)
            self.assertIn("No reads found", combined)


class NormalizeAssemblyRemapGfaTests(unittest.TestCase):
    """remap_gfa prunes dropped contigs, remaps P path names, and supports hyphenated segment IDs."""

    @classmethod
    def setUpClass(cls):
        cls.mod = load_normalize_assembly()

    def test_prunes_short_contigs_and_remaps_p_path_names(self):
        with tempfile.TemporaryDirectory(prefix="danaseq_f5_") as tmp:
            tdir = Path(tmp)
            in_fa = tdir / "in.fa"
            in_fa.write_text(
                ">contig_1\n" + "A" * 500 + "\n>contig_2\n" + "C" * 1500 + "\n",
                encoding="utf-8",
            )
            in_gfa = tdir / "in.gfa"
            in_gfa.write_text(
                "H\tVN:Z:1.0\n"
                f"S\tcontig_1\t{'A' * 500}\n"
                f"S\tcontig_2\t{'C' * 1500}\n"
                "L\tcontig_1\t+\tcontig_2\t+\t0M\n"
                "P\tcontig_1\tedge_1+\t*\n"
                "P\tcontig_2\tedge_2+\t*\n",
                encoding="utf-8",
            )
            out_fa = tdir / "out.fa"
            out_info = tdir / "info.txt"
            out_map = tdir / "map.tsv"
            out_gfa = tdir / "out.gfa"

            argv = [
                "normalize_assembly.py",
                "--input", str(in_fa),
                "--output", str(out_fa),
                "--output-info", str(out_info),
                "--output-map", str(out_map),
                "--header-format", "flye",
                "--gfa-in", str(in_gfa),
                "--gfa-out", str(out_gfa),
                "--min-len", "1000",
            ]
            with mock.patch.object(sys, "argv", argv), contextlib.redirect_stderr(io.StringIO()):
                self.mod.main()

            gfa_lines = out_gfa.read_text(encoding="utf-8").splitlines()
            s_names = [ln.split("\t")[1] for ln in gfa_lines if ln.startswith("S\t")]
            l_lines = [ln for ln in gfa_lines if ln.startswith("L\t")]
            p_lines = [ln for ln in gfa_lines if ln.startswith("P\t")]
            self.assertEqual(s_names, ["contig_1"])
            self.assertEqual(l_lines, [])
            self.assertEqual(p_lines, ["P\tcontig_1\tedge_2+\t*"])

    def test_remap_gfa_hyphenated_segments_and_path_id(self):
        with tempfile.TemporaryDirectory(prefix="danaseq_f5_hyph_") as tmp:
            tdir = Path(tmp)
            hyph_in = tdir / "in.gfa"
            hyph_out = tdir / "out.gfa"
            hyph_in.write_bytes(
                b"H\tVN:Z:1.0\n"
                b"P\tu1_len-1500\tu1_len-1500+,plain_seg-\t*\n"
            )
            self.mod.remap_gfa(
                str(hyph_in),
                str(hyph_out),
                {"u1_len-1500": "contig_1", "plain_seg": "contig_2"},
            )
            self.assertEqual(
                hyph_out.read_bytes().splitlines()[1],
                b"P\tcontig_1\tcontig_1+,contig_2-\t*",
            )

    def test_control_all_kept_simple_s_and_l_lines(self):
        with tempfile.TemporaryDirectory(prefix="danaseq_f5_ctrl_") as tmp:
            tdir = Path(tmp)
            in_fa = tdir / "in.fa"
            in_fa.write_text(
                ">ctgA length=1200 circular=no coverage=10\n" + "A" * 1200 + "\n"
                ">ctgB length=2000 circular=yes coverage=30\n" + "G" * 2000 + "\n",
                encoding="utf-8",
            )
            in_gfa = tdir / "in.gfa"
            in_gfa.write_text(
                "H\tVN:Z:1.0\n"
                f"S\tctgA\t{'A' * 1200}\n"
                f"S\tctgB\t{'G' * 2000}\n"
                "L\tctgB\t+\tctgA\t-\t50M\n",
                encoding="utf-8",
            )
            out_fa = tdir / "out.fa"
            out_info = tdir / "info.txt"
            out_map = tdir / "map.tsv"
            out_gfa = tdir / "out.gfa"
            argv = [
                "normalize_assembly.py",
                "--input", str(in_fa),
                "--output", str(out_fa),
                "--output-info", str(out_info),
                "--output-map", str(out_map),
                "--header-format", "metamdbg",
                "--gfa-in", str(in_gfa),
                "--gfa-out", str(out_gfa),
                "--min-len", "1000",
            ]
            with mock.patch.object(sys, "argv", argv), contextlib.redirect_stderr(io.StringIO()):
                self.mod.main()

            gfa_lines = out_gfa.read_text(encoding="utf-8").splitlines()
            self.assertIn("L\tcontig_1\t+\tcontig_2\t-\t50M", gfa_lines)


class TetramerFreqsMaskAmbiguousTests(unittest.TestCase):
    """mask_ambiguous preserves short unambiguous contigs and terminal sequence while masking short internal islands."""

    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory(prefix="danaseq_f7_")
        cls.bin_path = Path(cls._tmp.name) / "tetramer_freqs"
        subprocess.run(
            ["gcc", "-O2", "-o", str(cls.bin_path), str(TETRAMER_FREQS_C), "-lm"],
            check=True,
        )

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def _run_calculate_tnf(self, fasta_text: str) -> tuple[dict[str, list[float]], dict[str, float]]:
        with tempfile.TemporaryDirectory(prefix="danaseq_f7_run_") as tmp:
            tdir = Path(tmp)
            in_fa = tdir / "in.fa"
            in_fa.write_text(fasta_text, encoding="utf-8")
            rendered = extract_nf_process_script(
                ASSEMBLY_NF, "CALCULATE_TNF", {"assembly": str(in_fa)}
            ).replace("tetramer_freqs", f'"{self.bin_path}"')
            subprocess.run(
                ["bash", "-ue", "-c", rendered],
                cwd=tdir,
                capture_output=True,
                text=True,
                check=True,
            )
            tnf = {}
            for ln in (tdir / "tnf.tsv").read_text(encoding="utf-8").splitlines():
                cols = ln.split("\t")
                tnf[cols[0]] = [float(x) for x in cols[1:]]
            gc = {}
            for ln in (tdir / "gc.tsv").read_text(encoding="utf-8").splitlines()[1:]:
                cid, val = ln.split("\t")
                gc[cid] = float(val)
            return tnf, gc

    def test_short_unambiguous_and_terminal_ns_preserve_valid_bases(self):
        tnf, gc = self._run_calculate_tnf(
            f">short48\n{'ACGT' * 12}\n"
            f">leftN47\n{'A' * 47}N{'C' * 80}\n"
            f">leftN48\n{'A' * 48}N{'C' * 80}\n"
            f">rightN47\n{'C' * 80}N{'A' * 47}\n"
        )
        self.assertAlmostEqual(sum(tnf["short48"]), 1.0, places=6)
        self.assertAlmostEqual(gc["short48"], 50.0, places=2)
        for cid in ("leftN47", "leftN48", "rightN47"):
            self.assertGreater(tnf[cid][0], 0.0, f"Expected non-zero AAAA/TTTT bin for {cid}")
            self.assertLess(gc[cid], 100.0, f"Expected GC < 100% when A tract is retained for {cid}")

        c_src = TETRAMER_FREQS_C.read_text(encoding="utf-8")
        self.assertNotIn("(int)span", c_src)

    def test_control_long_contig_and_short_internal_island_between_ns(self):
        tnf, gc = self._run_calculate_tnf(
            f">clean100\n{'ACGT' * 25}\n"
            f">internal_island\n{'A' * 80}N{'C' * 20}N{'A' * 80}\n"
        )
        self.assertAlmostEqual(sum(tnf["clean100"]), 1.0, places=6)
        self.assertAlmostEqual(gc["clean100"], 50.0, places=2)
        # The 20 bp C island between two Ns (span 22 < 50) must be masked, leaving 0% GC.
        self.assertAlmostEqual(sum(tnf["internal_island"]), 1.0, places=6)
        self.assertEqual(gc["internal_island"], 0.0)


if __name__ == "__main__":
    unittest.main()
