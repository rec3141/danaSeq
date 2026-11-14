#!/usr/bin/env python3
"""Regression and control tests for nanopore_live, illumina_rna, illumina_amplicon, and illumina_assembly fixes."""

import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import textwrap
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
MINIFORGE_PYTHON = Path("/opt/miniforge3/bin/python3")
PYTHON_BIN = str(MINIFORGE_PYTHON if MINIFORGE_PYTHON.exists() else Path(sys.executable))
NXF_PKG = Path("/opt/miniforge3/pkgs/nextflow-25.10.4-h2a3209d_0")
JAVA_PKGS = sorted(Path("/opt/miniforge3/pkgs").glob("openjdk-*/bin"))


def _install_offline_nextflow(dest_bin: Path) -> Path:
    src = (NXF_PKG / "bin" / "nextflow").read_text()
    dist = NXF_PKG / "share" / "nextflow" / "dist"
    lines = [
        f"NXF_DIST={dist}" if ln.startswith("NXF_DIST=") else ln
        for ln in src.splitlines()
    ]
    dest_bin.write_text("\n".join(lines) + "\n")
    dest_bin.chmod(0o755)
    return dest_bin


def _nxf_env(nxf_home: Path, extra_path: str = "") -> dict:
    nxf_home.mkdir(parents=True, exist_ok=True)
    java_dir = str(JAVA_PKGS[-1])
    java_home = str(JAVA_PKGS[-1].parent)
    parts = [p for p in [extra_path, java_dir, "/opt/miniforge3/bin", os.environ.get("PATH", "")] if p]
    env = dict(os.environ)
    env["PATH"] = ":".join(parts)
    env["JAVA_HOME"] = java_home
    env["NXF_HOME"] = str(nxf_home)
    env["NXF_OFFLINE"] = "true"
    env["NXF_ANSI_LOG"] = "false"
    env["NXF_DISABLE_CHECK_LATEST"] = "true"
    return env


def _install_container_shims(fakebin: Path, captured_json: Path) -> None:
    fakebin.mkdir(parents=True, exist_ok=True)
    for runtime in ("docker", "apptainer"):
        shim = fakebin / runtime
        shim.write_text(
            f"#!{PYTHON_BIN}\n"
            "import json, os, sys\n"
            "args = sys.argv[1:]\n"
            "if args[:1] == ['run']:\n"
            f"    open({str(captured_json)!r}, 'w').write(json.dumps({{'runtime': {runtime!r}, 'argv': args}}))\n"
            "    if os.access('.', os.W_OK):\n"
            "        open('.nextflow.log', 'w').write('Session UUID: 00000000-0000-0000-0000-000000000000\\n')\n"
            "sys.exit(0)\n"
        )
        shim.chmod(0o755)


class LiveAndIlluminaHighsTest(unittest.TestCase):
    def test_f4_container_launchers_bind_external_paths_and_store_dir(self):
        """Container launchers mount external metadata, samplesheet, DBs, and writable store_dir."""
        amp_sh = REPO_ROOT / "illumina_amplicon" / "run-illumina-amplicon.sh"
        rna_sh = REPO_ROOT / "illumina_rna" / "run-illumina-rna.sh"
        live_sh = REPO_ROOT / "nanopore_live" / "run-realtime.sh"

        with tempfile.TemporaryDirectory(prefix="f4-reg-") as tmp:
            root = Path(tmp)
            fakebin = root / "fakebin"
            captured = root / "captured.json"
            _install_container_shims(fakebin, captured)
            env = dict(os.environ, PATH=f"{fakebin}:/usr/bin:/bin")

            inp = root / "reads_in"
            out = root / "out_dir"
            ext_meta_dir = root / "ext_meta"
            ext_sheet_dir = root / "ext_sheet"
            ext_sheet_reads = root / "sheet_reads"
            ext_primers_dir = root / "ext_primers"
            ext_refdb_dir = root / "ext_refdb"
            ext_earth_dir = root / "ext_earthdna"
            ext_bakta = root / "ext_bakta"
            ext_map_refs = root / "ext_mapping_refs"
            ext_dmnd_dir = root / "ext_dmnd"
            for d in (
                inp, out, ext_meta_dir, ext_sheet_dir, ext_sheet_reads,
                ext_primers_dir, ext_refdb_dir, ext_earth_dir,
                ext_bakta, ext_map_refs, ext_dmnd_dir,
            ):
                d.mkdir(parents=True)

            r1_ext = ext_sheet_reads / "S1_R1.fastq.gz"
            r2_ext = ext_sheet_reads / "S1_R2.fastq.gz"
            r1_ext.write_bytes(b"r1")
            r2_ext.write_bytes(b"r2")

            amp_sheet = ext_sheet_dir / "amplicon_sheet.csv"
            amp_sheet.write_text(
                "sample_id,run,plate_id,r1_path,r2_path,primer_pair\n"
                f"S1,run1,P1,{r1_ext},{r2_ext},EMP-16S\n"
            )
            rna_sheet = ext_sheet_dir / "rna_sheet.csv"
            rna_sheet.write_text(f"sample,r1,r2\nS1,{r1_ext},{r2_ext}\n")

            metadata = ext_meta_dir / "metadata.tsv"
            metadata.write_text("sample_name\tsite\nS1\tA\n")
            pfwd = ext_primers_dir / "fwd.fa"
            prev = ext_primers_dir / "rev.fa"
            pfwd.write_text(">f\nACGT\n")
            prev.write_text(">r\nTGCA\n")
            silva_db = ext_refdb_dir / "silva.fasta"
            silva_db.write_text(">ref1\nACGT\n")
            dmnd_file = ext_dmnd_dir / "HADEG.dmnd"
            dmnd_genes = ext_dmnd_dir / "HADEG.genes.tsv"
            dmnd_file.write_bytes(b"dmnd")
            dmnd_genes.write_text("target\tgene\nt1\tg1\n")

            # New store directories that do not exist prior to launcher invocation
            amp_store = root / "new_amp_store"
            live_store = root / "new_live_store"

            # 1. illumina_amplicon Docker mode
            res_amp_docker = subprocess.run(
                [
                    "bash", str(amp_sh),
                    "--runtime", "docker",
                    "--input", str(inp),
                    "--outdir", str(out / "amp_docker"),
                    "--samplesheet", str(amp_sheet),
                    "--metadata", str(metadata),
                    "--primers_fwd", str(pfwd),
                    "--primers_rev", str(prev),
                    "--ref_databases", f"silva:{silva_db}:Domain,Phylum",
                    "--earthdna_db", str(ext_earth_dir),
                    "--store_dir", str(amp_store),
                ],
                cwd=str(root),
                env=env,
                capture_output=True,
                text=True,
                timeout=15,
            )
            self.assertEqual(res_amp_docker.returncode, 0, res_amp_docker.stderr)
            self.assertTrue(amp_store.is_dir(), "Launcher must create --store_dir on host before docker run")
            amp_argv = json.loads(captured.read_text())["argv"]
            amp_mounts = [amp_argv[i + 1] for i, a in enumerate(amp_argv[:-1]) if a == "-v"]
            for expected_dir in (
                str(ext_sheet_dir), str(ext_sheet_reads), str(ext_meta_dir),
                str(ext_primers_dir), str(ext_refdb_dir), str(ext_earth_dir),
            ):
                self.assertTrue(
                    any(expected_dir in m for m in amp_mounts),
                    f"Missing bind mount covering {expected_dir} in {amp_mounts}",
                )
            self.assertIn(f"{amp_store}:{amp_store}", amp_mounts, "store_dir must be bound read-write")

            # 2. illumina_amplicon Apptainer mode
            sif_file = root / "test_amp.sif"
            sif_file.write_bytes(b"sif")
            res_amp_app = subprocess.run(
                [
                    "bash", str(amp_sh),
                    "--runtime=apptainer",
                    f"--sif={sif_file}",
                    f"--outdir={out / 'amp_app'}",
                    f"--samplesheet={amp_sheet}",
                    f"--metadata={metadata}",
                    f"--store_dir={amp_store}",
                ],
                cwd=str(root),
                env=env,
                capture_output=True,
                text=True,
                timeout=15,
            )
            self.assertEqual(res_amp_app.returncode, 0, res_amp_app.stderr)
            app_argv = json.loads(captured.read_text())["argv"]
            app_binds = [app_argv[i + 1] for i, a in enumerate(app_argv[:-1]) if a == "--bind"]
            app_bind_blob = ",".join(app_binds)
            for expected in (str(ext_sheet_dir), str(ext_sheet_reads), str(ext_meta_dir), f"{amp_store}:{amp_store}"):
                self.assertIn(expected, app_bind_blob)

            # 3. illumina_rna Docker mode (--metadata isolated from other paths)
            refs_dir = root / "rna_refs"
            refs_dir.mkdir()
            res_rna = subprocess.run(
                [
                    "bash", str(rna_sh),
                    "--docker",
                    "--samplesheet", str(rna_sheet),
                    "--references", str(refs_dir),
                    "--metadata", str(metadata),
                    "--outdir", str(out / "rna"),
                ],
                cwd=str(root),
                env=env,
                capture_output=True,
                text=True,
                timeout=15,
            )
            self.assertEqual(res_rna.returncode, 0, res_rna.stderr)
            rna_argv = json.loads(captured.read_text())["argv"]
            rna_mounts = [rna_argv[i + 1] for i, a in enumerate(rna_argv[:-1]) if a == "-v"]
            self.assertTrue(
                any(str(ext_meta_dir) in m or str(metadata) in m for m in rna_mounts),
                f"--metadata not bound in illumina_rna mounts: {rna_mounts}",
            )

            # 4. nanopore_live Docker mode (--bakta_db, --mapping_refs, --protein_dbs, --metadata, --store_dir)
            res_live = subprocess.run(
                [
                    "bash", str(live_sh),
                    "--docker",
                    "--input", str(inp),
                    "--outdir", str(out / "live"),
                    "--bakta_db", str(ext_bakta),
                    "--mapping_refs", str(ext_map_refs),
                    "--protein_dbs", str(dmnd_file),
                    "--metadata", str(metadata),
                    "--store_dir", str(live_store),
                ],
                cwd=str(root),
                env=env,
                capture_output=True,
                text=True,
                timeout=15,
            )
            self.assertEqual(res_live.returncode, 0, res_live.stderr)
            self.assertTrue(live_store.is_dir(), "run-realtime.sh must create --store_dir before docker run")
            live_argv = json.loads(captured.read_text())["argv"]
            live_mounts = [live_argv[i + 1] for i, a in enumerate(live_argv[:-1]) if a == "-v"]
            for expected_ro in (str(ext_bakta), str(ext_map_refs), str(dmnd_file), str(dmnd_genes), str(metadata)):
                self.assertTrue(
                    any(expected_ro in m and m.endswith(":ro") for m in live_mounts),
                    f"Missing read-only mount for {expected_ro}: {live_mounts}",
                )
            self.assertIn(f"{live_store}:{live_store}", live_mounts)

    def test_f4_control_basic_input_outdir_mounts_preserved(self):
        """Control: basic --input and --outdir container mounts continue to work across all three launchers."""
        amp_sh = REPO_ROOT / "illumina_amplicon" / "run-illumina-amplicon.sh"
        rna_sh = REPO_ROOT / "illumina_rna" / "run-illumina-rna.sh"
        live_sh = REPO_ROOT / "nanopore_live" / "run-realtime.sh"

        with tempfile.TemporaryDirectory(prefix="f4-ctrl-") as tmp:
            root = Path(tmp)
            fakebin = root / "fakebin"
            captured = root / "captured.json"
            _install_container_shims(fakebin, captured)
            env = dict(os.environ, PATH=f"{fakebin}:/usr/bin:/bin")

            inp = root / "inp"
            refs = root / "refs"
            kdb = root / "kdb"
            for d in (inp, refs, kdb):
                d.mkdir()
            (kdb / "hash.k2d").write_bytes(b"k2d")

            res_amp = subprocess.run(
                ["bash", str(amp_sh), "--runtime", "docker", "--input", str(inp), "--outdir", str(root / "out_amp")],
                cwd=str(root), env=env, capture_output=True, text=True, timeout=15,
            )
            self.assertEqual(res_amp.returncode, 0, res_amp.stderr)
            amp_mounts = [
                a for i, a in enumerate(json.loads(captured.read_text())["argv"])
                if i > 0 and json.loads(captured.read_text())["argv"][i - 1] == "-v"
            ]
            self.assertIn(f"{inp}:{inp}:ro", amp_mounts)
            self.assertIn(f"{root / 'out_amp'}:{root / 'out_amp'}", amp_mounts)

            res_rna = subprocess.run(
                ["bash", str(rna_sh), "--docker", "--input", str(inp), "--references", str(refs), "--outdir", str(root / "out_rna")],
                cwd=str(root), env=env, capture_output=True, text=True, timeout=15,
            )
            self.assertEqual(res_rna.returncode, 0, res_rna.stderr)
            rna_mounts = [
                a for i, a in enumerate(json.loads(captured.read_text())["argv"])
                if i > 0 and json.loads(captured.read_text())["argv"][i - 1] == "-v"
            ]
            self.assertTrue(any(str(inp) in m for m in rna_mounts))
            self.assertTrue(any(str(refs) in m for m in rna_mounts))

            res_live = subprocess.run(
                ["bash", str(live_sh), "--docker", "--input", str(inp), "--outdir", str(root / "out_live"), "--kraken_db", str(kdb)],
                cwd=str(root), env=env, capture_output=True, text=True, timeout=15,
            )
            self.assertEqual(res_live.returncode, 0, res_live.stderr)
            live_mounts = [
                a for i, a in enumerate(json.loads(captured.read_text())["argv"])
                if i > 0 and json.loads(captured.read_text())["argv"][i - 1] == "-v"
            ]
            self.assertIn(f"{inp}:/data/input:ro", live_mounts)
            self.assertIn(f"{kdb}:/kraken_db:ro", live_mounts)

    def _setup_sortmerna_stubs(self, fakebin: Path, capture_file: Path) -> None:
        fakebin.mkdir(parents=True, exist_ok=True)
        reformat = fakebin / "reformat.sh"
        reformat.write_text(
            textwrap.dedent(
                """\
                #!/usr/bin/env bash
                for a in "$@"; do
                    case "$a" in
                        out1=*|out2=*|out=*) printf "@r1\\nACGT\\n+\\nIIII\\n" | gzip -c > "${a#*=}" ;;
                    esac
                done
                echo "Input: 2 reads" >&2
                """
            )
        )
        reformat.chmod(0o755)

        sortmerna = fakebin / "sortmerna"
        sortmerna.write_text(
            f"#!{PYTHON_BIN}\n"
            "import json, os, sys\n"
            "args = sys.argv[1:]\n"
            "refs = [args[i + 1] for i, a in enumerate(args[:-1]) if a == '--ref']\n"
            "contents = [open(r).read() if os.path.isfile(r) else None for r in refs]\n"
            "rec = {'refs': refs, 'refs_exist': [os.path.isfile(r) for r in refs], 'contents': contents}\n"
            f"open({str(capture_file)!r}, 'w').write(json.dumps(rec))\n"
            "open('other_fwd.fq.gz', 'wb').write(b'x')\n"
            "open('other_rev.fq.gz', 'wb').write(b'x')\n"
            "open('aligned.log', 'w').write('aligned\\n')\n"
        )
        sortmerna.chmod(0o755)

    def test_f7_remove_rrna_sortmerna_stages_refs_and_handles_fna_and_gz(self):
        """REMOVE_RRNA_SORTMERNA stages relative --rrna_refs, expands .fna/.gz, converts U->T, and fails on empty dir."""
        main_nf = (REPO_ROOT / "illumina_rna" / "main.nf").read_text()
        rrna_nf_path = REPO_ROOT / "illumina_rna" / "modules" / "rrna.nf"
        rrna_nf = rrna_nf_path.read_text()
        self.assertRegex(
            main_nf,
            r"REMOVE_RRNA_SORTMERNA\s*\(\s*ch_reads\s*,\s*ch_sortmerna_refs\s*\)",
            "main.nf must pass staged reference channel into REMOVE_RRNA_SORTMERNA",
        )
        self.assertIn("path(refs)", rrna_nf.split("process REMOVE_RRNA_SORTMERNA", 1)[1].split("output:", 1)[0])

        with tempfile.TemporaryDirectory(prefix="f7-reg-") as tmp:
            root = Path(tmp)
            nf_bin = _install_offline_nextflow(root / "nextflow")
            fakebin = root / "fakebin"
            capture = root / "sortmerna_capture.json"
            self._setup_sortmerna_stubs(fakebin, capture)

            reads = root / "S1.fq.gz"
            with gzip.open(reads, "wt") as gz:
                gz.write("@r1\nACGT\n+\nIIII\n")

            # Relative directory with .fna, .fasta.gz (RNA U alphabet), and .fna.gz
            rel_dir = root / "rel_silva_refs"
            rel_dir.mkdir()
            (rel_dir / "1_ssu.fna").write_text(">ssu1\nACGUGCAU\n")
            with gzip.open(rel_dir / "2_lsu.fasta.gz", "wt") as gz:
                gz.write(">lsu1\nuuacgg\n")
            with gzip.open(rel_dir / "3_5s.fna.gz", "wt") as gz:
                gz.write(">5s1\nGGCCAA\n")

            cfg = root / "nf.config"
            cfg.write_text("conda.enabled = false\n")
            wf = root / "run_sortmerna.nf"
            wf.write_text(
                textwrap.dedent(
                    f"""\
                    params.outdir = "{root / 'out'}"
                    params.store_dir = null
                    params.read_len = 150
                    params.rrna_refs = null
                    include {{ REMOVE_RRNA_SORTMERNA }} from '{rrna_nf_path}'
                    workflow {{
                        ch_reads = Channel.of([[id: 'S1'], file("{reads}")])
                        ch_sortmerna_refs = Channel.value(file(params.rrna_refs, checkIfExists: true))
                        REMOVE_RRNA_SORTMERNA(ch_reads, ch_sortmerna_refs)
                    }}
                    """
                )
            )

            res = subprocess.run(
                [str(nf_bin), "-c", str(cfg), "run", str(wf), "--rrna_refs", "rel_silva_refs"],
                cwd=str(root),
                env=_nxf_env(root / "nxf_home_1", str(fakebin)),
                capture_output=True,
                text=True,
                timeout=40,
            )
            self.assertEqual(res.returncode, 0, res.stdout + "\n" + res.stderr)
            rec = json.loads(capture.read_text())
            self.assertEqual(len(rec["refs"]), 3)
            self.assertEqual(rec["refs_exist"], [True, True, True])
            joined_contents = "\n".join(rec["contents"])
            self.assertIn("ACGTGCAT", joined_contents)
            self.assertIn("TTACGG", joined_contents)
            self.assertIn("GGCCAA", joined_contents)

            # Empty directory must fail non-zero rather than invoking sortmerna with no --ref
            empty_dir = root / "empty_refs"
            empty_dir.mkdir()
            capture.unlink()
            res_empty = subprocess.run(
                [str(nf_bin), "-c", str(cfg), "run", str(wf), "--rrna_refs", str(empty_dir)],
                cwd=str(root),
                env=_nxf_env(root / "nxf_home_2", str(fakebin)),
                capture_output=True,
                text=True,
                timeout=40,
            )
            self.assertNotEqual(res_empty.returncode, 0)
            self.assertFalse(capture.exists(), "sortmerna must not run when zero reference FASTAs exist")

    def test_f7_control_absolute_fa_and_fasta_refs_passed_to_sortmerna(self):
        """Control: absolute directory containing .fa and .fasta references succeeds."""
        rrna_nf_path = REPO_ROOT / "illumina_rna" / "modules" / "rrna.nf"
        rrna_block = rrna_nf_path.read_text().split("process REMOVE_RRNA_SORTMERNA", 1)[1].split("output:", 1)[0]
        takes_staged_refs = "path(refs)" in rrna_block

        with tempfile.TemporaryDirectory(prefix="f7-ctrl-") as tmp:
            root = Path(tmp)
            nf_bin = _install_offline_nextflow(root / "nextflow")
            fakebin = root / "fakebin"
            capture = root / "sortmerna_capture.json"
            self._setup_sortmerna_stubs(fakebin, capture)

            reads = root / "S1.fq.gz"
            with gzip.open(reads, "wt") as gz:
                gz.write("@r1\nACGT\n+\nIIII\n")

            abs_refs = root / "abs_silva"
            abs_refs.mkdir()
            (abs_refs / "a.fa").write_text(">a\nACGTACGT\n")
            (abs_refs / "b.fasta").write_text(">b\nTGCATGCA\n")

            call_expr = (
                f'REMOVE_RRNA_SORTMERNA(ch_reads, Channel.value(file("{abs_refs}", checkIfExists: true)))'
                if takes_staged_refs
                else "REMOVE_RRNA_SORTMERNA(ch_reads)"
            )
            cfg = root / "nf.config"
            cfg.write_text("conda.enabled = false\n")
            wf = root / "ctrl_sortmerna.nf"
            wf.write_text(
                textwrap.dedent(
                    f"""\
                    params.outdir = "{root / 'out'}"
                    params.store_dir = null
                    params.read_len = 150
                    params.rrna_refs = "{abs_refs}"
                    include {{ REMOVE_RRNA_SORTMERNA }} from '{rrna_nf_path}'
                    workflow {{
                        ch_reads = Channel.of([[id: 'S1'], file("{reads}")])
                        {call_expr}
                    }}
                    """
                )
            )
            res = subprocess.run(
                [str(nf_bin), "-c", str(cfg), "run", str(wf)],
                cwd=str(root),
                env=_nxf_env(root / "nxf_home_ctrl", str(fakebin)),
                capture_output=True,
                text=True,
                timeout=40,
            )
            self.assertEqual(res.returncode, 0, res.stdout + "\n" + res.stderr)
            rec = json.loads(capture.read_text())
            self.assertEqual(len(rec["refs"]), 2)
            self.assertEqual(rec["refs_exist"], [True, True])
            self.assertTrue((root / "out" / "preprocess" / "S1" / "S1.norrna.fq.gz").is_file())

    def test_f8_augment_viz_functions_reads_gzipped_gene_counts_from_summarize_expression(self):
        """augment_viz_functions.py loads .gene_counts.tsv.gz produced by summarize_expression.py."""
        summarize_py = REPO_ROOT / "illumina_rna" / "bin" / "summarize_expression.py"
        augment_py = REPO_ROOT / "illumina_rna" / "bin" / "augment_viz_functions.py"

        with tempfile.TemporaryDirectory(prefix="f8-reg-") as tmp:
            root = Path(tmp)
            ref = "MAG_01"
            for sub in (
                "genes", "contigs", "summary", "flagstat", "readcounts", "viz",
                "mag/binning/checkm2", "mag/annotation/bakta/basic", "mag/binning/dastool",
            ):
                (root / sub).mkdir(parents=True, exist_ok=True)
            expr_ref_dir = root / "expression" / ref
            expr_ref_dir.mkdir(parents=True)

            (root / "mag/annotation/bakta/basic/annotation.tsv").write_text(
                "#Sequence Id\tType\tStart\tStop\tStrand\tLocus Tag\tGene\tProduct\tDbXrefs\n"
                + "".join(
                    f"c1\tCDS\t{i*100+1}\t{(i+1)*100}\t+\tg{i}\tgene{i}\tprod_{i}\tKEGG:K0000{i}\n"
                    for i in range(6)
                )
            )

            gene_annot = root / "genes.tsv.gz"
            with gzip.open(gene_annot, "wt") as gz:
                gz.write("gene_id\tcontig\tstart\tend\tstrand\tlength\ttype\tgene\tproduct\tdbxref\n")
                for i in range(6):
                    gz.write(f"g{i}\tc1\t{i*100+1}\t{(i+1)*100}\t+\t100\tCDS\tgene{i}\tprod_{i}\tKEGG:K0000{i}\n")

            contig_len = root / "contig_lengths.tsv"
            contig_len.write_text("contig\tlength\nc1\t1000\n")

            for s_idx, s_name in enumerate(("S1", "S2"), start=1):
                with gzip.open(root / "genes" / f"{s_name}_vs_{ref}.genes.tsv.gz", "wt") as gz:
                    gz.write("gene_id\tcount\n")
                    for i in range(6):
                        gz.write(f"g{i}\t{(i + 1) * 10 * s_idx}\n")
                with gzip.open(root / "contigs" / f"{s_name}_vs_{ref}.contigs.tsv.gz", "wt") as gz:
                    gz.write("contig\tcount\nc1\t100\n")
                (root / "summary" / f"{s_name}_vs_{ref}.summary").write_text("Status\tsample\nAssigned\t100\n")
                (root / "flagstat" / f"{s_name}_vs_{ref}.flagstat").write_text(
                    "200 + 0 primary mapped (100%)\n200 + 0 properly paired (100%)\n"
                )
                (root / "readcounts" / f"{s_name}.counts.tsv").write_text(f"{s_name}\traw\t300\n{s_name}\tfiltered\t250\n")

            res_sum = subprocess.run(
                [
                    PYTHON_BIN, "-B", str(summarize_py),
                    "--ref", ref,
                    "--genes-dir", str(root / "genes"),
                    "--contigs-dir", str(root / "contigs"),
                    "--summary-dir", str(root / "summary"),
                    "--flagstat-dir", str(root / "flagstat"),
                    "--readcounts-dir", str(root / "readcounts"),
                    "--gene-annot", str(gene_annot),
                    "--contig-lengths", str(contig_len),
                    "--out", str(expr_ref_dir),
                ],
                capture_output=True,
                text=True,
                timeout=30,
            )
            self.assertEqual(res_sum.returncode, 0, res_sum.stderr)
            self.assertTrue((expr_ref_dir / f"{ref}.gene_counts.tsv.gz").is_file())
            self.assertFalse((expr_ref_dir / f"{ref}.gene_counts.tsv").exists())

            (root / "viz" / "references.json").write_text(json.dumps([{"name": ref}]))
            res_aug = subprocess.run(
                [
                    PYTHON_BIN, "-B", str(augment_py),
                    "--viz-dir", str(root / "viz"),
                    "--mag-dir", str(root / "mag"),
                    "--expression-dir", str(root / "expression"),
                ],
                capture_output=True,
                text=True,
                timeout=30,
            )
            self.assertEqual(res_aug.returncode, 0, res_aug.stderr)

            with gzip.open(root / "viz" / "functions.json.gz", "rt") as gz:
                functions = json.load(gz)
            with gzip.open(root / "viz" / "tsne.json.gz", "rt") as gz:
                tsne_out = json.load(gz)

            self.assertIn(ref, functions)
            self.assertEqual(functions[ref]["samples"], ["S1", "S2"])
            self.assertEqual(functions[ref]["genes"], [f"g{i}" for i in range(6)])
            self.assertEqual(functions[ref]["annotations"][0]["product"], "prod_0")
            self.assertIn(ref, tsne_out)
            self.assertEqual(len(tsne_out[ref]["coords"]), 6)

    def test_f8_control_uncompressed_gene_counts_tsv_still_supported(self):
        """Control: augment_viz_functions.py continues to read uncompressed .gene_counts.tsv files."""
        augment_py = REPO_ROOT / "illumina_rna" / "bin" / "augment_viz_functions.py"

        with tempfile.TemporaryDirectory(prefix="f8-ctrl-") as tmp:
            root = Path(tmp)
            ref = "MAG_PLAIN"
            for sub in ("viz", f"expression/{ref}", "mag/binning/checkm2", "mag/annotation/bakta/basic", "mag/binning/dastool"):
                (root / sub).mkdir(parents=True, exist_ok=True)

            counts_tsv = root / "expression" / ref / f"{ref}.gene_counts.tsv"
            counts_tsv.write_text(
                "gene_id\tS1\tS2\n"
                + "".join(f"g{i}\t{(i+1)*5}\t{(i+1)*15}\n" for i in range(6))
            )
            (root / "viz" / "references.json").write_text(json.dumps([{"name": ref}]))

            res = subprocess.run(
                [
                    PYTHON_BIN, "-B", str(augment_py),
                    "--viz-dir", str(root / "viz"),
                    "--mag-dir", str(root / "mag"),
                    "--expression-dir", str(root / "expression"),
                ],
                capture_output=True,
                text=True,
                timeout=30,
            )
            self.assertEqual(res.returncode, 0, res.stderr)
            with gzip.open(root / "viz" / "functions.json.gz", "rt") as gz:
                functions = json.load(gz)
            self.assertIn(ref, functions)
            self.assertEqual(functions[ref]["genes"], [f"g{i}" for i in range(6)])

    def _extract_groovy_functions(self, nf_text: str) -> str:
        blocks = []
        m_derive = re.search(r"(def deriveSampleId\([^\)]*\)\s*\{.*?^\})", nf_text, flags=re.S | re.M)
        if m_derive:
            blocks.append(m_derive.group(1))
        m_disc = re.search(r"(def discoverReads\(\)\s*\{.*?^\})", nf_text, flags=re.S | re.M)
        if m_disc:
            blocks.append(m_disc.group(1))
        return "\n\n".join(blocks)

    def _build_asm_harness(self, inp_dir: Path) -> str:
        asm_nf = (REPO_ROOT / "illumina_assembly" / "main.nf").read_text()
        funcs = self._extract_groovy_functions(asm_nf)
        if "def discoverReads()" in funcs:
            return (
                f'params.input = "{inp_dir}"\n\n'
                + funcs
                + "\n\nworkflow {\n"
                + "    ch_raw_reads = Channel.fromList(discoverReads())\n"
                + '    ch_raw_reads.map { meta, r1, r2 -> "${meta.id}<-${r1.name}:${r2.name}" }\n'
                + '        .toSortedList().view { "ASM_IDS:" + it.join(",") }\n'
                + "}\n"
            )
        m_chan = re.search(r"(ch_raw_reads\s*=\s*Channel\.fromFilePairs\(.*?^\s*\})", asm_nf, flags=re.S | re.M)
        assert m_chan is not None
        return (
            f'params.input = "{inp_dir}"\n\nworkflow {{\n    '
            + m_chan.group(1)
            + '\n    ch_raw_reads.map { meta, r1, r2 -> "${meta.id}<-${r1.name}:${r2.name}" }\n'
            + '        .toSortedList().view { "ASM_IDS:" + it.join(",") }\n}\n'
        )

    def _build_rna_harness(self, inp_dir: Path) -> str:
        rna_nf = (REPO_ROOT / "illumina_rna" / "main.nf").read_text()
        funcs = self._extract_groovy_functions(rna_nf)
        return (
            f'params.samplesheet = null\nparams.input = "{inp_dir}"\n\n'
            + funcs
            + "\n\nworkflow {\n"
            + "    def s = discoverReads()\n"
            + '    println("RNA_IDS:" + s.collect { "${it[0]}<-${it[1].name}:${it[2].name}" }.join(","))\n'
            + "}\n"
        )

    def test_f9_multi_underscore_sample_ids_preserved_and_true_duplicates_fail(self):
        """Multi-underscore sample IDs stay distinct in illumina_assembly and illumina_rna; true lane collisions fail."""
        with tempfile.TemporaryDirectory(prefix="f9-reg-") as tmp:
            root = Path(tmp)
            nf_bin = _install_offline_nextflow(root / "nextflow")
            inp = root / "reads"
            inp.mkdir()
            for fn in (
                "siteA_rep1_S1_L001_R1_001.fastq.gz",
                "siteA_rep1_S1_L001_R2_001.fastq.gz",
                "siteA_rep2_S2_L001_R1_001.fastq.gz",
                "siteA_rep2_S2_L001_R2_001.fastq.gz",
                "cohort_R1_sub3_R1_001.fastq.gz",
                "cohort_R1_sub3_R2_001.fastq.gz",
            ):
                with gzip.open(inp / fn, "wt") as gz:
                    gz.write(f"@{fn}\nACGT\n+\nIIII\n")

            asm_wf = root / "asm.nf"
            asm_wf.write_text(self._build_asm_harness(inp))
            res_asm = subprocess.run(
                [str(nf_bin), "run", str(asm_wf)],
                cwd=str(root),
                env=_nxf_env(root / "nxf_asm"),
                capture_output=True,
                text=True,
                timeout=35,
            )
            self.assertEqual(res_asm.returncode, 0, res_asm.stdout + "\n" + res_asm.stderr)
            self.assertIn(
                "ASM_IDS:"
                "cohort_R1_sub3<-cohort_R1_sub3_R1_001.fastq.gz:cohort_R1_sub3_R2_001.fastq.gz,"
                "siteA_rep1<-siteA_rep1_S1_L001_R1_001.fastq.gz:siteA_rep1_S1_L001_R2_001.fastq.gz,"
                "siteA_rep2<-siteA_rep2_S2_L001_R1_001.fastq.gz:siteA_rep2_S2_L001_R2_001.fastq.gz",
                res_asm.stdout,
            )

            rna_wf = root / "rna.nf"
            rna_wf.write_text(self._build_rna_harness(inp))
            res_rna = subprocess.run(
                [str(nf_bin), "run", str(rna_wf)],
                cwd=str(root),
                env=_nxf_env(root / "nxf_rna"),
                capture_output=True,
                text=True,
                timeout=35,
            )
            self.assertEqual(res_rna.returncode, 0, res_rna.stdout + "\n" + res_rna.stderr)
            self.assertIn(
                "RNA_IDS:"
                "cohort_R1_sub3<-cohort_R1_sub3_R1_001.fastq.gz:cohort_R1_sub3_R2_001.fastq.gz,"
                "siteA_rep1<-siteA_rep1_S1_L001_R1_001.fastq.gz:siteA_rep1_S1_L001_R2_001.fastq.gz,"
                "siteA_rep2<-siteA_rep2_S2_L001_R1_001.fastq.gz:siteA_rep2_S2_L001_R2_001.fastq.gz",
                res_rna.stdout,
            )

            # Unconcatenated multi-lane files for the same sample must fail loudly in illumina_assembly
            dup_inp = root / "dup_reads"
            dup_inp.mkdir()
            for fn in (
                "siteA_rep1_S1_L001_R1_001.fastq.gz",
                "siteA_rep1_S1_L001_R2_001.fastq.gz",
                "siteA_rep1_S1_L002_R1_001.fastq.gz",
                "siteA_rep1_S1_L002_R2_001.fastq.gz",
            ):
                with gzip.open(dup_inp / fn, "wt") as gz:
                    gz.write(f"@{fn}\nACGT\n+\nIIII\n")
            asm_dup_wf = root / "asm_dup.nf"
            asm_dup_wf.write_text(self._build_asm_harness(dup_inp))
            res_dup = subprocess.run(
                [str(nf_bin), "run", str(asm_dup_wf)],
                cwd=str(root),
                env=_nxf_env(root / "nxf_asm_dup"),
                capture_output=True,
                text=True,
                timeout=35,
            )
            self.assertNotEqual(res_dup.returncode, 0)
            self.assertIn("Duplicate sample IDs: siteA_rep1", res_dup.stdout + "\n" + res_dup.stderr)

    def test_f9_control_single_token_illumina_and_bgi_ids_unchanged(self):
        """Control: single-token Illumina (_S1_L001_R1_001) and BGI (_1.fq.gz) names derive identical IDs."""
        with tempfile.TemporaryDirectory(prefix="f9-ctrl-") as tmp:
            root = Path(tmp)
            nf_bin = _install_offline_nextflow(root / "nextflow")
            inp = root / "reads"
            inp.mkdir()
            for fn in (
                "sample1_S1_L001_R1_001.fastq.gz",
                "sample1_S1_L001_R2_001.fastq.gz",
                "sample2_S2_L001_R1_001.fastq.gz",
                "sample2_S2_L001_R2_001.fastq.gz",
            ):
                with gzip.open(inp / fn, "wt") as gz:
                    gz.write(f"@{fn}\nACGT\n+\nIIII\n")

            asm_wf = root / "asm_ctrl.nf"
            asm_wf.write_text(self._build_asm_harness(inp))
            res_asm = subprocess.run(
                [str(nf_bin), "run", str(asm_wf)],
                cwd=str(root),
                env=_nxf_env(root / "nxf_asm_ctrl"),
                capture_output=True,
                text=True,
                timeout=35,
            )
            self.assertEqual(res_asm.returncode, 0, res_asm.stdout + "\n" + res_asm.stderr)
            self.assertIn(
                "ASM_IDS:"
                "sample1<-sample1_S1_L001_R1_001.fastq.gz:sample1_S1_L001_R2_001.fastq.gz,"
                "sample2<-sample2_S2_L001_R1_001.fastq.gz:sample2_S2_L001_R2_001.fastq.gz",
                res_asm.stdout,
            )

            # Add a BGI pair for illumina_rna control
            for fn in ("bgiSample_1.fq.gz", "bgiSample_2.fq.gz"):
                with gzip.open(inp / fn, "wt") as gz:
                    gz.write(f"@{fn}\nACGT\n+\nIIII\n")
            rna_wf = root / "rna_ctrl.nf"
            rna_wf.write_text(self._build_rna_harness(inp))
            res_rna = subprocess.run(
                [str(nf_bin), "run", str(rna_wf)],
                cwd=str(root),
                env=_nxf_env(root / "nxf_rna_ctrl"),
                capture_output=True,
                text=True,
                timeout=35,
            )
            self.assertEqual(res_rna.returncode, 0, res_rna.stdout + "\n" + res_rna.stderr)
            self.assertIn(
                "RNA_IDS:"
                "bgiSample<-bgiSample_1.fq.gz:bgiSample_2.fq.gz,"
                "sample1<-sample1_S1_L001_R1_001.fastq.gz:sample1_S1_L001_R2_001.fastq.gz,"
                "sample2<-sample2_S2_L001_R1_001.fastq.gz:sample2_S2_L001_R2_001.fastq.gz",
                res_rna.stdout,
            )


if __name__ == "__main__":
    unittest.main()
