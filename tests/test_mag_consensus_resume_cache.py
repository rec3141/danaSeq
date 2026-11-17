"""Consensus steps rerun on -resume whenever their inputs' content changes.

MAGSCOT_CONSENSUS runs from the real module under Nextflow with stub prodigal and
hmmsearch (as in test_mag_storedir_failures). Each case resumes after one change and
checks whether the task was cached or rerun, and what mapping it published:
  unchanged inputs            -> cached
  same labels, new assignments -> rerun, new mapping
  recovered binner (3rd input) -> rerun, mapping includes it
  assembly bytes changed, size and mtime kept -> rerun
  marker database bytes changed, size and mtime kept -> rerun
  process script changed      -> rerun
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from test_mag_storedir_failures import resolve_nextflow  # noqa: E402

MODULE = HERE.parent / "mag_analysis" / "modules" / "binning.nf"
BIN = HERE.parent / "mag_analysis" / "bin"


def same_size_edit(path: Path, old: bytes, new: bytes):
    """Replace bytes in place, keeping the file's size and mtime."""
    assert len(old) == len(new)
    st = path.stat()
    data = path.read_bytes()
    assert old in data
    path.write_bytes(data.replace(old, new, 1))
    os.utime(path, ns=(st.st_atime_ns, st.st_mtime_ns))


class ConsensusResumeCacheTest(unittest.TestCase):
    def setUp(self):
        self.d = Path(tempfile.mkdtemp(prefix="consensus-cache-"))
        self.nf, extra = resolve_nextflow(self.d)
        if not self.nf:
            self.skipTest("nextflow not available")
        d = self.d
        shutil.copy(MODULE, d / "binning.nf")  # a copy, so the script-change case can edit it
        (d / "bin").symlink_to(BIN, target_is_directory=True)
        files = {
            "asm.fa": ">c1\nACTG\n>c2\nACTG\n>c3\nACTG\n",
            "one.tsv": "c1\tb1\n",
            "two.tsv": "c2\tb2\n",
            "three.tsv": "c3\tb3\n",
            "db/gtdbtk_rel207_tigrfam.hmm": "stubA\n",
            "db/gtdbtk_rel207_Pfam-A.hmm": "stubB\n",
            "db/profiles/gtdb_rel207_default_markers.tsv": "marker\tset\nM1\tA\n",
            "stubs/prodigal": '#!/bin/bash\nwhile [[ "$1" != "-a" ]]; do shift; done\n'
                              'printf ">c1_1\\nMA\\n>c2_1\\nMA\\n>c3_1\\nMA\\n" > "$2"\n',
            "stubs/hmmsearch": '#!/bin/bash\nprintf "c1_1 - M1 M1 1e-10\\nc2_1 - M1 M1 1e-10\\nc3_1 - M1 M1 1e-10\\n" > "$2"\n',
        }
        for name, text in files.items():
            p = d / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(text)
            if name.startswith("stubs/"):
                p.chmod(0o755)
        self.env = {**os.environ, "PATH": f"{d}/stubs:/opt/miniforge3/envs/checkm2/bin:{extra}:{os.environ['PATH']}",
                    "NXF_OFFLINE": "true", "NXF_HOME": str(d / ".nxf"), "NXF_ANSI_LOG": "false"}
        self.inputs = ["one", "two"]

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def run_nf(self, resume=True):
        d = self.d
        files = ",".join(f"file('{d}/{n}.tsv')" for n in self.inputs)
        labels = ",".join(f"'{n}'" for n in self.inputs)
        (d / "main.nf").write_text(f"""
params.outdir = '{d}/out'
params.store_dir = null
params.magscot_hmm_dir = '{d}/db'
include {{ MAGSCOT_CONSENSUS }} from './binning'
workflow {{
    MAGSCOT_CONSENSUS(file('{d}/asm.fa'), [{files}], [{labels}],
        file('{d}/db/gtdbtk_rel207_tigrfam.hmm'), file('{d}/db/gtdbtk_rel207_Pfam-A.hmm'),
        file('{d}/db/profiles/gtdb_rel207_default_markers.tsv'))
}}
""")
        args = [str(self.nf), "run", str(d / "main.nf"), "-work-dir", str(d / "work")]
        if resume:
            args.append("-resume")
        r = subprocess.run(args, cwd=d, env=self.env, text=True, capture_output=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        cached = "Cached process > MAGSCOT_CONSENSUS" in r.stdout
        mapping = (d / "out/binning/magscot/magscot_bins.tsv").read_text()
        return cached, mapping

    def first_run(self):
        cached, mapping = self.run_nf(resume=False)
        self.assertFalse(cached)
        return mapping

    def test_unchanged_resume_is_cached(self):
        m0 = self.first_run()
        cached, m1 = self.run_nf()
        self.assertTrue(cached)
        self.assertEqual(m0, m1)
        self.assertFalse((self.d / "store").exists(), "consensus must not write a storeDir")

    def test_changed_assignments_same_labels_rerun(self):
        self.first_run()
        (self.d / "two.tsv").write_text("c2\tb2\nc3\tb2\n")
        cached, m = self.run_nf()
        self.assertFalse(cached)
        self.assertIn("c3", m)

    def test_recovered_binner_reruns(self):
        m0 = self.first_run()
        self.assertNotIn("c3", m0)
        self.inputs = ["one", "two", "three"]
        cached, m = self.run_nf()
        self.assertFalse(cached)
        self.assertIn("c3", m)

    def test_assembly_bytes_same_size_mtime_rerun(self):
        self.first_run()
        same_size_edit(self.d / "asm.fa", b">c3\nACTG", b">c3\nACTA")
        cached, _ = self.run_nf()
        self.assertFalse(cached)

    def test_marker_db_bytes_same_size_mtime_rerun(self):
        self.first_run()
        same_size_edit(self.d / "db/gtdbtk_rel207_tigrfam.hmm", b"stubA", b"stubC")
        cached, _ = self.run_nf()
        self.assertFalse(cached)

    def test_script_change_reruns(self):
        self.first_run()
        p = self.d / "binning.nf"
        text = p.read_text()
        i = text.index("process MAGSCOT_CONSENSUS")
        j = text.index("    mkdir -p bins\n", i)
        p.write_text(text[:j] + "    mkdir -p bins  # changed\n" + text[j + len("    mkdir -p bins\n"):])
        cached, _ = self.run_nf()
        self.assertFalse(cached)


MAIN_NF = HERE.parent / "mag_analysis" / "main.nf"


def main_nf_function(name: str) -> str:
    """A top-level def from main.nf, verbatim."""
    text = MAIN_NF.read_text(encoding="utf-8")
    i = text.index(f"def {name}(")
    j = text.index("\n}\n", i) + 3
    return text[i:j]


class BinetteResumeCacheTest(unittest.TestCase):
    """BINETTE_CONSENSUS (real module, stub binette) with main.nf's fastaDirsDigest."""

    def setUp(self):
        self.d = Path(tempfile.mkdtemp(prefix="binette-cache-"))
        self.nf, extra = resolve_nextflow(self.d)
        if not self.nf:
            self.skipTest("nextflow not available")
        d = self.d
        shutil.copy(MODULE, d / "binning.nf")
        (d / "asm.fa").write_text(">c1\nACGT\n>c2\nTGCA\n")
        for b, text in {"binA/bin.1.fa": ">c1\nACGT\n", "binB/bin.2.fa": ">c2\nTGCA\n",
                        "db/checkm2.dmnd": "dmndA\n"}.items():
            (d / b).parent.mkdir(parents=True, exist_ok=True)
            (d / b).write_text(text)
        stub = d / "stubs/binette"
        stub.parent.mkdir()
        stub.write_text('#!/bin/bash\nmkdir -p binette_out/final_bins\n'
                        'while [ "$1" != --bin_dirs ]; do shift; done; shift\n'
                        'cat "$1"/*.fa > binette_out/final_bins/b1.fa\n')
        stub.chmod(0o755)
        self.env = {**os.environ, "PATH": f"{d}/stubs:{extra}:{os.environ['PATH']}",
                    "NXF_OFFLINE": "true", "NXF_HOME": str(d / ".nxf"), "NXF_ANSI_LOG": "false"}
        (d / "main.nf").write_text(f"""
params.outdir = '{d}/out'
params.store_dir = null
params.binette_min_bin_bp = 1
params.binette_max_bin_bp = 1000000000
include {{ BINETTE_CONSENSUS }} from './binning'
{main_nf_function('fastaDirsDigest')}
workflow {{
    def dirs = [file('{d}/binA'), file('{d}/binB')]
    BINETTE_CONSENSUS(file('{d}/asm.fa'), Channel.of([dirs, fastaDirsDigest(dirs)]), file('{d}/db/checkm2.dmnd'))
}}
""")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def run_nf(self, resume=True):
        args = [str(self.nf), "run", str(self.d / "main.nf"), "-work-dir", str(self.d / "work")]
        if resume:
            args.append("-resume")
        r = subprocess.run(args, cwd=self.d, env=self.env, text=True, capture_output=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return "Cached process > BINETTE_CONSENSUS" in r.stdout

    def test_unchanged_resume_is_cached(self):
        self.assertFalse(self.run_nf(resume=False))
        self.assertTrue(self.run_nf())

    def test_fasta_bytes_in_unchanged_directory_rerun(self):
        self.assertFalse(self.run_nf(resume=False))
        same_size_edit(self.d / "binA/bin.1.fa", b"ACGT", b"ACGA")
        self.assertFalse(self.run_nf())
        self.assertIn("ACGA", (self.d / "out/binning/binette/bins/binette_001.fa").read_text())

    def test_checkm2_db_bytes_same_size_mtime_rerun(self):
        self.assertFalse(self.run_nf(resume=False))
        same_size_edit(self.d / "db/checkm2.dmnd", b"dmndA", b"dmndB")
        self.assertFalse(self.run_nf())


class UnusedConsensusDatabasesTest(unittest.TestCase):
    """Full main.nf runs: --binners_only never resolves the consensus databases."""

    def setUp(self):
        from test_mag_vamb_tax_and_gtdbtk_gate import BinnersOnlyGtdbtkGateTests
        self.fx = BinnersOnlyGtdbtkGateTests("test_requested_binners_report_missing_prerequisites")
        self.fx.setUp()
        self.addCleanup(self.fx.doCleanups)

    def test_binners_only_ignores_missing_consensus_databases(self):
        missing = str(self.fx.d / "does-not-exist")
        r = self.fx._run_mag("binners_only_missing_dbs", True, self.fx.d / "store",
                             "--run_magscot", "true", "--magscot_hmm_dir", missing,
                             "--run_binette", "true", "--checkm2_db", missing)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertNotIn("MAGSCOT_CONSENSUS", r.stdout)
        self.assertNotIn("BINETTE_CONSENSUS", r.stdout)

    def test_full_run_with_missing_magscot_database_stops_with_its_path(self):
        missing = str(self.fx.d / "no-magscot-db")
        r = self.fx._run_mag("full_missing_magscot_db", False, self.fx.d / "store2",
                             "--run_magscot", "true", "--magscot_hmm_dir", missing)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("no-magscot-db", r.stdout + r.stderr)


if __name__ == "__main__":
    unittest.main()
