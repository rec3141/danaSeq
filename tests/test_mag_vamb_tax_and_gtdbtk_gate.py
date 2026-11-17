#!/usr/bin/env python3
"""Regression and control tests for BIN_VAMB_TAX taxonomy formatting and --binners_only GTDB-Tk gating."""

import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from typing import Union
import unittest

REPO = Path(__file__).resolve().parents[1]
MAG_DIR = REPO / "mag_analysis"
MAIN_NF = MAG_DIR / "main.nf"
BINNING_NF = MAG_DIR / "modules" / "binning.nf"
TAXONOMY_NF = MAG_DIR / "modules" / "taxonomy.nf"

LOCAL_VAMB_TAX = Path("/opt/miniforge3/envs/vamb/lib/python3.11/site-packages/vamb/taxonomy.py")
NXF_TEMPLATE = Path("/opt/miniforge3/pkgs/nextflow-25.10.4-h2a3209d_0/bin/nextflow")
NXF_DIST = Path("/opt/miniforge3/pkgs/nextflow-25.10.4-h2a3209d_0/share/nextflow/dist")


def write_executable(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


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


def extract_nf_process_script(nf_path: Path, process_name: str, subs: dict[str, str]) -> str:
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


def extract_sendsketch_parser() -> str:
    rendered = extract_nf_process_script(TAXONOMY_NF, "SENDSKETCH_CLASSIFY", {})
    m = re.search(r"python3 - sendsketch_raw\.txt <<'PYEOF'\n(.*?)\nPYEOF", rendered, re.DOTALL)
    if not m:
        raise ValueError("Could not locate SENDSKETCH_CLASSIFY Python parser")
    return m.group(1)


def extract_vamb_tax_converter(sendsketch_tsv: Path, assembly_fa: Path) -> str:
    rendered = extract_nf_process_script(
        BINNING_NF,
        "BIN_VAMB_TAX",
        {"sendsketch_contigs": str(sendsketch_tsv), "assembly": str(assembly_fa)},
    )
    m = re.search(r'python3 -c "\n(.*?)\n"', rendered, re.DOTALL)
    if not m:
        raise ValueError("Could not locate BIN_VAMB_TAX taxonomy converter")
    return m.group(1)


class BinVambTaxConversionTests(unittest.TestCase):
    """BIN_VAMB_TAX maps unclassified contigs to empty root and truncates lineages at missing ranks."""

    def test_unclassified_maps_to_empty_and_missing_rank_truncates_unambiguously(self):
        with tempfile.TemporaryDirectory(prefix="danaseq_c02_f6_") as tmp:
            d = Path(tmp)
            raw_sketch = d / "sendsketch_raw.txt"
            cols_c2 = "\t".join(
                [
                    "90.0%", "85.0%", "95.00%", "0", "0", "0", "10", "10",
                    "1", "1000", "1", "Bac1", "d:Bacteria;p:Bacillota;o:Bacillales",
                ]
            )
            cols_c3 = "\t".join(
                [
                    "90.0%", "85.0%", "96.00%", "0", "0", "0", "10", "10",
                    "2", "1000", "1", "Bac2", "d:Bacteria;p:Bacillota;c:Bacilli;o:Bacillales",
                ]
            )
            raw_sketch.write_text(
                f"Query: c1\tSketchLen: 100\nNo hits.\n"
                f"Query: c2\tSketchLen: 100\nWKID\tKID\tANI\n{cols_c2}\n"
                f"Query: c3\tSketchLen: 100\nWKID\tKID\tANI\n{cols_c3}\n",
                encoding="utf-8",
            )
            r_sketch = subprocess.run(
                [sys.executable, "-B", "-", str(raw_sketch)],
                input=extract_sendsketch_parser(),
                cwd=d,
                capture_output=True,
                text=True,
            )
            self.assertEqual(r_sketch.returncode, 0, r_sketch.stderr)

            assembly = d / "assembly.fa"
            assembly.write_text(">c1\nACGT\n>c2\nACGT\n>c3\nACGT\n", encoding="utf-8")
            conv_py = extract_vamb_tax_converter(d / "sendsketch_contigs.tsv", assembly)
            r_vtax = subprocess.run(
                [sys.executable, "-B", "-c", conv_py],
                cwd=d,
                capture_output=True,
                text=True,
            )
            self.assertEqual(r_vtax.returncode, 0, r_vtax.stderr)
            self.assertIn("2/3 contigs classified", r_vtax.stderr)

            lines = (d / "vamb_taxonomy.tsv").read_text(encoding="utf-8").splitlines()
            self.assertEqual(
                lines,
                [
                    "contigs\tpredictions",
                    "c1\t",
                    "c2\tBacteria;Bacillota",
                    "c3\tBacteria;Bacillota;Bacilli;Bacillales",
                ],
            )

            # Verify VAMB's rank-tree invariant directly (no taxon appears at multiple depths).
            if LOCAL_VAMB_TAX.is_file():
                tax_lines = LOCAL_VAMB_TAX.read_text(encoding="utf-8").splitlines()
                idx = next(i for i, ln in enumerate(tax_lines) if ln.startswith("def assert_unambiguous_ranks"))
                ns: dict[str, object] = {"Union": Union}

                class DummyContigTax:
                    def __init__(self, ranks: list[str]) -> None:
                        self.ranks = ranks

                class DummyTax:
                    def __init__(self, items: list[DummyContigTax]) -> None:
                        self.contig_taxonomies = items

                ns["ContigTaxonomy"] = DummyContigTax
                ns["Taxonomy"] = DummyTax
                ns["PredictedTaxonomy"] = DummyTax
                exec("\n".join(tax_lines[idx:]), ns)
                parsed = []
                for ln in lines[1:]:
                    pred = ln.split("\t", 1)[1]
                    parsed.append(DummyContigTax([] if pred == "" else pred.split(";")))
                ns["assert_unambiguous_ranks"](DummyTax(parsed))

    def test_zero_classified_contigs_skips_vamb_invocation_cleanly(self):
        with tempfile.TemporaryDirectory(prefix="danaseq_c02_f6_allu_") as tmp:
            d = Path(tmp)
            sendsketch_tsv = d / "sendsketch_contigs.tsv"
            sendsketch_tsv.write_text(
                "contig_id\tstatus\tANI\tref_name\tlineage\n"
                "c1\tU\t0\tUnclassified\tUnclassified\n",
                encoding="utf-8",
            )
            assembly = d / "assembly.fa"
            assembly.write_text(">c1\n" + "ACGT" * 600 + "\n", encoding="utf-8")
            depths = d / "depths.txt"
            depths.write_text(
                "contigName\tcontigLen\ttotalAvgDepth\ts1.bam\ts1.bam-var\n"
                "c1\t2400\t10.0\t10.0\t1.0\n",
                encoding="utf-8",
            )
            stub_bin = d / "bin"
            stub_bin.mkdir()
            write_executable(
                stub_bin / "vamb",
                "#!/usr/bin/env bash\necho 'vamb should not be called when 0 contigs are classified' >&2\nexit 99\n",
            )
            write_executable(stub_bin / "torch_device.sh", "#!/usr/bin/env bash\necho cpu\n")
            rendered = extract_nf_process_script(
                BINNING_NF,
                "BIN_VAMB_TAX",
                {
                    "sendsketch_contigs": str(sendsketch_tsv),
                    "assembly": str(assembly),
                    "depths": str(depths),
                    "params.gpu": "false",
                    "task.cpus": "1",
                },
            )
            env = dict(os.environ, PATH=f"{stub_bin}:{os.environ.get('PATH', '')}")
            res = subprocess.run(["bash", "-ue", "-c", rendered], cwd=d, env=env, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, res.stderr)
            self.assertTrue((d / "vamb_tax_bins.tsv").is_file())
            self.assertEqual((d / "vamb_tax_bins.tsv").read_text(encoding="utf-8"), "")

    def test_control_complete_seven_rank_lineage_preserved(self):
        with tempfile.TemporaryDirectory(prefix="danaseq_c02_f6_ctrl_") as tmp:
            d = Path(tmp)
            sendsketch_tsv = d / "sendsketch_contigs.tsv"
            sendsketch_tsv.write_text(
                "contig_id\tstatus\tANI\tref_name\tlineage\n"
                "c1\tC\t98.50\tRef\td__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus subtilis\n",
                encoding="utf-8",
            )
            assembly = d / "assembly.fa"
            assembly.write_text(">c1\nACGT\n", encoding="utf-8")
            conv_py = extract_vamb_tax_converter(sendsketch_tsv, assembly)
            res = subprocess.run([sys.executable, "-B", "-c", conv_py], cwd=d, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, res.stderr)
            lines = (d / "vamb_taxonomy.tsv").read_text(encoding="utf-8").splitlines()
            self.assertEqual(
                lines,
                [
                    "contigs\tpredictions",
                    "c1\tBacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus subtilis",
                ],
            )


class BinnersOnlyGtdbtkGateTests(unittest.TestCase):
    """--binners_only skips GTDBTK_CLASSIFY so storeDir stays clean for the full consensus run."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="danaseq_c02_f7_")
        self.d = Path(self.tmp.name)
        self.bin_dir = self.d / "bin"
        self.bin_dir.mkdir()
        self.nxf_bin = resolve_nextflow(self.bin_dir)
        if self.nxf_bin is None:
            self.skipTest("Nextflow is not available in this environment")

        write_executable(
            self.bin_dir / "metabat2",
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            "out=''\n"
            "while [[ $# -gt 0 ]]; do\n"
            '  case "$1" in -o) out="$2"; shift 2 ;; *) shift ;; esac\n'
            "done\n"
            'mkdir -p "$(dirname "$out")"\n'
            'printf \'>c1\\nACGTACGTACGT\\n>c2\\nACGTACGTACGT\\n\' > "${out}.1.fa"\n',
        )
        write_executable(
            self.bin_dir / "DAS_Tool",
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            "out=''\n"
            "while [[ $# -gt 0 ]]; do\n"
            '  case "$1" in -o) out="$2"; shift 2 ;; *) shift ;; esac\n'
            "done\n"
            'mkdir -p "${out}_DASTool_bins"\n'
            'printf \'>c1\\nACGTACGTACGT\\n\' > "${out}_DASTool_bins/metabat_001_sub.fa"\n'
            'printf \'c1\\tmetabat_001_sub\\n\' > "${out}_DASTool_contig2bin.tsv"\n'
            'printf \'bin\\tSCG_completeness\\tSCG_redundancy\\nmetabat_001_sub\\t100\\t0\\n\' > "${out}_DASTool_summary.tsv"\n'
            'printf \'bin\\tSCG_completeness\\nmetabat_001\\t100\\n\' > "${out}_allBins.eval"\n',
        )
        self.gtdbtk_calls_log = self.d / "gtdbtk_invocations.log"
        write_executable(
            self.bin_dir / "gtdbtk",
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            f"echo invoked >> '{self.gtdbtk_calls_log}'\n"
            "gdir=''; odir=''\n"
            "while [[ $# -gt 0 ]]; do\n"
            '  case "$1" in --genome_dir) gdir="$2"; shift 2 ;; --out_dir) odir="$2"; shift 2 ;; *) shift ;; esac\n'
            "done\n"
            'mkdir -p "$odir"\n'
            "hdr='user_genome\\tclassification\\tclosest_genome_reference\\tclosest_genome_reference_radius\\tclosest_genome_taxonomy\\tclosest_genome_ani\\tclosest_genome_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings'\n"
            'printf "%b\\n" "$hdr" > "$odir/gtdbtk.bac120.summary.tsv"\n'
            'for f in "$gdir"/*.fa; do\n'
            '  [ -e "$f" ] || continue\n'
            '  stem="$(basename "$f" .fa)"\n'
            '  printf \'%s\\td__Bacteria;p__Bacillota\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\t11\\tNA\\tNA\\n\' "$stem" >> "$odir/gtdbtk.bac120.summary.tsv"\n'
            "done\n",
        )

        self.asm = self.d / "assembly.fa"
        self.asm.write_text(">c1\nACGTACGTACGT\n>c2\nACGTACGTACGT\n", encoding="utf-8")
        self.depths = self.d / "depths.txt"
        self.depths.write_text(
            "contigName\tcontigLen\ttotalAvgDepth\ts1.bam\ts1.bam-var\n"
            "c1\t12\t10.0\t10.0\t1.0\n"
            "c2\t12\t10.0\t10.0\t1.0\n",
            encoding="utf-8",
        )
        self.gtdb_db = self.d / "gtdb_db"
        self.gtdb_db.mkdir()
        self.override_cfg = self.d / "override.config"
        self.override_cfg.write_text(
            "conda.enabled = false\n"
            "timeline.enabled = false\n"
            "report.enabled = false\n"
            "trace.enabled = false\n"
            "dag.enabled = false\n"
            "params.gtdbtk_min_memory = '0 B'\n"
            "process {\n"
            "    cpus = 1\n"
            "    memory = 512.MB\n"
            "    maxRetries = 0\n"
            "    withLabel: 'process_low' { cpus = 1; memory = 512.MB }\n"
            "    withLabel: 'process_medium' { cpus = 1; memory = 512.MB }\n"
            "    withLabel: 'process_high' { cpus = 1; memory = 512.MB }\n"
            "    withName: 'GTDBTK_CLASSIFY' { cpus = 1; memory = 512.MB }\n"
            "}\n",
            encoding="utf-8",
        )

    def tearDown(self):
        self.tmp.cleanup()

    def _run_mag(self, name: str, binners_only: bool, store_dir: Path, extra=()) -> subprocess.CompletedProcess:
        rcwd = self.d / name
        rcwd.mkdir()
        env = dict(
            os.environ,
            PATH=f"{self.bin_dir}:{os.environ.get('PATH', '')}",
            NXF_HOME=str(self.d / "nxf_home"),
            NXF_OFFLINE="true",
            NXF_ANSI_LOG="false",
            NXF_JVM_ARGS="-Xms64m -Xmx512m",
        )
        return subprocess.run(
            [
                str(self.nxf_bin),
                "-log", str(rcwd / "nextflow.log"),
                "run", str(MAIN_NF),
                "-c", str(self.override_cfg),
                "-w", str(rcwd / "work"),
                "--assembly", str(self.asm),
                "--depths", str(self.depths),
                "--outdir", str(rcwd / "out"),
                "--store_dir", str(store_dir),
                "--annotator", "none",
                "--run_ecossdb", "false",
                "--run_viz", "false",
                "--gpu", "false",
                "--binners_only", "true" if binners_only else "false",
                "--run_gtdbtk", "true",
                "--gtdbtk_db", str(self.gtdb_db),
                *extra,
            ],
            cwd=rcwd,
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )

    def test_binners_only_skips_gtdbtk_so_consensus_run_classifies_dastool_bins(self):
        store_dir = self.d / "shared_store"
        r1 = self._run_mag("run1_binners", binners_only=True, store_dir=store_dir)
        self.assertEqual(r1.returncode, 0, r1.stdout + "\n" + r1.stderr)
        stored_tsv = store_dir / "taxonomy" / "gtdbtk" / "gtdbtk_taxonomy.tsv"
        self.assertFalse(stored_tsv.exists())
        self.assertFalse(self.gtdbtk_calls_log.exists())

        r2 = self._run_mag("run2_full", binners_only=False, store_dir=store_dir)
        self.assertEqual(r2.returncode, 0, r2.stdout + "\n" + r2.stderr)
        self.assertTrue(stored_tsv.is_file())
        genomes = {
            ln.split("\t")[0]
            for ln in stored_tsv.read_text(encoding="utf-8").splitlines()[1:]
            if ln.strip()
        }
        self.assertIn("metabat_001", genomes)
        self.assertIn("dastool-metabat_001_sub", genomes)
        self.assertEqual(len(self.gtdbtk_calls_log.read_text(encoding="utf-8").splitlines()), 1)

    def test_control_direct_full_run_classifies_consensus_and_raw_bins(self):
        store_dir = self.d / "fresh_store"
        r = self._run_mag("run_ctrl", binners_only=False, store_dir=store_dir)
        self.assertEqual(r.returncode, 0, r.stdout + "\n" + r.stderr)
        stored_tsv = store_dir / "taxonomy" / "gtdbtk" / "gtdbtk_taxonomy.tsv"
        genomes = {
            ln.split("\t")[0]
            for ln in stored_tsv.read_text(encoding="utf-8").splitlines()[1:]
            if ln.strip()
        }
        self.assertEqual(genomes, {"metabat_001", "dastool-metabat_001_sub"})

    def test_too_little_memory_fails_gtdbtk_fast_and_run_completes(self):
        store_dir = self.d / "lowmem_store"
        r = self._run_mag("run_lowmem", binners_only=False, store_dir=store_dir,
                          extra=("--gtdbtk_min_memory", "1 GB"))
        self.assertEqual(r.returncode, 0, r.stdout + "\n" + r.stderr)
        self.assertFalse((store_dir / "taxonomy" / "gtdbtk" / "gtdbtk_taxonomy.tsv").exists())
        self.assertFalse(self.gtdbtk_calls_log.exists())          # gtdbtk never started
        log = (self.d / "run_lowmem" / "nextflow.log").read_text(encoding="utf-8")
        self.assertIn("GTDBTK_CLASSIFY (gtdbtk)` terminated with an error exit status (1) -- Error is ignored", log)
        errs = [p.read_text(encoding="utf-8")
                for p in (self.d / "run_lowmem" / "work").glob("*/*/.command.err")]
        self.assertTrue(any("GTDB-Tk needs 1 GB of memory" in e for e in errs), errs)


if __name__ == "__main__":
    unittest.main()
