"""Tests for mge_bin_annotate.py (cases from the reviewed design, section 5). Run: python3 -m pytest -q"""
import csv
import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("mba", Path(__file__).with_name("mge_bin_annotate.py"))
mba = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mba)

SEQS = {"host1": 10000, "host2": 5000, "plas1": 3000, "amb1": 2000, "free": 4000}


def write(p, text):
    p.write_text(text)
    return str(p)


def status(p, s):
    return write(p, f"tool\tversion\texit_code\tn_input\tn_output\tstatus\ngenomad\t1.13.0\t0\t5\t3\t{s}\n")


def tsv(p, header, rows):
    return write(p, "\t".join(header) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))


@pytest.fixture
def toy(tmp_path):
    asm = write(tmp_path / "asm.fa", "".join(f">{c}\n{'A' * n}\n" for c, n in SEQS.items()))
    c2b = write(tmp_path / "c2b.tsv", "host1\tbinA\nhost2\tbinA\nplas1\tbinA\namb1\tbinB\n")
    files = dict(
        assembly=asm, contig2bin=c2b, binner="dastool",
        genomad_status=status(tmp_path / "g.status", "ok"),
        genomad_scores=tsv(tmp_path / "scores.tsv", ["seq_name", "chromosome_score", "plasmid_score", "virus_score"],
                           [["host1", .9, .05, .05], ["host2", .95, .03, .02], ["plas1", .02, .97, .01],
                            ["amb1", .3, .4, .3], ["free", .9, .05, .05]]),
        virus_summary=tsv(tmp_path / "v.tsv", ["seq_name", "length", "topology", "virus_score"],
                          [["host1|provirus_101_2100", 2000, "Provirus", .95]]),
        plasmid_summary=tsv(tmp_path / "p.tsv", ["seq_name", "length", "topology", "plasmid_score"],
                            [["plas1", 3000, "DTR", .97]]),
        provirus=tsv(tmp_path / "pv.tsv", ["seq_name", "source_seq", "start", "end", "v_vs_c_score"],
                     [["host1|provirus_101_2100", "host1", 101, 2100, 40.0]]),
        checkv_status=None, checkv_quality=None, alt_export=False, rescue_unbinned=False,
        min_plasmid_score=None, min_virus_score=None, outdir=str(tmp_path / "out"))
    return tmp_path, files


def run(files, **over):
    args = type("A", (), {**files, **over})()
    return mba.annotate(args)


def rows(out, name):
    with open(Path(out) / name) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def test_toy_case(toy):
    _, f = toy
    run(f)
    cls = {r["contig_id"]: r for r in rows(f["outdir"], "contig_class.tsv")}
    assert set(cls) == set(SEQS)                       # every contig exactly once
    assert cls["plas1"]["genomad_call"] == "plasmid"
    assert cls["amb1"]["genomad_call"] == "ambiguous"
    assert cls["host1"]["n_provirus_regions"] == "1"
    assert cls["plas1"]["topology"] == "DTR"           # observation only
    groups = rows(f["outdir"], "mge_groups.tsv")
    assert {g["member"] for g in groups} == {"plas1", "host1|provirus_101_2100"}
    assert all(g["checkv_completeness"] in ("unavailable", "unknown") for g in groups)
    links = {l["group_id"]: l for l in rows(f["outdir"], "host_links.tsv")}
    assert links["mge_host1|provirus_101_2100"]["evidence_type"] == "predicted_integration"
    assert links["mge_host1|provirus_101_2100"]["host_bin"] == "binA"
    assert links["mge_plas1"]["host_bin"] == "unknown"


def test_genomad_failed_is_unavailable_and_blocks_export(toy):
    tmp, f = toy
    f["genomad_status"] = status(tmp / "gf.status", "failed")
    run(f)
    assert {r["genomad_call"] for r in rows(f["outdir"], "contig_class.tsv")} == {"unavailable"}
    with pytest.raises(mba.InputError):
        run(f, alt_export=True, min_plasmid_score=.9, min_virus_score=.9)


def test_genomad_missing_status_is_not_run(toy):
    _, f = toy
    f["genomad_status"] = None
    run(f)
    assert {r["genomad_status"] for r in rows(f["outdir"], "contig_class.tsv")} == {"not_run"}
    assert {r["genomad_call"] for r in rows(f["outdir"], "contig_class.tsv")} == {"unavailable"}


def test_no_hits(toy):
    tmp, f = toy
    f.update(genomad_status=status(tmp / "gn.status", "no_hits"),
             genomad_scores=tsv(tmp / "ns.tsv", ["seq_name", "chromosome_score", "plasmid_score", "virus_score"], []),
             virus_summary=tsv(tmp / "nv.tsv", ["seq_name", "length", "topology", "virus_score"], []),
             plasmid_summary=tsv(tmp / "np.tsv", ["seq_name", "length", "topology", "plasmid_score"], []),
             provirus=tsv(tmp / "npv.tsv", ["seq_name", "source_seq", "start", "end", "v_vs_c_score"], []))
    run(f, alt_export=True, min_plasmid_score=.9, min_virus_score=.9)
    assert {r["genomad_call"] for r in rows(f["outdir"], "contig_class.tsv")} == {"not_scored"}
    assert rows(f["outdir"], "mge_groups.tsv") == []
    assert all(r["mge_operation"] == "kept" for r in rows(f["outdir"], "bins_with_mge_separate.tsv"))


def test_two_proviruses_one_contig(toy):
    tmp, f = toy
    f["provirus"] = tsv(tmp / "pv2.tsv", ["seq_name", "source_seq", "start", "end", "v_vs_c_score"],
                        [["host1|provirus_101_2100", "host1", 101, 2100, 40],
                         ["host1|provirus_5001_9000", "host1", 5001, 9000, 30]])
    run(f, alt_export=True, min_plasmid_score=.5, min_virus_score=.5)
    assert len(rows(f["outdir"], "mge_regions.tsv")) == 2
    assert sum(g["kind"] == "provirus" for g in rows(f["outdir"], "mge_groups.tsv")) == 2
    alt = {r["contig_id"]: r for r in rows(f["outdir"], "bins_with_mge_separate.tsv")}
    assert alt["host1"]["bin"] == "binA"               # region-bearing contig never moves


def test_contig_absent_from_genomad_is_not_scored(toy):
    tmp, f = toy
    f["genomad_scores"] = tsv(tmp / "s2.tsv", ["seq_name", "chromosome_score", "plasmid_score", "virus_score"],
                              [["host1", .9, .05, .05]])
    run(f)
    cls = {r["contig_id"]: r["genomad_call"] for r in rows(f["outdir"], "contig_class.tsv")}
    assert cls["free"] == "not_scored" and cls["host2"] == "not_scored"


def test_checkv_absent_leaves_classification(toy):
    tmp, f = toy
    f["checkv_status"] = status(tmp / "c.status", "failed")
    run(f)
    cls = {r["contig_id"]: r["genomad_call"] for r in rows(f["outdir"], "contig_class.tsv")}
    assert cls["plas1"] == "plasmid"
    prov = [g for g in rows(f["outdir"], "mge_groups.tsv") if g["kind"] == "provirus"][0]
    assert prov["checkv_status"] == "failed" and prov["checkv_completeness"] == "unavailable"


@pytest.mark.parametrize("start,end", [(0, 10), (50, 40), (9000, 10001)])
def test_invalid_coordinates(toy, start, end):
    tmp, f = toy
    f["provirus"] = tsv(tmp / "bad.tsv", ["seq_name", "source_seq", "start", "end", "v_vs_c_score"],
                        [[f"host1|provirus_{start}_{end}", "host1", start, end, 1]])
    with pytest.raises(mba.InputError):
        run(f)


def test_duplicate_region_ids(toy):
    tmp, f = toy
    f["provirus"] = tsv(tmp / "dup.tsv", ["seq_name", "source_seq", "start", "end", "v_vs_c_score"],
                        [["host1|provirus_101_2100", "host1", 101, 2100, 1]] * 2)
    with pytest.raises(mba.InputError):
        run(f)


def test_alt_export_needs_thresholds(toy):
    _, f = toy
    with pytest.raises(mba.InputError):
        run(f, alt_export=True)


def test_alt_export_peel_rescue_and_accounting(toy):
    tmp, f = toy
    f["plasmid_summary"] = tsv(tmp / "p2.tsv", ["seq_name", "length", "topology", "plasmid_score"],
                               [["plas1", 3000, "DTR", .97], ["free", 4000, "No terminal repeats", .95]])
    run(f, alt_export=True, min_plasmid_score=.9, min_virus_score=.9)
    alt = {r["contig_id"]: r for r in rows(f["outdir"], "bins_with_mge_separate.tsv")}
    assert alt["plas1"]["mge_operation"] == "peel" and alt["plas1"]["bin"] == "mge_plas1"
    assert "free" not in alt                          # unbinned, no rescue
    assert set(alt) == {"host1", "host2", "plas1", "amb1"}   # same contigs as the original binning
    run(f, alt_export=True, rescue_unbinned=True, min_plasmid_score=.9, min_virus_score=.9)
    alt = {r["contig_id"]: r for r in rows(f["outdir"], "bins_with_mge_separate.tsv")}
    assert alt["free"]["mge_operation"] == "rescue"


def test_sole_contig_bin_not_duplicated(toy):
    tmp, f = toy
    f["contig2bin"] = write(tmp / "c2b2.tsv", "host1\tbinA\nplas1\tbinP\n")
    run(f, alt_export=True, min_plasmid_score=.9, min_virus_score=.9)
    alt = {r["contig_id"]: r for r in rows(f["outdir"], "bins_with_mge_separate.tsv")}
    assert alt["plas1"]["bin"] == "binP" and alt["plas1"]["mge_operation"] == "kept_sole_contig"


def test_binned_contig_missing_from_assembly(toy):
    tmp, f = toy
    f["contig2bin"] = write(tmp / "c2b3.tsv", "ghost\tbinZ\n")
    with pytest.raises(mba.InputError):
        run(f)


def test_checkv_join_success(toy):
    tmp, f = toy
    f["checkv_status"] = status(tmp / "cok.status", "ok")
    f["checkv_quality"] = tsv(tmp / "cq.tsv", ["contig_id", "completeness", "completeness_method", "checkv_quality"],
                              [["host1|provirus_101_2100", 87.5, "AAI-based", "High-quality"]])
    run(f)
    prov = [g for g in rows(f["outdir"], "mge_groups.tsv") if g["kind"] == "provirus"][0]
    assert prov["checkv_completeness"] == "87.5" and prov["checkv_method"] == "AAI-based"


@pytest.mark.parametrize("key", ["genomad_scores", "virus_summary", "plasmid_summary", "provirus"])
def test_missing_required_genomad_file_is_error(toy, key):
    tmp, f = toy
    f[key] = write(tmp / f"empty_{key}.tsv", "")
    with pytest.raises(mba.InputError):
        run(f)


def test_provirus_summary_without_region_record_is_error(toy):
    # Kite's reproduction: whole-host virus call plus a provirus summary row, coordinates file emptied.
    tmp, f = toy
    f["virus_summary"] = tsv(tmp / "v3.tsv", ["seq_name", "length", "topology", "virus_score"],
                             [["host1", 10000, "No terminal repeats", .95],
                              ["host1|provirus_101_2100", 2000, "Provirus", .95]])
    f["provirus"] = tsv(tmp / "pv3.tsv", ["seq_name", "source_seq", "start", "end", "v_vs_c_score"], [])
    with pytest.raises(mba.InputError):
        run(f, alt_export=True, min_plasmid_score=.5, min_virus_score=.5)


def test_new_bin_ids_avoid_existing_bins(toy):
    tmp, f = toy
    f["contig2bin"] = write(tmp / "c2b4.tsv", "host1\tmge_plas1\nhost2\tmge_plas1\nplas1\tmge_plas1\n")
    run(f, alt_export=True, min_plasmid_score=.9, min_virus_score=.9)
    alt = {r["contig_id"]: r for r in rows(f["outdir"], "bins_with_mge_separate.tsv")}
    assert alt["plas1"]["mge_operation"] == "peel"
    assert alt["plas1"]["bin"] != "mge_plas1" and alt["host1"]["bin"] == "mge_plas1"


@pytest.mark.parametrize("bad", [-0.1, 1.5, float("nan"), float("inf")])
def test_threshold_must_be_unit_interval(toy, bad):
    _, f = toy
    with pytest.raises(mba.InputError):
        run(f, alt_export=True, min_plasmid_score=bad, min_virus_score=.9)


@pytest.mark.parametrize("bad", ["nan", "-1", "2", "x"])
def test_scores_must_be_unit_interval(toy, bad):
    tmp, f = toy
    f["genomad_scores"] = tsv(tmp / "bs.tsv", ["seq_name", "chromosome_score", "plasmid_score", "virus_score"],
                              [["host1", bad, .05, .05]])
    with pytest.raises(mba.InputError):
        run(f)


@pytest.mark.parametrize("key,header,row", [
    ("genomad_scores", ["seq_name", "chromosome_score", "plasmid_score", "virus_score"], ["host1", .9, .05, .05]),
    ("plasmid_summary", ["seq_name", "length", "topology", "plasmid_score"], ["plas1", 3000, "DTR", .97]),
])
def test_duplicate_evidence_ids_are_error(toy, key, header, row):
    tmp, f = toy
    f[key] = tsv(tmp / f"dup_{key}.tsv", header, [row, row])
    with pytest.raises(mba.InputError):
        run(f)


def test_exact_ids_not_prefix_split(toy):
    tmp, f = toy
    f["plasmid_summary"] = tsv(tmp / "p4.tsv", ["seq_name", "length", "topology", "plasmid_score"],
                               [["plas1|extra", 3000, "DTR", .97]])
    with pytest.raises(mba.InputError):
        run(f)
