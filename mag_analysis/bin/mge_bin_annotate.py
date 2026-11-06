#!/usr/bin/env python3
"""Annotate binned contigs with geNomad/CheckV mobile-element evidence, without changing the bins.

Reads the outputs of GENOMAD_CLASSIFY and CHECKV_QUALITY (modules/mge.nf), their status files, the
assembly (for contig lengths) and the final contig-to-bin table. Writes four tables that keep
classification, grouping and host attribution apart:

  contig_class.tsv  one row per assembly contig: geNomad scores and call, topology (an observation,
                    not completeness), number of provirus regions, original bin
  mge_regions.tsv   one row per geNomad provirus region (1-based inclusive coordinates on its parent)
  mge_groups.tsv    one candidate group per geNomad plasmid/virus contig and per provirus region,
                    with CheckV quality where an exact ID match exists; never assumed complete
  host_links.tsv    host hypotheses per group; v1 has one evidence type, predicted_integration
                    (a provirus region inside an assembled contig that a binner placed in a bin)

With --alt-export it also writes bins_with_mge_separate.tsv, an alternative binning for benchmarks
(the original bins remain the default output). It needs geNomad status ok/no_hits and explicit,
uncalibrated score thresholds. "peel" moves a qualifying whole contig out of a multi-contig bin;
"rescue" (only with --rescue-unbinned) gives a qualifying unbinned contig its own group. Contigs
carrying provirus regions are never moved, so host sequence never travels with a virus.

A missing or failed geNomad run makes every classification "unavailable"; a contig absent from a
successful run is "not_scored". Neither is ever reported as chromosome. Design and review:
/data/dev/agentchat/tessera-danaseq-mge-step-design.md (Kite, Argon).
"""
import argparse
import csv
import gzip
import sys
from pathlib import Path

OK_STATES = {"ok", "no_hits"}
STATUS_VALUES = {"ok", "no_hits", "failed", "not_run"}


class InputError(Exception):
    """A hard input error: the step stops rather than write a misleading table."""


def read_status(path):
    """Return the status word from a status.tsv (tool, version, exit_code, n_input, n_output, status)."""
    if path is None or not Path(path).is_file() or Path(path).stat().st_size == 0:
        return "not_run"
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) != 1 or rows[0].get("status") not in STATUS_VALUES:
        raise InputError(f"malformed status file {path}")
    return rows[0]["status"]


def fasta_lengths(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    lengths, name = {}, None
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if name in lengths:
                    raise InputError(f"duplicate contig id in assembly: {name}")
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def read_tsv(path):
    if path is None or not Path(path).is_file() or Path(path).stat().st_size == 0:
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_contig2bin(path):
    c2b = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("@", "#")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                raise InputError(f"bad contig2bin line: {line!r}")
            if parts[0] in c2b and c2b[parts[0]] != parts[1]:
                raise InputError(f"contig {parts[0]} assigned to two bins")
            c2b[parts[0]] = parts[1]
    return c2b


def annotate(a):
    lengths = fasta_lengths(a.assembly)
    c2b = read_contig2bin(a.contig2bin)
    unknown = sorted(set(c2b) - set(lengths))
    if unknown:
        raise InputError(f"{len(unknown)} binned contigs are not in the assembly, e.g. {unknown[0]}")

    g_status = read_status(a.genomad_status)
    c_status = read_status(a.checkv_status)
    g_ok = g_status in OK_STATES

    scores, virus, plasmid, regions = {}, {}, {}, []
    if g_ok:
        for r in read_tsv(a.genomad_scores):
            scores[r["seq_name"]] = r
        for r in read_tsv(a.virus_summary):
            if "|provirus_" not in r["seq_name"]:
                virus[r["seq_name"]] = r
        for r in read_tsv(a.plasmid_summary):
            plasmid[r["seq_name"]] = r
        seen = set()
        for r in read_tsv(a.provirus):
            rid, parent = r["seq_name"], r["source_seq"]
            start, end = int(r["start"]), int(r["end"])
            if rid in seen:
                raise InputError(f"duplicate provirus region id {rid}")
            seen.add(rid)
            if parent not in lengths:
                raise InputError(f"provirus region {rid} on unknown contig {parent}")
            if not (1 <= start <= end <= lengths[parent]):
                raise InputError(f"provirus region {rid} has invalid coordinates {start}-{end} "
                                 f"on {parent} (length {lengths[parent]})")
            regions.append({"region_id": rid, "parent_contig": parent, "start": start, "end": end,
                            "source": "genomad_provirus", "score": r.get("v_vs_c_score", "")})
        for name in list(scores) + list(virus) + list(plasmid):
            if name.split("|")[0] not in lengths:
                raise InputError(f"geNomad reports {name}, which is not in the assembly")

    checkv = {}
    if c_status in OK_STATES:
        for r in read_tsv(a.checkv_quality):
            checkv[r["contig_id"]] = r

    n_regions = {}
    for reg in regions:
        n_regions[reg["parent_contig"]] = n_regions.get(reg["parent_contig"], 0) + 1

    def call_of(c):
        if not g_ok:
            return "unavailable"
        if c in plasmid and c in virus:
            return "conflicting"
        if c in plasmid:
            return "plasmid"
        if c in virus:
            return "virus"
        s = scores.get(c)
        if s is None:
            return "not_scored"
        best = max(("chromosome", "plasmid", "virus"), key=lambda k: float(s[f"{k}_score"]))
        # Scored as mobile but not kept by geNomad's own filters: leave it unresolved.
        return "chromosome" if best == "chromosome" else "ambiguous"

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    calls = {}
    with open(out / "contig_class.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["contig_id", "length", "genomad_status", "chromosome_score", "plasmid_score",
                    "virus_score", "genomad_call", "topology", "n_provirus_regions", "original_bin"])
        for c, length in lengths.items():
            s = scores.get(c, {})
            calls[c] = call_of(c)
            topo = (plasmid.get(c) or virus.get(c) or {}).get("topology", "") or "unknown"
            w.writerow([c, length, g_status, s.get("chromosome_score", "NA"), s.get("plasmid_score", "NA"),
                        s.get("virus_score", "NA"), calls[c], topo, n_regions.get(c, 0),
                        c2b.get(c, "unbinned")])

    with open(out / "mge_regions.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["region_id", "parent_contig", "start", "end", "source", "score"])
        for reg in regions:
            w.writerow([reg[k] for k in ("region_id", "parent_contig", "start", "end", "source", "score")])

    def checkv_cols(seq_id):
        if c_status not in OK_STATES:
            return [c_status, "unavailable", "unavailable", "unavailable"]
        r = checkv.get(seq_id)
        if r is None:
            return [c_status, "unknown", "unknown", "unknown"]
        return [c_status, r.get("completeness", "NA"), r.get("completeness_method", "NA"),
                r.get("checkv_quality", "NA")]

    groups = []
    for c in lengths:
        if calls[c] in ("plasmid", "virus"):
            groups.append({"group_id": f"mge_{c}", "kind": calls[c], "member": c, "parent_contig": c,
                           "start": 1, "end": lengths[c], "score": (plasmid.get(c) or virus.get(c))
                           .get(f"{calls[c]}_score", "")})
    for reg in regions:
        groups.append({"group_id": f"mge_{reg['region_id']}", "kind": "provirus", "member": reg["region_id"],
                       "parent_contig": reg["parent_contig"], "start": reg["start"], "end": reg["end"],
                       "score": reg["score"]})
    with open(out / "mge_groups.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["group_id", "kind", "member", "parent_contig", "start", "end", "genomad_score",
                    "checkv_status", "checkv_completeness", "checkv_method", "checkv_quality"])
        for g in groups:
            cv = checkv_cols(g["member"]) if g["kind"] != "plasmid" else ["n/a", "unknown", "n/a", "n/a"]
            w.writerow([g["group_id"], g["kind"], g["member"], g["parent_contig"], g["start"], g["end"],
                        g["score"]] + cv)

    with open(out / "host_links.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["group_id", "host_bin", "evidence_type", "evidence_chain"])
        for g in groups:
            parent_bin = c2b.get(g["parent_contig"])
            if g["kind"] == "provirus" and parent_bin is not None:
                w.writerow([g["group_id"], parent_bin, "predicted_integration",
                            f"geNomad provirus region {g['member']} in assembled contig {g['parent_contig']}; "
                            f"{g['parent_contig']} assigned to bin {parent_bin} by {a.binner}"])
            else:
                w.writerow([g["group_id"], "unknown", "none", ""])

    if a.alt_export:
        alt_export(a, out, lengths, c2b, calls, n_regions, plasmid, virus, g_status)
    return {"contigs": len(lengths), "groups": len(groups), "regions": len(regions),
            "genomad_status": g_status, "checkv_status": c_status}


def alt_export(a, out, lengths, c2b, calls, n_regions, plasmid, virus, g_status):
    if g_status not in OK_STATES:
        raise InputError(f"--alt-export refused: geNomad status is {g_status}")
    if a.min_plasmid_score is None or a.min_virus_score is None:
        raise InputError("--alt-export needs explicit --min-plasmid-score and --min-virus-score")
    bin_size = {}
    for c, b in c2b.items():
        bin_size[b] = bin_size.get(b, 0) + 1

    def qualifies(c):
        if n_regions.get(c, 0) > 0:
            return False
        if calls[c] == "plasmid":
            return float(plasmid[c]["plasmid_score"]) >= a.min_plasmid_score
        if calls[c] == "virus":
            return float(virus[c]["virus_score"]) >= a.min_virus_score
        return False

    with open(out / "bins_with_mge_separate.tsv", "w", newline="") as f:
        f.write(f"# alternative binning for benchmarks; thresholds plasmid>={a.min_plasmid_score} "
                f"virus>={a.min_virus_score} are UNCALIBRATED; original bins remain the default output\n")
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["contig_id", "bin", "mge_operation"])
        for c in lengths:
            b = c2b.get(c)
            q = qualifies(c)
            if b is not None and q and bin_size[b] > 1:
                w.writerow([c, f"mge_{c}", "peel"])
            elif b is not None:
                w.writerow([c, b, "kept_sole_contig" if q else "kept"])
            elif q and a.rescue_unbinned:
                w.writerow([c, f"mge_{c}", "rescue"])


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--assembly", required=True)
    p.add_argument("--contig2bin", required=True)
    p.add_argument("--binner", default="consensus")
    p.add_argument("--genomad-status")
    p.add_argument("--genomad-scores", help="aggregated classification (genomad_summary.tsv)")
    p.add_argument("--virus-summary")
    p.add_argument("--plasmid-summary")
    p.add_argument("--provirus")
    p.add_argument("--checkv-status")
    p.add_argument("--checkv-quality")
    p.add_argument("--alt-export", action="store_true")
    p.add_argument("--rescue-unbinned", action="store_true")
    p.add_argument("--min-plasmid-score", type=float)
    p.add_argument("--min-virus-score", type=float)
    p.add_argument("--outdir", default=".")
    a = p.parse_args(argv)
    try:
        summary = annotate(a)
    except InputError as e:
        print(f"[ERROR] mge_bin_annotate: {e}", file=sys.stderr)
        return 2
    print("\t".join(f"{k}={v}" for k, v in summary.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
