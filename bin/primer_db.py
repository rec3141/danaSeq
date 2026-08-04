"""Curated primer database and read-based primer detection.

Ported from OMC's portal/app/primers.py so the pipeline can answer, from the
reads themselves, which primers a sample carries — and say which gene those
primers target. Keeping the table here is what lets samples.json report
`assay_gene`/`assay_lineage` directly, instead of emitting a bare primer name
that only a downstream service holding the table could interpret.

What deliberately did NOT come across: fetching reads from SRA, parsing primers
out of submission metadata, and the multi-run orchestration over accessions.
Those are a submission system's job, not a pipeline's.

Detection is two-tier:
  (a) match sampled 5' ends against the curated database;
  (b) failing that, build a degenerate IUPAC consensus of the conserved 5'
      prefix de novo.

Both are read sampling plus string matching over a few hundred reads — cheap
compared with DETECT_PRIMERS, which ran a full cutadapt pass per candidate
primer file per sample purely to count survivors.
"""
from __future__ import annotations

import csv
import gzip
import logging
import math
import os
import re
from collections import Counter

# IUPAC nucleotide codes → the set of bases each represents.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
# Reverse map: frozenset of bases → the tightest IUPAC code.
_IUPAC_REV = {frozenset(v): k for k, v in IUPAC.items()}

# Curated primer database — common 16S/18S/ITS amplicon primers (5'->3'),
# mirroring microscape's bundled sets plus a few widely used pairs. Sequences
# use IUPAC degeneracy; the reverse primer is written 5'->3' as synthesised.
# Curated core of named primer pairs, hand-verified against real reads. This is
# the canonical layer: 18S / protist primers (which FoodMicrobionet does not
# cover) plus the standard 16S/ITS pairs with clean names. The vendored
# FoodMicrobionet tables are merged on top for breadth (see _load_vendored_
# primers), deduped by sequence so these canonical entries win name resolution.
#
# Sequences are the biological primer only (5'->3', IUPAC, adapters stripped).
# Detection matches on SEQUENCE, never on the name in SRA metadata — EMP renamed
# 515FB->515F(Parada)/806R(Apprill), so submitter names are unreliable (exactly
# what mislabelled PRJNA1473294's 18S runs as "16S"). Sources: Herlemann 2011;
# Parada 2016 / Apprill 2015 / EMP; Caporaso 2011; Quince 2011; Lane 1991;
# Stoeck 2010; Amaral-Zettler 2009; Comeau 2011; White 1990; Gardes & Bruns
# 1993; Ihrmark 2012; UNITE; pr2-primers (Vaulot 2022).
_CORE_PRIMER_DB = [
    # ── 16S rRNA (bacteria / archaea) ──
    {"name": "341F", "rev_name": "805R", "region": "16S V3-V4",
     "fwd": "CCTACGGGNGGCWGCAG", "rev": "GACTACHVGGGTATCTAATCC"},
    {"name": "515F", "rev_name": "806R", "region": "16S V4",  # Parada/Apprill (EMP)
     "fwd": "GTGYCAGCMGCCGCGGTAA", "rev": "GGACTACNVGGGTWTCTAAT"},
    {"name": "515F", "rev_name": "806R", "region": "16S V4",  # Caporaso 2011 (original)
     "fwd": "GTGCCAGCMGCCGCGGTAA", "rev": "GGACTACHVGGGTWTCTAAT"},
    {"name": "515F", "rev_name": "926R", "region": "16S V4-V5",  # EMP long
     "fwd": "GTGYCAGCMGCCGCGGTAA", "rev": "CCGYCAATTYMTTTRAGTTT"},
    {"name": "27F", "rev_name": "1492R", "region": "16S (near full length)",
     "fwd": "AGAGTTTGATCMTGGCTCAG", "rev": "TACGGYTACCTTGTTACGACTT"},
    # ── 18S rRNA (eukaryotes / protists) ──
    {"name": "TAReuk454FWD1", "rev_name": "TAReukREV3", "region": "18S V4",
     "fwd": "CCAGCASCYGCGGTAATTCC", "rev": "ACTTTCGTTCTTGATYRA"},
    {"name": "E572F", "rev_name": "E1009R", "region": "18S V4",  # Comeau 2011
     "fwd": "CYGCGGTAATTCCAGCTC", "rev": "AYGGTATCTRATCRTCTTYG"},
    {"name": "Euk1391F", "rev_name": "EukBr", "region": "18S V9",  # EMP
     "fwd": "GTACACACCGCCCGTC", "rev": "TGATCCTTCTGCAGGTTCACCTAC"},
    {"name": "1389F", "rev_name": "1510R", "region": "18S V9",  # Amaral-Zettler 2009
     "fwd": "TTGTACACACCGCCC", "rev": "CCTTCYGCAGGTTCACCTAC"},
    # ── ITS (fungi) ──
    {"name": "ITS1F", "rev_name": "ITS2", "region": "fungal ITS1",
     "fwd": "CTTGGTCATTTAGAGGAAGTAA", "rev": "GCTGCGTTCTTCATCGATGC"},
    {"name": "ITS1", "rev_name": "ITS4", "region": "fungal ITS (full)",  # White 1990
     "fwd": "TCCGTAGGTGAACCTGCGG", "rev": "TCCTCCGCTTATTGATATGC"},
    {"name": "ITS3", "rev_name": "ITS4", "region": "fungal ITS2",  # White 1990
     "fwd": "GCATCGATGAAGAACGCAGC", "rev": "TCCTCCGCTTATTGATATGC"},
    {"name": "gITS7", "rev_name": "ITS4", "region": "fungal ITS2",  # Ihrmark 2012
     "fwd": "GTGARTCATCGARTCTTTG", "rev": "TCCTCCGCTTATTGATATGC"},
]

_IUPAC_CHARS = set("ACGTRYSWKMBDHVN")
_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                         "primers", "tables")


def _insert_from_expected(expected, fwd: str, rev: str) -> int | None:
    """Insert length after primer removal, or None if not usable.

    Rejects values that cannot be a real amplicon: non-numeric, or an insert
    outside 50..2000bp once the primers come off. A wrong length here would
    silence a genuine overlap warning or raise a false one, so a missing value is
    better than a bad one.
    """
    try:
        total = int(float(str(expected).strip()))
    except (TypeError, ValueError):
        return None
    insert = total - len(fwd) - len(rev)
    return insert if 50 <= insert <= 2000 else None


def insert_length_for(fwd: str, rev: str) -> int | None:
    """Expected merged-fragment length for a primer pair, by exact sequence match."""
    f, r = fwd.strip().upper(), rev.strip().upper()
    for rec in PRIMER_DB:
        if rec.get("insert_length") and rec["fwd"] == f and rec["rev"] == r:
            return rec["insert_length"]
    return None


def _load_vendored_primers() -> list[dict]:
    """Parse the vendored FoodMicrobionet primer tables (16S + ITS).

    MIT-licensed data from github.com/ep142/FoodMicrobionet — see data/README.md.
    Schema: Target_region, primer_f_name, primer_f_seq, primer_r_name,
    primer_r_seq, reference, expected_length|notes. Skips rows that are empty,
    contain non-IUPAC characters (a stray typo), or are adapter-laden (a real
    metabarcoding primer is <=30 bp; longer entries carry sequencing adapters
    that wouldn't match demultiplexed reads).
    """
    out = []
    for fname, marker in (("primer_pairs_bacteria.txt", "16S"),
                          ("primer_pairs_fungi.txt", "ITS"),
                          ("primer_pairs_18S_pr2.txt", "18S")):
        path = os.path.join(_DATA_DIR, fname)
        try:
            with open(path, encoding="latin-1") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    fwd = (row.get("primer_f_seq") or "").strip().upper()
                    rev = (row.get("primer_r_seq") or "").strip().upper()
                    if not fwd or not rev:
                        continue
                    if len(fwd) > 30 or len(rev) > 30:
                        continue  # adapter/pad-laden, not a clean primer
                    if set(fwd) - _IUPAC_CHARS or set(rev) - _IUPAC_CHARS:
                        continue  # stray non-nucleotide character
                    region = (row.get("Target_region") or "").strip()
                    # Some tables already prefix the marker in Target_region
                    # (pr2's "18S V4"); others give a bare region ("V3-V4").
                    if region.upper().startswith(marker):
                        label = region
                    elif region:
                        label = f"{marker} {region}"
                    else:
                        label = marker
                    rec = {
                        "name": (row.get("primer_f_name") or "?").strip(),
                        "rev_name": (row.get("primer_r_name") or "?").strip(),
                        "region": label,
                        "fwd": fwd, "rev": rev,
                    }
                    # expected_length is the amplicon INCLUDING both primers, so
                    # the insert a merged pair has to span is that minus the two
                    # primers cutadapt removed. Only the 16S table carries it.
                    ins = _insert_from_expected(row.get("expected_length"), fwd, rev)
                    if ins:
                        rec["insert_length"] = ins
                    out.append(rec)
        except OSError as e:
            logging.getLogger(__name__).warning("primer table %s unreadable: %s", fname, e)
    return out


# A bare marker name is ambiguous across research communities: "16S" is the
# prokaryotic SSU rRNA to a microbial ecologist and the mitochondrial LSU rRNA
# (rrnL / mt-rnr2) to a zoologist barcoding animals, and "ITS" is fungal here but
# plant elsewhere. Anything reporting an assay therefore has to say which gene it
# means, so every entry carries the lineage its primers actually target rather
# than leaving a reader — or a model writing Methods — to assume.
_GENE_LINEAGE = {
    "16S": ("16S rRNA", "Bacteria/Archaea"),
    "18S": ("18S rRNA", "Eukaryota"),
    "ITS": ("ITS", "Fungi"),
}


def _split_region(label: str | None) -> tuple[str | None, str | None]:
    """Split a region label into (marker, sub-region).

    "16S V3-V4" -> ("16S", "V3-V4");  "ITS1" -> ("ITS", "1");  "18S" -> ("18S", None)
    """
    text = (label or "").strip()
    if not text:
        return None, None
    upper = text.upper()
    if "ITS" in upper:
        # "ITS1", "ITS2", "fungal ITS2" — the digit is the spacer, not a V-region.
        token = next(t for t in text.split() if "ITS" in t.upper())
        return "ITS", (token.upper().replace("ITS", "").strip() or None)
    marker, _, rest = text.partition(" ")
    return marker, (rest.strip() or None)


def _enrich(entry: dict) -> dict:
    """Add gene / lineage / sub-region alongside the collapsed `region` label."""
    marker, sub = _split_region(entry.get("region"))
    gene, lineage = _GENE_LINEAGE.get((marker or "").upper(), (marker, None))
    out = dict(entry)
    if gene:
        out["gene"] = gene
    if lineage:
        out["lineage"] = lineage
    if sub:
        out["sub_region"] = sub
    return out


def _build_primer_db() -> list[dict]:
    """Core (canonical, verified) primers first, then vendored ones deduped by
    sequence — so a pair we curated keeps its clean name over any FMBN variant."""
    db, seen = [], {}
    for p in _CORE_PRIMER_DB + _load_vendored_primers():
        key = (p["fwd"], p["rev"])
        if key in seen:
            # Same pair already held. Keep the curated name, but take an
            # insert_length the duplicate has and the winner lacks — only the
            # vendored 16S table carries lengths, so otherwise a curated pair
            # would silently lose the one field it cannot supply itself.
            if p.get("insert_length") and not seen[key].get("insert_length"):
                seen[key]["insert_length"] = p["insert_length"]
            continue
        rec = _enrich(p)
        seen[key] = rec
        db.append(rec)
    return db


PRIMER_DB = _build_primer_db()


def describe_pair(fwd_name: str | None, rev_name: str | None = None) -> dict | None:
    """Interpret an observed primer pair: which gene, whose lineage, what region.

    This is the OMC half of assay provenance (issue #57). The pipeline reports
    only what it can observe — that adapter `341Fv3` matched these reads — because
    a primer FASTA carries nothing but the name (and cutadapt truncates headers at
    whitespace, so it cannot be annotated into the log either). The curated table
    here is what turns that name into "bacterial/archaeal 16S rRNA, V3-V4".

    Matching prefers the exact pair, then the forward name alone, then the
    reverse. Returns None when nothing matches, and omits `region` when the pair
    resolves to conflicting sub-regions — an unknown assay must read as unknown
    rather than as a confident guess.
    """
    fwd = (fwd_name or "").strip()
    rev = (rev_name or "").strip()
    if not fwd and not rev:
        return None

    def _match(pred):
        return [p for p in PRIMER_DB if pred(p)]

    hits = []
    if fwd and rev:
        hits = _match(lambda p: p.get("name") == fwd and p.get("rev_name") == rev)
    if not hits and fwd:
        hits = _match(lambda p: p.get("name") == fwd)
    if not hits and rev:
        hits = _match(lambda p: p.get("rev_name") == rev)
    if not hits:
        return None

    genes = {p["gene"] for p in hits if p.get("gene")}
    if len(genes) != 1:
        return None  # the name is ambiguous across genes — say nothing

    out: dict = {"gene": genes.pop()}
    lineages = {p["lineage"] for p in hits if p.get("lineage")}
    if len(lineages) == 1:
        out["lineage"] = lineages.pop()
    subs = {p["sub_region"] for p in hits if p.get("sub_region")}
    if len(subs) == 1:
        out["region"] = subs.pop()
    return out

_DB_MATCH_MIN = 0.6   # min fraction of reads whose 5' matches a DB forward primer
_CONSENSUS_STOP_ENTROPY = 1.7  # bits; above this a column is "biological", stop


def _iupac_regex(seq: str) -> re.Pattern:
    """Compile an IUPAC-aware regex matching `seq` against plain-ACGT reads."""
    return re.compile("".join(f"[{IUPAC.get(b, b)}]" for b in seq.upper()))


def sample_reads(fastq_path: str, n: int = 500) -> list[str]:
    """Return up to `n` sequence lines from a (optionally gzipped) FASTQ."""
    reads: list[str] = []
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    try:
        with opener(fastq_path, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    reads.append(line.strip().upper())
                    if len(reads) >= n:
                        break
    except (OSError, EOFError):
        pass
    return reads


def _match_fraction(reads: list[str], primer: str, max_offset: int = 3) -> float:
    """Fraction of reads whose 5' end matches `primer` (small offset allowed)."""
    if not reads:
        return 0.0
    rx = _iupac_regex(primer)
    L = len(primer)
    hits = 0
    for r in reads:
        for off in range(max_offset + 1):
            seg = r[off:off + L]
            if len(seg) == L and rx.match(seg):
                hits += 1
                break
    return hits / len(reads)


def _consensus_primer(reads: list[str], length: int = 25, cover: float = 0.9) -> str:
    """De-novo degenerate consensus of the first `length` bp across `reads`.

    At each position, pick the smallest set of bases covering `cover` of the
    reads and emit its IUPAC code. Stop at the first high-entropy column — that
    is where the conserved primer ends and the biological (variable) region
    begins.
    """
    out: list[str] = []
    for i in range(length):
        col = [r[i] for r in reads if len(r) > i and r[i] in "ACGT"]
        if len(col) < max(10, 0.5 * len(reads)):
            break
        counts = Counter(col)
        tot = len(col)
        ent = -sum((c / tot) * math.log2(c / tot) for c in counts.values())
        bases: list[str] = []
        acc = 0
        for b, c in counts.most_common():
            bases.append(b)
            acc += c
            if acc / tot >= cover:
                break
        if ent > _CONSENSUS_STOP_ENTROPY and len(bases) >= 3:
            break
        out.append(_IUPAC_REV.get(frozenset(bases), "N"))
    return "".join(out)


def detect_from_reads(r1_path: str, r2_path: str | None = None, n: int = 500) -> dict | None:
    """Tier 3. Infer primers from a sample of reads (given FASTQ paths)."""
    r1 = sample_reads(r1_path, n)
    r2 = sample_reads(r2_path, n) if r2_path else []
    return detect_from_read_lists(r1, r2)


def detect_from_read_lists(r1: list[str], r2: list[str] | None = None) -> dict | None:
    """Tier 3 core, operating on already-sampled reads.

    Split out from detect_from_reads so callers that have reads in hand (e.g.
    detect_primer_sets, which re-probes the same sample against several
    candidate sets) don't re-read the FASTQ each time.

    Returns a primer dict, or None if there aren't enough reads to try.
    """
    r2 = r2 or []
    if len(r1) < 20:
        return None

    # 3a — database match: score each pair by forward match on R1 (and, if we
    # have R2, reverse match on R2), keep the best.
    best = None
    for p in PRIMER_DB:
        fwd_frac = _match_fraction(r1, p["fwd"])
        rev_frac = _match_fraction(r2, p["rev"]) if r2 else None
        score = fwd_frac if rev_frac is None else (fwd_frac + rev_frac) / 2
        if best is None or score > best["score"]:
            best = {**p, "score": score, "fwd_frac": fwd_frac, "rev_frac": rev_frac}
    if best and best["fwd_frac"] >= _DB_MATCH_MIN:
        return {
            "fwd": best["fwd"], "rev": best["rev"],
            "fwd_name": best["name"], "rev_name": best["rev_name"],
            "region": best["region"], "source": "inferred-db",
            "confidence": round(best["score"], 2),
        }

    # 3b — de-novo consensus of the conserved 5' prefix.
    fwd = _consensus_primer(r1)
    rev = _consensus_primer(r2) if r2 else ""
    if len(fwd) < 8:
        return None
    return {
        "fwd": fwd, "rev": rev, "fwd_name": "inferred", "rev_name": "inferred",
        "region": "unknown", "source": "inferred-denovo", "confidence": None,
    }
