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
import subprocess
import tempfile
from collections import Counter
from functools import lru_cache

# IUPAC nucleotide codes → the set of bases each represents.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
# Reverse map: frozenset of bases → the tightest IUPAC code.
_IUPAC_REV = {frozenset(v): k for k, v in IUPAC.items()}
# Set form, for testing whether two IUPAC codes can encode a common base.
_SETS = {k: frozenset(v) for k, v in IUPAC.items()}

# Bases of each read the consensus is built from. Long enough to hold any
# catalogue primer plus a barcode or spacer in front of it, short enough that the
# alignment is dominated by the primer rather than by the amplicon behind it.
_CONSENSUS_READ_LEN = 50

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


_MIN_PRIMER_OVERLAP = 15


def _same_primer(a: str, b: str) -> bool:
    """Do these two strings denote the same primer?

    Exact match, or one is the other with bases trimmed off an end. The same
    primer is routinely written both ways: 341F is published as CCTACGGGNGGCWGCAG
    but is carried in submissions, and written into primer fastas, as the 16bp
    CTACGGGNGGCWGCAG. Requiring exact equality means the commonest 16S V3-V4
    primer on earth misses its own database entry.

    Floored at _MIN_PRIMER_OVERLAP so a short fragment cannot match half the table.
    """
    if a == b:
        return True
    if min(len(a), len(b)) < _MIN_PRIMER_OVERLAP:
        return False
    return (a.endswith(b) or b.endswith(a) or a.startswith(b) or b.startswith(a))


def insert_length_for(fwd: str, rev: str) -> int | None:
    """Expected merged-fragment length for a primer pair.

    Exact matches win. Failing that, end-trimmed variants are accepted — but only
    if every such match agrees on the length. Two entries offering different
    answers means the pair is ambiguous, and a wrong fragment length is worse than
    none: it would silence a real overlap warning or invent a false one.
    """
    f, r = fwd.strip().upper(), rev.strip().upper()
    with_len = [rec for rec in PRIMER_DB if rec.get("insert_length")]

    for rec in with_len:
        if rec["fwd"] == f and rec["rev"] == r:
            return rec["insert_length"]

    lengths = {rec["insert_length"] for rec in with_len
               if _same_primer(rec["fwd"], f) and _same_primer(rec["rev"], r)}
    return lengths.pop() if len(lengths) == 1 else None


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

# A primer we infer must look like one the field actually uses, so the bounds
# come from the catalogue rather than from taste.
_DEGENERATE = set("RYSWKMBDHVN")
MIN_PRIMER_LEN = min(len(s) for p in PRIMER_DB for s in (p["fwd"], p["rev"]) if s)
MAX_PRIMER_LEN = max(len(s) for p in PRIMER_DB for s in (p["fwd"], p["rev"]) if s)
MAX_PRIMER_DEGENERACY = max(
    sum(1 for c in s if c in _DEGENERATE) / len(s)
    for p in PRIMER_DB for s in (p["fwd"], p["rev"]) if s)

# A candidate is accepted on the evidence of the reads, not on its alignment
# score, so the score only has to be loose enough not to miss the real primer.
_CANDIDATE_MIN_ANCHORS = 10   # unambiguous consensus positions a candidate must span
_CANDIDATE_MIN_ID = 0.90      # agreement across those positions
_VERIFY_MIN = 0.30            # fraction of reads that must carry it to accept it
_PRIMER_MAX_START = 20        # a primer sits at the 5' end, behind at most a barcode

# Reads at the head of a FASTQ come from one corner of one tile and run roughly a
# Phred below the body of the file; the difference is gone by ~1000 reads.
_SAMPLE_SKIP = 1000
_SAMPLE_MIN_MEAN_Q = 30
_SAMPLE_MIN_BASE_Q = 15


def _degeneracy(seq: str) -> float:
    return sum(1 for c in seq if c in _DEGENERATE) / len(seq) if seq else 1.0


def sample_reads(fastq_path: str, n: int = 500, trim: int | None = None,
                 skip: int = _SAMPLE_SKIP, min_mean_q: int = _SAMPLE_MIN_MEAN_Q,
                 min_base_q: int = _SAMPLE_MIN_BASE_Q) -> list[str]:
    """Return up to `n` good reads from a (optionally gzipped) FASTQ.

    Skips the first `skip` reads and drops any read carrying an N, a mean quality
    below `min_mean_q`, or any base below `min_base_q`, so the consensus is built
    from the body of the file rather than its weakest corner. Falls back to
    whatever it found if the file is too small to satisfy the filters.
    """
    reads: list[str] = []
    fallback: list[str] = []
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    try:
        with opener(fastq_path, "rt") as f:
            seen = 0
            while len(reads) < n:
                if not f.readline():
                    break
                seq = f.readline().strip().upper()
                f.readline()
                qual = f.readline().strip()
                if not seq or not qual:
                    break
                seen += 1
                if trim:
                    seq, qual = seq[:trim], qual[:trim]
                    if len(seq) < trim:
                        continue
                if len(fallback) < n:
                    fallback.append(seq)
                if seen <= skip or "N" in seq:
                    continue
                phred = [ord(c) - 33 for c in qual]
                if not phred or sum(phred) / len(phred) < min_mean_q:
                    continue
                if min(phred) < min_base_q:
                    continue
                reads.append(seq)
    except (OSError, EOFError):
        pass
    return reads if len(reads) >= 20 else fallback


def derived_primer_name(seq: str) -> str:
    """Name for a primer we inferred rather than recognised: its own sequence.

    The name is the only field that survives the trip to samples.json — cutadapt
    reports the FASTA header it matched and nothing else — so whatever identity a
    de-novo primer has must live in the name. The sequence is that identity, and
    carrying it directly means a reader can see what was trimmed without
    resolving anything, while two runs of the same primer still line up because
    identical sequences give identical names.

    The name is also the grouping key: plateFor() folds it into the key for
    AUTO_TRIM and LEARN_ERRORS, so samples trimmed with the same sequence share
    an error model and samples trimmed with different ones do not.
    """
    return re.sub(r"\s+", "", (seq or "")).upper()


def _mafft(seqs: list[str], timeout: int = 120) -> list[str]:
    """Align `seqs` with mafft, or return [] if it is unavailable or fails."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for i, s in enumerate(seqs):
            fh.write(f">r{i}\n{s}\n")
        path = fh.name
    try:
        proc = subprocess.run(["mafft", "--quiet", "--retree", "1", path],
                              capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.SubprocessError):
        return []
    finally:
        os.unlink(path)
    if proc.returncode != 0:
        return []
    aln, cur = [], []
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            if cur:
                aln.append("".join(cur))
                cur = []
        else:
            cur.append(line.strip().upper())
    if cur:
        aln.append("".join(cur))
    return aln


def _consensus_primer(reads: list[str], cover: float = 0.9,
                      min_mean_info: float = 1.2) -> str:
    """De-novo degenerate consensus of the conserved 5' block across `reads`.

    Aligns the reads to each other first. A positional consensus assumes every
    read begins at the primer, which heterogeneity spacers and untrimmed inline
    barcodes both violate — out-of-phase reads then average into a string of
    degenerate codes that matches nothing. Under an alignment the conserved block
    finds its own frame, so the primer survives whatever precedes it.

    Scores a window by total information rather than mean: a barcode is short and
    perfectly conserved, and mean information prefers it to the longer primer
    behind it. The window must also start near the 5' end, since adapter
    read-through is every bit as conserved as a primer and lands mid-read on
    short inserts.
    """
    aln = _mafft(reads)
    if not aln:
        return ""
    width = len(aln[0])
    cons: list[str] = []
    info: list[float] = []
    for i in range(width):
        col = [s[i] for s in aln if s[i] in "ACGT"]
        if len(col) < 0.5 * len(aln):
            cons.append("-")
            info.append(-1.0)
            continue
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
        cons.append(_IUPAC_REV.get(frozenset(bases), "N"))
        info.append(2.0 - ent)

    best = None
    for length in range(MIN_PRIMER_LEN, MAX_PRIMER_LEN + 1):
        for start in range(min(_PRIMER_MAX_START, len(info) - length) + 1):
            window = info[start:start + length]
            if any(x < 0 for x in window):
                continue
            if sum(window) / length < min_mean_info:
                continue
            if best is None or sum(window) > best[0]:
                best = (sum(window), start, length)
    if not best:
        return ""
    return "".join(cons[best[1]:best[1] + best[2]])


def read_support(reads: list[str], primer: str, max_offset: int = _PRIMER_MAX_START,
                 min_id: float = 0.90) -> float:
    """Fraction of `reads` carrying `primer` within the 5' window."""
    if not reads or not primer:
        return 0.0
    length = len(primer)
    hits = 0
    for r in reads:
        for off in range(max_offset + 1):
            seg = r[off:off + length]
            if len(seg) < length:
                break
            ok = sum(1 for a, b in zip(seg, primer)
                     if _SETS.get(a, _SETS["N"]) & _SETS.get(b, _SETS["N"]))
            if ok / length >= min_id:
                hits += 1
                break
    return hits / len(reads)


def catalogue_candidates(seq: str) -> list[tuple[str, str, int, float]]:
    """Catalogue primers compatible with `seq`, best first.

    Scores only the positions where the consensus commits to a single base. A
    degenerate code is compatible with almost anything, so counting it as
    agreement lets a mostly-degenerate consensus match the whole catalogue.

    Deliberately permissive: every candidate is checked against the reads before
    it is accepted, so the cost of a loose score here is a wasted check, while the
    cost of a tight one is missing the primer that was actually used.
    """
    best: dict[tuple[str, str], tuple[int, float]] = {}
    for p in PRIMER_DB:
        for key, name in (("fwd", p["name"]), ("rev", p["rev_name"])):
            ref = p[key]
            if not ref or not name:
                continue
            for si in range(len(seq) - _CANDIDATE_MIN_ANCHORS + 1):
                for ri in range(len(ref) - _CANDIDATE_MIN_ANCHORS + 1):
                    n = min(len(seq) - si, len(ref) - ri)
                    idx = [k for k in range(n) if seq[si + k] not in _DEGENERATE]
                    if len(idx) < _CANDIDATE_MIN_ANCHORS:
                        continue
                    ok = sum(1 for k in idx
                             if seq[si + k] in IUPAC.get(ref[ri + k], ""))
                    score = ok / len(idx)
                    if score < _CANDIDATE_MIN_ID:
                        continue
                    cur = best.get((name, ref))
                    if cur is None or (len(idx), score) > cur:
                        best[(name, ref)] = (len(idx), score)
    return sorted(((name, ref, a, s) for (name, ref), (a, s) in best.items()),
                  key=lambda t: (-t[2], -t[3]))


def resolve_primer(consensus: str, reads: list[str]) -> tuple[str, str, dict]:
    """Settle on the sequence to trim with: (sequence, source, detail).

    The catalogue comes first — a recognised primer carries the right ends,
    whereas a consensus runs on into conserved amplicon and would trim real
    sequence away — but a name is not evidence. Candidates are tried best-first
    and the first the reads confirm wins; a primer that is not in the data would
    otherwise reach cutadapt, and with --discard-untrimmed that empties the run.

    Only when nothing in the catalogue survives is a primer invented, and then it
    must be at least as long and no more degenerate than the catalogue allows.
    """
    if not consensus:
        return "", "none", {"reason": "no conserved 5' block"}
    tried = []
    for name, ref, anchors, score in catalogue_candidates(consensus):
        obs = read_support(reads, ref)
        tried.append({"name": name, "anchors": anchors,
                      "id": round(score, 3), "read_support": round(obs, 3)})
        if obs >= _VERIFY_MIN:
            return ref, "catalogue", {"name": name, "anchors": anchors,
                                      "id": round(score, 3),
                                      "read_support": round(obs, 3),
                                      "rejected": tried[:-1]}
    if len(consensus) < MIN_PRIMER_LEN:
        return "", "none", {"reason": f"consensus {len(consensus)}bp is shorter than "
                                      f"the shortest catalogue primer ({MIN_PRIMER_LEN}bp)",
                            "rejected": tried}
    degen = _degeneracy(consensus)
    if degen > MAX_PRIMER_DEGENERACY:
        return "", "none", {"reason": f"consensus degeneracy {degen:.2f} exceeds the "
                                      f"most degenerate catalogue primer "
                                      f"({MAX_PRIMER_DEGENERACY:.2f})",
                            "rejected": tried}
    return consensus, "de-novo", {"degeneracy": round(degen, 3), "rejected": tried}


def detect_from_reads(r1_path: str, r2_path: str | None = None, n: int = 500) -> dict | None:
    """Infer the primers a sample carries, from the reads themselves."""
    r1 = sample_reads(r1_path, n, trim=_CONSENSUS_READ_LEN)
    r2 = sample_reads(r2_path, n, trim=_CONSENSUS_READ_LEN) if r2_path else []
    return detect_from_read_lists(r1, r2)


@lru_cache(maxsize=1)
def _ssu_references() -> dict[str, str]:
    """The E. coli 16S and yeast 18S sequences primer coordinates are quoted in."""
    path = os.path.join(_DATA_DIR, "ssu_references.fa")
    refs: dict[str, list[str]] = {}
    name = None
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    refs[name] = []
                elif name:
                    refs[name].append(line.strip().upper())
    except OSError:
        return {}
    return {k: "".join(v) for k, v in refs.items()}


def locate_on_ssu(primer: str, min_id: float = 0.80) -> list[dict]:
    """Where `primer` sits on the SSU references, best match first.

    A primer's name already is its coordinate — 515F is the forward primer at
    E. coli position 515 — so the position is the identity the catalogue only
    approximates, and it is derived rather than looked up. Matching both
    references also lets a universal primer say so instead of being forced into
    one gene: the SSU V4 region is homologous across domains, which is why one
    sequence is catalogued as 515F for 16S and V4_1f for 18S.

    An empty result is the useful case: a primer that lands on neither gene is
    not a ribosomal primer, and the run should not be handed an rRNA reference
    database.
    """
    if not primer:
        return []
    out = []
    for ref_name, ref in _ssu_references().items():
        length = len(primer)
        best = (0.0, None)
        for i in range(len(ref) - length + 1):
            ok = sum(1 for k in range(length)
                     if ref[i + k] in _SETS.get(primer[k], _SETS["N"]))
            score = ok / length
            if score > best[0]:
                best = (score, i + 1)
        if best[0] >= min_id:
            out.append({"reference": ref_name, "identity": round(best[0], 3),
                        "start": best[1], "end": best[1] + length - 1})
    return sorted(out, key=lambda d: -d["identity"])


def detect_from_read_lists(r1: list[str], r2: list[str] | None = None) -> dict | None:
    """Infer a sample's primers from reads already in hand.

    Split out from detect_from_reads so callers holding reads don't re-read the
    FASTQ. Returns a primer dict, or None when the reads yield nothing usable.
    """
    r2 = r2 or []
    if len(r1) < 20:
        return None

    fwd, fwd_src, fwd_detail = resolve_primer(_consensus_primer(r1), r1)
    if not fwd:
        return None
    rev, rev_src, rev_detail = ("", "none", {})
    if r2:
        rev, rev_src, rev_detail = resolve_primer(_consensus_primer(r2), r2)

    fwd_loc = locate_on_ssu(fwd)
    rev_loc = locate_on_ssu(rev) if rev else []
    # Both ends landing nowhere on either SSU gene means the assay is not
    # ribosomal, whatever the reads denoise into. Taxonomy against an rRNA
    # database would still return confident labels, and they would be noise.
    ribosomal = bool(fwd_loc or rev_loc)

    return {
        "fwd": fwd, "rev": rev,
        "fwd_name": fwd if fwd_src == "de-novo" else fwd_detail.get("name", fwd),
        "rev_name": rev if rev_src == "de-novo" else rev_detail.get("name", rev),
        "region": "unknown",
        "source": f"inferred-{fwd_src}",
        "confidence": fwd_detail.get("read_support"),
        "ribosomal": ribosomal,
        "fwd_detail": fwd_detail, "rev_detail": rev_detail,
        "fwd_location": fwd_loc, "rev_location": rev_loc,
    }
