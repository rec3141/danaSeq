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
import bisect
import math
import os
import random
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


# How much of two primers has to line up before they are called the same primer.
# The catalogue holds the same primer at several lengths — TAReukFWD1 is
# TAReuk454FWD1 three bases shorter — so agreement is judged over the shorter of
# the two, from the 5' end where both are anchored.
_MIN_SEQ_OVERLAP = 15


def _same_primer(observed: str, catalogued: str) -> bool:
    """True when an observed sequence is the catalogue's primer.

    Degenerate positions are compared by whether the two codes can stand for a
    common base, not by whether they are written the same way: a consensus taken
    from reads records the bases those reads actually carried, which is at most
    what the ordered primer allowed and is usually less.
    """
    if not observed or not catalogued:
        return False
    n = min(len(observed), len(catalogued))
    if n < _MIN_SEQ_OVERLAP:
        return False
    for i in range(n):
        if not (set(IUPAC.get(observed[i], "")) & set(IUPAC.get(catalogued[i], ""))):
            return False
    return True


def _end_matches(entry: dict, token: str, end: str, exact: bool) -> bool:
    """Does `token` identify this catalogue entry's forward or reverse primer?

    A token is whatever the run recorded, which is a catalogue name when the
    primers were supplied and the primer's own sequence when they were detected
    from the reads. Both are accepted because the caller cannot tell them apart
    and should not have to.

    `exact` distinguishes the two strengths of evidence. Compatibility is a wide
    net — one degenerate position is enough to make a primer for one gene answer
    to another's — and the catalogue holds near-duplicates that differ by a
    single code, so it can only be trusted once nothing exact has been found.
    """
    name_key, seq_key = ("name", "fwd") if end == "fwd" else ("rev_name", "rev")
    if entry.get(name_key) == token:
        return True
    seq = (entry.get(seq_key) or "").upper()
    token = token.upper()
    return seq == token if exact else _same_primer(token, seq)


def _end_gene(token: str, end: str) -> str | None:
    """The gene one end alone points at, or None if it is silent or ambiguous."""
    if not token:
        return None
    for exact in (True, False):
        hits = [p for p in PRIMER_DB if _end_matches(p, token, end, exact)]
        if hits:
            genes = {p["gene"] for p in hits if p.get("gene")}
            return genes.pop() if len(genes) == 1 else None
    return None


def describe_pair(fwd: str | None, rev: str | None = None) -> dict | None:
    """Interpret an observed primer pair: which gene, whose lineage, what region.

    This is the OMC half of assay provenance (issue #57). The pipeline reports
    only what it can observe, and what that is depends on how the run was set up:
    a supplied primer is known by the name in the cutadapt log, while a detected
    one is known only by the sequence the reads agreed on, since there was never
    a name to record. Both arrive here and both are looked up, because a detected
    primer is very often a catalogue primer that nobody typed the name of
    (danaSeq #53). The curated table is what turns either into "bacterial/archaeal
    16S rRNA, V3-V4".

    Matching prefers the exact pair, then the forward alone, then the reverse.
    Returns None when nothing matches, and omits `region` when the pair resolves
    to conflicting sub-regions — an unknown assay must read as unknown rather
    than as a confident guess.
    """
    fwd = (fwd or "").strip()
    rev = (rev or "").strip()
    if not fwd and not rev:
        return None

    def _match(pred):
        return [p for p in PRIMER_DB if pred(p)]

    # Strongest evidence first: both ends exactly is one assay, one end
    # compatibly is a family of them, and the checks below decide whether that
    # family still agrees on an answer.
    hits = []
    for exact in (True, False):
        if fwd and rev:
            hits = _match(lambda p: _end_matches(p, fwd, "fwd", exact)
                          and _end_matches(p, rev, "rev", exact))
        if hits:
            break

    # No pair matched, and the two ends may be why. A forward for one gene and a
    # reverse for another is a real thing to find in a deposited run, and
    # answering with the forward's assay would state an amplicon that was never
    # amplified. Resolve the ends separately and say so when they disagree.
    if not hits:
        fwd_gene = _end_gene(fwd, "fwd")
        rev_gene = _end_gene(rev, "rev")
        if fwd_gene and rev_gene and fwd_gene != rev_gene:
            return {"gene_conflict": f"{fwd_gene} forward with {rev_gene} reverse",
                    "fwd_gene": fwd_gene, "rev_gene": rev_gene}
        for exact in (True, False):
            if fwd:
                hits = _match(lambda p: _end_matches(p, fwd, "fwd", exact))
            if not hits and rev:
                hits = _match(lambda p: _end_matches(p, rev, "rev", exact))
            if hits:
                break
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
_CATALOGUE_LENS = sorted(
    len(s) for s in {s for p in PRIMER_DB for s in (p["fwd"], p["rev"]) if s})
MIN_PRIMER_LEN = _CATALOGUE_LENS[0]

# How long a primer we are willing to invent, which is not the same as the
# longest one on file. The catalogue's longest entry is a tagged ITS primer whose
# tag is counted in that length, and two thirds of the table sits between 17 and
# 22bp — searching out to the extreme means the window keeps extending past the
# primer into the amplicon, where the added columns are variable. That inflates
# the consensus' degeneracy until it is rejected outright, and where it is not,
# cutadapt trims the extra bases off the front of every read.
MAX_DENOVO_LEN = _CATALOGUE_LENS[int(0.95 * (len(_CATALOGUE_LENS) - 1))]

# A candidate is accepted on the evidence of the reads, not on its alignment
# score, so the score only has to be loose enough not to miss the real primer.
_CANDIDATE_MIN_ANCHORS = 10   # unambiguous consensus positions a candidate must span
_CANDIDATE_MIN_ID = 0.90      # agreement across those positions
_VERIFY_MIN = 0.30            # fraction of reads that must carry it to accept it
_PRIMER_MAX_START = 20        # a primer sits at the 5' end, behind at most a barcode

# Reads at the head of a FASTQ come from one corner of one tile and run roughly a
# Phred below the body of the file; the difference is gone by ~1000 reads.
_SAMPLE_SKIP = 1000
# Reservoir sampling needs a source of randomness; a fixed one keeps detection
# reproducible, which -resume and provenance both depend on.
_SAMPLE_SEED = 1337
_SAMPLE_MIN_MEAN_Q = 30
_SAMPLE_MIN_BASE_Q = 15


def _degeneracy(seq: str) -> float:
    """Bits of ambiguity per base: 0 for a fixed base, 1 for Y, 2 for N.

    Counting degenerate positions instead prices a Y the same as an N, which puts
    a two-fold 28-mer and a four-fold one on the same side of any cap that
    separates them by a single position — so which of two samples carrying the
    same primer is accepted comes down to sampling noise.
    """
    if not seq:
        return 2.0
    return sum(math.log2(len(IUPAC.get(c, "ACGT"))) for c in seq) / len(seq)


# A consensus is only worth trimming with if it is as specific as the primers the
# field actually uses, so the ceiling is the most ambiguous one in the catalogue.
MAX_PRIMER_DEGENERACY = max(
    _degeneracy(s) for p in PRIMER_DB for s in (p["fwd"], p["rev"]) if s)


def sample_reads(fastq_path: str, n: int = 500, trim: int | None = None,
                 skip: int = _SAMPLE_SKIP, min_mean_q: int = _SAMPLE_MIN_MEAN_Q,
                 min_base_q: int = _SAMPLE_MIN_BASE_Q) -> list[str]:
    """Return up to `n` good reads from a (optionally gzipped) FASTQ.

    Skips the first `skip` reads and drops any read carrying an N, a mean quality
    below `min_mean_q`, or any base below `min_base_q`, so the consensus is built
    from the body of the file rather than its weakest corner. Falls back to
    whatever it found if the file is too small to satisfy the filters.

    Sampled across the whole file rather than from the first `n` reads that pass.
    A FASTQ is not in random order, and a mixed library need not be interleaved:
    PRJNA779070 writes its 16S and 18S amplicons in contiguous runs, so a window
    taken near the start reports one of them and never learns the other exists —
    SRR16930999 reads as 16S from its first 1,500 reads and is 72% 18S.

    Reservoir sampling, on a fixed seed so the same file always yields the same
    reads: detection has to be reproducible across a -resume, and a consensus
    that moves between runs is not evidence of anything.
    """
    reads: list[str] = []
    fallback: list[str] = []
    rng = random.Random(_SAMPLE_SEED)
    kept = 0
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    try:
        with opener(fastq_path, "rt") as f:
            seen = 0
            while True:
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
                kept += 1
                if len(reads) < n:
                    reads.append(seq)
                else:
                    j = rng.randrange(kept)
                    if j < n:
                        reads[j] = seq
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
                      max_mean_ambiguity: float = MAX_PRIMER_DEGENERACY) -> str:
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

    Total information rises with every column added, so the window is held back
    by what it may not exceed rather than by what it maximises: the degeneracy of
    the codes it would emit, measured the way resolve_primer measures it. Past the
    primer the columns are amplicon and the mean climbs, which is where it stops.
    """
    aln = _mafft(reads)
    if not aln:
        return ""
    width = len(aln[0])
    cons: list[str] = []
    info: list[float] = []
    amb: list[float] = []
    # Columns most reads do not occupy are insertions carried by a few of them,
    # not part of the primer, and they grow with the sample: at 50 reads one
    # falls inside the primer, at 500 half a dozen do. Drop them and close the
    # gap, so the primer stays one run of columns however many reads are read.
    for i in range(width):
        col = [s[i] for s in aln if s[i] in "ACGT"]
        if len(col) < 0.5 * len(aln):
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
        code = _IUPAC_REV.get(frozenset(bases), "N")
        cons.append(code)
        info.append(2.0 - ent)
        amb.append(math.log2(len(IUPAC[code])))

    # Earliest workable start first, and only then the most informative length
    # from it. Ranking by information alone lets the window slide 3', because
    # trading a degenerate 5' column for a clean one further in always scores
    # better — and it slides by a different amount in each sample, so the primers
    # no longer line up and collapse_primers grinds the shared core down to
    # nothing. Where the start has to move it is the degeneracy bound that moves
    # it, past a barcode or a heterogeneity spacer, which is the only reason a
    # primer is not already at the 5' end.
    best = None
    for start in range(_PRIMER_MAX_START + 1):
        for length in range(MIN_PRIMER_LEN, MAX_DENOVO_LEN + 1):
            if start + length > len(info):
                break
            if sum(amb[start:start + length]) / length > max_mean_ambiguity:
                continue
            total = sum(info[start:start + length])
            if best is None or total > best[0]:
                best = (total, start, length)
        if best:
            break
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


def best_supported_catalogue_primer(reads: list[str], end: str = "fwd",
                                    min_support: float = _VERIFY_MIN) -> tuple[str, str, float]:
    """Catalogue primer the most reads actually carry: (sequence, name, support).

    A consensus describes the average read, which is only the primer when most
    reads carry one. A library that is half primer-bearing and half adapter
    read-through averages into something too degenerate to use, while the primer
    is still plainly present in half the reads — so ask each catalogue primer
    directly how many reads carry it, rather than asking what the reads agree on.
    """
    best = ("", "", 0.0)
    seen: set[str] = set()
    for p in PRIMER_DB:
        ref = p[end]
        name = p["name"] if end == "fwd" else p["rev_name"]
        if not ref or ref in seen:
            continue
        seen.add(ref)
        obs = read_support(reads, ref)
        if obs > best[2]:
            best = (ref, name or "", obs)
    return best if best[2] >= min_support else ("", "", best[2])


def _template_consumed(primer: str) -> int:
    """How much of the gene a primer eats, for choosing between candidates.

    Lower is better, and the number only means anything against other candidates
    for the same end. A forward primer consumes up to its last base, so its end
    coordinate is the measure. A reverse primer is written as the other strand
    and extends leftwards along the reference, so the base it consumes to is its
    *start* coordinate and a higher one leaves more amplicon — negated so that
    smaller still means better.

    A primer that places nowhere sorts last: nothing can be said about where it
    stops, and a placed candidate is the better evidence.
    """
    locs = locate_on_ssu(primer)
    if not locs:
        return 10 ** 6
    return min((l["end"] if l["strand"] == "+" else -l["start"]) for l in locs)

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
    tried, confirmed = [], []
    for name, ref, anchors, score in catalogue_candidates(consensus) if consensus else []:
        obs = read_support(reads, ref)
        record = {"name": name, "anchors": anchors,
                  "id": round(score, 3), "read_support": round(obs, 3)}
        tried.append(record)
        if obs >= _VERIFY_MIN:
            confirmed.append((len(ref), -score, name, ref, anchors, obs, record))

    # Among those the reads confirm, the one that cuts least far into the gene.
    #
    # The catalogue holds one primer site at several lengths: 341F ends at E. coli
    # 357, and 341Fmod is the same primer plus 8 bases, ending at 365. Reads
    # cannot separate them — the bases past a primer are conserved across the
    # community, so they are as invariant as the primer itself, and both clear
    # the threshold on the same sample (0.986 and 0.980).
    #
    # The two mistakes are not equal. Cut too far and the surplus comes out of
    # template: real variation there is deleted, and the reads carrying it fail
    # the match and are discarded — 11 points of yield on PRJNA661323. Cut too
    # short and a constant remnant stays on every read, which denoising collapses
    # to a single sequence and which costs nothing but where the amplicon starts.
    # One loses data, the other keeps it (danaSeq #59).
    #
    # It is the 3' end that decides this, not the length: 343F is shorter than
    # 341F only because it starts later, and starting later removes nothing extra
    # — cutadapt's -g takes everything ahead of the match regardless. So the
    # candidate eating least template wins, and among those reaching equally far
    # the longest, being the most specific match.
    if confirmed:
        confirmed.sort(key=lambda c: (_template_consumed(c[3]), -len(c[3]), c[1]))
        _len, _negscore, name, ref, anchors, obs, record = confirmed[0]
        return ref, "catalogue", {"name": name, "anchors": anchors,
                                  "id": record["id"],
                                  "read_support": record["read_support"],
                                  "rejected": [t for t in tried if t is not record]}

    # The consensus is the average read, and averages a mixed library into
    # something no candidate resembles. Ask the catalogue directly instead.
    seq, name, obs = best_supported_catalogue_primer(reads)
    if seq:
        return seq, "catalogue-scan", {"name": name, "read_support": round(obs, 3),
                                       "consensus": consensus, "rejected": tried}

    if not consensus:
        return "", "none", {"reason": "no conserved 5' block and no catalogue "
                                      "primer carried by enough reads",
                            "rejected": tried}
    if len(consensus) < MIN_PRIMER_LEN:
        return "", "none", {"reason": f"consensus {len(consensus)}bp is shorter than "
                                      f"the shortest catalogue primer ({MIN_PRIMER_LEN}bp)",
                            "rejected": tried}
    degen = _degeneracy(consensus)
    if degen > MAX_PRIMER_DEGENERACY:
        return "", "none", {"reason": f"consensus carries {degen:.2f} bits of "
                                      f"ambiguity per base, more than the most "
                                      f"degenerate catalogue primer "
                                      f"({MAX_PRIMER_DEGENERACY:.2f})",
                            "rejected": tried}
    # read_support on the de-novo branch too: it is what says how many reads
    # actually carry the block, and a run whose reads are half in one orientation
    # and half in the other sits near a half on both ends.
    return consensus, "de-novo", {"degeneracy": round(degen, 3),
                                  "read_support": round(read_support(reads, consensus), 3),
                                  "rejected": tried}


def detect_from_reads(r1_path: str, r2_path: str | None = None, n: int = 500) -> dict | None:
    """Infer the primers a sample carries, from the reads themselves."""
    r1 = sample_reads(r1_path, n, trim=_CONSENSUS_READ_LEN)
    r2 = sample_reads(r2_path, n, trim=_CONSENSUS_READ_LEN) if r2_path else []
    return detect_from_read_lists(r1, r2)


def shared_region(a: str, b: str, min_id: float = 0.90,
                  min_len: int = MIN_PRIMER_LEN) -> str:
    """The longest stretch two primers agree on, over every relative placement.

    Records of the same primer differ at the 5' end and agree after it, because
    what varies between samples is the barcode or spacer in front rather than the
    primer itself. Scanning every placement — not just a common prefix — finds
    the agreed stretch, and returning it is what drops the barcode.
    """
    best = ""
    for ai in range(len(a) - min_len + 1):
        for bi in range(len(b) - min_len + 1):
            n = min(len(a) - ai, len(b) - bi)
            if n < min_len or n <= len(best):
                continue
            ok = sum(1 for k in range(n)
                     if _SETS.get(a[ai + k], _SETS["N"]) & _SETS.get(b[bi + k], _SETS["N"]))
            if ok / n >= min_id:
                best = a[ai:ai + n]
    return best


def collapse_primers(primers: list[str]) -> list[tuple[str, int]]:
    """Reduce per-sample primers to the smallest set that still separates assays.

    Every extra entry in the FASTA is another adapter cutadapt can pick per
    sample, and plateFor() keys the error-model groups on which one it picked —
    so one entry per sample splits a single assay into as many groups as there
    are samples. Compatible primers are merged and each group keeps only what all
    its members share.
    """
    groups: list[list] = []          # [core, member count]
    for seq in primers:
        if not seq:
            continue
        for g in groups:
            core = shared_region(seq, g[0])
            if core:
                g[0] = core
                g[1] += 1
                break
        else:
            groups.append([seq, 1])
    # Merging shrinks a core, which can bring two groups within reach of each
    # other, so keep merging until nothing more comes together.
    merged = True
    while merged:
        merged = False
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                core = shared_region(groups[i][0], groups[j][0])
                if core:
                    groups[i][0] = core
                    groups[i][1] += groups[j][1]
                    del groups[j]
                    merged = True
                    break
            if merged:
                break
    return sorted(((g[0], g[1]) for g in groups), key=lambda t: -t[1])


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


_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}


def reverse_complement(seq: str) -> str:
    """Reverse complement, preserving IUPAC degeneracy."""
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


# A placement has to beat what the primer's own degeneracy would reach by
# chance. Identity alone cannot say that: scoring counts a match whenever the
# reference base falls inside the primer's IUPAC set, so an N matches
# everywhere and a B three times in four. nifH's PolF is 28 bases of which ten
# are degenerate, and the best of the ~6,000 offsets across both references and
# both strands reaches 0.82 on yeast 18S — enough to clear any fixed floor, on a
# gene that is not ribosomal at all.
#
# So compare the match to the primer's own null: each position matches by chance
# with probability |IUPAC set| / 4, which gives a mean and variance for the
# number of matching positions, and the placement is scored in standard
# deviations above that mean.
#
# Where the bar sits is decided by the two ends of the observed gap. Below it,
# the highest-scoring thing that is not a ribosomal primer: the truncated nifH
# forward at 4.33, with COI around 4.1. Above it, the weakest placement the
# catalogue still depends on: LABY-A's reverse end at 5.34, the only one of 259
# catalogue placements anywhere near the bar. Everything else that is genuinely
# SSU sits at 6.7-8.0.
#
# 12S is the case this cannot decide. Mitochondrial 12S is an SSU rRNA, so the
# teleo forward reads as one (0.88 on E. coli, z 6.0) and is right to. Telling
# it from 16S needs a 12S reference, not a threshold (#51).
_MIN_PLACEMENT_Z = 5.0


def _placement_z(primer: str, matched: int) -> float:
    """How many standard deviations `matched` positions beat chance for `primer`."""
    ps = [len(_SETS.get(b, _SETS["N"])) / 4 for b in primer]
    mu = sum(ps)
    var = sum(p * (1 - p) for p in ps)
    if var <= 0:
        # Every position degenerate: the primer matches everywhere and says
        # nothing about where it is.
        return 0.0
    return (matched - mu) / var ** 0.5


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
    database. A match counts only when it beats the primer's own degeneracy by
    _MIN_PLACEMENT_Z, so a degenerate primer for some other gene reports nothing
    rather than the best of several thousand accidents.
    """
    if not primer:
        return []
    length = len(primer)
    # A reverse primer is written as the reverse complement of the strand the
    # reference is on, so searching it as given can never match. Try both and
    # report which orientation placed it.
    forms = (("+", primer), ("-", reverse_complement(primer)))
    out = []
    for ref_name, ref in _ssu_references().items():
        best = (0.0, None, "+")
        for strand, seq in forms:
            for i in range(len(ref) - length + 1):
                ok = sum(1 for k in range(length)
                         if ref[i + k] in _SETS.get(seq[k], _SETS["N"]))
                score = ok / length
                if score > best[0]:
                    best = (score, i + 1, strand)
        z = _placement_z(primer, best[0] * length)
        if best[0] >= min_id and z >= _MIN_PLACEMENT_Z:
            out.append({"reference": ref_name, "identity": round(best[0], 3),
                        "start": best[1], "end": best[1] + length - 1,
                        "strand": best[2], "z": round(z, 2)})
    return sorted(out, key=lambda d: -d["identity"])


# How close a second reference has to score before a primer is called universal.
# Judged on what the sequence does rather than on what the primer is
# conventionally called: 1492R is spoken of as universal and is 1.0 on E. coli
# against 0.909 on yeast, which is a preference, while 515F is 1.0 on both and
# genuinely cannot tell the domains apart.
_UNIVERSAL_MARGIN = 0.05


def is_universal(locations: list[dict]) -> bool:
    """True when a primer places as well on a second reference as on the first.

    Worth knowing separately from where it placed. A universal primer still has
    to have its coordinates quoted against one reference, and which one wins is
    then decided by a tie — so calling the assay bacterial because E. coli was
    reached first states a domain the sequence never established.
    """
    if not locations or len(locations) < 2:
        return False
    return (locations[0]["identity"] - locations[1]["identity"]) <= _UNIVERSAL_MARGIN

@lru_cache(maxsize=1)
def _catalogue_spans() -> tuple:
    """Every catalogue primer's placement on the SSU references.

    (reference, start, end, name), one row per primer per reference it lands on.
    """
    rows, seen = [], set()
    for p in PRIMER_DB:
        for end in ("fwd", "rev"):
            seq = p[end]
            name = p["name"] if end == "fwd" else p["rev_name"]
            if not seq or seq in seen:
                continue
            seen.add(seq)
            for loc in locate_on_ssu(seq):
                rows.append((loc["reference"], loc["start"], loc["end"], name or seq))
    return tuple(rows)


# Where the variable regions sit on E. coli 16S, in the J01695 numbering every
# 16S primer name already refers to — 515F is the forward primer at position 515.
# Checked against the catalogue: of 108 amplicons whose label names a region, 102
# do contain the interval below. All six exceptions are the label being wrong
# rather than the coordinates — k1f/k1r is filed as V3-V4 and spans 515-805,
# which is V4 alone.
_VARIABLE_REGIONS_16S = {
    "V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682),
    "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173), "V8": (1243, 1294),
    "V9": (1435, 1465),
}

# How much of a variable region an amplicon has to span before it is said to
# cover it. Not all of it: a primer sitting a few bases inside a boundary still
# reads that region, and the boundaries themselves are a convention drawn on a
# continuum.
_REGION_COVER = 0.9


@lru_cache(maxsize=1)
def _reference_ladder() -> tuple:
    """Positions that are the same place on both references, ascending.

    The SSU is one molecule across the domains, so a primer that matches both
    references is matching the same site twice and its two coordinates are a
    correspondence. The catalogue supplies these for free: the universal primers
    in it place on both, closely enough together to interpolate between.
    """
    by_16s: dict[int, list[int]] = {}
    seen = set()
    for p in PRIMER_DB:
        for end in ("fwd", "rev"):
            seq = p.get(end)
            if not seq or seq in seen:
                continue
            seen.add(seq)
            locs = {loc["reference"]: loc for loc in locate_on_ssu(seq)}
            if "ecoli_16S" in locs and "yeast_18S" in locs:
                by_16s.setdefault(locs["ecoli_16S"]["start"], []).append(
                    locs["yeast_18S"]["start"])
    # One 18S position per 16S position, so the ladder cannot double back.
    return tuple(sorted((a, sorted(v)[len(v) // 2]) for a, v in by_16s.items()))


def _to_18s(pos: int) -> int:
    """An E. coli 16S coordinate at the corresponding place on yeast 18S."""
    ladder = _reference_ladder()
    if not ladder:
        return pos
    xs = [a for a, _ in ladder]
    ys = [b for _, b in ladder]
    if pos <= xs[0]:
        return ys[0] + (pos - xs[0])
    if pos >= xs[-1]:
        return ys[-1] + (pos - xs[-1])
    i = bisect.bisect_right(xs, pos) - 1
    x0, x1, y0, y1 = xs[i], xs[i + 1], ys[i], ys[i + 1]
    return round(y0 + (pos - x0) * (y1 - y0) / (x1 - x0))


@lru_cache(maxsize=1)
def _variable_regions() -> dict:
    """{reference: {region: (start, end)}} in each reference's own numbering.

    18S is carried over from 16S through the ladder rather than tabulated
    separately: the regions are the same features of the same molecule, and the
    eukaryotic expansions are exactly what the interpolation accounts for — V4
    is 107 bases on E. coli and 213 on yeast. Checked the same way as the 16S
    table, at 87 of 93.
    """
    return {
        "ecoli_16S": dict(_VARIABLE_REGIONS_16S),
        "yeast_18S": {r: (_to_18s(a), _to_18s(b))
                      for r, (a, b) in _VARIABLE_REGIONS_16S.items()},
    }


def regions_in_span(reference: str, start: int, end: int) -> str | None:
    """Which variable regions an amplicon covers: "V4", "V3-V4", or None.

    This is what a placement is for. An assay whose primers the catalogue does
    not hold still lands somewhere on the gene, and where it lands is the same
    fact the catalogue would have supplied — so a region can be named for an
    assay nobody has catalogued (danaSeq #53).
    """
    table = _variable_regions().get(reference)
    if not table or not start or not end or end <= start:
        return None
    covered = []
    for name, (a, b) in table.items():
        overlap = min(end, b) - max(start, a) + 1
        if overlap >= _REGION_COVER * (b - a + 1):
            covered.append((int(name[1:]), name))
    if not covered:
        return None
    covered.sort()
    # Contiguous or not, the span is what it reaches: an amplicon covering V3
    # and V4 is V3-V4, and one reaching V3 and V5 is still described by its ends.
    return covered[0][1] if len(covered) == 1 else f"{covered[0][1]}-{covered[-1][1]}"

def amplicon_boundary(consensus: str, slack: int = 1) -> dict:
    """Whether a 5' block is a primer, or the cut left where one was removed.

    A primer occupies its own coordinates: a block that begins where catalogue
    primers begin is one. A block that begins in the base immediately after where
    they end is not a primer at all — it is the amplicon, delivered with the
    primer already taken off, and trimming it removes real sequence.

    Which primer was removed is not recoverable and is not guessed: ten
    catalogued primers end against E. coli 534, and nothing in the reads says
    which of them made the cut. The boundary is the observation, so the boundary
    is what is reported.

    Measure this on a sample's own consensus, not on a set collapsed across
    samples: collapsing keeps only what the members share and moves the edge off
    the boundary, which is the one coordinate the question turns on.

    Returns state "primer", "pre-trimmed", or "unplaced" — the last for a block
    that sits on neither reference, which is every non-ribosomal assay, and where
    there is no gene to measure against and so nothing to conclude.
    """
    locs = locate_on_ssu(consensus)
    if not locs:
        return {"state": "unplaced", "location": []}
    spans = _catalogue_spans()
    starts_on, starts_after = [], []
    for loc in locs:
        # The amplicon's outer edge is the 5' end of the read, which is the low
        # coordinate for a forward block and the high one for a reverse block.
        forward = loc["strand"] == "+"
        edge = loc["start"] if forward else loc["end"]
        for ref, start, end, name in spans:
            if ref != loc["reference"]:
                continue
            own = start if forward else end
            adjacent = end + 1 if forward else start - 1
            if abs(own - edge) <= slack:
                starts_on.append(name)
            elif abs(adjacent - edge) <= slack:
                starts_after.append(name)
    if starts_on:
        state = "primer"
    elif starts_after:
        state = "pre-trimmed"
    else:
        # On the gene but at no catalogued boundary: a primer nobody has
        # catalogued is likelier than a cut nobody makes.
        state = "primer"
    return {"state": state, "location": locs,
            "abuts": sorted(set(starts_after)) if state == "pre-trimmed" else []}


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

    # Measured on the sample's own consensus, which still carries the edge the
    # boundary test needs — a set collapsed across samples no longer does.
    fwd_edge = amplicon_boundary(fwd)
    rev_edge = amplicon_boundary(rev) if rev else {"state": "unplaced", "location": []}
    fwd_loc = fwd_edge["location"]
    rev_loc = rev_edge["location"]
    # Both ends landing nowhere on either SSU gene means the assay is not
    # ribosomal, whatever the reads denoise into. Taxonomy against an rRNA
    # database would still return confident labels, and they would be noise.
    ribosomal = bool(fwd_loc or rev_loc)

    return {
        "fwd": fwd, "rev": rev,
        "fwd_name": fwd if fwd_src == "de-novo" else fwd_detail.get("name", fwd),
        "rev_name": rev if rev_src == "de-novo" else rev_detail.get("name", rev),
        "region": "unknown",
        # Per end, because the ends are resolved independently and routinely
        # disagree — a de-novo forward primer alongside a catalogue reverse one.
        # One source for the record can only carry one of those answers, and the
        # end that loses is then treated as de-novo, which strips a catalogue
        # sequence of the protection that keeps it whole through collapsing.
        "fwd_source": f"inferred-{fwd_src}",
        "rev_source": f"inferred-{rev_src}",
        # Each end's primer looked for in the other end's reads. Reads deposited
        # in a random orientation put both primers in both files at roughly equal
        # frequency; correctly oriented reads put each primer in one file only.
        "fwd_in_r2": round(read_support(r2, fwd), 3) if (r2 and fwd) else None,
        "rev_in_r1": round(read_support(r1, rev), 3) if rev else None,
        "fwd_state": fwd_edge["state"], "rev_state": rev_edge["state"],
        "fwd_abuts": fwd_edge.get("abuts", []),
        "rev_abuts": rev_edge.get("abuts", []),
        "confidence": fwd_detail.get("read_support"),
        "ribosomal": ribosomal,
        "fwd_detail": fwd_detail, "rev_detail": rev_detail,
        "fwd_location": fwd_loc, "rev_location": rev_loc,
    }
