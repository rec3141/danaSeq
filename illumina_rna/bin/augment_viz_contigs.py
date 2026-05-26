#!/usr/bin/env python3
"""
For each MAG reference + its contigs:
  - Compute 4-mer (tetranucleotide) frequencies, fourth-root normalized
  - Compute per-contig GC content
  - Sum RNA read counts per contig across all samples (from idxstats)
  - Run sendsketch.sh per MAG to get best-hit taxonomy
  - Compute one t-SNE over all contigs combined

Emits:
  contigs.json.gz   { contigs: [ {id, mag, length, gc, x, y, reads_by_sample, total_reads} ] }
  taxonomy.json     { <mag>: {lineage, top_hit, ani} }
"""
import argparse, gzip, json, os, re, subprocess, sys
from collections import defaultdict, OrderedDict

import numpy as np
from sklearn.manifold import TSNE

# Canonical 4-mers (collapse RC). 256 raw -> 136 canonical.
def revcomp(s):
    t = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(t[c] for c in reversed(s))

def build_canonical_kmers(k=4):
    seen, lst = set(), []
    bases = "ACGT"
    def gen(prefix):
        if len(prefix) == k:
            rc = revcomp(prefix)
            canon = prefix if prefix <= rc else rc
            if canon not in seen:
                seen.add(canon); lst.append(canon)
            return
        for b in bases: gen(prefix + b)
    gen("")
    return lst

CANON = build_canonical_kmers(4)
CANON_INDEX = {km: i for i, km in enumerate(CANON)}

def kmer_index(kmer):
    rc = revcomp(kmer)
    return CANON_INDEX.get(kmer if kmer <= rc else rc)

def tnf_vector(seq, k=4):
    """Fourth-root-normalized canonical 4-mer freq vector."""
    seq = seq.upper()
    counts = np.zeros(len(CANON), dtype=np.float32)
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        j = kmer_index(kmer)
        if j is not None:
            counts[j] += 1
    total = counts.sum()
    if total > 0:
        counts /= total
    return np.power(counts, 0.25)

def gc_content(seq):
    seq = seq.upper()
    gc = sum(1 for c in seq if c in "GC")
    at = sum(1 for c in seq if c in "AT")
    return gc / (gc + at) if (gc + at) else 0.0

def parse_fasta(path):
    """Yield (id, sequence) — strips trailing metadata after first whitespace."""
    name, parts = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            yield name, "".join(parts)

def parse_idxstats(path):
    """Return {contig: mapped_reads}."""
    out = {}
    if not os.path.exists(path): return out
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] != "*":
                try:
                    out[parts[0]] = int(parts[2])
                except ValueError: pass
    return out

def run_sendsketch(fasta, sendsketch_bin):
    """Query JGI sketch server; return (best_hit_name, ani) or (None, None).
    Format-3 columns: #Query Ref ANI QSize RefSize QBases RBases QTaxID RTaxID KID WKID SSU
    """
    try:
        r = subprocess.run(
            [sendsketch_bin, f"in={fasta}", "format=3", "records=1", "color=f"],
            capture_output=True, text=True, timeout=120,
        )
        for line in (r.stdout + "\n" + r.stderr).splitlines():
            if line.startswith("#") or not line.strip(): continue
            cols = line.split("\t")
            if len(cols) < 3: continue
            try:
                ani = float(cols[2])
            except ValueError: continue
            return cols[1].strip(), ani   # col 1 = Ref (e.g. "Colwellia sp. Arc7-D")
        return None, None
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError) as e:
        print(f"  sendsketch failed for {fasta}: {e}", file=sys.stderr)
        return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refs-dir",   required=True, help="references/ dir of <mag>.fasta")
    ap.add_argument("--mapping-dir", required=True, help="illumina_rna_out/mapping/")
    ap.add_argument("--viz-dir",    required=True)
    ap.add_argument("--sendsketch", default="/home/grid/apps/BBTools-39.52/sendsketch.sh")
    ap.add_argument("--skip-taxonomy", action="store_true")
    args = ap.parse_args()

    fastas = sorted(f for f in os.listdir(args.refs_dir) if f.endswith(".fasta"))
    mags = [os.path.splitext(f)[0] for f in fastas]

    # Sum idxstats across all samples per (mag, contig). The publishDir pattern
    # in alignment.nf only matches *.{bam,bai,covhist,covstats,flagstat,idxstats}
    # so the actual .idxstats.tsv files weren't published — but the BAMs were
    # symlinked, so we can follow the symlink target back to the work dir which
    # has the sibling idxstats.
    contig_reads = defaultdict(lambda: defaultdict(int))    # mag -> contig -> reads
    contig_reads_by_sample = defaultdict(lambda: defaultdict(dict))
    samples_seen = set()
    for sample_dir in sorted(os.listdir(args.mapping_dir)):
        sd = os.path.join(args.mapping_dir, sample_dir)
        if not os.path.isdir(sd): continue
        samples_seen.add(sample_dir)
        for mag in mags:
            tag = f"{sample_dir}_vs_{mag}"
            # First try the publishDir path:
            idx = os.path.join(sd, mag, f"{tag}.idxstats.tsv")
            if not os.path.exists(idx):
                # Fall back: resolve the published BAM symlink → its work dir → sibling
                bam_sym = os.path.join(sd, mag, f"{tag}.sorted.bam")
                if os.path.islink(bam_sym) or os.path.exists(bam_sym):
                    work_dir = os.path.dirname(os.path.realpath(bam_sym))
                    cand = os.path.join(work_dir, f"{tag}.idxstats.tsv")
                    if os.path.exists(cand):
                        idx = cand
            for contig, n in parse_idxstats(idx).items():
                contig_reads[mag][contig] += n
                contig_reads_by_sample[mag][contig][sample_dir] = n
    samples = sorted(samples_seen)
    print(f"[contigs] {len(mags)} MAGs, {len(samples)} samples", file=sys.stderr)

    # Compute TNF + GC per contig
    all_tnf, meta = [], []
    for mag, fa in zip(mags, fastas):
        path = os.path.join(args.refs_dir, fa)
        n_c = 0
        for cid, seq in parse_fasta(path):
            if len(seq) < 1000:   # skip tiny contigs
                continue
            all_tnf.append(tnf_vector(seq))
            meta.append({
                "id": cid,
                "mag": mag,
                "length": len(seq),
                "gc": gc_content(seq),
            })
            n_c += 1
        print(f"  {mag}: {n_c} contigs ≥1kb kept", file=sys.stderr)

    if not all_tnf:
        print("[contigs] no contigs to embed — aborting", file=sys.stderr)
        sys.exit(1)

    X = np.array(all_tnf, dtype=np.float32)
    print(f"[contigs] TNF matrix: {X.shape}", file=sys.stderr)

    perp = max(5, min(30, X.shape[0] // 20))
    print(f"[contigs] running t-SNE (perplexity={perp})...", file=sys.stderr)
    tsne = TSNE(n_components=2, perplexity=perp, learning_rate="auto",
                init="pca", max_iter=500, random_state=42)
    coords = tsne.fit_transform(X)

    contigs_out = []
    for i, m in enumerate(meta):
        cr = contig_reads.get(m["mag"], {}).get(m["id"], 0)
        rbs = contig_reads_by_sample.get(m["mag"], {}).get(m["id"], {})
        contigs_out.append({
            **m,
            "x": float(coords[i, 0]),
            "y": float(coords[i, 1]),
            "total_reads": cr,
            "reads_by_sample": [rbs.get(s, 0) for s in samples],
        })

    with gzip.open(os.path.join(args.viz_dir, "contigs.json.gz"), "wt") as f:
        json.dump({"samples": samples, "contigs": contigs_out}, f)
    print(f"[contigs] wrote contigs.json.gz ({len(contigs_out)} contigs)", file=sys.stderr)

    # Taxonomy via sendsketch
    if args.skip_taxonomy:
        print("[contigs] skipping taxonomy", file=sys.stderr)
    else:
        tax_out = {}
        for mag, fa in zip(mags, fastas):
            path = os.path.join(args.refs_dir, fa)
            print(f"  sendsketch {mag}", file=sys.stderr)
            lineage, ani = run_sendsketch(path, args.sendsketch)
            tax_out[mag] = {"lineage": lineage, "ani": ani}
        with open(os.path.join(args.viz_dir, "taxonomy.json"), "w") as f:
            json.dump(tax_out, f)
        print(f"[contigs] wrote taxonomy.json ({sum(1 for v in tax_out.values() if v['lineage'])} hits)", file=sys.stderr)


if __name__ == "__main__":
    main()
