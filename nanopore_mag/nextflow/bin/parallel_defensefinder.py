#!/usr/bin/env python3
"""parallel_defensefinder.py — Parallel DefenseFinder wrapper

Splits a protein FASTA by contig (using GFF3 mapping), rewrites headers for
MacSyFinder's gembase format, runs N defense-finder instances in parallel,
and merges the output TSVs with IDs mapped back to original Bakta locus tags.

This works around the single-threaded system detection bottleneck in MacSyFinder:
in ordered_replicon mode, all proteins are treated as one giant replicon, so the
nested loop (replicons × models) runs in serial. With gembase + contig splitting,
each contig is an independent replicon and chunks run in parallel.

Usage:
    parallel_defensefinder.py --proteins FILE --gff FILE --output-dir DIR
        --workers N [--models-dir DIR]
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Parallel DefenseFinder wrapper")
    p.add_argument("--proteins", required=True, help="Protein FASTA (.faa)")
    p.add_argument("--gff", required=True, help="GFF3 annotation file")
    p.add_argument("--output-dir", required=True, help="Output directory")
    p.add_argument("--workers", type=int, default=4, help="Number of parallel chunks")
    p.add_argument("--models-dir", default=None, help="DefenseFinder models directory")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Step 1: Parse GFF3 → protein-to-contig mapping
# ---------------------------------------------------------------------------

def parse_gff_protein_map(gff_path):
    """Parse GFF3 to build {protein_id: contig_name} for CDS features."""
    prot2contig = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue
            contig = cols[0]
            m = re.search(r"ID=([^;]+)", cols[8])
            if m:
                prot2contig[m.group(1)] = contig
    return prot2contig


# ---------------------------------------------------------------------------
# Step 2: Rewrite FASTA headers for gembase format, group by contig
# ---------------------------------------------------------------------------

def rewrite_and_group_fasta(fasta_path, prot2contig):
    """Read protein FASTA, rewrite headers for gembase, group by contig.

    Gembase format: >REPLICON_PROTEINID (MacSyFinder splits on last '_')
    We strip underscores from the protein ID so the last '_' cleanly separates
    the contig (replicon) from the protein identifier.

    Returns:
        by_contig: {contig: [(gembase_header_line, seq_text), ...]}
        rev_map:   {gembase_id: original_id}
    """
    by_contig = defaultdict(list)
    rev_map = {}

    def flush(hdr, seq_parts):
        if hdr is None:
            return
        tokens = hdr[1:].split(None, 1)
        orig_id = tokens[0]
        desc = tokens[1] if len(tokens) > 1 else ""

        contig = prot2contig.get(orig_id, "unknown")
        if contig == "unknown":
            print(
                f"[WARNING] Protein {orig_id} not found in GFF "
                f"— assigning to 'unknown' replicon",
                file=sys.stderr,
            )

        # Remove underscores from protein ID for gembase compatibility
        clean = orig_id.replace("_", "")
        gembase_id = f"{contig}_{clean}"
        rev_map[gembase_id] = orig_id
        rev_map[clean] = orig_id  # in case some fields use protein part only

        new_hdr = f">{gembase_id}"
        if desc:
            new_hdr += f" {desc}"
        by_contig[contig].append((new_hdr + "\n", "".join(seq_parts)))

    hdr = None
    seq = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush(hdr, seq)
                hdr = line.rstrip("\n")
                seq = []
            else:
                seq.append(line)
    flush(hdr, seq)

    return by_contig, rev_map


# ---------------------------------------------------------------------------
# Step 3: Split contigs into N balanced chunks (greedy bin packing)
# ---------------------------------------------------------------------------

def split_chunks(by_contig, n):
    """Distribute contigs across N chunks, balancing protein count."""
    items = sorted(by_contig.items(), key=lambda kv: -len(kv[1]))
    chunks = [dict() for _ in range(n)]
    totals = [0] * n
    for contig, recs in items:
        idx = totals.index(min(totals))
        chunks[idx][contig] = recs
        totals[idx] += len(recs)
    return chunks


def write_chunk(chunk, path):
    """Write chunk records to a FASTA file, preserving within-contig order."""
    with open(path, "w") as fh:
        for contig in sorted(chunk):
            for hdr_line, seq_str in chunk[contig]:
                fh.write(hdr_line)
                fh.write(seq_str)


# ---------------------------------------------------------------------------
# Step 4: Run defense-finder on each chunk
# ---------------------------------------------------------------------------

def run_chunk(chunk_fasta, out_dir, models_dir):
    """Run defense-finder on one chunk. Returns (out_dir, returncode, stderr)."""
    cmd = [
        "defense-finder", "run",
        "--db-type", "gembase",
        "-o", str(out_dir),
        "-w", "1",
    ]
    if models_dir:
        cmd += ["--models-dir", str(models_dir)]
    cmd.append(str(chunk_fasta))
    r = subprocess.run(cmd, capture_output=True, text=True)
    return (str(out_dir), r.returncode, r.stderr)


# ---------------------------------------------------------------------------
# Step 5: Merge per-chunk outputs
# ---------------------------------------------------------------------------

def find_tsv(directory, pattern):
    """Glob for a single file matching pattern in directory."""
    hits = list(Path(directory).glob(pattern))
    return hits[0] if hits else None


def col_index(header_line, name):
    """Return 0-based column index for a header name, or -1 if missing."""
    cols = header_line.rstrip("\n").split("\t")
    try:
        return cols.index(name)
    except ValueError:
        return -1


def map_ids(val, rev):
    """Map comma-separated gembase IDs back to originals."""
    if not val:
        return val
    return ",".join(rev.get(v.strip(), v.strip()) for v in val.split(","))


def merge_all(chunk_dirs, output_dir, rev_map):
    """Merge per-chunk DefenseFinder outputs into final TSVs."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Shared sys_id renumbering state (closures share these)
    sys_id_map = {}
    counter = [0]

    def renumber(old_id):
        """Renumber a sys_id globally across chunks."""
        if old_id in sys_id_map:
            return sys_id_map[old_id]
        counter[0] += 1
        # Preserve type info: "contig_1_AbiD_3" → "contig_1_AbiD_{global_N}"
        m = re.match(r"^(.+)_(\d+)$", old_id)
        if m:
            new_id = f"{m.group(1)}_{counter[0]}"
        else:
            new_id = f"{old_id}_{counter[0]}"
        sys_id_map[old_id] = new_id
        return new_id

    def map_sys_ids(val):
        """Map comma-separated sys_ids through renumbering."""
        if not val:
            return val
        return ",".join(
            renumber(v.strip()) if v.strip() else v for v in val.split(",")
        )

    # Merge in order: systems first (builds sys_id_map), then genes (uses it)
    _merge_systems(chunk_dirs, out, rev_map, renumber)
    _merge_genes(chunk_dirs, out, rev_map, map_sys_ids)
    _merge_hmmer(chunk_dirs, out, rev_map)


# Default headers for empty output files
SYSTEMS_HEADER = (
    "sys_id\ttype\tsubtype\tactivity\tsys_beg\tsys_end\t"
    "protein_in_syst\tgenes_count\tname_of_profiles_in_sys\n"
)
GENES_HEADER = (
    "replicon\thit_id\tgene_name\thit_pos\tmodel_fqn\tsys_id\t"
    "sys_loci\tlocus_num\tsys_wholeness\tsys_score\tsys_occ\t"
    "hit_gene_ref\thit_status\thit_seq_len\thit_i_eval\thit_score\t"
    "hit_profile_cov\thit_seq_cov\thit_begin_match\thit_end_match\t"
    "counterpart\tused_in\ttype\tsubtype\tactivity\n"
)
HMMER_HEADER = (
    "hit_id\treplicon\thit_pos\thit_sequence_length\tgene_name\t"
    "i_eval\thit_score\thit_profile_cov\thit_seq_cov\t"
    "hit_begin_match\thit_end_match\n"
)


def _merge_systems(chunk_dirs, out, rev_map, renumber_fn):
    path = out / "systems.tsv"
    header = None
    ci = {}

    with open(path, "w") as fh:
        for d in chunk_dirs:
            f = find_tsv(d, "*_defense_finder_systems.tsv")
            if not f:
                continue
            with open(f) as inp:
                for i, line in enumerate(inp):
                    if i == 0:
                        if header is None:
                            header = line
                            fh.write(line)
                            for name in ("sys_id", "sys_beg", "sys_end",
                                         "protein_in_syst"):
                                ci[name] = col_index(line, name)
                        continue
                    if not line.strip():
                        continue
                    cols = line.rstrip("\n").split("\t")

                    # Renumber sys_id
                    idx = ci.get("sys_id", -1)
                    if 0 <= idx < len(cols):
                        cols[idx] = renumber_fn(cols[idx])

                    # Map protein IDs back to originals
                    for name in ("sys_beg", "sys_end", "protein_in_syst"):
                        idx = ci.get(name, -1)
                        if 0 <= idx < len(cols):
                            cols[idx] = map_ids(cols[idx], rev_map)

                    fh.write("\t".join(cols) + "\n")

    if header is None:
        with open(path, "w") as fh:
            fh.write(SYSTEMS_HEADER)


def _merge_genes(chunk_dirs, out, rev_map, map_sys_ids_fn):
    path = out / "genes.tsv"
    header = None
    ci = {}

    with open(path, "w") as fh:
        for d in chunk_dirs:
            f = find_tsv(d, "*_defense_finder_genes.tsv")
            if not f:
                continue
            with open(f) as inp:
                for i, line in enumerate(inp):
                    if i == 0:
                        if header is None:
                            header = line
                            fh.write(line)
                            for name in ("hit_id", "sys_id", "counterpart",
                                         "used_in"):
                                ci[name] = col_index(line, name)
                        continue
                    if not line.strip():
                        continue
                    cols = line.rstrip("\n").split("\t")

                    # Map hit_id back to original
                    idx = ci.get("hit_id", -1)
                    if 0 <= idx < len(cols):
                        cols[idx] = rev_map.get(cols[idx], cols[idx])

                    # Renumber sys_id
                    idx = ci.get("sys_id", -1)
                    if 0 <= idx < len(cols):
                        cols[idx] = map_sys_ids_fn(cols[idx])

                    # Renumber sys_id references in counterpart and used_in
                    for name in ("counterpart", "used_in"):
                        idx = ci.get(name, -1)
                        if 0 <= idx < len(cols):
                            cols[idx] = map_sys_ids_fn(cols[idx])

                    fh.write("\t".join(cols) + "\n")

    if header is None:
        with open(path, "w") as fh:
            fh.write(GENES_HEADER)


def _merge_hmmer(chunk_dirs, out, rev_map):
    path = out / "hmmer.tsv"
    header = None
    hit_id_col = -1

    with open(path, "w") as fh:
        for d in chunk_dirs:
            f = find_tsv(d, "*_defense_finder_hmmer.tsv")
            if not f:
                continue
            with open(f) as inp:
                for i, line in enumerate(inp):
                    if i == 0:
                        if header is None:
                            header = line
                            fh.write(line)
                            hit_id_col = col_index(line, "hit_id")
                        continue
                    if not line.strip():
                        continue
                    cols = line.rstrip("\n").split("\t")

                    if 0 <= hit_id_col < len(cols):
                        cols[hit_id_col] = rev_map.get(
                            cols[hit_id_col], cols[hit_id_col]
                        )

                    fh.write("\t".join(cols) + "\n")

    if header is None:
        with open(path, "w") as fh:
            fh.write(HMMER_HEADER)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    proteins = Path(args.proteins)
    gff = Path(args.gff)
    out_dir = Path(args.output_dir)
    n_workers = max(1, args.workers)

    # --- Handle empty input ---
    if not proteins.exists() or proteins.stat().st_size == 0:
        print("[WARNING] Empty protein FASTA — writing empty output TSVs",
              file=sys.stderr)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "systems.tsv").write_text(SYSTEMS_HEADER)
        (out_dir / "genes.tsv").write_text(GENES_HEADER)
        (out_dir / "hmmer.tsv").write_text(HMMER_HEADER)
        return

    # --- Step 1: Parse GFF → protein-to-contig map ---
    print(f"[INFO] Parsing GFF: {gff}", file=sys.stderr)
    prot2contig = parse_gff_protein_map(str(gff))
    print(f"[INFO] Found {len(prot2contig)} CDS entries in GFF", file=sys.stderr)

    # --- Step 2: Rewrite FASTA headers and group by contig ---
    print("[INFO] Rewriting FASTA headers for gembase format", file=sys.stderr)
    by_contig, rev_map = rewrite_and_group_fasta(str(proteins), prot2contig)
    total_prots = sum(len(v) for v in by_contig.values())
    print(f"[INFO] {total_prots} proteins across {len(by_contig)} contigs",
          file=sys.stderr)

    # --- Step 3: Split into balanced chunks ---
    n_chunks = min(n_workers, len(by_contig))
    if n_chunks < 1:
        n_chunks = 1
    print(f"[INFO] Splitting into {n_chunks} chunks for {n_workers} workers",
          file=sys.stderr)
    chunks = split_chunks(by_contig, n_chunks)

    # Write chunk FASTAs to temp directory
    tmp_dir = Path(tempfile.mkdtemp(prefix="df_parallel_",
                                     dir=str(out_dir.parent)))
    chunk_fastas = []
    chunk_outdirs = []
    for i, chunk in enumerate(chunks):
        if not chunk:
            continue
        fa = tmp_dir / f"chunk_{i}.faa"
        write_chunk(chunk, fa)
        chunk_fastas.append(fa)
        chunk_outdirs.append(tmp_dir / f"chunk_{i}_out")

    # --- Step 4: Run in parallel ---
    print(
        f"[INFO] Running {len(chunk_fastas)} defense-finder instances in parallel",
        file=sys.stderr,
    )

    try:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {
                pool.submit(run_chunk, fa, od, args.models_dir): i
                for i, (fa, od) in enumerate(zip(chunk_fastas, chunk_outdirs))
            }
            for fut in as_completed(futures):
                idx = futures[fut]
                od, rc, stderr = fut.result()
                if rc != 0:
                    print(
                        f"[WARNING] Chunk {idx} defense-finder exited {rc}:\n{stderr}",
                        file=sys.stderr,
                    )
                else:
                    print(f"[INFO] Chunk {idx} completed successfully",
                          file=sys.stderr)

        # --- Step 5: Merge outputs ---
        valid_dirs = [str(d) for d in chunk_outdirs if d.exists()]
        print(f"[INFO] Merging outputs from {len(valid_dirs)} chunks",
              file=sys.stderr)
        merge_all(valid_dirs, str(out_dir), rev_map)

        # Report results
        sys_file = out_dir / "systems.tsv"
        if sys_file.exists():
            n_systems = sum(1 for _ in open(sys_file)) - 1
            print(f"[INFO] Total defense systems found: {n_systems}",
                  file=sys.stderr)

    finally:
        # --- Step 6: Cleanup temp files ---
        shutil.rmtree(tmp_dir, ignore_errors=True)

    print(f"[INFO] Done. Output in {out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
