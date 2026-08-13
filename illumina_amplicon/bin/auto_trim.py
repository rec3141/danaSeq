#!/usr/bin/env python3
"""Quality profiling and automated truncation length selection (Python engine).

Vendored from the microscape package (microscape/quality.py).

Must stay in step with bin/dada2_auto_trim.R, the R engine's implementation —
the two are compared byte-for-byte on real and synthetic libraries.
"""
import gzip
import glob
import os
from typing import List, Optional, Tuple

import numpy as np


def _read_quals(files: List[str], n_reads: int) -> List[List[int]]:
    """Read quality scores from FASTQ files."""
    all_quals = []
    count = 0
    for fpath in files:
        opener = gzip.open if fpath.endswith(".gz") else open
        with opener(fpath, "rt") as fh:
            while count < n_reads:
                header = fh.readline()
                if not header:
                    break
                fh.readline()  # seq
                fh.readline()  # +
                qual_str = fh.readline().rstrip()
                if not qual_str:
                    break
                all_quals.append([ord(c) - 33 for c in qual_str])
                count += 1
        if count >= n_reads:
            break
    return all_quals


def _amplicon_lengths(lengths: List[int], floor_frac: float = 0.5) -> List[int]:
    """Drop fragment reads before measuring the length distribution.

    Amplicon libraries routinely carry adapter-dimer and primer-only fragments
    an order of magnitude shorter than the product. Those are real reads but not
    reads *of the amplicon*, and they make the length distribution bimodal —
    which is what breaks a low percentile.

    Measured on a real 18S library: 325 of 1,426 reads (22.8%) fell under 50 bp
    against a modal length of 232, so the 10th percentile landed at 31 bp, deep
    inside the fragment mode. Truncating there left 150 bp per read against a
    277 bp amplicon, and 99.5% of pairs then failed to merge — the sample
    vanished from the final table.

    Keeping only reads of at least `floor_frac` of the median restores a
    unimodal distribution, so the percentile measures what it is meant to. A
    read under half the typical length cannot carry the expected product anyway,
    and filterAndTrim would discard it regardless.
    """
    if not lengths:
        return lengths
    median = float(np.median(lengths))
    usable = [n for n in lengths if n >= median * floor_frac]
    # If nearly everything is a fragment the library is broken rather than
    # bimodal, and measuring the remainder would be less honest than measuring
    # all of it.
    return usable if len(usable) >= max(10, 0.05 * len(lengths)) else lengths


def _find_trunc_pos(quals_list: List[List[int]], min_q: float,
                    window: int, min_length: int = 0) -> int:
    """Find position where rolling median quality drops below threshold.

    Returns the quality-driven position, floored at *min_length*, then capped
    at the 10th percentile read length so short reads aren't rejected — while a
    population of short reads can't drag the cap down and truncate the
    full-length reads into oblivion.

    *min_length* floors the QUALITY position and nothing else. It exists so a
    quality dip at ~30bp cannot truncate reads below the length they need to
    overlap (issue #4). Applying it to the length cap instead asks for reads
    that were never sequenced: dada2 discards reads shorter than truncLen, so a
    floor above the library's real length silently returns zero reads for every
    sample. A 2x150 run whose amplicon reads are 126bp after primer removal was
    truncated at 148 and kept 607 of 101194 reads.
    """
    if not quals_list:
        return 0
    lengths = [len(q) for q in quals_list]
    max_len = max(lengths)
    # Cap truncation at the 10th-percentile read length so we don't reject more
    # than the shortest reads for being too short — but measure that percentile
    # over amplicon-length reads only. Fragments form a separate mode, and once
    # they exceed 10% of the library the percentile falls *into* that mode and
    # truncates the real reads below the length they need to overlap (see
    # _amplicon_lengths).
    len_cap = min(int(np.percentile(_amplicon_lengths(lengths), 10)), max_len)
    mat = np.full((len(quals_list), max_len), np.nan)
    for i, q in enumerate(quals_list):
        mat[i, :len(q)] = q
    medians = np.nanmedian(mat, axis=0)
    if len(medians) < window:
        return min(max(int(len(medians)), int(min_length)), len_cap)
    rolling = np.convolve(medians, np.ones(window) / window, mode="valid")
    quality_pos = int(len(medians))
    for i, val in enumerate(rolling):
        if val < min_q:
            quality_pos = i + window // 2
            break
    return min(max(quality_pos, int(min_length)), len_cap)


def auto_trim(
    input_dir: str,
    *,
    min_quality: float = 25.0,
    window: int = 10,
    n_reads: int = 10000,
    n_files: int = 20,
    min_length: int = 0,
    verbose: bool = False,
) -> dict:
    """Analyze quality profiles and recommend truncation lengths.

    Samples reads from paired-end FASTQs in *input_dir*, computes
    per-position median quality, and finds where a rolling median drops
    below *min_quality*.

    Args:
        input_dir: Directory containing paired FASTQ files.
        min_quality: Minimum rolling median Q score to keep a position.
        window: Rolling window size for quality smoothing.
        n_reads: Number of reads to sample.
        n_files: Maximum number of file pairs to sample from.
        verbose: Print progress messages.

    Returns:
        Dict with keys: trunc_len_fwd, trunc_len_rev, fwd_read_len,
        rev_read_len, n_reads_sampled, min_quality.
    """
    fwd_files = sorted(
        glob.glob(os.path.join(input_dir, "*_1.fastq.gz")) +
        glob.glob(os.path.join(input_dir, "*_R1*.fastq.gz")) +
        glob.glob(os.path.join(input_dir, "*_R1.*.fastq.gz"))
    )
    rev_files = sorted(
        glob.glob(os.path.join(input_dir, "*_2.fastq.gz")) +
        glob.glob(os.path.join(input_dir, "*_R2*.fastq.gz")) +
        glob.glob(os.path.join(input_dir, "*_R2.*.fastq.gz"))
    )
    # Deduplicate (patterns may overlap)
    fwd_files = sorted(set(fwd_files))
    rev_files = sorted(set(rev_files))

    if not fwd_files or not rev_files:
        raise FileNotFoundError(
            f"No paired FASTQ files found in {input_dir}"
        )

    # Sample from a subset of files
    n_sample = min(n_files, len(fwd_files))
    step = max(1, len(fwd_files) // n_sample)
    fwd_sample = fwd_files[::step][:n_sample]
    rev_sample = rev_files[::step][:n_sample]

    if verbose:
        print(f"[INFO] Sampling {n_reads} reads from {n_sample} file pairs")

    fwd_quals = _read_quals(fwd_sample, n_reads)
    rev_quals = _read_quals(rev_sample, n_reads)

    trunc_fwd = _find_trunc_pos(fwd_quals, min_quality, window, min_length)
    trunc_rev = _find_trunc_pos(rev_quals, min_quality, window, min_length)

    fwd_len = max(len(q) for q in fwd_quals) if fwd_quals else 0
    rev_len = max(len(q) for q in rev_quals) if rev_quals else 0

    if verbose:
        fwd_end_q = float(np.median([q[-1] for q in fwd_quals if q])) if fwd_quals else 0
        rev_end_q = float(np.median([q[-1] for q in rev_quals if q])) if rev_quals else 0
        print(f"[INFO] Forward: len={fwd_len}, Q@end={fwd_end_q:.0f}, trunc={trunc_fwd}")
        print(f"[INFO] Reverse: len={rev_len}, Q@end={rev_end_q:.0f}, trunc={trunc_rev}")

    return {
        "trunc_len_fwd": trunc_fwd,
        "trunc_len_rev": trunc_rev,
        "fwd_read_len": fwd_len,
        "rev_read_len": rev_len,
        "n_reads_sampled": len(fwd_quals),
        "min_quality": min_quality,
        "min_length": min_length,
    }


def main(argv=None):
    import argparse

    p = argparse.ArgumentParser(
        prog="auto_trim.py",
        description="Recommend truncation lengths from paired-end FASTQ quality profiles.")
    p.add_argument("input_dir")
    p.add_argument("--min-quality", type=float, default=25.0)
    p.add_argument("--window", type=int, default=10)
    p.add_argument("--n-reads", type=int, default=10000)
    p.add_argument("--n-files", type=int, default=20)
    p.add_argument("--min-length", type=int, default=0)
    p.add_argument("--output")
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args(argv)

    result = auto_trim(
        args.input_dir,
        min_quality=args.min_quality,
        window=args.window,
        n_reads=args.n_reads,
        n_files=args.n_files,
        min_length=args.min_length,
        verbose=args.verbose,
    )

    print(f"trunc_len_fwd={result['trunc_len_fwd']}")
    print(f"trunc_len_rev={result['trunc_len_rev']}")
    print(f"fwd_read_len={result['fwd_read_len']}")
    print(f"rev_read_len={result['rev_read_len']}")
    print(f"n_reads_sampled={result['n_reads_sampled']}")

    if args.output:
        with open(args.output, "w") as f:
            for key, val in result.items():
                f.write(f"{key}\t{val}\n")


if __name__ == "__main__":
    main()
