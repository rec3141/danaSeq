#!/usr/bin/env python3
"""
fastq_filter.py — Single-pass streaming FASTQ dedup + quality/length filter.

Replaces: cat *.fastq.gz | dedup | filtlong -t TARGET

Single pass, constant memory, pipe-friendly. Reads are scored and immediately
accepted or rejected. A running score distribution drives a dynamic threshold
that tightens as accepted_bases approaches target_bases.

Scoring follows filtlong's algorithm (without global z-score normalization):
  length_score  = 100 * len / (len + 5000)
  mean_quality  = mean(1 - 10^(-Q/10)) * 100
  window_quality = min sliding-window avg * 100
  final = sqrt(length_score * mean_quality) * (0.5 + 0.5 * window_ratio)

Memory: O(1) beyond the dedup set (~64 bytes per unique read ID).
"""

import sys
import gzip
import math
import argparse
import os
import shutil
import subprocess
# Precompute Phred+33 -> P(correct) lookup indexed by byte value (0-255)
# This avoids ord(c)-33 in the hot loop — index directly by byte value
PHRED_LUT = [0.0] * 256
for _q in range(94):
    PHRED_LUT[_q + 33] = 1.0 - 10.0 ** (-_q / 10.0)


def score_read_bytes(qual_bytes, length, window_size):
    """Score a read from raw quality bytes. Returns score.

    qual_bytes: bytes object (ASCII quality line, NOT decoded to str)
    length: int, read length (len of sequence)
    """
    if length == 0:
        return 0.0

    lut = PHRED_LUT
    length_score = 100.0 * length / (length + 5000.0)

    # Mean quality — direct byte indexing, no str decode needed
    total_q = 0.0
    for b in qual_bytes:
        total_q += lut[b]
    mean_q = total_q / length * 100.0

    if mean_q <= 0:
        return 0.0

    # Window quality: min sliding-window average
    if length <= window_size:
        window_q = mean_q
    else:
        wsum = 0.0
        for i in range(window_size):
            wsum += lut[qual_bytes[i]]
        min_wsum = wsum
        for i in range(1, length - window_size + 1):
            wsum += lut[qual_bytes[i + window_size - 1]] - lut[qual_bytes[i - 1]]
            if wsum < min_wsum:
                min_wsum = wsum
        window_q = min_wsum / window_size * 100.0

    window_ratio = min(window_q / mean_q, 1.0)
    return math.sqrt(length_score * mean_q) * (0.5 + 0.5 * window_ratio)


def estimate_total_bases(filenames):
    """Estimate total bases from gzipped file sizes. Returns 0 if unknown."""
    total = 0
    for f in filenames:
        if f == '-':
            return 0
        try:
            total += os.path.getsize(f)
        except OSError:
            return 0
    # gzipped FASTQ: ~4x compression, ~50% of uncompressed is sequence
    return int(total * 2)


PIGZ = shutil.which('pigz')


def open_fastq(fname):
    """Open a FASTQ file. Uses pigz for parallel gzip decompression if available."""
    if fname == '-':
        return sys.stdin.buffer
    elif fname.endswith('.gz'):
        if PIGZ:
            proc = subprocess.Popen(
                [PIGZ, '-dc', '-p', str(getattr(sys.modules[__name__], 'PIGZ_THREADS', 4)), fname],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            # Attach proc to the stdout so we can wait() on cleanup
            proc.stdout._proc = proc
            return proc.stdout
        return gzip.open(fname, 'rb')
    else:
        return open(fname, 'rb')


def main():
    parser = argparse.ArgumentParser(
        description='Single-pass streaming FASTQ dedup + quality/length filter')
    parser.add_argument('-t', '--target_bases', type=int, default=0,
                        help='Target bases to keep (0 = keep all, just dedup)')
    parser.add_argument('--min_length', type=int, default=0,
                        help='Minimum read length')
    parser.add_argument('--window_size', type=int, default=250,
                        help='Quality scoring window (default: 250)')
    parser.add_argument('--no_dedupe', action='store_true',
                        help='Skip deduplication')
    parser.add_argument('-o', '--output', default='-',
                        help='Output file (- for stdout, .gz for gzipped)')
    parser.add_argument('-p', '--threads', type=int, default=4,
                        help='Threads for pigz decompression (default: 4)')
    parser.add_argument('input', nargs='*', default=['-'],
                        help='Input FASTQ files (gzipped ok, - for stdin)')
    args = parser.parse_args()

    global PIGZ_THREADS
    PIGZ_THREADS = args.threads

    target = args.target_bases
    filtering = target > 0
    dedupe = not args.no_dedupe
    seen_ids = set() if dedupe else None
    window_size = args.window_size
    min_length = args.min_length

    # Estimate total bases for threshold calibration
    est_total = estimate_total_bases(args.input) if filtering else 0

    # Running score statistics (Welford's online algorithm)
    n_scored = 0
    score_mean = 0.0
    score_m2 = 0.0
    scored_bases = 0  # total bases seen (after dedup/minlen)
    threshold = 0.0  # accept everything until we have stats (first 100 reads)

    # Output
    if args.output == '-':
        out = sys.stdout.buffer
    elif args.output.endswith('.gz'):
        out = gzip.open(args.output, 'wb', compresslevel=3)
    else:
        out = open(args.output, 'wb')

    # Stats
    total_reads = 0
    dup_reads = 0
    short_reads = 0
    accepted_reads = 0
    accepted_bases = 0

    for fname in args.input:
        fh = open_fastq(fname)
        try:
            while True:
                header = fh.readline()
                if not header:
                    break
                seq_line = fh.readline()
                fh.readline()  # +
                qual_line = fh.readline()
                if not seq_line or not qual_line:
                    break

                total_reads += 1
                seq_raw = seq_line.rstrip(b'\n')
                read_len = len(seq_raw)
                qual_raw = qual_line.rstrip(b'\n')

                # Dedup by read ID (work in bytes, avoid decode)
                if dedupe:
                    hdr_raw = header.rstrip(b'\n')
                    # Extract read ID: skip '@', take first space-delimited field
                    sp = hdr_raw.find(b' ', 1)
                    rid = hdr_raw[1:sp] if sp > 0 else hdr_raw[1:]
                    if rid in seen_ids:
                        dup_reads += 1
                        continue
                    seen_ids.add(rid)
                else:
                    hdr_raw = header.rstrip(b'\n')

                # Min length
                if read_len < min_length:
                    short_reads += 1
                    continue

                # If no filtering, just output (all bytes, no decode)
                if not filtering:
                    out.write(hdr_raw)
                    out.write(b'\n')
                    out.write(seq_raw)
                    out.write(b'\n+\n')
                    out.write(qual_raw)
                    out.write(b'\n')
                    accepted_reads += 1
                    accepted_bases += read_len
                    continue

                score = score_read_bytes(qual_raw, read_len, window_size)

                # Update running statistics (Welford)
                n_scored += 1
                scored_bases += read_len
                delta = score - score_mean
                score_mean += delta / n_scored
                score_m2 += delta * (score - score_mean)

                # Update threshold every 1000 reads once we have stats
                if n_scored >= 100 and n_scored % 100 == 0:
                    score_std = math.sqrt(score_m2 / n_scored)
                    if score_std > 0:
                        # Remaining bases we still need
                        remaining_need = target - accepted_bases
                        # Remaining input bases (estimate from what we've seen)
                        if est_total > 0:
                            remaining_input = max(est_total - scored_bases, 1)
                        else:
                            # No estimate: assume we're halfway through
                            remaining_input = scored_bases

                        if remaining_need <= 0:
                            # Already at target — only accept top reads
                            threshold = score_mean + 2.0 * score_std
                        else:
                            # What fraction of remaining input do we need?
                            need_frac = remaining_need / remaining_input
                            if need_frac >= 1.0:
                                # Need more than what's left — accept everything
                                threshold = 0.0
                                continue
                            # Convert to z-score: need_frac=0.5 -> z=0,
                            # need_frac=0.1 -> z≈1.28, need_frac=0.9 -> z≈-1.28
                            # Using inverse normal approximation
                            # Clamp to avoid extreme thresholds
                            need_frac = max(0.01, min(need_frac, 0.99))
                            # Rational approximation of inverse normal CDF
                            # (Beasley-Springer-Moro, good enough for this)
                            p = 1.0 - need_frac
                            if p <= 0.5:
                                t = math.sqrt(-2.0 * math.log(p))
                                z = t - (2.515517 + 0.802853*t + 0.010328*t*t) / \
                                    (1.0 + 1.432788*t + 0.189269*t*t + 0.001308*t*t*t)
                            else:
                                t = math.sqrt(-2.0 * math.log(1.0 - p))
                                z = -(t - (2.515517 + 0.802853*t + 0.010328*t*t) / \
                                    (1.0 + 1.432788*t + 0.189269*t*t + 0.001308*t*t*t))
                            threshold = score_mean + z * score_std

                # Accept/reject
                if score >= threshold:
                    out.write(hdr_raw)
                    out.write(b'\n')
                    out.write(seq_raw)
                    out.write(b'\n+\n')
                    out.write(qual_raw)
                    out.write(b'\n')
                    accepted_reads += 1
                    accepted_bases += read_len

        finally:
            if fname != '-':
                fh.close()
                if hasattr(fh, '_proc'):
                    fh._proc.wait()

    if args.output != '-':
        out.close()

    # Report
    print(f"[fastq_filter] Total reads:  {total_reads:,}", file=sys.stderr)
    if dup_reads:
        print(f"[fastq_filter] Dupes removed: {dup_reads:,}", file=sys.stderr)
    if short_reads:
        print(f"[fastq_filter] Below min_len: {short_reads:,}", file=sys.stderr)
    print(f"[fastq_filter] Accepted:     {accepted_reads:,} reads, "
          f"{accepted_bases:,} bases ({accepted_bases/1e9:.1f} Gbp)",
          file=sys.stderr)
    if filtering:
        print(f"[fastq_filter] Target was:   {target:,} bases ({target/1e9:.1f} Gbp)",
              file=sys.stderr)
        overshoot = (accepted_bases - target) / target * 100 if target > 0 else 0
        print(f"[fastq_filter] Overshoot:    {overshoot:+.1f}%", file=sys.stderr)
        print(f"[fastq_filter] Final thresh: {threshold:.2f} "
              f"(mean={score_mean:.2f})", file=sys.stderr)


if __name__ == '__main__':
    main()
