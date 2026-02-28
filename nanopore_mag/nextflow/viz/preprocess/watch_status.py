#!/usr/bin/env python3
"""
Pipeline status watcher daemon.

Polls .nextflow.log, trace.txt, and work directories to produce a live
pipeline_status.json file consumed by the Svelte dashboard.

Usage:
  python3 watch_status.py --results /path/to/outdir --output /path/to/viz/data \
      --work-dir /path/to/work --interval 30

  python3 watch_status.py --results /path/to/outdir --output /path/to/viz/data \
      --work-dir /path/to/work --once   # single snapshot, then exit
"""

import argparse
import glob
import json
import os
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# Reuse build_pipeline_status() from preprocess.py
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(__file__))
from preprocess import build_pipeline_status


# ---------------------------------------------------------------------------
# Work directory resolution from .nextflow.log
# ---------------------------------------------------------------------------

def parse_work_dirs(log_path):
    """Map process names to work directory hashes from .nextflow.log.

    Parses lines like:
        [3b/462b91] Submitted process > PROC_NAME (optional_tag)

    Returns: {proc_name: "3b/462b91"}
    """
    work_hashes = {}
    if not os.path.isfile(log_path):
        return work_hashes

    pattern = re.compile(r'\[([0-9a-f]{2}/[0-9a-f]{6})\] Submitted process > (\w+)')
    with open(log_path, 'r', errors='replace') as f:
        for line in f:
            m = pattern.search(line)
            if m:
                hash_prefix, proc_name = m.group(1), m.group(2)
                # Keep the LAST submission (retries overwrite)
                work_hashes[proc_name] = hash_prefix

    return work_hashes


def resolve_work_dir(work_base, hash_prefix):
    """Resolve a hash prefix like '3b/462b91' to a full work directory path."""
    pattern = os.path.join(work_base, hash_prefix + '*')
    matches = glob.glob(pattern)
    if matches:
        return matches[0]
    return None


# ---------------------------------------------------------------------------
# Timestamp parsing from .nextflow.log
# ---------------------------------------------------------------------------

def parse_submit_times(log_path):
    """Extract submission timestamps per process from .nextflow.log.

    Lines look like:
        Feb-27 15:41:58.123 [Task submitter] INFO  ... [3b/462b91] Submitted process > PROC_NAME

    Returns: {proc_name: datetime}
    """
    times = {}
    if not os.path.isfile(log_path):
        return times

    # Nextflow log timestamp format: Mon-DD HH:MM:SS.mmm
    pattern = re.compile(
        r'^([A-Z][a-z]{2}-\d{2} \d{2}:\d{2}:\d{2})\.\d+ .* Submitted process > (\w+)'
    )
    current_year = datetime.now().year

    with open(log_path, 'r', errors='replace') as f:
        for line in f:
            m = pattern.search(line)
            if m:
                ts_str, proc_name = m.group(1), m.group(2)
                try:
                    # Parse with year included to avoid Python 3.15 deprecation
                    dt = datetime.strptime(f"{current_year}-{ts_str}", '%Y-%b-%d %H:%M:%S')
                    times[proc_name] = dt
                except ValueError:
                    pass

    return times


# ---------------------------------------------------------------------------
# Error message extraction
# ---------------------------------------------------------------------------

# Known error signatures: (pattern, human-readable template)
ERROR_SIGNATURES = [
    (
        re.compile(r'torch\.(?:cuda\.)?OutOfMemoryError.*?CUDA out of memory.*?'
                   r'Tried to allocate ([\d.]+ [A-Za-z]+).*?'
                   r'([\d.]+ [A-Za-z]+) free', re.DOTALL),
        lambda m: f"CUDA out of memory (GPU had {m.group(2)} free, needed {m.group(1)})"
    ),
    (
        re.compile(r'CUDA out of memory', re.IGNORECASE),
        lambda m: "CUDA out of memory"
    ),
    (
        re.compile(r'Process exceeded running time limit'),
        lambda m: "Exceeded time limit. Process needs more time."
    ),
    (
        re.compile(r'MemoryError|std::bad_alloc'),
        lambda m: "Out of system memory"
    ),
    (
        re.compile(r'FileNotFoundError:.*?[\'"](.+?)[\'"]'),
        lambda m: f"Missing input file: {m.group(1)}"
    ),
    (
        re.compile(r'\[WARNING\].*produced no bins'),
        lambda m: "Produced no bins (check GPU/memory)"
    ),
]


def extract_error_message(stderr_path):
    """Read .command.err and extract a human-readable error message."""
    if not os.path.isfile(stderr_path):
        return None

    try:
        with open(stderr_path, 'r', errors='replace') as f:
            lines = f.readlines()
    except Exception:
        return None

    # Use last 50 lines for pattern matching
    tail = ''.join(lines[-50:])

    for pattern, formatter in ERROR_SIGNATURES:
        m = pattern.search(tail)
        if m:
            return formatter(m)

    # Generic fallback: last non-empty stderr line (truncated)
    for line in reversed(lines):
        stripped = line.strip()
        if stripped:
            return stripped[:200]

    return None


def check_soft_failure(stderr_path):
    """Check if process had a soft failure (exited 0 but with warnings)."""
    if not os.path.isfile(stderr_path):
        return None

    try:
        with open(stderr_path, 'r', errors='replace') as f:
            content = f.read()
    except Exception:
        return None

    if re.search(r'\[WARNING\].*produced no bins', content):
        return "Produced no bins (check GPU/memory)"
    if re.search(r'\[WARNING\].*exited with code (\d+)', content):
        m = re.search(r'\[WARNING\](.*)', content)
        if m:
            return m.group(1).strip()

    return None


def get_exit_code(work_dir):
    """Read exit code from .exitcode file in work directory."""
    exitcode_path = os.path.join(work_dir, '.exitcode')
    if os.path.isfile(exitcode_path):
        try:
            with open(exitcode_path, 'r') as f:
                return int(f.read().strip())
        except (ValueError, IOError):
            pass
    return None


# ---------------------------------------------------------------------------
# Process-specific progress probes
# ---------------------------------------------------------------------------

def _count_lines(path):
    """Count lines in a file efficiently."""
    if not os.path.isfile(path):
        return 0
    count = 0
    try:
        with open(path, 'rb') as f:
            for _ in f:
                count += 1
    except Exception:
        pass
    return count


def _count_files(pattern):
    """Count files matching a glob pattern."""
    try:
        return len(glob.glob(pattern))
    except Exception:
        return 0


def _format_count(n):
    """Format a number like 1940000 as '1.94M' or 703000 as '703K'."""
    if n >= 1e6:
        return f"{n/1e6:.2f}M"
    if n >= 1e3:
        return f"{n/1e3:.0f}K"
    return str(n)


def probe_calculate_gene_depths(work_dir):
    """Progress probe for CALCULATE_GENE_DEPTHS."""
    done_path = os.path.join(work_dir, 'bedcov_raw.tsv')
    total_path = os.path.join(work_dir, 'genes.bed')

    done = _count_lines(done_path)
    total = _count_lines(total_path)

    if total == 0:
        return None
    progress = min(done / total, 1.0)
    detail = f"{_format_count(done)} / {_format_count(total)} genes"
    return {'progress': round(progress, 3), 'detail': detail}


def probe_minpath(work_dir):
    """Progress probe for MINPATH."""
    done = _count_files(os.path.join(work_dir, 'details', '*.report.txt'))
    total = _count_files(os.path.join(work_dir, 'per_mag', '*.tsv'))

    if total == 0:
        return None
    progress = min(done / total, 1.0)
    detail = f"{done:,} / {total:,} MAGs"
    return {'progress': round(progress, 3), 'detail': detail}


def probe_kofamscan(work_dir):
    """Progress probe for KOFAMSCAN."""
    done = _count_files(os.path.join(work_dir, 'tmp', 'tabular', '*'))
    # KO profile count: check kofam_db profiles directory
    # The total is typically ~27,000 but we count from the work dir if possible
    total_path = os.path.join(work_dir, 'tmp', 'split', '*')
    total = _count_files(total_path)

    if total == 0:
        # Fallback: try to find from profile.desc or just use a reasonable estimate
        return None
    progress = min(done / total, 1.0)
    detail = f"{done:,} / {total:,} KO profiles"
    return {'progress': round(progress, 3), 'detail': detail}


def probe_binning_phase(work_dir, binner_name):
    """Phase-based probe for GPU binners (COMEBin, LorBin, SemiBin2)."""
    phases = []

    if binner_name == 'BIN_COMEBIN':
        phases = [
            ('comebin_out/data_augmentation/', 'Data augmentation'),
            ('comebin_out/train/', 'Training network'),
            ('comebin_out/', 'Binning'),
        ]
    elif binner_name == 'BIN_LORBIN':
        phases = [
            ('lorbin_tmp/', 'Processing contigs'),
            ('lorbin_out/', 'Binning'),
        ]
    elif binner_name == 'BIN_SEMIBIN2':
        phases = [
            ('output_recluster_bins/', 'Reclustering bins'),
            ('output_bins/', 'Generating bins'),
            ('data.csv', 'Preparing features'),
        ]
    elif binner_name == 'BAKTA_EXTRA':
        phases = [
            ('annotation.gff3', 'Writing output'),
            ('annotation.hypotheticals.tsv', 'Annotating hypotheticals'),
            ('annotation.json', 'Running annotation'),
        ]

    # Check phases in reverse order (later phases = more progress)
    for path_suffix, phase_label in phases:
        check_path = os.path.join(work_dir, path_suffix)
        if os.path.exists(check_path):
            return {'progress': None, 'detail': phase_label}

    return None


# Registry: process name -> probe function
PROGRESS_PROBES = {
    'CALCULATE_GENE_DEPTHS': probe_calculate_gene_depths,
    'MINPATH': probe_minpath,
    'KOFAMSCAN': probe_kofamscan,
}

# Phase-based probes (separate because they need the binner name)
PHASE_PROBES = {'BIN_COMEBIN', 'BIN_LORBIN', 'BIN_SEMIBIN2', 'BAKTA_EXTRA'}


def run_probe(proc_name, work_dir):
    """Run the appropriate probe for a process. Returns dict or None."""
    if proc_name in PROGRESS_PROBES:
        return PROGRESS_PROBES[proc_name](work_dir)
    if proc_name in PHASE_PROBES:
        return probe_binning_phase(work_dir, proc_name)
    return None


# ---------------------------------------------------------------------------
# Enhanced status builder
# ---------------------------------------------------------------------------

def build_enhanced_status(results_dir, work_base):
    """Build an enhanced pipeline status dict with timing, progress, and errors.

    Extends build_pipeline_status() with:
      - Per-process timing (submitted, elapsed_min)
      - Progress probes for long-running processes
      - ETA estimates
      - Human-readable error messages for failed processes
      - Work directory paths and key file references
    """
    base_status = build_pipeline_status(results_dir)
    if base_status is None:
        return {
            'processes': {},
            'pipeline_total': 0,
            'pipeline_completed': 0,
            'pipeline_running': 0,
            'pipeline_pending': 0,
            'pipeline_failed': 0,
            'timestamp': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
            'pipeline_active': False,
        }

    log_path = os.path.join(results_dir, 'pipeline_info', 'nextflow.log')

    # Parse work dirs and submit times
    work_hashes = parse_work_dirs(log_path)
    submit_times = parse_submit_times(log_path)

    now = datetime.now()
    enhanced_processes = {}

    for proc_name, status in base_status['processes'].items():
        proc_info = {'status': status}

        # Resolve work directory
        work_dir = None
        if proc_name in work_hashes and work_base:
            work_dir = resolve_work_dir(work_base, work_hashes[proc_name])

        if status == 'running':
            # Add timing info
            if proc_name in submit_times:
                submitted = submit_times[proc_name]
                proc_info['submitted'] = submitted.strftime('%Y-%m-%dT%H:%M:%S')
                elapsed = (now - submitted).total_seconds() / 60
                proc_info['elapsed_min'] = round(elapsed, 1)
            else:
                elapsed = None

            # Run progress probe
            if work_dir:
                proc_info['work_dir'] = work_dir
                proc_info['files'] = {
                    'stderr': '.command.err',
                    'stdout': '.command.log',
                    'script': '.command.sh',
                }
                probe_result = run_probe(proc_name, work_dir)
                if probe_result:
                    proc_info['progress'] = probe_result.get('progress')
                    proc_info['progress_detail'] = probe_result.get('detail')

                    # Estimate ETA from progress + elapsed
                    if (probe_result.get('progress') and elapsed
                            and probe_result['progress'] > 0.01):
                        remaining = elapsed * (1 - probe_result['progress']) / probe_result['progress']
                        proc_info['eta_min'] = round(remaining, 1)
                    else:
                        proc_info['eta_min'] = None
                else:
                    proc_info['progress'] = None
                    proc_info['progress_detail'] = None
                    proc_info['eta_min'] = None
            else:
                proc_info['progress'] = None
                proc_info['progress_detail'] = None
                proc_info['eta_min'] = None

        elif status == 'failed':
            # Extract error message
            if work_dir:
                proc_info['work_dir'] = work_dir
                proc_info['files'] = {
                    'stderr': '.command.err',
                    'stdout': '.command.log',
                    'script': '.command.sh',
                    'trace': '.command.trace',
                }
                stderr_path = os.path.join(work_dir, '.command.err')
                proc_info['error'] = extract_error_message(stderr_path)
                proc_info['exit_code'] = get_exit_code(work_dir)
            else:
                proc_info['error'] = None
                proc_info['exit_code'] = None

        elif status == 'completed':
            # Check for soft failures (warnings in completed processes)
            if work_dir:
                stderr_path = os.path.join(work_dir, '.command.err')
                warning = check_soft_failure(stderr_path)
                if warning:
                    proc_info['status'] = 'warning'
                    proc_info['error'] = warning
                    proc_info['work_dir'] = work_dir
                    proc_info['files'] = {
                        'stderr': '.command.err',
                        'stdout': '.command.log',
                        'script': '.command.sh',
                    }

            # Add duration from trace if available
            trace_path = os.path.join(results_dir, 'pipeline_info', 'trace.txt')
            if os.path.isfile(trace_path):
                try:
                    import pandas as pd
                    trace_df = pd.read_csv(trace_path, sep='\t')
                    proc_rows = trace_df[trace_df['process'] == proc_name]
                    if not proc_rows.empty:
                        duration = proc_rows.iloc[-1].get('duration', '')
                        if isinstance(duration, str):
                            proc_info['duration'] = duration
                except Exception:
                    pass

        enhanced_processes[proc_name] = proc_info

    # Recount with warning status
    status_values = [p['status'] for p in enhanced_processes.values()]
    counts = Counter(status_values)

    pipeline_active = counts.get('running', 0) > 0 or counts.get('pending', 0) > 0

    return {
        'processes': enhanced_processes,
        'pipeline_total': len(enhanced_processes),
        'pipeline_completed': counts.get('completed', 0) + counts.get('warning', 0),
        'pipeline_running': counts.get('running', 0),
        'pipeline_pending': counts.get('pending', 0),
        'pipeline_failed': counts.get('failed', 0),
        'pipeline_warning': counts.get('warning', 0),
        'timestamp': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'pipeline_active': pipeline_active,
    }


# ---------------------------------------------------------------------------
# JSON writer
# ---------------------------------------------------------------------------

def write_json(path, data):
    """Write JSON atomically (write to .tmp then rename)."""
    tmp_path = path + '.tmp'
    with open(tmp_path, 'w') as f:
        json.dump(data, f, separators=(',', ':'))
    os.replace(tmp_path, path)


# ---------------------------------------------------------------------------
# Daemon main loop
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Pipeline status watcher — writes pipeline_status.json for live dashboard')
    parser.add_argument('--results', required=True,
                        help='Pipeline output directory (contains pipeline_info/)')
    parser.add_argument('--output', required=True,
                        help='Viz data directory (where pipeline_status.json is written)')
    parser.add_argument('--work-dir', required=True,
                        help='Nextflow work directory (for progress probes)')
    parser.add_argument('--interval', type=int, default=30,
                        help='Polling interval in seconds (default: 30)')
    parser.add_argument('--once', action='store_true',
                        help='Run once and exit (for testing)')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    output_path = os.path.join(args.output, 'pipeline_status.json')

    print(f"[watch_status] Monitoring: {args.results}", file=sys.stderr)
    print(f"[watch_status] Work dir:   {args.work_dir}", file=sys.stderr)
    print(f"[watch_status] Output:     {output_path}", file=sys.stderr)
    print(f"[watch_status] Interval:   {args.interval}s", file=sys.stderr)

    while True:
        try:
            status = build_enhanced_status(args.results, args.work_dir)
            write_json(output_path, status)

            n_running = status.get('pipeline_running', 0)
            n_completed = status.get('pipeline_completed', 0)
            n_total = status.get('pipeline_total', 0)
            ts = datetime.now().strftime('%H:%M:%S')
            print(f"[watch_status] {ts}  {n_completed}/{n_total} completed, "
                  f"{n_running} running", file=sys.stderr)

            if args.once or not status.get('pipeline_active'):
                if not status.get('pipeline_active'):
                    print("[watch_status] Pipeline no longer active — exiting.",
                          file=sys.stderr)
                break

        except Exception as e:
            print(f"[watch_status] Error: {e}", file=sys.stderr)

        time.sleep(args.interval)


if __name__ == '__main__':
    main()
