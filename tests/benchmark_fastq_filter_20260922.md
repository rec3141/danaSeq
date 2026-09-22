# fastq_filter scaling benchmark, 2026-09-22

This measures the deterministic two-pass implementation committed accidentally
with `01ac977` and documented/tested by `eacc0c4`, against Filtlong v0.3.1.
It is a performance comparison, not a claim that the tools assign identical
scores or retain identical reads.

## Results

All inputs are real nanopore reads from the marine collection. Benchmark-only
unique prefixes were added to read IDs because Filtlong rejects duplicate IDs;
sequences and qualities were unchanged. Input was gzip, output was plain FASTQ,
and the target was 50% of input bases. Each tool ran twice, sequentially on the
same grex Skylake node, with execution order reversed for replicate 2.
`fastq_filter` used `--no_dedupe` for parity with Filtlong.

| input sequence | gzip bytes | fastq_filter mean | Filtlong mean | speedup | peak RSS fastq_filter / Filtlong |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100,003,614 bp | 101,119,977 | 2.70 s | 9.865 s | 3.65x | 6.0 / 243.8 MB |
| 1,000,018,739 bp | 1,009,749,806 | 26.82 s | 98.215 s | 3.66x | 6.4 / 289.4 MB |
| 10,000,001,035 bp | 10,094,857,157 | 270.61 s | 981.025 s | 3.63x | 6.5 / 771.4 MB |

Scaling from 1 to 10 Gbp was 10.09x for `fastq_filter` and 9.99x for
Filtlong. The planned 25 and 50 Gbp measurements were stopped after that
confirmed linear scaling. Dataset staging time is excluded.

Both repetitions of each tool produced identical output SHA256 hashes at every
size. At 10 Gbp, `fastq_filter` retained 5,000,000,478 bases against a target of
5,000,000,517; Filtlong retained 5,000,006,090. Different output read counts and
hashes between tools are expected because their scoring and budget semantics
differ.

## Provenance and limitations

- `fastq_filter.cpp` SHA256:
  `5651ab37198ca1107b298ba384e0eb3aed6841082d01d31ceb207dd00f78b0e7`
- compiled using g++ `-O2 -std=c++11`, matching `nanopore_assembly/install.sh`
- Filtlong v0.3.1, upstream commit `83893671ed26a271880973e87dfd3f3088ea9c09`,
  compiled with its release Makefile (`-O3 -mtune=native`)
- grex Slurm job 7590855; all recorded tool invocations exited 0
- kernel caches were not dropped; reversed execution order limits systematic
  cache bias, and replicate timings were very close
- benchmark harness staging was inefficient and took substantially longer than
  the measured tools; it is outside all reported timings

Raw per-run timing, memory, output count, base count and hashes are in
[`benchmark_fastq_filter_20260922.csv`](benchmark_fastq_filter_20260922.csv).

