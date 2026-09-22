# fastq_filter regressions

Run `python3 tests/test_fastq_filter.py -v` from the repository root. This
compiles the current C++ source into a temporary directory and runs small
fixtures locally; it does not replace the installed executable. Requires
Python 3, g++, zlib development headers, and Linux `/dev/full`/resource limits.
Set `FASTQ_FILTER_TEST_BINARY` to test a separately built binary (for example,
a sanitizer build).

## Two-pass contract

Eligible reads are visited in descending score order, then ascending read-ID
hash, then lexicographic full-record order. Identical records are
interchangeable. Every bucket is sorted before emission, not just the bucket
containing the target cutoff. Output records and order are independent of
input permutation and bucket count, for identical input records and options.

Selection is greedy: skip reads that exceed the remaining base budget and
continue into lower score buckets. It is not a knapsack optimizer, so the
target need not be reached exactly. With deduplication enabled, the first
record in this canonical order represents each read ID, even if that record
does not fit the remaining budget. This chooses the highest-scoring eligible
representative instead of the first record encountered in the input stream.
Duplicate representatives are resolved by full record content, not arrival
order. `--onepass` and percentage mode retain their streaming semantics.

Scratch files live in a private directory under `--spill_dir`; at most 512
spill writers are open concurrently, with a lower cap derived from the process
descriptor limit (16 writers with a limit of 32). Ordinary exceptions, input errors,
spill/output errors and broken pipes return nonzero and remove scratch files.
SIGKILL or machine failure can still leave the private directory behind.
On failure, output may be partial and must not be used.

Sorting every emitted bucket adds scoring/sorting and seek work compared
with copying whole buckets verbatim. Candidate indices use roughly 32 bytes
per read on 64-bit systems, one bucket at a time, plus the optional ID dedup
set. A concentrated score distribution can put most reads in one bucket;
bucket count is not a hard memory bound. Large-workload performance has not
been benchmarked for this correction.

Coverage includes complete-bucket permutations, lower-bucket budget filling,
64/512/4096-bucket equivalence, duplicate IDs with/without deduplication,
descriptor-limited LRU eviction/reopening, buffered and direct spill failures,
plain/gzip/stdout output failures, malformed/truncated input, and streaming
compatibility.
