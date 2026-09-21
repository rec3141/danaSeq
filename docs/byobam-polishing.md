# BYOBAM: supply the polishing alignment instead of letting Flye build it

Status: proposal. The Flye side already exists upstream; the pipeline side does not.

## Why

Timings from the marine co-assembly (8.9 Gbp, 1,176,983 contigs, 250 Gbp reads,
192-core node, work dir on Lustre):

| phase | wall-clock |
|---|---|
| minimap2 alignment + sort | 52 min |
| `samtools index -c -@ 4` | **5h53m** |
| separating alignment into bubbles | 2h04m |
| correcting bubbles | 2h53m |
| **polish total** | **11h42m** |

The alignment is not the problem. Indexing is half the wall-clock, at roughly
4 MB/s over an 85 GB BAM, on four threads hardcoded in
`flye/polishing/alignment.py`. Lustre contention makes it worse -- the
freshwater assembly was writing its own 96 GB BAM over the same window.

Separately, `MAP_READS` then maps every sample against the polished assembly:
276 tasks, each rebuilding a split minimap2 index over the same 8.9 Gbp
reference, ~1.5 h and 48 GB peak per task. That is a second full mapping pass
over reads Flye has already aligned.

## What Flye already supports

`flye --polish-target` accepts a BAM in place of reads:

- `flye/main.py` sets `bam_input` when a reads argument ends in `.bam`
- `flye/polishing/polish.py` then skips `make_alignment` entirely
- and, critically, does **not** delete the file:
  `if not bam_input: os.remove(alignment_file)`

Constraints, both enforced in `main.py` and both already satisfied by the
pipeline: a single BAM only, and `--iterations 1`.

So Flye does not need patching for this. The alignment simply becomes ours.

## Proposed shape

```
PREPARE_READS ──► all_reads.fastq + read_map.tsv.gz
                        │
                  MAP_TO_DRAFT            one minimap2 over the draft assembly,
                        │                 our threads, sort --write-index,
                        │                 work dir on node-local scratch
                        ├──► FLYE_POLISH --polish-target draft.fasta aln.bam
                        │                 (skips its own mapping and indexing)
                        └──► DEMUX_BAM    split by read_map into per-sample BAMs
                                          published to mapping/ as today
```

`MAP_READS` stays for metaMDBG, myloasm and `--polish false`, where no such
alignment exists.

## What it buys

- removes the 5h53m index pass: `samtools sort --write-index` writes the CSI
  during sorting (samtools >= 1.10; danaSeq pins >= 1.17 in the flye env, and
  Flye now prefers a PATH samtools over its vendored 1.9)
- removes the second mapping pass entirely (276 tasks, ~1.5 h, 48 GB each)
- removes the failure mode those tasks caused: on the marine run 116 task
  failures across 39 of 276 samples, each retried once then silently ignored,
  leaving 14% of samples absent from the depth matrix
- one alignment, one set of flags, used for both polishing and coverage

## What it costs / open questions

- **Coordinates are the draft, not the polished assembly.** Contig names are
  identical (`normalize_assembly.py` renumbers before polishing and Flye
  preserves names), so all binners are unaffected -- they work per contig and
  match by name. Only `CALCULATE_GENE_DEPTHS` is coordinate-sensitive, and its
  sole consumer is viz. Accepted 2026-09-20.
- **Read set differs from today's `MAP_READS`.** The polish alignment is
  post-dedupe, post-filtlong, post-REMOVE_HUMAN; `MAP_READS` maps the raw
  per-barcode reads. Arguably an improvement, but it is a change.
- **Flye's own flags differ from `MAP_READS`.** Flye emits `-N 10 -p 0.5
  --secondary-seq`; `MAP_READS` uses `--secondary=no`. If we own the alignment
  we choose, but the polisher wants secondaries and coverm does not -- so the
  demux must filter `-F 0x100`, as the prototype already does.
- **Disk.** The alignment must survive from mapping through polishing and
  demultiplexing; 85 GB for marine, 96 GB for freshwater. Node-local scratch is
  the right home (and would likely have avoided the 4 MB/s index), but it has to
  outlive a single process, not a single task.
- **Validation still outstanding.** Demuxed depths have not yet been compared
  against a full `MAP_READS` run. The marine run that would have provided the
  baseline was cancelled because 39 samples were missing; a re-run with the
  corrected memory label (danaSeq 2ad6fb7) would provide it.

## The other samtools cost: per-region spawns

BYOBAM removes the indexing pass, but not this one, because it lives inside
Flye's polishing loop regardless of who built the alignment.

Flye calls samtools in six places:

| site | call | frequency |
|---|---|---|
| `alignment.py:265` | `view -T ref -u -` | once (pipe stage) |
| `alignment.py:266` | `sort -@ 4 -l 1 -m MEM` | once |
| `alignment.py:277` | `index -c -@ 4` | once -- 5h53m on marine |
| `sam_parser.py:358` | `depth -r ctg:s-e -a -m 0 -Q 10 -l 100` | **per region** |
| `sam_parser.py:378` | `view bam ctg:s-e` | **per region** |
| `sam_parser.py:443` | `view bam` (whole file) | once |

`SynchonizedChunkManager` emits one region per contig, and both `make_bubbles`
and `consensus` drive it. On the marine assembly that is 1,176,983 subprocess
spawns; freshwater's BAM has 1,325,561 references.

Measured on fir, samtools 1.22, page-cache warm:

```
samtools view -H <96.7 GB bam>   ->  0.59 s, 364 MB RSS
header size                      ->  40,032,496 bytes / 1,325,561 @SQ records
```

Every region query pays that before reading a single alignment:

```
1,325,561 regions x 0.59 s ~= 217 core-hours of pure header parsing
```

which across 128 workers is ~1.7 h of wall-clock -- close to the 2h04m observed
for bubbles, and again for consensus (2h53m).

Newer samtools helps, but does not fix this:
- `samtools depth` was rewritten in 1.13 and is materially faster
- htslib >= 1.10 hashes reference names; older paths did linear scans, which is
  pathological at 1.3M refs
- BGZF threading improved, which matters for sort/index given more than 4 threads
- but the 0.59 s above IS 1.22. A newer binary cannot fix re-parsing a 40 MB
  header 1.3M times.

The fix is architectural: batch many regions per invocation (`samtools view`
accepts multiple regions, or `-M -L regions.bed`), or hold one reader open per
worker instead of forking per contig. Batching 1000 contigs per spawn would cut
header parses ~1000x. This only bites assemblies with very many contigs -- i.e.
exactly the metagenome case, and not the isolate case Flye was tuned for.

## Prototype status

`demux_bam_by_sample.py` splits a pooled BAM by `read_map.tsv.gz`. Verified
120/120 on synthetic data with known ground truth, and 108/108 correct with zero
misassignments on the real marine polishing BAM. Read maps for both habitats are
built and validated (marine 276 samples / 97,363,917 reads; freshwater 232 /
87,815,882) -- every read in every staged file is accounted for, because the
sample comes from the filename rather than the header.
