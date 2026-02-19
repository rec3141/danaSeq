# Nextflow Caching and Resume: How It Works

This document explains how Nextflow caching works in the Dana MAG pipeline,
why concurrent runs sometimes fail to share cached data, and how to work around it.

## The Core Concept

Every Nextflow task gets a **128-bit hash** computed from its inputs. When you
run with `-resume`, Nextflow checks if a task with that hash already completed
successfully. If so, the task is skipped ("cached"). If not, it re-executes.

The hash is computed from **all** of the following:

| Hash input           | Example                                    |
|----------------------|--------------------------------------------|
| Session ID (UUID)    | `701ae656-bd37-4b35-9418-8a32917137cd`     |
| Process name         | `ASSEMBLY_FLYE`                            |
| Script block content | The literal bash script after interpolation |
| Input files          | Path + size + mtime (default mode)         |
| Input values         | Channel values (strings, numbers, booleans)|
| Conda environment    | `envs/flye.yml` contents                   |
| Container image      | Docker/Singularity tag                     |
| `bin/` scripts       | Any bundled scripts used by the task       |
| Global variables     | Workflow-scope variables used in script    |

**The session ID is part of the hash.** This is the single most important thing
to understand. It means:

- Two runs with **different** session IDs will **never** share cached tasks
- To reuse cached results, you must `-resume` the **same** session

## Where Cache Data Lives

Nextflow stores cache data in two places:

```
nextflow/
├── .nextflow/
│   ├── history              # Log of all runs (session IDs, run names, timestamps)
│   └── cache/
│       └── <session-id>/
│           └── db/          # LevelDB: hash → task metadata
│               ├── LOCK     # Exclusive lock file (only one reader/writer)
│               ├── *.ldb    # Data files
│               └── ...
└── work/
    └── <hash-prefix>/<hash>/    # Actual task files
        ├── .command.sh          # The script that ran
        ├── .command.log         # stdout
        ├── .command.err         # stderr
        ├── .exitcode            # Exit status
        └── <outputs>            # Task output files
```

**Both** must be intact for resume to work:
- `.nextflow/cache/` tells Nextflow "this hash was previously executed"
- `work/` contains the actual output files to verify and reuse

## How `-resume` Behaves

### `nextflow run main.nf -resume` (no session specified)

Resumes the **most recent** run from `.nextflow/history`. Uses that run's session
ID, so all task hashes match the previous execution.

### `nextflow run main.nf -resume <session-id>`

Resumes a **specific** prior run. Useful when you want to resume an older run,
not the most recent one.

### `nextflow run main.nf` (no `-resume`)

Starts fresh with a **new** session ID. All tasks re-execute, even if identical
work directories exist from prior runs. The old cached data is not consulted.

### `run-mag.sh` behavior

`run-mag.sh` auto-detects the session ID from
`<outdir>/pipeline_info/run_command.sh` if it exists. This means re-running
with the same `--outdir` will automatically resume the correct session. Using a
**new** `--outdir` starts a new session (no cache reuse). You can override with
`--session <uuid>`.

## Why Concurrent Runs Fail

### Problem: LevelDB Lock

The `.nextflow/cache/<session>/db/LOCK` file is an **exclusive lock**. Only one
Nextflow process can open a given session's cache at a time. This means:

- You **cannot** run two pipelines from the same launch directory simultaneously
- You **cannot** resume the same session from two places at once
- Even `nextflow log` is blocked while a pipeline is running

### Problem: Different Sessions Don't Share Cache

If you launch a second run with a different session ID (including a brand-new
`--session $(uuidgen)`), it creates fresh hashes that don't match any existing
work directories. Everything re-executes from scratch.

### Problem: Resuming a Session That Never Ran Certain Tasks

If you resume session A but want to run tasks that session A never executed,
those tasks will run fresh. Only tasks that session A **successfully completed**
are cached.

## Practical Scenarios

### Scenario 1: Re-run with an extra module (the common case)

You ran the pipeline once and now want to add `--run_marferret true`.

**What works:**
```bash
# Resume from the same outdir (run-mag.sh auto-detects session)
./run-mag.sh --input /data/reads --outdir results_full --run_marferret true --marferret_db /path/to/db
```
All previously completed tasks are cached. Only MarFERReT (and any new
dependencies) actually executes.

**What fails:**
```bash
# New outdir = new session = no cache
./run-mag.sh --input /data/reads --outdir results_marferret_test --run_marferret true --marferret_db /path/to/db
```
Everything re-executes because the new session has no history.

### Scenario 2: Concurrent runs

You want to test a new module while a long-running pipeline is still going.

**What fails:**
- Running from the same directory → LevelDB lock conflict
- Resuming the running session → lock conflict
- New session → no cached data, full re-execution

**What works (but requires setup):**
1. Wait for the first run to complete
2. Resume that session with the extra flags added

**Alternative for concurrent testing:**
Use a separate launch directory with a shared work directory:
```bash
mkdir /tmp/marferret-test && cd /tmp/marferret-test
nextflow run /data/danav2/nanopore_mag/nextflow/main.nf \
    -w /data/danav2/nanopore_mag/nextflow/work \
    --input /data/reads --outdir results_test \
    --run_marferret true --marferret_db /path/to/db
```
This avoids the lock conflict but starts a **new session**, so no cache reuse
occurs. The shared `-w` means work directories aren't duplicated on disk, but
all tasks still re-execute.

**True concurrent caching** requires the cloud cache store
(`NXF_CLOUDCACHE_PATH`), which uses object storage instead of LevelDB and
supports multiple readers/writers. This is not practical for local deployments.

### Scenario 3: A process fails and kills unrelated tasks

If process X fails and Nextflow terminates the pipeline, any running process Y
that gets killed (exit 143 = SIGTERM) will **not** be cached. On resume, Y must
re-execute from scratch even though X was the root cause.

**Mitigation** (applied in this pipeline):
```groovy
errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore' }
```
This ensures failed processes are **ignored** rather than terminating the
pipeline, so unrelated processes continue to completion.

## What Invalidates the Cache

Any change to a hash input forces re-execution of the affected task **and all
downstream tasks**. Common causes:

| Change                              | Effect                              |
|-------------------------------------|-------------------------------------|
| Edit the process script block       | That task + all downstream re-run   |
| Edit a `bin/` script used by task   | That task + all downstream re-run   |
| Change `--assembly_cpus` or memory  | Processes using those params re-run |
| Touch/modify an input file          | That task + all downstream re-run   |
| Change conda YAML                   | That task + all downstream re-run   |
| Add a new `--run_X true` flag       | Only new tasks execute (if wired correctly) |
| Change `--outdir`                   | **No effect on caching** (publishDir only) |

**Debugging cache misses:** Run with `-dump-hashes` and compare logs:
```bash
nextflow -log run1.log run main.nf -dump-hashes -resume
# diff run1.log run2.log to see which hash inputs changed
```

## Persistent Caching with storeDir

### The Problem with `-resume`

`-resume` relies on the work directory and `.nextflow/cache/` being intact. If either
is deleted (e.g., disk cleanup, accidental `rm -rf work/`), all cached results are
lost and expensive processes like assembly (~2h), binning (~1h), and annotation (~1h)
must re-run.

### How storeDir Works

`storeDir` provides a separate, persistent cache layer. When enabled via
`--store_dir /path/to/store`:

1. Each process stores outputs directly in `store_dir/<category>/<subcategory>/`
2. On subsequent runs, Nextflow checks if all declared output files exist in storeDir
3. If all outputs are present, the process is **skipped entirely** — no work directory
   task is created, no resources are consumed
4. Nextflow shows "Stored" status for skipped processes in the log

### Key Differences from `-resume`

| Feature                | `-resume`                    | `storeDir`                    |
|------------------------|------------------------------|-------------------------------|
| Session-dependent?     | Yes (hash includes session)  | No (just checks file existence) |
| Survives `rm -rf work/`? | No                        | Yes                           |
| Survives new session?  | No                           | Yes                           |
| Output location        | work/ + publishDir copy      | storeDir (publishDir ignored) |
| Cache granularity      | Per-task hash                | Per-process output files      |
| Concurrent run support | No (LevelDB lock)            | Yes (file-based)              |

### publishDir Interaction

**When storeDir is active, publishDir is silently ignored.** Outputs go directly to the
storeDir path instead of being copied to the outdir. The storeDir structure mirrors the
publishDir layout so paths are consistent.

### Seeding storeDir from Existing Results

Use `seed-store-dir.sh` to populate a storeDir from a previous pipeline run without
re-executing:

```bash
# Auto-detect: hardlink (same filesystem) or symlink (cross-filesystem)
./seed-store-dir.sh <results_dir> <store_dir>

# Force a specific mode
./seed-store-dir.sh --mode hardlink <results_dir> <store_dir>
./seed-store-dir.sh --mode symlink  <results_dir> <store_dir>
./seed-store-dir.sh --mode copy     <results_dir> <store_dir>
```

The script handles legacy naming conventions:
- Binner TSVs: `contig_bins.tsv` → `{name}_bins.tsv` (required for storeDir)
- Prokka outputs: `PROKKA_<timestamp>.*` → `annotation.*`
- MAP_TO_BINS: Combines `per_mag/` + `community_annotations.tsv` into single storeDir

### When to Use storeDir vs. `-resume`

- **`-resume` only**: Sufficient for quick one-off runs where you don't need durability
- **`--store_dir`**: Recommended for both testing and production — survives work directory
  cleanup, shareable across sessions, and protects against accidental cache loss
- **Both together**: They complement each other. `-resume` handles the fast path
  (hash match → skip), storeDir handles the durable path (files exist → skip).

## Best Practices for This Pipeline

1. **Always use `run-mag.sh`** — it records the session ID for reliable resume
2. **Add flags, don't change them** — appending `--run_marferret true` to an
   existing run command preserves all cached tasks
3. **Never change `--input` or assembly params** on a resume — it invalidates
   the assembly and forces a full re-run (hours on real data)
4. **One run at a time** — Nextflow's local cache doesn't support concurrency;
   wait for the current run to finish before starting another
5. **Use `--outdir` consistently** — `run-mag.sh` binds session IDs to output
   directories; switching `--outdir` means losing cache
6. **Use `--store_dir`** — persistent caching that survives work directory cleanup;
   seed it from existing results with `seed-store-dir.sh`

## References

- [Nextflow: Caching and resuming](https://www.nextflow.io/docs/latest/cache-and-resume.html)
- [Seqera: Demystifying Nextflow resume](https://seqera.io/blog/demystifying-nextflow-resume)
- [Seqera: Troubleshooting resume](https://seqera.io/blog/troubleshooting-nextflow-resume)
- [Nextflow training: Cache and resume](https://training.nextflow.io/2.1/basic_training/cache_and_resume/)
