# Iterative Metagenomics Pipeline — New Tool Evaluation

**Date:** 2026-02-21
**Purpose:** Document tool selection, conda availability, dependency analysis, and environment isolation strategy for the iterative nanopore metagenomics monitoring pipeline.

## Motivation

The current danaSeq MAG pipeline performs single-pass co-assembly + binning. The iterative pipeline extends this to support recurrent monitoring of aquatic environments:

1. **Growing MAG reference collection** — each run adds new MAGs to a dereplicated catalog
2. **Read recruitment** — incoming reads are recruited against the catalog to track strain dynamics
3. **Strain-level resolution** — new long-read assemblers resolve strains to 99% ANI
4. **Dereplication** — fast ANI-based clustering prevents catalog bloat across time points

## Tool Summary

### Assembly

| Tool | Version | Package | Platform | Description |
|------|---------|---------|----------|-------------|
| metaMDBG | 1.3.1 | `metamdbg` | linux-64, osx-64 | De Bruijn graph assembler for ONT metagenomics. Includes nanoMDBG mode for simplex R10.4+ reads with minimizer-space error correction. |
| myloasm | 0.4.0 | `myloasm` | linux-64 only | High-resolution strain assembler for ONT reads. Resolves strains to 99% ANI in a single command. Rust binary with minimal deps. |

### Dereplication and ANI

| Tool | Version | Package | Platform | Description |
|------|---------|---------|----------|-------------|
| skani | 0.3.1 | `skani` | all platforms | Fast ANI calculation via sparse chaining of k-mer matches. Rust static binary. Required by galah and skDER. |
| galah | 0.4.2 | `galah` | all platforms | Fast MAG dereplication using skani for ANI clustering. Rust binary from Parks lab. |
| sourmash | 4.9.4 | `sourmash-minimal` | noarch | FracMinHash sketching for rapid composition screening, taxonomy (via `gather`), and containment analysis. |
| dRep | 3.6.2 | `drep` | noarch | Established MAG dereplication with quality-aware winner selection. Heavy deps (CheckM v1, mummer4, scikit-learn). |
| skDER | 1.3.4 | `skder` | all platforms | Fast skani-based dereplication + CiDDER protein-coding dereplication. |

### Strain-Level Analysis

| Tool | Version | Package | Platform | Description |
|------|---------|---------|----------|-------------|
| InStrain | 1.10.0 | `instrain` | noarch | Strain-level population genomics. popANI vs conANI for tracking strain dynamics across time points. |
| Floria | 0.0.2 | `floria` | linux-64, osx-64 | Strain-aware phasing from long reads. Produces phased haplotype blocks. Same author as skani/myloasm. |
| Strainy | 1.2 | `strainy` | noarch | Assembly-graph-based strain phasing. Modifies GFA graphs to separate strains. |

## Dependency Conflict Analysis

Several Python-heavy tools have conflicting version pins that make co-installation impossible:

### Critical Conflicts

| Tool | Conflict | Severity |
|------|----------|----------|
| **InStrain** | `biopython <=1.74` (current is 1.84+) | **CRITICAL** — incompatible with any tool needing modern biopython |
| **skDER** | `python >=3.10,<3.11.0a0` | **HIGH** — locks to Python 3.10 exactly |
| **Strainy** | `numpy <1.27`, `scipy <1.13`, `networkx <3.4`, `pysam <0.23` | **HIGH** — tight upper bounds conflict with newer tools |
| **dRep** | `checkm-genome` (pulls pplacer, prodigal, HMMER, mummer4) | **MEDIUM** — heavy transitive deps, solver conflicts |

### Safe Combinations (Rust binaries, minimal deps)

The following tools are compiled Rust/C++ binaries with no Python dependency conflicts:
- metaMDBG (C++ with minimap2)
- myloasm (Rust)
- skani (Rust)
- galah (Rust + skani)
- Floria (Rust)

These could theoretically share an environment, but we follow the pipeline convention of one env per logical function.

## Environment Isolation Strategy

Eight new conda environments, following the `dana-mag-*` naming convention:

| Environment | YAML | Tools | Rationale |
|-------------|------|-------|-----------|
| `dana-mag-metamdbg` | `metamdbg.yml` | metaMDBG 1.3.1 | ONT metagenome assembly; C++ binary, includes minimap2 |
| `dana-mag-myloasm` | `myloasm.yml` | myloasm 0.4.0 | Strain-resolution assembly; Rust binary, Linux-only |
| `dana-mag-derep` | `derep.yml` | galah 0.4.2, skani 0.3.1, sourmash-minimal | Fast dereplication + ANI + sketching; all Rust-based, safe to combine |
| `dana-mag-drep` | `drep.yml` | dRep 3.6.2 | Established dereplication; isolated due to heavy deps (CheckM v1, mummer4) |
| `dana-mag-instrain` | `instrain.yml` | InStrain 1.10.0 | Strain tracking; **must be isolated** (biopython <=1.74 pin) |
| `dana-mag-strainy` | `strainy.yml` | Strainy 1.2 | Strain phasing; isolated due to tight numpy/scipy/pysam pins |
| `dana-mag-floria` | `floria.yml` | Floria 0.0.2 | Strain phasing; Rust binary, pre-release versioning |
| `dana-mag-skder` | `skder.yml` | skDER 1.3.4 | Fast dereplication; isolated due to Python 3.10 pin |

## Installation

All environments are built by the existing `install.sh`:

```bash
cd /data/danav2/nanopore_mag/nextflow
./install.sh          # builds all environments (existing + new)
./install.sh --check  # verify installations
```

Individual environments can be built manually:

```bash
mamba env create -p conda-envs/dana-mag-metamdbg -f envs/metamdbg.yml
mamba env create -p conda-envs/dana-mag-derep -f envs/derep.yml
# etc.
```

## Database Requirements

| Tool | Database | Size | Download |
|------|----------|------|----------|
| sourmash | GTDB representative genomes | ~2 GB | `sourmash databases --list` or download from [sourmash.bio](https://sourmash.bio/databases/) |
| dRep | CheckM v1 data | ~1.4 GB | `checkm data setRoot /path/to/checkm_data` |
| InStrain | None (uses BAM + FASTA) | — | — |
| skani | None (self-contained) | — | — |
| galah | None (uses skani internally) | — | — |

## Intended Pipeline Integration

### Phase 1: Assembly + Dereplication (near-term)

```
Incoming reads → metaMDBG assembly → galah dereplication against catalog
                                   → sourmash gather for composition screening
                                   → skani dist for pairwise ANI
```

### Phase 2: Strain Resolution (medium-term)

```
Per-MAG read recruitment → myloasm strain assembly (99% ANI resolution)
                         → Floria haplotype phasing
                         → InStrain profile for population-level metrics
```

### Phase 3: Longitudinal Tracking (long-term)

```
Time series → InStrain compare across time points
            → popANI / conANI strain identity tracking
            → Strainy graph-based strain separation
            → dRep quality-aware catalog updates
```

## Tool Provenance

| Tool | Author(s) | Publication | License |
|------|-----------|-------------|---------|
| metaMDBG | Benoit et al. | Nature Biotechnology 2024 | MIT |
| myloasm | Shaw, Marin, Li | bioRxiv 2025 | MIT |
| galah | Woodcroft (Parks lab) | — | GPL-3.0 |
| skani | Shaw & Yu | Nature Methods 2023 | MIT |
| sourmash | Brown & Irber | JOSS 2016 | BSD-3 |
| dRep | Olm et al. | ISME Journal 2017 | MIT |
| InStrain | Olm et al. | Nature Biotechnology 2021 | MIT |
| Floria | Shaw et al. | Bioinformatics 2024 | MIT |
| Strainy | Kazantseva et al. | bioRxiv 2023 | GPL-3.0 |
| skDER | Rauf et al. | Microbial Genomics 2025 | BSD-3 |

## Notes

- **myloasm** is Linux-only (no macOS builds). This matches our shipboard deployment target.
- **Floria** v0.0.2 is early versioning despite being published. The tool is functional.
- **InStrain**'s `biopython <=1.74` pin is a known limitation. The authors have not updated it as of v1.10.0.
- **galah + skani + sourmash** are bundled in one env (`derep`) because they are all Rust-based with no Python conflicts and serve the same logical function (fast ANI/dereplication/sketching).
- **dRep** is kept separate from the galah/skani env despite overlapping function because dRep's CheckM v1 dependency is extremely heavy and conflict-prone. For routine dereplication, prefer galah; use dRep when quality-aware winner selection is needed.
