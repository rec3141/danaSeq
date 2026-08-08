# Parameters

All parameters can be set via command line (`--param value`) or in a params file.

```bash
nextflow run rec3141/microscape-nf --help
```

---

## Input / Output

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | Directory of paired-end `*.fastq.gz` files |
| `--outdir` | `results` | Output directory |
| `--store_dir` | off | Persistent cache directory (skip completed steps across runs) |

---

## Language and Engine

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--lang` | `python` | Script language: `python` (papa2 + microscape) or `R` (dada2 + microscapeR) |
| `--dada_engine` | follows `--lang` | Override DADA2 engine independently: `python` or `R` |

---

## Primer Removal

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_primer_removal` | `false` | Skip if input is already trimmed |
| `--primer_auto` | `true` | Auto-select primer file by filename prefix (16S/18S/ITS) |
| `--primers` | auto | Override with a specific primer FASTA |
| `--primer_error_rate` | `0.12` | Cutadapt max error rate |

---

## DADA2 Quality Filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--maxEE` | `2` | Max expected errors per read |
| `--truncQ` | `11` | Truncate at first base with quality <= Q |
| `--maxN` | `0` | Max Ns allowed |
| `--truncLen_fwd` | `0` | Truncate forward reads at position N (0 = off) |
| `--truncLen_rev` | `0` | Truncate reverse reads at position N (0 = off) |
| `--min_overlap` | `10` | Minimum overlap for paired-read merging |

---

## QC Filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_seq_length` | `50` | Remove ASVs shorter than N bp |
| `--min_samples` | `2` | Remove ASVs present in fewer than N samples |
| `--min_seqs` | `3` | Remove ASVs with fewer than N total reads |
| `--min_reads` | `100` | Remove samples with fewer than N reads |

---

## Taxonomy

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--ref_databases` | required | Reference databases, semicolon-separated: `"name:path:Level1,Level2,..."` |
| `--run_phylogeny` | `false` | Build MSA + NJ phylogenetic tree |

**Example:**
```bash
--ref_databases "silva:/db/silva_nr99_v138.1_train_set.fa.gz:Domain,Phylum,Class,Order,Family,Genus;pr2:/db/pr2_train.fasta:Domain,Supergroup,Division,Class,Order,Family,Genus,Species"
```

---

## Metadata

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--metadata` | off | Path to MIMARKS-compliant TSV/CSV metadata file |
| `--sample_id_column` | `sample_name` | Column name in metadata matching sample IDs |

---

## Network

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_prevalence` | `5` | Min samples an ASV must appear in for network analysis |

---

## Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--threads` | `8` | General thread count |
| `--dada_cpus` | `8` | CPUs for DADA2 processes |
| `--dada_memory` | `16 GB` | Memory for DADA2 processes |
