# BBSketch GTDB TaxServer — Setup & Usage Guide

A local sketch-based taxonomy server using BBTools and GTDB, providing millisecond-speed taxonomic classification of contigs, MAGs, or assemblies against the full GTDB database.

## Overview

The BBTools `TaxServer` loads a pre-built sketch database into memory and serves queries over HTTP. Combined with GTDB (Genome Taxonomy Database), this gives you fast, indexed taxonomy lookups against ~144K reference genomes — no internet required.

**Performance example:** 759 contigs classified in **5.5 seconds** against 143,614 GTDB genomes.

---

## Prerequisites

- **BBTools** (tested with v39.52): `conda install -c bioconda bbmap`
- **Java** (JDK 11+): comes with BBTools via conda
- **RAM**: The full GTDB R226 database requires ~40 GB RAM when indexed. Plan for at least 64 GB; 128 GB gives comfortable headroom.
- **Disk**: ~15 GB for the sketch file + taxonomy tree

---

## Part 1: Building a GTDB Sketch Database

### 1.1 Download GTDB Metadata and Taxonomy

```bash
# Get GTDB taxonomy and metadata
# Visit https://gtdb.ecogenomic.org/downloads for latest release
wget https://data.gtdb.ecogenomic.org/releases/latest/ar53_taxonomy.tsv.gz
wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz

# Combine into one taxonomy file
zcat ar53_taxonomy.tsv.gz bac120_taxonomy.tsv.gz > gtdb_taxonomy.tsv
```

### 1.2 Download GTDB Representative Genomes

```bash
# Get genome paths from GTDB metadata
# You'll need the genomic_files_all directory or individual genome FTPs
# For representative genomes (~144K), download from NCBI using accessions

# Example: download genomes listed in GTDB metadata
# (This can take a while — ~144K genomes)
mkdir -p genomes/
# Use your preferred method: datasets CLI, wget from NCBI FTP, etc.
```

### 1.3 Build the NCBI-style Taxonomy Tree

BBTools TaxServer requires an NCBI-format taxonomy tree. If your genomes use NCBI taxIDs, you need:

```bash
# Download NCBI taxonomy dump
mkdir -p taxonomy && cd taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xzf taxdump.tar.gz

# Build the BBTools taxtree
taxtree.sh names.txt nodes.txt merged.txt tree.taxtree.gz
cd ..
```

### 1.4 Sketch All Genomes

```bash
# Sketch all genomes into a single sketch file
# -Xmx sets Java heap; adjust based on available RAM
# k=31 is standard for nucleotide sketches
# sizemult=1 uses default autosize (~10K kmers per bacterial genome)

bbsketch.sh \
    in=genomes/*.fna.gz \
    out=gtdb_r226_k31.sketch \
    k=31 \
    autosize=t \
    processaliases=t \
    -Xmx16g

# Check the result
ls -lh gtdb_r226_k31.sketch
grep -c "^#SZ:" gtdb_r226_k31.sketch  # Count sketches
```

**Notes on sketch building:**
- `processaliases=t` — resolves NCBI accession→taxID mappings
- For large genome sets, you may need to sketch in batches and merge with `mergesketch.sh`
- The sketch file is plain text (one header + data block per genome)
- Our full GTDB R226 sketch: **11 GB, 143,614 sketches**

### 1.5 Building a Custom/Subset Database

To build a database from your own genomes or a subset:

```bash
# Sketch individual genomes
for genome in genomes/*.fasta; do
    name=$(basename $genome .fasta)
    sketch.sh in=$genome out=sketches/${name}.sketch k=31
done

# Merge into one file
cat sketches/*.sketch > my_custom.sketch

# Or use bbsketch.sh directly on multiple files
bbsketch.sh in=genomes/*.fasta out=my_custom.sketch k=31
```

---

## Part 2: Starting the TaxServer

### 2.1 Basic Server Launch

```bash
# Start the server with indexing enabled
java -da -Xmx110g \
    -cp /path/to/bbmap/current/ \
    tax.TaxServer \
    tree=taxonomy/tree.taxtree.gz \
    port=3068 \
    k=31 \
    index \
    /path/to/gtdb_r226_k31.sketch
```

**Parameter guide:**

| Parameter | Description |
|-----------|-------------|
| `-Xmx110g` | Java max heap — must be large enough for sketches + index |
| `-da` | Disable assertions (faster) |
| `tree=` | Path to taxonomy tree file |
| `port=` | HTTP port to serve on (default: 3068) |
| `k=31` | K-mer size (must match sketch) |
| `index` | Build a k-mer index for fast lookups (highly recommended) |

### 2.2 Memory Requirements

The index uses hash tables that double in size at load thresholds. Empirical measurements for GTDB R226:

| Sketches | Sketch File | RAM Used | Notes |
|----------|-------------|----------|-------|
| 100 | 7.8 MB | 2.3 GB | Fixed overhead dominates |
| 1,000 | 74 MB | 2.4 GB | Still mostly overhead |
| 5,000 | 360 MB | 4.1 GB | Index starting to grow |
| 10,000 | 723 MB | 4.5 GB | |
| 20,000 | 1.5 GB | 9.0 GB | |
| 30,000 | 2.2 GB | 9.9 GB | |
| 40,000 | 2.9 GB | 26.4 GB | Hash table resize here |
| 50,000 | 3.7 GB | 27.3 GB | |
| 75,000 | 5.5 GB | 29.6 GB | |
| 100,000 | 7.3 GB | 31.9 GB | |
| **143,614** | **11 GB** | **39.7 GB** | **Full GTDB R226** |

**Key insight:** Memory scaling is NOT linear — it follows a staircase pattern due to hash table doubling. Set `-Xmx` to at least 2.5× the "RAM Used" value to allow headroom during index construction.

### 2.3 Running as a Background Service

```bash
# Run detached with nohup
nohup java -da -Xmx110g \
    -cp /path/to/bbmap/current/ \
    tax.TaxServer \
    tree=taxonomy/tree.taxtree.gz \
    port=3068 k=31 index \
    /path/to/gtdb_r226_k31.sketch \
    > taxserver.log 2>&1 &

echo "TaxServer PID: $!"

# Monitor startup (takes ~5 minutes for full GTDB)
tail -f taxserver.log
# Wait for "Ready!" message
```

### 2.4 Running Without Index

If RAM is limited, you can skip the index:

```bash
# Omit the 'index' flag — queries will be slower (brute-force scan)
java -da -Xmx64g -cp /path/to/bbmap/current/ \
    tax.TaxServer tree=taxonomy/tree.taxtree.gz \
    port=3068 k=31 \
    /path/to/gtdb_r226_k31.sketch
```

**Trade-off:** Without index, `unique2`, `unique3`, `contam2`, and `refhits` output columns are unavailable, and query speed drops significantly.

---

## Part 3: Querying the Server

### 3.1 Basic Query with sendsketch.sh

```bash
# Classify a single assembly (all contigs merged into one sketch)
sendsketch.sh \
    in=my_assembly.fasta \
    address=http://HOSTNAME:3068/sketch \
    k=31 \
    records=5 \
    format=3

# Classify each contig separately
sendsketch.sh \
    in=my_assembly.fasta \
    address=http://HOSTNAME:3068/sketch \
    k=31 \
    records=1 \
    persequence \
    format=3 \
    out=results.txt
```

**Important:** The address must include `http://` and the `/sketch` path.

### 3.2 Network Access

The TaxServer binds to all interfaces (`0.0.0.0`) by default. Any machine on the network can query it:

```bash
# From another machine on the network
sendsketch.sh \
    in=sample.fasta \
    address=http://10.151.50.41:3068/sketch \
    k=31 records=3 format=3
```

For access from outside the local network, set up a firewall rule or SSH tunnel:

```bash
# SSH tunnel from remote machine
ssh -L 3068:localhost:3068 user@server_ip

# Then query locally
sendsketch.sh in=sample.fasta address=http://localhost:3068/sketch k=31
```

### 3.3 Output Formats

| `format=` | Description |
|-----------|-------------|
| `2` | Default human-readable (one header line per query, columnar) |
| `3` | Simple TSV (easy to parse, one line per hit) |
| `4` or `json` | JSON output |

### 3.4 Useful Query Flags

**Output columns:**

| Flag | Description | Needs Index? |
|------|-------------|:---:|
| `records=N` | Top N hits per query | No |
| `printall=t` | Enable all output columns | No |
| `printtaxa=t` | Full taxonomic lineage | No |
| `completeness=t` | Genome completeness estimate | No |
| `printcontam=t` | Contamination estimate | No |
| `printunique=t` | Matches unique to this reference | No |
| `printunique2=t` | Matches unique to this reference's taxon | **Yes** |
| `printunique3=t` | Query kmers unique to this taxon | **Yes** |
| `printcontam2=t` | Contamination using only kmer hits to other taxa | **Yes** |
| `printrefhits=t` | Avg reference sketches hit per shared kmer | **Yes** |
| `score=t` | Composite score (used for sorting) | No |
| `printdepth=t` | Average depth of sketch kmers | No |
| `printgsize=t` | Estimated genome size | No |

**Filtering:**

```bash
# Minimum thresholds
minani=0.90          # Only hits with ≥90% ANI
minwkid=0.01         # Minimum weighted kmer identity
minhits=5            # Minimum shared kmers

# Taxonomic filtering
include=Pseudomonadota           # Only show hits in this clade
exclude=Streptomycetaceae        # Exclude this clade
level=genus                      # Best hit per genus (instead of per species)

# Input handling
persequence              # One sketch per sequence (vs. one per file)
mode=sequence            # Same as persequence
translate=t              # Predict genes and use protein sketches
```

### 3.5 Batch Processing

For many samples, loop over input files:

```bash
# Process all assemblies in a directory
for sample in assemblies/*.fasta; do
    name=$(basename $sample .fasta)
    sendsketch.sh \
        in=$sample \
        address=http://localhost:3068/sketch \
        k=31 records=1 format=3 \
        persequence \
        out=results/${name}_gtdb.txt \
        -Xmx4g
done
```

### 3.6 Quick Taxonomy Lookups

The TaxServer also serves taxonomy queries (no sketch needed):

```bash
# By taxID
curl http://localhost:3068/id/2588536

# By organism name
curl http://localhost:3068/name/methylopumilus_universalis

# By accession
curl http://localhost:3068/accession/NZ_CP040977.1

# Common ancestor of multiple taxa
curl http://localhost:3068/id/ancestor/2588536,1581557
```

---

## Part 4: Our Setup on Ratnakara

| Item | Value |
|------|-------|
| Machine | Ratnakara (Mac Pro, 128 GB RAM) |
| IP | 10.151.50.41 |
| Port | 3068 |
| Database | GTDB R226 (143,614 representative genomes) |
| Sketch file | `/home/nakha/databases/bbsketch_gtdb/gtdb_r226_k31.sketch` (11 GB) |
| Taxonomy tree | `bbsketch_test/taxonomy/tree.taxtree.gz` |
| Java heap | `-Xmx110g` |
| RAM in use | ~40 GB |
| Index | Enabled (1.2 billion unique k-mer hashes) |
| Load time | ~5 minutes |
| Query speed | 759 contigs in 5.5 seconds |

### Start command

```bash
cd /home/nakha/.openclaw/workspace/bbsketch_test
nohup java -da -Xmx110g \
    -cp /home/nakha/miniconda3/opt/bbmap-39.52-0/current/ \
    tax.TaxServer \
    tree=taxonomy/tree.taxtree.gz \
    port=3068 k=31 index \
    /home/nakha/databases/bbsketch_gtdb/gtdb_r226_k31.sketch \
    > /tmp/taxserver_gtdb.log 2>&1 &
```

### Query from any machine on the network

```bash
sendsketch.sh \
    in=my_contigs.fasta \
    address=http://10.151.50.41:3068/sketch \
    k=31 records=3 format=3 \
    persequence printall=t color=f \
    out=results.txt
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| OOM during index build | Increase `-Xmx`. Full GTDB needs ≥80 GB heap; use 110 GB for safety |
| OOM during sketch loading (no index) | The sketch alone expands ~6× in memory. 11 GB sketch needs ~64 GB |
| Server returns help page | Make sure address ends with `/sketch` |
| "invalid URI scheme" error | Include `http://` in the address |
| Slow queries | Enable `index` flag when starting server |
| `unique2`/`contam2` missing | These require `index` to be enabled |
| Server dies when shell closes | Use `nohup ... &` and `disown`, or run via systemd |
