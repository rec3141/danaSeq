# GTDB sketch database for the sendsketch server

`SENDSKETCH` in nanopore_live (and `SENDSKETCH_CLASSIFY` in mag_analysis) query
a local BBTools `TaxServer`. These scripts build its database from a GTDB
release, so that the lineage it returns is the GTDB lineage, the same taxonomy
GTDB-Tk assigns.

| Script | Does |
|---|---|
| `build_sketch.sh` | Whole build: taxonomy tree, unpack genomes, sketch in batches, stamp taxids, concatenate |
| `build_gtdb_taxdump.py` | GTDB taxonomy tables -> NCBI-style `nodes.dmp`/`names.dmp` + `accession2taxid.tsv` |
| `set_sketch_taxids.py` | Sets each sketch's `ID:` to its GTDB species taxid and `NM0:` to its accession, matched by the `FN:` file name |

## Build

Inputs from `https://data.gtdb.ecogenomic.org/releases/release<NNN>/`:
`genomic_files_reps/gtdb_genomes_reps_r<NNN>.tar.gz` (192 GB for R232) and
`bac120_taxonomy_r<NNN>.tsv.gz` + `ar53_taxonomy_r<NNN>.tsv.gz`. Check them
against the release's `MD5SUM.txt`.

```bash
build_sketch.sh gtdb_genomes_reps_r232.tar.gz <dir with taxonomy tables> <workdir>
```

The genomes unpack to about the size of the tarball. `BATCH` and `JOBS` set the
batch size and concurrency (default 2000 genomes, 6 runs of 5 threads / 12 GB).

## Serve

```bash
java -da -Xmx<heap> -cp <bbtools>/current tax.TaxServer port=3068 k=32,24 \
    index=auto tree=gtdb_r232.taxtree.gz gtdb_r232_k2432.sketch
```

`k=32,24` is dual-k, which stock BBTools serves correctly; a single-k server
needs the `SketchTool.loadSketchesFromString` fix (bbushnell/BBTools#22).
The indexed server holds roughly 40 GB for 143,614 sketches, growing in steps
as its hash tables double; leave the heap well above that.

## Output

`printtaxa=t` returns `d:Bacteria;p:Bacillota;...;s:Natranaerobius thermophilus`,
and `taxName` is the GTDB species name. The pipeline parsers turn this into
`d__...;s__...`.

## Constraints the tree has to meet

- **Rank names.** The top rank is `domain`. BBTools treats `superkingdom` as a
  separate level printed with a different token, and the pipeline parsers keep
  only `d: p: c: o: f: g: s:`, so a `superkingdom` tree leaves every `d__` empty.
- **Reserved taxids.** BBTools looks up Bacteria as taxid 2 and Archaea as 2157,
  and treats ids such as 2759 (Eukaryota) and 10239 (Viruses) specially. The
  domains keep those two ids and every other GTDB taxon is numbered from
  3,000,000. Ids are stable for a rerun on the same release, not across releases.
- **`nodes.dmp` columns.** BBTools' parser needs at least four `|`-separated
  fields per row, so rows carry NCBI's empty trailing columns.

The sketch's `NM0:` is the genome accession, so a hit can be joined to the
release's `*_metadata_r<NNN>.tsv.gz` (which has `ncbi_taxid` and
`ncbi_taxonomy`) when an NCBI lineage is also wanted.
