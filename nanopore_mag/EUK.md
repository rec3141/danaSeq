# Eukaryotic Pipeline — Research & Planning

## Current Implementation (Phase 1)

Two classifiers run on assembled contigs, gated by `--run_eukaryotic`:

| Tool | Method | Input | Output | RAM | Contigs |
|------|--------|-------|--------|-----|---------|
| **Tiara** | Deep learning k-mer NN | Contigs FASTA | Per-contig class + probabilities | <2 GB | >= 3 kb |
| **Whokaryote** | Gene structure random forest | Contigs FASTA + GFF | Per-contig prediction | <1 GB | >= 5 kb |

Tiara uniquely detects organellar sequences (plastid vs mitochondrial). Whokaryote uses biological features (intergenic distance, gene density, gene length) orthogonal to Tiara's k-mer approach. Whokaryote also runs tiara internally as one of its features.

**Upstream fix**: Whokaryote's GFF parser crashed on non-Prodigal GFF files (Bakta, PGAP). We patched `calculate_features.py` to handle standard GFF3 and submitted [PR #13](https://github.com/LottePronk/whokaryote/pull/13) (fixes [issue #12](https://github.com/LottePronk/whokaryote/issues/12)). The patched code is installed in the local conda env.

**Test results** (24,585 contigs from Bakta-annotated marine metagenome):
- Tiara: 24,585 classified (11,830 bacteria, 357 eukarya, 106 archaea, 58 prokarya, 15 organelle, 64 unknown)
- Whokaryote: 12,130 classified (those >= 5 kb)

---

## Planned Additions (Phase 2)

### MetaEuk — Eukaryotic Gene Prediction

**What**: Predicts multi-exon protein-coding genes in eukaryotic contigs using homology to a reference protein database. Handles introns via dynamic programming over six-frame translated fragments. Developed by the Soeding Lab (MMseqs2 group).

**Why**: Prodigal/Bakta cannot predict introns. Eukaryotic genes with introns get fragmented/truncated. MetaEuk recovers full multi-exon genes. Also used internally by BUSCO v5/v6 as the default eukaryotic gene predictor.

**Install**: `mamba install -c bioconda metaeuk` (latest: v7-bba0d80, May 2024, GPL-3.0)

**Algorithm**:
1. Six-frame translation of contigs
2. Homology search against reference protein DB (MMseqs2 engine)
3. Dynamic programming to select optimal exon sets per target-contig-strand
4. Redundancy reduction of overlapping predictions
5. Output: protein FASTA, codon FASTA, GFF

**Output files**:
- `.fas` — predicted protein sequences
- `.codon.fas` — nucleotide coding sequences
- `.gff` — gene structure (non-standard attributes: `Target_ID`, `TCS_ID` instead of `ID`/`Parent`)
- `.headersMap.tsv` — internal ID mapping

**Usage**:
```bash
metaeuk easy-predict \
    contigs.fasta \
    protein_db \
    output_prefix \
    tmp_dir \
    --threads 16 \
    --split-memory-limit 50G \
    -e 100 \
    --metaeuk-eval 0.0001 \
    --metaeuk-tcov 0.6 \
    --min-length 20 \
    --max-intron 10000 \
    --remove-tmp-files 1
```

**Known issues**:
- GFF format is non-standard — may need reformatting for downstream tools
- Temp files can grow very large (use `--remove-tmp-files 1`)
- `--disk-space-limit` reportedly unreliable (GitHub #37)
- Cannot discover genes with no homology to known proteins

#### Database Options

MetaEuk requires a protein reference database. Options ranked by size:

| Database | Type | Compressed | Est. RAM | Download | Marine Relevance |
|----------|------|-----------|----------|----------|-----------------|
| SwissProt eukaryota | Sequence | ~300 MB | 4-8 GB | `metaeuk databases` | Low (curated but small) |
| MMETSP | Sequence | ~3 GB | 8-16 GB | Zenodo | Excellent (marine-specific) |
| OrthoDB v11 eukaryota | Sequence | ~8.5 GB | 15-20 GB | Greifswald mirror | Good (~2,000 euk genomes) |
| OrthoDB v12 eukaryota | Sequence | ~23 GB | 100-180 GB | Greifswald mirror | Best (~5,800 euk genomes, needs >128 GB) |
| UniRef50 eukaryota | Sequence | ~5-10 GB | 8-16 GB | `metaeuk databases` | Moderate (clustered) |
| Full 88M profiles | Profile | ~200 GB | 181+ GB | GWDG server | Maximum (needs --exhaustive-search) |

**MMETSP** (Marine Microbial Eukaryote Transcriptome Sequencing Project): 678 marine protist transcriptomes, ~18.5M proteins. Specifically relevant for marine metagenomics (diatoms, dinoflagellates, ciliates, haptophytes, etc.). Can be combined with OrthoDB via `mmseqs concatdbs`.

**OrthoDB is NOT in `metaeuk databases`**: OrthoDB must be downloaded manually from the Greifswald bioinformatics mirror and converted with `metaeuk createdb`. UniRef/SwissProt databases can be auto-downloaded via `metaeuk databases`.

**RAM management**: MetaEuk supports `--split-memory-limit XG` to cap RAM usage. This flag splits the database into chunks and processes them sequentially, trading runtime for memory. Total RAM usage is approximately `split-memory-limit / 0.8` (the remaining 20% goes to overhead).

**Recommendation for 64 GB RAM**: Use OrthoDB v11 Eukaryota (~8.5 GB compressed, ~15-20 GB as MMseqs2 DB). Fits within `--split-memory-limit 50G` (~62.5 GB total). OrthoDB v12 (~23 GB compressed, ~70-80 GB decompressed) requires >128 GB RAM even with splitting. The full 88M profile database (181 GB estimated) is not feasible.

**Estimated runtime**: 1-4 hours for ~24,000 contigs with 16 cores against a moderate sequence database.

---

### EukCC — Eukaryotic MAG Quality Assessment

**What**: Estimates completeness and contamination of eukaryotic MAGs using clade-specific single-copy marker genes + phylogenetic placement. The eukaryotic equivalent of CheckM2.

**Why**: CheckM2 is prokaryote-only. EukCC fills the gap for eukaryotic bins.

**Install**: `mamba install -c bioconda eukcc` (latest: v2.1.3, Jan 2025, GPL-3.0, actively maintained)

**Algorithm**:
1. Gene prediction with MetaEuk (bundled, no separate DB needed)
2. Phylogenetic placement with EPA-ng onto reference tree
3. Dynamic selection from 477 clade-specific single-copy marker gene sets (built from 734 eukaryotic genomes)
4. HMMER scan against chosen marker set
5. Completeness = fraction of expected markers found; contamination = fraction duplicated

**Database**:
- `eukcc2_db_ver_1.2` — 6 GB compressed, ~12 GB extracted
- URL: `http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.2.tar.gz`
- FTP confirmed accessible
- Includes MetaEuk internally (no separate MetaEuk DB needed for EukCC)

**Usage**:
```bash
# Single bin
eukcc single --db /path/to/eukcc2_db --out outfolder --threads 8 bin.fa

# Batch (all bins in a directory)
eukcc folder --db /path/to/eukcc2_db --out outfolder --threads 8 --suffix .fa bins/
```

**Output**: CSV with columns: `fasta`, `completeness`, `contamination`, `ncbi_lng`

**Resource requirements**: 4-8 GB RAM, minutes per bin with 8 threads.

**Comparison with BUSCO**:
- EukCC: automated marker set selection (477 sets), designed for MAGs
- BUSCO: manual lineage selection, designed for isolate assemblies
- EukCC more accurate for incomplete/contaminated genomes (MAGs)

**Limitations**:
- Works on binned MAGs, not individual unbinned contigs
- Performance varies by clade (fungi: excellent; ciliates, cryptophytes: weak)
- Needs >= ~50% completeness for reliable estimation
- Excludes Bilateria (animals) and vascular plants

---

## Pipeline Integration Plan

```
Contigs (from Flye)
    │
    ├── Tiara (all contigs >= 3 kb)
    ├── Whokaryote (contigs >= 5 kb + GFF)
    │
    ├── MetaEuk (euk gene prediction on all contigs)
    │     └── Output: euk proteins (.faa), CDS (.codon.fas), GFF
    │
    ├── Existing binning pipeline (unchanged)
    │     └── SemiBin2, MetaBAT2, ..., DAS Tool
    │
    └── Quality assessment
          ├── CheckM2 (prokaryotic bins)
          └── EukCC (eukaryotic bins, identified by Tiara/Whokaryote)
```

### Open Questions

1. **MetaEuk database**: Currently using OrthoDB v11 Eukaryota. Could add MMETSP via `mmseqs concatdbs` for better marine coverage.
2. **Which contigs to run MetaEuk on**: All contigs, or only those classified as eukaryotic by Tiara? Running on all is simpler but slower.
3. **EukCC bin selection**: How to route bins to EukCC vs CheckM2? Options:
   - Use Tiara/Whokaryote consensus on constituent contigs (majority vote)
   - Run both on all bins (wasteful but thorough)
   - Run EukCC only on bins with any eukaryotic contigs
4. **MetaEuk GFF reformatting**: Need a script to convert MetaEuk's non-standard GFF to standard GFF3 for downstream tools.

---

## References

- Levy Karin et al. (2020) MetaEuk — sensitive, high-throughput gene discovery and annotation for large-scale eukaryotic metagenomics. *Microbiome* 8:48. doi:10.1186/s40168-020-00808-x
- Saary et al. (2020) Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC. *Genome Biology* 21:244. doi:10.1186/s13059-020-02155-4
- Whokaryote: Pronk & Medema (2022) Whokaryote: distinguishing eukaryotic and prokaryotic contigs in metagenomes based on gene structure. *Microbial Genomics* 8(7). doi:10.1099/mgen.0.000823
- Tiara: Karlicki et al. (2022) Tiara: deep learning-based classification system for eukaryotic sequences. *Bioinformatics* 38(2):344-350. doi:10.1093/bioinformatics/btab672
