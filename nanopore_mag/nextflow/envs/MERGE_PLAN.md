# Container Image Slimming Plan

Current compressed image: 14.4 GB (47 layers, 35 conda envs, ~32 GB uncompressed)
Target: <10 GB compressed

## Root cause
Each conda env duplicates: Python runtime (30-35 MB), libopenblas (40 MB),
libicudata (32 MB), libstdc++ (21 MB), numpy/scipy/pandas (170+ MB each).
35 envs = massive duplication. Also: tiara pulls full GPU PyTorch (1.4 GB).

## Testing approach
Each proposed merge tested with `mamba create --dry-run` on HPC (2026-03-10).
Results recorded inline. Only merges that solve cleanly are included.

---

## Tested merges: 35 envs → 11 envs

### 1. dana-mag-assembly (flye + mapping + bbmap + metamdbg + myloasm)
- Deps: flye=2.9.5 filtlong samtools>=1.17 minimap2=2.28 coverm bbmap pigz
  metamdbg=1.3.1 myloasm=0.4.0 nextflow openjdk
- DRY-RUN: PASS (141 MB download)
- Current: 2763 MB (5 envs) → Est: ~1500 MB

### 2. dana-mag-binning (binning + semibin + comebin-deps)
- Deps: metabat2=2.17 maxbin2=2.2.7 das_tool=1.1.7 semibin=2.1.0 pytorch-cpu
  python>=3.10 + comebin deps (scanpy networkx numba leidenalg etc)
  + pip: lorbin comebin (forked for modern python)
- DRY-RUN: PASS (561 MB download)
- Comebin fork needed: bump python 3.7→3.10+, pytorch 1.10→cpu, numpy 1.19→modern
  Source code is standard torch/sklearn/pandas — no 3.7-specific syntax (verified)
- Current: 6645 MB (3 envs) → Est: ~3000 MB

### 3. dana-mag-quality (genomad + checkm2 + checkv + whokaryote)
- Deps: genomad checkm2 pyrodigal>=3.7 checkv pytorch-cpu python>=3.10 prodigal
  scikit-learn biopython + pip: whokaryote tiara
- DRY-RUN: PASS for conda deps (430 MB download)
- whokaryote/tiara must be pip-installed: bioconda tiara pins pytorch<1.8
  whokaryote setup.cfg only requires python>=3.8 (no pytorch pin in source)
- Tiara standalone env DELETED — whokaryote bundles it
- Note: for GPU systems, swap pytorch-cpu for pytorch-gpu
- Current: 8078 MB (4 envs, incl tiara) → Est: ~2500 MB

### 4. dana-mag-annotate (bakta + emapper + marferret)
- Deps: bakta>=1.9 eggnog-mapper>=2.1.12 diamond>=2.1.23 pandas
- DRY-RUN: PASS (49 MB download)
- Current: 2342 MB (3 envs) → Est: ~1200 MB

### 5. dana-mag-classify (kaiju + kraken2 + rrna + kofamscan + metaeuk)
- Deps: kaiju=1.10.1 kraken2=2.1.3 barrnap=0.9 vsearch>=2.28 aragorn>=1.2.41
  kofamscan hmmer>=3.3 metaeuk python>=3.8
- DRY-RUN: PASS (8 MB download — mostly cached)
- Current: 2261 MB (5 envs) → Est: ~900 MB

### 6. dana-mag-genomic (defensefinder + macsyfinder + integron + dbcan + islandpath)
- Deps: defense-finder macsyfinder>=2.1 integron_finder dbcan>=4.0 hmmer>=3.3
- DRY-RUN: PASS (96 MB download)
- Note: genomad moved to quality env (shares tensorflow with checkm2)
- Current: 3060 MB (5 envs) → Est: ~1200 MB

### 7. dana-mag-derep (derep + skder + drep)
- Deps: galah=0.4.2 skani=0.3.1 sourmash-minimal skder=1.3.4 drep=3.6.2
- DRY-RUN: PASS (96 MB download)
- Current: 2160 MB (3 envs) → Est: ~1200 MB

### 8. dana-mag-strain (instrain + strainy + floria)
- Deps: strainy=1.2 floria=0.0.2 samtools pysam numpy pandas matplotlib-base
  seaborn-base networkx h5py biopython lmfit tqdm psutil python>=3.10
  + pip: instrain
- DRY-RUN: PASS (199 MB download)
- Current: 2155 MB (3 envs) → Est: ~1500 MB

### 9. dana-mag-pathviz (pathway + viz)
- Deps: python>=3.10 pandas numpy scipy scikit-learn matplotlib seaborn
  biopython plotly glpk umap-learn nodejs>=18
- DRY-RUN: PASS (307 MB download)
- Current: 2257 MB (2 envs) → Est: ~1200 MB

### 10. prokka — DROPPED
- Replaced by bakta (now in annotate env)
- Saves 2741 MB (incl 821 MB HMM databases)

### 11. dana-bbmap — merged into assembly (env #1)

---

## Failed merges (do not attempt)

- **prokka + bakta**: FAIL — Perl 5.26 pin conflicts with bakta's zlib/tbl2asn
- **whokaryote via conda**: FAIL — bioconda tiara pins pytorch<1.8
  Workaround: pip install whokaryote+tiara into quality env

---

## Database externalization
Move databases out of container into download-databases.sh (mounted at runtime):
- integron_finder Functional_annotation HMMs: 83 MB
- tiara neural net models: 96 MB (x2 copies currently)
- semibin pretrained models: 36 MB
- das_tool SCG HMM profiles: 82 MB

Keep in container (small enough / required for solver):
- checkm2 data+models: 33 MB
- genomad data: 32 MB

## Dockerfile cleanup step (after all env installs)
Strip packaging cruft from all envs — saves ~500+ MB even after merges:
```dockerfile
RUN find /opt/conda/envs -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null; \
    rm -rf /opt/conda/envs/*/share/{locale,doc,info,man,cups,icons,gir-1.0,easel} \
           /opt/conda/envs/*/lib/python*/site-packages/{pip,setuptools} \
    && find /opt/conda/envs -name "*.pyc" -delete 2>/dev/null; \
    true
```

## Summary

| # | Merged env | Envs merged | Current (MB) | Est (MB) |
|---|-----------|-------------|-------------|---------|
| 1 | assembly | 5 | 2763 | 1500 |
| 2 | binning | 3 | 6645 | 3000 |
| 3 | quality | 4 | 8078 | 2500 |
| 4 | annotate | 3 | 2342 | 1200 |
| 5 | classify | 5 | 2261 | 900 |
| 6 | genomic | 5 | 3060 | 1200 |
| 7 | derep | 3 | 2160 | 1200 |
| 8 | strain | 3 | 2155 | 1500 |
| 9 | pathviz | 2 | 2257 | 1200 |
| | **TOTAL** | **33→9** | **31721** | **15200** |

Plus savings from:
- Prokka dropped: -2741 MB
- Cruft stripping: -500 MB
- DB externalization: -300 MB

Estimated uncompressed: ~15 GB → compressed: **~7-8 GB** (from 14.4 GB)

## Action items
1. [x] Fix tiara.yml: add pytorch-cpu (done)
2. [ ] Delete tiara env and prokka env from Dockerfile.base
3. [ ] Fork COMEBin: bump to python>=3.10, pytorch-cpu, modern numpy
4. [ ] Create 9 merged env YAML files
5. [ ] Update Dockerfile.base: 35 RUN lines → 9 + cleanup step
6. [ ] Move databases to download-databases.sh (integron, tiara, semibin, das_tool)
7. [ ] Update wrapper scripts in Dockerfile.base
8. [ ] Update Nextflow modules: conda env names in process directives
9. [ ] Remove prokka references from pipeline (main.nf, modules)
10. [ ] Rebuild base image, verify all tools work
11. [ ] Verify download-databases.sh fetches externalized DBs
