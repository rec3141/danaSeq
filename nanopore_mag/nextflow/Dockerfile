# ============================================================================
# Dana MAG Assembly Pipeline - Docker Container
# ============================================================================
#
# Self-contained image with all MAG assembly tools pre-installed in isolated
# conda environments. Uses Miniforge for conda package management.
#
# Thirty-five environments match the local install (install.sh):
#   dana-mag-flye          - Flye, Filtlong, Nextflow, OpenJDK
#   dana-mag-mapping       - minimap2, samtools, CoverM
#   dana-mag-semibin       - SemiBin2, LorBin (CPU-only PyTorch)
#   dana-mag-comebin       - COMEBin (CPU-only PyTorch)
#   dana-mag-binning       - MetaBAT2, MaxBin2, DAS_Tool
#   dana-mag-prokka        - Prokka gene annotation
#   dana-mag-bakta         - Bakta gene annotation (alternative to Prokka)
#   dana-mag-kaiju         - Kaiju (protein-level taxonomy)
#   dana-mag-kofamscan     - KofamScan (KEGG Orthology via adaptive HMM thresholds)
#   dana-mag-emapper       - eggNOG-mapper (COG/GO/EC/KEGG/Pfam annotation)
#   dana-mag-dbcan         - run_dbcan / dbCAN3 (CAZyme annotation)
#   dana-mag-genomad       - geNomad (virus + plasmid + provirus detection)
#   dana-mag-checkv        - CheckV (viral quality assessment)
#   dana-mag-integron      - IntegronFinder (integron + gene cassette detection)
#   dana-mag-islandpath    - IslandPath-DIMOB (genomic island detection)
#   dana-mag-macsyfinder   - MacSyFinder v2 (secretion systems + conjugation)
#   dana-mag-defensefinder - DefenseFinder (anti-phage defense systems)
#   dana-mag-checkm2       - CheckM2 (MAG quality assessment)
#   dana-mag-tiara         - Tiara (deep learning eukaryotic contig classification)
#   dana-mag-whokaryote    - Whokaryote (gene structure-based eukaryotic classification)
#   dana-mag-metaeuk       - MetaEuk (eukaryotic gene prediction, multi-exon)
#   dana-mag-kraken2       - Kraken2 (k-mer contig-level taxonomy)
#   dana-mag-rrna          - barrnap + vsearch (rRNA gene detection + SILVA classification)
#   dana-mag-pathway       - MinPath + KEGG-Decoder (parsimony pathways + biogeochemical heatmaps)
#   dana-mag-marferret     - DIAMOND + Python/pandas (MarFERReT marine eukaryotic classification)
#   dana-mag-viz           - Node.js + Python/pandas/scipy (viz dashboard preprocessing + build)
#   dana-bbmap             - BBMap (optional dedupe)
#   dana-mag-metamdbg      - metaMDBG / nanoMDBG (de Bruijn graph ONT metagenome assembler)
#   dana-mag-myloasm       - myloasm (high-resolution ONT strain assembler)
#   dana-mag-derep         - galah, skani, sourmash (fast MAG dereplication + ANI + sketching)
#   dana-mag-drep          - dRep (established MAG dereplication)
#   dana-mag-instrain      - InStrain (strain-level population genomics)
#   dana-mag-strainy       - Strainy (strain-aware assembly from long reads)
#   dana-mag-floria        - Floria (strain-aware phasing from long reads)
#   dana-mag-skder         - skDER (fast skani-based dereplication)
#
# Build:
#   docker build -t danaseq-mag .
#
# Run (as current user):
#   docker run --user $(id -u):$(id -g) \
#       -v /path/to/reads:/data/input:ro \
#       -v /path/to/output:/data/output \
#       danaseq-mag \
#       run /pipeline/main.nf \
#           --input /data/input --outdir /data/output \
#           -resume
#
# ============================================================================

FROM condaforge/miniforge3:latest AS base

LABEL maintainer="rec3141@gmail.com"
LABEL description="Dana MAG Assembly Pipeline - Metagenome-assembled genome reconstruction"
LABEL org.opencontainers.image.source="https://github.com/rec3141/danaSeq"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        curl \
        git \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

# ============================================================================
# Install bioinformatics tools from pinned YAML files
# ============================================================================

COPY envs/ /tmp/envs/

RUN mamba env create -y -p /opt/conda/envs/dana-mag-flye     -f /tmp/envs/flye.yml        && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-mapping  -f /tmp/envs/mapping.yml      && mamba clean -afy
# Use CPU-only PyTorch in Docker (~500MB vs ~5GB for GPU variant)
RUN mamba env create -y -p /opt/conda/envs/dana-mag-semibin  -f /tmp/envs/semibin-cpu.yml  && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-comebin  -f /tmp/envs/comebin-cpu.yml  && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-binning  -f /tmp/envs/binning.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-prokka   -f /tmp/envs/prokka.yml       && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-bakta    -f /tmp/envs/bakta.yml        && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-kaiju    -f /tmp/envs/kaiju.yml        && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-kofamscan -f /tmp/envs/kofamscan.yml    && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-emapper  -f /tmp/envs/emapper.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-dbcan    -f /tmp/envs/dbcan.yml         && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-genomad  -f /tmp/envs/genomad.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-checkv   -f /tmp/envs/checkv.yml       && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-integron   -f /tmp/envs/integron.yml     && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-islandpath   -f /tmp/envs/islandpath.yml    && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-macsyfinder    -f /tmp/envs/macsyfinder.yml    && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-defensefinder -f /tmp/envs/defensefinder.yml && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-checkm2       -f /tmp/envs/checkm2.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-bbmap        -f /tmp/envs/bbmap.yml         && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-tiara       -f /tmp/envs/tiara.yml        && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-whokaryote  -f /tmp/envs/whokaryote.yml   && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-metaeuk     -f /tmp/envs/metaeuk.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-kraken2     -f /tmp/envs/kraken2.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-rrna        -f /tmp/envs/rrna.yml         && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-pathway     -f /tmp/envs/pathway.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-marferret   -f /tmp/envs/marferret.yml    && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-viz         -f /tmp/envs/viz.yml          && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-metamdbg    -f /tmp/envs/metamdbg.yml     && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-myloasm     -f /tmp/envs/myloasm.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-derep       -f /tmp/envs/derep.yml        && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-drep        -f /tmp/envs/drep.yml         && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-instrain    -f /tmp/envs/instrain.yml     && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-strainy     -f /tmp/envs/strainy.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-floria      -f /tmp/envs/floria.yml       && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-skder       -f /tmp/envs/skder.yml        && mamba clean -afy

# Nextflow + OpenJDK are in the flye env (via flye.yml)

# ============================================================================
# Create wrapper scripts so tools are callable from PATH
# ============================================================================

RUN mkdir -p /usr/local/bin/dana-tools \
    && for tool in flye filtlong nextflow java; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-flye/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in minimap2 samtools coverm; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-mapping/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in SemiBin2 LorBin; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-semibin/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in run_comebin.sh gen_cov_file.sh; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-comebin/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in metabat2 run_MaxBin.pl DAS_Tool jgi_summarize_bam_contig_depths; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-binning/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in prokka; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-prokka/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in bakta; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-bakta/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in exec_annotation; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-kofamscan/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in emapper.py download_eggnog_data.py; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-emapper/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in run_dbcan dbcan_build; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-dbcan/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in genomad; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-genomad/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in checkv; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-checkv/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in integron_finder; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-integron/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in kaiju kaiju-addTaxonNames kaiju2table; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-kaiju/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in hmmscan; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-islandpath/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in macsyfinder macsydata; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-macsyfinder/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in defense-finder; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-defensefinder/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in checkm2; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-checkm2/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in bbduk.sh dedupe.sh; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-bbmap/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in tiara; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-tiara/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in whokaryote.py; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-whokaryote/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in metaeuk; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-metaeuk/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in kraken2; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-kraken2/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in barrnap vsearch; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-rrna/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in KEGG-decoder; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-pathway/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in diamond; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-marferret/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in node npm python3; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-viz/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in metaMDBG; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-metamdbg/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in myloasm; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-myloasm/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in galah skani sourmash; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-derep/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in dRep; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-drep/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in inStrain; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-instrain/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in strainy; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-strainy/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in floria; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-floria/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in skder; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-skder/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done

ENV PATH="/usr/local/bin/dana-tools:${PATH}"

# ============================================================================
# Writable home for non-root users (supports docker run --user)
# ============================================================================

RUN mkdir -p /home/dana && chmod 777 /home/dana
ENV HOME=/home/dana

# ============================================================================
# Copy pipeline files
# ============================================================================

COPY . /pipeline/

# Docker config: disable per-process conda (tools already on PATH via wrappers)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

WORKDIR /home/dana

ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
