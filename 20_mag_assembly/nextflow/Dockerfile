# ============================================================================
# Dana MAG Assembly Pipeline - Docker Container
# ============================================================================
#
# Self-contained image with all MAG assembly tools pre-installed in isolated
# conda environments. Uses Miniforge for conda package management.
#
# Twelve environments match the local install (install.sh):
#   dana-mag-flye     - Flye, Filtlong, Nextflow, OpenJDK
#   dana-mag-mapping  - minimap2, samtools, CoverM
#   dana-mag-semibin  - SemiBin2, LorBin (CPU-only PyTorch)
#   dana-mag-comebin  - COMEBin (CPU-only PyTorch)
#   dana-mag-binning  - MetaBAT2, MaxBin2, DAS_Tool
#   dana-mag-prokka   - Prokka gene annotation
#   dana-mag-genomad  - geNomad (virus + plasmid + provirus detection)
#   dana-mag-checkv   - CheckV (viral quality assessment)
#   dana-mag-integron   - IntegronFinder (integron + gene cassette detection)
#   dana-mag-islandpath   - IslandPath-DIMOB (genomic island detection)
#   dana-mag-macsyfinder  - MacSyFinder v2 (secretion systems + conjugation)
#   dana-mag-checkm2      - CheckM2 (MAG quality assessment)
#   dana-bbmap        - BBMap (optional dedupe)
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
RUN mamba env create -y -p /opt/conda/envs/dana-mag-genomad  -f /tmp/envs/genomad.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-checkv   -f /tmp/envs/checkv.yml       && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-integron   -f /tmp/envs/integron.yml     && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-islandpath   -f /tmp/envs/islandpath.yml    && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-macsyfinder -f /tmp/envs/macsyfinder.yml  && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-mag-checkm2     -f /tmp/envs/checkm2.yml      && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-bbmap        -f /tmp/envs/bbmap.yml         && mamba clean -afy

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
    && for tool in islandpath; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-islandpath/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in macsyfinder macsydata; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-macsyfinder/bin:$PATH\nexec %s "$@"\n' "$tool" \
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
