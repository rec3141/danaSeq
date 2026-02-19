# ============================================================================
# Dana Pipeline - Docker Container
# ============================================================================
#
# Self-contained image with all pipeline tools pre-installed in isolated
# conda environments. Uses Miniforge for conda package management.
#
# Four environments are needed due to dependency conflicts:
#   dana-bbmap   - BBMap (samtools/libdeflate conflicts with R)
#   dana-prokka  - Prokka (BioPerl pins perl to 5.26)
#   dana-bakta   - Bakta (modern alternative to Prokka)
#   dana-tools   - Everything else (nextflow, openjdk, filtlong, kraken2, hmmer, R/DuckDB, perl)
#
# Build:
#   docker build -t danaseq-realtime .
#
# Run (as current user):
#   docker run --user $(id -u):$(id -g) \
#       -v /path/to/data:/data -v /path/to/krakendb:/kraken_db \
#       danaseq-realtime \
#       run /pipeline/main.nf \
#           --input /data \
#           --kraken_db /kraken_db \
#           --run_kraken --run_prokka \
#           -resume
#
# ============================================================================

FROM condaforge/miniforge3:latest AS base

LABEL maintainer="rec3141@gmail.com"
LABEL description="Dana Pipeline - Real-time Nanopore metagenomic processing"
LABEL org.opencontainers.image.source="https://github.com/rec3141/danaSeq"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        curl \
    && rm -rf /var/lib/apt/lists/*

# ============================================================================
# Install bioinformatics tools from pinned YAML files
# ============================================================================

COPY envs/ /tmp/envs/

RUN mamba env create -y -p /opt/conda/envs/dana-bbmap    -f /tmp/envs/bbmap.yml    && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-prokka   -f /tmp/envs/prokka.yml   && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-bakta    -f /tmp/envs/bakta.yml    && mamba clean -afy
RUN mamba env create -y -p /opt/conda/envs/dana-tools    -f /tmp/envs/tools.yml    && mamba clean -afy

# Download tetramer_freqs_esom.pl (not available via conda)
RUN curl -fsSL \
        https://raw.githubusercontent.com/tetramerFreqs/Binning/master/tetramer_freqs_esom.pl \
        -o /opt/conda/envs/dana-tools/bin/tetramer_freqs_esom.pl \
    && chmod +x /opt/conda/envs/dana-tools/bin/tetramer_freqs_esom.pl

# ============================================================================
# Create wrapper scripts so tools are callable from PATH
# ============================================================================

RUN mkdir -p /usr/local/bin/dana-tools \
    && for tool in bbduk.sh reformat.sh sendsketch.sh; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-bbmap/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in prokka; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-prokka/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in bakta; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-bakta/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    && for tool in nextflow filtlong kraken2 gawk hmmsearch perl tetramer_freqs_esom.pl Rscript R; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-tools/bin:$PATH\nexec %s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done

ENV PATH="/usr/local/bin/dana-tools:${PATH}"
# No JAVA_HOME — each conda env has its own JDK; BBMap's javasetup.sh finds
# java from PATH which resolves to the correct env-specific binary.

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
