# ============================================================================
# Dana Pipeline - Docker Container
# ============================================================================
#
# Self-contained image with all pipeline tools pre-installed.
# Uses Miniforge for conda package management.
#
# Build:
#   docker build -t danaseq-realtime .
#
# Run:
#   docker run -v /path/to/data:/data -v /path/to/krakendb:/kraken_db \
#       danaseq-realtime \
#       nextflow run /pipeline/main.nf \
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
        gawk \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -s https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow -version

# ============================================================================
# Install bioinformatics tools into a single conda environment
# Inside Docker there's no host env to protect, so one env is fine.
# ============================================================================

RUN mamba create -y -n dana \
        -c conda-forge -c bioconda \
        bbmap \
        filtlong \
        kraken2 \
        prokka \
        hmmer \
        perl \
        r-base \
        r-dbi \
        r-duckdb \
        r-readr \
    && mamba clean -afy

# Make the dana environment the default for all commands
ENV PATH="/opt/conda/envs/dana/bin:${PATH}"

# Download tetramer_freqs_esom.pl (not available via conda)
RUN curl -fsSL \
        https://raw.githubusercontent.com/tetramerFreqs/Binning/master/tetramer_freqs_esom.pl \
        -o /opt/conda/envs/dana/bin/tetramer_freqs_esom.pl \
    && chmod +x /opt/conda/envs/dana/bin/tetramer_freqs_esom.pl

# ============================================================================
# Copy pipeline files
# ============================================================================

COPY . /pipeline/

# Docker config: disable per-process conda (tools already on PATH)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

WORKDIR /pipeline

ENTRYPOINT ["nextflow", "-c", "/pipeline/docker.config"]
CMD ["run", "/pipeline/main.nf", "--help"]
