# ============================================================================
# Illumina MAG Pipeline - Container Image
# ============================================================================
#
# Single-stage build containing all 4 conda environments and the pipeline code.
# Much smaller than the nanopore_mag image (~4 envs vs ~9).
#
# Build:
#   docker build -t danaseq-illumina-mag illumina_mag/nextflow/
#
# Run (as current user):
#   docker run --user $(id -u):$(id -g) \
#       -v /path/to/reads:/data/input:ro \
#       -v /path/to/output:/data/output \
#       -v /path/to/human_ref:/data/human_ref:ro \
#       danaseq-illumina-mag \
#       run /pipeline/main.nf \
#           --input /data/input --outdir /data/output \
#           --human_ref /data/human_ref \
#           -resume
#
# ============================================================================

FROM condaforge/miniforge3:latest

LABEL maintainer="rec3141@gmail.com"
LABEL description="Illumina MAG Pipeline - Metagenomic assembly + binning from Illumina paired-end reads"
LABEL org.opencontainers.image.source="https://github.com/rec3141/danaSeq"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        curl \
    && rm -rf /var/lib/apt/lists/*

# ============================================================================
# Install conda environments
# ============================================================================

COPY envs/ /tmp/envs/

# 1. BBTools + FastQC + samtools + Nextflow
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -y -p /opt/conda/envs/dana-illumina-mag-bbmap -f /tmp/envs/bbmap.yml

# 2. Megahit
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -y -p /opt/conda/envs/dana-illumina-mag-megahit -f /tmp/envs/megahit.yml

# 3. SPAdes
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -y -p /opt/conda/envs/dana-illumina-mag-spades -f /tmp/envs/spades.yml

# 4. MetaBAT2 + binning tools
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -y -p /opt/conda/envs/dana-illumina-mag-binning -f /tmp/envs/binning.yml

# ============================================================================
# Strip packaging cruft
# ============================================================================

RUN find /opt/conda/envs -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null; \
    rm -rf /opt/conda/envs/*/share/{locale,doc,info,man} ; \
    find /opt/conda/envs -name "*.pyc" -delete 2>/dev/null; \
    rm -rf /tmp/envs/ ; \
    true

# ============================================================================
# Create wrapper scripts so tools are callable from PATH
# ============================================================================

RUN mkdir -p /usr/local/bin/dana-tools \
    # --- bbmap env ---
    && for tool in bbduk.sh bbmap.sh bbmerge.sh bbmerge-auto.sh bbnorm.sh \
                   clumpify.sh dedupe.sh filterbytile.sh reformat.sh \
                   removehuman.sh sendsketch.sh statswrapper.sh tadpole.sh \
                   fastqc samtools nextflow java pigz; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-illumina-mag-bbmap/bin:$PATH\nexec /opt/conda/envs/dana-illumina-mag-bbmap/bin/%s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    # --- megahit env ---
    && for tool in megahit; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-illumina-mag-megahit/bin:$PATH\nexec /opt/conda/envs/dana-illumina-mag-megahit/bin/%s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    # --- spades env ---
    && for tool in spades.py metaspades.py; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-illumina-mag-spades/bin:$PATH\nexec /opt/conda/envs/dana-illumina-mag-spades/bin/%s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done \
    # --- binning env ---
    && for tool in metabat2 jgi_summarize_bam_contig_depths; do \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-illumina-mag-binning/bin:$PATH\nexec /opt/conda/envs/dana-illumina-mag-binning/bin/%s "$@"\n' "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done

ENV PATH="/usr/local/bin/dana-tools:${PATH}"

# ============================================================================
# Add pipeline code
# ============================================================================

COPY . /pipeline/

# Docker config: disable per-process conda (tools already on PATH via wrappers)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

# Create non-root user for Nextflow
RUN useradd -m -s /bin/bash dana || true
WORKDIR /home/dana

ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
