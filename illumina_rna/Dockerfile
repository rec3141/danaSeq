# ============================================================================
# Dana Illumina RNA Pipeline - Pipeline Image
# ============================================================================
#
# Thin layer on top of the shared assembly base image (Dockerfile.base) which
# contains conda envs for the danaSeq pipelines. The RNA-specific tools
# (SortMeRNA, subread/featureCounts) layer on top.
#
# Build:
#   docker build -t danaseq-illumina-rna .
#
# Run:
#   docker run --user $(id -u):$(id -g) \
#       -v /path/to/reads:/data/input:ro \
#       -v /path/to/refs:/data/references:ro \
#       -v /path/to/output:/data/output \
#       danaseq-illumina-rna \
#       run /pipeline/main.nf \
#           --input /data/input --references /data/references \
#           --outdir /data/output -resume
# ============================================================================

FROM ghcr.io/rec3141/danaseq-assembly-base:latest

LABEL maintainer="rec3141@gmail.com"
LABEL description="Dana Illumina RNA Pipeline — BBTools QC + SortMeRNA + BBmap-to-references + featureCounts"
LABEL org.opencontainers.image.source="https://github.com/rec3141/danaSeq"

# Add the RNA-seq env on top of the base image's bbmap env
COPY envs/rnaseq.yml /tmp/rnaseq.yml
RUN mamba env create -y -p /opt/conda/envs/dana-illumina-rna-rnaseq -f /tmp/rnaseq.yml \
    && rm -f /tmp/rnaseq.yml \
    && mamba clean -afy

# Copy pipeline code
COPY . /pipeline/

# Docker config: disable per-process conda (tools already on PATH via wrappers)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

WORKDIR /home/dana

ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
