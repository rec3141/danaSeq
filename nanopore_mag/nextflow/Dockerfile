# ============================================================================
# Dana MAG Assembly Pipeline - Pipeline Image
# ============================================================================
#
# Thin layer on top of the pre-built base image (Dockerfile.base) which
# contains all 35 conda environments and wrapper scripts.
#
# This image only adds the Nextflow pipeline code, so it rebuilds in seconds
# on every push to main.
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

FROM ghcr.io/rec3141/danaseq-mag-base:latest

LABEL maintainer="rec3141@gmail.com"
LABEL description="Dana MAG Assembly Pipeline - Metagenome-assembled genome reconstruction"
LABEL org.opencontainers.image.source="https://github.com/rec3141/danaSeq"
LABEL org.opencontainers.image.description="danaSeq MAG Assembly Pipeline — Metagenome-assembled genome reconstruction from Oxford Nanopore long reads"

# Copy pipeline code
COPY . /pipeline/

# Docker config: disable per-process conda (tools already on PATH via wrappers)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

WORKDIR /home/dana

ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
