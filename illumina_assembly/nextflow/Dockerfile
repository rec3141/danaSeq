# ============================================================================
# Dana Illumina Assembly Pipeline - Pipeline Image
# ============================================================================
#
# Thin layer on top of the shared assembly base image (Dockerfile.base)
# which contains conda envs for both nanopore and illumina assembly.
#
# Build:
#   docker build -t danaseq-illumina-assembly .
#
# Run:
#   docker run --user $(id -u):$(id -g) \
#       -v /path/to/reads:/data/input:ro \
#       -v /path/to/output:/data/output \
#       danaseq-illumina-assembly \
#       run /pipeline/main.nf \
#           --input /data/input --outdir /data/output -resume
#
# ============================================================================

FROM ghcr.io/rec3141/danaseq-assembly-base:latest

LABEL maintainer="rec3141@gmail.com"
LABEL description="Dana Illumina Assembly Pipeline — multi-assembler consensus + mapping + depth"
LABEL org.opencontainers.image.source="https://github.com/rec3141/danaSeq"

# Copy pipeline code
COPY . /pipeline/

# Docker config: disable per-process conda (tools already on PATH via wrappers)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

WORKDIR /home/dana

ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
