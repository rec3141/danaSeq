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
# Provenance for runs launched from the .sif. There is no git repository inside the
# image, so `workflow.commitId` is null and the run manifest would otherwise record
# no build at all. CI already tags the image with the SHA; this puts it somewhere
# the pipeline can read.
ARG GIT_SHA=unknown
ARG GIT_REF=unknown
ARG BUILD_DATE=unknown
ENV DANASEQ_GIT_SHA=${GIT_SHA} \
    DANASEQ_GIT_REF=${GIT_REF} \
    DANASEQ_BUILD_DATE=${BUILD_DATE}
LABEL org.opencontainers.image.revision="${GIT_SHA}" \
      org.opencontainers.image.version="${GIT_REF}" \
      org.opencontainers.image.created="${BUILD_DATE}"

COPY . /pipeline/

# Docker config: disable per-process conda (tools already on PATH via wrappers)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

WORKDIR /home/dana

ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
