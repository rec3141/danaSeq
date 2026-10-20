# ============================================================================
# Dana MAG Analysis Pipeline - Pipeline Image
# ============================================================================
#
# Thin layer on top of the pre-built base image (Dockerfile.base-slim), which
# contains all conda environments and wrapper scripts. This image adds Node for
# the viz build and the Nextflow pipeline code.
#
# Build:
#   docker build -t danaseq-mag-analysis .
#
# Run (as current user):
#   docker run --user $(id -u):$(id -g) \
#       -v /path/to/assembly.fasta:/data/assembly.fasta:ro \
#       -v /path/to/depths.txt:/data/depths.txt:ro \
#       -v /path/to/output:/data/output \
#       danaseq-mag-analysis \
#       run /pipeline/main.nf \
#           --assembly /data/assembly.fasta --depths /data/depths.txt \
#           --outdir /data/output -resume
#
# ============================================================================

# danaseq-mag-base-slim-gpu (CUDA PyTorch in the binning envs) builds the
# danaseq-mag-analysis-gpu variant.
ARG BASE_IMAGE=ghcr.io/rec3141/danaseq-mag-base-slim:latest
FROM ${BASE_IMAGE}

LABEL maintainer="rec3141@gmail.com"
LABEL description="Dana MAG Analysis Pipeline — Technology-agnostic downstream MAG analysis"
LABEL org.opencontainers.image.source="https://github.com/rec3141/danaSeq"
LABEL org.opencontainers.image.description="danaSeq MAG Analysis Pipeline — binning, annotation, taxonomy, metabolism, MGEs, eukaryotes, ecosystem services, viz"

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

# Node, for building the viz SPA out of the image. The site is rebuilt from the
# image whenever the front end moves (omc-refresh-sites.py), so the image has to
# be able to build it without reaching the network — which is also what lets it
# build on an offline compute node.
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - \
    && apt-get install -y --no-install-recommends nodejs \
    && rm -rf /var/lib/apt/lists/*

# Bakta from our fork: 1.12.2 plus a fix that lets pseudogene detection finish
# on whole co-assemblies (it ran for days on one core). The base image strips pip
# from its envs, so the base conda's pip installs into the env's interpreter.
ARG BAKTA_COMMIT=fc13deb142f85b6291ddce388b799733e5cdd8c6
RUN SP=/opt/conda/envs/dana-mag-annotate/lib/python3.11/site-packages \
    && rm -rf "$SP"/bakta "$SP"/bakta-*.dist-info \
    && /opt/conda/bin/pip --python /opt/conda/envs/dana-mag-annotate/bin/python install \
        --no-deps --no-cache-dir --root-user-action=ignore \
        "bakta @ git+https://github.com/rec3141/bakta.git@${BAKTA_COMMIT}" \
    && grep -q uniref90_counts "$SP"/bakta/features/cds.py

COPY . /pipeline/

# Pre-install the viz app's node_modules, so a later `vite build` needs nothing
# fetched.
RUN cd /pipeline/viz && npm ci --no-audit --no-fund

# Docker config: disable per-process conda (tools already on PATH via wrappers)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

WORKDIR /home/dana

ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
