# ============================================================================
# Dana MAG Analysis Pipeline - Pipeline Image
# ============================================================================
#
# Thin layer on top of the pre-built base image (Dockerfile.base-slim) which
# contains all conda environments and wrapper scripts.
#
# This image only adds the Nextflow pipeline code, so it rebuilds in seconds
# on every push to main.
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

FROM ghcr.io/rec3141/danaseq-mag-base-slim:latest

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

# COMEBin: the base image has its dependencies in dana-mag-binning and a
# dana-tools wrapper that execs run_comebin.sh from that env, but not the code.
# Install it the way install.sh does: the rec3141 fork under the env's share/,
# a run_comebin.sh that resolves -a/-o/-p against the caller's directory before
# changing into the source tree, and gen_cov_file.sh on PATH.
ARG COMEBIN_COMMIT=8fab9d73268d134b2d33b3d9ea5853a327a5ed70
RUN E=/opt/conda/envs/dana-mag-binning \
    && git init -q "$E/share/COMEBin" \
    && git -C "$E/share/COMEBin" fetch -q --depth 1 https://github.com/rec3141/COMEBin.git "${COMEBIN_COMMIT}" \
    && git -C "$E/share/COMEBin" checkout -q FETCH_HEAD \
    && rm -rf "$E/share/COMEBin/.git" \
    && printf '%s\n' \
        '#!/usr/bin/env bash' \
        'COMEBIN_ROOT="$(dirname "$(dirname "$(readlink -f "$0")")")/share/COMEBin/COMEBin"' \
        'ORIG_CWD="$(pwd)"' \
        'RESOLVED_ARGS=()' \
        'while (( $# )); do' \
        '    case "$1" in' \
        '        -a|-o|-p) RESOLVED_ARGS+=("$1"); shift' \
        '                  RESOLVED_ARGS+=("$(cd "${ORIG_CWD}" && realpath -m "$1")"); shift ;;' \
        '        *) RESOLVED_ARGS+=("$1"); shift ;;' \
        '    esac' \
        'done' \
        'cd "${COMEBIN_ROOT}" && exec bash run_comebin.sh "${RESOLVED_ARGS[@]}"' \
        > "$E/bin/run_comebin.sh" \
    && chmod +x "$E/bin/run_comebin.sh" \
    && test -f "$E/share/COMEBin/COMEBin/run_comebin.sh" \
    && ln -sf "$E/share/COMEBin/COMEBin/scripts/gen_cov_file.sh" "$E/bin/gen_cov_file.sh" \
    && chmod +x "$E/share/COMEBin/COMEBin/scripts/gen_cov_file.sh" \
    && chmod -R a+rX "$E/share/COMEBin"

COPY . /pipeline/

# Pre-install the viz app's node_modules, so a later `vite build` needs nothing
# fetched.
RUN cd /pipeline/viz && npm ci --no-audit --no-fund

# Docker config: disable per-process conda (tools already on PATH via wrappers)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

WORKDIR /home/dana

ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
