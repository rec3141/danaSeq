# ============================================================================
# Dana MAG Analysis Pipeline - Pipeline Image
# ============================================================================
#
# Layer on top of the pre-built base image (Dockerfile.base-slim), which holds
# most conda environments and wrapper scripts. This image adds the tools the
# base lacks (VAMB, Binette, GTDB-Tk and COMEBin envs) and the Nextflow
# pipeline code. Edits to envs/vamb.yml, binette.yml, gtdbtk.yml or
# comebin-cpu.yml do not trigger a rebuild on their own
# (build-mag-analysis.yml ignores envs/**).
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

# Tools the slim base never carried: each has its own env in the modules
# (conda "dana-mag-vamb" etc.), which a container run ignores, so without
# these the binners and GTDB-Tk were missing from PATH.
COPY envs/vamb.yml envs/binette.yml envs/gtdbtk.yml /tmp/envs-extra/
RUN for e in vamb binette gtdbtk; do \
        mamba env create -y -p /opt/conda/envs/dana-mag-$e -f /tmp/envs-extra/$e.yml || exit 1; \
    done \
    && mamba clean -afy \
    && rm -rf /tmp/envs-extra \
    && for pair in vamb:vamb binette:binette gtdbtk:gtdbtk; do \
        env="${pair%%:*}"; tool="${pair#*:}"; \
        test -x "/opt/conda/envs/dana-mag-$env/bin/$tool" || exit 1; \
        printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-%s/bin:$PATH\nexec /opt/conda/envs/dana-mag-%s/bin/%s "$@"\n' "$env" "$env" "$tool" \
            > "/usr/local/bin/dana-tools/$tool" \
        && chmod +x "/usr/local/bin/dana-tools/$tool"; \
    done

# COMEBin in its own env from the bioconda package, which pins the full
# dependency set (checkm-genome, biolib, atomicwrites, tensorboard, ...) that the
# merged binning env lacks. The base image's dana-tools wrapper pointed at a
# run_comebin.sh that was never installed; docker/run_comebin.sh replaces it.
# CheckM marker sets come from CHECKM_DATA_PATH, which the launcher sets from
# --db_dir.
COPY envs/comebin-cpu.yml /tmp/comebin.yml
RUN mamba env create -y -p /opt/conda/envs/dana-mag-comebin -f /tmp/comebin.yml \
    && mamba clean -afy \
    && rm -f /tmp/comebin.yml \
    && test -x /opt/conda/envs/dana-mag-comebin/bin/run_comebin.sh
COPY docker/run_comebin.sh /usr/local/bin/dana-tools/run_comebin.sh
RUN printf '#!/bin/sh\nexport PATH=/opt/conda/envs/dana-mag-comebin/bin:$PATH\nexec /opt/conda/envs/dana-mag-comebin/bin/gen_cov_file.sh "$@"\n' \
        > /usr/local/bin/dana-tools/gen_cov_file.sh \
    && chmod +x /usr/local/bin/dana-tools/run_comebin.sh /usr/local/bin/dana-tools/gen_cov_file.sh

COPY . /pipeline/

# Pre-install the viz app's node_modules, so a later `vite build` needs nothing
# fetched.
RUN cd /pipeline/viz && npm ci --no-audit --no-fund

# Docker config: disable per-process conda (tools already on PATH via wrappers)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

WORKDIR /home/dana

ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
