# syntax=docker/dockerfile:1

# ============================================================
# danaSeq illumina_amplicon container
# Pre-installs the Python engine (papa2, byte-identical to R dada2)
# environments so the pipeline runs without conda env creation.
# ============================================================
FROM condaforge/miniforge3:latest

LABEL org.opencontainers.image.source="https://github.com/rec3141/danaSeq"
LABEL org.opencontainers.image.description="Microscape amplicon sequencing pipeline"

# System dependencies
# - default-jre-headless: Nextflow needs a JVM at both build time (nextflow
#   -version below) and runtime; the miniforge3 base ships no Java.
# - build-essential + zlib1g-dev: papa2 compiles a C lib (`make libpapa2.so`)
#   during its wheel build; it needs make + gcc/g++ and the zlib dev headers.
# - which: GNU which. The base image defines a `which` shell FUNCTION that passes
#   --tty-only --read-alias --read-functions, but Debian 12+ ships debianutils' which,
#   which accepts only [-as]. Every `which java` in nextflow's launcher therefore printed
#   "Illegal option --" ahead of its real output — harmless to runs, but it polluted every
#   invocation and made the version probe record the error as the version.
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps curl default-jre-headless build-essential zlib1g-dev which \
    && rm -rf /var/lib/apt/lists/*

# Node.js (NodeSource LTS) for BUNDLE_VIZ_SITE, which builds the static viz SPA
# that microscape.app hosts (vite build). The viz app's deps are pre-installed
# below so the build runs offline on HPC compute nodes.
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - \
    && apt-get install -y --no-install-recommends nodejs \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow and bake its framework jar into the image.
# The launcher self-downloads the framework jar on first run; on offline HPC
# compute nodes that fails. Caching under $HOME/.nextflow is not portable
# because apptainer overrides HOME at runtime (to the host home), orphaning the
# cache — so pin NXF_HOME to a fixed in-image path and bake the jar there. The
# second, network-free invocation verifies the bake so the build fails loudly
# if the offline path ever regresses.
ENV NXF_HOME=/opt/nextflow
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow -version \
    && chmod -R a+rX "$NXF_HOME" \
    && NXF_OFFLINE=true nextflow -version

# Copy environment specs
COPY envs/ /tmp/envs/

# Upstream HEAD of the remaining git dependency, resolved by CI at build time.
#
# An unpinned URL does not change, so Docker restores the cached layer and pip
# never re-fetches: papa2's version PR merged three minutes before a build and
# the image still shipped a four-month-old commit. Naming the SHA here changes
# the layer key whenever upstream moves, so "unpinned" actually means latest.
# It is recorded rather than enforced: what was really installed is read back
# from pip's direct_url.json by bin/tool_versions.sh.
ARG PAPA2_SHA=main

# Create Python environment (papa2 from git, everything else from conda)
# No `mamba clean` here: /opt/conda/pkgs is a BuildKit cache mount (already
# excluded from the image layer), and cleaning it fails with "Device or
# resource busy" because it's a live mount point.
RUN --mount=type=cache,target=/opt/conda/pkgs \
    echo "upstream: papa2=${PAPA2_SHA}" \
    && mamba env create -y -p /opt/conda/envs/microscape-python \
        -f /tmp/envs/python.yml

# Put both envs on PATH so tools are available without conda activate
ENV PATH="/opt/conda/envs/microscape-python/bin:${PATH}"

# Provenance, baked at build time.
#
# There is no git repository inside the image, so `workflow.commitId` is null for every
# run launched from the .sif — which is how this pipeline actually runs. Three different
# main.nf versions were once in play at once (the one that produced a dataset, the one
# shipped in the container, and the one in git) with no way to tell them apart. CI
# already tags the image with the SHA; this puts it somewhere the pipeline can read.
ARG GIT_SHA=unknown
ARG GIT_REF=unknown
ARG BUILD_DATE=unknown
ENV MICROSCAPE_GIT_SHA=${GIT_SHA} \
    MICROSCAPE_GIT_REF=${GIT_REF} \
    MICROSCAPE_BUILD_DATE=${BUILD_DATE}
LABEL org.opencontainers.image.revision="${GIT_SHA}" \
      org.opencontainers.image.version="${GIT_REF}" \
      org.opencontainers.image.created="${BUILD_DATE}"

# Copy pipeline code
COPY . /pipeline/
RUN chmod +x /pipeline/bin/*.R /pipeline/bin/*.sh /pipeline/entrypoint.sh

# Pre-install the viz app's node_modules so BUNDLE_VIZ_SITE's `npm install
# --prefer-offline` + `vite build` succeed on offline compute nodes.
RUN cd /pipeline/viz && npm ci --no-audit --no-fund

# Create docker.config that disables per-process conda (tools are on PATH)
RUN printf 'conda.enabled = false\n' > /pipeline/docker.config

# Strip caches to reduce image size
RUN find /opt/conda -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true \
    && find /opt/conda -name "*.pyc" -delete 2>/dev/null || true \
    && rm -rf /tmp/envs

# Smoke test
RUN python3 -c "import papa2; print('papa2', papa2.__version__)" \
    && python3 -c "import numpy; print('numpy', numpy.__version__)" \
    && cutadapt --version \
    && /pipeline/bin/auto_trim.py --help > /dev/null

# Writable home for non-root users (HPC)
RUN mkdir -p /home/microscape && chmod 777 /home/microscape
ENV HOME=/home/microscape

WORKDIR /data
ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
