# ============================================================================
# Dana Nanopore Assembly Pipeline - Pipeline Image
# ============================================================================
#
# Thin layer on top of the shared assembly base image (Dockerfile.base)
# which contains conda envs for both nanopore and illumina assembly.
#
# Build:
#   docker build -t danaseq-nanopore-assembly .
#
# Run:
#   docker run --user $(id -u):$(id -g) \
#       -v /path/to/reads:/data/input:ro \
#       -v /path/to/output:/data/output \
#       danaseq-nanopore-assembly \
#       run /pipeline/main.nf \
#           --input /data/input --outdir /data/output -resume
#
# ============================================================================

FROM ghcr.io/rec3141/danaseq-assembly-base:latest

LABEL maintainer="rec3141@gmail.com"
LABEL description="Dana Nanopore Assembly Pipeline — preprocessing, co-assembly, mapping, depth calculation"
LABEL org.opencontainers.image.source="https://github.com/rec3141/danaSeq"

# Copy pipeline code
COPY . /pipeline/

# Compile C tools (fastq_filter, tetramer_freqs)
RUN g++ -O3 -o /pipeline/bin/fastq_filter /pipeline/bin/fastq_filter.cpp -lz -lpthread && \
    gcc -O3 -o /pipeline/bin/tetramer_freqs /pipeline/bin/tetramer_freqs.c -lz -lm

# Docker config: disable per-process conda (tools already on PATH via wrappers)
RUN echo 'conda.enabled = false' > /pipeline/docker.config

WORKDIR /home/dana

ENTRYPOINT ["/pipeline/entrypoint.sh"]
CMD ["run", "/pipeline/main.nf", "--help"]
