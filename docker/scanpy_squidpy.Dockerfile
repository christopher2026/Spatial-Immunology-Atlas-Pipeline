# Image for the QC, clustering, spatial clustering and spatial statistics modules.
#
# Design notes worth being able to defend in an interview:
#  * micromamba base rather than plain python-slim: the single-cell stack has heavy compiled
#    dependencies (igraph, numba/llvmlite, hdf5) that resolve far more reliably from conda-forge
#    than from source builds on pip.
#  * Every scientific dependency is pinned to an exact version. "It worked in August 2026" has to
#    still be true in a year, and an unpinned `pip install scanpy` guarantees it will not be.
#  * Non-root user (micromamba's default `mambauser`), which is what Nextflow and most HPC/cloud
#    executors expect.
#
# Build:  docker build -f docker/scanpy_squidpy.Dockerfile -t ghcr.io/OWNER/stpipe-scanpy:0.1.0 .
FROM mambaorg/micromamba:2.0-debian12-slim

LABEL org.opencontainers.image.title="stpipe-scanpy-squidpy" \
      org.opencontainers.image.description="scanpy + squidpy + scrublet for single-cell and spatial QC, clustering and spatial statistics" \
      org.opencontainers.image.licenses="MIT"

USER root
# procps supplies `ps`, which Nextflow shells out to in order to collect per-task resource metrics.
# build-essential supplies g++ for dependencies such as annoy when PyPI has no matching wheel for the
# selected Python/platform. Without it pip falls back to a source build and the image build fails.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps ca-certificates build-essential \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/env-scanpy.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml \
 && micromamba clean --all --yes \
 && rm /tmp/env.yml

# Make the env's interpreter the default for non-login shells, so Nextflow `script:` blocks
# (which do not source a profile) find the right python.
ENV PATH="/opt/conda/bin:$PATH" \
    MPLBACKEND=Agg \
    NUMBA_CACHE_DIR=/tmp/numba_cache \
    PYTHONDONTWRITEBYTECODE=1

RUN python -c "import scanpy, squidpy, scrublet, leidenalg; print('scanpy', scanpy.__version__); print('squidpy', squidpy.__version__)"

CMD ["/bin/bash"]
