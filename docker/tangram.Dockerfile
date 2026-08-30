# Image for the deconvolution module (Tangram).
#
# Separate from the scanpy image on purpose. Tangram pulls in torch, which is by far the largest and
# most version-sensitive dependency in the project; isolating it means a torch upgrade cannot break
# your QC modules, and the QC image stays small enough to pull quickly in CI. This container-per-tool
# boundary is standard nf-core practice and is a direct answer to "why not one big image?".
#
# CPU-only torch: the Visium lymph node sample is ~4k spots against a 73k-cell reference, which maps
# in minutes on CPU. Pulling CUDA wheels would add ~3 GB for no benefit.
#
# Build:  docker build -f docker/tangram.Dockerfile -t ghcr.io/OWNER/stpipe-tangram:0.1.0 .
FROM mambaorg/micromamba:2.0-debian12-slim

LABEL org.opencontainers.image.title="stpipe-tangram" \
      org.opencontainers.image.description="Tangram + CPU torch for mapping single-cell profiles onto Visium spots" \
      org.opencontainers.image.licenses="MIT"

USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER docker/env-tangram.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml \
 && micromamba clean --all --yes \
 && rm /tmp/env.yml

ENV PATH="/opt/conda/bin:$PATH" \
    MPLBACKEND=Agg \
    NUMBA_CACHE_DIR=/tmp/numba_cache \
    PYTHONDONTWRITEBYTECODE=1 \
    OMP_NUM_THREADS=1

# torch respects OMP_NUM_THREADS; Nextflow controls parallelism per task via `cpus`, and letting torch
# also spawn a thread per core causes oversubscription that makes tasks slower, not faster. The module
# overrides this with the task's actual cpu allocation.

RUN python -c "import torch, tangram, scanpy; print('torch', torch.__version__); print('tangram', tangram.__version__)"

CMD ["/bin/bash"]
