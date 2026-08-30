# Image for the report module: assembles figures and summary tables into a single standalone HTML file.
#
# Intentionally tiny. The report step does no scientific computation — it reads PNGs and TSVs produced
# upstream and renders a Jinja2 template. A slim python base (not micromamba) is right here because
# there are no compiled scientific dependencies to resolve.
#
# Build:  docker build -f docker/report.Dockerfile -t ghcr.io/OWNER/stpipe-report:0.1.0 .
FROM python:3.11.11-slim-bookworm

LABEL org.opencontainers.image.title="stpipe-report" \
      org.opencontainers.image.description="Jinja2 + pandas for assembling the pipeline HTML report" \
      org.opencontainers.image.licenses="MIT"

RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

# Hashless but exact pins. Figures are embedded as base64 so the HTML is a single shareable file
# with no relative-path dependencies.
RUN pip install --no-cache-dir \
      jinja2==3.1.5 \
      pandas==2.2.3 \
      tabulate==0.9.0 \
      markupsafe==3.0.2

ENV PYTHONDONTWRITEBYTECODE=1

RUN python -c "import jinja2, pandas; print('jinja2', jinja2.__version__); print('pandas', pandas.__version__)"

CMD ["/bin/bash"]
