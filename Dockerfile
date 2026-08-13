# BrainTrace Dockerfile
# =====================
# Build:   docker build -t braintrace:v0.1.12 .
# Run:     docker run -p 8501:8501 braintrace:v0.1.12
#              (Streamlit web interface at http://localhost:8501)
#          docker run --rm --entrypoint braintrace braintrace:v0.1.12 --help
#              (CLI help)
#
# For full reproducibility of manuscript results, mount the Bo2023 external
# data directory and use requirements_reproducible.txt instead.

FROM python:3.11.9-slim-bookworm

LABEL org.opencontainers.image.title="BrainTrace"
LABEL org.opencontainers.image.description="Hierarchical brain-origin candidate ranking from RNA expression profiles"
LABEL org.opencontainers.image.version="v0.1.12"
LABEL org.opencontainers.image.source="https://github.com/wz7717/cfrna-brain-tracing"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN useradd --create-home --shell /bin/bash braintrace
WORKDIR /home/braintrace/app

# Copy dependency files first (layer caching)
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY --chown=braintrace:braintrace . .

# Install the package in development mode for CLI entry
RUN pip install --no-cache-dir -e .

# Create data directories
RUN mkdir -p /home/braintrace/app/data/models

# Switch to non-root user
USER braintrace

# Expose Streamlit default port
EXPOSE 8501

# Default: Streamlit web interface
# Override with `--entrypoint braintrace` for CLI usage
ENTRYPOINT ["streamlit", "run", "app/main.py", "--server.address=0.0.0.0"]
CMD []
