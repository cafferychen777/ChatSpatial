FROM python:3.12-slim-bookworm

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    TQDM_DISABLE=1 \
    MPLBACKEND=Agg \
    CHATSPATIAL_OUTPUT_DIR=/outputs

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        git \
        libgdal-dev \
        libgeos-dev \
        libgl1 \
        libhdf5-dev \
        libspatialindex-dev \
        libxml2-dev \
        libxslt1-dev \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml README.md LICENSE ./
COPY chatspatial ./chatspatial

RUN python -m pip install --upgrade pip setuptools wheel \
    && python -m pip install .

RUN mkdir -p /data /outputs /workspace \
    && python -c "import chatspatial; print(f'ChatSpatial {chatspatial.__version__} ready')" \
    && python -m chatspatial --version \
    && python -m chatspatial server --help >/tmp/chatspatial-server-help.txt

WORKDIR /workspace
VOLUME ["/data", "/outputs"]
EXPOSE 8000

ENTRYPOINT ["python", "-m", "chatspatial"]
CMD ["server", "--transport", "stdio"]
