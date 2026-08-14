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

# Dependencies resolve from the metadata alone, so they are installed against a
# placeholder package first. Copying the source before this step made every
# edit to a single .py file rebuild the whole dependency tree — the layer that
# compiles scipy, h5py and friends, and the one step that is slow on arm64.
COPY pyproject.toml README.md LICENSE ./
RUN mkdir -p chatspatial \
    && printf '__version__ = "0.0.0"\n' > chatspatial/__init__.py \
    && python -m pip install --upgrade pip setuptools wheel \
    && python -m pip install .

# Now the real source, installed over the placeholder without touching the
# dependency layer above.
COPY chatspatial ./chatspatial
RUN python -m pip install --no-deps --force-reinstall .

RUN mkdir -p /data /outputs /workspace \
    && python -c "import chatspatial; print(f'ChatSpatial {chatspatial.__version__} ready')" \
    && python -m chatspatial --version \
    && python -m chatspatial server --help >/tmp/chatspatial-server-help.txt

WORKDIR /workspace
VOLUME ["/data", "/outputs"]
EXPOSE 8000

ENTRYPOINT ["python", "-m", "chatspatial"]
CMD ["server", "--transport", "stdio"]
