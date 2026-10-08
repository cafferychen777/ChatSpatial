# Reproduction image for the manuscript's key results.
#
# Built on top of the published ChatSpatial image (pass its digest as
# BASE_IMAGE) and adds only what the reproduction needs: R + CARD for the
# deconvolution case study, and the reproduction scripts with their committed
# reference values.
#
#   docker build -f reproducibility/Dockerfile \
#     --build-arg BASE_IMAGE=ghcr.io/cafferychen777/chatspatial@sha256:<digest> \
#     -t chatspatial:repro reproducibility
#
# Run:  docker run --rm -v "$PWD/repro:/outputs" <image> chatspatial-reproduce

ARG BASE_IMAGE=ghcr.io/cafferychen777/chatspatial:latest
FROM ${BASE_IMAGE}

# Default location of the reproduction input archive (prepared inputs and
# recorded tool-call logs); can be overridden at run time.
ARG REPRO_DATA_URL=""
ENV CHATSPATIAL_REPRO_DATA_URL=${REPRO_DATA_URL} \
    CHATSPATIAL_REPRO_DATA_DIR=/data/chatspatial-repro \
    CHATSPATIAL_CODE_DIR=/app \
    OMP_NUM_THREADS=1

# R and the CARD dependencies that Debian packages; the rest come from a CRAN
# snapshot and Bioconductor 3.16, both contemporary with Debian's R 4.2. CARD
# (1.1, the commit used for the manuscript) builds its single-cell reference
# with MuSiC, which is therefore installed from the same commit as well.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        cmake \
        r-base-dev \
        r-cran-rcpp r-cran-rcpparmadillo r-cran-mcmcpack r-cran-fields \
        r-cran-sp r-cran-sf r-cran-dplyr r-cran-tidyr r-cran-rlang \
        r-cran-rann r-cran-ggplot2 r-cran-ggforce r-cran-reshape2 \
        r-cran-rcolorbrewer r-cran-nnls r-cran-pbmcapply r-cran-rcppml \
        r-cran-nmf r-cran-spatstat.random r-cran-gtools r-cran-v8 \
        r-cran-jsonlite r-cran-remotes r-cran-quadprog \
        r-bioc-singlecellexperiment r-bioc-summarizedexperiment \
        r-bioc-biobase r-bioc-limma \
    && rm -rf /var/lib/apt/lists/*
RUN Rscript -e 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2023-06-01", \
                                  BioCsoft = "https://bioconductor.org/packages/3.16/bioc")); \
      install.packages(c("wrMisc", "concaveman", "scatterpie", "ggcorrplot", "TOAST"), Ncpus = 4); \
      remotes::install_github("xuranw/MuSiC@f21fe67f5670d5e9fca0ad7550abaae3423eb59c", \
                              dependencies = FALSE, upgrade = "never"); \
      remotes::install_github("YingMa0107/CARD@2d64b91abb5cdd0c7f576b1c5d4727c84e7c93a0", \
                              dependencies = FALSE, upgrade = "never"); \
      library(MuSiC); library(CARD); cat("CARD", as.character(packageVersion("CARD")), "\n")'

# Python side: rpy2 bridge and CPU-only torch. SpaGCN comes from the base image:
# spagcn-modern >= 1.4 runs the published SpaGCN 1.2.7 model by default and
# reproduces its labels exactly. rpy2 3.6 needs a newer R than the 4.2 that
# Debian bookworm ships.
RUN python -m pip install --index-url https://download.pytorch.org/whl/cpu torch \
    && python -m pip install "rpy2>=3.5.17,<3.6" "anndata2ri>=1.3.0" requests \
    && python -c "import SpaGCN, rpy2.robjects; print('SpaGCN', SpaGCN.__version__ if hasattr(SpaGCN, '__version__') else '', '+ rpy2 ready')"

COPY scripts /opt/chatspatial-repro/scripts
COPY data/dlpfc_benchmark /opt/chatspatial-repro/data/dlpfc_benchmark
COPY data/casestudy_reproducibility /opt/chatspatial-repro/data/casestudy_reproducibility
RUN printf '#!/bin/sh\nexec python /opt/chatspatial-repro/scripts/reproduce_key_results.py "$@"\n' \
        > /usr/local/bin/chatspatial-reproduce \
    && chmod +x /usr/local/bin/chatspatial-reproduce \
    && chatspatial-reproduce --help >/dev/null

WORKDIR /outputs
ENTRYPOINT []
CMD ["chatspatial-reproduce"]
