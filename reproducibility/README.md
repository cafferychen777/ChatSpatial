# LLM-free reproduction image

This directory builds the ChatSpatial reproduction image, which is published
with each release under the `-repro` tag. The image re-executes tool calls
that language models recorded during the manuscript experiments, without
calling any LLM, and compares the regenerated numbers with the reference
values in `data/`:

- **DLPFC spatial domains.** The 30 recorded `identify_spatial_domains` calls
  on the DLPFC sections 151673, 151507 and 151669 (10 calls per section) are
  validated against the current parameter schema and executed. The step
  reports the ARI against the expert annotation for each section and pooled.
  By default each distinct validated parameter set is executed once; pass
  `--every-trial` to execute all 30 calls.
- **OSCC CARD deconvolution.** The recorded `deconvolve_data` calls on the
  OSCC sections S1 and S9 are validated again and executed with CARD and the
  Puram HNSCC single-cell reference. The step reports the number of valid calls
  per schema condition, the mean pairwise Pearson correlation of the cell-type
  proportions across valid calls, and the largest deviation from the recorded
  proportion matrices.

The prepared inputs and the recorded tool-call logs are in the
[`repro-inputs-2026-10-08`](https://github.com/cafferychen777/ChatSpatial/releases/tag/repro-inputs-2026-10-08)
release archive. The image downloads this archive (about 480 MB) when it runs.

## Run

```bash
docker run --rm -v "$PWD/out:/outputs" \
  cafferyyang777/chatspatial@sha256:ac3ba1ed8235c3d53ca72005d8ad2556cac007f5ddb932ae3e1d1241fa2cf9d8 \
  chatspatial-reproduce
```

This digest is the v1.5.9 reproduction image (`cafferyyang777/chatspatial:v1.5.9-repro`).
The same image is available as `ghcr.io/cafferychen777/chatspatial@<digest>`.
The table of reference and reproduced values is written to
`out/reproduction/reproduction_summary.csv`, and the command exits with a
non-zero status if any value deviates from its reference by more than the
tolerance (`--tolerance`, default 1e-3). To reuse an unpacked copy of the input
archive, mount it at `/data/chatspatial-repro`.

## Contents

- `Dockerfile` adds R, CARD, MuSiC and SpaGCN to the published ChatSpatial
  image and installs the `chatspatial-reproduce` command.
- `scripts/reproduce_key_results.py` is the reproduction entry point. It
  imports helper functions from four experiment modules next to it
  (`ablation_e2e.py`, `ablation_invocation.py`, `casestudy_reproducibility.py`
  and `dlpfc_benchmark_analysis.py`), which use `paths.py` to locate the
  ChatSpatial source. `scripts/compat/` restores the SciPy sparse-matrix alias
  that SpaGCN 1.2.7 uses.
- `data/` holds the reference values: the DLPFC benchmark results and
  ground-truth labels, and the CARD case-study metrics.

The release workflow (`.github/workflows/docker.yml`) builds this image on top
of the release image digest, and `.github/workflows/reproduce.yml` runs it on
clean Linux x86-64 and arm64 runners.

## Manuscript analysis code

The experiment and analysis scripts for the manuscript, the aggregate result
tables and the supplementary tables are in
[ChatSpatial-Reproducibility](https://github.com/cafferychen777/ChatSpatial-Reproducibility).
