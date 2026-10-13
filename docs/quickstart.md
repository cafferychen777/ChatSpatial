# Quick Start

This page is for your **first successful analysis after installation**.

- Install first: [Installation](installation.md)
- Configure your client first: [Configuration Guide](advanced/configuration.md)
- If something fails: [Troubleshooting](advanced/troubleshooting.md)

---

## Download the Example

The first workflow below uses only the **standard installation**. It does not
require optional method extras or a single-cell reference dataset.

Download [card_spatial.h5ad](https://github.com/cafferychen777/ChatSpatial/releases/download/v0.3.0-data/card_spatial.h5ad)
from the [sample-data release](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data).
On macOS or Linux:

```bash
mkdir -p chatspatial-example
cd chatspatial-example
curl -fL -o card_spatial.h5ad \
  https://github.com/cafferychen777/ChatSpatial/releases/download/v0.3.0-data/card_spatial.h5ad
pwd
```

Use the directory printed by `pwd` to construct the absolute path in the prompts
below. The sample contains **428 spatial locations and 25,753 genes**. It has
spatial coordinates but **no histology image**: the plots show expression or
clusters at those coordinates, without an H&E background.

## Your First Analysis

Open your configured MCP client. Replace the example path with your actual
absolute file path, then ask:

```text
Load /absolute/path/to/card_spatial.h5ad as generic spatial data
and plot REG3A expression at the spatial coordinates.
```

For Docker, use the path inside the container instead, such as
`/data/card_spatial.h5ad`. For example,
`-v /Users/alice/chatspatial-example:/data:ro` makes the host directory available
at `/data`. See [Docker / GHCR](docker.md) for the full mount configuration.

Continue in the same conversation:

```text
Preprocess the loaded data with the default QC and normalization settings.
Then compute PCA, neighbors, UMAP, and Leiden clustering.
```

```text
Plot the Leiden clusters at the spatial coordinates, and show a UMAP
colored by the same clusters.
```

This follows `load_data` → `preprocess_data` → `compute_embeddings` →
`visualize_data`. Preprocessing alone does not compute clusters; the embedding
step creates the `leiden` labels used by both plots. Leiden here clusters the
expression-neighbor graph; spatial plotting shows where those clusters occur.

To save the result, choose a writable absolute path:

```text
Export the analyzed dataset to /absolute/path/to/results/pancreas.h5ad.
```

With Docker, export to a mounted output directory such as `/outputs/pancreas.h5ad`.

## What Success Looks Like

- Loading reports 428 locations, 25,753 genes, and available spatial coordinates.
- Preprocessing reports QC statistics and selected highly variable genes.
- Embedding computation reports PCA, neighbors, UMAP, and Leiden clusters.
- The spatial and UMAP plots use the same cluster labels.
- Export returns the absolute path of the saved H5AD file.

Filtered gene counts and cluster counts depend on parameters and runtime
versions. Plots may be displayed inline or returned as saved file paths,
depending on the client.

## Optional: SpaGCN Spatial Domains

SpaGCN is **not included in the standard installation**. First install the
`spatial-domains` extra in the environment that your MCP client actually starts.
For `uvx`, update the registered server arguments to use:

```text
--from chatspatial[spatial-domains] chatspatial server
```

For a persistent environment, run this with that environment active:

```bash
uv pip install 'chatspatial[spatial-domains]'
```

See [Installation](installation.md#optional-method-families) for platform
requirements and [Configuration](advanced/configuration.md) for client syntax.
Restart the client after updating its configuration. Loaded datasets are held
in server memory, so load the sample and preprocess it again after the restart.
Then ask:

```text
Identify 7 spatial domains with SpaGCN using expression and spatial
coordinates only; disable histology because this sample has no tissue image.
Then plot the resulting domains at the spatial coordinates.
```

The histology setting is `spagcn_use_histology=false`. An H&E overlay requires
a dataset that actually includes a tissue image and the corresponding spatial
alignment metadata.

## More Examples

For annotation, deconvolution, and other workflows, see [Examples](examples.md)
for their data requirements and [Installation](installation.md) for method
extras. The optional
[card_reference_filtered.h5ad](https://github.com/cafferychen777/ChatSpatial/releases/download/v0.3.0-data/card_reference_filtered.h5ad)
reference file is available from the same sample-data release, which is versioned
separately from the software package.

---

## Before You Troubleshoot

Check these basics first:

- Use an **absolute** data path, not `~/...` or `./...`
- If using Docker, use the container path from your mount, such as `/data/sample.h5ad`; full Docker mount rules live in [Docker / GHCR](docker.md)
- Restart your client after MCP configuration changes
- Run preprocessing before most downstream analyses
- Run `compute_embeddings` before tools that require PCA, UMAP, clusters, or neighbor graphs

If that does not fix the problem, go to the [Troubleshooting Guide](advanced/troubleshooting.md).

---

## Next Steps

- [Concepts](concepts.md) — choose the right analysis strategy
- [Examples](examples.md) — complete workflow recipes
- [Methods Reference](advanced/methods-reference.md) — exact tool parameters and defaults
