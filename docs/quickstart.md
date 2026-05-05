# Quick Start

This page is for your **first successful analysis after installation**.

- Install first: [Installation](installation.md)
- Configure your client first: [Configuration Guide](advanced/configuration.md)
- If something fails: [Troubleshooting](advanced/troubleshooting.md)

---

## Your First Analysis

Open your MCP client and type:

```text
Load /absolute/path/to/spatial_data.h5ad and show me the tissue structure
```

If you configured ChatSpatial through Docker, use the container path from your mount instead:

```text
Load /data/spatial_data.h5ad and show me the tissue structure
```

For example, `-v /Users/alice/spatial-data:/data:ro` makes host files under `/Users/alice/spatial-data` visible as `/data` inside ChatSpatial. For the full Docker path model, see [Docker / GHCR](docker.md).

Then continue with:

```text
Normalize and cluster the data
```

```text
Identify spatial domains with SpaGCN
```

---

## Common First Prompts

| Task | Say This |
|------|----------|
| Load data | `Load my Visium data from /absolute/path/to/data` |
| Preprocess | `Normalize and cluster the data` |
| Find domains | `Identify spatial domains with SpaGCN` |
| Annotate cells | `Annotate cell types using the reference` |
| Deconvolve | `Estimate cell type proportions` |
| Find spatial genes | `Find spatially variable genes` |
| Visualize | `Show expression of CD3D on the tissue` |

---

## Sample Data

Try these public test files from the sample-data release:

- [card_spatial.h5ad](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data) — pancreatic spatial data
- [card_reference_filtered.h5ad](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data) — reference dataset

The sample-data release is versioned separately from the ChatSpatial package and Docker image.

---

## What Success Looks Like

**After loading:**

```text
Dataset loaded successfully
- ID: spatial_data_abc123
- 3000 spots, 18000 genes
- Platform: Visium
```

**After preprocessing:**

```text
Preprocessing complete
- Filtered to 2800 spots, 2000 HVGs
- Computed 30 PCs, UMAP
- Found 8 clusters (Leiden)
```

Visualizations appear directly in the chat or client UI.

---

## Before You Troubleshoot

Check these basics first:

- Use an **absolute** data path, not `~/...` or `./...`
- If using Docker, use the container path from your mount, such as `/data/sample.h5ad`; full Docker mount rules live in [Docker / GHCR](docker.md)
- Restart your client after MCP configuration changes
- Run preprocessing before most downstream analyses

If that does not fix the problem, go to the [Troubleshooting Guide](advanced/troubleshooting.md).

---

## Next Steps

- [Concepts](concepts.md) — choose the right analysis strategy
- [Examples](examples.md) — complete workflow recipes
- [Methods Reference](advanced/methods-reference.md) — exact tool parameters and defaults
