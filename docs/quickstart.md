# Quick Start

First time? Complete [Installation](../INSTALLATION.md) first, then come back here.

---

## Your First Analysis

Open Claude and type:

```
Load /path/to/spatial_data.h5ad and show the tissue structure
```

ChatSpatial handles the rest.

---

## What You Can Do

| Task | Say This |
|------|----------|
| Load data | "Load my Visium data from /path/to/data" |
| Preprocess | "Normalize and cluster the data" |
| Find domains | "Identify spatial domains with SpaGCN" |
| Annotate cells | "Annotate cell types using the reference" |
| Deconvolve | "Estimate cell type proportions" |
| Analyze patterns | "Find spatially variable genes" |
| Visualize | "Show expression of CD3D on the tissue" |

---

## Sample Data

Test with our datasets:

- [card_spatial.h5ad](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data) (7.7 MB) - Pancreatic spatial
- [card_reference_filtered.h5ad](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data) (36 MB) - Reference

---

## What Success Looks Like

**After loading:**
```
Dataset loaded successfully
- ID: spatial_data_abc123
- 3000 spots, 18000 genes
- Platform: Visium
```

**After preprocessing:**
```
Preprocessing complete
- Filtered to 2800 spots, 2000 HVGs
- Computed 30 PCs, UMAP
- Found 8 clusters (Leiden)
```

**Visualizations** appear directly in the chat.

---

## Quick Troubleshooting

| Issue | Fix |
|-------|-----|
| Tools not showing | Restart Claude |
| "Dataset not found" | Use absolute path: `/Users/you/data.h5ad` |
| Analysis fails | Run preprocessing first |

More: [Troubleshooting Guide](advanced/troubleshooting.md)

---

## Next Steps

- [Concepts](concepts.md) — Understand what each analysis does
- [Examples](examples.md) — Complete workflows
- [Methods Reference](advanced/methods-reference.md) — All 20 tools with parameters
