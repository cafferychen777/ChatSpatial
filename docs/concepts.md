# Concepts

Understanding spatial transcriptomics analysis.

---

## What is Spatial Transcriptomics?

Spatial transcriptomics measures gene expression while preserving the physical location of cells in tissue. Unlike standard single-cell RNA sequencing, it tells you **where** cells are, not just what they express.

**Key insight**: Location matters. A tumor cell behaves differently depending on whether it's surrounded by immune cells or fibroblasts.

---

## Core Analysis Types

### Spatial Domains

**What it does**: Groups tissue regions based on similar gene expression patterns.

**When to use**: First step after preprocessing. Identifies tissue architecture like tumor regions, immune infiltrates, or tissue layers.

**Method selection**:

| Your Data | Recommended Method |
|-----------|-------------------|
| Visium with H&E image | SpaGCN (uses histology) |
| High-resolution (Xenium, MERFISH) | STAGATE or GraphST |
| Quick exploratory analysis | Leiden clustering |

---

### Cell Type Annotation vs Deconvolution

These two concepts are often confused. Here's the difference:

| | **Annotation** | **Deconvolution** |
|---|----------------|-------------------|
| **Output** | "This spot is T cells" | "This spot is 60% T cells, 30% macrophages, 10% fibroblasts" |
| **Best for** | Single-cell resolution data | Spot-based data (Visium) |
| **Assumption** | One cell type per spot | Multiple cell types per spot |

**Rule of thumb**:
- **Xenium, MERFISH, CosMx**: Use annotation (single-cell resolution)
- **Visium, Slide-seq**: Use deconvolution (multiple cells per spot)

---

### Cell Communication

**What it does**: Identifies which cell types are "talking" to each other through ligand-receptor interactions.

**Key concept**: Cell A expresses a ligand (signal molecule), Cell B expresses the receptor. If they're spatially close, they may be communicating.

**Species matters**: Use the correct database:
- Human: `liana_resource="consensus"`
- Mouse: `liana_resource="mouseconsensus"`

---

### RNA Velocity

**What it does**: Predicts future cell states by comparing spliced vs unspliced RNA.

**Key insight**: If a gene has more unspliced RNA, it's being upregulated. If more spliced, it's being downregulated. This tells you the "direction" cells are moving.

**Requirement**: Your data must have `spliced` and `unspliced` layers (from velocyto, kallisto, or STARsolo).

---

## Choosing Methods

### Deconvolution Methods

| Method | Speed | Accuracy | When to Use |
|--------|-------|----------|-------------|
| **FlashDeconv** | Fast | Good | Default choice, quick exploration |
| **Cell2location** | Slow | Excellent | Final analysis, publication |
| **RCTD** | Fast | Good | R users, batch processing |
| **CARD** | Medium | Good | Need spatial imputation |

**Accuracy vs Speed tradeoff**: Start with FlashDeconv for exploration, run Cell2location for final figures.

---

### Annotation Methods

| Method | Requires | Best For |
|--------|----------|----------|
| **Tangram** | Reference scRNA-seq | Most accurate when reference matches tissue |
| **scANVI** | Reference scRNA-seq | Large datasets, GPU available |
| **CellAssign** | Marker gene list | When you know cell type markers |
| **mLLMCelltype** | Nothing | Quick automated annotation |

---

### Spatial Statistics

| Analysis | Question It Answers |
|----------|-------------------|
| **Moran's I** | "Is this gene spatially clustered?" (global) |
| **Local Moran's I** | "Where are the clusters?" (local hotspots) |
| **Getis-Ord Gi*** | "Where are the high/low expression hotspots?" |
| **Neighborhood enrichment** | "Do these cell types co-localize?" |
| **Co-occurrence** | "How does co-localization change with distance?" |

---

## Understanding Results

### Interpreting Deconvolution

Good deconvolution results show:
- Cell type proportions sum to ~1.0 per spot
- Known tissue structure is visible (e.g., epithelium vs stroma)
- Proportions correlate with histology

Warning signs:
- One cell type dominates everywhere (>80%)
- Proportions don't match expected tissue composition
- Results change dramatically with different methods

---

### Interpreting Spatial Statistics

**Moran's I interpretation**:
- I > 0: Clustered (similar values near each other)
- I ~ 0: Random
- I < 0: Dispersed (dissimilar values near each other)

**p-value**: Tests if pattern is significant vs random.

---

## Common Pitfalls

### 1. Skipping Preprocessing

Most analyses fail because preprocessing wasn't run. Always preprocess first:
```
"Preprocess the data"
```

### 2. Wrong Species Parameter

Cell communication analysis requires correct species:
```
# Human data
species="human"

# Mouse data
species="mouse", liana_resource="mouseconsensus"
```

### 3. Expecting Single-Cell Resolution from Visium

Visium spots contain 1-10 cells. Use deconvolution to estimate proportions, not annotation to assign types.

### 4. Using GPU Methods Without GPU

Methods like Cell2location are 10-100x slower without GPU. Either:
- Set `use_gpu=False` explicitly
- Use CPU-friendly alternatives (FlashDeconv, RCTD)

---

## Workflow Patterns

### Standard Discovery Workflow

```
Load → Preprocess → Domains → Markers → Visualize
```

Best for: Initial exploration of new dataset.

### Reference-Based Workflow

```
Load spatial → Load reference → Preprocess both → Deconvolve → Communicate
```

Best for: When you have matching single-cell reference data.

### Publication Workflow

```
Load → Preprocess → Domains → Deconvolve → Statistics → Communication → Velocity
```

Best for: Comprehensive analysis for publication.

---

## Next Steps

- [Quick Start](quickstart.md) — Get running in 5 minutes
- [Examples](examples.md) — See all analysis types
- [Methods Reference](advanced/methods-reference.md) — Full parameter details
