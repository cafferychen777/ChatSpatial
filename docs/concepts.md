# Concepts

Understanding spatial transcriptomics analysis.

---

## What is Spatial Transcriptomics?

Spatial transcriptomics measures gene expression while preserving the physical location of cells in tissue. Unlike standard single-cell RNA sequencing, it tells you **where** cells are, not just what they express.

**Key insight**: Location matters. A tumor cell behaves differently depending on whether it's surrounded by immune cells or fibroblasts.

A spatial transcriptomics dataset usually combines four pieces of information:

1. **Expression** — genes measured in spots or cells.
2. **Coordinates** — where each spot or cell sits in the tissue.
3. **Metadata** — sample, condition, cluster, cell type, or patient labels.
4. **Histology image** — optional tissue morphology, most common in Visium-like data.

Most analyses ask one of six biological questions:

| Question | Analysis family |
|----------|-----------------|
| What tissue regions or niches exist? | Spatial domains |
| What cells are here? | Annotation or deconvolution |
| Which genes have spatial structure? | Spatially variable genes or spatial statistics |
| Which cells may communicate in space? | Cell communication |
| What differs between conditions or regions? | Differential or condition comparison |
| How do cell states change? | Trajectory or RNA velocity |

Choose methods in this order: first define the biological question, then check
platform and data requirements, then choose a method that fits your compute
budget and validation needs.

---

## Core Analysis Types

### Spatial Domains

**What it does**: Groups tissue regions based on similar gene expression patterns.

**When to use**: First step after preprocessing. Identifies tissue architecture like tumor regions, immune infiltrates, or tissue layers.

**How to choose**:

- If your data includes informative histology images, prefer a histology-aware method.
- If your data is high-resolution single-cell spatial data, prefer graph or deep-learning domain methods.
- If you want a quick baseline, use clustering-first exploration.

See [Methods Reference](advanced/methods-reference.md) for the full supported domain methods and defaults.

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

### Deconvolution

Use deconvolution when a spot contains multiple cell types and you want proportions rather than a single label.

**How to choose**:
- start with a fast method for exploration
- move to a slower, stronger method for final figures
- prefer methods with explicit spatial modeling if tissue structure matters to the question

See [Methods Reference](advanced/methods-reference.md) for the full deconvolution method list, defaults, and requirements.

---

### Annotation

Use annotation when your platform already has single-cell resolution or when you want one dominant label per cell/spot.

**How to choose**:
- use transfer methods when you have a strong matching reference
- use marker-based methods when marker genes are well established
- use automated methods for quick initial labeling, then validate biologically

See [Methods Reference](advanced/methods-reference.md) for supported annotation methods and exact requirements.

---

### Spatial Statistics

Spatial statistics answer different spatial questions.

**How to choose**:
- use global autocorrelation when you want one summary statistic for a gene
- use local hotspot methods when you want to locate spatially enriched regions
- use neighborhood or co-occurrence analyses when your question is about cell-type organization rather than gene-level spatial patterning

See [Methods Reference](advanced/methods-reference.md) for the full analysis-type matrix and required inputs.

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

Cell communication analysis depends on species-specific resources. Use the resource that matches the organism, especially for mouse data.

See [Methods Reference](advanced/methods-reference.md) for the canonical species and resource settings.

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

