---
name: cell-composition
description: |
  Determine cell type composition at each spatial location through deconvolution or annotation.
  Use when user wants to know what cell types exist, their proportions, or where specific cells are located.
  Triggers: "cell type composition", "deconvolution", "what cells are here", "cell type proportions",
  "estimate cell types", "annotate cells", "which cell types".
context: fork
agent: general-purpose
allowed-tools: Read, Grep, Glob
---

# Cell Composition Analysis

## Overview

This skill answers the fundamental question: **What cell types exist at each spatial location, and in what proportions?**

For spot-based data (Visium), this requires deconvolution to estimate cell type mixtures.
For single-cell resolution data (Xenium/MERFISH), this requires direct cell type annotation.

## Decision Tree: Which Approach?

```
START: User wants cell type information
    │
    ├─ Q: What is the data resolution?
    │   │
    │   ├─ Spot-based (Visium, Slide-seq)
    │   │   └─ Q: Do you have a reference scRNA-seq dataset?
    │   │       ├─ YES → Deconvolution (see below)
    │   │       └─ NO → Q: Do you have marker genes?
    │   │           ├─ YES → Marker-based annotation
    │   │           └─ NO → Use public atlas as reference
    │   │
    │   └─ Single-cell (Xenium, MERFISH, CosMx)
    │       └─ Q: Do you have a reference dataset?
    │           ├─ YES → Transfer learning (Tangram/scANVI)
    │           └─ NO → Marker-based or LLM annotation
    │
    └─ Execute appropriate workflow
```

## Deconvolution Method Selection

### Quick Reference Table

| Your Scenario | Recommended Method | Why |
|---------------|-------------------|-----|
| Quick exploration | **FlashDeconv** | Fastest, good accuracy |
| Publication quality | **RCTD** (doublet mode) | Gold standard, well-validated |
| Large dataset (>50k spots) | **Cell2location** | Scalable, GPU-accelerated |
| Need spatial imputation | **CARD** | Can impute cell-type-specific expression |
| No matched reference | **Tangram** | More flexible with reference |
| Deep learning preference | **DestVI/Stereoscope** | Variational inference |

### Detailed Method Guide

#### FlashDeconv (Recommended Default)
- **Speed**: Fastest (~seconds to minutes)
- **Accuracy**: Good for most applications
- **When to use**: Initial exploration, iterative analysis
- **Limitations**: Less accurate for rare cell types

#### RCTD (Publication Standard)
- **Speed**: Moderate (~minutes)
- **Accuracy**: Excellent, especially with doublet mode
- **When to use**: Final results, manuscript figures
- **Modes**:
  - `doublet`: High-resolution platforms (Visium HD, Slide-seq)
  - `full`: Standard Visium
  - `multi`: When spots may contain >2 cell types
- **Note**: R-based, requires rpy2

#### Cell2location (Large Scale)
- **Speed**: Slow but scalable (~hours)
- **Accuracy**: Excellent for complex tissues
- **When to use**: Large datasets, need uncertainty estimates
- **Requirements**: GPU recommended, scvi-tools

#### CARD (Spatial Imputation)
- **Speed**: Moderate
- **Accuracy**: Good
- **Unique feature**: Can impute cell-type-specific gene expression
- **When to use**: Need spatial expression patterns per cell type
- **Note**: R-based

## Workflow: Deconvolution

### Step 1: Prepare Reference Data

**Requirements for reference scRNA-seq:**
- Cell type annotations in `adata.obs`
- Sufficient cells per type (>50 recommended)
- Matching species and tissue type
- Quality-controlled and normalized

```
If no reference available:
1. Check CellxGene for public atlases
2. Use Tabula Sapiens (human) or Tabula Muris (mouse)
3. Consider marker-based annotation instead
```

### Step 2: Validate Reference

Before deconvolution, verify:
- [ ] Cell type labels are present and clean
- [ ] No ambiguous categories (remove "unknown", "doublet")
- [ ] Sufficient cells per type
- [ ] Gene overlap with spatial data

### Step 3: Run Deconvolution

Use `deconvolve_data` tool with:
- `method`: Selected method from decision tree
- `reference_data_id`: Loaded reference dataset
- `cell_type_key`: Column name with cell type labels

### Step 4: Validate Results

**Quality checks:**
- [ ] Proportions sum to approximately 1
- [ ] Spatial patterns make biological sense
- [ ] Known cell type distributions match literature

**Red flags:**
- One cell type dominates everywhere → check reference quality
- Proportions don't vary spatially → possible normalization issue
- Expected cell types missing → check gene overlap

### Step 5: Visualize

Essential visualizations:
1. **Spatial pie charts**: Overall composition per spot
2. **Individual cell type maps**: Spatial distribution of each type
3. **Dominant cell type**: Which type is most abundant at each spot
4. **Diversity index**: Cellular heterogeneity across tissue

Use `visualize_data` with `plot_type="deconvolution"`.

## Workflow: Direct Annotation (Single-Cell Resolution)

For Xenium/MERFISH data:

### Option A: Transfer from Reference
```
1. Load reference scRNA-seq with cell type labels
2. Use annotate_cell_types with method="tangram" or "scanvi"
3. Validate transferred labels
```

### Option B: Marker-Based
```
1. Define marker genes per cell type
2. Use annotate_cell_types with method="cellassign" or "marker"
3. Validate against known markers
```

### Option C: LLM-Assisted (Experimental)
```
1. Identify cluster marker genes
2. Use mLLMCelltype for LLM-based annotation
3. Verify with domain knowledge
```

## Common Issues and Solutions

| Problem | Likely Cause | Solution |
|---------|-------------|----------|
| All spots same composition | Poor reference match | Try different reference or method |
| Missing expected cell types | Low gene overlap | Check shared genes, try Tangram |
| Unrealistic proportions | Normalization mismatch | Ensure consistent normalization |
| Method fails | Dependency issue | Check R packages installed (rpy2) |

## Output Interpretation

Deconvolution produces:
- **Proportions matrix**: Cell type fractions per spot (stored in `adata.obsm`)
- **Spatial patterns**: Where each cell type is enriched
- **Biological context**: Cell type organization in tissue architecture

Use these results for:
- Understanding tissue composition
- Identifying cell type niches
- Downstream cell-cell communication analysis
