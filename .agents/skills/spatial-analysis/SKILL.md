---
name: spatial-analysis
description: |
  Perform complete spatial transcriptomics analysis workflow to understand tissue architecture.
  Use when user wants to analyze spatial data from Visium, Xenium, MERFISH, Slide-seq, or other platforms.
  Triggers: "analyze spatial data", "what's the tissue structure", "spatial transcriptomics analysis",
  "load and analyze", "identify spatial domains", "cluster spatial data".
context: fork
agent: general-purpose
allowed-tools: Read, Grep, Glob
---

# Spatial Analysis

## Overview

This skill guides complete spatial transcriptomics analysis from raw data to biological insights.
The fundamental question: **What is the spatial organization of this tissue?**

## Workflow Decision Tree

```
START: User provides spatial data
    │
    ├─ Q: What platform?
    │   ├─ Visium/Visium HD → spot-based, may have histology
    │   ├─ Xenium/MERFISH/CosMx → single-cell resolution, imaging-based
    │   ├─ Slide-seq → spot-based, no histology
    │   └─ Unknown → check adata.uns for platform info
    │
    ├─ Q: Single sample or multiple?
    │   ├─ Single → proceed to standard workflow
    │   └─ Multiple → include integration step (Harmony recommended)
    │
    └─ Execute workflow based on answers
```

## Standard Workflow

### Step 1: Data Loading and Understanding

```
1. Load data using load_data tool
2. Examine dataset profile:
   - n_cells, n_genes
   - Available annotations (adata.obs columns)
   - Spatial coordinates availability
   - Tissue image availability
3. Understand the biological context from user
```

**Key questions to ask user:**
- What tissue/organ is this?
- What is your biological question?
- Do you have a reference dataset for cell type annotation?

### Step 2: Quality Control and Preprocessing

```
Platform-specific QC thresholds:

| Metric | Visium | Xenium/MERFISH | Slide-seq |
|--------|--------|----------------|-----------|
| min_genes | 200 | 50 | 100 |
| min_cells | 3 | 3 | 3 |
| max_mito% | 20% | 10% | 15% |
| HVG count | 2000-3000 | 500-1000 | 1500-2000 |
```

Use `preprocess_data` tool with appropriate parameters.

### Step 3: Spatial Domain Identification

**Method Selection Guide:**

| Scenario | Recommended Method | Reason |
|----------|-------------------|--------|
| Visium with histology | SpaGCN | Uses image features |
| Visium without histology | STAGATE/Leiden | Graph-based |
| Single-cell resolution | GraphST | Handles high resolution |
| Quick exploration | Leiden/Louvain | Fast, interpretable |
| Need fine structure | Higher resolution (1.5-2.0) | More clusters |
| Need broad domains | Lower resolution (0.3-0.5) | Fewer clusters |

Use `identify_spatial_domains` tool.

### Step 4: Validation and Refinement

**Checklist:**
- [ ] Do spatial domains match histological regions (if available)?
- [ ] Are marker genes for each domain biologically meaningful?
- [ ] Is the number of domains reasonable for the tissue?

Use `find_markers` to identify domain-specific genes.
Use `visualize_data` with plot_type="feature" to examine spatial patterns.

### Step 5: Multi-Sample Integration (if applicable)

When analyzing multiple samples:

```
1. Load all samples individually
2. Use integrate_samples with method selection:
   - harmony: Fast, good for batch correction (recommended default)
   - scvi: Deep learning, better for complex batches
   - bbknn: K-nearest neighbor based
   - scanorama: Panoramic stitching approach

3. Re-run clustering on integrated data
4. Validate integration quality
```

### Step 6: Visualization and Reporting

**Essential visualizations:**
1. Spatial scatter colored by clusters/domains
2. UMAP embedding for global structure
3. Marker gene expression patterns
4. Violin plots for key genes

Use `visualize_data` tool with appropriate plot_type.

## Platform-Specific Guidance

### Visium
- Resolution: ~55μm spots, ~1-10 cells per spot
- Typical workflow: Full analysis → deconvolution for cell composition
- Histology integration: Use SpaGCN or include image features

### Xenium/MERFISH
- Resolution: Single-cell
- Typical workflow: Full analysis → direct cell type annotation
- Large datasets: May need subsampling or chunked processing

### Slide-seq
- Resolution: ~10μm beads
- No histology: Rely purely on expression patterns
- Bead-specific QC important

## Output Expectations

A complete spatial analysis should provide:
1. **Spatial domains** with biological interpretation
2. **Marker genes** for each domain
3. **Visualizations** showing spatial organization
4. **Quality metrics** confirming data validity

## Downstream Analysis Pointers

After establishing spatial structure, users may want:
- Cell type composition → use `/cell-composition` skill
- Cell-cell communication → use `/cell-interaction` skill
- Developmental trajectories → use `/cell-dynamics` skill
- Publication figures → use `/publication-ready` skill
