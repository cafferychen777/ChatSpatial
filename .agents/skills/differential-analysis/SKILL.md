---
name: differential-analysis
description: |
  Find differentially expressed genes between conditions, regions, or cell types.
  Use when user wants to compare gene expression, find markers, or identify condition-specific changes.
  Triggers: "differential expression", "DEG", "marker genes", "compare conditions",
  "what genes differ", "find markers", "condition comparison", "fold change".
context: fork
agent: general-purpose
allowed-tools: Read, Grep, Glob
---

# Differential Analysis

## Overview

This skill answers: **What genes/features differ between groups in this spatial data?**

Three main comparison types:
1. **Marker genes**: What defines each cluster/cell type?
2. **Condition comparison**: How does treatment/disease change expression?
3. **Regional comparison**: What differs between spatial domains?

## Decision Tree: Which Comparison?

```
Q: What are you comparing?
│
├─ Clusters or cell types
│   └─ Marker Finding
│       ├─ Wilcoxon rank-sum - Default, robust
│       ├─ t-test - Parametric, fast
│       └─ Logistic regression - Handles confounders
│
├─ Experimental conditions (treated vs control)
│   └─ Condition Comparison
│       ├─ Pseudobulk + DESeq2 - Gold standard for replicates
│       ├─ MAST - Single-cell aware
│       └─ Mixed models - Complex designs
│
└─ Spatial domains or regions
    └─ Regional Comparison
        ├─ Same as marker finding
        └─ Consider spatial autocorrelation
```

## Marker Gene Finding

### Purpose

Identify genes that distinguish one group from all others (1 vs rest) or from a specific group (pairwise).

### Method Selection

| Method | When to Use | Strengths |
|--------|-------------|-----------|
| Wilcoxon | Default choice | Non-parametric, robust |
| t-test | Large datasets | Fast, well-understood |
| Logreg | Confounders present | Adjusts for covariates |

### Workflow

#### Step 1: Define Groups

Ensure grouping column exists in `adata.obs`:
- Clusters: `leiden`, `louvain`
- Cell types: `cell_type`, `annotation`
- Domains: `spatial_domain`

#### Step 2: Run Marker Finding

Use `find_markers` tool with:
- `groupby`: Column name defining groups
- `method`: "wilcoxon", "t-test", or "logreg"
- `n_genes`: Number of top markers to return per group

#### Step 3: Interpret Results

**Key output columns:**
- `names`: Gene name
- `scores`: Test statistic
- `pvals_adj`: Adjusted p-value (FDR corrected)
- `logfoldchanges`: Log2 fold change
- `pct_nz_group`: Percent expressing in group
- `pct_nz_reference`: Percent expressing in other groups

**Quality markers criteria:**
```
Strong markers typically have:
- logFC > 1 (2-fold change) or < -1
- padj < 0.05
- pct_nz_group > 0.25 (expressed in >25% of group)
- pct_nz_group >> pct_nz_reference (specificity)
```

## Condition Comparison

### Purpose

Compare expression between experimental conditions (e.g., disease vs healthy, treated vs control).

### Considerations

**Replicates matter:**
- With replicates: Use pseudobulk + DESeq2/edgeR
- Without replicates: Treat as exploratory only

**Confounders:**
- Batch effects: Include batch as covariate
- Cell type composition: Control or stratify

### Workflow

#### Step 1: Set Up Comparison

Identify:
- Condition column in `adata.obs`
- Reference group (baseline/control)
- Confounding variables to adjust

#### Step 2: Run Comparison

Use `compare_conditions` tool with:
- `condition_key`: Column with conditions
- `reference`: Reference/control group name
- `method`: "pseudobulk", "mast", or "wilcoxon"
- `covariates`: Optional confounders to adjust

#### Step 3: Interpret Results

**Volcano plot interpretation:**
- X-axis: Log2 fold change (effect size)
- Y-axis: -log10(p-value) (significance)
- Upper right: Upregulated in condition
- Upper left: Downregulated in condition

**Thresholds:**
```
Significant DEGs typically:
- |logFC| > 0.5 (at least 1.4-fold change)
- padj < 0.05
- For stricter: |logFC| > 1, padj < 0.01
```

## Regional/Spatial Domain Comparison

### Purpose

Compare expression between different spatial regions or domains.

### Special Considerations

**Spatial autocorrelation:**
- Nearby spots are not independent
- Standard tests may be anti-conservative
- Consider spatial-aware methods or pseudobulk

### Workflow

#### Step 1: Define Regions

Regions can come from:
- Spatial domain identification
- Manual annotation
- Anatomical boundaries

#### Step 2: Run Comparison

Use `find_markers` with regions as groups:
- Pairwise: Compare specific domains
- 1 vs rest: Find domain-specific genes

#### Step 3: Validate Spatially

Visualize top DEGs on spatial map:
- Confirm pattern matches expected region
- Check for artifacts at boundaries

## Visualization

### Essential Plots

1. **Volcano plot**: Overview of all comparisons
2. **Heatmap**: Top markers across groups
3. **Violin plots**: Distribution of key genes
4. **Spatial maps**: Expression in tissue context

Use `visualize_data` with:
- `plot_type="volcano"` for volcano plots
- `plot_type="heatmap"` for expression heatmaps
- `plot_type="violin"` for distributions

## Multiple Testing Correction

### Why It Matters

Testing thousands of genes inflates false positives:
- 20,000 genes × 0.05 = 1,000 false positives expected

### Correction Methods

| Method | Stringency | When to Use |
|--------|------------|-------------|
| Benjamini-Hochberg (FDR) | Moderate | Default, controls false discovery rate |
| Bonferroni | Strict | When false positives are costly |
| None | None | Exploratory only |

**Always use adjusted p-values for final results!**

## Common Analysis Scenarios

### Scenario 1: Find Cell Type Markers
```
Goal: What genes define each annotated cell type?
1. find_markers with groupby="cell_type"
2. Filter: padj < 0.05, logFC > 1, pct_nz > 0.25
3. Validate known markers present
4. Visualize top 3-5 per type
```

### Scenario 2: Disease vs Healthy
```
Goal: What changes in disease state?
1. compare_conditions with condition_key="disease_state"
2. Set reference="healthy"
3. Add covariates if batch effects present
4. Volcano plot to visualize
5. Pathway analysis on significant DEGs
```

### Scenario 3: Tumor Margin Analysis
```
Goal: What genes differ at tumor-stroma boundary?
1. Define regions: "tumor_core", "margin", "stroma"
2. find_markers with groupby="region"
3. Focus on margin-specific genes
4. Cross-reference with communication analysis
```

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| No significant DEGs | Low power or no effect | Check sample size, relax threshold |
| Too many DEGs | Loose threshold or batch | Tighten threshold, check batches |
| Known markers missing | Expression too low | Lower min expression filter |
| Inconsistent results | Technical noise | Use pseudobulk, add replicates |

## Output Summary

Successful differential analysis provides:
- **DEG list**: Ranked genes with statistics
- **Effect sizes**: Log fold changes
- **Significance**: Adjusted p-values
- **Visualizations**: Volcano, heatmap, spatial maps
- **Downstream input**: Gene lists for enrichment analysis
