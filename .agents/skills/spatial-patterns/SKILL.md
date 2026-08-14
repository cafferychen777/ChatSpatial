---
name: spatial-patterns
description: |
  Identify genes and features with significant spatial patterns and organization.
  Use when user wants to find spatially variable genes, spatial autocorrelation, or spatial clustering patterns.
  Triggers: "spatially variable genes", "spatial pattern", "Moran's I", "spatial autocorrelation",
  "which genes vary spatially", "SVG", "SpatialDE", "SPARK", "hot spots", "Getis-Ord".
context: fork
agent: general-purpose
allowed-tools: Read, Grep, Glob
---

# Spatial Pattern Analysis

## Overview

This skill answers: **Which genes/features show significant spatial organization in this tissue?**

Two complementary analyses:
1. **Spatially Variable Genes (SVG)**: Find genes with non-random spatial expression
2. **Spatial Statistics**: Quantify and test spatial patterns

## Decision Tree: Which Analysis?

```
Q: What do you want to know?
│
├─ Which genes have spatial patterns?
│   └─ SVG Detection
│       ├─ SpatialDE - Gaussian Process, handles noise well
│       ├─ SPARK-X - Fast, scalable to large datasets
│       └─ Moran's I per gene - Classic, interpretable
│
├─ Are spatial patterns statistically significant?
│   └─ Spatial Statistics
│       ├─ Global Moran's I - Overall spatial autocorrelation
│       ├─ Getis-Ord Gi* - Local hot/cold spot detection
│       └─ Ripley's K/L - Point pattern analysis
│
└─ Where are the spatial hotspots?
    └─ Local Spatial Statistics
        ├─ LISA - Local Moran's I
        └─ Getis-Ord Gi* - Hot spot analysis
```

## Spatially Variable Gene Detection

### Method Selection

| Method | Speed | Best For | Output |
|--------|-------|----------|--------|
| SpatialDE | Moderate | Noise robustness | p-values, effect sizes |
| SPARK-X | Fast | Large datasets (>50k spots) | p-values, adjusted |
| Moran's I | Fast | Simple interpretation | I statistic, p-value |
| Sepal | Moderate | Multi-scale patterns | Scale-specific genes |

### Workflow

#### Step 1: Run SVG Detection

Use `find_spatial_genes` tool with:
- `method`: "spatialDE", "sparkx", or "moran"
- `n_top_genes`: Number of top genes to return (default 100)

#### Step 2: Interpret Results

**Key output columns:**
- `pval` / `padj`: Statistical significance
- `moranI` / `effect_size`: Strength of spatial pattern
- `gene`: Gene identifier

**Filtering criteria:**
```
Recommended thresholds:
- Adjusted p-value < 0.05 (stringent: < 0.01)
- Effect size > 0.1 (moderate spatial pattern)
- Expression level: Filter very low-expressed genes
```

#### Step 3: Visualize Top SVGs

Use `visualize_data` with `plot_type="spatial"`:
- Show expression of top SVGs spatially
- Overlay with tissue domains for context

## Spatial Statistics Analysis

### Global Statistics

**Moran's I**: Measures overall spatial autocorrelation
- Range: -1 (dispersed) to +1 (clustered)
- 0 = random distribution
- Use for: Overall pattern assessment

**Geary's C**: Alternative to Moran's I
- Range: 0 (clustered) to 2 (dispersed)
- 1 = random
- More sensitive to local variation

### Local Statistics

**LISA (Local Moran's I)**: Identifies local clusters
- High-High: Hot spots (high values surrounded by high)
- Low-Low: Cold spots (low values surrounded by low)
- High-Low / Low-High: Spatial outliers

**Getis-Ord Gi***: Hot spot analysis
- Positive z-score: Hot spot (high value cluster)
- Negative z-score: Cold spot (low value cluster)
- Use for: Precise hotspot identification

### Workflow

#### Step 1: Compute Spatial Statistics

Use `analyze_spatial_statistics` tool with:
- `statistic`: "moran", "geary", "getis_ord", or "ripley"
- `feature`: Gene or obs column to analyze

#### Step 2: Interpret Results

| Statistic | Interpretation |
|-----------|----------------|
| Moran's I > 0.3 | Strong positive spatial autocorrelation |
| Moran's I < -0.3 | Strong negative autocorrelation (checkerboard) |
| Getis-Ord z > 1.96 | Significant hot spot (p < 0.05) |
| Getis-Ord z < -1.96 | Significant cold spot (p < 0.05) |

## Common Analysis Scenarios

### Scenario 1: Find Tissue Boundary Genes
```
Goal: Genes marking transitions between regions
Approach:
1. Run SVG detection (SpatialDE recommended)
2. Overlay top SVGs with domain boundaries
3. Look for genes with gradient patterns
```

### Scenario 2: Identify Hot Spots of Activity
```
Goal: Where is gene X most active?
Approach:
1. Run Getis-Ord Gi* on gene expression
2. Map significant hot spots spatially
3. Cross-reference with cell type composition
```

### Scenario 3: Compare Spatial Organization
```
Goal: Do two conditions differ in spatial organization?
Approach:
1. Compute Moran's I for each condition separately
2. Compare I statistics
3. Test for significant difference
```

## Biological Interpretation

### What SVGs Tell You

- **Boundary markers**: Genes defining tissue compartments
- **Gradient genes**: Continuous spatial variation
- **Spot-specific**: Localized expression patterns
- **Tissue architecture**: Structural organization genes

### Combining with Other Analyses

| Combine With | Insight |
|--------------|---------|
| Spatial domains | Which genes drive domain identity |
| Cell composition | Cell type-specific spatial patterns |
| Cell communication | Spatially organized signaling |

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| No significant SVGs | Low spatial variation | Lower threshold, check data quality |
| Too many SVGs | Loose threshold | Use stricter adjusted p-value |
| Unexpected patterns | Batch effects | Check for technical artifacts |
| Slow computation | Large dataset | Use SPARK-X, subsample for exploration |

## Output Summary

Successful spatial pattern analysis provides:
- **SVG list**: Ranked spatially variable genes
- **Statistics**: Significance and effect sizes
- **Spatial maps**: Visualization of top patterns
- **Hotspots**: Localized areas of interest
