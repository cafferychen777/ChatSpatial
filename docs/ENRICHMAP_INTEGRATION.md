# EnrichMap Integration Guide

## Overview

ChatSpatial integrates EnrichMap, a sophisticated tool for spatially-aware gene set enrichment analysis in spatial transcriptomics data. EnrichMap computes enrichment scores with spatial smoothing and covariate correction, providing more accurate results than traditional enrichment methods.

## Features

### Core Capabilities

1. **Spatially-aware Enrichment**
   - K-nearest neighbor spatial smoothing
   - Spatial covariate correction using GAM (Generalized Additive Models)
   - Batch-wise normalization for multi-sample analysis

2. **Flexible Gene Set Input**
   - Single gene list analysis
   - Multiple gene sets simultaneously
   - Custom gene weights support

3. **Statistical Analysis**
   - Spatial autocorrelation metrics (Moran's I, Getis-Ord)
   - Permutation testing for significance
   - Cluster-specific gene correlations

4. **Visualization**
   - Spatial enrichment maps
   - Multi-panel visualizations for multiple signatures
   - Integration with ChatSpatial's visualization framework

## Installation

EnrichMap is included in ChatSpatial's third_party directory and is automatically available. To install its dependencies:

```bash
pip install -e .[enrichmap]
```

Or install all dependencies:

```bash
pip install -e .[all]
```

### Required Dependencies

- statannotations
- pygam
- scikit-gstat
- adjustText
- splot
- dask==2024.11.2

## Usage Examples

### Basic Enrichment Analysis

#### Single Gene Set
```python
# Analyze T cell enrichment
"Analyze enrichment for genes CD3D, CD3E, CD8A and name it T_cell"

# With spatial smoothing disabled
"Analyze enrichment for CD3D, CD3E, CD8A without spatial smoothing"
```

#### Multiple Gene Sets
```python
# Analyze multiple immune signatures
"Analyze enrichment for these gene sets: T_cell: CD3D, CD3E, CD8A; B_cell: CD19, MS4A1, CD79A; Macrophage: CD68, CD163"
```

### Advanced Options

#### Batch Correction
```python
# For multi-sample analysis
"Analyze T cell enrichment with batch correction using sample_id column"
```

#### Custom Parameters
```python
# Adjust spatial neighbors
"Analyze enrichment for gene set CD3D, CD3E with 10 spatial neighbors"

# Without spatial correction
"Analyze enrichment without spatial covariate correction"
```

### Visualization

ChatSpatial follows a separation of analysis and visualization. After running enrichment analysis, use the visualization tool:

```python
# Basic spatial visualization
"Visualize enrichment plot for T_cell"
"Show spatial plot for T_cell_score"

# Multiple scores
"Visualize enrichment for T_cell, B_cell, and Macrophage"

# Violin plots by cluster
"Create violin plot with plot_type enrichment for T_cell_score"

# Gene contribution heatmap (if supported in future)
# Would show the weight of each gene in the enrichment score

# Advanced visualizations (if EnrichMap plotting is available)
# - Moran's I correlogram
# - Variogram analysis
# - Cross-correlation between signatures
```

Note: All visualizations are handled by the `visualize_data` tool with `plot_type="enrichment"`

### Spatial Metrics

Compute spatial statistics for enrichment scores:

```python
# Not yet exposed via MCP interface, but available in the API
# Would compute Moran's I, Getis-Ord statistics
```

## Technical Details

### Algorithm Overview

1. **Gene Weight Inference**
   - Automatically infers gene weights based on expression patterns
   - Or uses provided custom weights

2. **Score Computation**
   - Weighted sum of gene expressions
   - Z-score normalization (optional batch-wise)

3. **Spatial Smoothing**
   - K-nearest neighbor graph construction
   - Weighted averaging with spatial neighbors

4. **Spatial Correction**
   - GAM fitting with spatial coordinates as smooth terms
   - Residuals represent spatially-corrected scores

### Output Structure

The enrichment analysis returns:

```python
{
    "data_id": "data_1",
    "signatures": ["T_cell", "B_cell"],
    "score_columns": ["T_cell_score", "B_cell_score"],
    "gene_contributions": {
        "T_cell": {"CD3D": 0.4, "CD3E": 0.35, "CD8A": 0.25}
    },
    "summary_stats": {
        "T_cell": {
            "mean": 0.15,
            "std": 0.8,
            "min": -2.1,
            "max": 3.2,
            "median": 0.1,
            "q25": -0.4,
            "q75": 0.6,
            "n_genes": 3
        }
    }
}
```

### Data Storage

- Enrichment scores are stored in `adata.obs['{signature_name}_score']`
- Gene contributions are stored in `adata.uns['gene_contributions']`

## Best Practices

1. **Gene Set Size**
   - Minimum 2 genes per signature
   - Optimal: 10-50 genes
   - Very large sets may dilute signal

2. **Spatial Parameters**
   - Default 6 neighbors works well for most datasets
   - Increase for smoother results
   - Decrease for more local patterns

3. **Batch Effects**
   - Always use batch correction for multi-sample analysis
   - Ensure batch labels are in `adata.obs`

4. **Interpretation**
   - Positive scores indicate enrichment
   - Scores are relative within dataset
   - Consider spatial patterns, not just magnitude

## Troubleshooting

### Common Issues

1. **"No common genes found"**
   - Check gene names match between data and gene set
   - Ensure genes are in `adata.var_names`

2. **Memory errors with large datasets**
   - Reduce number of spatial neighbors
   - Process in batches if needed

3. **GAM fitting errors**
   - May occur with very few spots
   - Disable spatial correction: `correct_spatial_covariates=False`

## Citation

If you use EnrichMap in your research, please cite:

```
Celik C & Secrier M (2025). EnrichMap: Spatially-informed enrichment analysis 
for functional interpretation of spatial transcriptomics. 
bioRxiv: https://www.biorxiv.org/content/10.1101/2025.05.30.656960v1
```

## Advanced Usage (API)

For developers wanting to use EnrichMap directly:

```python
from chatspatial.tools.enrichment_analysis import (
    perform_enrichment_analysis,
    compute_spatial_metrics,
    cluster_gene_correlation
)

# Perform enrichment
result = await perform_enrichment_analysis(
    data_id="data_1",
    data_store=data_store,
    gene_sets={"T_cell": ["CD3D", "CD3E", "CD8A"]},
    n_neighbors=6,
    smoothing=True,
    correct_spatial_covariates=True,
    batch_key="sample_id"
)

# Compute spatial metrics
metrics = await compute_spatial_metrics(
    data_id="data_1",
    data_store=data_store,
    score_key="T_cell_score",
    metrics=["morans_i", "getis_ord"],
    n_perms=999
)

# Analyze cluster-specific correlations
correlations = await cluster_gene_correlation(
    data_id="data_1",
    data_store=data_store,
    signature_name="T_cell",
    cluster_key="leiden",
    correlation_method="pearson"
)
```