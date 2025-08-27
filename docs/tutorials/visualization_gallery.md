# Visualization Gallery

Comprehensive guide to all visualization capabilities in ChatSpatial. This gallery showcases 20+ plot types with examples, parameters, and best practices for spatial transcriptomics data visualization.

## Overview

ChatSpatial provides rich visualization capabilities optimized for LLM agents and interactive analysis. All plots are returned as MCP Image objects for direct display in agent interfaces.

## Plot Categories

### 🗺️ Spatial Plots
- [Spatial Gene Expression](#spatial-gene-expression)
- [Spatial Domains](#spatial-domains) 
- [Cell Type Mapping](#cell-type-mapping)
- [Spatial Statistics](#spatial-statistics)

### 🧬 Expression Analysis
- [Gene Expression Heatmaps](#gene-expression-heatmaps)
- [Violin Plots](#violin-plots)
- [Dot Plots](#dot-plots)
- [Ridge Plots](#ridge-plots)

### 🔬 Dimensionality Reduction
- [UMAP Plots](#umap-plots)
- [t-SNE Plots](#t-sne-plots)
- [PCA Plots](#pca-plots)

### 🕸️ Network Analysis
- [Cell Communication Networks](#cell-communication-networks)
- [Gene Regulatory Networks](#gene-regulatory-networks)
- [Pathway Networks](#pathway-networks)

### 📊 Quality Control
- [QC Metrics](#qc-metrics)
- [Batch Effects](#batch-effects)
- [Technical Artifacts](#technical-artifacts)

---

## Spatial Plots

### Spatial Gene Expression

Visualize gene expression overlaid on spatial coordinates.

```python
# Single gene expression
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial",
    genes=["Gad1"],
    color_by="expression",
    save_path="spatial_gad1.pdf"
)

# Multiple genes in subplots
visualize_data(
    data_id="mouse_brain", 
    plot_type="spatial",
    genes=["Gad1", "Slc17a7", "Mbp", "Cx3cr1"],
    save_path="spatial_marker_genes.pdf"
)
```

**Best for**: 
- Identifying spatial expression patterns
- Comparing gene expression across regions
- Validating marker genes

**Parameters**:
- `genes`: List of gene names
- `color_by`: "expression" (default), "log_expression"
- `cmap`: Color palette ("viridis", "plasma", "magma")
- `spot_size`: Size of spatial spots (default: auto)

### Spatial Domains

Display identified spatial domains with distinct colors.

```python
# Basic spatial domains plot
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial_domains", 
    color_by="spatial_domains",
    save_path="spatial_domains.pdf"
)

# Domains with gene expression overlay
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial_domains",
    genes=["Gad1"],
    overlay="expression",
    save_path="domains_with_expression.pdf"
)
```

**Best for**:
- Visualizing tissue architecture
- Validating clustering results
- Comparing domain identification methods

**Parameters**:
- `domain_column`: Column containing domain labels (default: "spatial_domains")
- `overlay`: Gene expression overlay ("expression", "none")
- `legend`: Show legend (True/False)
- `alpha`: Transparency (0.0-1.0)

### Cell Type Mapping

Show cell type annotations on spatial coordinates.

```python
# Cell type spatial mapping
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial",
    color_by="cell_type",
    save_path="cell_type_spatial.pdf"
)

# Cell type proportions by region
visualize_data(
    data_id="mouse_brain", 
    plot_type="spatial_composition",
    color_by="cell_type",
    region_by="spatial_domains",
    save_path="cell_composition.pdf"
)
```

**Best for**:
- Validating cell type annotations
- Understanding cellular organization
- Comparing annotation methods

### Spatial Statistics

Visualize spatial statistics and patterns.

```python
# Spatial autocorrelation
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial_autocorr",
    genes=["Gad1", "Slc17a7"], 
    save_path="spatial_autocorr.pdf"
)

# Neighborhood analysis
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial_neighbors",
    method="delaunay",
    save_path="spatial_neighbors.pdf"
)
```

---

## Expression Analysis

### Gene Expression Heatmaps

Complex heatmaps for expression pattern analysis.

```python
# Marker gene heatmap
visualize_data(
    data_id="mouse_brain",
    plot_type="heatmap",
    genes=["Gad1", "Slc17a7", "Mbp", "Cx3cr1"],
    group_by="cell_type",
    save_path="marker_heatmap.pdf"
)

# Spatial domain markers
visualize_data(
    data_id="mouse_brain",
    plot_type="heatmap", 
    genes=top_markers,
    group_by="spatial_domains",
    cluster_genes=True,
    save_path="domain_markers.pdf"
)
```

**Best for**:
- Comparing expression across groups
- Identifying marker genes
- Clustering analysis validation

**Parameters**:
- `group_by`: Grouping variable
- `cluster_genes`: Cluster genes by expression
- `cluster_obs`: Cluster observations
- `cmap`: Color palette
- `standard_scale`: Standardization ("var", "obs", None)

### Violin Plots

Distribution plots for expression analysis.

```python
# Expression distribution by cell type
visualize_data(
    data_id="mouse_brain",
    plot_type="violin",
    genes=["Gad1"],
    group_by="cell_type",
    save_path="gad1_violin.pdf"
)

# Multiple genes violin plot
visualize_data(
    data_id="mouse_brain", 
    plot_type="violin",
    genes=["Gad1", "Slc17a7", "Mbp"],
    group_by="spatial_domains",
    save_path="multi_gene_violin.pdf"
)
```

**Best for**:
- Expression distribution comparison
- Statistical analysis visualization
- Outlier identification

### Dot Plots

Compact representation of expression patterns.

```python
# Marker gene dot plot
visualize_data(
    data_id="mouse_brain",
    plot_type="dotplot", 
    genes=marker_genes,
    group_by="cell_type",
    save_path="marker_dotplot.pdf"
)
```

**Best for**:
- Compact marker visualization
- Cross-group comparisons
- Publication-ready figures

---

## Dimensionality Reduction

### UMAP Plots

UMAP projections with various colorings.

```python
# Basic UMAP
visualize_data(
    data_id="mouse_brain",
    plot_type="umap",
    color_by="cell_type",
    save_path="umap_celltype.pdf"
)

# UMAP with gene expression
visualize_data(
    data_id="mouse_brain",
    plot_type="umap",
    genes=["Gad1"],
    color_by="expression",
    save_path="umap_gad1.pdf"
)

# Multiple coloring options
visualize_data(
    data_id="mouse_brain",
    plot_type="umap",
    color_by=["cell_type", "spatial_domains", "batch"],
    save_path="umap_multi_color.pdf"
)
```

**Best for**:
- Cell type visualization
- Batch effect assessment
- Clustering validation

**Parameters**:
- `color_by`: Coloring variable
- `legend_loc`: Legend position
- `size`: Point size
- `alpha`: Transparency

### t-SNE Plots

t-SNE projections for non-linear dimensionality reduction.

```python
# t-SNE visualization
visualize_data(
    data_id="mouse_brain",
    plot_type="tsne",
    color_by="cell_type", 
    perplexity=30,
    save_path="tsne_celltype.pdf"
)
```

### PCA Plots

Principal component analysis visualization.

```python
# PCA plot
visualize_data(
    data_id="mouse_brain",
    plot_type="pca",
    color_by="cell_type",
    components=[1, 2],
    save_path="pca_celltype.pdf"
)

# PCA loadings
visualize_data(
    data_id="mouse_brain",
    plot_type="pca_loadings",
    components=[1, 2], 
    n_genes=10,
    save_path="pca_loadings.pdf"
)
```

---

## Network Analysis

### Cell Communication Networks

Visualize cell-cell communication patterns.

```python
# Communication network
visualize_data(
    data_id="mouse_brain",
    plot_type="cell_communication",
    method="liana",
    save_path="comm_network.pdf"
)

# Pathway-specific networks
visualize_data(
    data_id="mouse_brain", 
    plot_type="communication_pathways",
    pathways=["NOTCH", "WNT", "TGFb"],
    save_path="pathway_networks.pdf"
)

# Spatial communication
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial_communication", 
    max_distance=100.0,
    save_path="spatial_comm.pdf"
)
```

**Best for**:
- Cell communication analysis
- Pathway visualization
- Spatial interaction patterns

### Gene Regulatory Networks

Network visualization for regulatory relationships.

```python
# TF-target networks
visualize_data(
    data_id="mouse_brain",
    plot_type="regulatory_network",
    transcription_factors=["Sox2", "Pax6"],
    save_path="tf_network.pdf"
)
```

---

## Quality Control

### QC Metrics

Essential quality control visualizations.

```python
# QC metrics overview
visualize_data(
    data_id="mouse_brain",
    plot_type="qc_metrics",
    metrics=["n_genes", "n_counts", "pct_mito"],
    save_path="qc_overview.pdf"
)

# Spatial QC
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial_qc", 
    save_path="spatial_qc.pdf"
)

# Gene detection rate
visualize_data(
    data_id="mouse_brain",
    plot_type="gene_detection",
    save_path="gene_detection.pdf"
)
```

**Best for**:
- Data quality assessment
- Outlier detection
- Preprocessing validation

### Batch Effects

Identify and visualize batch effects.

```python
# Batch effect visualization
visualize_data(
    data_id="mouse_brain",
    plot_type="batch_effects",
    batch_key="sample_id",
    save_path="batch_effects.pdf"
)

# Before/after correction
visualize_data(
    data_id="mouse_brain",
    plot_type="batch_correction",
    before_key="X_pca",
    after_key="X_pca_harmony", 
    save_path="batch_correction.pdf"
)
```

---

## Advanced Visualization Techniques

### Multi-Panel Figures

Create comprehensive multi-panel visualizations.

```python
# Multi-panel spatial analysis
visualize_data(
    data_id="mouse_brain",
    plot_type="multi_panel",
    panels=[
        {"type": "spatial", "color_by": "cell_type"},
        {"type": "spatial_domains"},
        {"type": "umap", "color_by": "cell_type"},
        {"type": "heatmap", "genes": marker_genes}
    ],
    layout="2x2",
    save_path="comprehensive_analysis.pdf"
)
```

### Interactive Plots

Generate interactive visualizations (when supported).

```python
# Interactive spatial plot
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial",
    interactive=True,
    save_path="interactive_spatial.html"
)
```

### Custom Styling

Apply custom styling and themes.

```python
# Custom styled plot
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial",
    style_config={
        "figure_size": (12, 8),
        "dpi": 300,
        "color_palette": "Set1",
        "background": "white",
        "grid": False
    },
    save_path="custom_styled.pdf"
)
```

## Plot Export Options

### File Formats

Supported export formats:
- **PDF**: Vector graphics, publication quality
- **PNG**: Raster format, web-friendly
- **SVG**: Scalable vector graphics
- **HTML**: Interactive plots
- **Base64**: MCP image objects for agents

### Resolution Settings

```python
# High resolution for publication
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial",
    dpi=300,
    figure_size=(10, 8),
    save_path="high_res_spatial.pdf"
)
```

## Best Practices

### 1. Plot Selection Guide

| Data Type | Recommended Plots |
|-----------|------------------|
| Expression patterns | Spatial, UMAP, Violin |
| Cell types | Spatial, UMAP, Dot plot |
| Spatial domains | Spatial domains, Heatmap |
| QC assessment | QC metrics, Spatial QC |
| Communication | Network, Spatial communication |

### 2. Color Palette Guidelines

- **Categorical data**: Use distinct colors (Set1, Set2)
- **Continuous data**: Use sequential palettes (viridis, plasma)
- **Diverging data**: Use diverging palettes (RdBu, RdYlBu)

### 3. Performance Tips

- Use appropriate resolution (150 DPI for screen, 300 DPI for print)
- Optimize point size for data density
- Consider subsampling for very large datasets
- Cache complex computations

### 4. Accessibility

- Use colorblind-friendly palettes
- Include clear legends and labels
- Provide alternative text descriptions
- Ensure sufficient contrast

This comprehensive visualization gallery provides all the tools needed for publication-quality spatial transcriptomics visualization. Choose the appropriate plot types based on your analysis goals and data characteristics.