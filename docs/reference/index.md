---
layout: default
title: Reference
nav_order: 4
has_children: true
permalink: /reference/
---

# Reference Documentation
{: .no_toc }

Comprehensive reference materials for ChatSpatial users and developers.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

This section provides detailed reference information for all aspects of ChatSpatial, from API documentation to troubleshooting guides.

## Quick References

### 📋 [Quick Reference Guides](quick-reference/)
Essential information at your fingertips:
- [All Tools](quick-reference/all-tools.html) - Complete list of analysis tools
- [Common Workflows](quick-reference/common-workflows.html) - Standard analysis patterns
- [Troubleshooting Index](quick-reference/troubleshooting-index.html) - Quick problem resolution

### 🔧 [API Documentation](api/)
Technical documentation for developers:
- [Data Models](api/data_models.html) - Input/output data structures
- [Error Handling](api/error_handling.html) - Error codes and handling
- [MCP Protocol](api/mcp_protocol.html) - Model Context Protocol details

### 📊 [Data and Configuration](/)
Format specifications and setup guides:
- [Data Formats](data-formats.html) - Supported spatial data formats
- [Configuration](configuration.html) - Advanced configuration options
- [Performance](performance.html) - Optimization and benchmarking

### 🛠️ [Troubleshooting](troubleshooting/)
Problem-solving resources:
- [Common Issues](troubleshooting/common_issues.html) - Frequently encountered problems
- [Installation Problems](troubleshooting/installation.html) - Setup and dependency issues
- [Analysis Errors](troubleshooting/analysis_errors.html) - Method-specific troubleshooting

## Analysis Methods Reference

### Core Analysis Tools

| Tool | Description | Key Parameters |
|------|-------------|----------------|
| `load_data` | Load spatial transcriptomics data | `data_path`, `data_type` |
| `preprocess_data` | Quality control and normalization | `normalization_method`, `n_top_genes` |
| `visualize_data` | Create spatial visualizations | `plot_type`, `feature`, `colormap` |
| `annotate_cells` | Cell type annotation | `method`, `reference_data_id` |
| `analyze_spatial_data` | Spatial pattern analysis | `analysis_type`, `n_neighbors` |

### Advanced Analysis Tools

| Tool | Description | Requirements |
|------|-------------|--------------|
| `identify_spatial_domains` | Spatial domain identification | SpaGCN, STAGATE, Leiden/Louvain |
| `analyze_cell_communication` | Cell-cell interaction analysis | LIANA or CellPhoneDB |
| `analyze_trajectory_data` | RNA velocity and pseudotime | CellRank or Palantir |
| `deconvolve_data` | Spatial deconvolution | Cell2location or RCTD |
| `find_spatial_genes` | Spatially variable genes | GASTON, SpatialDE, or SPARK |

### Integration and Utilities

| Tool | Description | Use Cases |
|------|-------------|-----------|
| `integrate_samples` | Multi-sample integration | Batch effect correction |
| `register_spatial_data` | Spatial alignment | Multi-section studies |
| `analyze_enrichment` | Pathway enrichment | Functional interpretation |
| `find_markers` | Differential expression | Cluster characterization |

## Data Format Specifications

### Supported Input Formats

#### H5AD (AnnData) - Recommended
```python
# Expected structure
adata.X              # Gene expression matrix (cells x genes)
adata.obs            # Cell metadata
adata.var            # Gene metadata  
adata.obsm['spatial'] # Spatial coordinates (cells x 2)
adata.uns            # Unstructured metadata
```

#### 10X Visium Format
```
sample_folder/
├── filtered_feature_bc_matrix.h5
└── spatial/
    ├── tissue_positions_list.csv
    ├── scalefactors_json.json
    ├── tissue_hires_image.png
    └── tissue_lowres_image.png
```

#### H5 (HDF5) Format
```python
# Gene expression matrix
/matrix/             # Sparse matrix data
/matrix/barcodes     # Cell barcodes
/matrix/features     # Gene information
```

### Output Formats

#### Analysis Results
- **Tabular data**: CSV files with statistics and annotations
- **Visualizations**: PNG/PDF plots and interactive HTML
- **Processed data**: H5AD files with analysis results

#### Metadata Structure
```python
# Analysis metadata stored in adata.uns
adata.uns['spatial_domains']    # Domain assignments
adata.uns['cell_communication'] # Interaction results  
adata.uns['trajectory']         # Pseudotime results
adata.uns['deconvolution']      # Cell type proportions
```

## Method Parameters Reference

### Spatial Domain Identification

#### SpaGCN Parameters
```python
params = {
    "method": "spagcn",
    "n_domains": 7,           # Number of domains
    "spagcn_p": 0.5,         # Smoothing parameter
    "spagcn_s": 1,           # Scaling factor
    "spagcn_b": 49,          # Bandwidth
    "resolution": 0.5        # Clustering resolution
}
```

#### STAGATE Parameters
```python
params = {
    "method": "stagate", 
    "stagate_epochs": 1000,      # Training epochs
    "stagate_learning_rate": 0.001, # Learning rate
    "stagate_dim_output": 512,   # Output dimensions
    "n_domains": 7               # Target domains
}
```

### Cell Communication Analysis

#### LIANA Parameters
```python
params = {
    "method": "liana",
    "liana_resource": "consensus",    # Database
    "liana_local_metric": "cosine",   # Local metric
    "liana_global_metric": "morans",  # Global metric
    "perform_spatial_analysis": True
}
```

#### CellPhoneDB Parameters
```python
params = {
    "method": "cellphonedb",
    "cellphonedb_iterations": 1000,     # Permutations
    "cellphonedb_pvalue": 0.05,         # P-value threshold
    "cellphonedb_use_microenvironments": True
}
```

### Trajectory Analysis

#### CellRank Parameters
```python
params = {
    "method": "cellrank",
    "cellrank_n_states": 5,              # Terminal states
    "cellrank_kernel_weights": [0.8, 0.2], # Velocity/connectivity weights
    "spatial_weight": 0.5                 # Spatial constraint
}
```

## Error Codes Reference

### Common Error Types

| Error Code | Description | Solution |
|------------|-------------|----------|
| `DATA_LOAD_ERROR` | Failed to load data | Check file path and format |
| `MISSING_SPATIAL` | No spatial coordinates | Ensure spatial coordinates in obsm['spatial'] |
| `DEPENDENCY_ERROR` | Missing required package | Install missing dependencies |
| `MEMORY_ERROR` | Insufficient memory | Reduce dataset size or increase memory |
| `PARAMETER_ERROR` | Invalid parameters | Check parameter types and ranges |

### Method-Specific Errors

| Method | Error | Cause | Solution |
|--------|-------|-------|----------|
| SpaGCN | `HISTOLOGY_ERROR` | Missing histology image | Provide tissue image or disable histology |
| Cell2location | `GPU_ERROR` | CUDA issues | Use CPU or fix GPU setup |
| RCTD | `R_ERROR` | R/rpy2 problems | Check R installation |
| LIANA | `DATABASE_ERROR` | Missing interaction database | Install LIANA databases |

## Performance Guidelines

### Memory Requirements

| Dataset Size | Minimal | Advanced | Experimental |
|--------------|---------|----------|--------------|
| < 5K cells | 2GB RAM | 4GB RAM | 8GB RAM |
| 5K-20K cells | 4GB RAM | 8GB RAM | 16GB RAM |
| 20K-100K cells | 8GB RAM | 16GB RAM | 32GB RAM |
| > 100K cells | 16GB RAM | 32GB RAM | 64GB RAM |

### Processing Time Estimates

| Analysis Type | Small Dataset | Large Dataset |
|---------------|---------------|---------------|
| Basic preprocessing | < 1 min | 5-15 min |
| Spatial domains (SpaGCN) | 2-5 min | 15-30 min |
| Cell communication | 5-10 min | 30-60 min |
| Trajectory analysis | 10-20 min | 1-2 hours |
| Deep learning methods | 15-30 min | 2-4 hours |

### Optimization Tips

1. **Use appropriate data types** - Sparse matrices for large datasets
2. **Chunk processing** - Process data in batches for memory efficiency
3. **GPU acceleration** - Use CUDA for compatible methods
4. **Parallel processing** - Set n_jobs parameter for multi-core usage

## Version Compatibility

### Python Version Support

| Python Version | Minimal | Advanced | Experimental |
|----------------|---------|----------|--------------|
| 3.10 | ✅ | ✅ | ✅ |
| 3.11 | ✅ | ✅ | ✅ |
| 3.12 | ✅ | ✅ | ⚠️ |

### Key Dependencies

| Package | Minimal Version | Recommended | Notes |
|---------|----------------|-------------|-------|
| scanpy | 1.9.0+ | 1.10.0+ | Core single-cell analysis |
| squidpy | 1.2.0+ | 1.4.0+ | Spatial analysis tools |
| torch | 2.0.0+ | 2.1.0+ | Deep learning (advanced) |
| scvi-tools | 1.0.0+ | 1.1.0+ | Variational inference |

---

Need help finding something? Try the [search function](/) or browse the [troubleshooting guides](troubleshooting/).