---
layout: default
title: API Reference
parent: Reference
nav_order: 1
has_children: true
description: "Complete API documentation for all ChatSpatial tools"
---

# ChatSpatial API Reference
{: .fs-7 }

Complete reference for all ChatSpatial MCP tools, parameters, and data models.

## Overview

ChatSpatial provides 16 MCP tools for spatial transcriptomics analysis. Each tool follows the Model Context Protocol specification with:

- JSON Schema validation for all inputs and outputs
- Structured error handling with detailed error messages
- Type-safe parameters with automatic validation
- Return types include images, data, and metadata

## Tool Categories

| Category | Tools | Description |
|----------|-------|-------------|
| **[Data Management](#data-management)** | `load_data`, `preprocess_data` | Data loading, QC, and preprocessing |
| **[Cell Annotation](#cell-annotation)** | `annotate_cells` | 7 annotation methods with reference data support |
| **[Spatial Analysis](#spatial-analysis)** | `analyze_spatial_data`, `identify_spatial_domains`, `register_spatial_data` | Comprehensive spatial pattern analysis, domain identification, and registration |
| **[Gene Analysis](#gene-analysis)** | `find_spatial_genes`, `find_markers`, `analyze_enrichment` | Spatial variable genes, differential expression, and enrichment |
| **[Cell Communication](#cell-communication)** | `analyze_cell_communication` | Ligand-receptor interaction analysis |
| **[Deconvolution](#deconvolution)** | `deconvolve_data` | Cell type proportion estimation |
| **[Integration](#integration)** | `integrate_samples` | Multi-modal and batch integration |
| **[Trajectory](#trajectory)** | `analyze_velocity_data`, `analyze_trajectory_data` | RNA velocity analysis and trajectory inference |
| **[Visualization](#visualization)** | `visualize_data` | 20 plot types with MCP image outputs |

## Quick Reference

### Essential Tools

```python
# Data loading and preprocessing
load_data(data_path="data.h5ad", name="dataset")
preprocess_data(data_id="dataset", normalize_total=True, log1p=True)

# Core analysis
identify_spatial_domains(data_id="dataset", method="spagcn")
annotate_cells(data_id="dataset", method="tangram")
analyze_cell_communication(data_id="dataset", method="liana")
analyze_enrichment(data_id="dataset", method="spatial_enrichmap")

# Advanced spatial analysis
register_spatial_data(source_id="section1", target_id="section2")
analyze_spatial_data(data_id="dataset", params={"analysis_type": "geary", "genes": ["gene"]})

# Visualization
visualize_data(data_id="dataset", plot_type="spatial_domains")
```

### Parameter Patterns

All tools follow consistent parameter patterns:

- **`data_path`**: Path to data file (load_data only)
- **`data_type`**: Data format specification (load_data only)  
- **`data_id`**: Required string identifier for loaded datasets
- **`method`**: Analysis method selection with fallbacks
- **`*_key`**: Keys for accessing data layers (e.g., `spatial_key`, `batch_key`)
- **`use_*`**: Boolean flags for optional features
- **`n_*`**: Numeric parameters (neighbors, components, etc.)
- **`*_threshold`**: Filtering and significance thresholds

## Data Management

### load_data

Load spatial transcriptomics data from various formats.

**Signature:**

```python
load_data(
    data_path: str,
    data_type: str = "auto", 
    name: Optional[str] = None,
    context: Context = None
) -> SpatialDataset
```

**Supported Formats:**
- **H5AD**: AnnData format with spatial coordinates
- **CSV**: Expression matrix with separate coordinate file
- **H5/HDF5**: Hierarchical data format
- **10x Visium**: Space Ranger outputs
- **Zarr**: Cloud-optimized arrays

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_path` | `str` | - | Path to the data file or directory |
| `data_type` | `str` | `"auto"` | Type of spatial data (auto, 10x_visium, slide_seq, merfish, seqfish, other, h5ad). If 'auto', will try to determine the type from the file extension or directory structure |
| `name` | `Optional[str]` | `None` | Optional name for the dataset |
| `context` | `Context` | `None` | Optional MCP context for logging |

**Example:**
```python
result = load_data(
    data_path="data/mouse_brain_visium.h5ad",
    data_type="h5ad",
    name="mouse_brain"
)
print(f"Loaded dataset: {result.id}")
```

### preprocess_data

Preprocessing pipeline for spatial transcriptomics data.

**Signature:**

```python
preprocess_data(
    data_id: str,
    min_genes: int = 200,
    min_cells: int = 3,
    normalize_total: bool = True,
    log1p: bool = True,
    highly_variable_genes: bool = True,
    n_top_genes: int = 2000,
    pca: bool = True,
    neighbors: bool = True,
    clustering: bool = True,
    umap: bool = True
) -> PreprocessingResult
```

**Features:**

- Quality control filtering
- Normalization and scaling
- Highly variable gene selection
- Dimensionality reduction (PCA, UMAP)
- Neighbor graph construction
- Leiden clustering

## Cell Annotation

### annotate_cells

Cell type annotation with multiple methods.

**Signature:**

```python
annotate_cells(
    data_id: str,
    method: str = "tangram",
    reference_data_id: Optional[str] = None,
    marker_genes: Optional[Dict] = None,
    confidence_threshold: float = 0.5
) -> AnnotationResult
```

**Available Methods:**

| Method | Description | Requirements |
|--------|-------------|--------------|
| `tangram` | Spatial mapping with reference | Single-cell reference data |
| `sctype` | Automated cell type identification | Tissue type specification |
| `cell2location` | Probabilistic deconvolution | Reference signatures |
| `scanvi` | Semi-supervised annotation | Reference data with labels |
| `cellassign` | Probabilistic assignment | Marker gene matrix |
| `mllmcelltype` | Multi-modal LLM classifier | Pre-trained model |

**Example:**
```python
# Reference-based annotation with Tangram
result = annotate_cells(
    data_id="spatial_dataset",
    method="tangram",
    reference_data_id="reference_scRNA_dataset"
)

# CellAssign with custom marker genes
markers = {
    "T_cells": ["CD3D", "CD3E", "CD3G"],
    "B_cells": ["CD19", "CD20", "MS4A1"],
    "Macrophages": ["CD68", "CD163", "CSF1R"]
}

result = annotate_cells(
    data_id="dataset",
    method="cellassign",
    marker_genes=markers
)
```

## Spatial Analysis

### identify_spatial_domains

Identify spatial domains and tissue architecture.

**Signature:**

```python
identify_spatial_domains(
    data_id: str,
    method: str = "spagcn",
    n_clusters: Optional[int] = None,
    resolution: float = 1.0,
    spatial_key: str = "spatial"
) -> SpatialDomainResult
```

**Available Methods:**

| Method | Description | Use Case |
|--------|-------------|----------|
| `spagcn` | Graph convolutional networks | General spatial domains |
| `stagate` | Spatial-temporal attention | Complex tissue architecture |
| `leiden` | Community detection | Quick clustering |
| `louvain` | Modularity optimization | Alternative clustering |

### analyze_spatial_data

Spatial pattern analysis.

**Signature:**

```python
analyze_spatial_data(
    data_id: str,
    analysis_type: str = "autocorrelation",
    genes: Optional[List[str]] = None,
    method: str = "moran"
) -> SpatialAnalysisResult
```

**Analysis Types:**

- **`autocorrelation`**: Spatial autocorrelation (Moran's I, Geary's C)
- **`hotspots`**: Hotspot detection (Getis-Ord Gi*)
- **`patterns`**: Spatial expression patterns
- **`neighborhoods`**: Neighborhood analysis

### register_spatial_data

Register and align spatial transcriptomics data across multiple tissue sections.

**Signature:**

```python
register_spatial_data(
    source_id: str,
    target_id: str,
    method: str = "paste",
    landmarks: Optional[List[Dict[str, Any]]] = None
) -> Dict[str, Any]
```

**Available Methods:**

| Method | Description | Use Case |
|--------|-------------|----------|
| `paste` | PASTE algorithm for spatial alignment | Multi-slice integration |

**Features:**

- Cross-section spatial alignment
- Transformation matrix computation
- Landmark-guided registration
- Batch correction integration
- Quality metrics for alignment assessment

**Example:**
```python
# Register consecutive tissue sections
result = register_spatial_data(
    source_id="section_1",
    target_id="section_2", 
    method="paste"
)

print(f"Registration successful with transformation matrix")
print(f"Alignment quality score: {result['alignment_score']:.3f}")
```

### analyze_spatial_data (Enhanced)

Unified spatial statistics analysis with support for 12 different analysis types.

**Signature:**

```python
analyze_spatial_data(
    data_id: str,
    params: Dict[str, Any]
) -> SpatialAnalysisResult
```

**Available Analysis Types:**

| Analysis Type | Description | Key Parameters |
|--------------|-------------|----------------|
| `moran` | Global Moran's I spatial autocorrelation | `genes`, `moran_n_perms` |
| `local_moran` | Local Moran's I (LISA) for hotspot detection | `genes`, `n_neighbors` |
| `geary` | Geary's C spatial autocorrelation | `genes`, `moran_n_perms` |
| `getis_ord` | Getis-Ord Gi* hot/cold spot analysis | `genes`, `n_neighbors` |
| `neighborhood` | Neighborhood enrichment analysis | `cluster_key`, `n_neighbors` |
| `co_occurrence` | Cell type co-occurrence patterns | `cluster_key`, `n_neighbors` |
| `ripley` | Ripley's K/L point pattern analysis | `cluster_key` |
| `centrality` | Graph centrality measures | `cluster_key` |
| `bivariate_moran` | Bivariate spatial correlation | `gene_pairs` |
| `join_count` | Join count for categorical data | `cluster_key` |
| `network_properties` | Spatial network analysis | `cluster_key` |
| `spatial_centrality` | Spatial-specific centrality | `cluster_key` |

**New Unified Gene Selection:**

The `genes` parameter now provides unified gene selection across all relevant analysis types:

```python
# Example: Local Moran's I analysis
result = analyze_spatial_data(
    data_id="tissue_dataset",
    params={
        "analysis_type": "local_moran",
        "genes": ["CD8A", "FOXP3"],  # Unified parameter
        "n_neighbors": 6
    }
)

# Example: Geary's C analysis  
result = analyze_spatial_data(
    data_id="tissue_dataset",
    params={
        "analysis_type": "geary", 
        "genes": ["GAPDH", "MKI67"],  # Same unified parameter
        "n_neighbors": 8
    }
)
```

## Gene Analysis

### find_spatial_genes

Identify spatially variable genes using multiple methods.

**Signature:**

```python
find_spatial_genes(
    data_id: str,
    method: str = "gaston",
    n_genes: int = 1000,
    alpha: float = 0.05
) -> SpatialVariableGenesResult
```

**Available Methods:**

| Method | Description | Strengths |
|--------|-------------|-----------|
| `gaston` | Poisson regression with spatial binning | Fast, reliable |
| `spatialde` | Gaussian process models | Variable patterns |
| `spark` | Generalized linear mixed models | Statistical rigor |
| `trendsceek` | Marked point processes | Spatial trends |

### find_markers

Find marker genes for cell types or spatial domains.

**Signature:**

```python
find_markers(
    data_id: str,
    groupby: str = "cell_type",
    method: str = "wilcoxon",
    n_genes: int = 100,
    logfc_threshold: float = 0.25
) -> DifferentialExpressionResult
```

### analyze_enrichment

Perform gene set enrichment analysis on spatial transcriptomics data.

**Signature:**

```python
analyze_enrichment(
    data_id: str,
    method: str = "spatial_enrichmap",
    gene_sets: Optional[Union[List[str], Dict[str, List[str]]]] = None,
    gene_set_database: str = "GO_Biological_Process",
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    smoothing: bool = True,
    min_genes: int = 10
) -> EnrichmentResult
```

**Available Methods:**

| Method | Description | Use Case |
|--------|-------------|----------|
| `spatial_enrichmap` | Spatially-aware enrichment mapping | Spatial pathway analysis |
| `pathway_gsea` | Gene Set Enrichment Analysis | Ranked gene lists |
| `pathway_ora` | Over-representation analysis | Discrete gene sets |
| `pathway_enrichr` | Enrichr web service integration | Online databases |
| `pathway_ssgsea` | Single-sample GSEA | Sample-level enrichment |

**Features:**

- Spatial awareness for tissue-specific pathways
- Multiple database support (GO, KEGG, Reactome)
- Custom gene set analysis
- Spatial smoothing and covariate correction
- Statistical significance testing with FDR correction

**Example:**
```python
# Custom gene set enrichment
custom_pathways = {
    "Neuronal_Signaling": ["SNAP25", "SYN1", "GRIN1", "GRIA1"],
    "Glial_Function": ["GFAP", "AQP4", "S100B", "ALDH1L1"]
}

result = analyze_enrichment(
    data_id="brain_dataset",
    method="spatial_enrichmap",
    gene_sets=custom_pathways,
    smoothing=True,
    n_neighbors=8
)

print(f"Found {result.n_significant} significant pathways")
```

## Cell Communication

### analyze_cell_communication

Analyze cell-cell communication using ligand-receptor interactions.

**Signature:**

```python
analyze_cell_communication(
    data_id: str,
    method: str = "liana",
    groupby: str = "cell_type",
    spatial_mode: str = "global",
    database: str = "consensus"
) -> CellCommunicationResult
```

**Available Methods:**

| Method | Description | Features |
|--------|-------------|----------|
| `liana` | Comprehensive LR analysis | Multiple databases, spatial modes |
| `cellphonedb` | Statistical interaction testing | Permutation testing |
| `cellchat_liana` | CellChat via LIANA | Pathway analysis |

**Spatial Modes:**

- **`global`**: Cell type-level interactions
- **`local`**: Spatially-aware interactions
- **`bivariate`**: Pairwise spatial analysis

## Deconvolution

### deconvolve_data

Estimate cell type proportions in spatial transcriptomics data.

**Signature:**

```python
deconvolve_data(
    data_id: str,
    method: str = "cell2location",
    reference_data_id: Optional[str] = None,
    n_factors: int = 50
) -> DeconvolutionResult
```

**Available Methods:**

| Method | Description | Requirements |
|--------|-------------|--------------|
| `cell2location` | Bayesian deconvolution | Reference single-cell data |
| `stereoscope` | Probabilistic deconvolution | Reference signatures |
| `rctd` | Robust cell type decomposition | Reference profiles |

*Full documentation will be added in future versions*

## Integration

### integrate_samples

Integrate multiple spatial transcriptomics datasets.

**Signature:**

```python
integrate_samples(
    data_ids: List[str],
    method: str = "harmony",
    batch_key: str = "batch"
) -> IntegrationResult
```

**Available Methods:**

| Method | Description | Use Case |
|--------|-------------|----------|
| `harmony` | Harmony batch correction | Simple batch effects |
| `scvi` | Variational integration | Complex batch effects |
| `combat` | ComBat batch correction | Gene expression normalization |

**Note:** Harmony parameters are hardcoded in the implementation for optimal performance:
- `sigma=0.1` (diversity clustering penalty parameter)
- `max_iter_harmony=10` (maximum iterations for convergence)
- `nclust=None` (automatic cluster number detection)
- `verbose=True` (progress display enabled)

*Full documentation will be added in future versions*

## Trajectory

### analyze_velocity_data

Analyze RNA velocity to understand cellular dynamics.

**Signature:**

```python
analyze_velocity_data(
    data_id: str,
    method: str = "scvelo",
    mode: str = "dynamical"
) -> VelocityResult
```

**Available Methods:**

| Method | Description | Features |
|--------|-------------|----------|
| `scvelo` | Standard RNA velocity analysis | Stochastic, deterministic, and dynamical models |
| `velovi` | Deep learning RNA velocity | More accurate velocity with uncertainty quantification (requires scvi-tools) |
| `sirv` | Reference-based velocity | Transfer velocity from reference dataset (not yet implemented) |

### analyze_trajectory_data

Infer cellular trajectories and pseudotime from spatial data.

**Signature:**

```python
analyze_trajectory_data(
    data_id: str,
    method: str = "cellrank",
    spatial_weight: float = 0.5
) -> TrajectoryResult
```

**Available Methods:**

| Method | Description | Features |
|--------|-------------|----------|
| `dpt` | Diffusion pseudotime | Classic pseudotime inference (no velocity needed) |
| `palantir` | Probabilistic trajectory inference | Branch probability analysis (no velocity needed) |
| `cellrank` | RNA velocity-based trajectory inference | Fate mapping and terminal states (requires velocity) |

**Important Note**: VELOVI is a velocity computation method (see `analyze_velocity_data` above), not a trajectory inference method. After computing velocity with VELOVI, use CellRank, Palantir, or DPT for trajectory inference.

*Full documentation will be added in future versions*

## Visualization

### visualize_data

Create visualizations with MCP image outputs.

**Signature:**

```python
visualize_data(
    data_id: str,
    params: VisualizationParameters
) -> Image
```

**Key Parameters in VisualizationParameters:**
- `plot_type`: str = "spatial" (visualization type)
- `feature`: Optional[Union[str, List[str]]] = None (gene/column to visualize)
- `colormap`: str = "viridis" (color scheme)  
- `figure_size`: Optional[Tuple[int, int]] = None (width, height)
- `dpi`: int = 100 (resolution)

**Plot Types (20 Total):**

| Type | Description | Use Case |
|------|-------------|----------|
| `spatial` | Spatial gene expression | Gene visualization |
| `spatial_domains` | Spatial domain overlay | Domain identification |
| `umap` | UMAP embedding | Dimensionality reduction |
| `heatmap` | Expression heatmap | Multi-gene comparison |
| `violin` | Distribution plots | Expression distributions |
| `deconvolution` | Cell type proportion maps | Deconvolution results |
| `cell_communication` | Communication networks | Interaction visualization |
| `multi_gene` | Multi-gene spatial panels | Gene comparison |
| `lr_pairs` | Ligand-receptor pairs | LR interaction analysis |
| `gene_correlation` | Gene correlation analysis | Co-expression patterns |
| `rna_velocity` | RNA velocity plots | Trajectory inference |
| `trajectory` | Developmental trajectories | Pseudotime analysis |
| `spatial_analysis` | Spatial statistics (6 subtypes) | Pattern analysis |
| `gaston_isodepth` | GASTON isodepth contours | Spatial gene patterns |
| `gaston_domains` | GASTON domain visualization | Spatial domains |
| `gaston_genes` | GASTON gene analysis | Spatially variable genes |
| `spatial_enrichment` | Spatial enrichment maps | Functional enrichment |
| `pathway_enrichment` | Pathway enrichment plots | GSEA visualization |
| `spatial_interaction` | Cell-cell interactions | Spatial communication |
| `integration_check` | Integration quality plots | Batch correction QC |

**MCP Integration:**

All visualizations return MCP Image objects for direct display in LLM agents like Claude Desktop.

## Error Handling

### Error Types

ChatSpatial implements error handling:

| Error Type | Description | Common Causes |
|------------|-------------|---------------|
| **ValidationError** | Invalid parameters or data format | Wrong parameter types, out-of-range values |
| **DataError** | Missing data or incompatible datasets | Missing required columns, incompatible data structures |
| **MethodError** | Algorithm-specific failures | Method not applicable to data type |
| **ResourceError** | Memory or computation limits | Insufficient memory, timeout exceeded |
| **SystemError** | File I/O or environment issues | File not found, permission denied |

### Error Response Format

```json
{
  "error": {
    "code": 1001,
    "message": "Invalid parameter: n_clusters must be positive",
    "type": "ValidationError",
    "details": {"parameter": "n_clusters", "value": -1},
    "suggestions": ["Use n_clusters > 0", "Set n_clusters=None for auto"]
  }
}
```

## Usage Examples

### Chaining Analysis

```python
# Complete workflow  
result = load_data(data_path="data.h5ad", name="sample")
preprocess_data(data_id=result.id)
identify_spatial_domains(data_id=result.id, method="spagcn")
annotate_cells(data_id=result.id, method="tangram", reference_data_id="ref")
analyze_cell_communication(data_id=result.id, method="liana")
analyze_enrichment(data_id=result.id, method="spatial_enrichmap")
visualize_data(data_id=result.id, plot_type="spatial_domains")
```

### Parameter Optimization

```python
# Test multiple resolutions
for res in [0.5, 1.0, 1.5, 2.0]:
    identify_spatial_domains(
        data_id="sample",
        resolution=res,
        method="spagcn"
    )
```

### Batch Processing

```python
# Process multiple samples
sample_files = ["sample1.h5ad", "sample2.h5ad", "sample3.h5ad"]
for sample_file in sample_files:
    result = load_data(data_path=f"data/{sample_file}", name=sample_file.replace(".h5ad", ""))
    preprocess_data(data_id=result.id)
    identify_spatial_domains(data_id=result.id)

# Register multiple tissue sections
sections = ["section_1", "section_2", "section_3"]
for i in range(len(sections)-1):
    register_spatial_data(
        source_id=sections[i+1], 
        target_id=sections[i],
        method="paste"
    )
```

## See Also

- **[Getting Started](../../getting-started/)**: Installation and setup
- **[Tutorials](../tutorials/README.html)**: Step-by-step guides
- **[GitHub Repository](https://github.com/cafferychen777/ChatSpatial)**: Source code and issues

