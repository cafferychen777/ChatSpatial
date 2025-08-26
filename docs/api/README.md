# ChatSpatial API Reference

Complete reference for all ChatSpatial MCP tools, parameters, and data models.

## Overview

ChatSpatial provides **19 standardized MCP tools** for spatial transcriptomics analysis. Each tool follows the Model Context Protocol specification with:

- **JSON Schema validation** for all inputs and outputs
- **Structured error handling** with detailed error messages
- **Type-safe parameters** with automatic validation
- **Rich return types** including images, data, and metadata

## Tool Categories

| Category | Tools | Description |
|----------|-------|-------------|
| **[Data Management](#data-management)** | `load_data`, `preprocess_data` | Data loading, QC, and preprocessing |
| **[Cell Annotation](#cell-annotation)** | `annotate_cells` | 7 annotation methods with reference data support |
| **[Spatial Analysis](#spatial-analysis)** | `analyze_spatial_data`, `identify_spatial_domains`, `register_spatial_data`, `calculate_spatial_statistics` | Pattern analysis, domain identification, registration, and advanced statistics |
| **[Gene Analysis](#gene-analysis)** | `identify_spatial_genes`, `find_markers`, `analyze_enrichment` | Spatial variable genes, differential expression, and enrichment |
| **[Cell Communication](#cell-communication)** | `analyze_cell_communication` | Ligand-receptor interaction analysis |
| **[Deconvolution](#deconvolution)** | `deconvolve_data` | Cell type proportion estimation |
| **[Integration](#integration)** | `integrate_data` | Multi-modal and batch integration |
| **[Trajectory](#trajectory)** | `analyze_rna_velocity` | RNA velocity and trajectory inference |
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
- **`data_path`**: Path to the data file or directory
- **`data_type`**: Type of spatial data (auto, 10x_visium, slide_seq, merfish, seqfish, other, h5ad). If 'auto', will try to determine the type from the file extension or directory structure.
- **`name`**: Optional name for the dataset
- **`context`**: Optional MCP context for logging

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

Comprehensive preprocessing pipeline for spatial transcriptomics data.

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

**Key Features:**
- Quality control filtering
- Normalization and scaling
- Highly variable gene selection
- Dimensionality reduction (PCA, UMAP)
- Neighbor graph construction
- Leiden clustering

## Cell Annotation

### annotate_cells

Comprehensive cell type annotation with multiple methods.

**Signature:**
```python
annotate_cells(
    data_id: str,
    method: str = "marker_genes",
    reference_data_id: Optional[str] = None,
    marker_genes: Optional[Dict] = None,
    confidence_threshold: float = 0.5
) -> AnnotationResult
```

**Available Methods:**

| Method | Description | Requirements |
|--------|-------------|--------------|
| `marker_genes` | Traditional marker-based annotation | Marker gene dictionary |
| `tangram` | Spatial mapping with reference | Single-cell reference data |
| `sctype` | Automated cell type identification | Tissue type specification |
| `cell2location` | Probabilistic deconvolution | Reference signatures |
| `scanvi` | Semi-supervised annotation | Reference data with labels |
| `cellassign` | Probabilistic assignment | Marker gene matrix |
| `mllmcelltype` | Multi-modal LLM classifier | Pre-trained model |

**Example:**
```python
# Marker-based annotation
markers = {
    "T_cells": ["CD3D", "CD3E", "CD3G"],
    "B_cells": ["CD19", "CD20", "MS4A1"],
    "Macrophages": ["CD68", "CD163", "CSF1R"]
}

result = annotate_cells(
    data_id="dataset",
    method="marker_genes",
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

| Method | Description | Best For |
|--------|-------------|----------|
| `spagcn` | Graph convolutional networks | General spatial domains |
| `stagate` | Spatial-temporal attention | Complex tissue architecture |
| `banksy` | Spatial clustering | Neighborhood-aware domains |
| `leiden` | Community detection | Quick clustering |
| `louvain` | Modularity optimization | Alternative clustering |

### analyze_spatial_data

Comprehensive spatial pattern analysis.

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

## Gene Analysis

### identify_spatial_genes

Identify spatially variable genes using multiple methods.

**Signature:**
```python
identify_spatial_genes(
    data_id: str,
    method: str = "gaston",
    n_genes: int = 1000,
    alpha: float = 0.05
) -> SpatialVariableGenesResult
```

**Available Methods:**

| Method | Description | Strengths |
|--------|-------------|-----------|
| `gaston` | Poisson regression with spatial binning | Fast, robust |
| `spatialde` | Gaussian process models | Flexible patterns |
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

| Method | Description | Best For |
|--------|-------------|----------|
| `spatial_enrichmap` | Spatially-aware enrichment mapping | Spatial pathway analysis |
| `pathway_gsea` | Gene Set Enrichment Analysis | Ranked gene lists |
| `pathway_ora` | Over-representation analysis | Discrete gene sets |
| `pathway_enrichr` | Enrichr web service integration | Online databases |
| `pathway_ssgsea` | Single-sample GSEA | Sample-level enrichment |

**Key Features:**
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

## Visualization

### visualize_data

Create comprehensive visualizations with MCP image outputs.

**Signature:**
```python
visualize_data(
    data_id: str,
    plot_type: str = "spatial",
    feature: Optional[str] = None,
    color_by: Optional[str] = None,
    size: Tuple[int, int] = (10, 8),
    dpi: int = 300
) -> VisualizationResult
```

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

ChatSpatial implements comprehensive error handling:

1. **ValidationError**: Invalid parameters or data format
2. **DataError**: Missing data or incompatible datasets
3. **MethodError**: Algorithm-specific failures
4. **ResourceError**: Memory or computation limits
5. **SystemError**: File I/O or environment issues

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

## Advanced Usage

### Chaining Analysis

```python
# Complete workflow  
result = load_data(data_path="data.h5ad", name="sample")
preprocess_data(data_id=result.id)
identify_spatial_domains(data_id=result.id, method="spagcn")
annotate_cells(data_id=result.id, method="tangram", reference_data_id="ref")
analyze_cell_communication(data_id=result.id, method="liana")
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
```

## See Also

- **[Getting Started](../getting_started.md)**: Installation and setup
- **[Tutorials](../tutorials/README.md)**: Step-by-step guides
- **[GitHub Repository](https://github.com/cafferychen777/ChatSpatial)**: Source code and issues

