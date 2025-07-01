# ChatSpatial API Reference

## Table of Contents

1. [Overview](#overview)
2. [MCP Protocol](#mcp-protocol)
3. [Core Server](#core-server)
4. [Analysis Tools](#analysis-tools)
5. [Data Models](#data-models)
6. [Utilities](#utilities)
7. [HTTP API](#http-api)
8. [Error Handling](#error-handling)
9. [Examples](#examples)

## Overview

ChatSpatial is a Model Context Protocol (MCP) server for spatial transcriptomics analysis. It provides a comprehensive suite of tools for analyzing spatial gene expression data through a standardized API.

### Architecture

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  MCP Client │────▶│ MCP Server  │────▶│   Analysis  │
│   (Claude)  │◀────│(ChatSpatial)│◀────│    Tools    │
└─────────────┘     └─────────────┘     └─────────────┘
                           │                     │
                           ▼                     ▼
                    ┌─────────────┐     ┌─────────────┐
                    │  Data Store │     │   Models    │
                    │  (AnnData)  │     │ (Pydantic) │
                    └─────────────┘     └─────────────┘
```

## MCP Protocol

### Transport Modes

**stdio (Default)**
```bash
python -m chatspatial.server
```

**HTTP/SSE**
```bash
python -m chatspatial.http_server --port 8000
```

### Tool Registration

All tools are automatically registered with the MCP server:

```python
@mcp.tool()
async def tool_name(params: ParameterModel) -> ResultModel:
    """Tool documentation"""
    pass
```

## Core Server

### Server Module (`server.py`)

Main MCP server implementation providing tool registration and execution.

#### Key Components

**Data Store**
```python
data_store: Dict[str, Any] = {}  # Global data storage
```

**Server Instance**
```python
mcp = FastMCP("ChatSpatial")
```

## Analysis Tools

### Data Loading

#### `load_data`
Load spatial transcriptomics data from various formats.

**Parameters:**
- `data_path` (str): Path to data file or directory
- `data_type` (str): Format type ["auto", "10x_visium", "slide_seq", "merfish", "seqfish", "h5ad"]
- `name` (str, optional): Dataset name

**Returns:** `SpatialDataset`

**Example:**
```python
result = await load_data(
    data_path="/path/to/data.h5ad",
    data_type="auto",
    name="Sample1"
)
```

### Preprocessing

#### `preprocess_data`
Preprocess spatial transcriptomics data with QC, normalization, and clustering.

**Parameters:**
- `data_id` (str): Dataset identifier
- `params` (AnalysisParameters): Preprocessing parameters

**Returns:** `PreprocessingResult`

**Example:**
```python
result = await preprocess_data(
    data_id="data_1",
    params=AnalysisParameters(
        normalization="log",
        n_hvgs=3000,
        scale=True
    )
)
```

### Visualization

#### `visualize_data`
Create various visualizations of spatial data.

**Parameters:**
- `data_id` (str): Dataset identifier
- `params` (VisualizationParameters): Visualization settings

**Returns:** Visualization resource URI

**Plot Types:**
- `spatial`: Spatial expression plots
- `umap`: Dimensionality reduction
- `heatmap`: Expression heatmaps
- `violin`: Distribution plots
- `spatial_domains`: Domain visualization
- `cell_communication`: LR networks
- `deconvolution`: Cell type proportions
- `trajectory`: Pseudotime/fate probabilities
- `spatial_analysis`: Analysis results
- `multi_gene`: Multiple gene comparison
- `lr_pairs`: Ligand-receptor pairs
- `gene_correlation`: Correlation plots
- `gaston_isodepth`: GASTON coordinates
- `gaston_domains`: GASTON domains
- `gaston_genes`: GASTON patterns

**Example:**
```python
result = await visualize_data(
    data_id="data_1",
    params=VisualizationParameters(
        plot_type="spatial",
        feature="CD3E",
        colormap="Reds"
    )
)
```

### Cell Type Annotation

#### `annotate_cells`
Annotate cell types using various methods.

**Parameters:**
- `data_id` (str): Dataset identifier
- `params` (AnnotationParameters): Annotation settings

**Methods:**
- `marker_genes`: Marker-based annotation
- `tangram`: Single-cell mapping
- `scanvi`: Semi-supervised learning
- `cellassign`: Probabilistic assignment
- `mllmcelltype`: LLM-based annotation

**Returns:** `AnnotationResult`

**Example:**
```python
result = await annotate_cells(
    data_id="data_1",
    params=AnnotationParameters(
        method="marker_genes",
        confidence_threshold=0.8
    )
)
```

### Spatial Analysis

#### `analyze_spatial_data`
Perform various spatial pattern analyses.

**Parameters:**
- `data_id` (str): Dataset identifier
- `params` (SpatialAnalysisParameters): Analysis settings

**Analysis Types:**
- `neighborhood`: Cell type enrichment
- `morans_i`: Spatial autocorrelation
- `getis_ord`: Hot/cold spots
- `ripley`: Point patterns
- `co_occurrence`: Cell co-localization
- `centrality`: Network importance

**Returns:** `SpatialAnalysisResult`

**Example:**
```python
result = await analyze_spatial_data(
    data_id="data_1",
    params=SpatialAnalysisParameters(
        analysis_type="neighborhood",
        n_neighbors=30
    )
)
```

### Differential Expression

#### `find_markers`
Find differentially expressed genes between groups.

**Parameters:**
- `data_id` (str): Dataset identifier
- `group_key` (str): Grouping column
- `group1` (str): First group
- `group2` (str): Second group or "rest"
- `n_top_genes` (int): Number of top genes
- `method` (str): Statistical method

**Returns:** `DifferentialExpressionResult`

**Example:**
```python
result = await find_markers(
    data_id="data_1",
    group_key="cell_type",
    group1="T_cells",
    group2="rest",
    n_top_genes=100,
    method="wilcoxon"
)
```

### RNA Velocity

#### `analyze_velocity_data`
Compute RNA velocity from spliced/unspliced counts.

**Parameters:**
- `data_id` (str): Dataset identifier
- `params` (RNAVelocityParameters): Velocity settings

**Modes:**
- `stochastic`: Original formulation
- `deterministic`: Without noise
- `dynamical`: Learn kinetics
- `sirv`: Spatial velocity

**Returns:** `RNAVelocityResult`

### Trajectory Analysis

#### `analyze_trajectory_data`
Infer cellular trajectories and pseudotime.

**Parameters:**
- `data_id` (str): Dataset identifier
- `params` (TrajectoryParameters): Trajectory settings

**Methods:**
- `cellrank`: Velocity-based
- `palantir`: Diffusion maps
- `dpt`: Diffusion pseudotime

**Returns:** `TrajectoryResult`

### Multi-Sample Integration

#### `integrate_samples`
Integrate multiple spatial samples.

**Parameters:**
- `data_ids` (List[str]): Dataset identifiers
- `params` (IntegrationParameters): Integration settings

**Methods:**
- `harmony`: Iterative correction
- `bbknn`: Graph-based
- `scanorama`: Panoramic stitching
- `mnn`: Mutual nearest neighbors

**Returns:** `IntegrationResult`

**Example:**
```python
result = await integrate_samples(
    data_ids=["data_1", "data_2", "data_3"],
    params=IntegrationParameters(
        method="harmony",
        theta=1.0
    )
)
```

### Spatial Deconvolution

#### `deconvolve_data`
Estimate cell type proportions in spots.

**Parameters:**
- `data_id` (str): Dataset identifier
- `params` (DeconvolutionParameters): Deconvolution settings

**Methods:**
- `cell2location`: Bayesian
- `spotiphy`: Fast PyTorch
- `rctd`: Robust decomposition
- `destvi`: Deep generative
- `stereoscope`: Probabilistic

**Returns:** `DeconvolutionResult`

### Spatial Domains

#### `identify_spatial_domains`
Identify spatial regions with similar expression.

**Parameters:**
- `data_id` (str): Dataset identifier
- `params` (SpatialDomainParameters): Domain settings

**Methods:**
- `stagate`: Graph attention
- `spagcn`: Graph convolution
- `banksy`: Multi-scale
- `leiden`: Standard clustering

**Returns:** `SpatialDomainResult`

### Cell Communication

#### `analyze_cell_communication`
Analyze ligand-receptor interactions.

**Parameters:**
- `data_id` (str): Dataset identifier
- `params` (CellCommunicationParameters): Communication settings

**Methods:**
- `liana`: LIANA+ framework
- Additional methods planned

**Returns:** `CellCommunicationResult`

### Enrichment Analysis

#### `analyze_enrichment`
Spatially-aware gene set enrichment.

**Parameters:**
- `data_id` (str): Dataset identifier
- `gene_sets` (Union[List[str], Dict[str, List[str]]]): Gene sets
- `spatial_key` (str): Spatial coordinates
- `n_neighbors` (int): Spatial neighbors
- `smoothing` (bool): Apply smoothing
- `correct_spatial_covariates` (bool): GAM correction

**Returns:** Enrichment scores and statistics

### Spatial Variable Genes

#### `find_spatial_genes`
Identify spatially variable genes using GASTON.

**Parameters:**
- `data_id` (str): Dataset identifier
- `params` (SpatialVariableGenesParameters): GASTON settings

**Returns:** `SpatialVariableGenesResult`

## Data Models

### Parameter Models

All parameter models inherit from Pydantic BaseModel for validation:

- `AnalysisParameters`: Preprocessing settings
- `VisualizationParameters`: Plot configuration
- `AnnotationParameters`: Cell type annotation
- `SpatialAnalysisParameters`: Spatial analysis
- `RNAVelocityParameters`: Velocity computation
- `TrajectoryParameters`: Trajectory inference
- `IntegrationParameters`: Sample integration
- `DeconvolutionParameters`: Spatial deconvolution
- `SpatialDomainParameters`: Domain identification
- `CellCommunicationParameters`: LR analysis
- `SpatialVariableGenesParameters`: GASTON settings

### Result Models

All result models provide structured outputs:

- `PreprocessingResult`: QC metrics, cell/gene counts
- `AnnotationResult`: Cell types, confidence scores
- `SpatialAnalysisResult`: Spatial statistics
- `DifferentialExpressionResult`: DE genes, fold changes
- `RNAVelocityResult`: Velocity computation status
- `TrajectoryResult`: Pseudotime, fate probabilities
- `IntegrationResult`: Integration metrics
- `DeconvolutionResult`: Cell type proportions
- `SpatialDomainResult`: Domain labels, quality
- `CellCommunicationResult`: LR interactions
- `SpatialVariableGenesResult`: Spatial genes

## Utilities

### Data Loading
```python
from chatspatial.utils import load_spatial_data

dataset = await load_spatial_data(
    data_path="/path/to/data",
    data_type="auto"
)
```

### Error Handling
```python
from chatspatial.utils import mcp_tool_error_handler

@mcp_tool_error_handler()
async def my_tool(params):
    # Tool implementation
    pass
```

### Image Processing
```python
from chatspatial.utils import fig_to_image

image = fig_to_image(
    matplotlib_figure,
    max_size_mb=3.0
)
```

### Parameter Validation
```python
from chatspatial.utils import manual_parameter_validation

@manual_parameter_validation(
    ("params", validate_params_func)
)
async def my_tool(params):
    pass
```

## HTTP API

### Endpoints

**Root Information**
```
GET /
```

**Health Check**
```
GET /health
```

**Create Session**
```
POST /sessions
```

**RPC Endpoint**
```
POST /rpc
Content-Type: application/json

{
    "jsonrpc": "2.0",
    "method": "tools/call",
    "params": {
        "name": "load_data",
        "arguments": {
            "data_path": "/path/to/data.h5ad"
        }
    },
    "id": "1"
}
```

**Server-Sent Events**
```
GET /sse
```

### Session Management

Sessions provide isolated data stores:

```javascript
// Create session
const response = await fetch('http://localhost:8000/sessions', {
    method: 'POST'
});
const { session_id } = await response.json();

// Use session
const headers = { 'X-Session-Id': session_id };
```

## Error Handling

### Error Types

**MCP Errors**
- Invalid parameters
- Tool not found
- Execution errors

**Analysis Errors**
- Data not found
- Invalid analysis parameters
- Computation failures

### Error Response Format

```json
{
    "jsonrpc": "2.0",
    "error": {
        "code": -32603,
        "message": "Error description",
        "data": {
            "details": "Additional information"
        }
    },
    "id": "request-id"
}
```

## Examples

### Complete Analysis Workflow

```python
# 1. Load data
dataset = await load_data(
    data_path="visium_sample.h5ad",
    data_type="10x_visium"
)

# 2. Preprocess
preprocess_result = await preprocess_data(
    data_id=dataset.id,
    params=AnalysisParameters(
        normalization="log",
        n_hvgs=3000,
        scale=True
    )
)

# 3. Annotate cells
annotation_result = await annotate_cells(
    data_id=dataset.id,
    params=AnnotationParameters(
        method="marker_genes"
    )
)

# 4. Spatial analysis
spatial_result = await analyze_spatial_data(
    data_id=dataset.id,
    params=SpatialAnalysisParameters(
        analysis_type="neighborhood",
        cluster_key="cell_type"
    )
)

# 5. Find markers
markers = await find_markers(
    data_id=dataset.id,
    group_key="cell_type",
    group1="T_cells",
    group2="rest"
)

# 6. Visualize
image = await visualize_data(
    data_id=dataset.id,
    params=VisualizationParameters(
        plot_type="spatial",
        feature="CD3E"
    )
)
```

### Multi-Sample Analysis

```python
# Load multiple samples
samples = []
for path in sample_paths:
    dataset = await load_data(path)
    samples.append(dataset.id)

# Integrate
integrated = await integrate_samples(
    data_ids=samples,
    params=IntegrationParameters(
        method="harmony"
    )
)

# Analyze integrated data
result = await analyze_spatial_data(
    data_id=integrated.integrated_data_id,
    params=SpatialAnalysisParameters(
        analysis_type="neighborhood"
    )
)
```

### Advanced Deconvolution

```python
# Load reference
ref = await load_data("scRNA_reference.h5ad")

# Deconvolve
deconv_result = await deconvolve_data(
    data_id="spatial_data",
    params=DeconvolutionParameters(
        method="cell2location",
        reference_data_id=ref.id,
        cell_type_key="cell_type",
        use_gpu=True
    )
)

# Visualize proportions
vis = await visualize_data(
    data_id="spatial_data",
    params=VisualizationParameters(
        plot_type="deconvolution",
        n_cell_types=6
    )
)
```

## Best Practices

1. **Error Handling**: Always handle tool errors appropriately
2. **Parameter Validation**: Use proper parameter models
3. **Session Management**: Create sessions for isolated analyses
4. **Resource Cleanup**: Close sessions when done
5. **Image Optimization**: Set appropriate size limits
6. **Batch Processing**: Use integration for multiple samples
7. **Method Selection**: Choose appropriate methods for data type

## Performance Tips

1. **Data Loading**: Use H5AD format for faster loading
2. **Preprocessing**: Subsample large datasets initially
3. **Visualization**: Limit figure size for MCP transport
4. **Integration**: Use approximation methods for large datasets
5. **GPU Usage**: Enable GPU for supported methods

## Troubleshooting

### Common Issues

**"Dataset not found"**
- Check data_id is correct
- Ensure data was loaded successfully

**"Memory error"**
- Reduce number of features/cells
- Use subsampling
- Enable sparse matrix operations

**"Tool timeout"**
- Increase timeout in HTTP transport
- Use smaller test datasets
- Check for infinite loops

**"Visualization too large"**
- Reduce DPI
- Use JPEG format
- Limit figure size

## Version Information

- ChatSpatial Version: 1.0.0
- MCP Protocol: 1.0
- Python: 3.8+
- AnnData: 0.8+