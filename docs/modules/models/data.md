# Data Models Documentation

## Overview

The `data.py` module defines all parameter and data model classes used throughout the ChatSpatial MCP system. These Pydantic models provide type safety, validation, and clear interfaces for all analysis operations.

## Purpose

Data models in the MCP architecture serve as:
- **Type Safety**: Ensure correct parameter types and values
- **Validation**: Automatic validation of inputs
- **Documentation**: Self-documenting parameter structures
- **Consistency**: Uniform interfaces across all tools
- **Extensibility**: Easy to add new parameters
- **Serialization**: JSON serialization for MCP protocol

## Core Models

### SpatialDataset
```python
class SpatialDataset(BaseModel):
    """Represents a loaded spatial transcriptomics dataset"""
    id: str  # Unique identifier (e.g., "data_1")
    name: str  # Human-readable name
    data_type: str  # Type: "10x_visium", "merfish", etc.
    description: str  # Dataset description
```

**Purpose**: Returned when loading data, provides dataset metadata and unique ID for subsequent operations.

### AnalysisParameters
```python
class AnalysisParameters(BaseModel):
    """Parameters for data preprocessing and analysis"""
    # Normalization
    normalization: str = "log"  # "log", "sct", "none"
    
    # Filtering
    filter_genes_min_cells: Optional[int] = None
    filter_cells_min_genes: Optional[int] = None
    
    # Subsampling
    subsample_spots: Optional[int] = None
    subsample_genes: Optional[int] = None
    subsample_random_seed: int = 42
    
    # Processing
    scale: bool = False
    n_hvgs: int = 2000
    n_pcs: int = 50
    
    # Advanced
    use_scvi_preprocessing: bool = False
```

**Key Features**:
- Adaptive defaults based on data type
- Supports multiple normalization methods
- Optional subsampling for large datasets
- scVI integration placeholder

## Visualization Parameters

### VisualizationParameters
```python
class VisualizationParameters(BaseModel):
    """Parameters for data visualization"""
    # Core parameters
    plot_type: str = "spatial"
    feature: Optional[str] = None
    features: Optional[List[str]] = None
    
    # Spatial parameters
    spatial_key: str = "spatial"
    img_key: Optional[str] = None
    crop_coord: Optional[Tuple[int, int, int, int]] = None
    size: Optional[float] = None
    alpha: float = 1.0
    
    # Visual customization
    colormap: str = "viridis"
    figure_size: Tuple[int, int] = (8, 8)
    
    # Analysis-specific
    groups: Optional[List[str]] = None
    groupby: Optional[str] = None
    analysis_key: Optional[str] = None
    analysis_sub_type: Optional[str] = None
    
    # Advanced options
    show_legend: bool = True
    show_title: bool = True
    save_path: Optional[str] = None
```

**Plot Types**:
- `spatial`: Gene/feature expression in tissue
- `heatmap`: Expression heatmaps
- `violin`: Distribution plots
- `umap`: Dimensionality reduction
- `spatial_domains`: Domain visualization
- `cell_communication`: LR networks
- `deconvolution`: Cell type proportions
- Plus 10+ specialized types

## Cell Type Annotation Parameters

### AnnotationParameters
```python
class AnnotationParameters(BaseModel):
    """Parameters for cell type annotation"""
    # Method selection
    method: Literal["marker_genes", "correlation", "tangram", 
                    "scanvi", "cellassign", "gptcelltype", 
                    "scrgcl", "popv", "supervised", 
                    "mllmcelltype"] = "marker_genes"
    
    # Common parameters
    marker_genes: Optional[Dict[str, List[str]]] = None
    reference_data_id: Optional[str] = None
    batch_key: Optional[str] = None
    layer: Optional[str] = None
    confidence_threshold: float = 0.7
    
    # Method-specific parameters
    n_top_genes: int = 2000
    
    # mLLMCellType parameters
    mllm_params: Dict[str, Any] = {
        "model": "gpt-4o-mini",
        "temperature": 0.1,
        "max_cells_per_request": 500,
        "use_cell_ontology": True,
        "custom_prompt": None,
        "parallel_requests": 5
    }
```

**Method Selection**:
- Traditional: `marker_genes`, `correlation`
- Transfer learning: `tangram`, `scanvi`
- Deep learning: `cellassign`, `scvi`
- LLM-based: `mllmcelltype`, `gptcelltype`

## Spatial Analysis Parameters

### SpatialAnalysisParameters
```python
class SpatialAnalysisParameters(BaseModel):
    """Parameters for spatial pattern analysis"""
    # Analysis type
    analysis_type: Literal["neighborhood", "co_occurrence", 
                          "ripley", "morans_i", "centrality", 
                          "getis_ord", "spatial_variability",
                          "interaction_matrix", 
                          "spatial_connectivity"] = "neighborhood"
    
    # Common parameters
    cluster_key: Optional[str] = None
    n_neighbors: int = 30
    spatial_key: str = "spatial"
    
    # Method-specific parameters
    morans_i_gene: Optional[str] = None
    getis_ord_genes: Optional[List[str]] = None
    fdr_correction: bool = True
    centrality_type: str = "degree"
    ripley_radii: Optional[List[float]] = None
    
    # Advanced
    random_seed: int = 42
```

**Analysis Types**:
- `neighborhood`: Cell type enrichment
- `morans_i`: Spatial autocorrelation
- `getis_ord`: Hot/cold spots
- `ripley`: Point pattern analysis
- `centrality`: Network importance

## Trajectory Analysis Parameters

### RNAVelocityParameters
```python
class RNAVelocityParameters(BaseModel):
    """Parameters for RNA velocity analysis"""
    # Mode selection
    mode: Literal["stochastic", "deterministic", 
                  "dynamical", "sirv"] = "dynamical"
    
    # Common parameters
    min_shared_counts: int = 30
    n_pcs: int = 30
    n_neighbors: int = 30
    
    # Advanced parameters
    filter_genes: bool = True
    filter_genes_dispersion: Optional[Tuple[float, float]] = None
    
    # Dynamical mode specific
    fit_basal_transcription: bool = True
    steady_state_prior: Optional[List[bool]] = None
    
    # SIRV (spatial) specific
    spatial_weight: float = 0.5
```

### TrajectoryParameters
```python
class TrajectoryParameters(BaseModel):
    """Parameters for trajectory inference"""
    # Method selection
    method: Literal["cellrank", "palantir", "dpt"] = "cellrank"
    
    # Start/end points
    start_cell: Optional[str] = None
    end_cells: Optional[List[str]] = None
    
    # Common parameters
    n_neighbors: int = 30
    n_pcs: int = 30
    
    # CellRank specific
    weight_connectivity: float = 0.8
    weight_velocity: float = 0.2
    softmax_scale: float = 4
    n_states: Optional[int] = None
    
    # Palantir specific
    n_diffusion_components: int = 10
    knn_palantir: int = 30
```

## Integration Parameters

### IntegrationParameters
```python
class IntegrationParameters(BaseModel):
    """Parameters for multi-sample integration"""
    # Method selection
    method: Literal["harmony", "bbknn", "scanorama", 
                    "mnn", "scvi", "multivi", 
                    "totalvi"] = "harmony"
    
    # Common parameters
    batch_key: str = "sample"
    n_pcs: int = 50
    
    # Harmony specific
    theta: float = 1.0  # Diversity penalty
    lambda_harmony: float = 1.0  # Ridge penalty
    sigma: float = 0.1  # Cluster width
    
    # BBKNN specific
    neighbors_within_batch: int = 3
    n_neighbors: int = 30
    
    # Additional options
    scale_data: bool = True
    copy: bool = False
```

## Deconvolution Parameters

### DeconvolutionParameters
```python
class DeconvolutionParameters(BaseModel):
    """Parameters for spatial deconvolution"""
    # Method selection
    method: Literal["cell2location", "spotiphy", "rctd", 
                    "destvi", "stereoscope", 
                    "spotlight"] = "cell2location"
    
    # Reference data
    reference_data_id: Optional[str] = None
    cell_type_key: str = "cell_type"
    
    # Common parameters
    use_gpu: bool = True
    verbose: bool = True
    
    # Cell2location specific
    detection_alpha: float = 20.0
    max_epochs: int = 30000
    batch_size: Optional[int] = None
    
    # Method-specific params stored in dicts
    cell2location_params: Dict[str, Any] = {}
    spotiphy_params: Dict[str, Any] = {}
    # ... etc
```

## Spatial Domain Parameters

### SpatialDomainParameters
```python
class SpatialDomainParameters(BaseModel):
    """Parameters for spatial domain identification"""
    # Method selection
    method: Literal["spagcn", "stagate", "bayesspace", 
                    "banksy", "leiden", 
                    "louvain"] = "stagate"
    
    # Common parameters
    n_domains: Optional[int] = None
    alpha: float = 0.5  # Spatial weight
    n_neighbors: int = 10
    random_seed: int = 42
    
    # SpaGCN specific
    beta: int = 49
    refine: bool = True
    spatial_radius: Optional[float] = None
    
    # STAGATE specific
    hidden_dims: List[int] = [512, 30]
    n_epochs: int = 1000
    
    # Method-specific params in dicts
```

## Spatial Variable Genes Parameters

### SpatialVariableGenesParameters
```python
class SpatialVariableGenesParameters(BaseModel):
    """Parameters for spatial gene identification"""
    # Preprocessing
    preprocessing_method: Literal["glmpca", 
                                  "pearson_residuals"] = "glmpca"
    num_hidden_genes: int = 200
    
    # Neural network architecture
    spatial_hidden_layers: List[int] = [32, 16]
    expression_hidden_layers: List[int] = [32, 16]
    
    # Training parameters
    epochs: int = 1000
    learning_rate: float = 0.001
    batch_size: int = 256
    validation_split: float = 0.1
    early_stopping_patience: int = 50
    
    # Analysis parameters
    n_domains: int = 5
    num_bins: int = 70
    continuous_quantile: float = 0.9
    discontinuous_quantile: float = 0.9
    umi_threshold: int = 500
    
    # Technical
    device: str = "auto"
    random_seed: int = 42
    verbose: bool = True
```

## Cell Communication Parameters

### CellCommunicationParameters
```python
class CellCommunicationParameters(BaseModel):
    """Parameters for cell-cell communication analysis"""
    # Core parameters
    method: Literal["liana", "cellphonedb", "cellchat", 
                    "nichenet", "commot"] = "liana"
    groupby: str = "cell_type"
    use_raw: bool = True
    
    # Resource selection
    resource_name: str = "consensus"
    species: Optional[Literal["human", "mouse"]] = None
    
    # Filtering parameters
    min_cells: int = 10
    min_expr: float = 0.1
    
    # Analysis parameters
    n_perms: int = 1000
    seed: int = 42
    return_all_lrs: bool = True
    
    # Spatial parameters
    spatial_key: str = "spatial"
    connectivity_key: Optional[str] = None
    bandwidth: float = 150
    cutoff: float = 0.01
    
    # Performance
    use_gpu: bool = False
    verbose: bool = True
    inplace: bool = True
```

## Data Validation

### Type Validation
Pydantic automatically validates:
- Type correctness (int, float, str, etc.)
- Literal values (exact string matches)
- Optional fields (None allowed)
- List/Dict structures

### Range Constraints
Many parameters have implicit constraints:
```python
# Examples
n_neighbors: int = 30  # Should be > 0
alpha: float = 0.5  # Should be 0-1
n_pcs: int = 50  # Should be > 0
```

### Custom Validation
Tools may add additional validation:
```python
# In tool code
if params.n_neighbors > n_cells:
    params.n_neighbors = min(30, n_cells - 1)
```

## Default Values and Rationale

### Common Defaults
- `n_neighbors = 30`: Balance between local and global
- `n_pcs = 50`: Captures most variation
- `random_seed = 42`: Reproducibility
- `n_hvgs = 2000`: Balance information vs efficiency

### Method-Specific Defaults
- `theta = 1.0` (Harmony): Moderate integration
- `learning_rate = 0.001`: Stable training
- `min_cells = 10`: Reliable statistics
- `alpha = 0.5`: Balance spatial vs expression

## Relationships Between Models

### Pipeline Flow
```
SpatialDataset (load)
    ↓
AnalysisParameters (preprocess)
    ↓
AnnotationParameters (annotate)
    ↓
Multiple analysis types:
- SpatialAnalysisParameters
- TrajectoryParameters
- DeconvolutionParameters
- etc.
    ↓
VisualizationParameters (visualize)
```

### Shared Parameters
Many models share common parameters:
- `n_neighbors`: Graph construction
- `spatial_key`: Coordinate location
- `random_seed`: Reproducibility
- `verbose`: Output control

## Usage by Tools

### Data Loading
```python
# Server receives parameters
dataset = await load_data(
    data_path="/path/to/data.h5ad",
    data_type="10x_visium",
    name="Sample1"
)
# Returns SpatialDataset model
```

### Analysis Tools
```python
# Tool receives typed parameters
async def preprocess_data(
    data_id: str,
    params: AnalysisParameters,  # Type-safe
    context: Context
) -> PreprocessingResult:
    # Access validated parameters
    if params.normalization == "log":
        # ...
```

### Visualization
```python
# Flexible parameter handling
if isinstance(params, str):
    # Convert string to model
    params = VisualizationParameters(feature=params)
elif isinstance(params, dict):
    # Create from dict
    params = VisualizationParameters(**params)
```

## Best Practices

### 1. Start with Defaults
```python
# Use defaults initially
params = AnalysisParameters()  # All defaults

# Then customize as needed
params = AnalysisParameters(
    n_hvgs=3000,  # More genes
    scale=True    # Scale data
)
```

### 2. Parameter Validation
```python
# Let Pydantic validate
try:
    params = VisualizationParameters(
        plot_type="invalid_type"  # Will raise error
    )
except ValidationError as e:
    print(f"Invalid parameters: {e}")
```

### 3. Method-Specific Parameters
```python
# Use nested parameters for methods
params = DeconvolutionParameters(
    method="cell2location",
    cell2location_params={
        "detection_alpha": 30.0,
        "max_epochs": 50000
    }
)
```

### 4. Incremental Refinement
```python
# Start simple
params = SpatialAnalysisParameters(
    analysis_type="neighborhood"
)

# Add specifics as needed
params = SpatialAnalysisParameters(
    analysis_type="neighborhood",
    cluster_key="refined_types",
    n_neighbors=50
)
```

## Examples

### Example 1: Complete Analysis Pipeline
```python
# 1. Load data
dataset = await load_data(
    data_path="visium_data.h5ad",
    data_type="10x_visium"
)

# 2. Preprocess
preprocess_params = AnalysisParameters(
    normalization="log",
    n_hvgs=3000,
    scale=True,
    n_pcs=50
)
preprocess_result = await preprocess_data(
    dataset.id, preprocess_params
)

# 3. Annotate cells
annotation_params = AnnotationParameters(
    method="marker_genes",
    confidence_threshold=0.8
)
annotation_result = await annotate_cells(
    dataset.id, annotation_params
)

# 4. Spatial analysis
spatial_params = SpatialAnalysisParameters(
    analysis_type="neighborhood",
    cluster_key="cell_type"
)
spatial_result = await analyze_spatial_data(
    dataset.id, spatial_params
)

# 5. Visualize
vis_params = VisualizationParameters(
    plot_type="spatial",
    feature="cell_type",
    colormap="tab20"
)
image = await visualize_data(dataset.id, vis_params)
```

### Example 2: Advanced LLM Annotation
```python
params = AnnotationParameters(
    method="mllmcelltype",
    n_top_genes=3000,
    mllm_params={
        "model": "gpt-4o",
        "temperature": 0.05,  # More deterministic
        "max_cells_per_request": 1000,
        "use_cell_ontology": True,
        "custom_prompt": "Focus on immune cell subtypes",
        "parallel_requests": 10
    }
)
```

### Example 3: Multi-Sample Integration
```python
params = IntegrationParameters(
    method="harmony",
    batch_key=["patient", "batch"],  # Multiple levels
    theta=[2.0, 1.0],  # Different penalties
    n_pcs=100,  # More components
    scale_data=True
)
```

### Example 4: Complex Deconvolution
```python
params = DeconvolutionParameters(
    method="cell2location",
    reference_data_id="scRNA_reference",
    cell_type_key="detailed_types",
    use_gpu=True,
    cell2location_params={
        "detection_alpha": 200,
        "max_epochs": 50000,
        "batch_size": 2500,
        "lr": 0.002,
        "use_raw": False
    }
)
```

## Extending Models

### Adding New Parameters
```python
class NewToolParameters(BaseModel):
    """Parameters for new tool"""
    # Required parameters
    required_param: str
    
    # Optional with default
    optional_param: int = 10
    
    # Complex types
    gene_lists: Dict[str, List[str]]
    thresholds: List[float]
    
    # Validation
    class Config:
        extra = "forbid"  # No extra fields
```

### Version Compatibility
Models support backward compatibility:
```python
# Old version
params = {"method": "old_name"}

# New version with alias
method: Literal["new_name", "old_name"] = "new_name"
```

## Performance Considerations

### Model Creation
- Pydantic validation has minimal overhead
- Models are created once per tool call
- Validation catches errors early

### Serialization
- Models serialize to JSON automatically
- Used for MCP protocol communication
- Efficient for network transport

### Memory Usage
- Models are lightweight
- Parameter dictionaries minimize memory
- No data duplication

## Future Enhancements

1. **Validation Extensions**
   - Custom validators for ranges
   - Cross-field validation
   - Dynamic defaults

2. **Type Improvements**
   - More specific types
   - Union types for flexibility
   - Generic models

3. **Documentation**
   - Auto-generated from models
   - Interactive parameter builders
   - Validation helpers