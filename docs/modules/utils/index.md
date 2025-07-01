# Utilities Module Documentation

## Overview

The `utils` directory contains essential utility modules that support the ChatSpatial MCP system. These utilities handle data loading, error management, image processing, parameter validation, and visualization, providing the foundational infrastructure for all analysis tools.

## Module Summary

| Module | Purpose | Key Functions |
|--------|---------|---------------|
| `data_loader.py` | Load spatial transcriptomics data | `load_spatial_data()`, format detection |
| `error_handling.py` | Enhanced error management | Custom exceptions, MCP formatting |
| `tool_error_handling.py` | Tool-specific error handling | `@mcp_tool_error_handler()` decorator |
| `image_utils.py` | Image processing for MCP | `fig_to_image()`, size optimization |
| `mcp_parameter_handler.py` | Parameter validation | Manual validation, user-friendly errors |
| `plotting.py` | Visualization utilities | Spatial plots, clusters, heatmaps |
| `pydantic_error_handler.py` | Pydantic error handling | `@mcp_pydantic_error_handler()` |
| `output_utils.py` | Output suppression | Clean MCP responses |

## Data Loader (`data_loader.py`)

### Purpose
Handles loading of various spatial transcriptomics data formats with automatic format detection and validation.

### Key Functions

#### `load_spatial_data()`
```python
async def load_spatial_data(
    data_path: str,
    data_type: str = "auto",
    name: Optional[str] = None
) -> Dict[str, Any]
```

**Supported Formats**:
- H5AD files (AnnData)
- 10x Visium
- Slide-seq
- MERFISH
- seqFISH
- Generic spatial data

**Auto-detection Logic**:
1. Check file extension (.h5ad)
2. Check directory structure (filtered_feature_bc_matrix)
3. Analyze gene count for platform inference

### Example Usage
```python
# Auto-detect format
dataset = await load_spatial_data("/path/to/data")

# Specify format
dataset = await load_spatial_data(
    "/path/to/visium",
    data_type="10x_visium"
)

# Load with custom name
dataset = await load_spatial_data(
    "data.h5ad",
    name="Patient_001"
)
```

## Error Handling (`error_handling.py`)

### Purpose
Provides enhanced error handling with custom exceptions and MCP-compliant error formatting.

### Custom Exceptions

```python
class ChatSpatialError(Exception):
    """Base exception for ChatSpatial"""
    
class DataError(ChatSpatialError):
    """Data-related errors"""
    
class AnalysisError(ChatSpatialError):
    """Analysis-related errors"""
    
class VisualizationError(ChatSpatialError):
    """Visualization-related errors"""
```

### Error Context Manager
```python
class ErrorHandler:
    @contextmanager
    def handle_errors(self, operation: str):
        """Context manager for error handling"""
        try:
            yield
        except Exception as e:
            formatted_error = self.format_error(e, operation)
            raise ChatSpatialError(formatted_error)
```

### Usage
```python
error_handler = ErrorHandler()

with error_handler.handle_errors("data loading"):
    # Risky operation
    data = load_complex_data()
```

## Tool Error Handling (`tool_error_handling.py`)

### Purpose
Provides MCP-specific error handling for tool functions, ensuring errors are returned as proper result objects.

### Main Decorator

#### `@mcp_tool_error_handler()`
```python
def mcp_tool_error_handler(
    result_class=None,
    error_prefix="Error in tool execution"
):
    """Decorator for MCP tool error handling"""
```

### Features
- Catches all exceptions in tool functions
- Returns errors as result objects
- Preserves error context
- Logs errors with context

### Example Usage
```python
@mcp_tool_error_handler()
async def analyze_data(data_id: str, params: Any) -> AnalysisResult:
    # Tool implementation
    if not data_exists(data_id):
        raise ValueError(f"Dataset {data_id} not found")
    
    # Analysis code...
    return AnalysisResult(...)
```

## Image Utilities (`image_utils.py`)

### Purpose
Handles image processing and optimization for MCP transport, ensuring images are properly sized and formatted.

### Key Functions

#### `fig_to_image()`
```python
def fig_to_image(
    fig,
    max_size_mb: float = 5.0,
    format: str = "PNG",
    dpi: int = 100
) -> Image:
    """Convert matplotlib figure to MCP Image"""
```

**Optimization Strategy**:
1. Try original quality
2. Reduce DPI if too large
3. Try JPEG with quality reduction
4. Resize if still too large

#### `create_placeholder_image()`
```python
def create_placeholder_image(
    message: str = "Visualization not available",
    width: int = 800,
    height: int = 600
) -> Image:
    """Create placeholder when visualization fails"""
```

### Size Optimization
```python
# Automatic size optimization
image = fig_to_image(plt.gcf(), max_size_mb=3.0)

# Manual DPI control
image = fig_to_image(fig, dpi=150)

# Force JPEG for smaller size
image = fig_to_image(fig, format="JPEG")
```

## Parameter Handler (`mcp_parameter_handler.py`)

### Purpose
Provides manual parameter validation with user-friendly error messages for MCP tools.

### Validation Functions

#### `manual_parameter_validation()`
```python
def manual_parameter_validation(*validators):
    """Decorator for parameter validation"""
    
# Usage
@manual_parameter_validation(
    ("params", validate_analysis_params),
    ("data_id", validate_data_id)
)
async def tool_function(data_id: str, params: Any):
    pass
```

#### Specific Validators
```python
def validate_analysis_params(params: Any) -> Any:
    """Validate and convert analysis parameters"""
    
def validate_visualization_params(params: Any) -> Any:
    """Validate visualization parameters"""
    
def validate_spatial_analysis_params(params: Any) -> Any:
    """Validate spatial analysis parameters"""
```

### Features
- Type conversion
- Default value application
- User-friendly error messages
- Backward compatibility

## Plotting Utilities (`plotting.py`)

### Purpose
Provides high-level plotting functions for common spatial transcriptomics visualizations.

### Key Functions

#### `plot_spatial()`
```python
def plot_spatial(
    adata: AnnData,
    color: str,
    spatial_key: str = "spatial",
    spot_size: float = None,
    cmap: str = "viridis",
    title: str = None,
    figsize: Tuple[int, int] = (8, 8),
    **kwargs
) -> plt.Figure:
    """Create spatial plot of gene expression or annotations"""
```

#### `plot_clusters()`
```python
def plot_clusters(
    adata: AnnData,
    cluster_key: str = "leiden",
    basis: str = "umap",
    **kwargs
) -> plt.Figure:
    """Plot clusters in reduced dimensions"""
```

#### `plot_heatmap()`
```python
def plot_heatmap(
    adata: AnnData,
    genes: List[str],
    groupby: str = None,
    **kwargs
) -> plt.Figure:
    """Create expression heatmap"""
```

### Advanced Features
- Automatic color palette selection
- Size optimization for spatial plots
- Flexible layout options
- Integration with scanpy plotting

## Pydantic Error Handler (`pydantic_error_handler.py`)

### Purpose
Handles Pydantic validation errors with clear, user-friendly messages.

### Main Decorator

#### `@mcp_pydantic_error_handler()`
```python
def mcp_pydantic_error_handler(
    custom_messages: Dict[str, str] = None
):
    """Decorator for Pydantic error handling"""
```

### Error Formatting
```python
# Converts Pydantic errors like:
# "field required"
# To user-friendly:
# "The 'data_id' parameter is required but was not provided"

# Custom messages
@mcp_pydantic_error_handler({
    "n_neighbors": "Number of neighbors must be a positive integer"
})
async def analyze_spatial(params: SpatialAnalysisParameters):
    pass
```

## Output Utilities (`output_utils.py`)

### Purpose
Manages output suppression and formatting for clean MCP responses.

### Key Functions

#### `suppress_output()`
```python
@contextmanager
def suppress_output():
    """Suppress stdout/stderr during execution"""
    with suppress_output():
        # Noisy operation
        verbose_function()
```

#### `clean_output()`
```python
def clean_output(text: str) -> str:
    """Clean text output for MCP"""
    # Removes ANSI codes
    # Strips excessive whitespace
    # Formats for readability
```

## Best Practices

### 1. Error Handling
- Always use appropriate error handlers
- Provide context in error messages
- Return errors as result objects in tools
- Log errors for debugging

### 2. Image Optimization
- Set appropriate size limits
- Use JPEG for photos, PNG for plots
- Include placeholder for failures
- Test with various figure sizes

### 3. Parameter Validation
- Use decorators for consistency
- Provide helpful error messages
- Support backward compatibility
- Validate early in execution

### 4. Data Loading
- Use auto-detection when possible
- Validate data format
- Handle missing files gracefully
- Support multiple formats

### 5. Visualization
- Use high-level functions
- Apply consistent styling
- Optimize for MCP transport
- Include error handling

## Integration Example

```python
from chatspatial.utils import (
    load_spatial_data,
    mcp_tool_error_handler,
    manual_parameter_validation,
    validate_analysis_params,
    fig_to_image,
    plot_spatial
)

@mcp_tool_error_handler()
@manual_parameter_validation(
    ("params", validate_analysis_params)
)
async def analyze_and_visualize(
    data_path: str,
    params: AnalysisParameters,
    context: Context = None
) -> AnalysisResult:
    """Complete analysis with visualization"""
    
    # Load data with auto-detection
    dataset = await load_spatial_data(data_path)
    
    if context:
        await context.info(f"Loaded {dataset['n_cells']} cells")
    
    # Perform analysis
    adata = dataset["adata"]
    # ... analysis code ...
    
    # Create visualization
    fig = plot_spatial(
        adata,
        color="total_counts",
        title="Spatial Distribution"
    )
    
    # Convert to MCP Image
    image = fig_to_image(fig, max_size_mb=3.0)
    
    return AnalysisResult(
        data_id=dataset["id"],
        visualization=image,
        statistics={...}
    )
```

## Performance Considerations

### Memory Management
- Utilities use generators where possible
- Image optimization reduces memory usage
- Error handlers clean up resources

### Speed Optimization
- Cached imports for heavy libraries
- Vectorized operations in plotting
- Efficient error formatting

### Scalability
- Handles large datasets gracefully
- Progressive image quality reduction
- Chunked processing support

## Future Enhancements

1. **Data Loading**
   - Additional format support
   - Cloud storage integration
   - Streaming for large files

2. **Error Handling**
   - Error recovery strategies
   - Detailed error tracking
   - User-friendly suggestions

3. **Visualization**
   - Interactive plots
   - 3D visualization
   - Animation support

4. **Performance**
   - Parallel processing
   - GPU acceleration
   - Caching mechanisms