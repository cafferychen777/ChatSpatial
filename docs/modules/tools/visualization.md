# Visualization Tool Documentation

## Overview

The `visualization.py` module is the primary interface for creating publication-quality visualizations of spatial transcriptomics data. It provides a unified API for generating various plot types, from basic spatial expression plots to complex multi-panel figures for advanced analyses.

## Purpose

This module serves as the visualization engine for ChatSpatial, handling:
- Spatial expression mapping
- Dimensionality reduction plots
- Statistical visualizations
- Analysis-specific figures
- Multi-panel compositions
- Image optimization for MCP transport

## Architecture

### Core Components

1. **Main Entry Point**: `visualize_data()` function
2. **Plot Type Handlers**: Specialized functions for each visualization type
3. **Helper Functions**: Utilities for common operations
4. **Image Processing**: Conversion and optimization for MCP

### Key Design Principles

- **Modular**: Each plot type has its own handler
- **Flexible**: Extensive customization options
- **Robust**: Graceful error handling with fallbacks
- **Efficient**: Optimized for large datasets
- **Consistent**: Unified styling across plot types

## Available Plot Types

### Basic Visualizations

#### 1. Spatial (`spatial`)
- **Purpose**: Visualize gene expression or features in spatial context
- **Use cases**: Gene expression, cell types, continuous scores
- **Features**: Automatic color scaling, size adjustment, coordinate handling

#### 2. UMAP (`umap`)
- **Purpose**: Show dimensionality reduction embeddings
- **Use cases**: Cell clustering, batch effects, trajectories
- **Features**: Multiple coloring options, subset highlighting

#### 3. Heatmap (`heatmap`)
- **Purpose**: Display expression patterns across genes/cells
- **Use cases**: Marker genes, differential expression, gene modules
- **Features**: Hierarchical clustering, annotations, scaling options

#### 4. Violin (`violin`)
- **Purpose**: Show distribution of expression values
- **Use cases**: Gene expression by group, QC metrics
- **Features**: Statistical annotations, multi-gene support

### Analysis-Specific Visualizations

#### 5. Deconvolution (`deconvolution`)
- **Purpose**: Display cell type proportions from deconvolution
- **Use cases**: Cell2location, Stereoscope, DestVI results
- **Features**: Multi-panel layout, top cell types, spatial mapping

#### 6. Spatial Domains (`spatial_domains`)
- **Purpose**: Show identified spatial regions
- **Use cases**: STAGATE, SpaGCN, BayesSpace results
- **Features**: Discrete coloring, boundary visualization

#### 7. Cell Communication (`cell_communication`)
- **Purpose**: Visualize ligand-receptor interactions
- **Use cases**: LIANA, CellPhoneDB results
- **Features**: Network plots, spatial expression of LR pairs

#### 8. RNA Velocity (`velocity`)
- **Purpose**: Show RNA velocity vectors
- **Use cases**: scVelo, velocyto results
- **Features**: Stream plots, grid visualization

#### 9. Trajectory (`trajectory`)
- **Purpose**: Display pseudotime and lineages
- **Use cases**: CellRank, Palantir results
- **Features**: Pseudotime coloring, fate probabilities

#### 10. Spatial Analysis (`spatial_analysis`)
- **Purpose**: Visualize spatial statistics results
- **Sub-types**:
  - `morans_i`: Spatial autocorrelation
  - `neighborhood`: Neighborhood enrichment
  - `getis_ord`: Hot/cold spots
  - `ripley`: Spatial point patterns
  - `co_occurrence`: Cell type co-occurrence

### Advanced Visualizations

#### 11. Multi-gene (`multi_gene`)
- **Purpose**: Compare multiple genes spatially
- **Features**: Grid layout, consistent scaling

#### 12. LR Pairs (`lr_pairs`)
- **Purpose**: Ligand-receptor spatial patterns
- **Features**: Side-by-side comparison, interaction strength

#### 13. Gene Correlation (`gene_correlation`)
- **Purpose**: Spatial correlation between genes
- **Features**: Scatter plots, regression lines

#### 14. GASTON Visualizations
- **Sub-types**:
  - `gaston_isodepth`: Isodepth coordinates
  - `gaston_domains`: Spatial domains
  - `gaston_genes`: Gene expression patterns

#### 15. Enrichment (`enrichment`)
- **Purpose**: Gene set enrichment visualization
- **Features**: Spatial scores, smoothing effects

#### 16. GSEA (`gsea`)
- **Purpose**: Gene set enrichment analysis plots
- **Features**: Enrichment plots, leading edge

## Input Parameters

### VisualizationParameters
```python
class VisualizationParameters:
    # Core parameters
    plot_type: str = "spatial"  # Type of visualization
    feature: Optional[str] = None  # Gene or feature to plot
    features: Optional[List[str]] = None  # Multiple features
    
    # Spatial plot parameters
    spatial_key: str = "spatial"  # Coordinates key
    img_key: Optional[str] = None  # Background image
    crop_coord: Optional[Tuple] = None  # Crop coordinates
    
    # Visual customization
    colormap: str = "viridis"  # Color scheme
    figure_size: Tuple[int, int] = (8, 8)  # Figure dimensions
    point_size: Optional[float] = None  # Point/spot size
    alpha: float = 1.0  # Transparency
    
    # Analysis-specific
    analysis_key: Optional[str] = None  # Analysis result key
    cell_type_key: str = "cell_type"  # Cell type column
    domain_key: str = "spatial_domain"  # Domain column
    
    # Advanced options
    normalize: bool = True  # Normalize values
    log_scale: bool = False  # Log transform
    show_legend: bool = True  # Display legend
    save_path: Optional[str] = None  # Save location
```

## Implementation Details

### Plot Generation Workflow

1. **Parameter Validation**: Check required parameters for plot type
2. **Data Preparation**: Extract and transform relevant data
3. **Plot Creation**: Generate base figure with matplotlib/seaborn
4. **Customization**: Apply user-specified options
5. **Optimization**: Adjust for clarity and aesthetics
6. **Conversion**: Convert to MCP Image format

### Key Helper Functions

#### `_prepare_spatial_plot()`
- Handles coordinate systems
- Manages image backgrounds
- Adjusts for tissue boundaries

#### `_get_color_palette()`
- Selects appropriate colors
- Handles discrete vs continuous
- Ensures color-blind friendly options

#### `_optimize_figure_size()`
- Calculates optimal dimensions
- Considers spot density
- Maintains aspect ratios

### Error Handling

1. **Missing Data**: Clear error messages with suggestions
2. **Invalid Parameters**: Validation with helpful feedback
3. **Plotting Failures**: Fallback visualizations
4. **Memory Issues**: Automatic downsampling for large data

## Usage Examples

### Example 1: Basic Gene Expression
```python
# Single gene spatial expression
params = VisualizationParameters(
    plot_type="spatial",
    feature="CD3E",
    colormap="Reds",
    figure_size=(10, 10)
)
image = await visualize_data("data_1", params)
```

### Example 2: Cell Type Visualization
```python
# Show annotated cell types
params = VisualizationParameters(
    plot_type="spatial",
    feature="cell_type",
    colormap="tab20",
    point_size=50
)
image = await visualize_data("data_1", params)
```

### Example 3: Deconvolution Results
```python
# Cell type proportions from Cell2location
params = VisualizationParameters(
    plot_type="deconvolution",
    n_cell_types=6,  # Top 6 cell types
    colormap="Blues",
    figure_size=(15, 10)
)
image = await visualize_data("data_1", params)
```

### Example 4: Spatial Domains
```python
# STAGATE domains with boundaries
params = VisualizationParameters(
    plot_type="spatial_domains",
    domain_key="STAGATE_domains",
    show_boundaries=True,
    alpha=0.8
)
image = await visualize_data("data_1", params)
```

### Example 5: RNA Velocity
```python
# Velocity stream plot
params = VisualizationParameters(
    plot_type="velocity",
    velocity_key="velocity",
    basis="spatial",
    arrow_size=3
)
image = await visualize_data("data_1", params)
```

### Example 6: Multi-Gene Comparison
```python
# Compare T cell markers
params = VisualizationParameters(
    plot_type="multi_gene",
    features=["CD3E", "CD4", "CD8A", "FOXP3"],
    colormap="viridis",
    normalize=True,
    figure_size=(16, 12)
)
image = await visualize_data("data_1", params)
```

### Example 7: Trajectory Analysis
```python
# Pseudotime visualization
params = VisualizationParameters(
    plot_type="trajectory",
    pseudotime_key="palantir_pseudotime",
    fate_key="palantir_fate_probabilities",
    colormap="plasma"
)
image = await visualize_data("data_1", params)
```

### Example 8: Spatial Statistics
```python
# Moran's I results
params = VisualizationParameters(
    plot_type="spatial_analysis",
    analysis_sub_type="morans_i",
    feature="CCL21",
    show_statistics=True
)
image = await visualize_data("data_1", params)
```

## Best Practices

### 1. Plot Type Selection
- Use `spatial` for single features in tissue context
- Use `heatmap` for multiple genes across cells
- Use analysis-specific plots for method results
- Use `multi_gene` for comparative analysis

### 2. Color Selection
- Sequential colormaps for continuous data (viridis, Blues)
- Diverging colormaps for centered data (RdBu, coolwarm)
- Qualitative colormaps for categories (tab20, Set3)
- Consider color-blind friendly options

### 3. Figure Optimization
- Adjust `figure_size` based on spot density
- Use `point_size` to prevent overlapping
- Set appropriate `alpha` for dense regions
- Enable `normalize` for cross-sample comparison

### 4. Publication Quality
- Use high DPI settings (300+)
- Include scale bars and labels
- Maintain consistent styling
- Export in vector formats when possible

### 5. Performance Tips
- Downsample for exploratory analysis
- Use `crop_coord` for regions of interest
- Disable `show_legend` for cleaner export
- Pre-filter data before visualization

## Integration with Analyses

The visualization module automatically detects and visualizes results from:

1. **Preprocessing**: QC metrics, dimensionality reduction
2. **Annotation**: Cell types, confidence scores
3. **Spatial Analysis**: All spatial statistics methods
4. **Deconvolution**: All supported methods
5. **Domain Detection**: STAGATE, SpaGCN, etc.
6. **Trajectory**: CellRank, Palantir, VeloVI
7. **Communication**: LIANA, CellPhoneDB
8. **Enrichment**: EnrichMap, GSEA

## Error Messages and Troubleshooting

### Common Issues

1. **"Feature not found"**: Check feature name and available features
2. **"No spatial coordinates"**: Ensure spatial data is loaded correctly
3. **"Analysis not run"**: Run the analysis before visualization
4. **"Memory error"**: Reduce figure size or downsample data

### Debugging Tips

1. Check data structure with `data_store[data_id]['adata'].obs.columns`
2. Verify spatial coordinates in `.obsm['spatial']`
3. Confirm analysis results in `.uns` keys
4. Use smaller test datasets first

## Advanced Features

### 1. Custom Colormaps
```python
import matplotlib.pyplot as plt
custom_cmap = plt.cm.colors.LinearSegmentedColormap.from_list(
    'custom', ['white', 'blue', 'red']
)
params.colormap = custom_cmap
```

### 2. Multi-Panel Figures
The module automatically creates multi-panel layouts for:
- Deconvolution (top N cell types)
- Multi-gene comparisons
- Before/after comparisons

### 3. Interactive Elements
While primarily static, the module prepares data compatible with:
- Plotly conversion
- Bokeh integration
- Web-based viewers

## Future Enhancements

1. **3D Visualization**: For 3D spatial data
2. **Animation**: Time-series and dynamics
3. **Interactive Plots**: Zoom, pan, hover
4. **Custom Layouts**: User-defined arrangements
5. **Export Options**: SVG, PDF, high-res PNG
6. **Style Templates**: Journal-specific formatting