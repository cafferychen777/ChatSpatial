# Preprocessing Tool Documentation

## Overview

The `preprocessing.py` module provides a comprehensive preprocessing pipeline for spatial transcriptomics data. It implements adaptive preprocessing strategies that automatically adjust parameters based on data type (Visium, MERFISH, etc.), ensuring optimal processing for different experimental platforms.

## Purpose

Preprocessing is the critical first step in spatial transcriptomics analysis that prepares raw data for downstream analysis. This module handles:
- Quality control and filtering
- Normalization and transformation
- Feature selection
- Dimensionality reduction
- Clustering
- Spatial graph construction

## Architecture

### Adaptive Data Type Detection

The module automatically detects data type based on gene count:
- **MERFISH**: < 200 genes (targeted panels)
- **Visium**: > 10,000 genes (whole transcriptome)
- **Other**: 200-10,000 genes (custom platforms)

### Processing Pipeline

1. **QC Metrics Calculation** → 2. **Filtering** → 3. **Normalization** → 4. **HVG Selection** → 5. **Batch Correction** → 6. **Scaling** → 7. **PCA** → 8. **Neighborhood Graph** → 9. **UMAP** → 10. **Clustering** → 11. **Spatial Neighbors**

## Quality Control and Filtering

### QC Metrics Computed

```python
# Standard metrics
- total_counts: Total UMI/transcript counts per cell
- n_genes_by_counts: Number of genes detected per cell
- pct_counts_mt: Mitochondrial gene percentage (0 for spatial data)
```

### Filtering Parameters

#### Gene Filtering (`filter_genes_min_cells`)
Removes genes expressed in fewer than specified cells.

**Adaptive Defaults**:
- MERFISH: `max(1, n_cells/100)` - Keep most genes
- Other: `3` - Standard filtering

#### Cell Filtering (`filter_cells_min_genes`)
Removes cells with fewer than specified genes.

**Adaptive Defaults**:
- MERFISH: `min(50, n_genes/2)` - Relaxed for targeted panels
- Other: `200` - Standard quality threshold

## Normalization Methods

### 1. Log Normalization (`log`)
**Process**:
1. Normalize total counts to 10,000 per cell
2. Apply log1p transformation

**When to use**: Default choice for most data
```python
params = AnalysisParameters(normalization="log")
```

### 2. SCTransform (`sct`)
**Process**:
- Variance-stabilizing transformation
- Currently simplified implementation

**When to use**: High technical variability
```python
params = AnalysisParameters(normalization="sct")
```

### 3. No Normalization (`none`)
**When to use**: Pre-normalized data or custom workflows
```python
params = AnalysisParameters(normalization="none")
```

### 4. scVI Normalization (Future)
**Note**: Parameter exists but not yet implemented
```python
params = AnalysisParameters(use_scvi_preprocessing=True)
```

## Feature Selection

### Highly Variable Genes (HVGs)

**Parameters**:
- `n_hvgs`: Number of HVGs to select (default: 2000)

**Adaptive Behavior**:
- MERFISH: Uses all genes (bypass HVG selection)
- Visium: Selects top variable genes
- Automatic adjustment if n_hvgs > n_genes

## Dimensionality Reduction

### Principal Component Analysis (PCA)

**Parameters**:
- `n_pcs`: Number of components (default: 50)

**Adaptive Features**:
- Automatically reduces if n_pcs > min(n_cells, n_genes)
- Fallback to random PCA if standard fails
- Handles both dense and sparse matrices

### UMAP Embedding

**Purpose**: 2D visualization of high-dimensional data

**Fallbacks**:
1. Try standard UMAP
2. Fall back to t-SNE if UMAP fails
3. Use random coordinates as last resort

## Batch Effect Correction

### Combat Integration

**Automatic Trigger**: Applied when 'batch' column exists with multiple batches

**Process**:
1. Detect batch column
2. Apply Combat correction
3. Continue without correction if fails

**Example**:
```python
# Ensure data has 'batch' column
adata.obs['batch'] = batch_labels
# Correction applied automatically during preprocessing
```

## Clustering

### Leiden Algorithm

**Parameters**:
- Resolution: Automatically determined
- Based on neighborhood graph

**Fallbacks**:
1. Try Leiden clustering
2. Create dummy clusters if fails

## Spatial Components

### Spatial Neighbor Graph

**Computation Methods**:
1. Delaunay triangulation (preferred)
2. K-nearest neighbors (fallback)

**Adaptive Parameters**:
- Coordinate handling for different platforms
- Edge effect management

## Input Parameters

```python
class AnalysisParameters:
    # Normalization
    normalization: str = "log"  # Method: "log", "sct", "none"
    
    # Filtering
    filter_genes_min_cells: Optional[int] = None
    filter_cells_min_genes: Optional[int] = None
    
    # Subsampling
    subsample_spots: Optional[int] = None
    subsample_genes: Optional[int] = None
    subsample_random_seed: int = 42
    
    # Processing
    scale: bool = False  # Scale to unit variance
    n_hvgs: int = 2000  # Number of HVGs
    n_pcs: int = 50  # Number of PCs
    
    # Advanced (future)
    use_scvi_preprocessing: bool = False
```

## Usage Examples

### Example 1: Standard Visium Processing
```python
# Default preprocessing for 10x Visium
result = await preprocess_data(
    data_id="visium_sample",
    params=AnalysisParameters(
        normalization="log",
        filter_cells_min_genes=200,
        filter_genes_min_cells=3,
        n_hvgs=3000,
        n_pcs=50,
        scale=True
    )
)
```

### Example 2: MERFISH Preprocessing
```python
# Adapted for targeted panel
result = await preprocess_data(
    data_id="merfish_data",
    params=AnalysisParameters(
        normalization="log",
        # Use adaptive defaults
        filter_cells_min_genes=None,
        filter_genes_min_cells=None,
        scale=False  # Often not needed for MERFISH
    )
)
```

### Example 3: Large Dataset with Subsampling
```python
# Efficient processing for >100k cells
result = await preprocess_data(
    data_id="large_dataset",
    params=AnalysisParameters(
        subsample_spots=50000,
        subsample_genes=5000,
        normalization="log",
        n_hvgs=5000,
        n_pcs=100,
        scale=True
    )
)
```

### Example 4: Custom Filtering
```python
# Strict quality control
result = await preprocess_data(
    data_id="high_quality_data",
    params=AnalysisParameters(
        filter_cells_min_genes=500,
        filter_genes_min_cells=10,
        normalization="log",
        n_hvgs=4000
    )
)
```

### Example 5: Minimal Processing
```python
# For pre-processed data
result = await preprocess_data(
    data_id="processed_data",
    params=AnalysisParameters(
        normalization="none",
        filter_cells_min_genes=None,
        filter_genes_min_cells=None,
        n_hvgs=0  # Skip HVG selection
    )
)
```

### Example 6: Multi-Batch Dataset
```python
# Automatic batch correction
# First, ensure batch information exists
adata.obs['batch'] = ['batch1', 'batch2', ...]

result = await preprocess_data(
    data_id="multi_batch",
    params=AnalysisParameters(
        normalization="log",
        scale=True
    )
)
# Combat correction applied automatically
```

## Output Format

### PreprocessingResult
```python
class PreprocessingResult:
    data_id: str  # Dataset identifier
    n_cells: int  # Cells after filtering
    n_genes: int  # Genes after filtering
    n_hvgs: int  # Number of HVGs
    clusters: int  # Number of clusters
    qc_metrics: Dict[str, Any]  # Detailed QC info
```

### QC Metrics Dictionary
```python
{
    "n_cells_before": 5000,
    "n_cells_after": 4800,
    "n_genes_before": 20000,
    "n_genes_after": 15000,
    "mean_genes_per_cell": 3500,
    "mean_counts_per_cell": 12000,
    "median_genes_per_cell": 3200,
    "median_counts_per_cell": 11000
}
```

## Best Practices

### 1. Data Type Considerations
- **Visium**: Use standard parameters with scaling
- **MERFISH**: Relax filtering, skip scaling
- **Slide-seq**: Similar to Visium but check spot density
- **Custom**: Start with defaults, adjust based on QC

### 2. Parameter Selection
- Start with defaults for initial exploration
- Adjust filtering based on QC plots
- Increase n_hvgs for complex tissues
- Scale n_pcs with dataset size

### 3. Quality Control
- Always check QC metrics after preprocessing
- Visualize filtered vs unfiltered distributions
- Ensure biological variation is preserved
- Monitor batch effects

### 4. Performance Optimization
- Use subsampling for initial parameter tuning
- Process in batches for very large datasets
- Consider memory constraints for dense operations
- Cache preprocessed data for iterative analysis

## Integration with Downstream Tools

### Visualization
```python
# After preprocessing
vis_params = VisualizationParameters(
    plot_type="umap",
    feature="leiden",  # Show clusters
    colormap="tab20"
)
```

### Analysis
```python
# Ready for spatial analysis
spatial_params = SpatialAnalysisParameters(
    analysis_type="neighborhood",
    cluster_key="leiden"  # Use preprocessing clusters
)
```

## Error Handling and Robustness

### Fallback Mechanisms
1. **PCA**: Random PCA if standard fails
2. **UMAP**: t-SNE → random coordinates
3. **Clustering**: Dummy clusters if needed
4. **Spatial**: KNN if Delaunay fails

### Common Issues and Solutions

1. **"Matrix is singular"**
   - Reduce n_pcs
   - Check for zero-variance genes
   - Use scale=True

2. **"Memory error"**
   - Use subsampling
   - Process in chunks
   - Use sparse matrices

3. **"No HVGs found"**
   - Check normalization
   - Adjust n_hvgs
   - Verify gene filtering

4. **"UMAP failed"**
   - Check PCA results
   - Verify neighborhood graph
   - Try different n_neighbors

## Performance Considerations

### Computational Complexity
- Filtering: O(n)
- Normalization: O(n)
- HVG selection: O(n log n)
- PCA: O(min(n², p²))
- UMAP: O(n^1.14)

### Memory Usage
- Scales with n_cells × n_genes
- Sparse matrices reduce memory
- Subsampling for large datasets

### Optimization Tips
- Use HDF5 backing for huge datasets
- Enable chunked processing
- Leverage GPU for applicable steps
- Parallelize independent operations

## Future Enhancements

1. **scVI Integration**
   - Deep learning normalization
   - Advanced batch correction
   - Uncertainty quantification

2. **Additional Methods**
   - Seurat-style SCTransform
   - GLM-PCA for count data
   - Factor analysis

3. **Performance**
   - GPU acceleration
   - Distributed processing
   - Incremental PCA

4. **Automation**
   - Auto-parameter tuning
   - Quality thresholds
   - Optimal HVG selection