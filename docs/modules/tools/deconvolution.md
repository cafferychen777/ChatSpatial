# Deconvolution Tool Documentation

## Overview

The `deconvolution.py` module provides spatial deconvolution capabilities to estimate cell type proportions in spatial transcriptomics spots. Since each spot typically contains multiple cells, deconvolution is essential for understanding cellular composition and spatial organization of tissues.

## Purpose

Spatial deconvolution addresses a fundamental challenge in spatial transcriptomics: each spot captures RNA from multiple cells. This module:
- Estimates cell type proportions per spot
- Leverages single-cell reference data
- Provides multiple state-of-the-art methods
- Enables cellular composition mapping

## Available Methods

### 1. Cell2location (`cell2location`)

**Description**: Bayesian hierarchical model that maps single-cell transcriptomes to spatial locations.

**Key Features**:
- Accounts for technical differences between technologies
- Models cell-type specific gene expression
- Provides uncertainty estimates
- Handles batch effects

**Requirements**:
- Single-cell reference with cell type annotations
- GPU recommended for large datasets

### 2. Spotiphy (`spotiphy`)

**Description**: Fast PyTorch-based deconvolution optimized for spatial data.

**Key Features**:
- Extremely fast (seconds to minutes)
- Graph-based regularization
- Lightweight implementation
- Good accuracy-speed tradeoff

**Requirements**:
- Reference data or marker genes
- Minimal computational resources

### 3. RCTD (`rctd`)

**Description**: Robust Cell Type Decomposition using statistical modeling.

**Key Features**:
- Doublet/multi-cell type detection
- Statistical confidence scores
- Platform effect correction
- Rigorous hypothesis testing

**Requirements**:
- R installation with spacexr package
- High-quality reference data

### 4. DestVI (`destvi`)

**Description**: Deep generative model from scvi-tools for joint analysis.

**Key Features**:
- Multi-resolution deconvolution
- Continuous cell state modeling
- Uncertainty quantification
- Batch correction built-in

**Requirements**:
- scvi-tools installation
- GPU strongly recommended

### 5. Stereoscope (`stereoscope`)

**Description**: Probabilistic model using negative binomial distributions.

**Key Features**:
- Handles sparse data well
- Cell state proportion inference
- Platform-agnostic
- Good for rare cell types

**Requirements**:
- scvi-tools installation
- Moderate computational resources

### 6. SPOTlight (`spotlight`)

**Description**: NMF-based deconvolution with marker gene selection.

**Key Features**:
- Seeded NMF approach
- Automatic marker selection
- Integration with Seurat
- Good interpretability

**Requirements**:
- R installation with SPOTlight package
- Seurat-compatible data

## Input Parameters

### DeconvolutionParameters
```python
class DeconvolutionParameters:
    # Method selection
    method: str = "cell2location"
    
    # Reference data
    reference_data_id: Optional[str] = None
    cell_type_key: str = "cell_type"
    
    # Common parameters
    use_gpu: bool = True
    verbose: bool = True
    
    # Method-specific parameters
    cell2location_params: Dict = {
        "detection_alpha": 20.0,  # Prior for detection rate
        "max_epochs": 30000,      # Training epochs
        "batch_size": None,       # Auto-determined
        "lr": 0.002              # Learning rate
    }
    
    spotiphy_params: Dict = {
        "n_top_markers": 50,      # Markers per cell type
        "alpha": 0.5,             # Spatial regularization
        "max_iter": 1000,         # Optimization iterations
        "lambda_spatial": 0.1     # Spatial smoothing
    }
    
    rctd_params: Dict = {
        "doublet_mode": "full",   # full/doublet/multi
        "max_cores": 8,           # Parallel cores
        "confidence_threshold": 5  # Min confidence
    }
    
    destvi_params: Dict = {
        "n_latent": 10,          # Latent dimensions
        "n_layers": 2,           # Neural network layers
        "dropout_rate": 0.1      # Dropout rate
    }
    
    stereoscope_params: Dict = {
        "st_batch_size": 100,    # Spatial batch size
        "sc_batch_size": 100,    # Single-cell batch size
        "max_epochs": 5000       # Training epochs
    }
    
    spotlight_params: Dict = {
        "cluster_resolution": 1,  # Clustering resolution
        "min_cells": 10,         # Min cells per type
        "hvg": 3000             # Variable genes
    }
```

## Output Format

### DeconvolutionResult
```python
class DeconvolutionResult:
    # Core results
    proportions: pd.DataFrame  # Cell type proportions
    cell_types: List[str]      # Detected cell types
    
    # Method metadata
    method_used: str
    parameters_used: Dict[str, Any]
    
    # Quality metrics
    statistics: Dict[str, Any] = {
        "mean_entropy": float,  # Diversity measure
        "dominant_cell_type_pct": float,
        "n_cell_types": int,
        "sparsity": float
    }
    
    # Storage keys
    proportions_key: str  # Key in adata.obsm
    
    # Visualization hints
    visualization_params: Dict = {
        "plot_type": "deconvolution",
        "n_cell_types": 6
    }
    
    # Optional method-specific
    confidence_scores: Optional[pd.DataFrame]
    convergence_info: Optional[Dict]
```

## Implementation Details

### Common Preprocessing

All methods share common preprocessing steps:

1. **Reference Data Validation**
   - Check cell type annotations
   - Verify gene overlap
   - Handle gene name conversions

2. **Gene Selection**
   - Intersection of spatial and reference genes
   - Optional marker gene prioritization
   - Minimum expression filtering

3. **Normalization Alignment**
   - Match normalization between datasets
   - Scale factor computation
   - Platform effect estimation

### Helper Functions

#### `prepare_reference_data()`
- Filters low-quality cells
- Computes cell type signatures
- Selects informative genes

#### `validate_cell_types()`
- Ensures consistent naming
- Removes rare cell types
- Handles hierarchical annotations

#### `compute_proportions_summary()`
- Calculates diversity metrics
- Identifies dominant cell types
- Computes spatial statistics

## Usage Examples

### Example 1: Basic Cell2location Deconvolution
```python
# Load reference single-cell data
ref_result = await load_data(
    "reference_scRNA.h5ad",
    data_type="h5ad",
    name="immune_reference"
)

# Run deconvolution
deconv_result = await deconvolve_data(
    data_id="spatial_data",
    params=DeconvolutionParameters(
        method="cell2location",
        reference_data_id=ref_result.id,
        cell_type_key="cell_type",
        use_gpu=True
    )
)

# Results in deconv_result.proportions DataFrame
```

### Example 2: Fast Deconvolution with Spotiphy
```python
# Quick analysis for exploration
result = await deconvolve_data(
    data_id="spatial_data",
    params=DeconvolutionParameters(
        method="spotiphy",
        reference_data_id="ref_data",
        spotiphy_params={
            "n_top_markers": 30,  # Fewer markers for speed
            "max_iter": 500,      # Fewer iterations
            "alpha": 0.3          # Less regularization
        }
    )
)
```

### Example 3: High-Confidence with RCTD
```python
# Rigorous analysis with doublet detection
result = await deconvolve_data(
    data_id="spatial_data",
    params=DeconvolutionParameters(
        method="rctd",
        reference_data_id="ref_data",
        rctd_params={
            "doublet_mode": "doublet",  # Detect doublets
            "confidence_threshold": 10,   # High confidence
            "max_cores": 16              # Use parallelization
        }
    )
)
```

### Example 4: Comparing Methods
```python
methods = ["cell2location", "spotiphy", "stereoscope"]
results = {}

for method in methods:
    results[method] = await deconvolve_data(
        data_id="spatial_data",
        params=DeconvolutionParameters(
            method=method,
            reference_data_id="ref_data",
            use_gpu=True
        )
    )
    
# Compare results
comparison_df = pd.DataFrame({
    method: result.proportions["T_cells"] 
    for method, result in results.items()
})
```

### Example 5: Region-Specific Analysis
```python
# Focus on tumor region
tumor_mask = adata.obs["region"] == "tumor"
tumor_spots = adata[tumor_mask].obs_names

result = await deconvolve_data(
    data_id="spatial_data",
    params=DeconvolutionParameters(
        method="destvi",
        reference_data_id="ref_data",
        spot_subset=tumor_spots,  # Analyze subset
        destvi_params={
            "n_latent": 15,      # More complex model
            "n_layers": 3
        }
    )
)
```

### Example 6: Custom Marker Genes
```python
# Use known markers instead of full reference
marker_genes = {
    "T_cells": ["CD3D", "CD3E", "CD8A"],
    "B_cells": ["CD19", "MS4A1", "CD79A"],
    "Myeloid": ["CD14", "CD68", "FCGR3A"]
}

result = await deconvolve_data(
    data_id="spatial_data",
    params=DeconvolutionParameters(
        method="spotiphy",
        marker_genes=marker_genes,  # Skip reference
        spotiphy_params={
            "n_top_markers": None  # Use all provided
        }
    )
)
```

## Best Practices

### 1. Reference Data Preparation
- **Quality Control**: Remove low-quality cells
- **Annotation**: Ensure accurate cell type labels
- **Balance**: Avoid extreme imbalances in cell type proportions
- **Platform**: Match experimental platforms when possible

### 2. Method Selection
- **Cell2location**: Gold standard, use when accuracy is critical
- **Spotiphy**: Initial exploration, interactive analysis
- **RCTD**: When doublet detection is important
- **DestVI**: Complex tissues with continuous states
- **Stereoscope**: Rare cell type detection
- **SPOTlight**: When interpretability is key

### 3. Parameter Optimization
- Start with defaults
- Adjust based on convergence
- Validate with known markers
- Compare multiple methods

### 4. Validation Strategies
- Visualize spatial patterns
- Check marker gene expression
- Compare with histology
- Use orthogonal methods

## Performance Considerations

### Speed Comparison
```
Spotiphy     < 1 minute    ⭐⭐⭐⭐⭐
SPOTlight    ~ 5 minutes   ⭐⭐⭐⭐
RCTD         ~ 10 minutes  ⭐⭐⭐
Stereoscope  ~ 30 minutes  ⭐⭐
DestVI       ~ 45 minutes  ⭐⭐
Cell2location ~ 1 hour     ⭐
```

### Memory Usage
- Scales with n_spots × n_cell_types × n_genes
- GPU methods more memory efficient
- Consider subsampling for testing

### GPU Acceleration
- Cell2location: 5-10x speedup
- DestVI: 3-5x speedup
- Stereoscope: 2-3x speedup
- Spotiphy: Minimal benefit

## Troubleshooting

### Common Issues

1. **"No overlapping genes"**
   - Check gene naming (symbols vs IDs)
   - Verify species compatibility
   - Use `.var_names_make_unique()`

2. **"Memory error"**
   - Reduce batch size
   - Subsample spots/genes
   - Use GPU if available

3. **"Convergence failed"**
   - Increase max_epochs
   - Adjust learning rate
   - Check data normalization

4. **"All zeros in output"**
   - Verify reference quality
   - Check normalization
   - Increase detection_alpha

### Debug Mode
```python
# Enable verbose output
params = DeconvolutionParameters(
    method="cell2location",
    verbose=True,
    cell2location_params={
        "train_size": 0.9,  # Use validation
        "early_stopping": True
    }
)
```

## Integration with Other Modules

### Visualization
```python
# Automatic visualization after deconvolution
vis_params = VisualizationParameters(
    plot_type="deconvolution",
    n_cell_types=8,
    colormap="Blues"
)
```

### Spatial Analysis
```python
# Analyze deconvolved cell type patterns
spatial_params = SpatialAnalysisParameters(
    analysis_type="morans_i",
    feature="T_cells_proportion"
)
```

### Cell Communication
```python
# Use proportions for weighted analysis
comm_params = CellCommunicationParameters(
    method="liana",
    use_cell_proportions=True,
    proportions_key="deconvolution_proportions"
)
```

## Quality Control Metrics

### Entropy Score
Measures diversity of cell types per spot:
```
H = -Σ(p_i × log(p_i))
```
- High entropy: Mixed composition
- Low entropy: Dominated by one type

### Sparsity Score
Fraction of near-zero proportions:
- High sparsity: Clear segregation
- Low sparsity: Mixed throughout

### Spatial Coherence
Autocorrelation of proportions:
- High coherence: Spatially organized
- Low coherence: Random distribution

## Future Enhancements

1. **Additional Methods**
   - CARD (Conditional Autoregressive)
   - SpatialDWLS
   - Tangram integration

2. **Multi-modal Integration**
   - Protein co-detection
   - Morphology integration
   - Prior knowledge incorporation

3. **Uncertainty Quantification**
   - Confidence intervals
   - Bootstrap estimates
   - Bayesian posteriors

4. **Optimization**
   - Auto-hyperparameter tuning
   - Ensemble methods
   - Transfer learning