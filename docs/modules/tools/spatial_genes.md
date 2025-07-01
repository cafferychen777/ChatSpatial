# Spatial Variable Genes Tool Documentation

## Overview

The `spatial_genes.py` module implements spatial variable gene identification using GASTON (Generative Adversarial Spatial Transcriptomics Optimization Network) and SpatialDE methods. It identifies genes with significant spatial expression patterns, including continuous gradients and discontinuous transitions.

## Purpose

Spatial variable gene identification is crucial for:
- Discovering spatially organized gene expression patterns
- Identifying tissue organization principles
- Finding genes that define anatomical structures
- Understanding gradients in development and disease
- Detecting boundaries and transition zones
- Revealing hidden spatial axes in tissues

## GASTON Method

### Core Concept

GASTON learns a one-dimensional "isodepth" coordinate that varies smoothly across tissue, revealing the primary axis of spatial variation. This enables identification of genes that follow spatial patterns along this learned coordinate.

### Mathematical Principles

**Dual Neural Network Architecture**:
1. **Spatial Embedding Network (f)**: Maps 2D coordinates to 1D isodepth
   ```
   s = f(x, y) where s ∈ [0, 1]
   ```

2. **Expression Function Network (g)**: Predicts expression from isodepth
   ```
   ê_ij = g_j(s_i) for gene j at location i
   ```

**Joint Optimization**:
```
L = Σ_ij |e_ij - g_j(f(x_i, y_i))|² + λ·R(f)
```
Where R(f) is a regularization term for spatial smoothness

### Isodepth Interpretation

The isodepth coordinate represents:
- A continuous spatial axis (e.g., superficial to deep)
- The primary direction of gene expression variation
- A topographic map of the tissue
- Values range from 0 to 1 across the tissue

## Input Parameters

### SpatialVariableGenesParameters
```python
class SpatialVariableGenesParameters:
    # Preprocessing
    preprocessing_method: str = "glmpca"  # or "pearson_residuals"
    num_hidden_genes: int = 200  # Latent genes for GLM-PCA
    
    # Neural network architecture
    spatial_hidden_layers: List[int] = [32, 16]  # Spatial embedding
    expression_hidden_layers: List[int] = [32, 16]  # Expression prediction
    
    # Training parameters
    epochs: int = 1000
    learning_rate: float = 0.001
    batch_size: int = 256
    validation_split: float = 0.1
    early_stopping_patience: int = 50
    
    # Spatial encoding
    use_positional_encoding: bool = True
    positional_encoding_dim: int = 4
    
    # Analysis parameters
    n_domains: int = 5  # Spatial domains to identify
    num_bins: int = 70  # Bins for isodepth analysis
    continuous_quantile: float = 0.9  # Threshold for continuous
    discontinuous_quantile: float = 0.9  # Threshold for discontinuous
    umi_threshold: int = 500  # Min UMI for gene selection
    
    # Technical parameters
    device: str = "auto"  # cuda/cpu/auto
    random_seed: int = 42
    verbose: bool = True
```

## Preprocessing Methods

### 1. GLM-PCA (`glmpca`)
**Description**: Generalized linear model PCA for count data

**Advantages**:
- Handles overdispersion
- Appropriate for UMI counts
- Reduces technical noise

**Process**:
1. Fit GLM-PCA model on gene expression
2. Extract latent factors
3. Use residuals for spatial analysis

### 2. Pearson Residuals (`pearson_residuals`)
**Description**: Variance-stabilizing transformation

**Advantages**:
- Simple and fast
- Well-established method
- Good for initial exploration

**Process**:
1. Compute expected counts under null model
2. Calculate Pearson residuals
3. Clip extreme values

## Output Format

### SpatialVariableGenesResult
```python
class SpatialVariableGenesResult:
    # Core results
    n_continuous_genes: int  # Genes with gradients
    n_discontinuous_genes: int  # Genes with jumps
    continuous_genes: List[str]  # Gene names
    discontinuous_genes: List[str]  # Gene names
    
    # Spatial information
    n_spatial_domains: int  # Identified domains
    isodepth_key: str = "gaston_isodepth"  # In adata.obs
    spatial_domains_key: str = "gaston_domains"  # In adata.obs
    
    # Gene analysis
    gene_r2_values: Dict[str, float]  # Fit quality
    gene_max_lfc: Dict[str, float]  # Max fold change
    gene_patterns: Dict[str, str]  # Pattern type
    
    # Model performance
    final_loss: float
    training_history: Dict[str, List[float]]
    model_performance: Dict[str, float] = {
        "r2": float,  # Overall R-squared
        "pearson_correlation": float
    }
    
    # Storage
    model_key: str = "gaston_model"  # In adata.uns
    predictions_key: str = "gaston_predictions"  # In adata.obsm
```

## Implementation Details

### Neural Network Architecture

#### Spatial Embedding Network
```python
Input: 2D coordinates (x, y)
↓
Positional Encoding (optional)
↓
Hidden Layer 1 (32 units, ReLU)
↓
Hidden Layer 2 (16 units, ReLU)
↓
Output: 1D isodepth (sigmoid activation)
```

#### Expression Prediction Network
```python
Input: Isodepth value (s)
↓
Hidden Layer 1 (32 units, ReLU)
↓
Hidden Layer 2 (16 units, ReLU)
↓
Output: Gene expression (linear)
```

### Training Process

1. **Data Preparation**
   - Normalize spatial coordinates
   - Preprocess gene expression
   - Create train/validation split

2. **Model Training**
   - Joint optimization of both networks
   - Early stopping on validation loss
   - Learning rate scheduling

3. **Post-Training Analysis**
   - Bin spots by isodepth
   - Fit piecewise linear models
   - Classify gene patterns

### Gene Pattern Classification

#### Continuous Genes
- Show smooth gradients along isodepth
- High R² with isodepth
- Examples: morphogen gradients, metabolic zones

#### Discontinuous Genes
- Show sharp transitions
- High fold change between adjacent bins
- Examples: layer markers, boundary genes

## Usage Examples

### Example 1: Basic GASTON Analysis
```python
# Standard spatial gene discovery
result = await find_spatial_genes(
    data_id="spatial_data",
    params=SpatialVariableGenesParameters(
        preprocessing_method="glmpca",
        epochs=1000,
        n_domains=5
    )
)

print(f"Found {result.n_continuous_genes} continuous genes")
print(f"Found {result.n_discontinuous_genes} discontinuous genes")
```

### Example 2: High-Resolution Analysis
```python
# For Visium HD or dense data
result = await find_spatial_genes(
    data_id="hd_data",
    params=SpatialVariableGenesParameters(
        spatial_hidden_layers=[64, 32, 16],  # Deeper network
        expression_hidden_layers=[64, 32, 16],
        num_bins=100,  # More bins
        epochs=2000,
        batch_size=512
    )
)
```

### Example 3: Developmental Gradient Analysis
```python
# Identify developmental genes
result = await find_spatial_genes(
    data_id="embryo_data",
    params=SpatialVariableGenesParameters(
        preprocessing_method="glmpca",
        continuous_quantile=0.95,  # Strict threshold
        n_domains=3,  # Major regions
        positional_encoding_dim=8  # Higher resolution
    )
)

# Get top gradient genes
gradient_genes = sorted(
    result.continuous_genes,
    key=lambda g: result.gene_r2_values[g],
    reverse=True
)[:20]
```

### Example 4: Layer-specific Analysis
```python
# For layered tissues (e.g., cortex)
result = await find_spatial_genes(
    data_id="cortex_data",
    params=SpatialVariableGenesParameters(
        discontinuous_quantile=0.95,  # Find layer markers
        num_bins=50,
        n_domains=6,  # Cortical layers
        umi_threshold=1000  # Well-expressed genes
    )
)

# Layer marker genes
layer_markers = result.discontinuous_genes
```

### Example 5: Focused Gene Discovery
```python
# Find specific pattern types
result = await find_spatial_genes(
    data_id="data_1",
    params=SpatialVariableGenesParameters(
        preprocessing_method="pearson_residuals",  # Faster
        epochs=500,  # Quick analysis
        continuous_quantile=0.8,  # Relaxed
        discontinuous_quantile=0.99  # Very strict
    )
)

# Strong discontinuous patterns only
strong_boundaries = [
    g for g in result.discontinuous_genes
    if result.gene_max_lfc[g] > 2
]
```

## Best Practices

### 1. Data Preparation
- Ensure sufficient spatial coverage
- Remove low-quality spots
- Consider batch effects
- Check spatial coordinate scaling

### 2. Parameter Selection

#### Network Architecture
- Deeper networks for complex patterns
- Smaller networks for simple gradients
- Match complexity to data size

#### Training Parameters
- More epochs for larger datasets
- Early stopping prevents overfitting
- Adjust batch size for GPU memory

#### Analysis Thresholds
- Higher quantiles for stringent selection
- Lower UMI threshold for rare genes
- More bins for fine-grained analysis

### 3. Interpretation

#### Isodepth Values
- 0 to 1 represents spatial progression
- Smooth transitions indicate gradients
- Sharp changes suggest boundaries

#### Gene Classifications
- Validate with known spatial markers
- Visualize patterns for confirmation
- Consider biological plausibility

### 4. Validation
```python
# Visualize isodepth
vis_params = VisualizationParameters(
    plot_type="gaston_isodepth"
)

# Visualize top genes
vis_params = VisualizationParameters(
    plot_type="gaston_genes",
    features=result.continuous_genes[:6]
)
```

## Advanced Features

### 1. Custom Preprocessing
```python
# Use subset of genes
params = SpatialVariableGenesParameters(
    gene_subset=marker_genes,  # Focus on known markers
    preprocessing_method="custom"
)
```

### 2. Multi-scale Analysis
```python
# Hierarchical spatial patterns
scales = [30, 50, 100]  # Different bin numbers
multi_scale_results = {}

for scale in scales:
    result = await find_spatial_genes(
        data_id="data_1",
        params=SpatialVariableGenesParameters(
            num_bins=scale
        )
    )
    multi_scale_results[scale] = result
```

### 3. Spatial Constraints
```python
# Focus on specific regions
params = SpatialVariableGenesParameters(
    spatial_subset=region_mask,  # Boolean mask
    n_domains=3
)
```

## Troubleshooting

### Common Issues

1. **"No spatial patterns found"**
   - Check data quality
   - Reduce quantile thresholds
   - Increase UMI threshold
   - Try different preprocessing

2. **"Training loss not decreasing"**
   - Reduce learning rate
   - Check data normalization
   - Simplify network architecture
   - Increase batch size

3. **"Memory error"**
   - Reduce batch size
   - Use fewer genes
   - Simplify network
   - Use CPU if GPU limited

4. **"Poor isodepth learning"**
   - Check spatial coordinates
   - Adjust positional encoding
   - Try different architectures
   - Verify tissue has structure

### Debug Mode
```python
params = SpatialVariableGenesParameters(
    verbose=True,
    save_model=True,
    plot_training=True,
    log_interval=100
)
```

## Performance Considerations

### Computational Requirements
- GPU recommended for large datasets
- Training time: 5-30 minutes typical
- Memory: Scales with n_spots × n_genes

### Optimization Tips
- Start with subset for tuning
- Use GLM-PCA for better results
- Cache preprocessing results
- Parallelize gene analysis

## Integration with Other Modules

### SpatialDE Integration
The module also supports SpatialDE for comparison:
```python
# Run both methods
gaston_result = await find_spatial_genes(
    data_id="data_1",
    params=SpatialVariableGenesParameters(method="gaston")
)

# Future: SpatialDE comparison
# spatialde_result = await find_spatial_genes(
#     data_id="data_1", 
#     params=SpatialVariableGenesParameters(method="spatialde")
# )
```

### Visualization Integration
```python
# Three visualization types available
# 1. Isodepth coordinate
vis_params = VisualizationParameters(
    plot_type="gaston_isodepth"
)

# 2. Spatial domains
vis_params = VisualizationParameters(
    plot_type="gaston_domains"
)

# 3. Gene patterns
vis_params = VisualizationParameters(
    plot_type="gaston_genes",
    features=result.continuous_genes[:9]
)
```

### Downstream Analysis
```python
# Use spatial genes for enrichment
await analyze_enrichment(
    data_id="data_1",
    gene_sets={
        "continuous_spatial": result.continuous_genes,
        "discontinuous_spatial": result.discontinuous_genes
    }
)
```

## Biological Interpretation

### Continuous Patterns
- Morphogen gradients
- Metabolic zones
- Developmental axes
- Microenvironment gradients

### Discontinuous Patterns
- Tissue boundaries
- Cell type transitions
- Anatomical layers
- Disease margins

### Domain Interpretation
- Functional zones
- Cellular neighborhoods  
- Microanatomical structures
- Pathological regions

## Future Enhancements

1. **Additional Methods**
   - SpatialDE integration
   - Hotspot implementation
   - SPARK integration

2. **Advanced Features**
   - Multi-dimensional isodepth
   - Temporal spatial patterns
   - 3D tissue support

3. **Performance**
   - Distributed training
   - Approximate methods
   - Transfer learning

4. **Interpretation**
   - Automatic annotation
   - Pattern databases
   - Interactive exploration