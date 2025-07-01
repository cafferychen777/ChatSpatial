# Spatial Domains Tool Documentation

## Overview

The `spatial_domains.py` module provides advanced methods for identifying spatial domains in spatial transcriptomics data. Spatial domains are regions of tissue with similar gene expression patterns and spatial proximity, representing functional units like anatomical layers, microenvironments, or cell type niches.

## Purpose

Spatial domain identification is crucial for:
- Discovering tissue architecture and organization
- Identifying functional regions and microenvironments
- Detecting spatial patterns beyond individual genes
- Enabling region-specific downstream analysis
- Understanding tissue heterogeneity at scale

## Available Methods

### 1. SpaGCN (`spagcn`)

**Description**: Spatial domain detection using graph convolutional networks with optional histology integration.

**Key Features**:
- Integrates gene expression with spatial location
- Optional histology image incorporation
- Adaptive spatial weight calculation
- Domain refinement capabilities

**Mathematical Principle**:
- Constructs spatial neighbor graph
- Applies graph convolutional layers
- Uses unsupervised clustering on embeddings
- Refines boundaries using spatial constraints

### 2. STAGATE (`stagate`)

**Description**: Spatial domain identification via graph attention auto-encoder.

**Key Features**:
- Self-supervised learning approach
- Attention mechanism for neighbor importance
- Preserves both local and global structure
- GPU acceleration support

**Mathematical Principle**:
- Graph attention networks (GAT) for feature learning
- Autoencoder for dimensionality reduction
- Adaptive attention weights for spatial neighbors
- KL divergence optimization

### 3. BANKSY (`banksy`)

**Description**: Building Aggregated Neighborhood Kernel matriceS methodologY for multi-scale analysis.

**Key Features**:
- Multi-scale neighborhood aggregation
- Balances cellular and neighborhood expression
- Azimuthal and radial kernels
- Gradient detection capabilities

**Mathematical Principle**:
- Augments expression matrix with neighborhood features
- Uses Gaussian/exponential spatial kernels
- Combines multiple spatial scales
- Standard clustering on augmented features

### 4. Leiden/Louvain (`leiden`/`louvain`)

**Description**: Standard graph-based clustering with spatial constraints.

**Key Features**:
- Fast and scalable
- Resolution parameter control
- Works with spatial neighborhood graph
- No additional dependencies

**Mathematical Principle**:
- Modularity optimization
- Spatial k-NN graph construction
- Community detection algorithms

## Input Parameters

### SpatialDomainParameters
```python
class SpatialDomainParameters:
    # Method selection
    method: str = "stagate"  # spagcn, stagate, banksy, leiden, louvain
    
    # Common parameters
    n_domains: Optional[int] = None  # Number of domains
    alpha: float = 0.5  # Spatial weight (0-1)
    n_neighbors: int = 10  # Spatial neighbors
    random_seed: int = 42  # Reproducibility
    
    # SpaGCN specific
    spagcn_params: Dict = {
        "beta": 49,  # Spatial penalty
        "refine": True,  # Refine domains
        "spatial_radius": None,  # Auto-calculated
        "histology_image_path": None,  # Optional H&E
        "histology_seg_method": "louvain",
        "n_clusters": None  # Override n_domains
    }
    
    # STAGATE specific
    stagate_params: Dict = {
        "hidden_dims": [512, 30],  # Network architecture
        "n_epochs": 1000,  # Training epochs
        "learning_rate": 0.001,
        "weight_decay": 0.0001,
        "gradient_clipping": 5.0,
        "use_gpu": True,
        "verbose": True
    }
    
    # BANKSY specific
    banksy_params: Dict = {
        "k_geom": 15,  # Neighbors for BANKSY
        "lambda_param": 0.2,  # Cell vs neighbor weight
        "max_m": 1,  # Max order of kernels
        "sigma": 1.5,  # Spatial decay
        "pca_dims": 20,  # PCA components
        "annotation_key": None  # Cell types
    }
    
    # Graph clustering specific
    graph_params: Dict = {
        "resolution": 1.0,  # Clustering resolution
        "n_iterations": -1,  # Auto-determine
        "use_weights": True  # Edge weights
    }
```

## Output Format

### SpatialDomainResult
```python
class SpatialDomainResult:
    # Core results
    n_domains: int  # Number of identified domains
    domain_key: str  # Key in adata.obs
    
    # Optional refined domains
    refined_domain_key: Optional[str]
    
    # Embeddings
    embedding_key: Optional[str]  # Key in adata.obsm
    
    # Quality metrics
    statistics: Dict[str, Any] = {
        "silhouette_score": float,  # Clustering quality
        "calinski_harabasz": float,  # Separation metric
        "davies_bouldin": float,  # Compactness
        "spatial_coherence": float,  # Custom metric
        "domain_sizes": Dict[str, int]
    }
    
    # Method metadata
    method_used: str
    parameters_used: Dict[str, Any]
    
    # Convergence info
    convergence_info: Optional[Dict]
```

## Implementation Details

### Preprocessing Pipeline

1. **Data Validation**
   - Check spatial coordinates
   - Verify gene expression matrix
   - Handle missing values

2. **Feature Selection**
   - Use highly variable genes
   - Optional PCA reduction
   - Normalize if needed

3. **Spatial Graph Construction**
   - K-nearest neighbors
   - Radius-based neighbors
   - Edge weight calculation

4. **Method-Specific Processing**
   - Initialize parameters
   - Run core algorithm
   - Post-process results

### Method-Specific Implementations

#### SpaGCN Implementation
```python
# 1. Calculate spatial weight matrix
# 2. Integrate histology if provided
# 3. Run SpaGCN clustering
# 4. Optional domain refinement
# 5. Compute spatial metrics
```

#### STAGATE Implementation
```python
# 1. Build spatial graph
# 2. Initialize GAT autoencoder
# 3. Train with reconstruction loss
# 4. Extract embeddings
# 5. Cluster in latent space
```

#### BANKSY Implementation
```python
# 1. Compute neighborhood matrices
# 2. Create augmented feature matrix
# 3. Apply dimensionality reduction
# 4. Perform clustering
# 5. Multi-scale analysis
```

## Usage Examples

### Example 1: Brain Tissue with Layers
```python
# Identify cortical layers in brain
result = await identify_spatial_domains(
    data_id="brain_data",
    params=SpatialDomainParameters(
        method="stagate",
        n_domains=6,  # Expected layers
        alpha=0.8,  # High spatial weight
        stagate_params={
            "hidden_dims": [512, 30],
            "n_epochs": 2000
        }
    )
)
```

### Example 2: Tumor Microenvironment
```python
# Detect tumor regions and boundaries
result = await identify_spatial_domains(
    data_id="tumor_data",
    params=SpatialDomainParameters(
        method="spagcn",
        n_domains=5,  # Tumor, stroma, immune, etc.
        spagcn_params={
            "beta": 100,  # Strong spatial penalty
            "refine": True,  # Refine boundaries
            "histology_image_path": "/path/to/HE.tif"
        }
    )
)
```

### Example 3: Large-Scale Atlas
```python
# Efficient processing for atlas data
result = await identify_spatial_domains(
    data_id="atlas_data",
    params=SpatialDomainParameters(
        method="banksy",
        n_domains=20,  # Many regions
        banksy_params={
            "lambda_param": 0.3,
            "pca_dims": 50,  # More components
            "k_geom": 30  # Larger neighborhoods
        }
    )
)
```

### Example 4: Quick Exploration
```python
# Fast leiden clustering
result = await identify_spatial_domains(
    data_id="data_1",
    params=SpatialDomainParameters(
        method="leiden",
        graph_params={
            "resolution": 0.8,
            "n_iterations": 2
        }
    )
)
```

### Example 5: Multi-Method Comparison
```python
methods = ["stagate", "spagcn", "banksy"]
results = {}

for method in methods:
    results[method] = await identify_spatial_domains(
        data_id="data_1",
        params=SpatialDomainParameters(
            method=method,
            n_domains=7,
            alpha=0.5
        )
    )
    
# Compare spatial coherence scores
for method, result in results.items():
    print(f"{method}: {result.statistics['spatial_coherence']}")
```

## Best Practices

### 1. Method Selection

- **STAGATE**: Best overall performance, good for most tissues
- **SpaGCN**: When histology is available, boundary refinement needed
- **BANKSY**: Multi-scale patterns, gradient detection
- **Leiden**: Quick exploration, baseline comparison

### 2. Parameter Tuning

#### Number of Domains
- Start with expected anatomical regions
- Use silhouette analysis to optimize
- Visualize and adjust based on biology

#### Spatial Weight (alpha)
- 0.0-0.3: Expression-driven (functional domains)
- 0.4-0.6: Balanced (default)
- 0.7-1.0: Spatially-driven (anatomical domains)

#### Resolution (Graph methods)
- Lower (0.1-0.5): Fewer, larger domains
- Medium (0.8-1.2): Standard
- Higher (1.5-3.0): Many small domains

### 3. Validation Strategies

- Compare with histology annotations
- Check marker gene enrichment
- Validate spatial coherence
- Test parameter sensitivity

### 4. Common Workflows

```python
# 1. Exploratory analysis
leiden_result = await identify_spatial_domains(
    data_id, SpatialDomainParameters(method="leiden")
)

# 2. Refined analysis
stagate_result = await identify_spatial_domains(
    data_id, 
    SpatialDomainParameters(
        method="stagate",
        n_domains=leiden_result.n_domains
    )
)

# 3. Visualization
vis_params = VisualizationParameters(
    plot_type="spatial_domains",
    domain_key=stagate_result.domain_key
)
```

## Performance Considerations

### Time Complexity
- Leiden: O(n log n) - fastest
- BANKSY: O(n × k × m) - moderate
- SpaGCN: O(n² × iterations) - slower
- STAGATE: O(n × epochs × k) - GPU helps

### Memory Usage
- Scales with n_spots × n_features
- Spatial graphs can be large
- GPU methods need additional memory

### Optimization Tips
- Subsample for parameter tuning
- Use GPU for STAGATE
- Reduce features with PCA
- Adjust n_neighbors for density

## Troubleshooting

### Common Issues

1. **"Import error for method"**
   - Install method-specific package
   - Check Python version compatibility
   - Use conda environments

2. **"No spatial coordinates"**
   - Verify spatial_key parameter
   - Check coordinate format
   - Ensure data has spatial info

3. **"Timeout/Memory error"**
   - Reduce n_neighbors
   - Subsample spots
   - Use fewer features
   - Try faster method

4. **"Poor domain separation"**
   - Adjust alpha parameter
   - Try different methods
   - Check data quality
   - Increase n_epochs

### Debug Mode
```python
# Enable detailed logging
params = SpatialDomainParameters(
    method="stagate",
    stagate_params={
        "verbose": True,
        "save_model": True,
        "log_interval": 100
    }
)
```

## Quality Metrics

### Silhouette Score
- Measures cluster separation
- Range: [-1, 1], higher is better
- > 0.5: Good separation
- < 0.2: Poor separation

### Spatial Coherence
- Custom metric for spatial continuity
- Measures neighbor similarity
- Range: [0, 1], higher is better

### Domain Size Distribution
- Check for extreme imbalances
- Small domains may be artifacts
- Large domains may need splitting

## Integration with Other Modules

### Visualization
```python
# Visualize domains
vis_params = VisualizationParameters(
    plot_type="spatial_domains",
    domain_key=result.domain_key,
    show_boundaries=True
)
```

### Differential Expression
```python
# Find domain markers
de_result = await find_markers(
    data_id="data_1",
    group_key=result.domain_key,
    group1="Domain_1",
    group2="Domain_2"
)
```

### Cell Communication
```python
# Domain-specific communication
comm_params = CellCommunicationParameters(
    method="liana",
    groupby=result.domain_key,
    spatial_key="spatial"
)
```

## Advanced Features

### 1. Hierarchical Domains
```python
# Multi-resolution analysis
resolutions = [0.5, 1.0, 2.0]
hierarchical_domains = {}

for res in resolutions:
    result = await identify_spatial_domains(
        data_id="data_1",
        params=SpatialDomainParameters(
            method="leiden",
            graph_params={"resolution": res}
        )
    )
    hierarchical_domains[res] = result
```

### 2. Constrained Clustering
```python
# Use prior knowledge
params = SpatialDomainParameters(
    method="banksy",
    banksy_params={
        "annotation_key": "cell_type",  # Guide by cell types
        "lambda_param": 0.4
    }
)
```

### 3. Domain Refinement
```python
# Post-process domains
params = SpatialDomainParameters(
    method="spagcn",
    spagcn_params={
        "refine": True,
        "beta": 150,  # Strong smoothing
        "n_iterations": 10
    }
)
```

## Future Enhancements

1. **Additional Methods**
   - BayesSpace integration
   - Hidden Markov Random Fields
   - Deep learning approaches

2. **Multi-modal Integration**
   - Morphology features
   - Protein co-detection
   - Chromatin accessibility

3. **Temporal Domains**
   - Time-series analysis
   - Dynamic domain tracking
   - Developmental trajectories

4. **3D Domains**
   - 3D tissue reconstruction
   - Volume-based clustering
   - Cross-section integration