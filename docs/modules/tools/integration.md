# Integration Tool Documentation

## Overview

The `integration.py` module provides multi-sample integration capabilities for spatial transcriptomics data. It enables the joint analysis of multiple samples while preserving biological variation and removing technical batch effects.

## Purpose

Multi-sample integration is essential for:
- Combining data from multiple experiments or conditions
- Removing batch effects while preserving biological signals
- Enabling comparative analysis across samples
- Building comprehensive spatial atlases
- Increasing statistical power through data aggregation
- Harmonizing datasets from different technologies or protocols

## Available Integration Methods

### 1. Harmony (`harmony`)

**Description**: Fast and scalable integration using iterative clustering.

**Key Features**:
- Iterative linear correction
- Preserves local structure
- Handles multiple batch variables
- Memory efficient

**Algorithm**:
1. Soft clustering of cells
2. Within-cluster linear correction
3. Iterative refinement
4. PCA space transformation

**Best for**: Large datasets, multiple batches, initial exploration

### 2. BBKNN (`bbknn`)

**Description**: Batch-balanced k-nearest neighbors graph construction.

**Key Features**:
- Graph-based integration
- No explicit correction
- Preserves cell type connectivity
- Fast for downstream clustering

**Algorithm**:
1. Identify k-nearest neighbors per batch
2. Connect cells across batches
3. Create unified graph
4. No embedding modification

**Best for**: Discrete cell types, clustering-focused analysis

### 3. Scanorama (`scanorama`)

**Description**: Mutual nearest neighbors with panoramic stitching.

**Key Features**:
- Panoramic alignment
- Gene-level correction
- Handles missing cell types
- Robust to outliers

**Algorithm**:
1. Find mutual nearest neighbors
2. Calculate transformation vectors
3. Apply panoramic corrections
4. Generate integrated embedding

**Best for**: Heterogeneous datasets, missing cell types

### 4. MNN/Combat (`mnn`)

**Description**: Mutual nearest neighbors or Combat correction.

**Key Features**:
- Linear model correction
- Parametric/non-parametric options
- Well-established methods
- Flexible implementation

**Algorithm**:
- Combat: Empirical Bayes batch correction
- MNN: Mutual nearest neighbor matching

**Best for**: Traditional batch correction, linear effects

## Input Parameters

### IntegrationParameters
```python
class IntegrationParameters:
    # Method selection
    method: str = "harmony"  # harmony/bbknn/scanorama/mnn
    
    # Common parameters
    batch_key: str = "sample"  # Batch identifier column
    n_pcs: int = 50  # Number of PCs to use
    
    # Harmony specific
    harmony_params: Dict = {
        "max_iter_harmony": 10,  # Max iterations
        "sigma": 0.1,  # Width of soft clusters
        "theta": 1.0,  # Diversity penalty (0-2)
        "lambda": 1.0,  # Ridge regression penalty
        "block_size": 0.05,  # Proportion per block
        "max_iter_kmeans": 20,  # Kmeans iterations
        "epsilon_cluster": 1e-5,  # Convergence tolerance
        "epsilon_harmony": 1e-4  # Harmony tolerance
    }
    
    # BBKNN specific
    bbknn_params: Dict = {
        "neighbors_within_batch": 3,  # Per-batch neighbors
        "n_pcs": 50,  # PCs to use
        "trim": None,  # Outlier trimming
        "approx": True,  # Approximate search
        "metric": "euclidean",  # Distance metric
        "bandwidth": 1.0,  # Gaussian kernel width
        "local_connectivity": 1  # UMAP connectivity
    }
    
    # Scanorama specific
    scanorama_params: Dict = {
        "approx": True,  # Approximate alignment
        "sigma": 15,  # Alignment sigma
        "alpha": 0.1,  # Alignment strength
        "knn": 20,  # k-nearest neighbors
        "verbose": True  # Progress output
    }
    
    # MNN/Combat specific
    mnn_params: Dict = {
        "k": 20,  # Nearest neighbors
        "sigma": 1.0,  # Gaussian bandwidth
        "cos_norm_in": True,  # Cosine normalize input
        "cos_norm_out": True,  # Cosine normalize output
        "svd_dim": 0,  # SVD dimensions (0=auto)
        "var_adj": True,  # Adjust variance
        "compute_angle": False,  # Compute angles
        "mnn_order": None,  # Batch order (auto)
        "svd_mode": "rsvd"  # SVD algorithm
    }
    
    # Additional options
    scale_data: bool = True  # Scale before integration
    copy: bool = False  # Return copy of data
```

## Output Format

### IntegrationResult
```python
class IntegrationResult:
    # Core results
    integrated_data_id: str  # New dataset ID
    n_samples_integrated: int  # Number of samples
    n_cells_total: int  # Total cells
    integration_key: str  # Key for embeddings
    
    # Method info
    method_used: str  # Integration method
    parameters_used: Dict[str, Any]
    
    # Quality metrics
    integration_metrics: Dict[str, float] = {
        "silhouette_score": float,  # Batch mixing
        "ari_score": float,  # Cell type preservation
        "batch_entropy": float,  # Mixing entropy
        "kbet_rejection_rate": float  # Batch effect test
    }
    
    # Sample information
    sample_cell_counts: Dict[str, int]
    common_genes: int  # Shared genes
    
    # Embeddings
    corrected_pca_key: str = "X_pca_integrated"
    umap_key: str = "X_umap_integrated"
```

### Integrated AnnData Structure
```python
adata_integrated = {
    # Observations
    .obs: {
        "sample": Original sample IDs,
        "batch": Batch identifiers,
        "cell_type": Preserved annotations,
        "_indices": Original cell indices
    },
    
    # Embeddings
    .obsm: {
        "X_pca_integrated": Corrected PCA,
        "X_umap_integrated": Integrated UMAP,
        "X_spatial": Combined spatial coords
    },
    
    # Integration metadata
    .uns: {
        "integration": {
            "method": Method used,
            "params": Parameters,
            "metrics": Quality scores,
            "sample_info": Sample metadata
        }
    }
}
```

## Implementation Details

### Data Preprocessing Pipeline

1. **Sample Loading and Validation**
   ```python
   # Load each sample
   # Find common genes
   # Validate spatial coordinates
   # Check batch assignments
   ```

2. **Gene Selection**
   ```python
   # Highly variable genes per sample
   # Intersection or union approach
   # Remove batch-specific genes
   ```

3. **Normalization Alignment**
   ```python
   # Ensure consistent normalization
   # Scale if requested
   # Handle sparse matrices
   ```

### Integration Workflow

1. **Data Concatenation**
   - Vertical stacking of samples
   - Preserve sample metadata
   - Handle spatial coordinates

2. **Batch Effect Removal**
   - Method-specific correction
   - Preserve biological variation
   - Update embeddings

3. **Post-processing**
   - Recompute UMAP
   - Update clustering
   - Calculate metrics

### Spatial Coordinate Handling

```python
# Option 1: Keep original coordinates
spatial_coords = np.vstack([
    adata.obsm['spatial'] for adata in adatas
])

# Option 2: Align coordinates
aligned_coords = align_spatial_coordinates(
    coords_list, 
    method='procrustes'
)

# Option 3: Add sample offset
offset_coords = add_spatial_offset(
    coords_list,
    spacing=1000
)
```

## Usage Examples

### Example 1: Basic Multi-Sample Integration
```python
# Integrate multiple Visium samples
sample_ids = ["sample1", "sample2", "sample3"]
result = await integrate_samples(
    data_ids=sample_ids,
    params=IntegrationParameters(
        method="harmony",
        batch_key="sample",
        harmony_params={"theta": 1.0}
    )
)

print(f"Integrated {result.n_samples_integrated} samples")
print(f"Total cells: {result.n_cells_total}")
```

### Example 2: Complex Batch Structure
```python
# Multiple batch variables
result = await integrate_samples(
    data_ids=sample_ids,
    params=IntegrationParameters(
        method="harmony",
        batch_key=["patient", "batch", "technology"],
        harmony_params={
            "theta": [2.0, 1.0, 1.0],  # Per-variable penalty
            "lambda": [1.0, 1.0, 1.0],
            "max_iter_harmony": 20
        }
    )
)
```

### Example 3: Graph-based Integration
```python
# BBKNN for discrete cell types
result = await integrate_samples(
    data_ids=sample_ids,
    params=IntegrationParameters(
        method="bbknn",
        bbknn_params={
            "neighbors_within_batch": 5,
            "trim": 10,  # Remove outliers
            "approx": False  # Exact computation
        }
    )
)
```

### Example 4: Handling Missing Cell Types
```python
# Scanorama for heterogeneous samples
result = await integrate_samples(
    data_ids=["young_sample", "aged_sample"],
    params=IntegrationParameters(
        method="scanorama",
        scanorama_params={
            "sigma": 20,  # Looser alignment
            "alpha": 0.05,  # Gentle correction
            "knn": 30
        }
    )
)
```

## Best Practices

### 1. Data Preparation
- Ensure consistent preprocessing
- Use same gene set across samples
- Check annotation consistency
- Validate spatial coordinates

### 2. Method Selection

#### Dataset Size
- **Small (<10k cells)**: Any method
- **Medium (10k-100k)**: Harmony, BBKNN
- **Large (>100k)**: Harmony, BBKNN with approximation

#### Batch Structure
- **Simple batches**: Harmony, Combat
- **Complex batches**: Harmony with multiple keys
- **Technical replicates**: Combat
- **Biological replicates**: Gentle Harmony (theta=0.5)

#### Cell Type Distribution
- **All types present**: Any method
- **Missing types**: Scanorama, MNN
- **Rare populations**: BBKNN, careful parameters

### 3. Parameter Tuning

#### Harmony Theta
- 0: No integration
- 0.5-1: Gentle integration (biological replicates)
- 1-2: Standard integration
- 2-4: Strong integration (technical batches)

#### BBKNN Neighbors
- 3-5: Standard
- 10+: Strong integration
- 1-2: Minimal integration

### 4. Quality Control

```python
# Check integration quality
metrics = result.integration_metrics

# Good integration indicators:
# - High silhouette score (>0.5)
# - High ARI score (>0.7)
# - Low KBET rejection (<0.1)
# - High batch entropy (>0.8)
```

## Spatial-Specific Considerations

### 1. Coordinate Alignment
```python
# Align tissue sections
params = IntegrationParameters(
    method="harmony",
    align_spatial=True,
    spatial_alignment_method="procrustes"
)
```

### 2. Regional Integration
```python
# Integrate specific regions
# First subset to regions
region_masks = [adata.obs['region'] == 'cortex' 
                for adata in adatas]
```

### 3. Multi-modal Integration
```python
# Combine different technologies
params = IntegrationParameters(
    method="scanorama",
    handle_missing_genes=True,
    gene_selection="union"
)
```

## Troubleshooting

### Common Issues

1. **"No common genes"**
   - Check gene naming (symbols vs IDs)
   - Verify species compatibility
   - Use `.var_names_make_unique()`

2. **"Memory error"**
   - Reduce n_pcs
   - Use approximation methods
   - Subsample for testing

3. **"Over-correction"**
   - Reduce theta (Harmony)
   - Fewer neighbors (BBKNN)
   - Check biological replicates

4. **"Under-correction"**
   - Increase theta
   - More neighbors
   - Try different method

### Validation Strategies

1. **Biological Validation**
   - Known markers preserved
   - Cell types identifiable
   - Spatial patterns maintained

2. **Technical Validation**
   - Batch mixing in UMAP
   - Consistent clustering
   - Reproducible results

## Performance Considerations

### Time Complexity
- Harmony: O(n × k × iterations)
- BBKNN: O(n × k × b) for b batches
- Scanorama: O(n² × b) worst case
- Combat: O(n × g) for g genes

### Memory Usage
- Scales with total cells
- PCA reduction helps
- Sparse matrices beneficial

### Optimization Tips
- Start with subset
- Use approximations
- Reduce dimensions
- Parallel processing

## Advanced Features

### 1. Hierarchical Integration
```python
# Integrate in stages
# Stage 1: Technical replicates
tech_integrated = await integrate_samples(
    technical_replicates,
    params=IntegrationParameters(
        method="combat",
        batch_key="run"
    )
)

# Stage 2: Biological samples
final_integrated = await integrate_samples(
    [tech_integrated] + other_samples,
    params=IntegrationParameters(
        method="harmony",
        batch_key="condition"
    )
)
```

### 2. Reference-based Integration
```python
# Integrate to reference atlas
params = IntegrationParameters(
    method="scanorama",
    reference_idx=0,  # First sample is reference
    scanorama_params={
        "reference_only": True
    }
)
```

### 3. Trajectory Preservation
```python
# Maintain trajectories
params = IntegrationParameters(
    method="harmony",
    preserve_trajectories=True,
    trajectory_key="dpt_pseudotime"
)
```

## Integration with Other Modules

### Downstream Analysis
```python
# Use integrated data
integrated_id = result.integrated_data_id

# Spatial analysis on integrated data
spatial_result = await analyze_spatial_data(
    data_id=integrated_id,
    params=SpatialAnalysisParameters(
        analysis_type="neighborhood"
    )
)

# Cell communication across samples
comm_result = await analyze_cell_communication(
    data_id=integrated_id,
    params=CellCommunicationParameters(
        groupby="cell_type"
    )
)
```

### Visualization
```python
# Visualize integration
vis_params = VisualizationParameters(
    plot_type="umap",
    color="sample",  # Color by batch
    embedding="X_umap_integrated"
)

# Check cell types
vis_params = VisualizationParameters(
    plot_type="umap", 
    color="cell_type",
    embedding="X_umap_integrated"
)
```

## Future Enhancements

1. **Additional Methods**
   - LIGER integration
   - scVI-based methods
   - Seurat v5 bridge

2. **Spatial Features**
   - Spatial batch effects
   - Image-guided integration
   - 3D alignment

3. **Advanced Options**
   - Online integration
   - Incremental updates
   - Cross-species integration

4. **Performance**
   - GPU acceleration
   - Distributed computing
   - Approximate methods