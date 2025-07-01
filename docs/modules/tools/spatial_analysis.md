# Spatial Analysis Tool Documentation

## Overview

The `spatial_analysis.py` module provides comprehensive spatial statistics and pattern analysis capabilities for spatial transcriptomics data. It implements multiple spatial analysis methods to identify patterns, detect hotspots, analyze neighborhoods, and quantify spatial relationships in gene expression and cell type distributions.

## Purpose

Spatial analysis is crucial for understanding the spatial organization of tissues. This module enables researchers to:
- Detect spatial patterns in gene expression
- Identify cell type neighborhoods and interactions
- Find hotspots and cold spots of activity
- Quantify spatial autocorrelation
- Analyze spatial point patterns
- Measure network centrality in spatial graphs

## Architecture

### Core Components

1. **Main Entry Point**: `analyze_spatial_patterns()` function
2. **Method Handlers**: Individual functions for each analysis type
3. **Helper Functions**: Utilities for spatial computations
4. **Error Handling**: Robust error management with informative messages

### Key Features

- **Unified Interface**: Single entry point for all spatial analyses
- **Flexible Parameters**: Extensive customization options
- **Automatic Detection**: Intelligent handling of data structures
- **Comprehensive Output**: Detailed statistics and visualization-ready results

## Available Analysis Methods

### 1. Neighborhood Enrichment Analysis (`neighborhood`)

**Purpose**: Identifies which cell types tend to co-localize in spatial neighborhoods.

**Mathematical Basis**: 
- Computes enrichment score as observed/expected frequency ratio
- Uses permutation testing for statistical significance

**Key Features**:
- Analyzes cell type co-occurrence patterns
- Provides Z-scores and p-values
- Creates enrichment matrix for visualization

**Output**:
- Enrichment scores matrix
- Statistical significance
- Neighborhood composition statistics

### 2. Co-occurrence Analysis (`co_occurrence`)

**Purpose**: Measures spatial association between cell types at different distances.

**Mathematical Basis**:
- Analyzes occurrence and co-occurrence frequencies
- Tests against spatial randomness null hypothesis
- Interval-based distance analysis

**Key Features**:
- Multi-scale analysis (different distance intervals)
- Handles both occurrence and co-occurrence
- Statistical testing with multiple testing correction

**Output**:
- Co-occurrence scores by distance
- P-values for each cell type pair
- Distance interval results

### 3. Ripley's L Function (`ripley`)

**Purpose**: Tests whether spatial point patterns show clustering, randomness, or regularity.

**Mathematical Basis**:
```
L(r) = √(K(r)/π)
K(r) = A/n² × Σᵢⱼ I(dᵢⱼ ≤ r)
```
Where K(r) is Ripley's K function, A is area, n is number of points

**Key Features**:
- Multi-distance analysis
- Edge correction for boundary effects
- Comparison to complete spatial randomness (CSR)

**Output**:
- L-statistic values at different radii
- Expected values under CSR
- Clustering strength indicators

### 4. Moran's I Spatial Autocorrelation (`morans_i`)

**Purpose**: Quantifies global spatial autocorrelation in gene expression.

**Mathematical Basis**:
```
I = (n/W) × Σᵢⱼ wᵢⱼ(xᵢ - x̄)(xⱼ - x̄) / Σᵢ(xᵢ - x̄)²
```
Where w is spatial weight, x is expression value, W is sum of weights

**Key Features**:
- Global autocorrelation measure
- Statistical significance testing
- Handles continuous variables

**Output**:
- Moran's I statistic
- P-value
- Interpretation (clustered/dispersed/random)

### 5. Centrality Scores (`centrality`)

**Purpose**: Identifies spatially important cells based on network topology.

**Types**:
- **Degree Centrality**: Number of spatial neighbors
- **Closeness Centrality**: Average distance to all other cells
- **Betweenness Centrality**: Frequency on shortest paths
- **Eigenvector Centrality**: Importance based on neighbor importance

**Output**:
- Centrality scores per cell
- Summary statistics
- Highly central cells

### 6. Getis-Ord Gi* Hot Spot Analysis (`getis_ord`)

**Purpose**: Identifies statistically significant hot spots and cold spots of gene expression.

**Mathematical Basis**:
```
Gi* = (Σⱼ wᵢⱼxⱼ - X̄Σⱼ wᵢⱼ) / (S√[(nΣⱼ wᵢⱼ² - (Σⱼ wᵢⱼ)²)/(n-1)])
```

**Key Features**:
- Local spatial clustering detection
- Statistical significance with FDR correction
- Works with gene expression values

**Output**:
- Z-scores indicating hot/cold spots
- P-values with multiple testing correction
- Binary hot/cold spot classification

## Input Parameters

### SpatialAnalysisParameters
```python
class SpatialAnalysisParameters:
    # Core parameters
    analysis_type: str  # Type of analysis
    cluster_key: Optional[str] = None  # Cell type/cluster column
    n_neighbors: int = 30  # Number of spatial neighbors
    
    # Method-specific parameters
    # Moran's I
    morans_i_gene: Optional[str] = None  # Gene for autocorrelation
    
    # Getis-Ord
    getis_ord_genes: Optional[List[str]] = None  # Genes to analyze
    fdr_correction: bool = True  # Apply FDR correction
    
    # Ripley's L
    ripley_radii: Optional[List[float]] = None  # Distance radii
    
    # Centrality
    centrality_type: str = "degree"  # Type of centrality
    
    # General
    spatial_key: str = "spatial"  # Coordinates key
    random_seed: int = 42  # For reproducibility
```

## Implementation Details

### Spatial Neighbor Computation
1. Uses k-nearest neighbors or radius-based methods
2. Constructs spatial connectivity graph
3. Handles edge effects and boundary corrections

### Cluster Key Detection
The module intelligently detects cluster/cell type keys:
1. Checks user-specified `cluster_key`
2. Falls back to common keys: 'cell_type', 'celltype', 'leiden', 'louvain'
3. Creates default clustering if none found

### Error Handling
- Validates required parameters for each method
- Checks data availability (coordinates, features)
- Provides informative error messages
- Implements fallback strategies

### Performance Optimization
- Efficient neighbor search using spatial indices
- Vectorized computations where possible
- Memory-efficient sparse matrix operations
- Parallel processing for permutation tests

## Usage Examples

### Example 1: Basic Neighborhood Analysis
```python
# Analyze cell type neighborhoods
params = SpatialAnalysisParameters(
    analysis_type="neighborhood",
    cluster_key="cell_type",
    n_neighbors=30
)
result = await analyze_spatial_data("data_1", params)

# Enrichment scores available in result.statistics
```

### Example 2: Gene Expression Hotspots
```python
# Find hotspots for immune genes
params = SpatialAnalysisParameters(
    analysis_type="getis_ord",
    getis_ord_genes=["CD3E", "CD8A", "GZMB"],
    fdr_correction=True
)
result = await analyze_spatial_data("data_1", params)

# Hot/cold spots stored in adata.obs
```

### Example 3: Spatial Autocorrelation
```python
# Test spatial clustering of a marker gene
params = SpatialAnalysisParameters(
    analysis_type="morans_i",
    morans_i_gene="MKI67",
    n_neighbors=50
)
result = await analyze_spatial_data("data_1", params)

# Moran's I statistic and p-value in result.statistics
```

### Example 4: Multi-Scale Pattern Analysis
```python
# Ripley's L at multiple scales
params = SpatialAnalysisParameters(
    analysis_type="ripley",
    cluster_key="cell_type",
    ripley_radii=[50, 100, 200, 500]
)
result = await analyze_spatial_data("data_1", params)
```

### Example 5: Network Centrality
```python
# Find spatially important cells
params = SpatialAnalysisParameters(
    analysis_type="centrality",
    centrality_type="betweenness",
    n_neighbors=20
)
result = await analyze_spatial_data("data_1", params)
```

### Example 6: Comprehensive Spatial Analysis
```python
# Run multiple analyses in sequence
analyses = [
    ("neighborhood", {"cluster_key": "cell_type"}),
    ("morans_i", {"morans_i_gene": "TNF"}),
    ("getis_ord", {"getis_ord_genes": ["IL6", "IL1B"]}),
    ("ripley", {"cluster_key": "cell_type"})
]

results = {}
for analysis_type, extra_params in analyses:
    params = SpatialAnalysisParameters(
        analysis_type=analysis_type,
        **extra_params
    )
    results[analysis_type] = await analyze_spatial_data("data_1", params)
```

## Output Interpretation

### Neighborhood Enrichment
- **Positive scores**: Cell types co-occur more than expected
- **Negative scores**: Cell types avoid each other
- **Near zero**: Random spatial distribution

### Moran's I
- **I > 0**: Positive autocorrelation (clustering)
- **I < 0**: Negative autocorrelation (dispersion)
- **I ≈ 0**: Random spatial pattern
- **P-value**: Statistical significance

### Getis-Ord Gi*
- **High positive Z-score**: Hot spot
- **High negative Z-score**: Cold spot
- **Z-score near 0**: No significant clustering
- **FDR < 0.05**: Statistically significant

### Ripley's L
- **L(r) > r**: Clustering at distance r
- **L(r) < r**: Regularity at distance r
- **L(r) ≈ r**: Complete spatial randomness

## Best Practices

### 1. Method Selection
- Use `neighborhood` for cell type interaction analysis
- Use `getis_ord` for gene expression hotspot detection
- Use `morans_i` for global pattern testing
- Use `ripley` for multi-scale pattern analysis

### 2. Parameter Tuning
- Adjust `n_neighbors` based on tissue density
- Use multiple radii for comprehensive analysis
- Apply FDR correction for multiple genes
- Set appropriate spatial keys

### 3. Validation
- Visualize results spatially
- Compare with biological expectations
- Test parameter sensitivity
- Check statistical assumptions

### 4. Performance
- Start with subset for parameter tuning
- Use appropriate n_neighbors (typically 20-50)
- Consider computational cost for large datasets
- Cache results for iterative analysis

## Integration with Other Tools

### Visualization
Results integrate seamlessly with visualization tool:
```python
# Visualize neighborhood enrichment
vis_params = VisualizationParameters(
    plot_type="spatial_analysis",
    analysis_sub_type="neighborhood"
)

# Visualize hotspots
vis_params = VisualizationParameters(
    plot_type="spatial_analysis",
    analysis_sub_type="getis_ord",
    feature="CD8A"
)
```

### Downstream Analysis
- Use centrality scores for cell selection
- Apply hotspot results for region definition
- Incorporate neighborhood info in cell communication
- Combine with trajectory analysis

## Troubleshooting

### Common Issues

1. **"No cluster key found"**
   - Run cell type annotation first
   - Specify cluster_key parameter
   - Check adata.obs columns

2. **"Spatial coordinates not found"**
   - Ensure spatial data is loaded correctly
   - Check spatial_key parameter
   - Verify coordinate format

3. **"Memory error"**
   - Reduce n_neighbors
   - Subsample data
   - Use sparse matrix operations

4. **"No significant results"**
   - Check data quality
   - Adjust parameters
   - Try different methods

## Performance Considerations

### Computational Complexity
- Neighborhood: O(n × k)
- Moran's I: O(n²) for full matrix
- Getis-Ord: O(n × k × g) for g genes
- Ripley's: O(n² × r) for r radii

### Memory Usage
- Scales with n_neighbors
- Sparse matrices for efficiency
- Chunked processing available

### Optimization Tips
- Pre-compute spatial graphs
- Use spatial indices
- Parallelize gene-wise operations
- Cache intermediate results

## Future Enhancements

1. **Additional Methods**
   - Lee's L (bivariate association)
   - Spatial scan statistics
   - Kernel density estimation
   - Spatial regression models

2. **Performance**
   - GPU acceleration
   - Distributed computing
   - Approximate algorithms

3. **Integration**
   - 3D spatial analysis
   - Multi-scale frameworks
   - Spatial-temporal analysis

4. **Visualization**
   - Interactive spatial plots
   - 3D visualization
   - Animation support