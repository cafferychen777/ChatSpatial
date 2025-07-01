# Enrichment Analysis Tool Documentation

## Overview

The `enrichment_analysis.py` module provides spatially-aware gene set enrichment analysis using the EnrichMap approach. It enables the identification of biological pathways and gene signatures that are enriched in specific spatial regions, incorporating spatial smoothing and technical covariate correction.

## Purpose

Spatially-aware enrichment analysis is essential for:
- Understanding spatial distribution of biological pathways
- Identifying region-specific functional programs
- Discovering spatial gradients in cellular processes
- Mapping disease-related signatures in tissue
- Validating spatial domains with functional annotation
- Revealing microenvironment-specific activities

## EnrichMap Method

### Core Concept

EnrichMap extends traditional gene set enrichment by incorporating spatial context through three key steps:

1. **Gene Set Scoring**: Calculate enrichment scores for each spot
2. **Spatial Smoothing**: Average scores across spatial neighbors
3. **Covariate Correction**: Remove technical spatial biases using GAM

### Key Advantages

- Reduces noise through spatial smoothing
- Accounts for technical artifacts (edge effects, gradients)
- Preserves biological spatial patterns
- Enables discovery of subtle spatial enrichments
- Supports custom gene weights

## Mathematical Principles

### 1. Basic Enrichment Scoring

For a gene set G with genes {g₁, g₂, ..., gₙ}, the enrichment score for spot i:

```
Score_i = Σⱼ (w_j × Expression_{ij}) / Σⱼ w_j
```

Where:
- w_j is the weight for gene j (default: 1)
- Expression_{ij} is the expression of gene j in spot i

### 2. Spatial Smoothing

The spatially smoothed score incorporates information from k-nearest neighbors:

```
SmoothedScore_i = (Score_i + Σₖ Score_k) / (1 + |N_i|)
```

Where N_i is the set of spatial neighbors for spot i

### 3. GAM Correction

Generalized Additive Model removes technical covariates:

```
CorrectedScore = Score - f(spatial_x, spatial_y) - Σ g(covariate)
```

Where f and g are smooth functions learned by the GAM

## Input Parameters

### Main Function: `perform_enrichment_analysis`

```python
async def perform_enrichment_analysis(
    data_id: str,                          # Dataset identifier
    data_store: Dict[str, Any],            # Data storage
    gene_sets: Union[List[str], Dict[str, List[str]]],  # Gene sets
    score_keys: Optional[Union[str, List[str]]] = None,  # Names
    spatial_key: str = "spatial",          # Spatial coordinates
    n_neighbors: int = 6,                  # Spatial neighbors
    smoothing: bool = True,                # Apply smoothing
    correct_spatial_covariates: bool = True,  # GAM correction
    batch_key: Optional[str] = None,       # Batch column
    context: Context = None                # MCP context
) -> Dict[str, Any]
```

### Helper Functions

#### `compute_enrichment_scores`
```python
def compute_enrichment_scores(
    adata: AnnData,                        # Annotated data
    gene_set: List[str],                   # Genes to score
    weights: Optional[List[float]] = None,  # Gene weights
    method: str = "mean"                   # Scoring method
) -> np.ndarray
```

#### `spatial_smoothing`
```python
def spatial_smoothing(
    scores: np.ndarray,                    # Original scores
    spatial_coords: np.ndarray,            # Coordinates
    n_neighbors: int = 6,                  # Neighbors
    method: str = "knn"                    # Neighbor method
) -> np.ndarray
```

#### `correct_spatial_covariates`
```python
def correct_spatial_covariates(
    scores: np.ndarray,                    # Scores to correct
    spatial_coords: np.ndarray,            # Coordinates
    covariates: Optional[pd.DataFrame] = None,  # Additional
    spline_df: int = 5                    # Spline flexibility
) -> np.ndarray
```

## Output Format

### Result Structure
```json
{
    "signatures": ["signature_names"],
    "enrichment_scores": {
        "signature_name": {
            "raw": [array of raw scores],
            "smoothed": [array of smoothed scores],
            "corrected": [array of corrected scores]
        }
    },
    "summary_stats": {
        "signature_name": {
            "mean": 0.45,
            "std": 0.12,
            "min": -0.23,
            "max": 0.89,
            "n_genes": 25,
            "genes_found": 23,
            "spatial_autocorrelation": 0.65
        }
    },
    "spatial_metrics": {
        "signature_name": {
            "morans_i": 0.42,
            "morans_p": 0.001,
            "hotspots": [spot indices],
            "coldspots": [spot indices]
        }
    },
    "parameters": {
        "n_neighbors": 6,
        "smoothing": true,
        "correction": true
    }
}
```

### Score Interpretation

- **Raw scores**: Direct gene set expression average
- **Smoothed scores**: Spatially averaged values
- **Corrected scores**: Technical artifacts removed

Score ranges:
- Positive: Enrichment above background
- Negative: Depletion below background
- Magnitude: Strength of enrichment/depletion

## Implementation Details

### Data Flow

1. **Gene Validation**
   - Check gene presence in dataset
   - Handle missing genes gracefully
   - Report coverage statistics

2. **Score Computation**
   - Extract expression matrix
   - Apply weights if provided
   - Calculate mean/sum scores

3. **Spatial Processing**
   - Build spatial neighbor graph
   - Apply smoothing kernel
   - Preserve edge spots

4. **Statistical Correction**
   - Fit GAM model
   - Remove spatial trends
   - Preserve biological signal

### Performance Optimizations

- Vectorized operations for score calculation
- KD-tree for efficient neighbor search
- Sparse matrix operations where applicable
- Parallel processing for multiple gene sets
- Caching of spatial graphs

### Error Handling

- Validates gene set input format
- Checks spatial coordinate availability
- Handles empty gene sets
- Reports partial matches
- Graceful degradation without optional features

## Usage Examples

### Example 1: Single Pathway Analysis
```python
# Analyze inflammatory response
result = await analyze_enrichment(
    data_id="tissue_data",
    gene_sets=["TNF", "IL6", "IL1B", "CXCL8", "CCL2"],
    score_keys="Inflammation",
    smoothing=True,
    correct_spatial_covariates=True
)

# Extract results
inflammation_scores = result["enrichment_scores"]["Inflammation"]["corrected"]
hotspots = result["spatial_metrics"]["Inflammation"]["hotspots"]
```

### Example 2: Multiple Pathway Comparison
```python
# Compare multiple pathways
pathways = {
    "Inflammation": ["TNF", "IL6", "IL1B", "CXCL8"],
    "Proliferation": ["MKI67", "TOP2A", "PCNA", "MCM2"],
    "Apoptosis": ["CASP3", "CASP8", "BAX", "BCL2"],
    "Angiogenesis": ["VEGFA", "FLT1", "KDR", "ANGPT1"]
}

result = await analyze_enrichment(
    data_id="tumor_data",
    gene_sets=pathways,
    n_neighbors=10,  # Larger smoothing
    batch_key="patient"  # Correct for patient
)

# Compare pathway activities
for pathway in pathways:
    stats = result["summary_stats"][pathway]
    print(f"{pathway}: mean={stats['mean']:.3f}, Moran's I={stats['spatial_autocorrelation']:.3f}")
```

### Example 3: Cell Type Signatures
```python
# Score cell type presence
cell_markers = {
    "T_cells": ["CD3D", "CD3E", "CD8A", "CD4"],
    "B_cells": ["CD19", "MS4A1", "CD79A"],
    "Macrophages": ["CD68", "CD163", "MRC1"],
    "Fibroblasts": ["COL1A1", "PDGFRA", "THY1"]
}

result = await analyze_enrichment(
    data_id="immune_data",
    gene_sets=cell_markers,
    smoothing=True,
    n_neighbors=8
)
```

### Example 4: Custom Weighted Signatures
```python
# Use gene weights from differential expression
de_genes = {
    "gene": ["GENE1", "GENE2", "GENE3"],
    "logFC": [2.5, 1.8, 1.2],
    "pvalue": [0.001, 0.01, 0.05]
}

# Create weighted signature
genes = de_genes["gene"]
weights = de_genes["logFC"]  # Use fold change as weight

adata = data_store["data_id"]["adata"]
scores = compute_enrichment_scores(
    adata, 
    genes, 
    weights=weights,
    method="weighted_mean"
)
```

### Example 5: Spatial Pattern Analysis
```python
# Focus on spatial metrics
result = await analyze_enrichment(
    data_id="brain_data",
    gene_sets={
        "Layer_markers": layer_specific_genes,
        "Gradient_genes": gradient_genes
    },
    smoothing=False,  # No smoothing for sharp boundaries
    correct_spatial_covariates=True
)

# Analyze spatial patterns
for sig in result["signatures"]:
    metrics = result["spatial_metrics"][sig]
    if metrics["morans_i"] > 0.5 and metrics["morans_p"] < 0.01:
        print(f"{sig} shows significant spatial clustering")
```

## Best Practices

### 1. Gene Set Selection

- **Size**: 10-200 genes optimal
- **Quality**: Use curated databases (MSigDB, GO, KEGG)
- **Specificity**: Avoid overly broad sets
- **Validation**: Check gene presence first

### 2. Parameter Tuning

#### Neighbor Count
- 4-6: Sharp boundaries, local patterns
- 8-12: Moderate smoothing
- 15-20: Strong smoothing, regional patterns

#### Smoothing Decision
- Use for: Noisy data, gradients, large-scale patterns
- Skip for: Sharp boundaries, single-cell resolution

#### Covariate Correction
- Always use for: Technical replicates, batch effects
- Consider skipping: Small datasets, uniform processing

### 3. Quality Control

```python
# Check gene coverage
for sig, stats in result["summary_stats"].items():
    coverage = stats["genes_found"] / stats["n_genes"]
    if coverage < 0.5:
        print(f"Warning: {sig} has low gene coverage ({coverage:.1%})")
```

### 4. Interpretation

- **Spatial autocorrelation > 0.3**: Significant spatial pattern
- **Hotspots**: Regions of high enrichment
- **Gradients**: Progressive changes across tissue
- **Patches**: Discrete enriched regions

## Spatial Considerations

### Edge Effects
- GAM correction handles edge artifacts
- Consider masking edge spots
- Validate with interior regions

### Tissue Geometry
- Adjust neighbors for tissue shape
- Consider anisotropic smoothing
- Account for holes/gaps

### Multi-region Analysis
```python
# Analyze regions separately
regions = adata.obs["region"].unique()
regional_results = {}

for region in regions:
    mask = adata.obs["region"] == region
    regional_data = adata[mask]
    
    scores = compute_enrichment_scores(
        regional_data,
        gene_set,
        method="mean"
    )
    regional_results[region] = scores
```

## Advanced Features

### 1. Temporal Analysis
```python
# Track enrichment over time
timepoints = ["0h", "6h", "24h"]
temporal_enrichment = {}

for tp in timepoints:
    result = await analyze_enrichment(
        data_id=f"data_{tp}",
        gene_sets=pathway_dict
    )
    temporal_enrichment[tp] = result
```

### 2. Multi-scale Analysis
```python
# Different smoothing scales
scales = [4, 8, 16]
multiscale_results = {}

for scale in scales:
    result = await analyze_enrichment(
        data_id="data_1",
        gene_sets=gene_sets,
        n_neighbors=scale
    )
    multiscale_results[scale] = result
```

### 3. Comparative Enrichment
```python
# Compare conditions
conditions = ["control", "treated"]
comparative_results = {}

for condition in conditions:
    result = await analyze_enrichment(
        data_id=f"{condition}_data",
        gene_sets=pathways
    )
    comparative_results[condition] = result

# Calculate differences
for pathway in pathways:
    diff = (comparative_results["treated"]["enrichment_scores"][pathway]["corrected"] -
            comparative_results["control"]["enrichment_scores"][pathway]["corrected"])
```

## Troubleshooting

### Common Issues

1. **"No genes found in dataset"**
   - Check gene naming (symbols vs IDs)
   - Verify species compatibility
   - Try case variations

2. **"Spatial coordinates not found"**
   - Check spatial_key parameter
   - Verify coordinate format
   - Ensure spatial data loaded

3. **"GAM fitting failed"**
   - Reduce spline complexity
   - Check for constant values
   - Remove outliers

4. **"Memory error"**
   - Process gene sets in batches
   - Reduce neighbor count
   - Use sparse operations

### Validation Strategies

1. **Known Patterns**: Test with validated spatial markers
2. **Permutation**: Shuffle genes to test significance
3. **Cross-validation**: Split spots for robustness
4. **Visual Inspection**: Always visualize results

## Performance Considerations

### Computational Complexity
- Score calculation: O(n_spots × n_genes)
- Spatial smoothing: O(n_spots × k)
- GAM correction: O(n_spots × spline_df²)

### Memory Usage
- Scales with n_spots × n_signatures
- Neighbor graph can be large
- Consider chunking for many signatures

### Optimization Tips
- Pre-compute neighbor graphs
- Batch process signatures
- Use appropriate data types
- Cache intermediate results

## Integration with Visualization

```python
# Visualize enrichment scores
vis_params = VisualizationParameters(
    plot_type="spatial",
    feature="Inflammation_enrichment",  # Stored score
    colormap="RdBu_r",
    center=0  # Center colormap at 0
)

# Multi-signature heatmap
vis_params = VisualizationParameters(
    plot_type="enrichment",
    signatures=list(pathways.keys()),
    cluster_spots=True
)
```

## Future Enhancements

1. **Additional Methods**
   - GSEA-style running scores
   - AUCell integration
   - ssGSEA implementation

2. **Statistical Extensions**
   - Permutation testing
   - FDR correction
   - Spatial significance

3. **Advanced Features**
   - Gene set overlap handling
   - Hierarchical signatures
   - Dynamic gene sets

4. **Performance**
   - GPU acceleration
   - Distributed computation
   - Incremental updates