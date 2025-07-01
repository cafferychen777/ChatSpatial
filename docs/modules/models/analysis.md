# Analysis Models Documentation

## Overview

The `analysis.py` module defines all result models returned by analysis tools in the ChatSpatial MCP system. These Pydantic models provide structured, type-safe representations of analysis outputs, ensuring consistent interfaces for result handling and interpretation.

## Purpose

Result models in the MCP architecture serve to:
- **Standardize Outputs**: Uniform result format across all tools
- **Type Safety**: Ensure result fields are properly typed
- **Documentation**: Self-documenting result structures
- **Data Access**: Provide keys to access stored data in AnnData
- **Metadata**: Include analysis parameters and quality metrics
- **Serialization**: Enable JSON serialization for MCP protocol

## Core Result Models

### PreprocessingResult
```python
class PreprocessingResult(BaseModel):
    """Result of data preprocessing"""
    data_id: str  # Dataset identifier
    n_cells: int  # Number of cells after filtering
    n_genes: int  # Number of genes after filtering
    n_hvgs: int   # Number of highly variable genes
    clusters: int  # Number of clusters identified
    qc_metrics: Dict[str, Any]  # Quality control statistics
```

**Stored in AnnData**:
- `adata.obs`: QC metrics (n_genes_by_counts, total_counts)
- `adata.var`: HVG flags (highly_variable)
- `adata.obsm`: PCA (X_pca), UMAP (X_umap)
- `adata.uns`: Preprocessing parameters

**Quality Metrics Example**:
```python
qc_metrics = {
    "n_cells_before": 5000,
    "n_cells_after": 4500,
    "n_genes_before": 20000,
    "n_genes_after": 15000,
    "mean_genes_per_cell": 3200,
    "mean_counts_per_cell": 12000,
    "pct_cells_filtered": 10.0,
    "pct_genes_filtered": 25.0
}
```

### AnnotationResult
```python
class AnnotationResult(BaseModel):
    """Result of cell type annotation"""
    cell_types: List[str]  # Unique cell type labels
    cell_type_counts: Dict[str, int]  # Count per type
    confidence_scores: Optional[List[float]]  # Per-cell confidence
    method_used: str  # Annotation method
    parameters_used: Dict[str, Any]  # Method parameters
    annotation_key: str = "cell_type"  # Key in adata.obs
    
    # Method-specific results
    tangram_mapping_score: Optional[float] = None
    scanvi_latent_key: Optional[str] = None
    cellassign_probability_key: Optional[str] = None
```

**Stored in AnnData**:
- `adata.obs[annotation_key]`: Cell type labels
- `adata.obs[f"{annotation_key}_confidence"]`: Confidence scores
- `adata.obsm`: Method-specific embeddings (e.g., scANVI latent)
- `adata.uns`: Detailed results and parameters

**Interpretation**:
- `cell_types`: Ordered list of discovered cell types
- `confidence_scores`: Higher values indicate more reliable annotations
- Method-specific scores provide quality assessment

### SpatialAnalysisResult
```python
class SpatialAnalysisResult(BaseModel):
    """Result of spatial pattern analysis"""
    analysis_type: str  # Type of analysis performed
    statistics: Dict[str, Any]  # Analysis-specific statistics
    significant_features: Optional[List[str]]  # Significant genes/types
    spatial_key: str  # Spatial coordinates used
    parameters_used: Dict[str, Any]  # Analysis parameters
```

**Analysis-Specific Statistics**:

**Neighborhood Enrichment**:
```python
statistics = {
    "enrichment_matrix": np.ndarray,  # Cell type enrichment
    "zscore_matrix": np.ndarray,      # Statistical significance
    "analysis_key_in_adata": "neighborhood_enrichment"
}
```

**Moran's I**:
```python
statistics = {
    "morans_i": 0.42,  # Global autocorrelation
    "p_value": 0.001,  # Significance
    "gene": "CD3E",    # Analyzed gene
    "analysis_key_in_adata": "morans_i_CD3E"
}
```

**Getis-Ord Gi***:
```python
statistics = {
    "n_hotspots": 150,     # Hot spot count
    "n_coldspots": 75,     # Cold spot count
    "genes_analyzed": [...],  # Gene list
    "analysis_key_in_adata": "getis_ord_stats"
}
```

### DifferentialExpressionResult
```python
class DifferentialExpressionResult(BaseModel):
    """Result of differential expression analysis"""
    top_genes: List[str]  # Top DE genes
    log_fold_changes: List[float]  # Log2 fold changes
    p_values: List[float]  # Raw p-values
    adjusted_p_values: List[float]  # FDR-corrected
    comparison: str  # e.g., "T_cells_vs_B_cells"
    method: str  # Statistical method used
    n_cells_group1: int  # Cells in first group
    n_cells_group2: int  # Cells in second group
    full_results_key: str = "de_results"  # Key in adata.uns
```

**Stored in AnnData**:
- `adata.uns[full_results_key]`: Complete DE results DataFrame
- Columns: gene, logFC, pval, pval_adj, pct_1, pct_2

**Interpretation**:
- Positive logFC: Higher in group1
- Negative logFC: Higher in group2
- adjusted_p_values < 0.05: Statistically significant

### DeconvolutionResult
```python
class DeconvolutionResult(BaseModel):
    """Result of spatial deconvolution"""
    proportions: pd.DataFrame  # Cell type proportions per spot
    cell_types: List[str]  # Deconvolved cell types
    method_used: str  # Deconvolution method
    parameters_used: Dict[str, Any]  # Method parameters
    
    # Quality metrics
    statistics: Dict[str, Any] = {
        "mean_entropy": float,  # Diversity measure
        "dominant_cell_type_pct": float,  # Dominance
        "n_cell_types": int,  # Types found
        "sparsity": float  # Proportion of zeros
    }
    
    # Storage
    proportions_key: str = "deconvolution_proportions"
    
    # Visualization hints
    visualization_params: Optional[Dict[str, Any]] = None
    
    # Method-specific
    confidence_scores: Optional[pd.DataFrame] = None
    convergence_info: Optional[Dict[str, Any]] = None
```

**Stored in AnnData**:
- `adata.obsm[proportions_key]`: Proportion matrix
- `adata.uns["deconvolution"]`: Method info and stats

**Quality Indicators**:
- High entropy: Mixed composition
- Low entropy: Dominated by single type
- Convergence info: Model reliability

### SpatialDomainResult
```python
class SpatialDomainResult(BaseModel):
    """Result of spatial domain identification"""
    n_domains: int  # Number of domains found
    domain_key: str = "spatial_domain"  # Key in adata.obs
    refined_domain_key: Optional[str] = None  # Refined labels
    
    # Embeddings
    embedding_key: Optional[str] = None  # Key in adata.obsm
    
    # Quality metrics
    statistics: Dict[str, Any] = {
        "silhouette_score": float,  # Clustering quality
        "calinski_harabasz": float,  # Separation
        "davies_bouldin": float,    # Compactness
        "spatial_coherence": float,  # Spatial continuity
        "domain_sizes": Dict[str, int]  # Size per domain
    }
    
    # Method info
    method_used: str  # STAGATE, SpaGCN, etc.
    parameters_used: Dict[str, Any]
    convergence_info: Optional[Dict[str, Any]] = None
```

**Stored in AnnData**:
- `adata.obs[domain_key]`: Domain labels
- `adata.obs[refined_domain_key]`: Refined labels (if applicable)
- `adata.obsm[embedding_key]`: Method embeddings

**Quality Assessment**:
- Silhouette > 0.5: Good separation
- High spatial coherence: Domains are contiguous
- Balanced domain sizes: No extreme outliers

### TrajectoryResult
```python
class TrajectoryResult(BaseModel):
    """Result of trajectory analysis"""
    # Pseudotime
    pseudotime_key: str = "pseudotime"
    pseudotime_computed: bool
    
    # Fate probabilities
    fate_probabilities_key: Optional[str] = None
    terminal_states: Optional[List[str]] = None
    n_fates: Optional[int] = None
    
    # Lineage information
    lineages: Optional[List[str]] = None
    branch_points: Optional[List[str]] = None
    
    # Method metadata
    method_used: str  # cellrank, palantir, dpt
    parameters_used: Dict[str, Any]
    
    # Quality metrics
    fate_correlation: Optional[float] = None
```

**Stored in AnnData**:
- `adata.obs[pseudotime_key]`: Pseudotime values
- `adata.obsm[fate_probabilities_key]`: Fate probability matrix
- `adata.uns["trajectory"]`: Detailed trajectory info

### RNAVelocityResult
```python
class RNAVelocityResult(BaseModel):
    """Result of RNA velocity analysis"""
    velocity_computed: bool
    velocity_key: str = "velocity"  # Key in adata.layers
    
    # Quality metrics
    velocity_confidence: Optional[float] = None
    velocity_coherence: Optional[float] = None
    
    # Method info
    method: str  # stochastic, dynamical, etc.
    parameters_used: Dict[str, Any]
    
    # Additional outputs
    velocity_genes: Optional[List[str]] = None
    kinetic_params: Optional[Dict[str, Any]] = None
```

**Stored in AnnData**:
- `adata.layers[velocity_key]`: Velocity matrix
- `adata.var["velocity_genes"]`: Genes used
- `adata.uns["velocity_params"]`: Parameters

### IntegrationResult
```python
class IntegrationResult(BaseModel):
    """Result of multi-sample integration"""
    integrated_data_id: str  # New dataset ID
    n_samples_integrated: int  # Number of samples
    n_cells_total: int  # Total cells
    integration_key: str = "batch"  # Batch key
    
    # Method info
    method_used: str  # harmony, bbknn, etc.
    parameters_used: Dict[str, Any]
    
    # Quality metrics
    integration_metrics: Dict[str, float] = {
        "silhouette_score": float,  # Batch mixing
        "ari_score": float,         # Type preservation
        "batch_entropy": float,     # Mixing entropy
        "kbet_rejection_rate": float  # Batch effect
    }
    
    # Sample information
    sample_cell_counts: Dict[str, int]
    common_genes: int
    
    # Embeddings
    corrected_pca_key: str = "X_pca_integrated"
    umap_key: str = "X_umap_integrated"
```

### CellCommunicationResult
```python
class CellCommunicationResult(BaseModel):
    """Result of cell-cell communication analysis"""
    n_interactions: int  # Total tested
    n_significant_pairs: int  # Significant LR pairs
    
    # Top interactions
    top_lr_pairs: List[str]  # Top ranked pairs
    top_source_target: List[Tuple[str, str]]  # Cell pairs
    
    # Method info
    method_used: str = "liana"
    analysis_mode: str  # "global" or "local"
    
    # Result storage
    global_results_key: Optional[str] = None
    local_results_key: Optional[str] = None
    
    # Spatial analysis
    local_analysis_performed: bool = False
    spatial_interactions: Optional[Dict[str, Any]] = None
    
    # Filtering info
    n_cell_types: int
    filtered_interactions: int
```

**Stored in AnnData**:
- `adata.uns[global_results_key]`: LR interaction DataFrame
- `adata.obsm[local_results_key]`: Spatial LR patterns

### SpatialVariableGenesResult
```python
class SpatialVariableGenesResult(BaseModel):
    """Result of spatial variable gene identification"""
    # Gene counts
    n_continuous_genes: int  # Gradient genes
    n_discontinuous_genes: int  # Jump genes
    continuous_genes: List[str]
    discontinuous_genes: List[str]
    
    # Spatial information
    n_spatial_domains: int
    isodepth_key: str = "gaston_isodepth"
    spatial_domains_key: str = "gaston_domains"
    
    # Gene analysis
    gene_r2_values: Dict[str, float]  # Fit quality
    gene_max_lfc: Dict[str, float]  # Max fold change
    gene_patterns: Dict[str, str]  # Pattern type
    
    # Model performance
    final_loss: float
    training_history: Dict[str, List[float]]
    model_performance: Dict[str, float]
    
    # Storage
    model_key: str = "gaston_model"
    predictions_key: str = "gaston_predictions"
```

## Data Storage Patterns

### Standard AnnData Locations

**Cell/Spot Annotations** (`adata.obs`):
- Cell types
- Clusters
- Pseudotime
- Domain labels
- Quality metrics

**Gene Annotations** (`adata.var`):
- Highly variable genes
- Velocity genes
- Spatial patterns

**Embeddings** (`adata.obsm`):
- PCA: X_pca
- UMAP: X_umap
- Spatial: spatial
- Method-specific: X_method

**Pairwise Matrices** (`adata.obsp`):
- Spatial connectivity
- Cell-cell distances

**Unstructured Data** (`adata.uns`):
- Analysis results
- Parameters
- Metadata
- Complex objects

### Naming Conventions

**Keys**:
- Use underscores: `cell_type`, not `cellType`
- Prefix with method: `stagate_domains`
- Include analysis type: `morans_i_CD3E`

**Suffixes**:
- `_key`: References to data locations
- `_scores`: Numerical scores
- `_params`: Parameters used

## Result Model Relationships

### Tool-Result Mapping

| Tool | Primary Result | Secondary Results |
|------|---------------|-------------------|
| preprocess_data | PreprocessingResult | - |
| annotate_cells | AnnotationResult | - |
| analyze_spatial_data | SpatialAnalysisResult | - |
| find_markers | DifferentialExpressionResult | - |
| deconvolve_data | DeconvolutionResult | - |
| identify_spatial_domains | SpatialDomainResult | - |
| analyze_velocity_data | RNAVelocityResult | - |
| analyze_trajectory_data | TrajectoryResult | RNAVelocityResult |
| integrate_samples | IntegrationResult | - |
| analyze_cell_communication | CellCommunicationResult | - |
| find_spatial_genes | SpatialVariableGenesResult | - |

### Result Dependencies

Some analyses depend on previous results:
```
PreprocessingResult
    ↓
AnnotationResult → SpatialAnalysisResult
    ↓               ↓
DeconvolutionResult CellCommunicationResult
```

## Result Interpretation

### Quality Indicators

**Good Quality**:
- High confidence scores (>0.7)
- Low p-values (<0.05)
- High silhouette scores (>0.5)
- Convergence achieved

**Poor Quality**:
- Low confidence (<0.5)
- High p-values (>0.1)
- Negative silhouette scores
- Convergence warnings

### Key Metrics

**Spatial Coherence**: How well spatial patterns align
- 0.8-1.0: Very coherent
- 0.5-0.8: Moderately coherent
- <0.5: Poor coherence

**Entropy**: Diversity measure
- High: Mixed composition
- Low: Homogeneous

**Fold Change**: Expression differences
- >2: Strong difference
- 1-2: Moderate
- <1: Weak

## Best Practices

### 1. Result Validation
```python
# Check key fields exist
if result.confidence_scores and min(result.confidence_scores) < 0.5:
    print("Warning: Low confidence annotations detected")

# Verify data storage
if result.annotation_key not in adata.obs:
    print("Error: Annotations not stored properly")
```

### 2. Accessing Stored Data
```python
# Use result keys to access data
adata = data_store[result.data_id]["adata"]
annotations = adata.obs[result.annotation_key]
embeddings = adata.obsm[result.embedding_key]
```

### 3. Result Chaining
```python
# Use results from one analysis in another
preprocess_result = await preprocess_data(...)
annotation_result = await annotate_cells(...)

# Access preprocessed data for further analysis
spatial_result = await analyze_spatial_data(
    data_id=preprocess_result.data_id,
    params=SpatialAnalysisParameters(
        cluster_key=annotation_result.annotation_key
    )
)
```

### 4. Error Handling
```python
# Check for optional fields
if result.convergence_info:
    if not result.convergence_info["converged"]:
        print("Warning: Analysis did not converge")

# Handle missing data gracefully
fate_probs = getattr(result, "fate_probabilities_key", None)
if fate_probs:
    probs = adata.obsm[fate_probs]
```

## Examples

### Example 1: Processing Preprocessing Results
```python
# Run preprocessing
preprocess_result = await preprocess_data(
    data_id="sample1",
    params=AnalysisParameters()
)

# Examine results
print(f"Filtered to {preprocess_result.n_cells} cells")
print(f"Found {preprocess_result.n_hvgs} variable genes")
print(f"Identified {preprocess_result.clusters} clusters")

# Check quality
qc = preprocess_result.qc_metrics
if qc["pct_cells_filtered"] > 50:
    print("Warning: Filtered >50% of cells")
```

### Example 2: Working with Spatial Domains
```python
# Identify domains
domain_result = await identify_spatial_domains(
    data_id="spatial_data",
    params=SpatialDomainParameters()
)

# Access domains
adata = data_store["spatial_data"]["adata"]
domains = adata.obs[domain_result.domain_key]

# Check quality
stats = domain_result.statistics
print(f"Silhouette score: {stats['silhouette_score']:.3f}")
print(f"Spatial coherence: {stats['spatial_coherence']:.3f}")

# Use refined domains if available
if domain_result.refined_domain_key:
    refined = adata.obs[domain_result.refined_domain_key]
    print(f"Refined domains available: {refined.nunique()} domains")
```

### Example 3: Analyzing Communication Results
```python
# Run communication analysis
comm_result = await analyze_cell_communication(
    data_id="immune_data",
    params=CellCommunicationParameters()
)

# Summary
print(f"Tested {comm_result.n_interactions} interactions")
print(f"Found {comm_result.n_significant_pairs} significant")

# Top interactions
for lr_pair in comm_result.top_lr_pairs[:5]:
    print(f"Top LR: {lr_pair}")

# Access full results
if comm_result.global_results_key:
    full_results = adata.uns[comm_result.global_results_key]
    # Filter for specific cell types
    tcell_interactions = full_results[
        full_results["source"] == "T_cells"
    ]
```

### Example 4: Complex Result Handling
```python
# GASTON spatial genes
spatial_genes_result = await find_spatial_genes(
    data_id="tissue_data",
    params=SpatialVariableGenesParameters()
)

# Analyze patterns
print(f"Continuous gradients: {spatial_genes_result.n_continuous_genes}")
print(f"Discontinuous patterns: {spatial_genes_result.n_discontinuous_genes}")

# Get top gradient genes
gradient_genes = sorted(
    spatial_genes_result.continuous_genes,
    key=lambda g: spatial_genes_result.gene_r2_values[g],
    reverse=True
)[:10]

# Check model performance
perf = spatial_genes_result.model_performance
if perf["r2"] < 0.5:
    print("Warning: Low model fit")
```

### Example 5: Result Visualization
```python
# Many results include visualization hints
if hasattr(result, "visualization_params") and result.visualization_params:
    # Use suggested parameters
    vis_result = await visualize_data(
        data_id=result.data_id,
        params=VisualizationParameters(**result.visualization_params)
    )
else:
    # Create custom visualization
    vis_result = await visualize_data(
        data_id=result.data_id,
        params=VisualizationParameters(
            plot_type="spatial",
            feature=result.annotation_key
        )
    )
```

## Serialization and Storage

### JSON Serialization
All result models support JSON serialization:
```python
# Serialize to JSON
result_json = result.json()

# Deserialize from JSON
result = AnnotationResult.parse_raw(result_json)
```

### Persistent Storage
Results can be saved with the dataset:
```python
# Store in AnnData
adata.uns["analysis_results"] = {
    "preprocessing": preprocess_result.dict(),
    "annotation": annotation_result.dict(),
    "spatial": spatial_result.dict()
}

# Save to disk
adata.write("analyzed_data.h5ad")
```

## Performance Considerations

### Memory Efficiency
- Results store references (keys) not data
- Large matrices remain in AnnData
- Minimal duplication

### Access Patterns
```python
# Efficient: Access via keys
domains = adata.obs[result.domain_key]

# Inefficient: Storing large data in result
# result.large_matrix = huge_array  # Don't do this
```

## Future Enhancements

1. **Result Validation**
   - Automatic quality checks
   - Standardized metrics
   - Warning thresholds

2. **Result Comparison**
   - Compare across methods
   - Ensemble results
   - Concordance metrics

3. **Visualization Integration**
   - Auto-generate plots
   - Interactive reports
   - Result dashboards

4. **Caching**
   - Result caching
   - Incremental updates
   - Checkpointing