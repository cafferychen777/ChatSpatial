# Differential Expression Tool Documentation

## Overview

The `differential.py` module provides differential expression (DE) analysis capabilities for spatial transcriptomics data. It enables identification of genes that are differentially expressed between groups of cells, such as different cell types, spatial regions, or experimental conditions.

## Purpose

Differential expression analysis is fundamental for:
- Identifying marker genes for cell types or states
- Understanding spatial gene expression patterns
- Comparing experimental conditions
- Discovering region-specific transcriptional programs
- Validating cell type annotations
- Characterizing disease states

## Core Functionality

The module performs pairwise differential expression testing between two groups of cells, with support for multiple statistical methods and comprehensive result reporting.

## Available Statistical Methods

### 1. Wilcoxon Rank-Sum Test (`wilcoxon`)
**Description**: Non-parametric test comparing gene expression ranks between groups.

**Advantages**:
- Robust to outliers
- No distributional assumptions
- Works well with zero-inflated data
- Default method in most single-cell tools

**Best for**: Most spatial transcriptomics comparisons

### 2. T-test (`t-test`)
**Description**: Parametric test assuming normal distribution.

**Advantages**:
- More powerful when assumptions are met
- Fast computation
- Well-understood statistics

**Best for**: Well-normalized data with sufficient samples

### 3. T-test with Overestimated Variance (`t-test_overestim_var`)
**Description**: T-test with variance overestimation for better FDR control.

**Advantages**:
- More conservative
- Better false discovery control
- Reduces false positives

**Best for**: Exploratory analysis, small sample sizes

### 4. Logistic Regression (`logreg`)
**Description**: Models gene expression as binary presence/absence.

**Advantages**:
- Handles sparse data well
- Can include covariates
- Provides effect sizes

**Best for**: Highly sparse data, marker gene identification

## Input Parameters

### Function Signature
```python
async def differential_expression(
    data_id: str,                    # Dataset identifier
    data_store: Dict[str, Any],      # Data storage
    group_key: str,                  # Column in adata.obs for grouping
    group1: str,                     # First group name
    group2: str,                     # Second group name
    n_top_genes: int = 50,          # Number of top genes to return
    method: str = "wilcoxon",       # Statistical method
    context: Context = None         # MCP context
) -> DifferentialExpressionResult
```

### Parameters Explained

- **group_key**: The column in `adata.obs` containing group labels (e.g., "cell_type", "leiden", "spatial_domain")
- **group1**: Name of the first group to compare
- **group2**: Name of the second group to compare (or "rest" for one-vs-rest)
- **n_top_genes**: Limits output to top differentially expressed genes
- **method**: Statistical test to use

## Output Format

### DifferentialExpressionResult
```python
class DifferentialExpressionResult:
    # Top genes information
    top_genes: List[str]          # Gene names
    log_fold_changes: List[float] # Log2 fold changes
    p_values: List[float]         # Raw p-values
    adjusted_p_values: List[float] # FDR-corrected p-values
    
    # Comparison metadata
    comparison: str               # "group1_vs_group2"
    method: str                   # Statistical method used
    n_cells_group1: int          # Cells in group 1
    n_cells_group2: int          # Cells in group 2
    
    # Full results
    full_results_key: str        # Key in adata.uns for all genes
```

### Result Interpretation

- **log_fold_changes**: 
  - Positive: Higher in group1
  - Negative: Higher in group2
  - Magnitude: Effect size

- **adjusted_p_values**: 
  - < 0.05: Statistically significant
  - < 0.01: Highly significant
  - Use these instead of raw p-values

## Implementation Details

### Preprocessing Steps

1. **Group Validation**: Ensures both groups exist and have sufficient cells
2. **Gene Filtering**: Removes genes with too few expressing cells
3. **Data Preparation**: Extracts expression matrix for comparison

### Statistical Pipeline

1. **Test Execution**: Runs selected statistical test
2. **Multiple Testing Correction**: Benjamini-Hochberg FDR correction
3. **Ranking**: Orders genes by significance and effect size
4. **Filtering**: Selects top genes based on adjusted p-value and fold change

### Error Handling

- Validates group existence
- Checks minimum cell requirements
- Handles missing data gracefully
- Provides informative error messages

## Usage Examples

### Example 1: Finding Cluster Markers
```python
# Identify markers for cluster 3
result = await find_markers(
    data_id="spatial_data",
    group_key="leiden",
    group1="3",
    group2="rest",  # Compare against all other clusters
    n_top_genes=100,
    method="wilcoxon"
)

# Top marker genes for cluster 3
print(f"Top markers: {result.top_genes[:10]}")
```

### Example 2: Comparing Spatial Regions
```python
# Compare tumor vs normal regions
result = await find_markers(
    data_id="cancer_data",
    group_key="tissue_region",
    group1="tumor",
    group2="normal",
    n_top_genes=200,
    method="t-test"
)

# Tumor-enriched genes
tumor_genes = [g for g, lfc in zip(result.top_genes, result.log_fold_changes) if lfc > 1]
```

### Example 3: Cell Type Markers with Conservative Testing
```python
# Use overestimated variance for robust markers
result = await find_markers(
    data_id="data_1",
    group_key="cell_type",
    group1="T_cells",
    group2="B_cells",
    n_top_genes=50,
    method="t-test_overestim_var"
)
```

### Example 4: Comparing Treatment Conditions
```python
# Treatment effect analysis
result = await find_markers(
    data_id="treatment_data",
    group_key="condition",
    group1="treated",
    group2="control",
    n_top_genes=500,  # Get more genes for pathway analysis
    method="wilcoxon"
)

# Filter for strong effects
significant_genes = [
    g for g, p, lfc in zip(result.top_genes, result.adjusted_p_values, result.log_fold_changes)
    if p < 0.01 and abs(lfc) > 1
]
```

### Example 5: Fine-grained Spatial Analysis
```python
# Compare specific spatial coordinates
# First, create groups based on spatial location
adata = data_store["data_1"]["adata"]
adata.obs["spatial_group"] = ["region_A" if x < 5000 else "region_B" 
                               for x in adata.obsm["spatial"][:, 0]]

# Then run DE
result = await find_markers(
    data_id="data_1",
    group_key="spatial_group",
    group1="region_A",
    group2="region_B",
    method="logreg"  # Good for spatial patterns
)
```

## Best Practices

### 1. Group Selection
- Ensure adequate sample size (>10 cells per group)
- Check group balance (avoid 10 vs 10000 comparisons)
- Validate biological relevance

### 2. Method Selection
- Default to Wilcoxon for general use
- Use t-test for well-powered comparisons
- Use logreg for marker discovery
- Use t-test_overestim_var for conservative analysis

### 3. Result Validation
- Check both p-value and fold change
- Visualize top genes spatially
- Validate with known markers
- Consider biological plausibility

### 4. Multiple Comparisons
```python
# Run multiple comparisons systematically
cell_types = ["T_cells", "B_cells", "Myeloid", "Fibroblasts"]
all_markers = {}

for ct in cell_types:
    result = await find_markers(
        data_id="data_1",
        group_key="cell_type",
        group1=ct,
        group2="rest",
        n_top_genes=50
    )
    all_markers[ct] = result.top_genes[:10]
```

## Spatial Considerations

### Spatial Autocorrelation
Genes with spatial patterns may show false positive DE results. Consider:
- Checking spatial distribution of DE genes
- Using spatial statistics alongside DE
- Accounting for spatial structure in analysis

### Region-based Analysis
```python
# Define regions based on domains
result = await find_markers(
    data_id="data_1",
    group_key="spatial_domain",  # From spatial domain analysis
    group1="Domain_1",
    group2="Domain_2"
)
```

### Distance-based Groups
```python
# Compare center vs periphery
adata = data_store["data_1"]["adata"]
center = np.mean(adata.obsm["spatial"], axis=0)
distances = np.linalg.norm(adata.obsm["spatial"] - center, axis=1)
adata.obs["location"] = ["center" if d < np.median(distances) else "periphery" 
                          for d in distances]
```

## Performance Considerations

### Computational Complexity
- Wilcoxon: O(n × m × log(m)) for n genes, m cells
- T-test: O(n × m)
- Logreg: O(n × m × iterations)

### Memory Usage
- Scales with number of genes and cells
- Sparse matrix handling for efficiency
- Results storage proportional to n_top_genes

### Optimization Tips
- Pre-filter genes with low expression
- Use chunking for very large datasets
- Consider subsampling for exploration
- Cache results for repeated queries

## Integration with Other Modules

### Visualization
```python
# Visualize top DE genes
top_gene = result.top_genes[0]
vis_params = VisualizationParameters(
    plot_type="spatial",
    feature=top_gene,
    groups_to_show=[group1, group2]  # Highlight compared groups
)
```

### Enrichment Analysis
```python
# Use DE genes for enrichment
enriched_genes = result.top_genes[:100]
enrichment_result = await analyze_enrichment(
    data_id="data_1",
    gene_sets={"DE_genes": enriched_genes}
)
```

### Validation
```python
# Validate with spatial analysis
for gene in result.top_genes[:5]:
    spatial_result = await analyze_spatial_data(
        data_id="data_1",
        params=SpatialAnalysisParameters(
            analysis_type="morans_i",
            morans_i_gene=gene
        )
    )
```

## Common Issues and Solutions

### 1. "No cells in group"
- Check group names match exactly
- Verify group_key exists in adata.obs
- Print unique values: `adata.obs[group_key].unique()`

### 2. "Too few cells for comparison"
- Increase group size by merging similar groups
- Use less stringent filtering
- Consider different grouping strategy

### 3. "No significant genes"
- Check data normalization
- Verify groups are truly different
- Increase n_top_genes
- Try different statistical method

### 4. "Memory error"
- Reduce number of genes tested
- Use sparse matrix operations
- Process in chunks

## Advanced Features

### Custom Gene Lists
```python
# Test specific genes only
genes_of_interest = ["CD3E", "CD4", "CD8A", "FOXP3"]
adata_subset = adata[:, genes_of_interest]
# Run DE on subset
```

### Covariate Adjustment
```python
# For future implementation
# Account for batch effects, spatial location, etc.
```

### Paired Comparisons
```python
# For matched samples (future feature)
# Useful for before/after treatment
```

## Quality Control

### Result Validation Checklist
1. ✓ Check sample sizes for both groups
2. ✓ Verify fold change direction makes sense
3. ✓ Examine p-value distribution
4. ✓ Visualize top genes spatially
5. ✓ Compare with known biology
6. ✓ Test method sensitivity

### Diagnostic Plots
- Volcano plots (fold change vs p-value)
- MA plots (mean vs fold change)
- P-value histograms
- Spatial expression patterns

## Future Enhancements

1. **Additional Methods**
   - MAST (for zero-inflation)
   - DESeq2 adaptation
   - Spatial-aware DE tests

2. **Advanced Features**
   - Covariate adjustment
   - Paired sample support
   - Time-series DE

3. **Performance**
   - GPU acceleration
   - Parallel processing
   - Approximate methods

4. **Integration**
   - Direct pathway analysis
   - Network enrichment
   - Multi-condition comparisons