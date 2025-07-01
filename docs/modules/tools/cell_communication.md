# Cell Communication Tool Documentation

## Overview

The `cell_communication.py` module analyzes cell-cell communication patterns in spatial transcriptomics data through ligand-receptor (LR) interactions. It leverages the LIANA+ framework to identify significant communication events and their spatial organization.

## Purpose

Cell communication analysis is essential for:
- Understanding tissue organization and function
- Identifying signaling networks in microenvironments
- Discovering disease-relevant interactions
- Mapping spatial communication patterns
- Revealing cell type crosstalk
- Characterizing developmental signaling

## Core Functionality

The module provides two main approaches:
1. **Cluster-based analysis**: Traditional cell type communication
2. **Spatial analysis**: Location-aware interaction patterns

## LIANA+ Framework

### Key Features
- Consensus approach combining multiple methods
- Comprehensive LR databases
- Statistical significance testing
- Spatial autocorrelation metrics
- Multi-species support

### Analysis Methods

#### 1. Rank Aggregate (Cluster-based)
Combines multiple scoring methods:
- **CellPhoneDB**: Statistical enrichment
- **Connectome**: Expression strength
- **log2FC**: Differential expression
- **NATMI**: Specificity scores
- **SingleCellSignalR**: Interaction probability
- **CellChat**: Probability and communication strength

#### 2. Spatial Bivariate (Location-based)
Analyzes spatial co-expression:
- **Local Moran's I**: Spatial autocorrelation
- **Lee's L**: Bivariate spatial association
- **Getis-Ord Gi***: Hot spot analysis

## Input Parameters

### CellCommunicationParameters
```python
class CellCommunicationParameters:
    # Core parameters
    method: str = "liana"  # Currently only LIANA+
    groupby: str = "cell_type"  # Cell grouping
    use_raw: bool = True  # Use raw counts
    
    # Resource selection
    resource_name: str = "consensus"  # LR database
    species: Optional[str] = None  # Auto-detected
    
    # Filtering parameters
    min_cells: int = 10  # Min cells per group
    min_expr: float = 0.1  # Min expression threshold
    
    # Analysis parameters
    n_perms: int = 1000  # Permutations
    seed: int = 42  # Random seed
    
    # Method-specific parameters
    return_all_lrs: bool = True  # All interactions
    aggregate_method: str = "rra"  # Rank aggregate
    
    # Spatial parameters
    spatial_key: str = "spatial"  # Coordinates
    connectivity_key: Optional[str] = None
    bandwidth: float = 150  # Spatial kernel
    cutoff: float = 0.01  # Distance cutoff
    
    # Performance parameters
    use_gpu: bool = False  # GPU acceleration
    verbose: bool = True  # Progress updates
    inplace: bool = True  # Modify adata
```

## Ligand-Receptor Databases

### Available Resources
1. **consensus** (default): Aggregated from multiple databases
2. **cellphonedb**: CellPhoneDB interactions
3. **cellchat**: CellChat database
4. **connectome**: Connectome database
5. **mouseconsensus**: Mouse-specific consensus

### Database Content
- Protein-protein interactions
- Secreted signaling molecules
- Cell surface receptors
- ECM interactions
- Multi-subunit complexes

## Output Format

### CellCommunicationResult
```python
class CellCommunicationResult:
    # Core results
    n_interactions: int  # Total interactions tested
    n_significant_pairs: int  # Significant LR pairs
    
    # Top interactions
    top_lr_pairs: List[str]  # Top ranked pairs
    top_source_target: List[Tuple[str, str]]  # Cell type pairs
    
    # Method info
    method_used: str = "liana"
    analysis_mode: str  # "global" or "local"
    
    # Result storage
    global_results_key: Optional[str]  # For cluster analysis
    local_results_key: Optional[str]  # For spatial analysis
    
    # Spatial analysis
    local_analysis_performed: bool
    spatial_interactions: Optional[Dict]  # Spatial patterns
    
    # Filtering info
    n_cell_types: int
    filtered_interactions: int
```

### Result Structure

#### Global Analysis (Cluster-based)
```python
adata.uns['liana_res'] = {
    'liana_res': DataFrame with columns:
        - source: Source cell type
        - target: Target cell type  
        - ligand_complex: Ligand name
        - receptor_complex: Receptor name
        - magnitude_rank: Combined score
        - specificity_rank: Specificity score
        - p_values: Individual method p-values
}
```

#### Local Analysis (Spatial)
```python
adata.obsm['local_liana'] = {
    'ligand-receptor': Array with spatial scores
}

adata.uns['local_liana'] = {
    'local_stats': Spatial statistics
    'categories': Spatial categories per spot
}
```

## Implementation Details

### Preprocessing Steps
1. **Species Detection**: Infers from gene names (human vs mouse)
2. **Data Validation**: Checks expression data format
3. **Group Filtering**: Removes small groups
4. **Gene Filtering**: Ensures LR genes are present

### Analysis Workflow

#### Cluster-based Analysis
1. Calculate mean expression per cell type
2. Score each LR pair with multiple methods
3. Aggregate scores using RRA
4. Filter by significance thresholds
5. Rank interactions

#### Spatial Analysis
1. Build spatial connectivity matrix
2. Calculate local expression products
3. Compute spatial autocorrelation
4. Identify spatial patterns
5. Categorize interaction zones

## Usage Examples

### Example 1: Basic Cell Type Communication
```python
# Standard analysis
result = await analyze_cell_communication(
    data_id="data_1",
    params=CellCommunicationParameters(
        groupby="cell_type",
        min_cells=20,
        return_all_lrs=False  # Only significant
    )
)

print(f"Found {result.n_significant_pairs} significant interactions")
print(f"Top pair: {result.top_lr_pairs[0]}")
```

### Example 2: Tumor Microenvironment
```python
# Immune-tumor interactions
result = await analyze_cell_communication(
    data_id="tumor_data",
    params=CellCommunicationParameters(
        groupby="cell_type_detailed",
        resource_name="cellphonedb",
        min_expr=0.2,  # Higher threshold
        n_perms=5000  # More permutations
    )
)

# Extract T cell - Tumor interactions
tcell_tumor = [
    (s, t, l, r) for s, t, l, r in result.top_interactions
    if "T_cell" in s and "Tumor" in t
]
```

### Example 3: Spatial Communication Patterns
```python
# Spatial LR analysis
result = await analyze_cell_communication(
    data_id="spatial_data",
    params=CellCommunicationParameters(
        groupby="cell_type",
        spatial_key="spatial",
        bandwidth=200,  # Larger neighborhoods
        local_analysis=True
    )
)

# Significant spatial patterns
if result.local_analysis_performed:
    print(f"Spatial patterns stored in {result.local_results_key}")
```

### Example 4: Developmental Gradients
```python
# Morphogen signaling
result = await analyze_cell_communication(
    data_id="embryo_data",
    params=CellCommunicationParameters(
        groupby="developmental_stage",
        resource_name="consensus",
        focus_ligands=["WNT3A", "SHH", "BMP4"],  # Key morphogens
        spatial_analysis=True
    )
)
```

### Example 5: Neural Circuits
```python
# Synaptic communication
result = await analyze_cell_communication(
    data_id="brain_data",
    params=CellCommunicationParameters(
        groupby="neuron_subtype",
        resource_name="consensus",
        focus_receptors=["GRIA1", "GRIN1", "GABRA1"],
        min_cells=5  # Rare subtypes
    )
)
```

## Best Practices

### 1. Data Preparation
- Ensure accurate cell type annotations
- Use raw counts when possible
- Check species compatibility
- Validate marker expression

### 2. Parameter Selection

#### Minimum Cells
- Standard: 10-20 cells
- Rare populations: 5-10 cells
- Major types: 50+ cells

#### Expression Threshold
- Default: 0.1 (10% of cells)
- Stringent: 0.2-0.3
- Permissive: 0.05

#### Permutations
- Quick analysis: 100-1000
- Publication: 5000-10000
- Spatial: 1000 minimum

### 3. Result Interpretation

#### Significance Levels
- magnitude_rank < 0.01: Very strong
- magnitude_rank < 0.05: Strong
- specificity_rank < 0.05: Cell type specific

#### Spatial Categories
- High-High: Co-localized high expression
- Low-Low: Co-localized low expression
- High-Low: Spatial segregation
- Not Significant: No pattern

### 4. Validation
- Check known interactions
- Validate with IF/FISH
- Compare across samples
- Test spatial proximity

## Spatial Analysis Details

### Connectivity Methods
1. **KNN**: K-nearest neighbors
2. **Radius**: Fixed distance
3. **Delaunay**: Triangulation

### Spatial Metrics

#### Local Moran's I
```
I_i = (x_i - x̄) × Σ_j w_ij(x_j - x̄)
```
Measures local spatial autocorrelation

#### Lee's L
```
L_ij = Σ_k w_ik × (x_k - x̄) × (y_k - ȳ)
```
Bivariate spatial association

### Interpretation
- Positive I: Clustering
- Negative I: Dispersion
- High L: Co-expression
- Low L: Segregation

## Troubleshooting

### Common Issues

1. **"No cell type annotations"**
   - Run annotation first
   - Check groupby column exists
   - Verify column values

2. **"No interactions found"**
   - Lower min_expr threshold
   - Check species setting
   - Verify gene names

3. **"Memory error"**
   - Reduce return_all_lrs
   - Filter interactions
   - Use subset of cells

4. **"LIANA not installed"**
   - Install: `pip install liana`
   - Check dependencies
   - Use conda environment

### Debug Mode
```python
params = CellCommunicationParameters(
    verbose=True,
    debug=True,
    save_intermediates=True
)
```

## Performance Optimization

### Computational Complexity
- Cluster: O(n_types² × n_lr_pairs)
- Spatial: O(n_spots × n_neighbors × n_lr_pairs)

### Memory Usage
- Scales with interactions × cell types
- Spatial needs connectivity matrix
- GPU can reduce memory pressure

### Speed Tips
- Pre-filter LR pairs
- Use fewer permutations initially
- Subset to regions of interest
- Parallelize when possible

## Visualization Integration

### Cluster Communication
```python
# Network visualization
vis_params = VisualizationParameters(
    plot_type="cell_communication",
    interaction_type="network",
    top_n=20
)
```

### Spatial Patterns
```python
# LR pair spatial expression
vis_params = VisualizationParameters(
    plot_type="lr_pairs",
    ligand="CXCL12",
    receptor="CXCR4"
)
```

### Dot Plots
```python
# Communication matrix
vis_params = VisualizationParameters(
    plot_type="cell_communication",
    interaction_type="dotplot"
)
```

## Advanced Features

### 1. Focused Analysis
```python
# Specific pathways
params = CellCommunicationParameters(
    focus_ligands=["TNF", "IL6", "IL1B"],
    focus_receptors=["TNFRSF1A", "IL6R", "IL1R1"],
    pathway_name="inflammation"
)
```

### 2. Multi-sample Comparison
```python
# Compare conditions
conditions = ["control", "treated"]
results = {}

for condition in conditions:
    result = await analyze_cell_communication(
        data_id=f"data_{condition}",
        params=params
    )
    results[condition] = result
```

### 3. Temporal Analysis
```python
# Time series communication
timepoints = ["0h", "6h", "24h"]
temporal_results = []

for tp in timepoints:
    result = await analyze_cell_communication(
        data_id=f"data_{tp}",
        params=params
    )
    temporal_results.append(result)
```

## Biological Applications

### 1. Cancer Research
- Tumor-immune crosstalk
- Metastatic niche formation
- Therapy resistance mechanisms
- Angiogenesis signaling

### 2. Neuroscience
- Synaptic networks
- Glial-neuron interaction
- Neuroinflammation
- Development guidance

### 3. Immunology
- Cytokine networks
- Checkpoint interactions
- Inflammatory cascades
- Tissue resident programs

### 4. Development
- Morphogen gradients
- Cell fate signaling
- Organogenesis
- Stem cell niches

## Quality Control

### Validation Checklist
1. ✓ Cell type annotations accurate
2. ✓ Sufficient cells per type
3. ✓ Expression thresholds appropriate
4. ✓ Known interactions detected
5. ✓ Spatial patterns make sense
6. ✓ Results reproducible

### Positive Controls
- Known interactions (e.g., PD1-PDL1)
- Housekeeping communication
- Validated pathways

## Future Enhancements

1. **Additional Methods**
   - CellChat integration
   - NicheNet downstream
   - Commot implementation

2. **Advanced Features**
   - Multi-modal integration
   - Trajectory-based communication
   - Cross-species analysis

3. **Spatial Extensions**
   - 3D communication
   - Dynamic signaling
   - Gradient modeling

4. **Machine Learning**
   - Interaction prediction
   - Pattern classification
   - Network inference