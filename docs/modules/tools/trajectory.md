# Trajectory Tool Documentation

## Overview

The `trajectory.py` module provides comprehensive trajectory analysis capabilities for spatial transcriptomics data, including RNA velocity computation and pseudotime inference. It integrates state-of-the-art methods to reveal cellular dynamics, differentiation trajectories, and fate decisions in spatial context.

## Purpose

Trajectory analysis in spatial transcriptomics enables:
- Understanding cellular differentiation processes
- Mapping developmental trajectories in tissue
- Identifying cell fate decisions
- Revealing dynamic processes in static snapshots
- Connecting spatial organization with temporal dynamics

## Components

### 1. RNA Velocity Analysis

RNA velocity leverages the ratio of spliced to unspliced RNA to infer future cell states and directionality of differentiation.

### 2. Trajectory Inference

Pseudotime ordering and fate probability estimation to understand cellular progression and lineage relationships.

## RNA Velocity Methods

### scVelo Modes

#### 1. Stochastic Mode (`stochastic`)
- **Description**: Original RNA velocity formulation
- **Assumptions**: Steady-state kinetics
- **Best for**: Quick analysis, stable systems

#### 2. Deterministic Mode (`deterministic`)
- **Description**: Deterministic model without noise
- **Assumptions**: Deterministic dynamics
- **Best for**: Clear trajectories, less noisy data

#### 3. Dynamical Mode (`dynamical`)
- **Description**: Learns gene-specific kinetics
- **Key Features**:
  - Infers splicing kinetics
  - Handles transient states
  - Most accurate but slower
- **Best for**: Complex dynamics, high-quality data

### SIRV Method (`sirv`)

**Description**: Spatial Inference of RNA Velocity - spatially-aware velocity computation.

**Key Features**:
- Incorporates spatial information
- Smooths velocities across neighbors
- Reduces noise in spatial context

**Mathematical Basis**:
```
v_spatial = αv_local + (1-α)Σ(w_ij × v_j) / Σw_ij
```
Where w_ij are spatial weights

## Trajectory Inference Methods

### 1. CellRank (`cellrank`)

**Description**: Combines RNA velocity with gene expression similarity for robust trajectory inference.

**Key Features**:
- Multi-kernel approach
- Markov chain analysis
- Automatic fate detection
- Handles uncertainty

**Mathematical Foundation**:
- Velocity kernel: Direction from RNA velocity
- Connectivity kernel: Expression similarity
- Combined kernel: K = α·K_velocity + (1-α)·K_connectivity
- GPCCA for fate identification

### 2. Palantir (`palantir`)

**Description**: Diffusion map-based trajectory inference with entropy-based pseudotime.

**Key Features**:
- No velocity requirement
- Branching trajectory support
- Fate probability computation
- Start cell specification

**Mathematical Foundation**:
- Diffusion maps for manifold learning
- Markov chain for fate probabilities
- Entropy for pseudotime: S(t) = -Σp_i·log(p_i)

### 3. DPT Fallback (`dpt`)

**Description**: Diffusion pseudotime as robust fallback method.

**Key Features**:
- Always works
- Simple and fast
- Single trajectory
- No fate probabilities

## Input Parameters

### RNAVelocityParameters
```python
class RNAVelocityParameters:
    # Core parameters
    mode: str = "dynamical"  # Velocity mode
    min_shared_counts: int = 30  # Min counts
    n_pcs: int = 30  # PCA components
    n_neighbors: int = 30  # KNN neighbors
    
    # Advanced parameters
    filter_genes: bool = True  # Gene filtering
    filter_genes_dispersion: Optional[Tuple] = None
    
    # Dynamical mode specific
    dynamical_params: Dict = {
        "fit_basal_transcription": True,
        "steady_state_prior": True,
        "n_top_genes": 2000
    }
    
    # SIRV specific
    sirv_params: Dict = {
        "spatial_weight": 0.5,  # Spatial influence
        "normalize_spatial": True
    }
```

### TrajectoryParameters
```python
class TrajectoryParameters:
    # Method selection
    method: str = "cellrank"  # cellrank/palantir
    
    # Start/end points
    start_cell: Optional[str] = None  # Start cell ID
    end_cells: Optional[List[str]] = None  # Terminal states
    
    # Common parameters
    n_neighbors: int = 30
    n_pcs: int = 30
    
    # CellRank specific
    cellrank_params: Dict = {
        "weight_velocity": 0.2,  # Velocity weight
        "softmax_scale": 4,  # Kernel scaling
        "n_states": None,  # Auto-detect
        "method": "GPCCA",  # Terminal state method
        "use_spatial": True  # Spatial kernel
    }
    
    # Palantir specific
    palantir_params: Dict = {
        "n_diffusion_components": 10,
        "knn": 30,
        "n_jobs": -1,
        "use_early_cell": True
    }
```

## Output Format

### RNAVelocityResult
```python
class RNAVelocityResult:
    # Core results
    velocity_computed: bool
    velocity_key: str = "velocity"  # In .layers
    
    # Quality metrics
    velocity_confidence: Optional[float]
    velocity_coherence: Optional[float]
    
    # Method info
    method: str
    parameters_used: Dict[str, Any]
    
    # Additional outputs
    velocity_genes: Optional[List[str]]
    kinetic_params: Optional[Dict]  # For dynamical
```

### TrajectoryResult
```python
class TrajectoryResult:
    # Pseudotime
    pseudotime_key: str = "pseudotime"  # In .obs
    pseudotime_computed: bool
    
    # Fate probabilities
    fate_probabilities_key: Optional[str]  # In .obsm
    terminal_states: Optional[List[str]]
    n_fates: Optional[int]
    
    # Lineage information
    lineages: Optional[List[str]]
    branch_points: Optional[List[str]]
    
    # Method metadata
    method_used: str
    parameters_used: Dict[str, Any]
    
    # Quality metrics
    fate_correlation: Optional[float]
```

## Implementation Details

### RNA Velocity Workflow

1. **Data Preparation**
   ```python
   # Check for spliced/unspliced layers
   # Filter genes and cells
   # Normalize if needed
   ```

2. **Velocity Computation**
   ```python
   # Mode-specific velocity calculation
   # Spatial smoothing for SIRV
   # Project to embedding
   ```

3. **Quality Control**
   ```python
   # Compute confidence scores
   # Check velocity coherence
   # Validate against known markers
   ```

### Trajectory Inference Workflow

1. **Preprocessing**
   ```python
   # PCA/dimensionality reduction
   # Neighbor graph construction
   # Start cell identification
   ```

2. **Method Execution**
   ```python
   # CellRank: Kernel combination → GPCCA
   # Palantir: Diffusion maps → Markov chain
   # DPT: Diffusion distance computation
   ```

3. **Post-processing**
   ```python
   # Extract pseudotime
   # Compute fate probabilities
   # Identify branch points
   ```

### Fallback Mechanism

The module implements a robust fallback strategy:
1. Try primary method (CellRank/Palantir)
2. If fails, try alternate method
3. If both fail, use DPT (always works)
4. Log warnings and continue

## Usage Examples

### Example 1: Standard RNA Velocity
```python
# Compute RNA velocity with dynamical model
velocity_result = await analyze_velocity_data(
    data_id="data_1",
    params=RNAVelocityParameters(
        mode="dynamical",
        n_neighbors=30,
        dynamical_params={
            "n_top_genes": 3000
        }
    )
)
```

### Example 2: Spatial RNA Velocity
```python
# SIRV for spatial smoothing
velocity_result = await analyze_velocity_data(
    data_id="spatial_data",
    params=RNAVelocityParameters(
        mode="sirv",
        sirv_params={
            "spatial_weight": 0.7,  # High spatial influence
            "normalize_spatial": True
        }
    )
)
```

### Example 3: CellRank Trajectory
```python
# Full trajectory analysis with CellRank
traj_result = await analyze_trajectory_data(
    data_id="data_1",
    params=TrajectoryParameters(
        method="cellrank",
        cellrank_params={
            "weight_velocity": 0.3,
            "n_states": 3,  # Expected endpoints
            "use_spatial": True
        }
    )
)
```

### Example 4: Palantir with Start Cell
```python
# Developmental trajectory from known start
traj_result = await analyze_trajectory_data(
    data_id="dev_data",
    params=TrajectoryParameters(
        method="palantir",
        start_cell="CELL_12345",  # Stem cell
        palantir_params={
            "n_diffusion_components": 15,
            "knn": 50
        }
    )
)
```

### Example 5: Combined Analysis
```python
# First velocity, then trajectory
velocity_result = await analyze_velocity_data(
    data_id="data_1",
    params=RNAVelocityParameters(mode="dynamical")
)

traj_result = await analyze_trajectory_data(
    data_id="data_1",
    params=TrajectoryParameters(
        method="cellrank",
        cellrank_params={"weight_velocity": 0.5}
    )
)
```

## Best Practices

### 1. Data Preparation

- **Quality Control**: Ensure sufficient unspliced counts
- **Gene Selection**: Use dynamical genes for velocity
- **Normalization**: Match methods between analyses
- **Batch Effects**: Correct before trajectory analysis

### 2. Method Selection

#### RNA Velocity
- **Stochastic**: Quick exploration
- **Deterministic**: Stable systems
- **Dynamical**: Publication quality
- **SIRV**: Spatial datasets

#### Trajectory Inference
- **CellRank**: When velocity is available
- **Palantir**: Complex branching without velocity
- **DPT**: Simple trajectories, fallback

### 3. Parameter Tuning

- **n_neighbors**: 20-50 based on data density
- **velocity weight**: 0.1-0.5 in CellRank
- **spatial weight**: 0.3-0.8 for SIRV
- **n_states**: Based on biological knowledge

### 4. Validation

- Check velocity coherence
- Validate against known markers
- Compare methods
- Visualize spatially

## Biological Scenarios

### Embryonic Development
```python
params = TrajectoryParameters(
    method="palantir",
    start_cell=early_embryo_cell,
    palantir_params={
        "n_diffusion_components": 20  # Complex landscape
    }
)
```

### Tumor Progression
```python
params = TrajectoryParameters(
    method="cellrank",
    cellrank_params={
        "n_states": 2,  # Normal vs malignant
        "weight_velocity": 0.4
    }
)
```

### Immune Response
```python
params = RNAVelocityParameters(
    mode="dynamical",
    dynamical_params={
        "fit_basal_transcription": False  # Activated cells
    }
)
```

### Tissue Regeneration
```python
params = TrajectoryParameters(
    method="cellrank",
    start_cell=stem_cell_id,
    cellrank_params={
        "use_spatial": True,  # Spatial gradients matter
        "weight_velocity": 0.6
    }
)
```

## Troubleshooting

### Common Issues

1. **"No velocity could be computed"**
   - Check spliced/unspliced layers
   - Increase min_shared_counts
   - Try different mode

2. **"No terminal states found"**
   - Adjust n_states parameter
   - Check data quality
   - Try different method

3. **"Memory error"**
   - Reduce n_neighbors
   - Subsample cells
   - Use fewer genes

4. **"Convergence failed"**
   - Increase iterations
   - Adjust parameters
   - Check for batch effects

### Debug Mode
```python
# Enable verbose output
params = RNAVelocityParameters(
    mode="dynamical",
    dynamical_params={
        "verbose": True,
        "plot_results": True
    }
)
```

## Performance Considerations

### Time Complexity
- Stochastic velocity: O(n × g)
- Dynamical velocity: O(n × g × iterations)
- CellRank: O(n² × iterations)
- Palantir: O(n² × components)

### Memory Usage
- Velocity: ~2-3x expression matrix
- Trajectory: ~n² for transition matrix
- Use sparse matrices when possible

### Optimization Tips
- Pre-filter genes aggressively
- Use GPU for CellRank
- Parallelize Palantir
- Cache intermediate results

## Integration with Other Modules

### Visualization
```python
# Visualize velocity
vis_params = VisualizationParameters(
    plot_type="velocity",
    velocity_key="velocity",
    basis="spatial"
)

# Visualize trajectory
vis_params = VisualizationParameters(
    plot_type="trajectory",
    pseudotime_key=result.pseudotime_key,
    fate_key=result.fate_probabilities_key
)
```

### Spatial Analysis
```python
# Spatial patterns of pseudotime
spatial_params = SpatialAnalysisParameters(
    analysis_type="morans_i",
    feature=result.pseudotime_key
)
```

### Differential Expression
```python
# Genes changing along trajectory
de_params = {
    "groupby": result.pseudotime_key,
    "method": "logreg"
}
```

## Advanced Features

### 1. Multi-Modal Integration
```python
# Combine RNA velocity with protein
params = RNAVelocityParameters(
    mode="dynamical",
    use_protein_velocity=True,  # If available
    protein_names=["CD4", "CD8"]
)
```

### 2. Spatial Constraints
```python
# Force trajectories to respect anatomy
params = TrajectoryParameters(
    method="cellrank",
    cellrank_params={
        "spatial_constraint": "soft",
        "anatomy_mask": tissue_regions
    }
)
```

### 3. Time Series
```python
# Multiple timepoints
params = TrajectoryParameters(
    method="cellrank",
    timepoint_key="day",
    cellrank_params={
        "time_forward_only": True
    }
)
```

## Future Enhancements

1. **VeloVI Integration**: Deep learning velocity
2. **MultiVelo**: Multi-modal velocity
3. **Spatial Trajectories**: True spatial paths
4. **Uncertainty Quantification**: Bayesian approaches
5. **Interactive Selection**: Web-based start/end cell picking
6. **Trajectory Alignment**: Cross-sample comparison