# Server Module Documentation

## Overview

The `server.py` module is the core implementation of the ChatSpatial MCP (Model Context Protocol) server. It provides the main FastMCP server instance that handles all spatial transcriptomics analysis operations through a set of well-defined tools.

## Architecture

### Core Components

1. **FastMCP Server Instance**
   - Uses FastMCP framework for MCP protocol implementation
   - Supports both stdio and SSE (Server-Sent Events) transport protocols
   - Instance name: "ChatSpatial"

2. **Data Store**
   - Global dictionary `data_store` that maintains loaded datasets
   - Each dataset is assigned a unique ID (e.g., "data_1", "data_2")
   - Stores AnnData objects and metadata

3. **Visualization Resources**
   - Global storage for generated visualization images
   - Implements MCP resource protocol for image access
   - Resources are saved with unique URIs (e.g., "visualization://viz_data_1_spatial_1234567890")

## Main Functions

### Data Management

#### `validate_dataset(data_id: str) -> None`
Validates that a dataset exists in the data store. Raises ValueError if dataset not found.

#### `load_data()`
- **Purpose**: Load spatial transcriptomics data from various formats
- **Parameters**:
  - `data_path`: Path to data file or directory
  - `data_type`: Format type (auto, 10x_visium, slide_seq, merfish, seqfish, other, h5ad)
  - `name`: Optional dataset name
- **Returns**: `SpatialDataset` object with dataset metadata
- **Supported Formats**: H5AD, 10x Visium, Slide-seq, MERFISH, seqFISH

### Preprocessing

#### `preprocess_data()`
- **Purpose**: Preprocess spatial transcriptomics data
- **Parameters**:
  - `data_id`: Dataset identifier
  - `params`: `AnalysisParameters` object
- **Returns**: `PreprocessingResult` object
- **Features**:
  - Multiple normalization methods (log, sct, none, scvi)
  - Quality control filtering
  - Dimensionality reduction
  - Optional scVI preprocessing for advanced denoising

### Visualization

#### `visualize_data()`
- **Purpose**: Generate various visualizations for spatial data
- **Parameters**:
  - `data_id`: Dataset identifier
  - `params`: `VisualizationParameters` object
- **Returns**: Formatted message with visualization resource URI
- **Plot Types**:
  - spatial: Spatial expression plots
  - heatmap: Gene expression heatmaps
  - violin: Violin plots
  - umap: UMAP embeddings
  - spatial_domains: Domain visualization
  - cell_communication: LR pair networks
  - deconvolution: Cell type proportions
  - trajectory: Pseudotime visualization
  - spatial_analysis: Analysis results
  - multi_gene: Multiple genes
  - lr_pairs: Ligand-receptor pairs
  - gene_correlation: Gene correlations
  - gaston_isodepth: GASTON isodepth coordinates
  - gaston_domains: GASTON spatial domains
  - gaston_genes: GASTON gene patterns

### Analysis Tools

#### `annotate_cells()`
- **Purpose**: Annotate cell types in spatial data
- **Parameters**:
  - `data_id`: Dataset identifier
  - `params`: `AnnotationParameters` object
- **Returns**: `AnnotationResult` object
- **Methods**:
  - marker_genes: Known marker gene annotation
  - tangram: Single-cell to spatial mapping
  - scanvi: Semi-supervised scANVI
  - cellassign: Probabilistic assignment
  - correlation, supervised, popv, gptcelltype, scrgcl

#### `analyze_spatial_data()`
- **Purpose**: Perform spatial pattern analysis
- **Parameters**:
  - `data_id`: Dataset identifier
  - `params`: `SpatialAnalysisParameters` object
- **Returns**: `SpatialAnalysisResult` object
- **Analysis Types**:
  - neighborhood: Spatial neighborhoods
  - morans_i: Spatial autocorrelation
  - getis_ord: Hot/cold spot analysis
  - ripley: Spatial point patterns
  - co_occurrence: Cell type co-occurrence
  - interaction_matrix: Cell interactions
  - spatial_connectivity: Connectivity analysis

#### `find_markers()`
- **Purpose**: Find differentially expressed genes between groups
- **Parameters**:
  - `data_id`: Dataset identifier
  - `group_key`: Observation column for grouping
  - `group1`, `group2`: Groups to compare
  - `n_top_genes`: Number of top genes
  - `method`: Statistical method (wilcoxon, t-test, etc.)
- **Returns**: `DifferentialExpressionResult` object

#### `analyze_velocity_data()`
- **Purpose**: Analyze RNA velocity in spatial context
- **Parameters**:
  - `data_id`: Dataset identifier
  - `params`: `RNAVelocityParameters` object
- **Returns**: `RNAVelocityResult` object
- **Modes**: stochastic, deterministic, dynamical

#### `analyze_trajectory_data()`
- **Purpose**: Infer cell trajectories and pseudotime
- **Parameters**:
  - `data_id`: Dataset identifier
  - `params`: `TrajectoryParameters` object
- **Returns**: `TrajectoryResult` object
- **Methods**:
  - cellrank: CellRank for trajectory inference
  - palantir: Palantir for branching trajectories
  - velovi: VeloVI for velocity-based trajectories

#### `integrate_samples()`
- **Purpose**: Integrate multiple spatial samples
- **Parameters**:
  - `data_ids`: List of dataset identifiers
  - `params`: `IntegrationParameters` object
- **Returns**: `IntegrationResult` object
- **Methods**:
  - harmony: Batch effect correction
  - bbknn: Batch-aware k-NN
  - scanorama: Panoramic integration
  - mnn: Mutual nearest neighbors
  - scvi: Probabilistic integration
  - multivi: Multi-modal integration
  - totalvi: Protein+RNA integration

#### `deconvolve_data()`
- **Purpose**: Estimate cell type proportions in spatial spots
- **Parameters**:
  - `data_id`: Dataset identifier
  - `params`: `DeconvolutionParameters` object
- **Returns**: `DeconvolutionResult` object
- **Methods**:
  - cell2location: Bayesian deconvolution
  - spotiphy: Graph-based deconvolution
  - rctd: Robust cell type decomposition
  - destvi: DestVI from scvi-tools
  - stereoscope: Stereoscope from scvi-tools

#### `identify_spatial_domains()`
- **Purpose**: Identify spatial domains/regions
- **Parameters**:
  - `data_id`: Dataset identifier
  - `params`: `SpatialDomainParameters` object
- **Returns**: `SpatialDomainResult` object
- **Methods**:
  - stagate: Graph attention network
  - spagcn: Graph convolutional network
  - bayesspace: Bayesian spatial clustering
  - banksy: Neighborhood-based clustering

#### `analyze_cell_communication()`
- **Purpose**: Analyze cell-cell communication via ligand-receptor pairs
- **Parameters**:
  - `data_id`: Dataset identifier
  - `params`: `CellCommunicationParameters` object
- **Returns**: `CellCommunicationResult` object
- **Methods**:
  - liana: LIANA+ framework
  - cellphonedb: CellPhoneDB
  - nichenet: NicheNet
  - commot: COMMOT

#### `analyze_enrichment()`
- **Purpose**: Spatially-aware gene set enrichment using EnrichMap
- **Parameters**:
  - `data_id`: Dataset identifier
  - `gene_sets`: Gene lists or dictionary of gene sets
  - `score_keys`: Names for gene signatures
  - `spatial_key`: Spatial coordinates key
  - `n_neighbors`: Neighbors for smoothing
  - `smoothing`: Enable spatial smoothing
  - `correct_spatial_covariates`: GAM correction
  - `batch_key`: Batch normalization key
- **Returns**: Enrichment scores and statistics

#### `find_spatial_genes()`
- **Purpose**: Identify spatial variable genes using GASTON
- **Parameters**:
  - `data_id`: Dataset identifier
  - `params`: `SpatialVariableGenesParameters` object
- **Returns**: `SpatialVariableGenesResult` object
- **Features**:
  - Deep learning-based spatial pattern detection
  - Isodepth coordinate learning
  - Continuous/discontinuous gene classification
  - Spatial domain identification

## Resource Management

### Visualization Resources

The server implements MCP resource protocol for visualization outputs:

1. **Resource Creation**: When visualizations are generated, they are saved as PNG files and registered as MCP resources
2. **Resource URIs**: Follow pattern `visualization://viz_{data_id}_{plot_type}_{timestamp}`
3. **Resource Storage**: Files saved in `visualization_resources/` directory
4. **Resource Access**: Available through MCP resource handlers

### Resource Handlers

- `list_visualization_resources()`: Lists all available visualization resources
- `read_visualization_resource(uri)`: Reads and returns visualization data as base64-encoded blob

## Error Handling

The server implements comprehensive error handling:

1. **Tool Error Handler**: `@mcp_tool_error_handler()` decorator catches and formats tool errors
2. **Pydantic Error Handler**: `@mcp_pydantic_error_handler()` handles parameter validation errors
3. **Manual Parameter Validation**: `@manual_parameter_validation()` provides fallback validation
4. **Formatted Error Messages**: All errors formatted according to MCP error protocol

## Transport Modes

### stdio Mode (Default)
- Standard input/output communication
- Used for command-line interface
- Single session per process

### SSE Mode
- Server-Sent Events for HTTP transport
- Supports multiple concurrent sessions
- Web browser compatibility

## Usage

### Command Line
```bash
# Default stdio mode
python -m chatspatial.server

# SSE mode for HTTP transport
python -m chatspatial.server --transport sse
```

### Integration
The server can be integrated with:
- Claude Desktop (via MCP configuration)
- Custom MCP clients
- Web applications (using SSE mode)

## Data Flow

1. **Data Loading**: User loads spatial data → Stored in `data_store` with unique ID
2. **Preprocessing**: Applied to stored data → Updates AnnData object in place
3. **Analysis**: Various analysis tools operate on preprocessed data → Results stored in AnnData
4. **Visualization**: Generate plots from analysis results → Saved as MCP resources
5. **Resource Access**: Clients retrieve visualizations through resource URIs

## Best Practices

1. **Session Management**: Each dataset gets unique ID for session isolation
2. **Memory Efficiency**: Large datasets handled through efficient AnnData operations
3. **Error Recovery**: Comprehensive error handling prevents server crashes
4. **Resource Cleanup**: Visualization resources persisted for session duration
5. **Parameter Validation**: Multiple validation layers ensure robust operation

## Extensions

The server is designed for extensibility:
- New tools can be added as decorated functions
- Additional data formats supported through data loader
- Custom analysis methods integrated through tool modules
- Visualization types extended in visualization module