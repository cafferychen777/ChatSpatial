# ChatSpatial - Interactive Spatial Transcriptomics Assistant

ChatSpatial is an interactive spatial transcriptomics data analysis assistant based on the Model Context Protocol (MCP), providing a suite of tools for spatial transcriptomics data processing, visualization, and analysis. It helps researchers analyze spatial transcriptomics data through natural language dialogue.

## Features

- **Data Loading**: Support for various spatial transcriptomics data formats (10x Visium, Slide-seq, MERFISH, seqFISH, etc.)
- **Enhanced Data Preprocessing**: User-controlled filtering, subsampling, normalization, and dimensionality reduction with intelligent defaults for different data types
- **Spatial Visualization**: Spatial distribution of gene expression, visualization of clustering results, etc.
- **Differential Expression Analysis**: Identification of differentially expressed genes between cell populations
- **Cell Type Annotation**: Multiple methods including marker-based, CellAssign, scANVI deep learning annotation
- **Spatial Analysis**: Spatial autocorrelation (Moran's I, Getis-Ord Gi*), neighborhood analysis, spatial trajectories, etc.
- **Spatial Domain Identification**: STAGATE, SpaGCN, and clustering-based methods for identifying spatial domains
- **Advanced Spatial Variable Genes**: GASTON (Graph Attention Spatial Transcriptomics Organizer Network) for learning tissue topology and identifying spatial gene patterns through deep learning
- **Cell Communication Analysis**: LIANA+ integration for fast and comprehensive ligand-receptor interaction analysis with spatial bivariate metrics
- **Advanced Deconvolution**: Complete scvi-tools integration with DestVI, Stereoscope, Cell2location, and traditional methods
- **Standardized Image Processing**: Unified image processing module ensuring all visualization functions return standardized Image objects

## Installation

```bash
# Clone the repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial

# Install basic dependencies
pip install -e .

# Install all optional dependencies
pip install -e .[all]

# Or install specific optional dependencies
pip install -e .[deconvolution]  # Install deconvolution-related dependencies
pip install -e .[spatial_domains]  # Install spatial domain identification dependencies
pip install -e .[spatial_genes]  # Install GASTON and spatial variable genes identification dependencies
pip install -e .[cell_communication]  # Install cell communication analysis dependencies
```

## Usage

### Starting the Server

```bash
# Using stdio transport (default)
chatspatial

# Using SSE transport with a specified port
chatspatial --transport sse --port 8000
```

### Using with Claude Desktop

To use ChatSpatial with Claude Desktop:

1. **Install Claude Desktop**: Download and install Claude Desktop from [Anthropic's website](https://claude.ai/desktop).

2. **Edit Claude Desktop Configuration**:
   - Open the Claude menu on your computer and select "Settings..."
   - Click on "Developer" in the left-hand bar of the Settings pane
   - Click on "Edit Config"
   - This will open the configuration file (`claude_desktop_config.json`) in your text editor

3. **Add ChatSpatial to the Configuration**:
   - Add the following to your configuration file:

   ```json
   {
     "mcpServers": {
       "chatspatial": {
         "command": "/path/to/your/venv/bin/chatspatial",
         "args": [],
         "env": {}
       }
     }
   }
   ```

   - Replace `/path/to/your/venv/bin/chatspatial` with the actual full path to your chatspatial executable
   - Save the file and close the editor

4. **Restart Claude Desktop**:
   - Completely close and restart Claude Desktop
   - After restarting, you should see a hammer icon in the bottom right corner of the input box
   - Click on the hammer icon to see the available tools, including ChatSpatial

5. **Use ChatSpatial in Claude**:
   - Start a new conversation in Claude
   - Click the hammer icon and select "ChatSpatial" from the tools list
   - Claude will connect to your ChatSpatial server
   - You can now interact with Claude to analyze spatial transcriptomics data

#### Example Prompts

Here are some example prompts to get started with ChatSpatial:

- "Load my 10x Visium dataset from `/path/to/data.h5ad`"
- "Preprocess my data with custom filtering: filter genes in less than 10 cells and subsample to 1000 spots"
- "Visualize the spatial expression of gene Cd8a"
- "Perform cell type annotation using marker genes"
- "Identify spatial domains using STAGATE method"
- "Find spatial variable genes using GASTON with GLM-PCA preprocessing"
- "Analyze cell communication using LIANA+"
- "Run spatial trajectory analysis"
- "Deconvolve my spatial data using the NNLS method"
- "Integrate multiple spatial samples"

## Available Tools

The server provides the following tools:

1. `load_data` - Load spatial transcriptomics data
2. `preprocess` - Preprocess data
3. `visualize` - Visualize data
4. `annotate` - Cell type annotation
5. `identify_domains` - Spatial domain identification
6. `find_spatial_genes_gaston` - Advanced spatial variable genes identification using GASTON
7. `analyze_communication` - Cell-cell communication analysis with LIANA+
8. `analyze_spatial_data` - Spatial analysis (Moran's I, Getis-Ord Gi*, neighborhood, co-occurrence, Ripley's K, centrality)
9. `find_markers` - Differential expression analysis
10. `integrate_samples` - Multi-sample integration
11. `analyze_trajectory_pseudotime` - Trajectory pseudotime analysis
12. `analyze_velocity` - RNA velocity analysis
13. `deconvolve` - Spatial transcriptomics deconvolution

### Cell Communication Analysis with LIANA+

ChatSpatial uses LIANA+ (Ligand-receptor Analysis) for fast and comprehensive cell communication analysis:

**Features:**

- **Spatial Bivariate Analysis**: Analyze ligand-receptor pairs using spatial bivariate metrics (cosine, pearson, spearman, jaccard)
- **Global Spatial Metrics**: Moran's I and Lee's L for spatial autocorrelation analysis
- **Multiple LR Databases**: Support for consensus, CellChat, CellPhoneDB, Connectome, and OmniPath databases
- **Fast Performance**: Optimized for interactive usage (1-2 minutes vs 10-30 minutes for other methods)
- **Comprehensive Visualization**: Spatial plots showing communication patterns and significance

**Example Usage:**

```text
"Analyze cell communication using LIANA+ with cosine similarity"
"Find ligand-receptor pairs using Moran's I global metric"
"Visualize spatial communication patterns for VEGFA-KDR interaction"
```

### Spatial Statistics Analysis with Getis-Ord Gi*

ChatSpatial now includes **Getis-Ord Gi*** local spatial autocorrelation analysis for detecting spatial hot spots and cold spots in gene expression:

**Features:**

- **Local Hot/Cold Spot Detection**: Identifies regions with significantly high (hot spots) or low (cold spots) gene expression
- **Multiple Testing Correction**: Supports Bonferroni, FDR-BH, or no correction
- **Flexible Gene Selection**: Analyze specific genes or top highly variable genes
- **Comprehensive Visualization**: Multi-panel spatial plots showing Z-scores with hot/cold spot counts
- **Statistical Summary**: Detailed statistics including total hot/cold spots and significant genes
- **PySAL Implementation**: Uses robust PySAL library for accurate Gi* calculations

**Example Usage:**

```text
"Find spatial hot spots for CCL21 using Getis-Ord Gi* with FDR correction"
"Analyze top 10 highly variable genes for local spatial autocorrelation"
"Detect immune infiltration hot spots using Getis-Ord analysis"
```

### Advanced Spatial Variable Genes with GASTON

ChatSpatial integrates GASTON (Graph Attention Spatial Transcriptomics Organizer Network) for state-of-the-art spatial variable gene identification:

**Features:**

- **Deep Learning Approach**: Neural network-based topology learning for tissue organization
- **Isodepth Mapping**: Learn continuous topographic coordinates that capture tissue structure
- **Spatial Domain Identification**: Automatic identification of spatial domains based on learned topology
- **Flexible Preprocessing**: Support for GLM-PCA and Pearson residuals preprocessing methods
- **Comprehensive Visualization**: Isodepth maps, spatial domains, and model performance metrics
- **Gene Pattern Classification**: Identify genes with continuous gradients vs. discontinuous patterns

**Example Usage:**

```text
"Find spatial variable genes using GASTON with GLM-PCA preprocessing"
"Identify spatial domains using GASTON with 1000 training epochs"
"Visualize tissue topology using GASTON isodepth mapping"
```

### Deep Learning Integration with scvi-tools

ChatSpatial provides complete integration with scvi-tools for state-of-the-art deep learning analysis:

**Cell Type Annotation:**
- **CellAssign**: Probabilistic cell type assignment using marker genes with confidence scores
- **scANVI**: Semi-supervised annotation with reference data transfer learning

**Spatial Deconvolution:**
- **DestVI**: Multi-resolution deconvolution with continuous sub-cell-type variation modeling
- **Stereoscope**: Spatial deconvolution using RNAStereoscope workflow

**Example Usage:**

```text
"Annotate cell types using CellAssign with marker genes"
"Use scANVI to transfer annotations from reference dataset"
"Deconvolve spatial data using DestVI with VampPrior"
"Run Stereoscope deconvolution for cell type proportions"
```

### Enhanced Data Preprocessing

ChatSpatial provides comprehensive, user-controlled data preprocessing with intelligent defaults:

**Features:**

- **User-Controlled Filtering**: Specify exact thresholds for gene and cell filtering
- **Smart Subsampling**: Subsample spots and genes for performance optimization
- **Data Type Adaptation**: Automatic detection and handling of different data types (Visium, MERFISH, etc.)
- **Flexible Parameters**: Full control over normalization, scaling, and dimensionality reduction
- **Quality Control**: Comprehensive QC metrics and filtering efficiency reporting

**Example Usage:**

```text
"Preprocess data with genes expressed in at least 15 cells and cells expressing at least 500 genes"
"Subsample my dataset to 800 spots and 2000 most variable genes for faster analysis"
"Use conservative filtering with 10 minimum cells per gene"
```

## Resources

The server provides the following resources:

- `dataset://{data_id}` - Get dataset information

## Dependencies

### Core Dependencies
- mcp - Model Context Protocol Python SDK
- numpy, pandas - Data processing
- matplotlib - Visualization
- scanpy - Single-cell data analysis
- squidpy - Spatial transcriptomics analysis
- anndata - AnnData data structure
- scikit-learn - Machine learning algorithms

### Optional Dependencies
- **Deep Learning & Advanced Methods**:
  - scvi-tools - Deep learning methods (CellAssign, scANVI, DestVI, Stereoscope)
  - torch, pyro-ppl - PyTorch and Pyro for deep learning models
  - GASTON - Advanced spatial variable genes identification with deep learning
  - glmpca - GLM-PCA preprocessing for GASTON

- **Specialized Analysis**:
  - liana - Cell communication analysis with LIANA+
  - cell2location - Spatial transcriptomics deconvolution
  - Spotiphy - Spatial transcriptomics deconvolution
  - cellrank - Trajectory analysis
  - scvelo - RNA velocity analysis

## License

MIT

## Contributing

Pull Requests and Issues are welcome!
