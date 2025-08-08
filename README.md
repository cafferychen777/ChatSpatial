# ChatSpatial MCP Server

**Interactive Spatial Transcriptomics Analysis via Model Context Protocol**

ChatSpatial is a **Model Context Protocol (MCP) server** that provides AI assistants with comprehensive spatial transcriptomics analysis capabilities. It enables natural language interaction with spatial data analysis through 32+ standardized MCP tools.

🔗 **MCP Protocol**: Seamless integration with Claude Desktop, Continue, and other MCP-compatible AI tools  
🧬 **Spatial Focus**: Specialized for spatial transcriptomics (10x Visium, Slide-seq, MERFISH, etc.)  
⚡ **Production Ready**: All core methods tested and validated for real-world usage

## MCP Tools & Capabilities

### 🔄 **Data Management** 
- `load_data` - Multi-format spatial data loading (10x Visium, Slide-seq, MERFISH, seqFISH)
- `preprocess_data` - QC, normalization, dimensionality reduction with intelligent defaults

### 📊 **Visualization**
- `visualize_data` - 15+ plot types (spatial, UMAP, violin, heatmap, trajectory, communication)
- Returns standard MCP `Image` objects for seamless AI assistant display

### 🔬 **Cell Analysis** 
- `annotate_cells` - Multiple annotation methods (marker-based, Tangram, scANVI, CellAssign)
- `find_markers` - Differential expression analysis between cell populations

### 🧬 **Spatial Analysis**
- `find_spatial_genes` - **GASTON** (deep learning), **SpatialDE**, **SPARK** methods ✅ All working
- `identify_spatial_domains` - SpaGCN, STAGATE, BANKSY spatial clustering
- `analyze_spatial_data` - Moran's I, Geary's C, spatial autocorrelation

### 💬 **Cell Communication**
- `analyze_cell_communication` - **LIANA+** ligand-receptor analysis with spatial bivariate metrics

### 🔗 **Integration & Trajectory**
- `integrate_samples` - Multi-sample integration (Harmony, scVI, BBKNN)
- `analyze_trajectory_data` - Pseudotime inference (Palantir, CellRank, DPT)
- `analyze_velocity_data` - RNA velocity analysis

### 🧮 **Advanced Methods**
- `deconvolve_data` - Cell type deconvolution (DestVI, Cell2location, RCTD, SPOTlight)
- `analyze_enrichment` - Pathway analysis (GSEA, ORA, Enrichr) with spatial smoothing

## Installation

### Python Version Requirements

- **Core ChatSpatial**: Python 3.8+
- **Recommended**: Python 3.10 or 3.11 for best compatibility

### Setting up a Dedicated Environment

Due to ChatSpatial's complex dependencies (PyTorch, scvi-tools, spatial analysis packages), we **strongly recommend** creating a dedicated Python environment:

```bash
# Create a dedicated environment for ChatSpatial
conda create -n chatspatial_env python=3.10
conda activate chatspatial_env

# Or using virtualenv
python3.10 -m venv chatspatial_env
source chatspatial_env/bin/activate  # On Windows: chatspatial_env\Scripts\activate
```

### Install ChatSpatial

```bash
# Clone the repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial

# Install basic dependencies
pip install -e .

# Install all optional dependencies (recommended for full functionality)
pip install -e .[all]

# Or install specific optional dependencies
pip install -e .[enrichmap]  # Install EnrichMap enrichment analysis dependencies
pip install -e .[deconvolution]  # Install deconvolution-related dependencies
pip install -e .[spatial_domains]  # Install spatial domain identification dependencies
pip install -e .[spatial_genes]  # Install GASTON and spatial variable genes identification dependencies
pip install -e .[cell_communication]  # Install cell communication analysis dependencies
```

### Verify Installation

```bash
# Test the installation
chatspatial --help

# Find the exact path to your chatspatial executable (needed for MCP configuration)
which chatspatial
# Or if using conda/virtual environment:
which python
```

### GASTON Installation (Optional)

GASTON is a deep learning method for spatial gene analysis. To use GASTON features:

#### Method 1: Install from PyPI (Recommended)

```bash
pip install gaston-spatial
```

#### Method 2: Install from Source

```bash
# Clone the official repository
git clone https://github.com/Arashz/GASTON.git
cd GASTON
pip install -e .
```

#### Development Workflow (Recommended for Contributors)

For development and comparison with the original implementation:

```bash
# Clone original repository for comparison (not committed to version control)
git clone https://github.com/Arashz/GASTON.git gaston_dev
cd gaston_dev
pip install -e .

# Run comparison tests to ensure compatibility
python -c "import gaston; print('GASTON installed successfully')"
```

**Note**: GASTON is optional but recommended for advanced spatial gene analysis. ChatSpatial also supports SpatialDE and SPARK methods for comprehensive spatial variable gene identification.

## Usage

### Starting the Server

```bash
# Using stdio transport (default)
chatspatial

# Using SSE transport with a specified port
chatspatial --transport sse --port 8000
```

## MCP Server Configuration

ChatSpatial implements the **Model Context Protocol (MCP)** standard for seamless AI assistant integration.

### Quick MCP Setup

1. **Find your Python path**:
```bash
# Activate your ChatSpatial environment
conda activate chatspatial_env

# Get the Python executable path 
which python
# Example output: /opt/anaconda3/envs/chatspatial_env/bin/python
```

2. **Configure your MCP client** (Claude Desktop, Continue, etc.):
```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/your/chatspatial_env/bin/python", 
      "args": ["-m", "chatspatial"],
      "env": {}
    }
  }
}
```

### MCP Protocol Details

**Server Information:**
- **Protocol**: MCP v2024-11-05
- **Tools**: 32 spatial transcriptomics analysis tools
- **Resources**: Automatic dataset and result management via `spatial://` URIs
- **Transport**: stdio (standard input/output)

**Data Flow:**
```
AI Assistant → MCP tools/call → ChatSpatial Tool → Analysis Result → MCP ToolResult
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
         "command": "/path/to/your/chatspatial_env/bin/python",
         "args": ["-m", "chatspatial"],
         "env": {}
       }
     }
   }
   ```

   - Replace `/path/to/your/chatspatial_env/bin/python` with the actual path to your ChatSpatial environment's Python executable
   - Use `python -m chatspatial` instead of direct executable for better compatibility
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

### Using with Cherry Studio

To use ChatSpatial with Cherry Studio (recommended for tasks requiring longer processing time):

1. **Install Cherry Studio**: Download and install Cherry Studio from the [official website](https://cherrystudio.ai/).

2. **Configure MCP Server**:
   - Open Cherry Studio Settings
   - Navigate to the MCP Servers section
   - Click "Add Server" and configure as follows:

   **Basic Configuration:**

   ```json
   {
     "name": "ChatSpatial",
     "command": "/path/to/your/chatspatial_env/bin/python",
     "args": ["-m", "chatspatial"],
     "env": {}
   }
   ```

   **Configuration Details:**
   - **Name**: `ChatSpatial` (or any name you prefer)
   - **Command**: `/path/to/your/chatspatial_env/bin/python` (path to your environment's Python)
   - **Args**: `["-m", "chatspatial"]` (run as Python module)
   - **Timeout**: Set to `3600` seconds in the separate timeout field
   - **Transport**: Select `stdio` (standard input/output)

   - Replace the command path with your actual ChatSpatial environment Python path
   - **Important**: Set the timeout to `3600` seconds to allow sufficient time for complex spatial analysis tasks
   - **Visualization**: All visualization functions are fully verified and production-ready

   **Finding Your Environment Path:**

   ```bash
   # Activate your ChatSpatial environment first
   conda activate chatspatial_env  # or: source chatspatial_env/bin/activate

   # Then find the Python path
   which python
   # Example output: /opt/anaconda3/envs/chatspatial_env/bin/python
   ```

3. **Save and Enable**:
   - Save the MCP server configuration
   - Enable the ChatSpatial MCP server
   - Restart Cherry Studio if required

4. **Use ChatSpatial in Cherry Studio**:
   - Start a new conversation
   - ChatSpatial tools will be automatically available
   - Cherry Studio's configurable timeout makes it ideal for computationally intensive spatial analysis

**Advantages of Cherry Studio over Claude Desktop:**

- **Configurable Timeout**: Set custom timeout (recommended: 3600s) for long-running analysis
- **Better Performance**: More suitable for computationally intensive spatial transcriptomics workflows
- **Stable Processing**: No interruptions during complex analysis tasks like GASTON, deconvolution, or large dataset processing

#### Example Prompts

Here are some example prompts to get started with ChatSpatial:

**Basic Workflow:**

- "Load my 10x Visium dataset from `/path/to/data.h5ad`"
- "Preprocess my data with custom filtering: filter genes in less than 10 cells and subsample to 1000 spots"
- "Visualize the spatial expression of gene Cd8a"

**Analysis and Visualization (Two-Step Process):**

- "Perform cell type annotation using marker genes" → "Visualize the cell type annotations"
- "Identify spatial domains using STAGATE method" → "Visualize the spatial domains"
- "Find spatial variable genes using GASTON with GLM-PCA preprocessing" → "Visualize GASTON isodepth map"
- "Analyze cell communication using LIANA+" → "Visualize cell communication results"
- "Integrate multiple spatial samples" → "Visualize the integrated UMAP colored by batch"
- "Deconvolve my spatial data using the NNLS method" → "Visualize deconvolution results"

**Advanced Analysis:**

- "Run spatial trajectory analysis"
- "Analyze spatial hot spots for immune genes using Getis-Ord Gi*"

## Available Tools

The server provides the following tools:

1. `load_data` - Load spatial transcriptomics data
2. `preprocess` - Preprocess data
3. `visualize` - Visualize data (supports all analysis results)
4. `annotate` - Cell type annotation
5. `identify_domains` - Spatial domain identification
6. `find_spatial_genes` - Spatial variable genes identification using GASTON, SpatialDE, or SPARK methods
7. `analyze_communication` - Cell-cell communication analysis with LIANA+
8. `analyze_spatial_data` - Spatial analysis (Moran's I, Getis-Ord Gi*, neighborhood, co-occurrence, Ripley's K, centrality)
9. `find_markers` - Differential expression analysis
10. `integrate_samples` - Multi-sample integration
11. `analyze_trajectory_pseudotime` - Trajectory pseudotime analysis
12. `analyze_velocity` - RNA velocity analysis
13. `deconvolve` - Spatial transcriptomics deconvolution

### Visualization Architecture

ChatSpatial uses a clean separation between analysis and visualization:

- All analysis tools focus solely on computation and return analysis results
- The `visualize` tool handles all visualization needs with various plot types
- This design ensures modularity and allows for flexible visualization options

To visualize analysis results, use the `visualize` tool with appropriate `plot_type`:

- `spatial` - Spatial gene expression or annotations
- `umap` - UMAP embeddings (e.g., for integrated data)
- `cell_communication` - Cell communication analysis results
- `spatial_domains` - Spatial domain identification results
- `gaston_isodepth`, `gaston_domains`, `gaston_genes` - GASTON analysis results
- And many more plot types for different analyses

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

### Spatially-aware Gene Set Enrichment with EnrichMap

ChatSpatial integrates EnrichMap for sophisticated gene set enrichment analysis with spatial awareness:

**Features:**

- **Spatial Smoothing**: Apply k-nearest neighbor smoothing to enrichment scores
- **Spatial Covariate Correction**: Use GAM to correct for spatial bias in enrichment
- **Batch Correction**: Support for multi-sample analysis with batch normalization
- **Gene Weighting**: Automatic or custom gene weights for enrichment calculation
- **Multi-signature Analysis**: Analyze multiple gene sets simultaneously
- **Comprehensive Visualization**: Spatial enrichment maps with customizable parameters

**Example Usage:**

```text
# Analysis (separate from visualization)
"Analyze T cell enrichment with gene set CD3D, CD3E, CD8A using spatial smoothing"
"Perform enrichment analysis for immune signatures with batch correction"

# Visualization (after running analyze_enrichment)
"Visualize enrichment plot for T_cell"
"Show spatial plot for T_cell_score"
"Create violin plot with plot_type enrichment for T_cell_score grouped by leiden"
```

### Spatial Variable Gene Analysis

ChatSpatial provides multiple methods for identifying spatial variable genes with different maturity levels:

#### GASTON (Recommended for Production) 🟢

**Graph Attention Spatial Transcriptomics Organizer Network** - Fully tested deep learning approach:

- **Production Ready**: Comprehensive testing and validation
- **Deep Learning Topology**: Neural network-based tissue organization learning
- **Isodepth Mapping**: Continuous topographic coordinates capturing tissue structure
- **Spatial Domain Detection**: Automatic identification of spatial domains
- **Easy Installation**: Optional pip installation with clear instructions

#### SpatialDE (Experimental Support) 🟡

**Spatial Differential Expression** - Requires compatibility fixes:

- **Implementation Status**: Code complete but needs environment setup
- **Compatibility Issues**: Requires fixes for modern scipy/pandas versions
- **Statistical Framework**: Gaussian process model for spatial gene expression
- **Installation**: Use `scripts/fixes/fix_spatialDE.py` for compatibility
- **Best For**: Users with specific research needs and technical expertise

#### SPARK (Advanced Users Only) 🔴

**Spatial Pattern Recognition** - Requires R environment:

- **Complex Dependencies**: Needs R environment + rpy2 + SPARK R package
- **Installation Challenge**: Multi-step setup process required
- **Count-Based Analysis**: Optimized for count-based spatial transcriptomics
- **Best For**: Users already working in R-based workflows

**Example Usage:**

```text
# GASTON (recommended for most users)
"Find spatial variable genes using GASTON with GLM-PCA preprocessing"
"Visualize GASTON isodepth map and spatial domains"

# SpatialDE (requires setup)
# First: python scripts/fixes/fix_spatialDE.py
"Find spatial variable genes using SpatialDE method"
"Visualize SpatialDE results with significance thresholds"

# SPARK (R environment required)
# First: pip install rpy2 && R -e "install.packages('SPARK')"
"Find spatial variable genes using SPARK with multiple kernels"
"Visualize SPARK spatial patterns"
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

### Optional Dependencies Mapping

Install optional packages to enable additional methods:

- SpaGCN (spatial domains): `pip install SpaGCN`
- STAGATE (spatial domains): `pip install STAGATE`
- BANKSY (spatial domains): `pip install banksy-utils`
- LIANA+ (cell communication): `pip install liana`
- esda/PySAL (Getis-Ord Gi*): `pip install esda libpysal`
- scvi-tools family (CellAssign, scANVI, DestVI, Stereoscope, VeloVI): `pip install scvi-tools`
- R-based methods (RCTD, SPOTlight): `pip install rpy2` and install R packages in your R environment
- Enrichment (GSEA/ssGSEA/EnrichMap): `pip install gseapy enrichmap`

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
  - GASTON - Deep learning spatial gene analysis (recommended)
  - glmpca - GLM-PCA preprocessing for GASTON
  - SpatialDE - Statistical spatial gene analysis
  - SPARK - Count-based spatial pattern recognition (requires R environment)

- **Specialized Analysis**:
  - liana - Cell communication analysis with LIANA+
  - cell2location - Spatial transcriptomics deconvolution
  - cellrank - Trajectory analysis
  - scvelo - RNA velocity analysis

## License

MIT

## Contributing

Pull Requests and Issues are welcome!
