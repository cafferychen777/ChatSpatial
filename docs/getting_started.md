# Getting Started with ChatSpatial

This comprehensive guide will help you install, configure, and run your first spatial transcriptomics analysis with ChatSpatial.

## Prerequisites

Before installing ChatSpatial, ensure you have:

- **Python 3.10+** (recommended: 3.10 or 3.11)
- **Conda** or **Miniconda** for environment management
- **Git** for cloning the repository
- **8GB+ RAM** for typical spatial datasets
- **Claude Desktop** or another MCP-compatible client

### 📥 Download Claude Desktop

If you don't have Claude Desktop yet:

1. **Visit**: [claude.ai](https://claude.ai)
2. **Download**: Claude Desktop for your operating system
3. **Sign up**: Create an Anthropic account if needed

### 🎥 Learn About MCP First

New to Model Context Protocol? Watch this quick introduction:

[![MCP Introduction](https://img.shields.io/badge/▶️_Watch-What_is_MCP?-blue?style=flat-square&logo=youtube)](https://www.youtube.com/watch?v=sfCBCyNyw7U)

**Additional Resources:**
- 📚 [Official MCP Documentation](https://modelcontextprotocol.io)
- 📰 [Anthropic's MCP Announcement](https://www.anthropic.com/news/model-context-protocol)
- 💻 [MCP GitHub Repository](https://github.com/modelcontextprotocol)

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial
```

### Step 2: Create Python Environment

```bash
# Create a new conda environment
conda create -n chatspatial python=3.10
conda activate chatspatial
```

### Step 3: Install Dependencies

```bash
# Install ChatSpatial in development mode
pip install -e .

# Verify installation
chatspatial --help
```

### Step 4: Install Optional Dependencies

For full functionality, install additional packages:

```bash
# R dependencies for scType
conda install -c conda-forge r-base r-essentials

# Additional spatial analysis tools
pip install SpaGCN STAGATE liana cellphonedb

# GPU support (optional, for faster processing)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

## MCP Configuration

### Claude Desktop Setup

1. **Locate Configuration File**:
   - **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
   - **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
   - **Linux**: `~/.config/Claude/claude_desktop_config.json`

2. **Add ChatSpatial Server**:

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/your/conda/envs/chatspatial/bin/python",
      "args": ["-m", "chatspatial"],
      "env": {
        "PYTHONPATH": "/path/to/ChatSpatial"
      }
    }
  }
}
```

3. **Find Your Python Path**:

```bash
# Activate your environment
conda activate chatspatial

# Get the Python path
which python
# Example output: /Users/username/miniconda3/envs/chatspatial/bin/python
```

### Alternative MCP Clients

ChatSpatial works with any MCP-compatible client:

- **MCP Inspector**: For development and testing
- **Custom Applications**: Using the MCP SDK
- **Other AI Assistants**: With MCP support

## Data Preparation

### Download Sample Data

ChatSpatial includes scripts to download standard datasets:

```bash
# Download demo datasets
python data/scripts/download_standard_datasets.py

# Verify data
ls data/demo_datasets/
```

### Supported Data Formats

- **AnnData (H5AD)**: Primary format, includes spatial coordinates
- **10x Visium**: Space Ranger outputs
- **CSV + Coordinates**: Expression matrix with spatial coordinates
- **H5/HDF5**: Hierarchical data format
- **Zarr**: Cloud-optimized format

### Data Requirements

Your spatial transcriptomics data should include:

1. **Gene Expression Matrix**: Genes × Spots/Cells
2. **Spatial Coordinates**: X, Y positions for each spot/cell
3. **Metadata** (optional): Cell types, batch information, etc.

## First Analysis Walkthrough

### Step 1: Start Claude Desktop

1. Open Claude Desktop
2. Verify ChatSpatial appears in the MCP servers list
3. Look for the 🔬 icon indicating spatial analysis tools

### Step 2: Load Your Data

```python
# Load a Visium dataset
load_data(
    file_path="data/demo_datasets/mouse_brain_visium.h5ad",
    data_id="mouse_brain"
)
```

### Step 3: Explore Data Structure

```python
# Get basic information about your dataset
analyze_spatial_data(
    data_id="mouse_brain",
    analysis_type="basic_stats"
)
```

### Step 4: Preprocessing

```python
# Standard preprocessing pipeline
preprocess_data(
    data_id="mouse_brain",
    normalize_total=True,
    log1p=True,
    highly_variable_genes=True,
    n_top_genes=2000
)
```

### Step 5: Spatial Domain Identification

```python
# Identify spatial domains using SpaGCN
identify_spatial_domains(
    data_id="mouse_brain",
    method="spagcn",
    n_clusters=7,
    resolution=1.0
)
```

### Step 6: Visualization

```python
# Create spatial domain visualization
visualize_data(
    data_id="mouse_brain",
    plot_type="spatial_domains",
    color_by="spatial_domains",
    title="Mouse Brain Spatial Domains"
)
```

## Understanding the Results

### Data Structure

After loading, your data contains:

- **`.X`**: Gene expression matrix
- **`.obs`**: Cell/spot metadata (including spatial domains)
- **`.var`**: Gene metadata
- **`.obsm['spatial']`**: Spatial coordinates
- **`.uns`**: Analysis results and parameters

### Spatial Domains

The `identify_spatial_domains` tool adds:

- **`spatial_domains`**: Cluster assignments in `.obs`
- **`spatial_domain_stats`**: Cluster statistics in `.uns`
- **`spatial_embeddings`**: Low-dimensional representations

### Visualization Outputs

ChatSpatial returns MCP Image objects that display directly in Claude Desktop:

- **High-resolution plots**: 300 DPI for publication quality
- **Interactive elements**: Hover information and zoom
- **Multiple formats**: PNG, SVG, PDF support
- **Customizable styling**: Colors, themes, annotations

## Next Steps

### Advanced Analysis

1. **Cell Type Annotation**:
   ```python
   annotate_cells(data_id="mouse_brain", method="tangram")
   ```

2. **Cell Communication Analysis**:
   ```python
   analyze_cell_communication(data_id="mouse_brain", method="liana")
   ```

3. **Spatial Variable Genes**:
   ```python
   identify_spatial_genes(data_id="mouse_brain", method="gaston")
   ```

### Explore Tutorials

- **[Basic Spatial Analysis](tutorials/basic_spatial_analysis.md)**: Complete workflow
- **[Cell Annotation Guide](tutorials/cell_annotation.md)**: Multiple annotation methods
- **[Visualization Gallery](tutorials/visualization_gallery.md)**: All plot types

### API Reference

- **[Tool Reference](api/README.md)**: Complete MCP tool documentation
- **[Parameter Guide](api/parameters.md)**: Detailed parameter descriptions
- **[Error Handling](api/error_handling.md)**: Troubleshooting guide

## Troubleshooting

### Common Issues

1. **Import Errors**: Ensure all dependencies are installed in the correct environment
2. **Memory Issues**: Use data subsampling for large datasets
3. **MCP Connection**: Verify Python path and environment variables
4. **R Dependencies**: Install R packages for scType functionality

### Getting Help

- **Documentation**: Browse this site for detailed guides
- **GitHub Issues**: Report bugs and request features
- **Discussions**: Ask questions and share experiences

### Performance Tips

- **Use SSD storage** for faster data loading
- **Increase memory** for large datasets (>50K cells)
- **Enable GPU** for deep learning methods
- **Subsample data** for initial exploration

