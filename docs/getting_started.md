# Getting Started with ChatSpatial

Transform your spatial transcriptomics analysis from complex coding to natural conversation! This guide shows you how to set up ChatSpatial and start analyzing your data through simple questions in Claude Desktop.

## 🎯 What You'll Achieve

By the end of this guide, you'll be able to:
- 💬 **Ask questions** about your spatial data in plain English
- 🧬 **Analyze tissue architecture** without writing code
- 🎨 **Generate beautiful visualizations** automatically
- 🔬 **Discover biological insights** through conversation

## 🚀 Quick Start (5 Minutes)

### Step 1: Get Claude Desktop

**New to Claude?** No problem!

1. 🌐 **Visit**: [claude.ai](https://claude.ai)
2. 📱 **Download**: Claude Desktop for your computer
3. 👤 **Sign up**: Create your free Anthropic account

### Step 2: Install ChatSpatial

**Don't worry - this is easier than it looks!**

Open your terminal/command prompt and run these commands:

```bash
# Create a new environment (like a clean workspace)
conda create -n chatspatial python=3.11
conda activate chatspatial

# Get ChatSpatial
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial
pip install -e .
```

### Step 3: Connect to Claude Desktop

**This is where the magic happens!**

1. **Find your Python path**:
   ```bash
   conda activate chatspatial
   which python
   ```
   Copy this path (something like `/Users/yourname/miniconda3/envs/chatspatial/bin/python`)

2. **Configure Claude Desktop**:
   - **Mac**: Open `~/Library/Application Support/Claude/claude_desktop_config.json`
   - **Windows**: Open `%APPDATA%\Claude\claude_desktop_config.json`
   - **Linux**: Open `~/.config/Claude/claude_desktop_config.json`

3. **Add this configuration** (replace the Python path with yours):
   ```json
   {
     "mcpServers": {
       "chatspatial": {
         "command": "/your/python/path/here",
         "args": ["-m", "chatspatial"],
         "env": {}
       }
     }
   }
   ```

4. **Restart Claude Desktop**

### 🎉 You're Ready!

Look for the 🔬 ChatSpatial tools in Claude Desktop. You should see spatial analysis tools available!

## 🎥 New to MCP? Watch This First!

[![MCP Introduction](https://img.shields.io/badge/▶️_Watch-What_is_MCP?-blue?style=flat-square&logo=youtube)](https://www.youtube.com/watch?v=sfCBCyNyw7U)

**Learn more:**
- 📚 [What is MCP?](what_is_mcp.md) - Detailed explanation
- 📰 [Anthropic's MCP Announcement](https://www.anthropic.com/news/model-context-protocol)
- 💻 [MCP GitHub Repository](https://github.com/modelcontextprotocol)

## 💬 Your First Conversation

Once ChatSpatial is connected, here's how easy it is to analyze your data:

---

**👤 You:** "Hi! I'm new to spatial transcriptomics. Can you help me analyze my Visium data?"

**🤖 Claude:** "I'd love to help you explore your spatial transcriptomics data! ChatSpatial makes it easy - no coding required. Just tell me what you'd like to discover about your tissue.

What kind of sample do you have, and where is your data file located?"

**👤 You:** "I have a mouse brain sample. The file is at `/Users/myname/data/brain_sample.h5ad`"

**🤖 Claude:** "Perfect! Let me load your mouse brain data and take a look at what we're working with."

*[ChatSpatial loads and analyzes the data]*

**🤖 Claude:** "Excellent! Your mouse brain sample looks great:
- 2,845 spots across the tissue
- 28,691 genes measured
- High-quality spatial coordinates included

What would you like to explore? I can help you:
- 🧠 Identify different brain regions
- 🔬 Find cell types in each area
- 📊 Discover spatially variable genes
- 🎨 Create beautiful visualizations
- 💬 Analyze cell-cell communication

What interests you most?"

**👤 You:** "I'd love to see the different brain regions!"

**🤖 Claude:** "Great choice! Let me identify the spatial domains in your brain tissue..."

*[Analysis happens automatically]*

**🤖 Claude:** "Amazing! I found 6 distinct brain regions:
- Cortical layers (outer brain)
- Hippocampus (memory center)
- White matter (connecting fibers)
- Ventricular zone (stem cells)
- Corpus callosum (brain bridge)
- Boundary regions

Here's a beautiful spatial map showing each region in different colors! 🧠✨"

---

### 🎯 See How Natural It Is?

- ❌ **No complex commands** to remember
- ❌ **No parameter tuning** required
- ❌ **No coding experience** needed
- ✅ **Just ask questions** in plain English!
- ✅ **Get instant insights** about your tissue
- ✅ **Beautiful visualizations** automatically generated

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

