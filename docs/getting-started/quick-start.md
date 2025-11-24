
title: Quick Start
---

# Quick Start Guide

Get ChatSpatial running in minutes and perform your first spatial analysis.

## Table of contents

## Prerequisites

Before starting, ensure you have:
- Python 3.10+ installed (required for MCP)
- Claude Desktop or compatible MCP client
- Git for cloning the repository

## Step 1: Installation

### Create Virtual Environment First

```bash
# Clone the repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd chatspatial

# Create and activate virtual environment
python3 -m venv chatspatial_env
source chatspatial_env/bin/activate  # macOS/Linux
# chatspatial_env\Scripts\activate   # Windows

# Install ChatSpatial with all features
pip install -e ".[full]"
```

💡 **Virtual environments prevent conflicts:** Always use a virtual environment for Python projects.

💡 For faster installation (80% features): `pip install -e .`

## Step 2: Configure MCP Client

### For Claude Desktop

1. **Find your virtual environment Python path:**
```bash
# In your activated virtual environment
which python
# Copy this path - you'll need it next
```

2. **Edit Claude Desktop configuration:**
   - Location: `~/Library/Application Support/Claude/claude_desktop_config.json`
   - Add ChatSpatial with your virtual environment path:

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial"]
    }
  }
}
```

**Example with real path:**
```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/Users/yourname/Projects/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial"]
    }
  }
}
```

### For Other MCP Clients

ChatSpatial supports the standard MCP protocol. Refer to your client's documentation for specific configuration steps.

## Step 3: Start the Server

### Option A: Automatic (with Claude Desktop)
The server starts automatically when you interact with ChatSpatial in Claude Desktop.

### Option B: Manual Start
```bash
# Make sure virtual environment is activated
source chatspatial_env/bin/activate  # macOS/Linux

# Start the MCP server manually
python -m chatspatial

# Or start with HTTP transport
python -m chatspatial server --transport sse --port 8000
```

## Step 4: Your First Analysis

Once ChatSpatial is configured, try these example analyses:

> 💬 **How to use ChatSpatial**: All examples below are **natural language commands** you type in your Claude chat (Desktop or Code).
>
> **First**, tell Claude: "Use ChatSpatial MCP for all my spatial transcriptomics analysis" - this ensures Claude uses the MCP tools instead of writing scripts.
>
> **Then**, simply chat in natural language! ChatSpatial MCP server interprets your requests and automatically runs the appropriate analysis tools. **No Python code needed!**

### Download Sample Data First

Before running analyses, download sample datasets from [ChatSpatial Releases](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data):
- `card_reference_filtered.h5ad` (36MB - pancreatic reference with 9 cell types)
- `card_spatial.h5ad` (7.7MB - spatial data with clear tissue domains)

### Load and Explore Data

**In your Claude chat, type:**

```
Load /Users/yourname/Downloads/card_spatial.h5ad and show me basic information about it
```

> ⚠️ **IMPORTANT**: Always use **absolute paths** (starting with `/`) when loading data files.

**Expected response**: ChatSpatial MCP will load the dataset and provide summary statistics including number of spots, genes, and quality metrics.

### Spatial Visualization

**In your Claude chat, type:**

```
Create a spatial plot showing the tissue structure and gene expression patterns
```

**What happens**: ChatSpatial MCP will generate spatial plots of your data with tissue coordinates and expression overlays.

### Identify Spatial Domains

**In your Claude chat, type:**

```
Identify spatial domains using spagcn method and visualize the results
```

**What happens**: ChatSpatial MCP will perform spatial domain identification and create visualization plots showing distinct tissue regions.

### Cell Type Annotation

**In your Claude chat, type:**

```
Annotate cell types using marker genes and show the spatial distribution
```

**What happens**: ChatSpatial MCP performs automated cell type annotation and creates spatial maps showing cell type distributions.

## Example Workflows

> 💬 **Conversational Analysis**: Copy and paste these natural language requests into your Claude chat. ChatSpatial MCP will execute each step automatically.

### Basic Spatial Analysis Workflow

**In your Claude chat, have this conversation:**

```
1. Load /Users/yourname/Downloads/card_spatial.h5ad

2. Show me quality control metrics and filter low-quality spots

3. Perform PCA and UMAP embedding on the data

4. Cluster the spots using Leiden clustering

5. Create spatial plots colored by clusters

6. Find marker genes for each cluster

7. Identify spatial domains using spagcn method
```

**What happens**: ChatSpatial MCP guides you through a complete basic analysis workflow, from data loading to spatial domain identification.

### Advanced Analysis Workflow

**In your Claude chat, have this conversation:**

```
1. Load /Users/yourname/Downloads/card_spatial.h5ad and /Users/yourname/Downloads/card_reference_filtered.h5ad and perform quality control

2. Annotate cell types using the reference data

3. Identify spatial domains and visualize them

4. Analyze cell-cell communication patterns using liana method

5. Find spatially variable genes using sparkx method

6. Deconvolve the spatial data using the reference with Cell2location
```

**What happens**: ChatSpatial MCP performs an advanced multi-method analysis including deconvolution, cell communication, and spatial statistics.

## Understanding Results

ChatSpatial returns results in several formats:

### Visualizations
- **Spatial plots**: Show gene expression or annotations overlaid on tissue coordinates
- **UMAP plots**: Display dimensionality reduction results
- **Heatmaps**: Show expression patterns across cell types or spatial regions

### Data Tables
- **Gene lists**: Differentially expressed or spatially variable genes
- **Statistics**: Spatial analysis results and significance tests
- **Annotations**: Cell type predictions and confidence scores

### Analysis Reports
- **Summary statistics**: Dataset characteristics and analysis parameters
- **Method details**: Algorithms used and parameter settings
- **Quality metrics**: Assessment of analysis quality and reliability

## Common Issues and Solutions

### Server Won't Start

**Problem:** MCP server fails to start

**Solutions:**
1. Check Python environment: `python --version`
2. Verify installation: `pip list | grep chatspatial`
3. Check dependencies: `python -c "import chatspatial; print('OK')"`

### Data Loading Errors

**Problem:** Cannot load spatial data

**Solutions:**
1. Verify file paths are correct
2. Check supported formats (H5AD, H5, MTX, Visium)
3. Ensure data contains spatial coordinates

### Memory Issues

**Problem:** Analysis fails with memory errors

**Solutions:**
1. Use smaller datasets for testing
2. Increase system memory
3. Use chunked processing for large datasets

### Missing Dependencies

**Problem:** Method not available errors

**Solutions:**
1. Check installed packages: `pip list | grep chatspatial`
2. Install missing dependencies: `pip install -e ".[advanced]"`
3. Use alternative methods when dependencies are missing

## Validation and Testing

### Test Your Setup

```bash
# Test basic functionality
python -c "import chatspatial; print('✅ Installation OK')"

# Test with demo data
python -c "
from chatspatial.utils.data_loader import load_demo_data
data = load_demo_data('visium_demo')
print(f'Loaded data with {data.n_obs} spots and {data.n_vars} genes')
"
```

### Performance Benchmarks

Run benchmark tests to ensure optimal performance:

```bash
# Quick performance test
python scripts/benchmark_basic.py

# Comprehensive benchmark
python scripts/benchmark_full.py
```

## Next Steps

Now that you have ChatSpatial running:

1. **Explore Tutorials**: Check out the [tutorials section](../tutorials/index.md) for detailed workflows
2. **Learn Analysis Methods**: Dive into specific analysis types in the [reference guide](../reference/index.md)
3. **Try Advanced Features**: Experiment with [advanced analysis methods](../tutorials/advanced/index.md)
4. **Join the Community**: Participate in [discussions](https://github.com/cafferychen777/ChatSpatial/discussions)

## Getting Help

If you encounter issues:

1. **Check the logs**: Look for error messages in the console output
2. **Review documentation**: Browse the [troubleshooting guide](../reference/troubleshooting/index.md)
3. **Search issues**: Check [existing GitHub issues](https://github.com/cafferychen777/ChatSpatial/issues)
4. **Ask for help**: Create a new issue with detailed error information

---

**Next:** [Configuration Guide](configuration.md) to customize your setup