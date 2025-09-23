---
layout: default
title: Quick Start
parent: Getting Started
nav_order: 2
---

# Quick Start Guide
{: .no_toc }

Get ChatSpatial running in minutes and perform your first spatial analysis.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

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

{: .highlight }
💡 **Virtual environments prevent conflicts:** Always use a virtual environment for Python projects.

{: .highlight }
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
      "command": "/Users/apple/Projects/chatspatial_env/bin/python",
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

### Download Sample Data First

Before running analyses, download sample datasets from [ChatSpatial Releases](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.2.0-data):
- `destvi_reference_small.h5ad` (7.8MB)
- `destvi_spatial_small.h5ad` (10MB)

### Load and Explore Data

```
Load /Users/yourname/Downloads/destvi_spatial_small.h5ad and show me basic information about it
```

> ⚠️ **IMPORTANT**: Always use **absolute paths** (starting with `/`) when loading data files.

Expected response: ChatSpatial will load the dataset and provide summary statistics.

### Spatial Visualization

```
Create a spatial plot showing the tissue structure and gene expression patterns
```

This will generate spatial plots of your data.

### Identify Spatial Domains

```
Identify spatial domains using SpaGCN and visualize the results
```

ChatSpatial will perform spatial domain identification and create visualization plots.

### Cell Type Annotation

```
Annotate cell types using marker genes and show the spatial distribution
```

This performs automated cell type annotation and spatial mapping.

## Example Workflows

### Basic Spatial Analysis Workflow

```
# 1. Load data (use absolute path!)
Load /Users/yourname/Downloads/destvi_spatial_small.h5ad

# 2. Quality control  
Show me quality control metrics and filter low-quality spots

# 3. Dimensionality reduction
Perform PCA and UMAP embedding on the data

# 4. Clustering
Cluster the spots using Leiden clustering

# 5. Spatial visualization
Create spatial plots colored by clusters

# 6. Find marker genes
Find marker genes for each cluster

# 7. Spatial domains
Identify spatial domains using SpaGCN
```

### Advanced Analysis Workflow

```
# 1. Load and preprocess (use absolute paths!)
Load /Users/yourname/Downloads/destvi_spatial_small.h5ad and /Users/yourname/Downloads/destvi_reference_small.h5ad and perform quality control

# 2. Cell type annotation
Annotate cell types using the marker gene approach

# 3. Spatial domains
Identify spatial domains and visualize them

# 4. Cell communication
Analyze cell-cell communication patterns using LIANA

# 5. Spatial variable genes
Find spatially variable genes using SPARK

# 6. Trajectory analysis
Perform RNA velocity analysis for developmental trajectories
```

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

{: .note }
**Problem:** MCP server fails to start

**Solutions:**
1. Check Python environment: `python --version`
2. Verify installation: `pip list | grep chatspatial`
3. Check dependencies: `python -c "import chatspatial; print('OK')"`

### Data Loading Errors

{: .note }
**Problem:** Cannot load spatial data

**Solutions:**
1. Verify file paths are correct
2. Check supported formats (H5AD, H5, MTX, Visium)
3. Ensure data contains spatial coordinates

### Memory Issues

{: .note }
**Problem:** Analysis fails with memory errors

**Solutions:**
1. Use smaller datasets for testing
2. Increase system memory
3. Use chunked processing for large datasets

### Missing Dependencies

{: .note }
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

1. **Explore Tutorials**: Check out the [tutorials section](../tutorials/) for detailed workflows
2. **Learn Analysis Methods**: Dive into specific analysis types in the [reference guide](../reference/)
3. **Try Advanced Features**: Experiment with [advanced analysis methods](../tutorials/advanced/)
4. **Join the Community**: Participate in [discussions](https://github.com/cafferychen777/ChatSpatial/discussions)

## Getting Help

If you encounter issues:

1. **Check the logs**: Look for error messages in the console output
2. **Review documentation**: Browse the [troubleshooting guide](../reference/troubleshooting/)
3. **Search issues**: Check [existing GitHub issues](https://github.com/cafferychen777/ChatSpatial/issues)
4. **Ask for help**: Create a new issue with detailed error information

---

**Next:** [Configuration Guide](configuration.html) to customize your setup