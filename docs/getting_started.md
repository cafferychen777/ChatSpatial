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

**This is the technical part, but we'll walk through it step by step!**

#### Option 1: Quick Install (Recommended)

Open your terminal/command prompt and run these commands:

```bash
# Create a new environment (like a clean workspace)
conda create -n chatspatial python=3.11
conda activate chatspatial

# Get ChatSpatial
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial

# Install ChatSpatial
pip install -e .

# Verify installation
chatspatial --help
```

#### Option 2: Step-by-Step Install (If you're new to conda)

**Step 2.1: Install Miniconda (if you don't have it)**

- **Windows**: Download from [miniconda.io](https://docs.conda.io/en/latest/miniconda.html)
- **Mac**: `brew install miniconda` or download from website
- **Linux**: `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh`

**Step 2.2: Create Environment**
```bash
# Open terminal/command prompt
conda create -n chatspatial python=3.11
# Say 'y' when prompted
```

**Step 2.3: Activate Environment**
```bash
conda activate chatspatial
# You should see (chatspatial) in your prompt
```

**Step 2.4: Install Git (if needed)**
```bash
# Windows: Download from git-scm.com
# Mac: brew install git
# Linux: sudo apt install git (Ubuntu) or sudo yum install git (CentOS)
```

**Step 2.5: Clone and Install ChatSpatial**
```bash
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial
pip install -e .
```

**Step 2.6: Test Installation**
```bash
chatspatial --help
# Should show help information
```

#### Troubleshooting Installation

**Common Issues:**

1. **"conda not found"**
   - Solution: Install Miniconda first, restart terminal

2. **"git not found"**
   - Solution: Install Git for your operating system

3. **Permission errors**
   - Solution: Don't use `sudo` with conda/pip

4. **Network issues**
   - Solution: Check internet connection, try again

5. **Python version errors**
   - Solution: Make sure you're using Python 3.10 or 3.11

### Step 3: Connect to Claude Desktop

**This is the most important step - connecting ChatSpatial to Claude Desktop!**

#### Step 3.1: Find Your Python Path

```bash
# Make sure you're in the right environment
conda activate chatspatial

# Find your Python path
which python
# On Windows, use: where python
```

**Copy the full path!** It should look like:
- **Mac/Linux**: `/Users/yourname/miniconda3/envs/chatspatial/bin/python`
- **Windows**: `C:\Users\yourname\miniconda3\envs\chatspatial\python.exe`

#### Step 3.2: Locate Claude Desktop Config File

**Find the configuration file for your operating system:**

**Mac:**
```bash
# Open Finder, press Cmd+Shift+G, paste this path:
~/Library/Application Support/Claude/claude_desktop_config.json
```

**Windows:**
```bash
# Open File Explorer, paste this in the address bar:
%APPDATA%\Claude\claude_desktop_config.json
```

**Linux:**
```bash
# Use your file manager or terminal:
~/.config/Claude/claude_desktop_config.json
```

#### Step 3.3: Edit Configuration File

**If the file doesn't exist, create it!**

Open the file in any text editor and add this configuration:

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

**⚠️ Important:** Replace `/your/python/path/here` with the actual path you copied in Step 3.1!

**Example configurations:**

**Mac example:**
```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/Users/john/miniconda3/envs/chatspatial/bin/python",
      "args": ["-m", "chatspatial"],
      "env": {}
    }
  }
}
```

**Windows example:**
```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "C:\\Users\\john\\miniconda3\\envs\\chatspatial\\python.exe",
      "args": ["-m", "chatspatial"],
      "env": {}
    }
  }
}
```

#### Step 3.4: Restart Claude Desktop

1. **Close Claude Desktop completely**
2. **Reopen Claude Desktop**
3. **Look for the 🔬 icon** - this means ChatSpatial is connected!

#### Troubleshooting MCP Connection

**Common Issues:**

1. **"No MCP servers found"**
   - Check the config file path is correct
   - Verify JSON syntax (use a JSON validator online)
   - Make sure Python path is correct

2. **"Server failed to start"**
   - Test Python path in terminal: `/your/python/path -m chatspatial --help`
   - Check environment is activated
   - Verify ChatSpatial is installed

3. **"Permission denied"**
   - Make sure Python path is executable
   - Don't use spaces in paths (use quotes if needed)

4. **JSON syntax errors**
   - Use a JSON validator to check syntax
   - Make sure all brackets and quotes match
   - No trailing commas

### Step 4: Verify Everything Works

**Let's make sure everything is set up correctly!**

#### Test 1: Check ChatSpatial Installation
```bash
# Activate environment
conda activate chatspatial

# Test ChatSpatial
python -m chatspatial --help
# Should show help information without errors
```

#### Test 2: Check MCP Connection
1. **Open Claude Desktop**
2. **Look for the 🔬 icon** in the tool panel
3. **Start a new conversation**
4. **Type**: "Do you have access to ChatSpatial tools?"

**If successful, Claude should respond with something like:**
> "Yes! I have access to ChatSpatial's spatial transcriptomics analysis tools. I can help you load data, identify spatial domains, annotate cell types, and much more!"

#### Test 3: Quick Function Test
**Ask Claude**: "Can you list the available ChatSpatial tools?"

**You should see tools like:**
- load_data
- preprocess_data
- identify_spatial_domains
- annotate_cells
- visualize_data
- And more!

### 🎉 Success! You're Ready to Analyze!

If all tests pass, you're ready to start analyzing spatial transcriptomics data through natural conversation!

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

