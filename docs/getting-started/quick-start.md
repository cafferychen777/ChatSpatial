# Getting Started with ChatSpatial

This guide demonstrates how to set up ChatSpatial and analyze spatial transcriptomics data through natural language queries in Claude Desktop, eliminating the need for complex coding.

## What You Will Achieve

By the end of this guide, you will be able to:
- **Ask questions** about your spatial data in plain English
- **Analyze tissue architecture** without writing code
- **Generate visualizations** automatically
- **Discover biological insights** through conversation

## 🚀 Quick Start (5 Minutes)

### Step 1: Get Claude Desktop

**New to Claude?** The following steps will guide you through the setup process.

1. 🌐 **Visit**: [claude.ai](https://claude.ai)
2. 📱 **Download**: Claude Desktop for your computer
3. 👤 **Sign up**: Create your free Anthropic account

### Step 2: Install ChatSpatial

**This section covers the technical installation steps, presented in detail.**

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

#### Option 2: Step-by-Step Install (If you are new to conda)

**Step 2.1: Install Miniconda (if you do not have it)**

- **Windows**: Download from [miniconda.io](https://docs.conda.io/en/latest/miniconda.html)
- **Mac**: `brew install miniconda` or download from website
- **Linux**: `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh`

**Step 2.2: Create Environment**

```bash
# Open terminal/command prompt
conda create -n chatspatial python=3.11
# Say 'y' when prompted
```

**Step 2.3: Initialize Conda (First Time Only)**

```bash
# Initialize conda for your shell (only needed once)
conda init bash  # For bash users
conda init zsh   # For zsh users (Mac default)

# Restart your terminal or run:
source ~/.bashrc  # Linux
source ~/.zshrc   # Mac
```

**Step 2.4: Activate Environment**

```bash
conda activate chatspatial
# You should see (chatspatial) in your prompt
```

**Step 2.5: Verify Environment Activation**

```bash
# Test that you are in the right environment
which python
# Should show: /path/to/miniconda3/envs/chatspatial/bin/python

python --version
# Should show: Python 3.11.x
```

**⚠️ If conda activate does not work:**

```bash
# Alternative method - use full path directly
/Users/yourname/miniconda3/envs/chatspatial/bin/python --version
# Replace 'yourname' with your actual username
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

**Step 2.6: Test ChatSpatial Installation**

**Test 1: Basic Import (may take 10-30 seconds)**

```bash
# Test ChatSpatial import - this loads scientific libraries
python -c 'import chatspatial; print("ChatSpatial import successful")'
# Note: Use single quotes to avoid bash issues with special characters
```

**Test 2: Command Line Interface**

```bash
# Test the command line interface
python -m chatspatial --help
# Should show help information immediately
```

**Expected Output:**

```text
Usage: python -m chatspatial server [OPTIONS]

  Start the ChatSpatial server.
  ...
```

**⚠️ Important Notes:**

- **Import time**: First `import chatspatial` may take 10-30 seconds (loading numpy, scipy, etc.)
- **Command help**: `python -m chatspatial --help` should be fast (< 5 seconds)
- **Use single quotes**: Avoid bash interpretation issues with `'` instead of `"`

#### Troubleshooting Installation

**Common Issues Based on Real Testing:**

1. **"conda activate doesn't work"**

   ```bash
   # Solution: Initialize conda first
   conda init bash  # or zsh
   source ~/.bashrc  # or ~/.zshrc
   # Or use full path: /path/to/miniconda3/envs/chatspatial/bin/python
   ```

2. **"bash: !': event not found"**

   ```bash
   # Problem: Using double quotes with ! character
   # Wrong: python -c "print('Hello!')"
   # Right: python -c 'print("Hello!")'
   ```

3. **"Import takes forever"**
   - **Normal**: First `import chatspatial` takes 10-30 seconds
   - **Wait patiently**: Scientific libraries (numpy, scipy) are loading
   - **Subsequent imports**: Will be faster

4. **"conda not found"**
   - Install Miniconda first, restart terminal
   - Make sure conda is in your PATH

5. **"Python version errors"**
   - Verify: `python --version` shows Python 3.11.x
   - Check environment: `which python` points to chatspatial env

### Step 3: Connect to Claude Desktop

**This is the most important step - connecting ChatSpatial to Claude Desktop!**

#### Step 3.1: Find Your Python Path

```bash
# Make sure you are in the right environment
conda activate chatspatial

# Verify you are in the correct environment
which python
# On Windows, use: where python
```

**Copy the COMPLETE path!** Real examples:

- **Mac/Linux**: `/Users/john/miniconda3/envs/chatspatial/bin/python`
- **Windows**: `C:\Users\john\miniconda3\envs\chatspatial\python.exe`

**⚠️ Important**:

- Copy the **entire path** including `/bin/python` or `\python.exe`
- Replace `john` with your actual username
- The path must point to the chatspatial environment, not base conda

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

**If the file does not exist, create it.**

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

**Common Issues Based on Real Testing:**

1. **"No MCP servers found"**
   ```bash
   # Test your Python path first
   /your/python/path -c 'import chatspatial; print("Works")'
   /your/python/path -m chatspatial --help
   ```
   - Verify JSON syntax (use jsonlint.com)
   - Check config file location is correct
   - Ensure Python path is complete (includes `/bin/python`)

2. **"Server failed to start"**
   ```bash
   # Debug step by step
   /your/python/path -c 'import chatspatial'  # Should work in 10-30s
   /your/python/path -m chatspatial --help    # Should be fast
   /your/python/path -m chatspatial server --transport stdio  # Should start
   ```
   - **First import**: Takes 10-30 seconds (normal!)
   - **Command help**: Should be fast (< 5 seconds)
   - **Server start**: Should show "Starting ChatSpatial server..."

3. **"JSON syntax errors"**
   - Use single quotes in test commands: `'import chatspatial'`
   - Validate JSON at jsonlint.com
   - No trailing commas in JSON
   - Match all brackets and quotes

4. **"Path issues"**
   - Use complete path: `/Users/john/miniconda3/envs/chatspatial/bin/python`
   - Not just: `/Users/john/miniconda3/envs/chatspatial/`
   - Test path works: `which python` when environment is active

### Step 4: Verify Everything Works

**The following tests verify that everything is set up correctly.**

#### Test 1: Verify ChatSpatial Installation

```bash
# Activate environment
conda activate chatspatial

# Test 1a: Basic import (may take 10-30 seconds)
python -c 'import chatspatial; print("Import successful")'

# Test 1b: Command line interface (should be fast)
python -m chatspatial --help

# Test 1c: Server startup test
python -m chatspatial server --transport stdio &
# You should see: "Starting ChatSpatial server with stdio transport..."
# Press Ctrl+C to stop
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

### Success! You Are Ready to Analyze

If all tests pass, you are ready to start analyzing spatial transcriptomics data through natural conversation.

## 🎥 New to MCP? Watch This First!

[![MCP Introduction](https://img.shields.io/badge/▶️_Watch-What_is_MCP?-blue?style=flat-square&logo=youtube)](https://www.youtube.com/watch?v=sfCBCyNyw7U)

**Learn more:**
- 📚 [What is MCP?](what_is_mcp.md) - Detailed explanation
- 📰 [Anthropic's MCP Announcement](https://www.anthropic.com/news/model-context-protocol)
- 💻 [MCP GitHub Repository](https://github.com/modelcontextprotocol)

## 💬 Your First Conversation

Once ChatSpatial is connected, here's how easy it is to analyze your data:

---

**👤 You:** "Hi! I am new to spatial transcriptomics. Can you help me analyze my Visium data?"

**🤖 Claude:** "I can help you explore your spatial transcriptomics data. ChatSpatial makes it easy - no coding required. Just tell me what you would like to discover about your tissue.

What kind of sample do you have, and where is your data file located?"

**👤 You:** "I have a mouse brain sample. The file is at `/Users/myname/data/brain_sample.h5ad`"

**🤖 Claude:** "Perfect! Let me load your mouse brain data and take a look at what we are working with."

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

**👤 You:** "I would like to see the different brain regions!"

**🤖 Claude:** "Great choice! Let me identify the spatial domains in your brain tissue..."

*[Analysis happens automatically]*

**🤖 Claude:** "Good! I found 6 distinct brain regions:
- Cortical layers (outer brain)
- Hippocampus (memory center)
- White matter (connecting fibers)
- Ventricular zone (stem cells)
- Corpus callosum (brain bridge)
- Boundary regions

Here is a spatial map showing each region in different colors."

---

### See How Natural It Is

- **No complex commands** to remember
- **No parameter tuning** required
- **No coding experience** needed
- **Ask questions** in plain English
- **Get instant insights** about your tissue
- **Visualizations** automatically generated

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
- **Other LLM Agents**: With MCP support

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
result = load_data(
    data_path="data/demo_datasets/mouse_brain_visium.h5ad",
    name="mouse_brain"
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

- **[Basic Spatial Analysis](../tutorials/core/basic_spatial_analysis.md)**: Complete workflow
- **[Cell Annotation Guide](../tutorials/cell_type_annotation.md)**: Multiple annotation methods
- **[Visualization Gallery](../tutorials/visualization_gallery.md)**: All plot types

### API Reference

- **[Tool Reference](../reference/api/README.md)**: Complete MCP tool documentation
- **[Data Models](../reference/api/data_models.md)**: Parameter schemas and data structures
- **[Error Handling](../reference/api/error_handling.md)**: Troubleshooting guide

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

