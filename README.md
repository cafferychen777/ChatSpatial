<div align="center">

<img src="https://raw.githubusercontent.com/cafferychen777/ChatSpatial/main/assets/images/chatspatial-logo.png" alt="ChatSpatial Logo" width="300"/>

# ChatSpatial

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/) [![MCP Protocol](https://img.shields.io/badge/MCP-v2025.03.26-green.svg)](https://modelcontextprotocol.io) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Docs](https://img.shields.io/badge/docs-available-blue)](https://cafferychen777.github.io/ChatSpatial/)

### Agentic Workflow Orchestration for Spatial Transcriptomics Analysis

</div>

**Eliminate the implementation tax. Focus on biological insight.**

ChatSpatial is an agentic workflow orchestration platform that integrates 60 state-of-the-art methods from fragmented Python and R ecosystems into a unified conversational interface. Built on the Model Context Protocol (MCP), it enables human-steered discovery through natural language in Claude Desktop or Claude Code, eliminating the need for manual data conversion and complex programming.

**🎯 Example: Analyze spatial transcriptomics data through conversation with Claude**

```text
👤 "Load my 10x Visium dataset and identify spatial domains"
🤖 ✅ Loaded 3,456 spots, 18,078 genes
    ✅ Identified 7 spatial domains using SpaGCN
    ✅ Generated spatial domain visualization

👤 "Find marker genes for domain 3 and create a heatmap"
🤖 ✅ Found 23 significant markers (adj. p < 0.05)
    ✅ Top markers: GFAP, S100B, AQP4 (astrocyte signature)
    ✅ Generated expression heatmap
```
*👤 = You chatting with Claude | 🤖 = ChatSpatial MCP executing analysis*


---

## 🚀 Why Researchers Choose ChatSpatial

<table>
<tr>
<td width="50%" valign="top">

### Before: Traditional Analysis
```python
# 50+ lines of code for basic analysis
import scanpy as sc
import squidpy as sq
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read_h5ad("data.h5ad")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
# ... 40+ more lines ...
```
❌ Hours of coding  
❌ Complex syntax  
❌ Error-prone  
❌ Steep learning curve  

</td>
<td width="50%" valign="top">

### After: ChatSpatial
**Chat with Claude using natural language:**
```text
"Analyze my Visium data and find
 spatially variable genes"
```

✅ **5 seconds to results**
✅ **Plain English in Claude chat**
✅ **Zero programming required**
✅ **Publication-ready output**
✅ **60 methods via MCP**  

</td>
</tr>
</table>

---

## ⚡ Quick Start (2 Minutes)

### 1. Create Virtual Environment & Install
```bash
# Clone and enter directory
git clone https://github.com/cafferychen777/ChatSpatial.git
cd chatspatial

# Create virtual environment (strongly recommended)
# For macOS with Homebrew Python:
/opt/homebrew/bin/python3.10 -m venv chatspatial_env  # macOS Homebrew
# For other systems:
# python3 -m venv chatspatial_env  # Linux/other systems
source chatspatial_env/bin/activate  # macOS/Linux
# chatspatial_env\Scripts\activate  # Windows

# Verify Python version and install
python --version  # Should be 3.10+
pip install --upgrade pip
pip install -e ".[full]"  # Recommended: All features included
```

> 💡 **Windows Users:** SingleR and PETSc acceleration are not available on Windows due to C++ compilation limitations. Use alternative cell type annotation methods (Tangram, scANVI, CellAssign). All R-based methods (RCTD, SPOTlight, Numbat) work on Windows. See [INSTALLATION.md](INSTALLATION.md#windows) for details.

### 2. Configure Your Client

<details>
<summary><strong>Option A: Claude Desktop</strong> (GUI Application)</summary>

> 💡 **New to Claude Desktop?** [Download Claude Desktop](https://claude.com/download) from Anthropic's official site (available for Mac & Windows)

Add to `~/Library/Application Support/Claude/claude_desktop_config.json`:
```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/chatspatial_env/bin/python",
      "args": ["-m", "chatspatial", "server"]
    }
  }
}
```
</details>

<details>
<summary><strong>Option B: Claude Code</strong> (Terminal/IDE)</summary>

**Step 1: Install Claude Code CLI**
```bash
npm install -g @anthropic-ai/claude-code
```

**Step 2: Find Your Virtual Environment Path**
```bash
# First, activate your virtual environment
source chatspatial_env/bin/activate

# Then get the exact Python path (copy this output!)
which python
```
**Example output:** `/Users/yourname/Research/chatspatial_env/bin/python`

**Step 3: Add ChatSpatial MCP Server**
```bash
# Replace <PYTHON_PATH_FROM_STEP_2> with the actual path from step 2
claude mcp add chatspatial <PYTHON_PATH_FROM_STEP_2> -- -m chatspatial server

# Example with real path:
# claude mcp add chatspatial /Users/yourname/Research/chatspatial_env/bin/python -- -m chatspatial server
```

**Step 4: Verify Installation**
```bash
claude mcp list
```

**Step 5: Restart Claude Code** ⚠️
```bash
# IMPORTANT: Exit and restart Claude Code for changes to take effect
/quit

# Then re-enter Claude Code
claude
```
> 💡 **Note**: After adding the MCP server, you MUST restart Claude Code using `/quit` and then re-enter with `claude` for ChatSpatial to appear in the available tools.

</details>

### 3. Download Sample Data & Start Analyzing

**🎯 Quick Test (1 minute):**

1. **Download sample datasets**: Go to [**ChatSpatial Releases**](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data)
2. **Download these files**:
   - `card_reference_filtered.h5ad` (36MB - pancreatic reference with 9 cell types)
   - `card_spatial.h5ad` (7.7MB - spatial data with clear tissue domains)

**3. Open Claude Desktop or Claude Code and chat with ChatSpatial:**

> 💡 **Important**: Tell Claude to use ChatSpatial MCP for your spatial analysis (e.g., "Use ChatSpatial MCP for all my spatial transcriptomics analysis"). This ensures Claude calls the MCP tools instead of writing scripts.

Simply type this natural language request in your Claude chat:

```text
Load /Users/yourname/Downloads/card_reference_filtered.h5ad and /Users/yourname/Downloads/card_spatial.h5ad, then show me the tissue structure
```

> 💡 **How it works**: ChatSpatial MCP server interprets your natural language and automatically calls the appropriate analysis tools. No coding required!

> ⚠️ **IMPORTANT**: Use **absolute paths** when loading data (e.g., `/Users/yourname/Downloads/card_reference_filtered.h5ad`)

**🎯 That's it!** ChatSpatial will load your data and create a beautiful spatial visualization - all through natural language conversation.

> **📚 Detailed Setup Guides**: [Claude Desktop](INSTALLATION.md#claude-desktop) | [Claude Code](INSTALLATION.md#claude-code)

---

## 🧬 What You Can Do

> 💬 **All examples below are natural language commands you type in Claude chat with ChatSpatial MCP**

### 📊 **Try These Examples** (After Loading Sample Data)

Chat with Claude using these natural language requests:

```text
Identify spatial domains using SpaGCN

Deconvolve the spatial data using the reference with Cell2location

Analyze cell communication between spatial regions

Find spatially variable genes and create heatmaps
```

**What happens**: ChatSpatial MCP automatically selects the right tools, runs the analysis, and returns publication-ready visualizations.

### 🔍 **Spatial Analysis - Just Ask in Plain English**

| Your Question | What ChatSpatial Does |
|---------------|----------------------|
| "Find spatial domains" | Runs SpaGCN/STAGATE/Leiden clustering |
| "Detect hotspots" | Applies Getis-Ord Gi* spatial statistics |
| "Map cell territories" | Performs spatial neighborhood analysis |

### 🧮 **Advanced Methods - No Coding Knowledge Required**

| Your Request | Methods Used |
|-------------|--------------|
| "Deconvolve this spatial data" | Cell2location + scvi-tools |
| "Analyze cell communication" | LIANA+ + CellPhoneDB |
| "Find developmental trajectories" | CellRank + Palantir |
| "Run pathway enrichment" | GSEA + spatial smoothing |

### 🎨 **Instant Visualizations**
- **Spatial plots** with tissue overlays
- **Interactive heatmaps** for gene expression
- **Communication networks** between cell types
- **Trajectory flow maps** for development
- **Domain boundary visualizations**

---

## 🎯 Choose Your Path

<table>
<tr>
<td width="33%" align="center">

### 🚀 **Researchers**
**Quick start?**

```bash
pip install -e .
```
✅ 80% of features  
✅ Most common methods  
✅ 6-minute install  

**→ [Research Quick Start](docs/getting-started/quick-start.md)**

</td>
<td width="33%" align="center">

### 🧠 **Power Users**
**Want everything?**

```bash
pip install -e ".[full]"
```
✅ 100% of features  
✅ All 16+ methods  
✅ Deep learning included  

**→ [Advanced Setup Guide](docs/tutorials/learning-paths/advanced.md)**

</td>
<td width="33%" align="center">

### 👩‍💻 **Developers**
**Want to contribute?**

```bash
pip install -e ".[dev]"
```
✅ Development tools  
✅ Testing framework  
✅ Documentation  

**→ [Contributor Guide](CONTRIBUTING.md)**

</td>
</tr>
</table>

---

## 🛠️ Technical Capabilities

<details>
<summary><strong>📊 Data Formats Supported</strong></summary>

- **10x Genomics**: Visium, Xenium
- **Spatial Technologies**: Slide-seq v2  
- **Multiplexed Imaging**: MERFISH, seqFISH
- **Standard Formats**: H5AD, H5, MTX, CSV

</details>

<details>
<summary><strong>🔬 Analysis Methods (12 Categories, 75+ Methods)</strong></summary>

| Category | Methods |
|----------|---------|
| **Cell Type Annotation** | Tangram, scANVI, CellAssign, mLLMCellType, sc-type, SingleR |
| **Spatial Domains** | SpaGCN, STAGATE, Leiden clustering |
| **Cell Communication** | LIANA+, CellPhoneDB, CellChat (via LIANA) |
| **Deconvolution** | Cell2location, DestVI, RCTD, Tangram, Stereoscope, SPOTlight |
| **CNV Analysis** | infercnvpy, Numbat (haplotype-aware CNV analysis) |
| **Spatial Variable Genes** | SpatialDE, SPARK-X |
| **Trajectory & Velocity** | CellRank, Palantir, DPT, scVelo, VeloVI |
| **Sample Integration** | Harmony, BBKNN, Scanorama, scVI |
| **Differential Expression** | Wilcoxon, t-test, Logistic Regression (scanpy methods) |
| **Gene Set Enrichment** | GSEA, ORA, ssGSEA, Enrichr, Spatial EnrichMap |
| **Spatial Statistics** | Moran's I, Local Moran's I (LISA), Geary's C, Getis-Ord Gi*, Neighborhood Enrichment, Co-occurrence, Ripley's K/L, Bivariate Moran's I, Join Count, Network Properties, Spatial Centrality |
| **Spatial Registration** | PASTE, STalign |

</details>

<details>
<summary><strong>⚙️ System Requirements</strong></summary>

- **Python**: 3.10+ (required for MCP)
- **Memory**: 8GB+ RAM (16GB+ for large datasets)  
- **Storage**: 5GB+ for dependencies
- **OS**: Linux, macOS, Windows (WSL recommended)
- **GPU**: Optional (speeds up deep learning methods)

</details>

---

## 🤝 Community & Support

- **📧 Issues**: [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues) for bug reports  
- **📖 Docs**: [Complete documentation](docs/) with tutorials
- **⭐ Star this repo** to follow development

---

## 📄 License & Citation

**MIT License** - Free for academic and commercial use.

If ChatSpatial helps your research, please cite:
```bibtex
@software{chatspatial2025,
  title={ChatSpatial: An Agentic Framework for Reproducible Cross-Platform Spatial Transcriptomics Analysis},
  author={Chen Yang and Xianyang Zhang and Jun Chen},
  year={2025},
  url={https://github.com/cafferychen777/ChatSpatial},
  note={Manuscript in preparation}
}
```

## 🙏 Built With

ChatSpatial stands on the shoulders of giants:
[**Scanpy**](https://scanpy.readthedocs.io/) • [**Squidpy**](https://squidpy.readthedocs.io/) • [**scvi-tools**](https://scvi-tools.org/) • [**LIANA**](https://liana-py.readthedocs.io/) • [**Anthropic MCP**](https://modelcontextprotocol.io/)

---

<div align="center">

**Made with ❤️ for the spatial transcriptomics community**

[⭐ **Star us on GitHub**](https://github.com/cafferychen777/ChatSpatial) if this project helps you!

</div>