# ChatSpatial 🧬

[Docs](docs/index.md) | [Tutorials](docs/tutorials/) | [API](docs/api/) | [Code of Conduct](CODE_OF_CONDUCT.md)

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![MCP Protocol](https://img.shields.io/badge/MCP-v2024.11.05-green.svg)](https://modelcontextprotocol.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/cafferychen777/ChatSpatial/actions/workflows/ci.yml/badge.svg)](https://github.com/cafferychen777/ChatSpatial/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-available-blue)](https://cafferychen777.github.io/ChatSpatial/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Code of Conduct](https://img.shields.io/badge/code%20of%20conduct-Contributor%20Covenant-ff69b4)](CODE_OF_CONDUCT.md)

## Interactive Spatial Transcriptomics Analysis via Model Context Protocol

ChatSpatial is a production-ready **Model Context Protocol (MCP) server** that provides LLM agents with comprehensive spatial transcriptomics analysis capabilities. It enables natural language interaction with complex spatial data analysis through 16 standardized MCP tools.

> **📁 Data Note**: Datasets are not included in the repository. Use the provided download scripts in `data/scripts/` or see [FINAL_MCP_DATASETS_REPORT.md](FINAL_MCP_DATASETS_REPORT.md) for dataset acquisition instructions.

## 🎯 Why ChatSpatial?

- **🔗 Universal Agent Integration**: Works seamlessly with Claude Desktop, Cherry Studio, Continue, and any MCP-compatible LLM agent
- **🧬 Spatial-First Design**: Purpose-built for spatial transcriptomics (10x Visium, Slide-seq, MERFISH, seqFISH)
- **⚡ Production Ready**: All core methods tested and validated with real-world datasets
- **🎛️ Comprehensive**: 16 tools covering the entire spatial analysis workflow
- **🛡️ Robust**: Advanced error handling and compatibility management

## 🚀 Quick Start

### Prerequisites

- Python 3.8+ (Python 3.10-3.11 recommended for best compatibility)
- An MCP-compatible LLM agent client (Claude Desktop, Cherry Studio, etc.)

### Installation

```bash
# 1. Clone repository
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial

# 2. Install ChatSpatial (core features)
pip install -e .

# 3. Or install with advanced features
pip install -e ".[advanced]"

# 4. Verify installation
python -c "import chatspatial; print('✅ Installation successful')"
```

### MCP Setup

```bash
# Find your Python path
which python
# Example: /opt/anaconda3/envs/chatspatial/bin/python
```

Add to your MCP client configuration:

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/your/chatspatial/bin/python",
      "args": ["-m", "chatspatial"],
      "env": {}
    }
  }
}
```

### First Analysis

```text
Load my 10x Visium dataset from /path/to/data.h5ad
Preprocess data with standard filtering
Annotate cell types using marker genes
Visualize spatial cell type distribution
```

For detailed installation instructions, see **[INSTALLATION.md](INSTALLATION.md)**.

## 🛠️ Core Capabilities

### 📊 **Data Management & Preprocessing**

- **Multi-format Loading**: 10x Visium, Slide-seq, MERFISH, seqFISH, H5AD
- **Intelligent Preprocessing**: QC, normalization, dimensionality reduction with smart defaults
- **Quality Control**: Comprehensive filtering and validation


### 🔬 **Cell Analysis**

- **Cell Type Annotation**: Marker-based, Tangram, scANVI, CellAssign, mLLMCellType, scType
- **Differential Expression**: Advanced marker discovery between cell populations
- **Data Integration**: Multi-sample integration (Harmony, scVI, BBKNN)

### 🧬 **Spatial Analysis**

- **Spatial Variable Genes**: GASTON (deep learning), SpatialDE, SPARK-X (non-parametric) methods
- **Spatial Domains**: SpaGCN, STAGATE, BANKSY, Leiden/Louvain clustering
- **Spatial Statistics**: Moran's I, Geary's C, Getis-Ord Gi*, spatial autocorrelation

### 💬 **Cell Communication**

- **LIANA+**: Fast spatial bivariate analysis (cosine, pearson, spearman, jaccard)
- **CellPhoneDB**: Statistical permutation testing with spatial microenvironments
- **CellChat via LIANA**: Advanced pattern recognition via LIANA integration

### 🧮 **Advanced Methods**

- **Spatial Deconvolution**: Cell2location, DestVI, RCTD, Stereoscope, Tangram, SPOTlight
- **Trajectory Analysis**: Palantir, CellRank, DPT pseudotime inference
- **RNA Velocity**: RNA velocity analysis with spatial context
- **Pathway Enrichment**: GSEA, ORA, Enrichr with spatial smoothing

### 📈 **Visualization**

- **20 Plot Types**: Spatial, UMAP, violin, heatmap, trajectory, communication plots
- **MCP Image Objects**: Seamless display in LLM agent clients
- **Interactive Support**: Plotly and Bokeh integration

## 🏗️ Architecture

ChatSpatial implements a clean **Model Context Protocol** architecture:

```text
LLM Agent Client → MCP Protocol → ChatSpatial Server → Analysis Tools → Results
                     ↓
           MCP ToolResult ← Visualization ← Data Processing
```

**Key Features:**
- **Protocol**: MCP v2024-11-05 compliance
- **Transport**: stdio (standard input/output) and SSE (Server-Sent Events)
- **Tools**: 16 spatial transcriptomics analysis tools
- **Resources**: Automatic dataset management via `spatial://` URIs
- **Error Handling**: Robust two-layer error management system

## 🔧 Tool Categories

| Category | Tools | Examples |
|----------|-------|----------|
| **Data** | `load_data`, `preprocess_data` | Load 10x Visium, quality control |
| **Visualization** | `visualize_data` | Spatial plots, UMAP, heatmaps |
| **Cell Analysis** | `annotate_cells`, `find_markers` | Cell typing, differential expression |
| **Spatial** | `find_spatial_genes`, `identify_spatial_domains` | GASTON, SpaGCN, spatial statistics |
| **Communication** | `analyze_cell_communication` | LIANA, CellPhoneDB, ligand-receptor |
| **Integration** | `integrate_samples`, `analyze_trajectory_data` | Harmony, pseudotime |
| **Advanced** | `deconvolve_data`, `analyze_enrichment` | Cell2location, GSEA |

## 🌟 Example Workflows

### Basic Spatial Analysis

```text
1. "Load my 10x Visium dataset from /path/to/data.h5ad"
2. "Preprocess with genes in ≥10 cells and cells with ≥500 genes"
3. "Identify spatial domains using SpaGCN"
4. "Visualize spatial domains"
```

### Cell Communication Analysis

```text
1. "Annotate cell types using marker genes"
2. "Analyze cell communication using LIANA with cosine similarity"
3. "Visualize communication for VEGFA-KDR interaction"
```

### Advanced Deep Learning

```text
1. "Find spatial variable genes using GASTON with GLM-PCA"
2. "Deconvolve spatial data using Cell2location"
3. "Visualize deconvolution results and GASTON isodepth map"
```

## 🔌 MCP Client Integration

### Claude Desktop
1. Download [Claude Desktop](https://claude.ai/desktop)
2. Edit configuration: Settings → Developer → Edit Config
3. Add ChatSpatial server configuration
4. Restart Claude Desktop
5. Look for hammer icon to access tools

**Recommended**: Use **Cherry Studio** for computationally intensive tasks due to configurable timeouts.

### Cherry Studio (Recommended for Heavy Analysis)
- **Configurable Timeout**: Set to 3600s for long-running analysis
- **Stable Processing**: No interruptions during complex tasks
- **Better Performance**: Ideal for GASTON, deconvolution, large datasets

See [INSTALLATION.md](INSTALLATION.md) for detailed client setup.

## 📦 Installation Options

ChatSpatial offers two installation levels:

```bash
# Core installation (recommended for most users)
pip install -e .

# Advanced installation (deep learning, specialized methods)
pip install -e ".[advanced]"

# Development installation (for contributors)
pip install -e ".[dev]"
```

## 📋 System Requirements

| Component | Requirement |
|-----------|-------------|
| **Python** | 3.8+ (3.10+ recommended) |
| **Memory** | 8GB+ RAM (16GB+ for large datasets) |
| **Storage** | 5GB+ for dependencies |
| **OS** | Linux, macOS, Windows (with WSL recommended) |

## 🧪 Production Status

ChatSpatial is production-ready with comprehensive testing:

- ✅ **16 MCP Tools**: All core tools tested and validated
- ✅ **Error Handling**: Robust two-layer error management
- ✅ **Visualization**: All plot types verified with MCP Image objects
- ✅ **Compatibility**: Tested with Claude Desktop, Cherry Studio LLM agent clients
- ✅ **Real Data**: Validated with 10x Visium, Slide-seq, MERFISH datasets

## 🔍 Troubleshooting

**Common Issues:**
- **Import Errors**: Check optional dependencies for specific methods
- **Memory Issues**: Use data subsampling for large datasets
- **Long Processing**: Use Cherry Studio with increased timeout
- **Visualization**: Ensure MCP client supports Image objects

See [UNIFIED_ERROR_HANDLING_MIGRATION_GUIDE.md](UNIFIED_ERROR_HANDLING_MIGRATION_GUIDE.md) for detailed troubleshooting.

## 📚 Documentation

- Docs Site: [docs/index.md](docs/index.md)
- Getting Started: [docs/getting_started.md](docs/getting_started.md)
- Tutorials: [docs/tutorials/](docs/tutorials/)
- API Reference: [docs/api/](docs/api/)

Additional guides:
- **[INSTALLATION.md](INSTALLATION.md)**
- **[Error Handling Guide](UNIFIED_ERROR_HANDLING_MIGRATION_GUIDE.md)**
- **[Dataset Guide](FINAL_MCP_DATASETS_REPORT.md)**
- **[Project Structure](PROJECT_STRUCTURE.md)**

## 🤝 Contributing

We welcome contributions! Please see:
- [CONTRIBUTING.md](CONTRIBUTING.md) - Contribution guidelines
- [SECURITY.md](SECURITY.md) - Security policy
- [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) - Technical documentation

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

ChatSpatial integrates and builds upon many excellent open-source projects:
- **MCP**: Model Context Protocol by Anthropic
- **Scanpy/Squidpy**: Single-cell and spatial analysis ecosystem
- **scvi-tools**: Deep learning for single-cell genomics
- **LIANA**: Ligand-receptor interaction analysis
- **GASTON**: Graph attention for spatial transcriptomics

## 📊 Citation

If you use ChatSpatial in your research, please cite:

```bibtex
@software{chatspatial2024,
  title={ChatSpatial: Interactive Spatial Transcriptomics Analysis via Model Context Protocol},
  author={ChatSpatial Development Team},
  year={2024},
  url={https://github.com/cafferychen777/ChatSpatial}
}
```

---

**Ready to analyze spatial data with LLM agent clients?** Install ChatSpatial and start exploring! 🚀