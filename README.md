# ChatSpatial 🧬 

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![MCP Protocol](https://img.shields.io/badge/MCP-v2024.11.05-green.svg)](https://modelcontextprotocol.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Tests](https://img.shields.io/badge/tests-passing-green.svg)]()

**Interactive Spatial Transcriptomics Analysis via Model Context Protocol**

ChatSpatial is a production-ready **Model Context Protocol (MCP) server** that provides AI assistants with comprehensive spatial transcriptomics analysis capabilities. It enables natural language interaction with complex spatial data analysis through 16 standardized MCP tools.

## 🎯 Why ChatSpatial?

- **🔗 Universal AI Integration**: Works seamlessly with Claude Desktop, Cherry Studio, Continue, and any MCP-compatible AI tool
- **🧬 Spatial-First Design**: Purpose-built for spatial transcriptomics (10x Visium, Slide-seq, MERFISH, seqFISH)
- **⚡ Production Ready**: All core methods tested and validated with real-world datasets
- **🎛️ Comprehensive**: 16 tools covering the entire spatial analysis workflow
- **🛡️ Robust**: Advanced error handling and compatibility management

## 🚀 Quick Start

### Prerequisites
- Python 3.8+ (Python 3.10-3.11 recommended for best compatibility)
- An MCP-compatible AI assistant (Claude Desktop, Cherry Studio, etc.)

### Installation
```bash
# 1. Create a dedicated environment (highly recommended)
conda create -n chatspatial python=3.10
conda activate chatspatial

# 2. Install ChatSpatial
git clone https://github.com/cafferychen777/ChatSpatial.git
cd ChatSpatial
pip install -e .

# 3. Verify installation
chatspatial --help
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
```
Load my 10x Visium dataset from /path/to/data.h5ad
Preprocess data with standard filtering
Annotate cell types using marker genes
Visualize spatial cell type distribution
```

For detailed installation instructions, see **[docs/INSTALLATION.md](docs/INSTALLATION.md)**.

## 🛠️ Core Capabilities

### 📊 **Data Management & Preprocessing**
- **Multi-format Loading**: 10x Visium, Slide-seq, MERFISH, seqFISH, H5AD
- **Intelligent Preprocessing**: QC, normalization, dimensionality reduction with smart defaults
- **Quality Control**: Comprehensive filtering and validation

### 🔬 **Cell Analysis**
- **Cell Type Annotation**: Marker-based, Tangram, scANVI, CellAssign, mLLMCellType, sc-type
- **Differential Expression**: Advanced marker discovery between cell populations
- **Data Integration**: Multi-sample integration (Harmony, scVI, BBKNN)

### 🧬 **Spatial Analysis**
- **Spatial Variable Genes**: GASTON (deep learning), SpatialDE, SPARK methods
- **Spatial Domains**: SpaGCN, STAGATE, BANKSY clustering
- **Spatial Statistics**: Moran's I, Geary's C, Getis-Ord Gi*, spatial autocorrelation

### 💬 **Cell Communication**
- **LIANA+**: Fast spatial bivariate analysis (cosine, pearson, spearman, jaccard)
- **CellPhoneDB v3**: Statistical permutation testing with spatial microenvironments
- **CellChat v2**: Advanced pattern recognition via LIANA integration

### 🧮 **Advanced Methods**
- **Spatial Deconvolution**: DestVI, Cell2location, RCTD, SPOTlight
- **Trajectory Analysis**: Palantir, CellRank, DPT pseudotime inference
- **RNA Velocity**: RNA velocity analysis with spatial context
- **Pathway Enrichment**: GSEA, ORA, Enrichr with spatial smoothing

### 📈 **Visualization**
- **15+ Plot Types**: Spatial, UMAP, violin, heatmap, trajectory, communication plots
- **MCP Image Objects**: Seamless display in AI assistants
- **Interactive Support**: Plotly and Bokeh integration

## 🏗️ Architecture

ChatSpatial implements a clean **Model Context Protocol** architecture:

```
AI Assistant → MCP Client → ChatSpatial Server → Analysis Tools → Results
                ↓
        MCP ToolResult ← Visualization ← Data Processing
```

**Key Features:**
- **Protocol**: MCP v2024-11-05 compliance
- **Transport**: stdio (standard input/output) and SSE (Server-Sent Events)
- **Tools**: 32 spatial transcriptomics analysis tools
- **Resources**: Automatic dataset management via `spatial://` URIs
- **Error Handling**: Robust two-layer error management system

## 🔧 Tool Categories

| Category | Tools | Examples |
|----------|-------|----------|
| **Data** | `load_data`, `preprocess_data` | Load 10x Visium, quality control |
| **Visualization** | `visualize_data` | Spatial plots, UMAP, heatmaps |
| **Cell Analysis** | `annotate_cells`, `find_markers` | Cell typing, differential expression |
| **Spatial** | `find_spatial_genes`, `identify_spatial_domains` | GASTON, SpaGCN, spatial statistics |
| **Communication** | `analyze_cell_communication` | LIANA+, CellPhoneDB, ligand-receptor |
| **Integration** | `integrate_samples`, `analyze_trajectory_data` | Harmony, pseudotime |
| **Advanced** | `deconvolve_data`, `analyze_enrichment` | DestVI, GSEA |

## 🌟 Example Workflows

### Basic Spatial Analysis
```
1. "Load my 10x Visium dataset from /path/to/data.h5ad"
2. "Preprocess with genes in ≥10 cells and cells with ≥500 genes"
3. "Identify spatial domains using STAGATE"
4. "Visualize spatial domains"
```

### Cell Communication Analysis
```
1. "Annotate cell types using marker genes"
2. "Analyze cell communication using LIANA+ with cosine similarity"
3. "Visualize communication for VEGFA-KDR interaction"
```

### Advanced Deep Learning
```
1. "Find spatial variable genes using GASTON with GLM-PCA"
2. "Deconvolve spatial data using DestVI"
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

See [docs/INSTALLATION.md](docs/INSTALLATION.md) for detailed client setup.

## 📦 Optional Dependencies

ChatSpatial uses a modular dependency system. Install specific modules as needed:

```bash
# Full installation (recommended)
pip install -e .[all]

# Or install specific modules
pip install -e .[enrichmap]        # EnrichMap spatial enrichment
pip install -e .[deconvolution]    # DestVI, Cell2location
pip install -e .[spatial_domains]  # SpaGCN, STAGATE, BANKSY  
pip install -e .[spatial_genes]    # GASTON, SpatialDE, SPARK
pip install -e .[cell_communication] # LIANA+, CellPhoneDB
```

## 📋 System Requirements

| Component | Requirement |
|-----------|-------------|
| **Python** | 3.8+ (3.10-3.11 recommended) |
| **Memory** | 8GB+ RAM (16GB+ for large datasets) |
| **Storage** | 5GB+ for dependencies |
| **OS** | Linux, macOS, Windows (with WSL recommended) |

## 🧪 Production Status

ChatSpatial is production-ready with comprehensive testing:

- ✅ **16 MCP Tools**: All core tools tested and validated
- ✅ **Error Handling**: Robust two-layer error management
- ✅ **Visualization**: All plot types verified with MCP Image objects
- ✅ **Compatibility**: Tested with Claude Desktop, Cherry Studio
- ✅ **Real Data**: Validated with 10x Visium, Slide-seq, MERFISH datasets

## 🔍 Troubleshooting

**Common Issues:**
- **Import Errors**: Check optional dependencies for specific methods
- **Memory Issues**: Use data subsampling for large datasets
- **Long Processing**: Use Cherry Studio with increased timeout
- **Visualization**: Ensure MCP client supports Image objects

See [docs/user_guides/ERROR_HANDLING_GUIDE.md](docs/user_guides/ERROR_HANDLING_GUIDE.md) for detailed troubleshooting.

## 📚 Documentation

| Document | Description |
|----------|-------------|
| **[INSTALLATION.md](docs/INSTALLATION.md)** | Detailed installation and setup guide |
| **[Error Handling Guide](docs/user_guides/ERROR_HANDLING_GUIDE.md)** | Troubleshooting and error resolution |
| **[MCP Tools Reference](docs/technical_docs/MCP_TOOLS_QUICK_REFERENCE.md)** | Complete tool documentation |
| **[Module Documentation](docs/modules/)** | API reference for all modules |

## 🤝 Contributing

We welcome contributions! Please see:
- [CONTRIBUTING.md](CONTRIBUTING.md) - Contribution guidelines
- [SECURITY.md](SECURITY.md) - Security policy
- [docs/technical_docs/](docs/technical_docs/) - Technical documentation

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

**Ready to analyze spatial data with AI?** Install ChatSpatial and start exploring! 🚀