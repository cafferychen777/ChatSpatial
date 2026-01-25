<div align="center">

# ChatSpatial

**MCP server for spatial transcriptomics analysis via natural language**

[![PyPI](https://img.shields.io/pypi/v/chatspatial)](https://pypi.org/project/chatspatial/)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Docs](https://img.shields.io/badge/docs-available-blue)](https://cafferychen777.github.io/ChatSpatial/)

</div>

---

<table>
<tr>
<td width="50%">

### ❌ Before
```python
import scanpy as sc
import squidpy as sq
adata = sc.read_h5ad("data.h5ad")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
# ... 40 more lines
```

</td>
<td width="50%">

### ✅ After
```text
"Load my Visium data and identify
 spatial domains"
```

```
✓ Loaded 3,456 spots, 18,078 genes
✓ Identified 7 spatial domains
✓ Generated visualization
```

</td>
</tr>
</table>

---

> **Note**: This is a **demo version**. If you encounter any issues or have feedback, please [open an issue](https://github.com/cafferychen777/ChatSpatial/issues) or contact us anytime. Your feedback helps us improve!

---

## Install

```bash
pip install chatspatial
```

## Configure

> **Important:** Use your virtual environment's Python path. Run `which python` to find it.

**Claude Code** (Recommended):

```bash
claude mcp add chatspatial /path/to/venv/bin/python -- -m chatspatial server
```

**Codex** (CLI or IDE extension):

```bash
codex mcp add chatspatial -- /path/to/venv/bin/python -m chatspatial server
```

Or edit `~/.codex/config.toml`:

```toml
[mcp_servers.chatspatial]
command = "/path/to/venv/bin/python"
args = ["-m", "chatspatial", "server"]
```

**Claude Desktop** — add to config file:

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/venv/bin/python",
      "args": ["-m", "chatspatial", "server"]
    }
  }
}
```

> Config locations: macOS `~/Library/Application Support/Claude/claude_desktop_config.json` · Windows `%APPDATA%\Claude\claude_desktop_config.json` · Linux `~/.config/Claude/claude_desktop_config.json`

Restart your client after configuration.


---

## Use

Open Claude and chat:

```text
Load /path/to/spatial_data.h5ad and show me the tissue structure
```

```text
Identify spatial domains using SpaGCN
```

```text
Find spatially variable genes and create a heatmap
```

---

## Capabilities

| Category | Methods |
|----------|---------|
| **Spatial Domains** | SpaGCN, STAGATE, GraphST, Leiden, Louvain |
| **Deconvolution** | FlashDeconv, Cell2location, RCTD, DestVI, Stereoscope, SPOTlight, Tangram, CARD |
| **Cell Communication** | LIANA+, CellPhoneDB, CellChat, FastCCC |
| **Cell Type Annotation** | Tangram, scANVI, CellAssign, mLLMCelltype, scType, SingleR |
| **Trajectory & Velocity** | CellRank, Palantir, DPT, scVelo, VeloVI |
| **Spatial Statistics** | Moran's I, Local Moran, Geary's C, Getis-Ord Gi*, Ripley's K, Neighborhood Enrichment |
| **Enrichment** | GSEA, ORA, Enrichr, ssGSEA, Spatial EnrichMap |
| **Spatial Genes** | SpatialDE, SPARK-X |
| **Integration** | Harmony, BBKNN, Scanorama, scVI |
| **Other** | CNV Analysis (inferCNVpy, Numbat), Spatial Registration (PASTE, STalign) |

**60+ methods** across 15 categories. **Supports** 10x Visium, Xenium, Slide-seq v2, MERFISH, seqFISH.

---

## Docs

- [Installation Guide](INSTALLATION.md) — detailed setup for all platforms
- [Examples](docs/examples.md) — step-by-step workflows
- [Methods Reference](docs/advanced/methods-reference.md) — all 20 tools documented
- [Full Documentation](https://cafferychen777.github.io/ChatSpatial/) — complete reference

---

## Citation

```bibtex
@software{chatspatial2025,
  title={ChatSpatial: Agentic Workflow for Spatial Transcriptomics},
  author={Chen Yang and Xianyang Zhang and Jun Chen},
  year={2025},
  url={https://github.com/cafferychen777/ChatSpatial}
}
```

<div align="center">

**MIT License** · [GitHub](https://github.com/cafferychen777/ChatSpatial) · [Issues](https://github.com/cafferychen777/ChatSpatial/issues)

</div>

<!-- mcp-name: io.github.cafferychen777/chatspatial -->
