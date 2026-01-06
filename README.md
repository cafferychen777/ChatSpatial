<div align="center">

# ChatSpatial

**MCP server for spatial transcriptomics analysis via natural language**

[![PyPI](https://img.shields.io/pypi/v/chatspatial)](https://pypi.org/project/chatspatial/)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
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

## Install

```bash
pip install chatspatial
```

## Configure

**Claude Desktop** — add to `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "python",
      "args": ["-m", "chatspatial", "server"]
    }
  }
}
```

**Claude Code**:

```bash
claude mcp add chatspatial python -- -m chatspatial server
```

> Restart Claude after configuration.

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
| **Spatial Domains** | SpaGCN, STAGATE, Leiden |
| **Deconvolution** | Cell2location, RCTD, Tangram, SPOTlight |
| **Cell Communication** | LIANA+, CellPhoneDB |
| **Cell Type Annotation** | Tangram, scANVI, CellAssign, mLLMCelltype |
| **Trajectory** | CellRank, Palantir, scVelo |
| **Spatial Statistics** | Moran's I, Getis-Ord Gi*, Ripley's K |
| **Gene Set Enrichment** | GSEA, ORA, Enrichr |

**60+ methods** across 12 categories. **Supports** 10x Visium, Xenium, Slide-seq v2, MERFISH, seqFISH.

---

## Docs

- [Installation Guide](INSTALLATION.md) — detailed setup for all platforms
- [Examples](docs/examples/) — step-by-step workflows
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
