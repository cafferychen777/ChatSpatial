# Quick Start

Get started in 3 steps.

---

## 1. Install

```bash
pip install chatspatial
```

Requires Python 3.11+ ([full installation options](advanced/installation.md))

---

## 2. Configure

Choose your client:

**Claude Desktop**

Add to config file (`~/Library/Application Support/Claude/claude_desktop_config.json`):

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

**Claude Code**

```bash
claude mcp add chatspatial python -- -m chatspatial server
```

Restart Claude after configuring.

---

## 3. Start Analyzing

Open Claude and type:

```
Load /path/to/spatial_data.h5ad and show the tissue structure
```

ChatSpatial handles the rest.

---

## What You Can Do

| Task | Example |
|------|---------|
| Load data | "Load my Visium data from /path/to/data" |
| Preprocess | "Normalize and cluster the data" |
| Find domains | "Identify spatial domains with SpaGCN" |
| Annotate cells | "Annotate cell types using the reference" |
| Deconvolve | "Estimate cell type proportions" |
| Analyze patterns | "Find spatially variable genes" |
| Visualize | "Show expression of CD3D on the tissue" |

---

## Sample Data

Test with our datasets:

- [card_spatial.h5ad](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data) (7.7 MB) - Pancreatic spatial data
- [card_reference_filtered.h5ad](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data) (36 MB) - Reference for deconvolution

---

## Troubleshooting

| Issue | Fix |
|-------|-----|
| Tools not showing | Restart Claude |
| "Dataset not found" | Use absolute path: `/Users/you/data.h5ad` |
| Analysis fails | Run preprocessing first |

More help: [Troubleshooting Guide](advanced/troubleshooting.md)

---

## Next Steps

- [Examples](examples.md) - Complete workflows
- [Methods Reference](advanced/methods-reference.md) - All 20 tools
