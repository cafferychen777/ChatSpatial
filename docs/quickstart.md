# Quick Start

Get started in 5 minutes. By the end, you'll have loaded data and created your first spatial visualization.

---

## 1. Install

```bash
pip install chatspatial
```

Requires Python 3.11+ ([full installation options](advanced/installation.md))

---

## 2. Configure

Choose your client. **Important:** Use your virtual environment's Python path.

```bash
# Find your Python path first
which python  # e.g., /Users/you/chatspatial_env/bin/python
```

**Claude Code** (Recommended)

```bash
claude mcp add chatspatial /path/to/chatspatial_env/bin/python -- -m chatspatial server
```

**Codex** (CLI or IDE)

```bash
codex mcp add chatspatial -- /path/to/chatspatial_env/bin/python -m chatspatial server
```

Or edit `~/.codex/config.toml`:

```toml
[mcp_servers.chatspatial]
command = "/path/to/chatspatial_env/bin/python"
args = ["-m", "chatspatial", "server"]
```

**Claude Desktop**

Add to config file (`~/Library/Application Support/Claude/claude_desktop_config.json`):

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

Restart your client after configuring. See [Configuration Guide](advanced/configuration.md) for detailed setup.

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

## What Success Looks Like

When ChatSpatial works correctly, you'll see:

**After loading data:**
```
Dataset loaded successfully
- ID: spatial_data_abc123
- 3000 spots, 18000 genes
- Platform: Visium
- Spatial coordinates: available
```

**After preprocessing:**
```
Preprocessing complete
- Filtered to 2800 spots, 2000 HVGs
- Normalized with log1p
- Computed 30 PCs, UMAP, 15 neighbors
- Found 8 clusters (Leiden)
```

**After visualization:**

ChatSpatial returns images directly in the chat. You'll see spatial plots showing gene expression or clusters overlaid on tissue coordinates.

---

## Next Steps

- [Concepts](concepts.md) - Understand what each analysis does
- [Examples](examples.md) - Complete workflows for every analysis type
- [Methods Reference](advanced/methods-reference.md) - All 20 tools with parameters
