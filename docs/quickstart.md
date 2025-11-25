# Quick Start (5 Minutes)

Get ChatSpatial running and analyze your first spatial dataset.

---

## Step 1: Install (2 minutes)

```bash
# Clone and create virtual environment
git clone https://github.com/cafferychen777/ChatSpatial.git
cd chatspatial
python3 -m venv chatspatial_env
source chatspatial_env/bin/activate  # macOS/Linux
# chatspatial_env\Scripts\activate  # Windows

# Install ChatSpatial
pip install -e ".[full]"  # All features (recommended)
# pip install -e .         # Core features only (faster)
```

**That's it for installation!** ✅

---

## Step 2: Configure Claude (1 minute)

### Find Your Python Path
```bash
# In your activated virtual environment
which python
# Copy this path (e.g., /Users/yourname/chatspatial_env/bin/python)
```

### Configure Your Client

**For Claude Desktop:**
Edit `~/Library/Application Support/Claude/claude_desktop_config.json`:
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

**For Claude Code:**
```bash
claude mcp add chatspatial /path/to/chatspatial_env/bin/python -- -m chatspatial server
```

Restart your client for changes to take effect.

---

## Step 3: Try It! (2 minutes)

### Download Sample Data
Get example datasets from [ChatSpatial Releases](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data):
- `card_spatial.h5ad` (7.7MB - pancreatic spatial data)
- `card_reference_filtered.h5ad` (36MB - reference with 9 cell types)

### Chat with Claude

Open Claude and start chatting naturally:

```text
👤 "Load /Users/yourname/Downloads/card_spatial.h5ad and show me what's in it"

🤖 ✅ Loaded 428 spots, 1,926 genes
    Tissue: Mouse pancreas
    Platform: Spatial transcriptomics

👤 "Identify spatial domains"

🤖 ✅ Found 7 distinct spatial domains
    ✅ Generated visualization

👤 "Find marker genes for domain 3"

🤖 ✅ Found 15 significant markers
    Top genes: Ins1, Ins2, Gcg (pancreatic islet signature)
```

**That's it! You're analyzing spatial data through conversation.** 🎉

---

## What to Try Next

**Basic Analysis:**
```text
"Show me gene expression patterns"
"Cluster the spots and visualize them"
"Find spatially variable genes"
```

**Advanced Analysis:**
```text
"Annotate cell types using the reference data"
"Analyze cell-cell communication patterns"
"Deconvolve spatial data with Cell2location"
```

**Custom Requests:**
```text
"Find marker genes enriched in region X"
"Compare gene expression between domains"
"Create a heatmap of top variable genes"
```

---

## Common Issues

**Problem: ChatSpatial tools not appearing**
- ✅ Restart Claude after configuration
- ✅ Verify Python path in config is correct
- ✅ Check virtual environment is activated

**Problem: Cannot load data**
- ✅ Use absolute paths (e.g., `/Users/name/data.h5ad`)
- ✅ Verify file exists at the path
- ✅ Check file format is supported (H5AD, H5, MTX, Visium)

**Problem: Analysis fails**
- ✅ Make sure you ran preprocessing first
- ✅ Check data quality (minimum cells/genes)
- ✅ See [Troubleshooting Guide](advanced/troubleshooting.md) for details

---

## Next Steps

- **[Examples](examples.md)** - More analysis examples and workflows
- **[Methods Reference](advanced/methods-reference.md)** - Complete list of all tools and parameters
- **[FAQ](advanced/faq.md)** - Frequently asked questions

**Need help?** [Open an issue](https://github.com/cafferychen777/ChatSpatial/issues)
