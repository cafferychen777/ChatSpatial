# Quick Start

Get ChatSpatial running in 5 minutes.

---

## 1. Install

```bash
pip install chatspatial
```

> For development or all optional features, see [Installation Guide](advanced/installation.md).

---

## 2. Configure

**Claude Desktop** — edit `~/Library/Application Support/Claude/claude_desktop_config.json`:

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

Restart Claude after configuration.

---

## 3. Use

Open Claude and chat:

```text
Load /path/to/spatial_data.h5ad and show me the tissue structure
```

```text
Identify spatial domains using SpaGCN
```

```text
Find marker genes for domain 3
```

That's it. You're analyzing spatial data through conversation.

---

## Sample Data

Download test datasets from [Releases](https://github.com/cafferychen777/ChatSpatial/releases/tag/v0.3.0-data):

- `card_spatial.h5ad` (7.7MB) — pancreatic spatial data
- `card_reference_filtered.h5ad` (36MB) — reference for deconvolution

---

## Common Issues

| Problem | Solution |
|---------|----------|
| Tools not appearing | Restart Claude after config |
| Cannot load data | Use absolute paths |
| Analysis fails | Run preprocessing first |

See [Troubleshooting](advanced/troubleshooting.md) for more.

---

## Next

- [Examples](examples.md) — Analysis workflows for every use case
- [Methods Reference](advanced/methods-reference.md) — All tools and parameters
