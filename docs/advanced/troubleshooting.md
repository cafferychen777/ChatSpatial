# Troubleshooting

This page is the canonical **symptom → fix** guide.

- For installation steps, see [Installation](../installation.md).
- For exact MCP client syntax, see [Configuration Guide](configuration.md).
- For a first-run workflow, see [Quick Start](../quickstart.md).

---

## MCP Connection Problems

### Tools not showing in the client

1. Confirm you used the correct config file for your client.
2. Confirm the Python path is an **absolute** path from `which python`.
3. Check the config file for JSON/TOML syntax errors.
4. Restart the client after configuration changes.
5. Test the server directly:

```bash
python -m chatspatial server --help
```

If you need the exact config file format, go back to the [Configuration Guide](configuration.md).

### "python not found" or "module not found"

- Make sure ChatSpatial is installed inside the environment you configured
- Re-run `which python` inside the activated environment
- Update the MCP config to use that exact path

---

## Docker / GHCR Problems

### `docker: command not found`

Install Docker Desktop or Docker Engine, confirm `docker --version` works, then restart your MCP client.

### Pull fails for the GHCR image

Check the image name and network access:

```bash
docker pull ghcr.io/cafferychen777/chatspatial:v1.2.6
```

### MCP tools do not appear when using Docker

- Use `--rm -i`, not `-it`, in MCP stdio configuration
- Use absolute host paths in `-v` mounts
- Restart the MCP client after changing configuration

### Dataset not found in Docker

Mount the host data directory and use the container path in prompts:

```bash
-v /Users/alice/spatial-data:/data:ro
```

```text
Load /data/sample.h5ad
```

Do not prompt with `/Users/alice/spatial-data/sample.h5ad`; that path exists on the host, not inside the container.

### Permission denied on mounted outputs

Confirm the host output directory exists and Docker has permission to write there. On Docker Desktop, also check file-sharing permissions for the mounted parent directory.

---

## Data Loading Problems

### "Dataset not found"

Use an **absolute** path:

```text
❌ ~/data/sample.h5ad
❌ ./data/sample.h5ad
✅ /Users/yourname/data/sample.h5ad
```

### File format not recognized

- **H5AD:** verify with `python -c "import scanpy as sc; sc.read_h5ad('file.h5ad')"`
- **Visium:** point to the directory containing the `spatial/` folder
- **HDF5 check:** `file yourdata.h5ad`

---

## Analysis Problems

### "Run preprocessing first"

Most analyses require preprocessing first.

```text
Preprocess the data
```

### "No significant results"

- check data quality (>500 spots, >1000 genes)
- lower significance thresholds
- try a different analysis method

### Cell communication fails

Use species/resource pairs that match the dataset:

```text
For mouse: species="mouse", liana_resource="mouseconsensus"
For human: species="human", liana_resource="consensus"
```

---

## Resource Problems

### System freezes / MemoryError

- subsample data for testing
- reduce batch sizes
- monitor memory with `top`
- use 32GB+ RAM or cloud resources for large datasets

### CUDA out of memory

- set `use_gpu=False`
- reduce batch size
- clear cached GPU memory if your workflow allows it

---

## Quick Fix Table

| Problem | First fix |
|---------|-----------|
| Import errors | Reinstall with `uv pip install chatspatial[full]` |
| `resolution-too-deep` | Use `uv` instead of `pip` |
| Client not connecting | Re-check config and restart the client |
| Docker pull fails | Run `docker pull ghcr.io/cafferychen777/chatspatial:v1.2.6` and check network access |
| Docker dataset not found | Mount the host data directory and prompt with `/data/...` |
| Path errors | Use absolute paths |
| Analysis fails immediately | Run preprocessing first |
| R methods fail | Install R and the required R packages |

---

## Still Stuck?

- [FAQ](faq.md) — short answers and pointers
- [Configuration Guide](configuration.md) — exact client syntax
- [Methods Reference](methods-reference.md) — tool parameters and defaults
- [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues) — report reproducible bugs
