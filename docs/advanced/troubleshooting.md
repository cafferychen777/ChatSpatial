# Troubleshooting

Quick solutions to common problems.

---

## MCP Connection

### Tools not showing in Claude

1. **Check config path:**
   - macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
   - Windows: `%APPDATA%\Claude\claude_desktop_config.json`

2. **Verify JSON syntax** (missing commas, brackets)

3. **Restart Claude** after config changes

4. **Test server:**
   ```bash
   python -m chatspatial server --help
   ```

---

## Data Loading

### "Dataset not found"

**Use absolute paths:**

```text
❌ ~/data/sample.h5ad
❌ ./data/sample.h5ad
✅ /Users/yourname/data/sample.h5ad
```

### File format not recognized

- **H5AD:** Verify with `python -c "import scanpy as sc; sc.read_h5ad('file.h5ad')"`
- **Visium:** Point to directory containing `spatial/` folder
- **Check file:** `file yourdata.h5ad` should show "HDF5"

---

## Analysis Errors

### "Run preprocessing first"

Most analyses require preprocessing. Ask Claude:
```text
"Preprocess the data"
```

### "No significant results"

- Check data quality (>500 spots, >1000 genes)
- Lower significance thresholds
- Try different methods

### Cell communication fails

```text
For mouse: species="mouse", liana_resource="mouseconsensus"
For human: species="human", liana_resource="consensus"
```

---

## Memory Issues

### System freezes / MemoryError

- Subsample data for testing
- Use smaller batch sizes
- Monitor with `top` command
- For large datasets: use 32GB+ RAM or cloud

### CUDA out of memory

- Set `use_gpu=False`
- Reduce batch size
- Run `torch.cuda.empty_cache()`

---

## Quick Fixes

| Problem | Solution |
|---------|----------|
| Import errors | `pip install --upgrade chatspatial[full]` |
| Claude not connecting | Restart Claude, check JSON config |
| Path errors | Use absolute paths |
| Analysis fails | Run preprocessing first |
| R methods fail | Install R + packages |

---

## Get Help

- [GitHub Issues](https://github.com/cafferychen777/ChatSpatial/issues) — Report bugs
- [FAQ](faq.md) — Common questions
