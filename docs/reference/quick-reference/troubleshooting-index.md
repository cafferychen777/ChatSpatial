# 🚨 Troubleshooting Quick Index

> **Fast Problem Resolution**: Use this index to jump to specific solutions in the [comprehensive troubleshooting guide](../troubleshooting/common_issues.md).

## 🔥 Emergency Quick Fixes

| Problem | One-Line Solution | Full Guide |
|---------|-------------------|------------|
| **MCP Server won't start** | Check Python path in MCP config | [MCP Connection Issues](../troubleshooting/common_issues.md#mcp-connection-issues) |
| **Out of memory** | `"Subsample data to 5000 cells"` | [Memory Problems](../troubleshooting/common_issues.md#large-dataset-handling) |
| **Analysis too slow** | Switch to Cherry Studio (3600s timeout) | [Performance Issues](../troubleshooting/common_issues.md#slow-performance) |
| **Data won't load** | Use absolute path: `"/full/path/to/file.h5ad"` | [Data Loading](../troubleshooting/common_issues.md#file-format-issues) |
| **No spatial coordinates** | Check `adata.obsm` keys: `"spatial"`, `"X_spatial"` | [Spatial Info Missing](../troubleshooting/common_issues.md#missing-spatial-information) |
| **Cell types not found** | Try `"Use scType for automatic annotation"` | [Annotation Fails](../troubleshooting/common_issues.md#cell-type-annotation-fails) |
| **Visualization blank** | `"Create simple scatter plot first"` | [Plot Problems](../troubleshooting/common_issues.md#plots-not-showing-or-corrupted) |
| **CUDA error** | `"Use CPU instead of GPU"` | [GPU Issues](../troubleshooting/common_issues.md#cuda-gpu-not-available) |

---

## 📋 Problem Categories

### 🔧 Installation & Setup
- **Import Errors**: Module not found, dependency conflicts
- **Environment Issues**: Conda vs pip, version mismatches  
- **MCP Configuration**: Server startup, connection problems
- **R Dependencies**: rpy2 installation for RCTD

**→ [Full Installation Guide](../troubleshooting/common_issues.md#installation-and-dependency-issues)**

### 💾 Performance & Memory
- **Memory Errors**: Out of RAM, CUDA memory issues
- **Speed Problems**: Analysis timeouts, slow performance
- **GPU Issues**: CUDA not available, GPU memory full
- **Large Datasets**: Handling 50K+ spots efficiently

**→ [Full Performance Guide](../troubleshooting/common_issues.md#memory-and-performance-problems)**

### 🔬 Analysis Issues  
- **Clustering Problems**: Too many/few clusters, poor domains
- **Annotation Failures**: No cell types found, marker issues
- **Spatial Analysis**: Missing coordinates, invalid patterns
- **Method Errors**: Algorithm-specific problems

**→ [Full Analysis Guide](../troubleshooting/common_issues.md#analysis-specific-issues)**

### 📊 Data & Visualization
- **Loading Errors**: File format issues, corrupted data
- **10X Visium**: Spatial folder structure problems
- **Plot Issues**: Blank/corrupted visualizations
- **Format Support**: Unsupported file types

**→ [Full Data Guide](../troubleshooting/common_issues.md#data-format-and-loading-problems)**

### 💬 Agent Communication
- **Unclear Responses**: Generic answers, wrong suggestions
- **Parameter Issues**: Inappropriate method choices
- **Progress Monitoring**: Long analyses without updates
- **Method Compatibility**: Incompatible analysis combinations

**→ [Full Communication Guide](../troubleshooting/common_issues.md#conversation-troubleshooting)**

---

## 🎯 Quick Diagnostic Questions

### "Is it working at all?"
```bash
# Test basic functionality
python -c "import chatspatial; print('✅ ChatSpatial installed')"
python -m chatspatial --help
npx @modelcontextprotocol/inspector python -m chatspatial
```

### "What's actually wrong?"
```text
"Show me system information and current dataset status"
"List what analyses have been completed"
"Display error details and suggest solutions"
```

### "How do I work around this?"
```text
"Suggest alternative methods for my analysis"
"What's the simplest way to get basic results?"
"Use CPU-only methods for this analysis"
```

---

## ⚡ Common Error Messages & Quick Fixes

| Error Message | Likely Cause | Quick Fix |
|---------------|--------------|-----------|
| `Dataset not found in data store` | Data not loaded | `"Load my data from /path/to/file"` |
| `Missing required dependencies` | Optional package missing | Check with dependency report |
| `Memory allocation failed` | Insufficient RAM | `"Process in smaller chunks"` |
| `CUDA out of memory` | GPU memory full | `"Use CPU mode"` |
| `No significant results found` | Parameters too strict | `"Use less stringent thresholds"` |
| `Visualization failed to render` | Data/parameter issue | `"Create basic plot first"` |
| `File not found` | Path error | Use absolute paths |
| `Unsupported file format` | Wrong file type | Check supported formats |
| `No spatial coordinates found` | Missing spatial data | Verify spatial info in data |
| `Annotation method failed` | Incompatible data/method | Try alternative method |

---

## 🏥 Emergency Protocols

### 🚑 **Critical Failure Recovery**
1. **Complete System Reset**:
   ```bash
   # Nuclear option: fresh environment
   conda create -n chatspatial-emergency python=3.10
   conda activate chatspatial-emergency
   pip install -e "."  # Minimal installation
   ```

2. **Minimal Working Analysis**:
   ```text
   "Load data with basic parameters only"
   "Use simple Leiden clustering"
   "Create basic spatial visualization" 
   "Export results immediately"
   ```

### 🔧 **Partial Recovery Strategies**
1. **Data Recovery**: Use subsets, simpler formats
2. **Method Fallbacks**: Basic algorithms, CPU-only
3. **Progressive Analysis**: Step-by-step validation
4. **Alternative Tools**: Python scripts, manual analysis

---

## 📞 When to Get Help

### ✅ **Try These First**
- [ ] Checked this troubleshooting index
- [ ] Read the [full troubleshooting guide](../troubleshooting/common_issues.md)
- [ ] Tested with minimal example
- [ ] Tried alternative methods
- [ ] Checked system resources

### 🆘 **Escalate If**
- Multiple error types occurring
- System crashes or freezes
- Data corruption suspected  
- Critical deadline pressure
- Novel error not in documentation

### 📋 **Information to Gather**
```bash
# System info
python --version
pip list | grep chatspatial
uname -a

# Error details  
- Full error message
- Steps to reproduce
- Data characteristics
- System specifications
```

---

## 💡 Prevention Tips

### 🎯 **Best Practices**
- **Test with small data first** before full analysis
- **Use absolute file paths** to avoid path errors
- **Monitor system resources** during analysis
- **Save intermediate results** to avoid recomputation
- **Keep backup environments** for critical work

### ⚠️ **Common Pitfalls**
- Using `git add -A` (never do this!)
- Mixing conda and pip installations  
- Running large analyses without testing
- Ignoring memory/timeout warnings
- Not checking data quality first

---

**🚨 Remember**: Most ChatSpatial issues have quick solutions. Start with the one-line fixes above, then escalate to the comprehensive troubleshooting guide if needed. When in doubt, ask the agent: *"Diagnose my current issue and suggest the fastest solution"*.

---

## 📚 Related Resources

- **[Complete Troubleshooting Guide](../troubleshooting/common_issues.md)** - Comprehensive solutions
- **[Installation Guide](../../getting-started/installation.md)** - Setup instructions
- **[API Reference](../api/)** - Technical documentation
- **[Common Workflows](./common-workflows.md)** - Standard analysis patterns
- **[All Tools Reference](./all-tools.md)** - Complete tool documentation