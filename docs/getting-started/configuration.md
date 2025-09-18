---
layout: default
title: Configuration
parent: Getting Started
nav_order: 3
---

# Configuration Guide
{: .no_toc }

Configure ChatSpatial for your specific needs and environment.
{: .fs-6 .fw-300 }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## MCP Client Configuration

### Claude Desktop

Configure ChatSpatial in Claude Desktop by editing the MCP configuration file:

**Location:**
- macOS: `~/Library/Application Support/Claude/claude_desktop_config.json`
- Windows: `%APPDATA%/Claude/claude_desktop_config.json`
- Linux: `~/.config/Claude/claude_desktop_config.json`

**Basic Configuration:**

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "python",
      "args": ["-m", "chatspatial"],
      "cwd": "/path/to/your/chatspatial"
    }
  }
}
```

**Advanced Configuration:**

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "python",
      "args": [
        "-m", "chatspatial",
        "--log-level", "INFO",
        "--max-memory", "8GB"
      ],
      "cwd": "/path/to/your/chatspatial",
      "env": {
        "CHATSPATIAL_DATA_DIR": "/path/to/your/data",
        "CHATSPATIAL_CACHE_DIR": "/path/to/cache"
      }
    }
  }
}
```

### Other MCP Clients

For other MCP-compatible clients, refer to their specific configuration documentation. ChatSpatial follows the standard MCP protocol.

## Environment Variables

Configure ChatSpatial behavior using environment variables:

### Data and Cache Directories

```bash
# Set custom data directory
export CHATSPATIAL_DATA_DIR="/path/to/your/spatial/data"

# Set custom cache directory
export CHATSPATIAL_CACHE_DIR="/path/to/cache"

# Set temporary directory
export CHATSPATIAL_TEMP_DIR="/path/to/temp"
```

### Performance Settings

```bash
# Set maximum memory usage
export CHATSPATIAL_MAX_MEMORY="8GB"

# Set number of CPU cores to use
export CHATSPATIAL_N_JOBS="4"

# Enable GPU acceleration (if available)
export CHATSPATIAL_USE_GPU="true"
```

### Logging Configuration

```bash
# Set log level (DEBUG, INFO, WARNING, ERROR)
export CHATSPATIAL_LOG_LEVEL="INFO"

# Set log file location
export CHATSPATIAL_LOG_FILE="/path/to/chatspatial.log"

# Enable verbose logging
export CHATSPATIAL_VERBOSE="true"
```

## Configuration File

Create a configuration file for persistent settings:

**File:** `~/.chatspatial/config.yml`

```yaml
# ChatSpatial Configuration File

# Data settings
data:
  default_data_dir: "/path/to/your/spatial/data"
  cache_dir: "/path/to/cache"
  temp_dir: "/tmp/chatspatial"
  
# Performance settings
performance:
  max_memory: "8GB"
  n_jobs: 4
  use_gpu: false
  chunk_size: 1000

# Logging settings
logging:
  level: "INFO"
  file: "/var/log/chatspatial.log"
  verbose: false

# Analysis defaults
analysis:
  default_clustering_resolution: 0.5
  default_n_neighbors: 15
  default_min_cells: 3
  
# Visualization settings
visualization:
  default_figsize: [8, 6]
  default_dpi: 300
  color_palette: "viridis"
  
# Method preferences
methods:
  preferred_deconvolution: "cell2location"
  preferred_spatial_domains: "spagcn"
  preferred_trajectory: "cellrank"
```

## Method-Specific Configuration

### Deep Learning Methods

Configure PyTorch and GPU settings:

```yaml
# Advanced configuration for deep learning
deep_learning:
  device: "auto"  # auto, cpu, cuda, mps
  precision: "float32"  # float32, float16
  batch_size: 128
  max_epochs: 1000
  early_stopping: true
  
# Cell2location specific settings
cell2location:
  max_epochs: 30000
  batch_size: null  # auto-determined
  use_gpu: true
  learning_rate: 0.002
```

### R Integration

Configure R-based methods:

```yaml
# R integration settings
r_integration:
  r_home: "/usr/local/lib/R"  # Path to R installation
  r_libs: "/usr/local/lib/R/library"  # R library path
  rpy2_timeout: 300  # Timeout for R operations
  
# RCTD specific settings
rctd:
  doublet_mode: "doublet"
  confidence_threshold: 25
  max_cores: 4
```

## Security Settings

Configure security and privacy settings:

```yaml
# Security configuration
security:
  # Data privacy
  allow_external_data: false
  mask_sensitive_info: true
  
  # Network settings
  allow_internet_access: false
  proxy_url: null
  
  # File system access
  restrict_file_access: true
  allowed_directories:
    - "/path/to/safe/data"
    - "/path/to/safe/output"
```

## Development Configuration

Settings for development and debugging:

```yaml
# Development settings
development:
  debug_mode: false
  profile_performance: false
  save_intermediate_results: false
  
# Testing configuration
testing:
  use_small_datasets: true
  skip_slow_tests: false
  test_data_dir: "/path/to/test/data"
```

## Validation and Debugging

### Validate Configuration

```bash
# Check configuration syntax
python -m chatspatial validate-config

# Show current configuration
python -m chatspatial show-config

# Test configuration with demo data
python -m chatspatial test-config
```

### Debug Configuration Issues

```bash
# Verbose configuration loading
python -m chatspatial --debug show-config

# Check environment variables
python -c "
import os
from chatspatial.config import get_config
config = get_config()
print('Configuration loaded successfully')
print(f'Data dir: {config.data.default_data_dir}')
"
```

## Advanced Configuration

### Custom Plugin Configuration

```yaml
# Plugin configuration
plugins:
  enabled_plugins:
    - "custom_deconvolution"
    - "experimental_spatial"
  
  plugin_directories:
    - "/path/to/custom/plugins"
    
  plugin_settings:
    custom_deconvolution:
      algorithm: "advanced_nnls"
      regularization: 0.01
```

### Resource Limits

```yaml
# Resource management
resources:
  # Memory limits
  max_memory_per_dataset: "4GB"
  max_total_memory: "16GB"
  
  # CPU limits
  max_cpu_percent: 80
  max_concurrent_analyses: 2
  
  # Storage limits
  max_cache_size: "10GB"
  max_temp_size: "5GB"
```

### Integration Settings

```yaml
# External tool integration
integrations:
  # Scanpy settings
  scanpy:
    settings:
      verbosity: 1
      set_figure_params:
        dpi: 80
        facecolor: "white"
  
  # Squidpy settings
  squidpy:
    settings:
      verbosity: 1
```

## Configuration Examples

### High-Performance Setup

```yaml
# Configuration for high-performance computing
performance:
  max_memory: "64GB"
  n_jobs: 16
  use_gpu: true
  chunk_size: 5000

deep_learning:
  device: "cuda"
  precision: "float16"
  batch_size: 512
```

### Memory-Constrained Setup

```yaml
# Configuration for limited memory environments
performance:
  max_memory: "2GB"
  n_jobs: 2
  use_gpu: false
  chunk_size: 100

analysis:
  use_sparse_matrices: true
  enable_memory_optimization: true
```

### Production Setup

```yaml
# Production environment configuration
logging:
  level: "WARNING"
  file: "/var/log/chatspatial/production.log"
  
security:
  allow_external_data: false
  restrict_file_access: true
  
performance:
  max_memory: "16GB"
  n_jobs: 8
```

## Troubleshooting Configuration

### Common Issues

{: .note }
**Problem:** Configuration file not found

**Solution:** Create the configuration directory and file:
```bash
mkdir -p ~/.chatspatial
cp /path/to/chatspatial/config/default.yml ~/.chatspatial/config.yml
```

{: .note }
**Problem:** Invalid YAML syntax

**Solution:** Validate your configuration file:
```bash
python -c "import yaml; yaml.safe_load(open('~/.chatspatial/config.yml'))"
```

{: .note }
**Problem:** Environment variables not recognized

**Solution:** Check environment variable names and restart the MCP server.

### Reset Configuration

To reset to default configuration:

```bash
# Backup current configuration
cp ~/.chatspatial/config.yml ~/.chatspatial/config.yml.backup

# Reset to defaults
python -m chatspatial reset-config

# Or manually remove configuration
rm ~/.chatspatial/config.yml
```

---

**Next:** Explore the [Tutorials](../tutorials/) to start analyzing spatial data