# Configuration Guide

Complete guide for configuring ChatSpatial with various MCP clients and customizing settings.

## MCP Client Configuration

### Claude Desktop

#### macOS Configuration

1. **Locate config file**:
   ```bash
   ~/Library/Application Support/Claude/claude_desktop_config.json
   ```

2. **Basic configuration**:
   ```json
   {
     "mcpServers": {
       "chatspatial": {
         "command": "/Users/username/miniconda3/envs/chatspatial/bin/python",
         "args": ["-m", "chatspatial"],
         "env": {
           "PYTHONPATH": "/path/to/ChatSpatial"
         }
       }
     }
   }
   ```

3. **Find your Python path**:
   ```bash
   conda activate chatspatial
   which python
   ```

#### Windows Configuration

1. **Config file location**:
   ```
   %APPDATA%\Claude\claude_desktop_config.json
   ```

2. **Configuration example**:
   ```json
   {
     "mcpServers": {
       "chatspatial": {
         "command": "C:\\Users\\username\\miniconda3\\envs\\chatspatial\\python.exe",
         "args": ["-m", "chatspatial"],
         "env": {
           "PYTHONPATH": "C:\\path\\to\\ChatSpatial"
         }
       }
     }
   }
   ```

#### Linux Configuration

1. **Config file location**:
   ```bash
   ~/.config/Claude/claude_desktop_config.json
   ```

2. **Configuration example**:
   ```json
   {
     "mcpServers": {
       "chatspatial": {
         "command": "/home/username/miniconda3/envs/chatspatial/bin/python",
         "args": ["-m", "chatspatial"],
         "env": {
           "PYTHONPATH": "/home/username/ChatSpatial"
         }
       }
     }
   }
   ```

### Advanced MCP Configuration

#### Custom Server Settings

```json
{
  "mcpServers": {
    "chatspatial": {
      "command": "/path/to/python",
      "args": ["-m", "chatspatial"],
      "env": {
        "PYTHONPATH": "/path/to/ChatSpatial",
        "CHATSPATIAL_DATA_DIR": "/path/to/data",
        "CHATSPATIAL_CACHE_DIR": "/path/to/cache",
        "CHATSPATIAL_LOG_LEVEL": "INFO",
        "CHATSPATIAL_MAX_MEMORY": "16GB",
        "CHATSPATIAL_TIMEOUT": "300"
      }
    }
  }
}
```

#### Multiple Server Instances

```json
{
  "mcpServers": {
    "chatspatial-dev": {
      "command": "/path/to/dev/python",
      "args": ["-m", "chatspatial", "--dev"],
      "env": {
        "CHATSPATIAL_LOG_LEVEL": "DEBUG"
      }
    },
    "chatspatial-prod": {
      "command": "/path/to/prod/python",
      "args": ["-m", "chatspatial", "--prod"],
      "env": {
        "CHATSPATIAL_LOG_LEVEL": "WARNING"
      }
    }
  }
}
```

## ChatSpatial Configuration

### Configuration File

Create `~/.chatspatial/config.yaml`:

```yaml
# Data settings
data:
  default_data_dir: "~/chatspatial_data"
  cache_dir: "~/chatspatial_cache"
  max_file_size: "10GB"
  supported_formats: ["h5ad", "csv", "h5", "zarr"]

# Analysis settings
analysis:
  default_n_neighbors: 15
  default_resolution: 1.0
  max_genes: 50000
  max_spots: 100000
  
# Visualization settings
visualization:
  default_dpi: 300
  default_format: "png"
  max_image_size: "50MB"
  color_palettes:
    - "viridis"
    - "plasma"
    - "tab20"

# Performance settings
performance:
  n_jobs: -1
  memory_limit: "16GB"
  timeout: 300
  use_gpu: true
  
# Logging settings
logging:
  level: "INFO"
  file: "~/.chatspatial/logs/chatspatial.log"
  max_size: "100MB"
  backup_count: 5
```

### Environment Variables

```bash
# Data directories
export CHATSPATIAL_DATA_DIR="/path/to/data"
export CHATSPATIAL_CACHE_DIR="/path/to/cache"

# Performance settings
export CHATSPATIAL_N_JOBS="8"
export CHATSPATIAL_MEMORY_LIMIT="16GB"
export CHATSPATIAL_TIMEOUT="600"

# GPU settings
export CHATSPATIAL_USE_GPU="true"
export CUDA_VISIBLE_DEVICES="0"

# Logging
export CHATSPATIAL_LOG_LEVEL="INFO"
export CHATSPATIAL_LOG_FILE="/path/to/logs/chatspatial.log"
```

### Command Line Options

```bash
# Start with custom settings
python -m chatspatial \
  --data-dir /path/to/data \
  --cache-dir /path/to/cache \
  --log-level DEBUG \
  --timeout 600 \
  --max-memory 32GB

# Development mode
python -m chatspatial --dev

# Production mode
python -m chatspatial --prod --quiet

# Custom config file
python -m chatspatial --config /path/to/config.yaml
```

## Tool-Specific Configuration

### Spatial Analysis Tools

```yaml
spatial_analysis:
  spagcn:
    default_beta: 49
    default_alpha: 0.1
    max_iterations: 1000
  
  stagate:
    hidden_dims: [512, 30]
    learning_rate: 0.001
    max_epochs: 1000
    
  gaston:
    n_bins: 10
    regularization: 0.0
```

### Visualization Settings

```yaml
visualization:
  spatial_plots:
    default_spot_size: 1.0
    default_alpha: 0.8
    default_colormap: "viridis"
    
  umap_plots:
    default_n_neighbors: 15
    default_min_dist: 0.1
    
  heatmaps:
    default_cmap: "RdBu_r"
    cluster_method: "ward"
```

### Data Processing

```yaml
preprocessing:
  quality_control:
    min_genes_per_spot: 200
    min_spots_per_gene: 3
    max_genes_per_spot: 5000
    mt_threshold: 20.0
    
  normalization:
    target_sum: 10000
    highly_variable_genes: 2000
    
  dimensionality_reduction:
    n_pcs: 50
    n_neighbors: 15
```

## Security Configuration

### Access Control

```yaml
security:
  allowed_paths:
    - "/home/user/data"
    - "/shared/spatial_data"
  
  blocked_paths:
    - "/system"
    - "/etc"
    
  max_file_size: "10GB"
  allowed_extensions:
    - ".h5ad"
    - ".csv"
    - ".h5"
    - ".zarr"
```

### Resource Limits

```yaml
limits:
  max_memory_per_analysis: "8GB"
  max_execution_time: 600
  max_concurrent_analyses: 3
  max_file_uploads: 10
```

## Troubleshooting Configuration

### Verify MCP Connection

```bash
# Test MCP server
python -m chatspatial --test-mcp

# Check configuration
python -m chatspatial --check-config

# Validate environment
python -m chatspatial --validate-env
```

### Debug Configuration Issues

```bash
# Enable debug logging
export CHATSPATIAL_LOG_LEVEL=DEBUG

# Check Python path
python -c "import sys; print('\n'.join(sys.path))"

# Verify imports
python -c "import chatspatial; print('OK')"

# Test tools
python -c "
from chatspatial.tools.data_management import load_data
print('Tools import OK')
"
```

### Common Configuration Problems

#### 1. Python Path Issues
```bash
# Solution: Use absolute paths
which python  # Get full path
# Use this path in MCP configuration
```

#### 2. Permission Errors
```bash
# Solution: Check file permissions
chmod +x /path/to/python
chmod -R 755 /path/to/ChatSpatial
```

#### 3. Environment Variables
```bash
# Solution: Set in shell profile
echo 'export CHATSPATIAL_DATA_DIR="/path/to/data"' >> ~/.bashrc
source ~/.bashrc
```

#### 4. Memory Issues
```bash
# Solution: Adjust memory limits
export CHATSPATIAL_MEMORY_LIMIT="8GB"
# Or modify config.yaml
```

## Performance Optimization

### Memory Optimization

```yaml
performance:
  memory_optimization:
    use_sparse_matrices: true
    chunk_size: 1000
    lazy_loading: true
    garbage_collection: true
```

### CPU Optimization

```yaml
performance:
  cpu_optimization:
    n_jobs: -1  # Use all cores
    parallel_backend: "loky"
    thread_pool_size: 8
```

### GPU Configuration

```yaml
performance:
  gpu:
    enabled: true
    device: "cuda:0"
    memory_fraction: 0.8
    allow_growth: true
```

## Next Steps

After configuration:

1. **Test Connection**: Verify MCP client can connect
2. **Run Examples**: Test with sample data
3. **Monitor Performance**: Check logs and resource usage
4. **Customize Settings**: Adjust based on your needs

See [Getting Started](../getting-started/README.md) for usage examples!
