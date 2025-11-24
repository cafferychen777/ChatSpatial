
title: Performance
description: Performance optimization tips and hardware recommendations
---

# Performance Tips

Optimize ChatSpatial performance for faster analysis and better resource utilization.

## Hardware Recommendations

### Minimum Requirements
- **CPU**: 4 cores, 2.5GHz
- **RAM**: 8GB
- **Storage**: 100GB SSD
- **Network**: Stable internet connection

### Recommended Setup
- **CPU**: 8+ cores, 3.0GHz+ (Intel i7/AMD Ryzen 7)
- **RAM**: 32GB+ DDR4
- **Storage**: 500GB+ NVMe SSD
- **GPU**: NVIDIA RTX 3060+ (for deep learning methods)
- **Network**: High-speed internet for large data downloads

### High-Performance Setup
- **CPU**: 16+ cores, 3.5GHz+ (Intel i9/AMD Ryzen 9)
- **RAM**: 64GB+ DDR4/DDR5
- **Storage**: 1TB+ NVMe SSD (PCIe 4.0)
- **GPU**: NVIDIA RTX 4080+ or A100
- **Network**: Gigabit ethernet

## Memory Optimization

### Data Loading Optimization

```python
# Use sparse matrices for count data
import scipy.sparse as sp
adata.X = sp.csr_matrix(adata.X)

# Load only necessary data
adata = sc.read_h5ad('data.h5ad', backed='r')  # Read-only mode

# Subsample large datasets
adata = adata[::10, :].copy()  # Every 10th spot
adata = adata[:, adata.var.highly_variable].copy()  # HVGs only
```

### Memory-Efficient Processing

```python
# Process in chunks
def process_in_chunks(adata, chunk_size=1000):
    n_chunks = (adata.n_obs + chunk_size - 1) // chunk_size
    results = []
    
    for i in range(n_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, adata.n_obs)
        chunk = adata[start:end, :].copy()
        
        # Process chunk
        result = analyze_chunk(chunk)
        results.append(result)
        
        # Free memory
        del chunk
        gc.collect()
    
    return combine_results(results)
```

### Memory Monitoring

```python
import psutil
import gc

def monitor_memory():
    """Monitor memory usage."""
    process = psutil.Process()
    memory_info = process.memory_info()
    
    print(f"RSS: {memory_info.rss / 1024**3:.2f} GB")
    print(f"VMS: {memory_info.vms / 1024**3:.2f} GB")
    print(f"Available: {psutil.virtual_memory().available / 1024**3:.2f} GB")

# Use throughout analysis
monitor_memory()
# ... run analysis ...
gc.collect()  # Force garbage collection
monitor_memory()
```

## CPU Optimization

### Parallel Processing

```python
# Set number of jobs
import os
os.environ['NUMBA_NUM_THREADS'] = '8'
os.environ['OMP_NUM_THREADS'] = '8'

# Preprocessing is automatically optimized
# Note: Parallel processing controlled at system level (NUMBA_NUM_THREADS, OMP_NUM_THREADS)
from chatspatial.tools import preprocess_data
result = preprocess_data(
    data_id="sample"
    # Parallelization handled internally via environment variables
)
```

### Efficient Algorithms

```python
# Use approximate algorithms for large datasets
from chatspatial.tools import identify_spatial_domains

# Fast approximate clustering
result = identify_spatial_domains(
    data_id="sample",
    method="leiden",  # Faster than spagcn method
    resolution=1.0
    # Leiden is fast by design, no iteration limiting needed
)

# For large datasets, subsample during preprocessing first
# Then run spatial domain identification
# preprocess_data(data_id="sample", subsample_spots=5000)
result = identify_spatial_domains(
    data_id="sample",
    method="spagcn"
    # Note: Subsample data during preprocessing, not here
)
```

### Threading Configuration

```python
# Optimize threading
import threadpoolctl

# Limit threads for specific operations
with threadpoolctl.threadpool_limits(limits=4, user_api='blas'):
    # CPU-intensive operation
    result = analyze_spatial_data(data_id="sample")
```

## GPU Acceleration

### CUDA Setup

```bash
# Install CUDA-enabled PyTorch
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# Verify CUDA
python -c "import torch; print(torch.cuda.is_available())"
```

### GPU-Accelerated Methods

```python
# Enable GPU for deep learning methods (deconvolution)
from chatspatial.tools import deconvolve_data

result = deconvolve_data(
    data_id="sample",
    method="cell2location",
    cell_type_key="cell_type",  # Required parameter
    reference_data_id="reference",  # Required parameter
    use_gpu=True
    # Note: GPU device selection handled automatically by scvi-tools
)

# GPU-accelerated spatial analysis (GraphST method)
result = identify_spatial_domains(
    data_id="sample",
    method="graphst",  # GraphST supports GPU
    graphst_use_gpu=True  # Enable GPU for GraphST
    # Note: STAGATE does not have GPU support, use GraphST instead
)
```

### GPU Memory Management

```python
import torch

# Clear GPU cache
torch.cuda.empty_cache()

# Monitor GPU memory
def gpu_memory_usage():
    if torch.cuda.is_available():
        allocated = torch.cuda.memory_allocated() / 1024**3
        cached = torch.cuda.memory_reserved() / 1024**3
        print(f"GPU Memory - Allocated: {allocated:.2f} GB, Cached: {cached:.2f} GB")

# Use mixed precision for memory efficiency
from torch.cuda.amp import autocast

with autocast():
    # GPU operations with reduced memory usage
    result = gpu_intensive_analysis()
```

## Storage Optimization

### SSD Configuration

```bash
# Check storage type
lsblk -d -o name,rota
# 0 = SSD, 1 = HDD

# Optimize SSD (Linux)
sudo fstrim -av  # TRIM unused blocks
```

### File Format Optimization

```python
# Use compressed HDF5
adata.write('data.h5ad', compression='gzip', compression_opts=9)

# Use Zarr for large datasets
adata.write_zarr('data.zarr', chunks=(1000, 1000))

# Optimize chunk sizes
import h5py
with h5py.File('data.h5ad', 'r') as f:
    print(f"Chunk size: {f['X'].chunks}")
```

### Data Caching

```python
# Enable caching
import os
os.environ['CHATSPATIAL_CACHE_DIR'] = '/fast/cache/directory'

# Use memory mapping
adata = sc.read_h5ad('data.h5ad', backed='r')
# Data stays on disk, loaded as needed
```

## Network Optimization

### Data Transfer

```python
# Compress data for transfer
import gzip
import pickle

# Compress results
with gzip.open('results.pkl.gz', 'wb') as f:
    pickle.dump(results, f)

# Use efficient serialization
import joblib
joblib.dump(results, 'results.joblib', compress=3)
```

### Batch Processing

```python
# Process multiple samples efficiently
samples = ['sample1', 'sample2', 'sample3']

# Batch load
for sample in samples:
    load_data(f"data/{sample}.h5ad", sample)

# Batch process
results = {}
for sample in samples:
    results[sample] = preprocess_data(data_id=sample)
```

## Analysis-Specific Optimizations

### Spatial Domain Identification

```python
# Fast methods for exploration
quick_result = identify_spatial_domains(
    data_id="sample",
    method="leiden",  # Fastest
    resolution=1.0
)

# High-quality methods for final analysis
final_result = identify_spatial_domains(
    data_id="sample",
    method="spagcn",  # More accurate
    n_domains=7,  # Number of spatial domains
    refine_domains=True  # Enable domain refinement
)
```

### Cell Communication Analysis

```python
# Optimize for large datasets
result = analyze_cell_communication(
    data_id="sample",
    method="liana",
    cell_type_column="cell_type",  # Required parameter
    perform_spatial_analysis=True,  # Spatial bivariate analysis
    liana_n_perms=100  # Reduce permutations for speed
    # Note: For speed, reduce n_perms or subsample data during preprocessing
)
```

### Visualization Optimization

```python
# Reduce image resolution for exploration
plot = visualize_data(
    data_id="sample",
    plot_type="spatial",
    dpi=150,  # Lower DPI
    size=(8, 6)  # Smaller size
)

# High-resolution for publication
final_plot = visualize_data(
    data_id="sample",
    plot_type="spatial",
    dpi=300,  # High DPI
    size=(12, 10),
    format="svg"  # Vector format
)
```

## Monitoring and Profiling

### Performance Monitoring

```python
import time
import psutil

class PerformanceMonitor:
    def __init__(self):
        self.start_time = time.time()
        self.start_memory = psutil.virtual_memory().used
    
    def checkpoint(self, name):
        current_time = time.time()
        current_memory = psutil.virtual_memory().used
        
        elapsed = current_time - self.start_time
        memory_diff = (current_memory - self.start_memory) / 1024**3
        
        print(f"{name}: {elapsed:.2f}s, Memory: {memory_diff:+.2f} GB")

# Usage
monitor = PerformanceMonitor()
preprocess_data(data_id="sample")
monitor.checkpoint("Preprocessing")
identify_spatial_domains(data_id="sample")
monitor.checkpoint("Spatial domains")
```

### Profiling Code

```python
# Profile with cProfile
import cProfile
import pstats

profiler = cProfile.Profile()
profiler.enable()

# Your analysis code here
result = analyze_spatial_data(data_id="sample")

profiler.disable()
stats = pstats.Stats(profiler)
stats.sort_stats('cumulative').print_stats(10)
```

### Memory Profiling

```python
# Install memory_profiler
# pip install memory-profiler

from memory_profiler import profile

@profile
def memory_intensive_function():
    # Your analysis code
    return result

# Run with: python -m memory_profiler script.py
```

## Troubleshooting Performance Issues

### Common Bottlenecks

1. **Memory Issues**
   - Reduce dataset size
   - Use sparse matrices
   - Process in chunks

2. **CPU Bottlenecks**
   - Increase parallelization
   - Use faster algorithms
   - Optimize threading

3. **I/O Bottlenecks**
   - Use SSD storage
   - Optimize file formats
   - Enable caching

4. **Network Issues**
   - Use local data
   - Compress transfers
   - Batch operations

### Performance Tuning Checklist

- [ ] Use appropriate hardware
- [ ] Optimize memory usage
- [ ] Enable parallel processing
- [ ] Use GPU acceleration when available
- [ ] Choose efficient algorithms
- [ ] Monitor resource usage
- [ ] Profile bottlenecks
- [ ] Cache intermediate results

## Benchmarking

### Performance Benchmarks

```python
def benchmark_methods():
    """Benchmark different spatial domain methods."""
    methods = ['leiden', 'spagcn', 'stagate']
    results = {}
    
    for method in methods:
        start_time = time.time()
        result = identify_spatial_domains(
            data_id="sample",
            method=method
        )
        elapsed = time.time() - start_time
        results[method] = elapsed
        print(f"{method}: {elapsed:.2f}s")
    
    return results
```

### Dataset Size Guidelines

| Dataset Size | Recommended Method | Expected Time | Memory Usage |
|--------------|-------------------|---------------|--------------|
| < 1K spots | Any method | < 1 min | < 2GB |
| 1K-10K spots | SpaGCN, STAGATE | 1-10 min | 2-8GB |
| 10K-50K spots | Leiden, subsampling | 10-30 min | 8-16GB |
| > 50K spots | Chunked processing | 30+ min | 16+ GB |

## Next Steps

After optimizing performance:

1. **Monitor Usage**: Track resource consumption
2. **Benchmark Methods**: Compare different approaches
3. **Scale Analysis**: Handle larger datasets
4. **Automate Workflows**: Create efficient pipelines

See [Configuration Guide](configuration.md) for system setup!