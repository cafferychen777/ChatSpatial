---
layout: default
title: Advanced Analysis
parent: Tutorials
nav_order: 6
description: "Advanced multi-modal spatial transcriptomics workflows"
---

# Advanced Analysis

This tutorial covers advanced multi-modal spatial transcriptomics analysis workflows using ChatSpatial. Learn how to combine multiple analysis methods and handle complex datasets.

## Overview

Advanced analysis in ChatSpatial involves:
- Multi-modal data integration
- Complex analytical workflows
- Cross-technology comparisons
- Advanced visualization techniques
- Performance optimization strategies

## Prerequisites

- Complete the [Basic Workflow](core/basic_spatial_analysis.md) tutorial
- Understand spatial transcriptomics technologies
- Familiarity with Python and data analysis concepts

## Multi-Modal Analysis Workflow

### 1. Loading Multiple Datasets

```python
# Load different spatial datasets
visium_result = load_data(
    data_path="data/mouse_brain_visium.h5ad", 
    name="visium_brain",
    data_type="visium"
)

merfish_result = load_data(
    data_path="data/mouse_brain_merfish.h5ad",
    name="merfish_brain", 
    data_type="merfish"
)
```

### 2. Cross-Platform Data Integration

```python
# Preprocess both datasets with consistent parameters
for data_id in [visium_result.id, merfish_result.id]:
    preprocess_data(
        data_id=data_id,
        normalization_method="sctransform",
        filter_genes=True,
        min_cells=5
    )
```

### 3. Comparative Cell Type Annotation

Use multiple annotation methods to cross-validate results:

```python
# Method 1: Reference-based annotation
annotate_cells(
    data_id=visium_result.id,
    method="tangram",
    reference_path="reference/mouse_brain_sc.h5ad"
)

# Method 2: Marker-based annotation  
annotate_cells(
    data_id=visium_result.id,
    method="sctype",
)

# Method 3: Probabilistic deconvolution
annotate_cells(
    data_id=visium_result.id,
    method="cell2location",
    reference_path="reference/mouse_brain_sc.h5ad"
)
```

### 4. Advanced Spatial Domain Analysis

Compare spatial organization across technologies:

```python
# High-resolution method for MERFISH
identify_spatial_domains(
    data_id=merfish_result.id,
    method="leiden", 
    n_domains=15,
    resolution=0.5
)

# Graph-based method for Visium
identify_spatial_domains(
    data_id=visium_result.id,
    method="spagcn",
    n_domains=10,
    resolution=1.0
)
```

## Advanced Visualization Strategies

### Multi-Panel Comparative Plots

```python
# Create comprehensive visualization
visualize_data(
    data_id=visium_result.id,
    plot_type="spatial_domains",
    save_path="results/visium_domains.pdf"
)

visualize_data(
    data_id=merfish_result.id, 
    plot_type="spatial_domains",
    save_path="results/merfish_domains.pdf"
)
```

### Custom Gene Expression Analysis

```python
# Identify technology-specific spatial genes
spatial_genes_visium = identify_spatial_genes(
    data_id=visium_result.id,
    method="spatialde",
    n_genes=100
)

spatial_genes_merfish = identify_spatial_genes(
    data_id=merfish_result.id,
    method="spatialde", 
    n_genes=100
)

# Visualize shared spatial genes
shared_genes = list(set(spatial_genes_visium.genes) & set(spatial_genes_merfish.genes))

visualize_data(
    data_id=visium_result.id,
    plot_type="gene_expression",
    genes=shared_genes[:6],
    save_path="results/shared_spatial_genes.pdf"
)
```

## Advanced Cell Communication Analysis

### Multi-Method Communication Analysis

```python
# LIANA+ comprehensive analysis
comm_liana = analyze_cell_communication(
    data_id=visium_result.id,
    method="liana",
    cell_type_column="cell_type"
)

# CellChat pathway analysis
comm_cellchat = analyze_cell_communication(
    data_id=visium_result.id,
    method="cellchat",
    cell_type_column="cell_type"
)

# Cross-validate results
visualize_data(
    data_id=visium_result.id,
    plot_type="cell_communication",
    save_path="results/communication_networks.pdf"
)
```

### Spatial Communication Patterns

```python
# Analyze distance-dependent communication
analyze_spatial_data(
    data_id=visium_result.id,
    analysis_type="spatial_communication",
    max_distance=100.0
)
```

## Performance Optimization

### Memory-Efficient Processing

```python
# Process large datasets in chunks
def process_large_dataset(data_path, chunk_size=5000):
    # Load with subsampling
    result = load_data(
        data_path=data_path,
        name="large_dataset",
        subsample=chunk_size
    )
    
    # Efficient preprocessing
    preprocess_data(
        data_id=result.id,
        normalization_method="pearson_residuals",  # Memory efficient
        filter_genes=True
    )
    
    return result
```

### Batch Processing Workflows

```python
# Process multiple samples efficiently
sample_paths = [
    "data/sample_1.h5ad",
    "data/sample_2.h5ad", 
    "data/sample_3.h5ad"
]

results = []
for i, path in enumerate(sample_paths):
    result = load_data(path, f"sample_{i+1}")
    preprocess_data(result.id)
    annotate_cells(result.id, method="sctype")
    results.append(result)
```

## Quality Control and Validation

### Cross-Method Validation

```python
# Validate spatial domains with multiple methods
methods = ["spagcn", "stagate", "leiden"]
domain_results = {}

for method in methods:
    domain_results[method] = identify_spatial_domains(
        data_id=visium_result.id,
        method=method,
        n_domains=8
    )

# Compare consistency across methods
analyze_spatial_data(
    data_id=visium_result.id,
    analysis_type="domain_consistency",
    methods=methods
)
```

### Statistical Validation

```python
# Marker gene validation
markers = find_markers(
    data_id=visium_result.id,
    group_by="spatial_domains",
    method="wilcoxon"
)

# Spatial autocorrelation analysis
analyze_spatial_data(
    data_id=visium_result.id,
    analysis_type="spatial_autocorrelation",
    genes=markers.genes[:20]
)
```

## Advanced Use Cases

### 1. Developmental Trajectory Analysis

```python
# Identify developmental trajectories in spatial context
analyze_spatial_data(
    data_id=visium_result.id,
    analysis_type="trajectory_inference",
    method="slingshot",
    root_cell_type="stem_cells"
)

# Visualize trajectories
visualize_data(
    data_id=visium_result.id,
    plot_type="trajectory",
    save_path="results/developmental_trajectories.pdf"
)
```

### 2. Disease vs Healthy Comparison

```python
# Load disease and control samples
disease_data = load_data("data/disease_sample.h5ad", "disease")
control_data = load_data("data/control_sample.h5ad", "control")

# Comparative analysis
for data_id in [disease_data.id, control_data.id]:
    preprocess_data(data_id)
    identify_spatial_domains(data_id, method="spagcn")
    
# Find disease-specific patterns
analyze_spatial_data(
    data_id=disease_data.id,
    analysis_type="differential_spatial",
    reference_id=control_data.id
)
```

### 3. Time-Series Spatial Analysis

```python
# Load time course data
timepoints = ["t0", "t6", "t12", "t24"]
time_data = {}

for tp in timepoints:
    time_data[tp] = load_data(f"data/timecourse_{tp}.h5ad", tp)
    preprocess_data(time_data[tp].id)

# Analyze temporal changes
analyze_spatial_data(
    data_ids=[data.id for data in time_data.values()],
    analysis_type="temporal_dynamics",
    timepoints=timepoints
)
```

## Best Practices

### 1. Parameter Selection

- **Resolution tuning**: Start with default parameters, then fine-tune
- **Method selection**: Consider data characteristics and biological questions
- **Validation**: Always cross-validate with multiple approaches

### 2. Computational Efficiency

- **Memory management**: Monitor memory usage for large datasets
- **Parallel processing**: Utilize multiple cores when available
- **Caching**: Save intermediate results to avoid recomputation

### 3. Result Interpretation

- **Biological validation**: Confirm results with known biology
- **Statistical testing**: Apply appropriate statistical tests
- **Visualization**: Use multiple visualization approaches

## Troubleshooting

### Common Issues

1. **Memory errors with large datasets**
   - Use subsampling or chunked processing
   - Increase available RAM
   - Use more memory-efficient algorithms

2. **Method convergence failures**
   - Try different initialization parameters
   - Use alternative methods
   - Check data quality

3. **Inconsistent results across methods**
   - Validate with known markers
   - Check preprocessing consistency
   - Consider method-specific biases

### Performance Tips

- Use SSD storage for large datasets
- Optimize Python environment (conda/mamba)
- Consider GPU acceleration where available
- Monitor system resources during analysis

## Next Steps

- Explore [Visualization Gallery](visualization_gallery.md) for advanced plotting
- Learn about [Cell Communication Analysis](analysis/cell_communication_analysis.md)
- Check out [Case Studies](../case_studies/mouse_brain_visium.md) for real-world examples

This advanced tutorial provides the foundation for sophisticated spatial transcriptomics analysis workflows. Combine these techniques based on your specific research questions and data characteristics.