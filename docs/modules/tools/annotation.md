# Annotation Tool Documentation

## Overview

The `annotation.py` module provides comprehensive cell type annotation functionality for spatial transcriptomics data. It supports multiple annotation methods ranging from simple marker-based approaches to sophisticated machine learning techniques, including integration with large language models.

## Purpose

Cell type annotation is a critical step in spatial transcriptomics analysis that assigns cell type labels to each spot or cell based on gene expression profiles. This module provides a unified interface to multiple annotation methods, allowing users to choose the most appropriate approach for their data.

## Available Methods

### 1. Marker Genes (`marker_genes`)
- **Description**: Simple and fast annotation using known marker genes
- **Best for**: Quick annotation when marker genes are well-established
- **Requirements**: List of marker genes per cell type

### 2. Correlation (`correlation`)
- **Description**: Assigns cell types based on expression correlation with reference data
- **Best for**: When you have well-annotated reference data
- **Requirements**: Reference dataset with cell type labels

### 3. Tangram (`tangram`)
- **Description**: Maps single-cell reference data to spatial data using optimal transport
- **Best for**: High-resolution spatial mapping with single-cell reference
- **Requirements**: Single-cell reference dataset
- **Key features**: Preserves spatial organization, handles batch effects

### 4. scANVI (`scanvi`)
- **Description**: Semi-supervised annotation using deep learning (from scvi-tools)
- **Best for**: Large datasets with partial annotations
- **Requirements**: Reference data or partial annotations
- **Key features**: Handles batch effects, provides uncertainty estimates

### 5. CellAssign (`cellassign`)
- **Description**: Probabilistic cell type assignment using marker genes
- **Best for**: When marker genes are known but expression is variable
- **Requirements**: Marker gene matrix
- **Key features**: Accounts for technical variability

### 6. mLLMCellType (`mllmcelltype`)
- **Description**: LLM-based cell type annotation using GPT models
- **Best for**: Exploratory analysis or when traditional methods fail
- **Requirements**: OpenAI API key
- **Key features**: Natural language understanding, handles novel cell types

## Input Parameters

### AnnotationParameters
```python
class AnnotationParameters:
    method: str = "marker_genes"  # Annotation method
    marker_genes: Dict[str, List[str]] = None  # Cell type markers
    reference_data_id: str = None  # Reference dataset ID
    batch_key: str = None  # Batch effect key
    n_top_genes: int = 2000  # Number of variable genes
    layer: str = None  # Expression layer to use
    confidence_threshold: float = 0.7  # Min confidence score
    
    # Method-specific parameters
    scanvi_params: Dict = {
        "n_latent": 30,
        "n_layers": 2,
        "dropout_rate": 0.1,
        "gene_likelihood": "zinb"
    }
    
    tangram_params: Dict = {
        "mode": "cells",
        "density_prior": "uniform",
        "num_epochs": 1000,
        "learning_rate": 0.1
    }
    
    cellassign_params: Dict = {
        "shrinkage": 3,
        "n_batches": 32
    }
    
    mllm_params: Dict = {
        "model": "gpt-4o-mini",
        "temperature": 0.1,
        "batch_size": 500,
        "use_cell_ontology": True
    }
```

## Output Format

### AnnotationResult
```python
class AnnotationResult:
    cell_types: List[str]  # Assigned cell types
    cell_type_counts: Dict[str, int]  # Count per type
    confidence_scores: Optional[List[float]]  # Confidence per cell
    method_used: str  # Annotation method
    parameters_used: Dict[str, Any]  # Parameters
    annotation_key: str  # Key in adata.obs
    
    # Method-specific outputs
    tangram_mapping_score: Optional[float]  # Tangram quality
    scanvi_latent_key: Optional[str]  # scANVI embeddings
    cellassign_probability_key: Optional[str]  # Probabilities
```

## Implementation Details

### Marker Gene Annotation
1. Calculates mean expression of marker genes per cell
2. Assigns cell type with highest mean expression
3. Applies confidence threshold
4. Fast but requires good markers

### Tangram Mapping
1. Aligns single-cell reference to spatial data
2. Uses optimal transport algorithm
3. Preserves spatial constraints
4. Returns mapping quality score

### scANVI Integration
1. Creates semi-supervised VAE model
2. Leverages labeled reference data
3. Provides uncertainty quantification
4. Handles batch effects naturally

### mLLMCellType Workflow
1. Selects top variable genes
2. Batches cells for API efficiency
3. Constructs prompts with expression data
4. Parses LLM responses
5. Handles retries and errors

## Default Marker Genes

The module includes predefined markers for common cell types:

```python
DEFAULT_MARKERS = {
    "T_cells": ["CD3D", "CD3E", "CD8A", "CD4"],
    "B_cells": ["CD19", "MS4A1", "CD79A"],
    "Myeloid": ["CD14", "FCGR3A", "CD68", "ITGAM"],
    "Fibroblasts": ["COL1A1", "COL3A1", "PDGFRA"],
    "Endothelial": ["PECAM1", "VWF", "CDH5"],
    "Epithelial": ["EPCAM", "KRT19", "KRT18"],
    "NK_cells": ["NCAM1", "NKG7", "GNLY"],
    "Plasma_cells": ["SDC1", "IGKC", "JCHAIN"],
    "Neutrophils": ["S100A8", "S100A9", "FCGR3B"]
}
```

## Usage Examples

### Example 1: Simple Marker-Based Annotation
```python
# Using default markers
result = await annotate_cells(
    data_id="data_1",
    params=AnnotationParameters(method="marker_genes")
)

# Using custom markers
custom_markers = {
    "Tumor": ["MKI67", "TOP2A", "PCNA"],
    "Stroma": ["COL1A1", "VIM", "ACTA2"]
}
result = await annotate_cells(
    data_id="data_1",
    params=AnnotationParameters(
        method="marker_genes",
        marker_genes=custom_markers
    )
)
```

### Example 2: Reference-Based Annotation with Tangram
```python
# Load reference single-cell data
ref_data = await load_data("reference.h5ad", name="scRNA_ref")

# Map to spatial data
result = await annotate_cells(
    data_id="spatial_data",
    params=AnnotationParameters(
        method="tangram",
        reference_data_id=ref_data.id,
        tangram_params={
            "mode": "cells",
            "num_epochs": 2000
        }
    )
)
```

### Example 3: Semi-Supervised with scANVI
```python
result = await annotate_cells(
    data_id="data_1",
    params=AnnotationParameters(
        method="scanvi",
        reference_data_id="annotated_ref",
        batch_key="batch",
        scanvi_params={
            "n_latent": 50,
            "gene_likelihood": "nb"
        }
    )
)
```

### Example 4: LLM-Based Annotation
```python
import os
os.environ["OPENAI_API_KEY"] = "your-key"

result = await annotate_cells(
    data_id="data_1",
    params=AnnotationParameters(
        method="mllmcelltype",
        mllm_params={
            "model": "gpt-4o",
            "temperature": 0.1,
            "use_cell_ontology": True
        }
    )
)
```

### Example 5: Ensemble Approach
```python
# First pass with markers
marker_result = await annotate_cells(
    data_id="data_1",
    params=AnnotationParameters(method="marker_genes")
)

# Refine uncertain cells with LLM
uncertain_mask = marker_result.confidence_scores < 0.5
if uncertain_mask.any():
    llm_result = await annotate_cells(
        data_id="data_1",
        params=AnnotationParameters(
            method="mllmcelltype",
            cell_subset=uncertain_mask
        )
    )
```

## Error Handling

The module implements comprehensive error handling:

1. **Method Validation**: Checks if method is supported
2. **Parameter Validation**: Validates method-specific requirements
3. **Reference Data Validation**: Ensures reference data exists and is compatible
4. **API Error Handling**: Retries and graceful degradation for LLM methods
5. **Memory Management**: Handles large datasets efficiently

## Best Practices

1. **Method Selection**:
   - Use marker genes for quick initial annotation
   - Use Tangram/scANVI for high-quality results with reference
   - Use mLLMCellType for exploration or difficult cases

2. **Parameter Tuning**:
   - Adjust confidence_threshold based on data quality
   - Increase n_top_genes for complex tissues
   - Tune method-specific parameters iteratively

3. **Validation**:
   - Always visualize results spatially
   - Check cell type proportions
   - Validate against known spatial patterns

4. **Performance**:
   - Marker genes: Fastest, scales linearly
   - Tangram/scANVI: Slower, GPU acceleration recommended
   - mLLMCellType: Rate-limited by API, use batching

## Integration Requirements

### Python Packages
- **Core**: scanpy, anndata, pandas, numpy
- **Tangram**: tangram-sc (requires separate installation)
- **scvi-tools**: scvi-tools>=0.20.0
- **mLLMCellType**: openai, scikit-learn

### External Services
- **mLLMCellType**: Requires OpenAI API key
- **GPU**: Recommended for deep learning methods

## Limitations and Considerations

1. **Marker Genes**: Quality depends on marker specificity
2. **Reference-Based**: Results limited by reference quality
3. **Tangram**: Computationally intensive for large datasets
4. **mLLMCellType**: API costs and rate limits apply
5. **Batch Effects**: Only some methods handle batch effects

## Future Enhancements

1. Support for additional methods (SingleR, Seurat)
2. Multi-modal integration (protein + RNA)
3. Spatial-aware annotation methods
4. Active learning for iterative refinement
5. Confidence calibration improvements