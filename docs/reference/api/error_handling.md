---
layout: default
title: Error Handling
parent: API Reference
grand_parent: Reference
nav_order: 2
description: "Error handling system and troubleshooting guide"
---

# Error Handling

This document describes the error handling system in ChatSpatial, including error codes, troubleshooting guides, and best practices for robust error recovery.

## Overview

ChatSpatial implements comprehensive error handling to provide clear, actionable feedback when operations fail. All errors follow a standardized format and include specific error codes for programmatic handling.

## Error Response Format

All ChatSpatial tools return errors in a consistent format:

```json
{
  "success": false,
  "error": {
    "type": "ErrorType",
    "code": "ERROR_CODE",
    "message": "Human-readable error description",
    "details": {
      "field": "specific_parameter",
      "value": "problematic_value",
      "suggestion": "How to fix the issue"
    },
    "traceback": "Full Python traceback (debug mode only)"
  }
}
```

## Error Categories

### 1. Validation Errors (VALIDATION_*)

Errors related to input parameter validation.

#### VALIDATION_SCHEMA
Invalid parameter schema or missing required fields.

```json
{
  "type": "ValidationError",
  "code": "VALIDATION_SCHEMA",
  "message": "Required parameter 'data_id' is missing",
  "details": {
    "required_fields": ["data_id"],
    "provided_fields": ["method"]
  }
}
```

#### VALIDATION_VALUE
Parameter values outside allowed ranges or types.

```json
{
  "type": "ValidationError", 
  "code": "VALIDATION_VALUE",
  "message": "Parameter 'method' has invalid value 'invalid_method'",
  "details": {
    "field": "method",
    "value": "invalid_method",
    "allowed_values": ["spagcn", "stagate", "banksy"],
    "suggestion": "Use one of the supported methods"
  }
}
```

### 2. Data Errors (DATA_*)

Errors related to data loading, access, and format issues.

#### DATA_NOT_FOUND
Referenced dataset or file not found.

```json
{
  "type": "DataError",
  "code": "DATA_NOT_FOUND", 
  "message": "Dataset with ID 'nonexistent_data' not found",
  "details": {
    "requested_id": "nonexistent_data",
    "available_ids": ["mouse_brain", "human_heart"],
    "suggestion": "Check the dataset ID or load the data first"
  }
}
```

#### DATA_FORMAT_ERROR
Data file format issues or corruption.

```json
{
  "type": "DataError",
  "code": "DATA_FORMAT_ERROR",
  "message": "Cannot read h5ad file: corrupted or invalid format",
  "details": {
    "file_path": "/path/to/data.h5ad",
    "error_details": "HDF5 file structure is invalid",
    "suggestion": "Verify the file is a valid AnnData h5ad file"
  }
}
```

#### DATA_MISSING_OBSM
Required spatial coordinate data missing.

```json
{
  "type": "DataError",
  "code": "DATA_MISSING_OBSM",
  "message": "Spatial coordinates not found in dataset",
  "details": {
    "required_key": "spatial",
    "available_keys": ["X_pca", "X_umap"],
    "suggestion": "Ensure the dataset contains spatial coordinate information"
  }
}
```

### 3. Analysis Errors (ANALYSIS_*)

Errors during analysis execution.

#### ANALYSIS_INSUFFICIENT_DATA
Not enough data for the requested analysis.

```json
{
  "type": "AnalysisError",
  "code": "ANALYSIS_INSUFFICIENT_DATA", 
  "message": "Insufficient cells for cell communication analysis",
  "details": {
    "min_required": 50,
    "actual_count": 23,
    "suggestion": "Use a larger dataset or adjust filtering parameters"
  }
}
```

#### ANALYSIS_CONVERGENCE_FAILURE
Algorithm failed to converge.

```json
{
  "type": "AnalysisError",
  "code": "ANALYSIS_CONVERGENCE_FAILURE",
  "message": "SpaGCN clustering failed to converge",
  "details": {
    "max_iterations": 1000,
    "final_iteration": 1000,
    "suggestion": "Try different parameters or use an alternative method"
  }
}
```

### 4. Resource Errors (RESOURCE_*)

System resource and dependency issues.

#### RESOURCE_MEMORY_ERROR
Insufficient memory for operation.

```json
{
  "type": "ResourceError",
  "code": "RESOURCE_MEMORY_ERROR",
  "message": "Insufficient memory for large dataset processing",
  "details": {
    "required_memory_gb": 32,
    "available_memory_gb": 8,
    "dataset_size": "50000 cells x 30000 genes",
    "suggestion": "Use data subsampling or increase available memory"
  }
}
```

#### RESOURCE_DEPENDENCY_ERROR
Missing or incompatible dependencies.

```json
{
  "type": "ResourceError",
  "code": "RESOURCE_DEPENDENCY_ERROR", 
  "message": "Required package 'tangram-sc' not installed",
  "details": {
    "missing_package": "tangram-sc",
    "install_command": "pip install tangram-sc",
    "suggestion": "Install the required package for Tangram annotation"
  }
}
```

### 5. Visualization Errors (VIZ_*)

Plotting and visualization specific errors.

#### VIZ_INVALID_PLOT_TYPE
Unsupported plot type requested.

```json
{
  "type": "VisualizationError",
  "code": "VIZ_INVALID_PLOT_TYPE",
  "message": "Plot type 'invalid_plot' is not supported",
  "details": {
    "requested_type": "invalid_plot", 
    "available_types": ["spatial", "umap", "heatmap", "violin"],
    "suggestion": "Use one of the supported plot types"
  }
}
```

#### VIZ_MISSING_DATA
Required data for visualization not available.

```json
{
  "type": "VisualizationError",
  "code": "VIZ_MISSING_DATA",
  "message": "Cannot create spatial plot: missing spatial coordinates",
  "details": {
    "plot_type": "spatial",
    "missing_data": "spatial coordinates",
    "suggestion": "Ensure the dataset has spatial coordinate information"
  }
}
```

## Error Recovery Strategies

### Automatic Recovery

ChatSpatial implements automatic recovery for certain error conditions:

1. **Parameter Auto-correction**: Automatically fix common parameter issues
2. **Fallback Methods**: Try alternative methods when the primary fails
3. **Data Repair**: Fix minor data format issues automatically

### Manual Recovery Guidelines

#### For Validation Errors
1. Check the parameter documentation
2. Verify required fields are provided
3. Ensure parameter values are within allowed ranges

#### For Data Errors
1. Verify file paths and permissions
2. Check data format compatibility
3. Ensure required data components are present

#### For Analysis Errors
1. Try alternative methods or parameters
2. Check data quality and preprocessing
3. Consider data subsampling for large datasets

#### For Resource Errors
1. Install missing dependencies
2. Increase available memory/compute resources  
3. Use more efficient analysis parameters

## Debugging and Logging

### Debug Mode

Enable detailed error reporting:

```python
import os
os.environ['CHATSPATIAL_DEBUG'] = '1'
```

This provides:
- Full Python tracebacks
- Detailed parameter information
- Internal state information
- Performance metrics

### Log Levels

- `ERROR`: Critical failures requiring user action
- `WARNING`: Issues that might affect results
- `INFO`: General operation information
- `DEBUG`: Detailed debugging information

## Common Error Scenarios

### Scenario 1: Dataset Not Preprocessed

**Error**: Trying to run analysis on raw data
```bash
ANALYSIS_MISSING_PREPROCESSING: Dataset must be preprocessed before analysis
```

**Solution**: Run preprocessing first
```python
preprocess_data(data_id="my_data")
```

### Scenario 2: Missing Cell Type Annotations

**Error**: Cell communication analysis without annotations
```bash
DATA_MISSING_ANNOTATIONS: Cell type annotations required for communication analysis
```

**Solution**: Annotate cells first
```python
annotate_cells(data_id="my_data", method="sctype")
```

### Scenario 3: Incompatible Analysis Method

**Error**: Method not suitable for data type
```bash
ANALYSIS_METHOD_INCOMPATIBLE: Tangram requires reference dataset
```

**Solution**: Provide reference or use different method
```python
annotate_cells(data_id="my_data", method="marker_based")
```

## Best Practices

1. **Always Check Success**: Verify `success: true` before using results
2. **Handle Errors Gracefully**: Implement proper error handling in agent workflows
3. **Use Specific Error Codes**: Match against error codes for programmatic handling
4. **Provide User Feedback**: Surface error messages to users with suggestions
5. **Log for Debugging**: Keep error logs for troubleshooting complex workflows

## Error Code Reference

| Code | Category | Description |
|------|----------|-------------|
| VALIDATION_SCHEMA | Validation | Missing required parameters |
| VALIDATION_VALUE | Validation | Invalid parameter values |
| DATA_NOT_FOUND | Data | Dataset/file not found |
| DATA_FORMAT_ERROR | Data | Invalid file format |
| DATA_MISSING_OBSM | Data | Missing spatial coordinates |
| DATA_MISSING_ANNOTATIONS | Data | Missing cell annotations |
| ANALYSIS_INSUFFICIENT_DATA | Analysis | Not enough data |
| ANALYSIS_CONVERGENCE_FAILURE | Analysis | Algorithm convergence failure |
| ANALYSIS_METHOD_INCOMPATIBLE | Analysis | Method not suitable for data |
| RESOURCE_MEMORY_ERROR | Resource | Insufficient memory |
| RESOURCE_DEPENDENCY_ERROR | Resource | Missing dependencies |
| VIZ_INVALID_PLOT_TYPE | Visualization | Unsupported plot type |
| VIZ_MISSING_DATA | Visualization | Missing required data |

This comprehensive error handling system ensures that ChatSpatial provides clear, actionable feedback to help users and agents recover from failures effectively.