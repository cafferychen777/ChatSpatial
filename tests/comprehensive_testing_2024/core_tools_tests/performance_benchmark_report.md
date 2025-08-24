# ChatSpatial Core Tools - Performance Benchmark Report
Generated: 2025-08-24 04:33:29.341850
Total measurements: 21

## Overall Performance Summary
- **Successful tests**: 17/21 (81.0%)
- **Average execution time**: 0.805s
- **Average memory usage**: 732.6MB
- **Average memory delta**: 11.3MB

## Import Performance
- Tests: 5/5 successful
- Avg time: 1.744s
- Max time: 5.570s
- Avg memory: 38.4MB
  - `import_cell_communication`: 2.052s, 5.9MB
  - `import_annotation`: 0.320s, 0.1MB
  - `import_deconvolution`: 5.570s, 185.5MB
  - `import_visualization`: 0.382s, 0.2MB
  - `import_spatial_analysis`: 0.395s, 0.0MB

## Create Performance
- Tests: 4/8 successful
- Avg time: 0.296s
- Max time: 0.384s
- Avg memory: 0.0MB
  - `create_params_CellCommunicationParameters`: 0.272s, 0.0MB
  - `create_params_AnnotationParameters`: 0.384s, 0.0MB
  - `create_params_DeconvolutionParameters`: 0.270s, 0.0MB
  - `create_params_SpatialAnalysisParameters`: 0.256s, 0.0MB
  **Failures**:
  - `create_synthetic_data_tiny`: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
  - `create_synthetic_data_small`: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
  - `create_synthetic_data_medium`: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
  - `create_synthetic_data_large`: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'

## Check Performance
- Tests: 8/8 successful
- Avg time: 0.473s
- Max time: 0.703s
- Avg memory: -0.0MB
  - `check_dependency_liana`: 0.398s, 0.0MB
  - `check_dependency_scvi`: 0.588s, 0.0MB
  - `check_dependency_tangram`: 0.703s, -0.3MB
  - `check_dependency_cellphonedb`: 0.484s, 0.0MB
  - `check_dependency_squidpy`: 0.412s, 0.0MB
  - `check_dependency_scanpy`: 0.397s, 0.0MB
  - `check_dependency_pandas`: 0.426s, 0.0MB
  - `check_dependency_numpy`: 0.374s, 0.0MB

## Data Performance
- Tests: 0/4 successful
  **Failures**:
  - `create_synthetic_data_tiny`: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
  - `create_synthetic_data_small`: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
  - `create_synthetic_data_medium`: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
  - `create_synthetic_data_large`: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'

## Dataset Scalability
## Performance Recommendations
### High Memory Usage Detected
- `import_deconvolution`: 185.5MB - Consider optimization
### Slow Operations Detected
- `import_deconvolution`: 5.57s - Consider caching or optimization