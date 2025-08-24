# ChatSpatial Core Tools Testing Report

Generated: 2025-08-24 04:19:04

## Executive Summary

- **Total Tests**: 13
- **Passed Tests**: 9
- **Failed Tests**: 4
- **Success Rate**: 69.2%
- **Overall Result**: ❌ FAIL

## Dataset Testing Results

### quick_test

**Description**: Quick Visium demo dataset (1K cells)
**Cells**: 1000
**Genes**: 18078

#### Preprocessing Results
✅ **Status**: PASSED
**Operations Completed**:
- qc_metrics: {'original_cells': 1000, 'original_genes': 18078, 'mt_genes_found': 0}
- normalization: {'status': 'completed'}
- highly_variable_genes: {'n_hvg': 2000, 'n_top_genes_requested': 2000}
- pca: {'n_components': 50, 'explained_variance_ratio': 0.21686838567256927}
- neighbors_umap: {'n_neighbors': 10, 'umap_computed': True}
- clustering: {'resolution': 0.4, 'n_clusters': 8}

#### Visualization Results
❌ **Status**: FAILED
**Error**: 1 validation error for VisualizationParameters
save_format
  Extra inputs are not permitted [type=extra_forbidden, input_value='png', input_type=str]
    For further information visit https://errors.pydantic.dev/2.10/v/extra_forbidden

#### Differential Analysis Results
✅ **Status**: PASSED
**Clusters Analyzed**: 8

### performance

**Description**: Performance test synthetic dataset (5K cells)
**Cells**: 5000
**Genes**: 3000

#### Preprocessing Results
✅ **Status**: PASSED
**Operations Completed**:
- qc_metrics: {'original_cells': 5000, 'original_genes': 3000, 'mt_genes_found': 0}
- normalization: {'status': 'completed'}
- highly_variable_genes: {'n_hvg': 1500, 'n_top_genes_requested': 1500}
- pca: {'n_components': 50, 'explained_variance_ratio': 0.01563967391848564}
- neighbors_umap: {'n_neighbors': 15, 'umap_computed': True}
- clustering: {'resolution': 0.6, 'n_clusters': 2}

#### Visualization Results
❌ **Status**: FAILED
**Error**: 1 validation error for VisualizationParameters
save_format
  Extra inputs are not permitted [type=extra_forbidden, input_value='png', input_type=str]
    For further information visit https://errors.pydantic.dev/2.10/v/extra_forbidden

#### Differential Analysis Results
✅ **Status**: PASSED
**Clusters Analyzed**: 2

### standard

**Description**: Standard Slide-seqV2 dataset (42K cells)
**Cells**: 41786
**Genes**: 4000

#### Preprocessing Results
✅ **Status**: PASSED
**Operations Completed**:
- qc_metrics: {'original_cells': 41786, 'original_genes': 4000, 'mt_genes_found': 0}
- normalization: {'status': 'completed'}
- highly_variable_genes: {'n_hvg': 2000, 'n_top_genes_requested': 2000}
- pca: {'n_components': 50, 'explained_variance_ratio': 0.06417469680309296}
- neighbors_umap: {'n_neighbors': 15, 'umap_computed': True}
- clustering: {'resolution': 0.6, 'n_clusters': 12}

#### Visualization Results
❌ **Status**: FAILED
**Error**: 1 validation error for VisualizationParameters
save_format
  Extra inputs are not permitted [type=extra_forbidden, input_value='png', input_type=str]
    For further information visit https://errors.pydantic.dev/2.10/v/extra_forbidden

#### Differential Analysis Results
✅ **Status**: PASSED
**Clusters Analyzed**: 12

### large_scale

**Description**: Large MERFISH dataset (73K cells)
**Cells**: 73655
**Genes**: 161

#### Preprocessing Results
✅ **Status**: PASSED
**Operations Completed**:
- qc_metrics: {'original_cells': 73655, 'original_genes': 161, 'mt_genes_found': 0}
- normalization: {'status': 'completed'}
- highly_variable_genes: {'n_hvg': 80, 'n_top_genes_requested': 80}
- pca: {'n_components': 50, 'explained_variance_ratio': 0.32880889215707754}
- neighbors_umap: {'n_neighbors': 15, 'umap_computed': True}
- clustering: {'resolution': 0.6, 'n_clusters': 16}

#### Visualization Results
❌ **Status**: FAILED
**Error**: 1 validation error for VisualizationParameters
save_format
  Extra inputs are not permitted [type=extra_forbidden, input_value='png', input_type=str]
    For further information visit https://errors.pydantic.dev/2.10/v/extra_forbidden

#### Differential Analysis Results
✅ **Status**: PASSED
**Clusters Analyzed**: 16

## Integration Testing Results

✅ **Status**: PASSED
**Datasets Integrated**: 2
**Common Genes**: 0
**Total Cells**: 6000

## Performance Metrics

| Dataset | Operation | Duration (s) | Memory Used (MB) | Status |
|---------|-----------|--------------|------------------|--------|
| quick_test | data_loading | 0.26 | 181.2 | ✅ |
| quick_test | preprocessing | 9.57 | -708.2 | ✅ |
| quick_test | visualization | 0.00 | 0.0 | ❌ |
| quick_test | differential | 2.74 | 519.2 | ✅ |
| performance | data_loading | 0.03 | 0.0 | ✅ |
| performance | preprocessing | 8.98 | -537.7 | ✅ |
| performance | visualization | 0.00 | 0.0 | ❌ |
| performance | differential | 2.88 | -220.5 | ✅ |
| standard | data_loading | 0.54 | 231.4 | ✅ |
| standard | preprocessing | 66.67 | -899.5 | ✅ |
| standard | visualization | 0.00 | 0.0 | ❌ |
| standard | differential | 7.83 | 374.1 | ✅ |
| large_scale | data_loading | 0.22 | 61.7 | ✅ |
| large_scale | preprocessing | 91.48 | -1628.4 | ✅ |
| large_scale | visualization | 0.00 | 0.0 | ❌ |
| large_scale | differential | 2.93 | 359.3 | ✅ |
| quick_test_performance | integration | 0.81 | 899.6 | ✅ |

## Technical Details

### Environment
- Python version: 3.13.2
- Platform: darwin
- Test output directory: /Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/core_tools_tests/test_results

### Dependencies
- scanpy: 1.11.0
- squidpy: 1.6.2
- matplotlib: 3.10.1

## Recommendations

⚠️ **Needs Attention**: Multiple failures detected, investigation required.

**Failed Operations**: quick_test-visualization, performance-visualization, standard-visualization, large_scale-visualization

---
*Report generated by ChatSpatial Core Tools Tester*
