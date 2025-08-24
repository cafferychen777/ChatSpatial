# ChatSpatial Core Tools - Simplified Test Report
Generated: 2025-08-24 04:28:01.691416

## Summary
- **Total Tests**: 5
- **Passed**: 3 ✅
- **Failed**: 2 ❌
- **Skipped**: 0 ⏭️
- **Success Rate**: 60.0%
- **Execution Time**: 4.85s

## Test Details
### imports_and_structure ✅
- **Status**: PASSED
- **Time**: 4.85s
- **message**: All core modules imported successfully

### parameter_validation ✅
- **Status**: PASSED
- **Time**: 4.85s
- **comm_method**: liana
- **annot_method**: marker_genes
- **deconv_method**: cell2location

### data_validation_system ❌
- **Status**: FAILED
- **Time**: 4.85s
- **error**: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'

### dependency_detection ✅
- **Status**: PASSED
- **Time**: 4.85s
- **dependencies**:
  - liana: {'available': True, 'version': '0.1.9', 'api_modules': ['fun', 'funcomics', 'method', 'mt', 'mu', 'multi', 'pl', 'plotting', 'resource', 'rs', 'testing']}
  - scvi: {'available': True, 'version': '1.3.2'}
  - cellphonedb: {'available': False}
  - tangram: {'available': True}

### basic_cell_communication ❌
- **Status**: FAILED
- **Time**: 4.85s
- **error**: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'

## Analysis & Recommendations
### Issues Found
- **data_validation_system**: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
- **basic_cell_communication**: Cannot cast ufunc 'add' output from dtype('float64') to dtype('int64') with casting rule 'same_kind'
### Dependency Status
- **liana**: ✅ Available (v0.1.9)
- **scvi**: ✅ Available (v1.3.2)
- **cellphonedb**: ❌ Not installed
- **tangram**: ✅ Available (vunknown)