# ChatSpatial Core Tools Test Report
Generated: 2025-08-24 04:23:59.464224

## Test Summary
- **Total Tests**: 3
- **Passed**: 0 ✅
- **Failed**: 3 ❌
- **Skipped**: 0 ⏭️
- **Success Rate**: 0.0%
- **Total Execution Time**: 4.61 seconds

## Cell Communication Results
- **liana_integration**: ❌ FAILED
  - Dataset: squidpy_merfish
  - Error: Error in cell communication analysis: LIANA+ analysis failed: module 'liana' has no attribute 'ut'
- **cellphonedb_integration**: ❌ FAILED
  - Dataset: squidpy_slideseqv2
  - Error: Test exceeded memory limit: 1126.16MB > 1000MB

## Performance Notes
- Tests use reduced parameters for speed (fewer epochs, permutations)
- Large datasets are automatically subsetted for testing
- GPU usage is disabled to ensure reproducibility
- Memory usage is monitored to prevent system issues

## Recommendations
### Failed Tests
- **liana_integration**: Check dependencies and data format
- **cellphonedb_integration**: Check dependencies and data format