# ChatSpatial Comprehensive Error Handling & Data Compatibility Tests

This test suite provides comprehensive testing for error handling and data compatibility in ChatSpatial. It follows Linus Torvalds' engineering principles: robust, practical, and focused on real-world scenarios.

## Overview

The test suite is divided into two main categories:

### 🚨 Error Handling Tests (`error_handling_tests/`)
Tests how the system handles various failure scenarios:

1. **Exceptional Input Handling** - Malformed data, invalid parameters, corrupted files
2. **Boundary Conditions** - Empty data, single cell, minimal/maximal datasets
3. **Resource Limitations** - Memory constraints, disk space issues
4. **Dependency Failures** - Third-party library failures with graceful fallbacks
5. **Network Exceptions** - HTTP mode connectivity issues

### 🔗 Compatibility Tests (`data_compatibility_tests/`)
Tests compatibility across different environments:

1. **Data Format Compatibility** - Different h5ad versions, CSV, HDF5 formats
2. **Parameter Version Compatibility** - Backward compatibility of tool parameters
3. **Platform Compatibility** - macOS, Linux, Windows behavior differences
4. **Dependency Version Compatibility** - Different versions of scanpy, pandas, numpy, etc.

## Test Data

The suite includes special test datasets designed to trigger specific error conditions:

- `empty_dataset.h5ad` - Completely empty dataset (0 cells)
- `single_cell.h5ad` - Dataset with only one cell
- `high_sparsity.h5ad` - Extremely sparse data (>99% zeros)
- `no_spatial.h5ad` - Dataset without spatial coordinates
- `nan_coordinates.h5ad` - Dataset with NaN spatial coordinates
- `identical_coordinates.h5ad` - All cells at same location
- `extreme_coordinates.h5ad` - Infinite/very large coordinate values
- `unicode_characters.h5ad` - Unicode text in metadata
- `mixed_data_types.h5ad` - Mixed data types in annotations

## Quick Start

### Run All Tests
```bash
python run_comprehensive_error_and_compatibility_tests.py
```

### Run Specific Test Categories
```python
# Error handling tests only
from error_handling_tests import ExceptionalInputTests, ErrorTestRunner
test = ExceptionalInputTests()
runner = ErrorTestRunner()
await runner.run_test_suite(test.get_test_cases())

# Compatibility tests only  
from data_compatibility_tests import DataFormatCompatibilityTests, CompatibilityTestRunner
test = DataFormatCompatibilityTests()
runner = CompatibilityTestRunner()
await runner.run_test_suite(test.get_test_cases())
```

### Generate Special Test Data
```bash
python generate_special_test_data.py
```

## Test Framework Architecture

### Error Handling Framework
- **`ErrorTestCase`** - Individual test definition
- **`ErrorCategory`** - Test categorization (exceptional_input, boundary_condition, etc.)
- **`ErrorSeverity`** - Impact level (critical, high, medium, low)
- **`ErrorTestRunner`** - Test execution engine

### Compatibility Framework
- **`CompatibilityTestCase`** - Individual compatibility test
- **`CompatibilityCategory`** - Test type (data_format, platform, etc.)
- **`CompatibilityMatrix`** - Cross-reference compatibility results
- **`CompatibilityTestRunner`** - Compatibility test execution

## Reports Generated

The test suite generates several types of reports:

### 📊 Comprehensive Report (`comprehensive_test_report_TIMESTAMP.json`)
- Overall test summary
- Error classification matrix
- Compatibility matrix  
- Critical issues identification
- Actionable recommendations

### 🔍 Detailed Reports
- `error_handling_detailed_TIMESTAMP.json` - Complete error test results
- `compatibility_detailed_TIMESTAMP.json` - Complete compatibility results
- `test_summary_TIMESTAMP.md` - Human-readable summary

### 📈 Key Metrics
- Overall success rate
- Error handling success rate
- Compatibility success rate
- Critical issues count
- Platform-specific results
- Dependency version compatibility

## Design Philosophy

Following Linus Torvalds' engineering principles:

### 🎯 "Good Taste"
- Eliminate special cases through good data structure design
- Consistent error handling patterns across all components
- Simple, direct test implementations

### 🛡️ "Never Break Userspace"
- Backward compatibility testing for all parameter changes
- Legacy parameter name support verification
- Data format compatibility across versions

### 🔧 Practical Engineering
- Test real-world failure scenarios, not theoretical edge cases
- Focus on production deployment readiness
- Clear, actionable error messages and recommendations

### 🚀 Simplicity
- Minimal abstractions - data and functions, not complex hierarchies
- Direct test implementations without over-engineering
- Clear separation of concerns

## Example Usage

```python
import asyncio
from run_comprehensive_error_and_compatibility_tests import ComprehensiveTestSuite

async def run_tests():
    suite = ComprehensiveTestSuite()
    report = await suite.run_comprehensive_tests()
    
    print(f"Success Rate: {report['overall_summary']['overall_success_rate']}")
    print(f"Critical Issues: {len(report['critical_issues'])}")
    
    return report

# Run tests
report = asyncio.run(run_tests())
```

## Contributing

When adding new tests:

1. **Follow the naming convention**: `test_<specific_scenario>`
2. **Use appropriate categories and severities**
3. **Include clear descriptions** of what each test validates
4. **Test real failure modes**, not artificial scenarios
5. **Generate actionable error messages**

### Adding Error Tests
```python
ErrorTestCase(
    name="test_your_scenario",
    category=ErrorCategory.EXCEPTIONAL_INPUT,
    severity=ErrorSeverity.HIGH,
    description="Clear description of what this tests",
    test_function=your_test_function,
    expected_exception=ExpectedException  # or None
)
```

### Adding Compatibility Tests
```python
CompatibilityTestCase(
    name="test_your_compatibility",
    category=CompatibilityCategory.DATA_FORMAT,
    description="What compatibility aspect this tests",
    test_function=your_test_function,
    version_info={"version": "specific_version"}
)
```

## Dependencies

- `scanpy` - Single-cell analysis
- `pandas` - Data manipulation
- `numpy` - Numerical computing
- `psutil` - System monitoring (for resource tests)
- `requests` - Network testing
- `packaging` - Version comparison

## Output Location

All reports are saved to `./reports/` directory with timestamps to avoid conflicts.

---

*This test suite ensures ChatSpatial is production-ready with robust error handling and broad compatibility.*