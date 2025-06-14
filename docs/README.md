# ChatSpatial Documentation

## Directory Structure

### `/docs/`
Main documentation directory

- **`analysis_reports/`** - Detailed analysis reports and bug investigation
  - `SPAGCN_DEEP_ANALYSIS.md` - Deep dive into SpaGCN code and performance issues
  - `SPATIAL_DOMAINS_FIXES.md` - Summary of fixes applied to spatial_domains.py
  - `SPATIAL_DOMAINS_FINAL_REPORT.md` - Final report on spatial domains module fixes
  
- **`image_utils.md`** - Documentation for image utility functions
- **`spatial_analysis_enhancement_plan.md`** - Enhancement plans for spatial analysis tools

### `/debug_scripts/`
Debug and validation scripts for development

- `debug_spatial_domains.py` - Debug script for spatial domains functionality
- `validate_fix.py` - Validation script for testing fixes

### `/tests/`
Test suite directory

- **`integration/`** - Integration tests
  - `final_test_spatial_domains.py` - Comprehensive spatial domains test suite
  
- Individual unit test files for each module

### `/examples/`
Example usage scripts demonstrating various functionalities

### `/data/`
Sample data files for testing and examples

### `/debug_scripts/`
Debug and validation scripts for development

### `/third_party/`
Third-party libraries and tools included in the project

## Key Files

- `pyproject.toml` - Project configuration and dependencies
- `README.md` - Main project README
- `chatspatial/` - Main source code directory

## Development Guidelines

1. **Analysis Reports**: Document major bug investigations in `/docs/analysis_reports/`
2. **Debug Scripts**: Place debugging utilities in `/debug_scripts/`
3. **Tests**: Add unit tests to `/tests/` and integration tests to `/tests/integration/`
4. **Examples**: Provide usage examples in `/examples/`

## Recent Major Fixes

- **Spatial Domains Module**: Comprehensive fixes for performance and reliability issues
  - Resolved "mysteriously fails to run" problems
  - Added timeout protection and smart preprocessing
  - Improved error handling and environment compatibility