# Scripts Directory

This directory contains utility scripts for the ChatSpatial project.

## Directory Structure

### `/fixes`
Scripts for fixing compatibility issues with external packages:
- `fix_paste_compatibility.py` - Fixes PASTE package compatibility with POT 0.9.5
- `fix_spatialDE.py` - Fixes SpatialDE compatibility with newer scipy/pandas versions

### `/tests`
Test scripts for validating functionality:
- `test_all_fixed_methods.py` - Tests all three fixed methods (PASTE, SPOTlight, SpatialDE)
- `test_methods_detailed.py` - Detailed testing of newly added spatial methods
- `test_spatialDE_fixed.py` - Specific tests for SpatialDE after fixes
- `test_spotlight_minimal.py` - Minimal test for SPOTlight functionality
- `test_summary.py` - Quick summary test of all methods
- `comprehensive_tests.py` - Comprehensive test suite with various scenarios
- `stress_tests.py` - Stress tests for edge cases and robustness

## Usage

To run fix scripts:
```bash
python scripts/fixes/fix_spatialDE.py
```

To run tests:
```bash
python scripts/tests/test_all_fixed_methods.py
```