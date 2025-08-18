# Archived Tests Directory

This directory contains historical test scripts that were previously scattered in the root directory. These tests are preserved for reference and historical context.

## Directory Structure

### `/scanvi/` - scANVI Method Tests
- **`test_scanvi_comprehensive.py`**: Comprehensive testing of the scANVI annotation method with simulated data
  - Tests basic functionality, error handling, and edge cases
  - Generated the report: `archived_docs/reports/SCANVI_TEST_REPORT.md`
  - Status: Completed with 80% success rate

### `/dependency_validation/` - Dependency System Tests  
- **`test_dependency_validation.py`**: Tests for the improved dependency validation system
  - Validates the new runtime dependency checking mechanism
  - Tests user-friendly error messages and installation guides
  - Status: Validates the dependency detection fixes implemented

### `/tangram/` - Tangram Method Tests
- **`test_tangram_fix.py`**: Tests for Tangram spatial mapping method
- **`debug_tangram.py`**: Debug utilities for Tangram integration issues
  - Related reports in: `archived_docs/reports/TANGRAM_*`
  - Status: Method improvements completed

### `/cellassign/` - CellAssign Method Tests
- **`test_cellassign_comprehensive.py`**: Comprehensive testing of CellAssign annotation method
  - Tests probabilistic cell type assignment
  - Generated results: `archived_docs/reports/cellassign_test_results.json`
  - Status: Method validation completed

## Purpose of Archiving

These test scripts were moved to maintain a clean root directory while preserving:
1. **Historical Context**: Understanding of past testing approaches
2. **Reference Material**: Examples of comprehensive testing strategies
3. **Bug Reports**: Documentation of issues encountered and resolved
4. **Method Validation**: Evidence of successful method implementations

## Current Testing

Active testing should use the organized test suite in `/tests/` directory:
- `/tests/unit/`: Unit tests for individual components
- `/tests/integration/`: Integration tests for workflows
- `/tests/e2e/`: End-to-end MCP integration tests

## Usage Notes

- These archived tests may reference outdated file paths
- Dependencies and imports may need adjustment if running
- For current testing, prefer the organized test suite in `/tests/`
- These files are primarily for reference and documentation purposes