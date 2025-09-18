# ChatSpatial Test Suite

This directory contains the official test suite for ChatSpatial MCP server.

## Directory Structure

```text
tests/
├── README.md          # This document
├── e2e/              # End-to-end test plans and documentation
├── fixtures/         # Test fixtures and mock data
└── integration/      # Integration test documentation
```

## Running Tests

Use pytest to run the test suite:

```bash
cd /Users/apple/Research/SpatialTrans_MCP/chatspatial
pytest tests/
```

## Directory Details

### e2e/
Contains comprehensive end-to-end test plans covering:
- Basic workflow validation
- Claude Desktop integration scenarios  
- Progress reporting functionality
- Error handling verification
- Performance benchmarks
- Cross-platform compatibility

### fixtures/
Provides standardized test data and utilities:
- Synthetic spatial transcriptomics datasets
- Mock data stores for testing
- Utility functions for test validation
- Various data quality scenarios

### integration/
Integration test documentation for:
- Trajectory analysis methods
- Enrichment analysis workflows
- Multi-component system testing

## Dependencies

Ensure required dependencies are installed:

```bash
# Base dependencies
pip install -r requirements.txt

# Full dependencies (including advanced features)  
pip install -r requirements-full.txt
```