# Integration Tests

This directory contains integration tests for the ChatSpatial MCP server functionality.

## Test Files

### test_trajectory_analysis.py
Tests for trajectory inference methods:
- CellRank (with and without PETSc/SLEPc)
- Palantir
- DPT (Diffusion Pseudotime)
- RNA velocity integration

### test_enrichment_analysis.py
Tests for enrichment analysis methods:
- GSEA (Gene Set Enrichment Analysis)
- ORA (Over-Representation Analysis)
- ssGSEA (Single Sample GSEA)
- EnrichMap (spatial enrichment)

## Running Tests

To run all integration tests:
```bash
cd /Users/apple/Research/SpatialTrans_MCP/chatspatial
python -m pytest tests/integration/
```

To run a specific test:
```bash
python tests/integration/test_trajectory_analysis.py
```

## Notes

- These tests use the pancreas velocity dataset
- Some tests require optional dependencies (e.g., gseapy for enrichment)
- Tests are designed to be run independently