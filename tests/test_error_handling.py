"""
Test error handling and user feedback enhancements
"""

import sys
import os
import asyncio
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
import matplotlib.pyplot as plt

# Add project root directory to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from spatial_transcriptomics_mcp.tools.spatial_analysis import analyze_spatial_unified
from spatial_transcriptomics_mcp.tools.visualization import visualize_data
from spatial_transcriptomics_mcp.models.data import SpatialAnalysisParameters, VisualizationParameters
from mcp.server.fastmcp.utilities.types import Image
from spatial_transcriptomics_mcp.models.analysis import SpatialAnalysisResult
from spatial_transcriptomics_mcp.utils.error_handling import (
    SpatialMCPError, DataNotFoundError, InvalidParameterError, 
    ProcessingError, DataCompatibilityError
)


class MockContext:
    """Mock MCP context"""
    def __init__(self):
        self.logs = []
        self.warnings = []
        self.infos = []
    
    async def info(self, message):
        print(f"INFO: {message}")
        self.infos.append(message)
    
    async def warning(self, message):
        print(f"WARNING: {message}")
        self.warnings.append(message)


async def test_error_handling():
    """Test error handling and user feedback enhancements"""
    print("Creating test data...")

    # Create a simple AnnData object for testing
    np.random.seed(42)
    n_obs = 100
    n_vars = 50

    X = np.random.rand(n_obs, n_vars)
    obs = {"leiden": np.random.choice(["0", "1", "2"], size=n_obs)}
    var = {"gene_name": [f"gene_{i}" for i in range(n_vars)]}

    # Create AnnData object
    adata = ad.AnnData(X, obs=obs, var=var)

    # Convert leiden column to categorical type
    adata.obs["leiden"] = adata.obs["leiden"].astype("category")

    # Add spatial coordinates
    adata.obsm["spatial"] = np.random.rand(n_obs, 2) * 100

    # Add highly variable genes
    adata.var["highly_variable"] = np.random.choice([True, False], size=n_vars, p=[0.2, 0.8])

    # Create data store
    data_store = {"test_data": {"adata": adata, "name": "Test Data"}}

    # Create mock context
    context = MockContext()
    
    # Test 1: Data not found error
    print("\nTest 1: Data not found error")
    try:
        result = await analyze_spatial_unified("non_existent_data", data_store,
                                              SpatialAnalysisParameters(), context)
        print("Error: Should have thrown DataNotFoundError exception")
    except DataNotFoundError as e:
        print(f"Successfully caught DataNotFoundError: {str(e)}")
        print(f"Number of warning messages: {len(context.warnings)}")
    except Exception as e:
        print(f"Error: Caught unexpected exception: {type(e).__name__}: {str(e)}")

    # Test 2: Invalid parameter error
    print("\nTest 2: Invalid parameter error")
    try:
        params = SpatialAnalysisParameters(analysis_type="invalid_type")
        result = await analyze_spatial_unified("test_data", data_store, params, context)
        print("Error: Should have thrown InvalidParameterError exception")
    except InvalidParameterError as e:
        print(f"Successfully caught InvalidParameterError: {str(e)}")
        print(f"Number of warning messages: {len(context.warnings)}")
    except Exception as e:
        print(f"Error: Caught unexpected exception: {type(e).__name__}: {str(e)}")

    # Test 3: Visualization error handling
    print("\nTest 3: Visualization error handling")
    try:
        params = VisualizationParameters(plot_type="invalid_plot")
        result = await visualize_data("test_data", data_store, params, context)
        print("Error: Should have thrown InvalidParameterError exception")
    except InvalidParameterError as e:
        print(f"Successfully caught InvalidParameterError: {str(e)}")
        print(f"Number of warning messages: {len(context.warnings)}")
    except Exception as e:
        print(f"Error: Caught unexpected exception: {type(e).__name__}: {str(e)}")
    
    # Test 4: Handle error but return placeholder image
    print("\nTest 4: Handle error but return placeholder image")

    # Create an AnnData object without spatial coordinates
    adata_no_spatial = ad.AnnData(X, obs=obs, var=var)
    data_store["no_spatial_data"] = {"adata": adata_no_spatial, "name": "No Spatial Data"}

    try:
        params = SpatialAnalysisParameters(analysis_type="neighborhood")
        result = await analyze_spatial_unified("no_spatial_data", data_store, params, context, return_type="image")
        print("Error: Should have thrown DataCompatibilityError exception")
    except DataCompatibilityError as e:
        print(f"Successfully caught DataCompatibilityError: {str(e)}")
        print(f"Number of warning messages: {len(context.warnings)}")
    except Exception as e:
        print(f"Error: Caught unexpected exception: {type(e).__name__}: {str(e)}")

    # Test 5: Test UMAP calculation failure but use PCA as fallback
    print("\nTest 5: Test UMAP calculation failure but use PCA as fallback")

    # Create a scenario that simulates UMAP calculation failure
    # Modify adata to make it fail during UMAP calculation
    adata_bad = adata.copy()
    adata_bad.X = np.zeros((n_obs, n_vars))  # All-zero matrix will cause UMAP calculation to fail
    data_store["bad_data"] = {"adata": adata_bad, "name": "Bad Data"}

    try:
        params = VisualizationParameters(plot_type="umap")
        result = await visualize_data("bad_data", data_store, params, context)
        print(f"Success: Returned placeholder image")
        print(f"Number of warning messages: {len(context.warnings)}")
        print(f"Number of info messages: {len(context.infos)}")

        # Save image to file
        if isinstance(result, Image):
            with open("error_handling_umap_fallback.png", "wb") as f:
                f.write(result.data)
            print("Image saved to error_handling_umap_fallback.png")
    except Exception as e:
        print(f"Error: Caught unexpected exception: {type(e).__name__}: {str(e)}")

    print("\nTesting completed!")


if __name__ == "__main__":
    asyncio.run(test_error_handling())
