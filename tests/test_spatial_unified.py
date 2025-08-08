"""
Test unified spatial analysis functionality
"""

import sys
import os
import asyncio
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path

# Add project root directory to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from spatial_transcriptomics_mcp.tools.spatial_analysis import analyze_spatial_unified
from spatial_transcriptomics_mcp.models.data import SpatialAnalysisParameters
from mcp.server.fastmcp.utilities.types import Image
from spatial_transcriptomics_mcp.models.analysis import SpatialAnalysisResult


class MockContext:
    """Mock MCP context"""
    async def info(self, message):
        print(f"INFO: {message}")

    async def warning(self, message):
        print(f"WARNING: {message}")


async def test_spatial_unified():
    """Test unified spatial analysis functionality"""
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

    # Test different spatial analysis types
    analysis_types = ["neighborhood", "co_occurrence"]

    for analysis_type in analysis_types:
        print(f"\nTesting {analysis_type} spatial analysis...")

        # Create spatial analysis parameters
        params = SpatialAnalysisParameters(
            analysis_type=analysis_type,
            n_neighbors=15,
            include_image=True
        )

        # Test returning Image object
        print(f"Testing Image object return...")
        try:
            result = await analyze_spatial_unified("test_data", data_store, params, context, return_type="image")
            print(f"Result type: {type(result)}")
            assert isinstance(result, Image), f"Expected Image, got {type(result)}"
            print(f"Image size: {len(result.data)/1024:.2f} KB")

            # Save image to file
            with open(f"unified_{analysis_type}_image.png", "wb") as f:
                f.write(result.data)
            print(f"Image saved to unified_{analysis_type}_image.png")

        except Exception as e:
            print(f"Error: {str(e)}")

        # Test returning SpatialAnalysisResult object
        print(f"Testing SpatialAnalysisResult object return...")
        try:
            result = await analyze_spatial_unified("test_data", data_store, params, context, return_type="result")
            print(f"Result type: {type(result)}")
            assert isinstance(result, SpatialAnalysisResult), f"Expected SpatialAnalysisResult, got {type(result)}"
            print(f"Analysis type: {result.analysis_type}")
            print(f"Statistics: {result.statistics}")
            print(f"Image: {'Yes' if result.result_image else 'No'}")

            # If there's an image, save to file
            if result.result_image:
                import base64
                img_data = base64.b64decode(result.result_image)
                with open(f"unified_{analysis_type}_result.png", "wb") as f:
                    f.write(img_data)
                print(f"Image saved to unified_{analysis_type}_result.png")

        except Exception as e:
            print(f"Error: {str(e)}")

    print("\nTesting completed!")


if __name__ == "__main__":
    asyncio.run(test_spatial_unified())
