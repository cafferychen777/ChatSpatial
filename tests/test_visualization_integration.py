"""
Test visualization and image processing integration
"""

import sys
import os
import unittest
import asyncio
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path

# Add project root directory to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from spatial_transcriptomics_mcp.tools.visualization import visualize_data
from spatial_transcriptomics_mcp.models.data import VisualizationParameters
from mcp.server.fastmcp.utilities.types import Image


class MockContext:
    """Mock MCP context"""
    async def info(self, message):
        print(f"INFO: {message}")
    
    async def warning(self, message):
        print(f"WARNING: {message}")


class TestVisualizationIntegration(unittest.TestCase):
    """Testvisualization功能与imageprocessing的集成"""

    def setUp(self):
        """setupTestdata"""
        # create一个简单的 AnnData object用于Test
        np.random.seed(42)
        n_obs = 100
        n_vars = 50
        
        X = np.random.rand(n_obs, n_vars)
        obs = {"cluster": np.random.choice(["A", "B", "C"], size=n_obs)}
        var = {"gene_name": [f"gene_{i}" for i in range(n_vars)]}
        
        # create AnnData object
        self.adata = ad.AnnData(X, obs=obs, var=var)
        
        # 添加spatialcoordinate
        self.adata.obsm["spatial"] = np.random.rand(n_obs, 2) * 100
        
        # 添加 UMAP coordinate
        self.adata.obsm["X_umap"] = np.random.rand(n_obs, 2) * 10
        
        # 添加 leiden 聚class
        self.adata.obs["leiden"] = np.random.choice(["0", "1", "2"], size=n_obs)
        
        # createdatastorage
        self.data_store = {"test_data": {"adata": self.adata, "name": "Test Data"}}

    def test_visualize_data_returns_image(self):
        """Test visualize_data functionreturn Image object"""
        # run异步Test
        result = asyncio.run(self._test_visualize_data())
        
        # validationresult是 Image object
        self.assertIsInstance(result, Image)
        
        # validationimagedata不为空
        self.assertIsNotNone(result.data)
        self.assertGreater(len(result.data), 0)

    async def _test_visualize_data(self):
        """异步Test visualize_data function"""
        # createparameter
        params = VisualizationParameters(
            feature="gene_0",
            plot_type="spatial",
            colormap="viridis"
        )
        
        # create模拟上下文
        context = MockContext()
        
        # callfunction
        result = await visualize_data("test_data", self.data_store, params, context)
        return result


if __name__ == "__main__":
    unittest.main()
