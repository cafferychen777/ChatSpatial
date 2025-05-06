"""
测试可视化和图像处理集成
"""

import sys
import os
import unittest
import asyncio
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path

# 添加项目根目录到 Python 路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from spatial_transcriptomics_mcp.tools.visualization import visualize_data
from spatial_transcriptomics_mcp.models.data import VisualizationParameters
from mcp.server.fastmcp.utilities.types import Image


class MockContext:
    """模拟 MCP 上下文"""
    async def info(self, message):
        print(f"INFO: {message}")
    
    async def warning(self, message):
        print(f"WARNING: {message}")


class TestVisualizationIntegration(unittest.TestCase):
    """测试可视化功能与图像处理的集成"""

    def setUp(self):
        """设置测试数据"""
        # 创建一个简单的 AnnData 对象用于测试
        np.random.seed(42)
        n_obs = 100
        n_vars = 50
        
        X = np.random.rand(n_obs, n_vars)
        obs = {"cluster": np.random.choice(["A", "B", "C"], size=n_obs)}
        var = {"gene_name": [f"gene_{i}" for i in range(n_vars)]}
        
        # 创建 AnnData 对象
        self.adata = ad.AnnData(X, obs=obs, var=var)
        
        # 添加空间坐标
        self.adata.obsm["spatial"] = np.random.rand(n_obs, 2) * 100
        
        # 添加 UMAP 坐标
        self.adata.obsm["X_umap"] = np.random.rand(n_obs, 2) * 10
        
        # 添加 leiden 聚类
        self.adata.obs["leiden"] = np.random.choice(["0", "1", "2"], size=n_obs)
        
        # 创建数据存储
        self.data_store = {"test_data": {"adata": self.adata, "name": "Test Data"}}

    def test_visualize_data_returns_image(self):
        """测试 visualize_data 函数返回 Image 对象"""
        # 运行异步测试
        result = asyncio.run(self._test_visualize_data())
        
        # 验证结果是 Image 对象
        self.assertIsInstance(result, Image)
        
        # 验证图像数据不为空
        self.assertIsNotNone(result.data)
        self.assertGreater(len(result.data), 0)

    async def _test_visualize_data(self):
        """异步测试 visualize_data 函数"""
        # 创建参数
        params = VisualizationParameters(
            feature="gene_0",
            plot_type="spatial",
            colormap="viridis"
        )
        
        # 创建模拟上下文
        context = MockContext()
        
        # 调用函数
        result = await visualize_data("test_data", self.data_store, params, context)
        return result


if __name__ == "__main__":
    unittest.main()
