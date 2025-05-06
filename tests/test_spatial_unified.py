"""
测试统一的空间分析功能
"""

import sys
import os
import asyncio
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path

# 添加项目根目录到 Python 路径
sys.path.insert(0, str(Path(__file__).parent.parent))

from spatial_transcriptomics_mcp.tools.spatial_analysis import analyze_spatial_unified
from spatial_transcriptomics_mcp.models.data import SpatialAnalysisParameters
from mcp.server.fastmcp.utilities.types import Image
from spatial_transcriptomics_mcp.models.analysis import SpatialAnalysisResult


class MockContext:
    """模拟 MCP 上下文"""
    async def info(self, message):
        print(f"INFO: {message}")

    async def warning(self, message):
        print(f"WARNING: {message}")


async def test_spatial_unified():
    """测试统一的空间分析功能"""
    print("创建测试数据...")

    # 创建一个简单的 AnnData 对象用于测试
    np.random.seed(42)
    n_obs = 100
    n_vars = 50

    X = np.random.rand(n_obs, n_vars)
    obs = {"leiden": np.random.choice(["0", "1", "2"], size=n_obs)}
    var = {"gene_name": [f"gene_{i}" for i in range(n_vars)]}

    # 创建 AnnData 对象
    adata = ad.AnnData(X, obs=obs, var=var)

    # 将leiden列转换为分类类型
    adata.obs["leiden"] = adata.obs["leiden"].astype("category")

    # 添加空间坐标
    adata.obsm["spatial"] = np.random.rand(n_obs, 2) * 100

    # 添加高变基因
    adata.var["highly_variable"] = np.random.choice([True, False], size=n_vars, p=[0.2, 0.8])

    # 创建数据存储
    data_store = {"test_data": {"adata": adata, "name": "Test Data"}}

    # 创建模拟上下文
    context = MockContext()

    # 测试不同的空间分析类型
    analysis_types = ["neighborhood", "co_occurrence"]

    for analysis_type in analysis_types:
        print(f"\n测试 {analysis_type} 空间分析...")

        # 创建空间分析参数
        params = SpatialAnalysisParameters(
            analysis_type=analysis_type,
            n_neighbors=15,
            include_image=True
        )

        # 测试返回 Image 对象
        print(f"测试返回 Image 对象...")
        try:
            result = await analyze_spatial_unified("test_data", data_store, params, context, return_type="image")
            print(f"结果类型: {type(result)}")
            assert isinstance(result, Image), f"Expected Image, got {type(result)}"
            print(f"图像大小: {len(result.data)/1024:.2f} KB")

            # 保存图像到文件
            with open(f"unified_{analysis_type}_image.png", "wb") as f:
                f.write(result.data)
            print(f"图像已保存到 unified_{analysis_type}_image.png")

        except Exception as e:
            print(f"错误: {str(e)}")

        # 测试返回 SpatialAnalysisResult 对象
        print(f"测试返回 SpatialAnalysisResult 对象...")
        try:
            result = await analyze_spatial_unified("test_data", data_store, params, context, return_type="result")
            print(f"结果类型: {type(result)}")
            assert isinstance(result, SpatialAnalysisResult), f"Expected SpatialAnalysisResult, got {type(result)}"
            print(f"分析类型: {result.analysis_type}")
            print(f"统计信息: {result.statistics}")
            print(f"图像: {'有' if result.result_image else '无'}")

            # 如果有图像，保存到文件
            if result.result_image:
                import base64
                img_data = base64.b64decode(result.result_image)
                with open(f"unified_{analysis_type}_result.png", "wb") as f:
                    f.write(img_data)
                print(f"图像已保存到 unified_{analysis_type}_result.png")

        except Exception as e:
            print(f"错误: {str(e)}")

    print("\n测试完成！")


if __name__ == "__main__":
    asyncio.run(test_spatial_unified())
