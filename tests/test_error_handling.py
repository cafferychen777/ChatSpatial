"""
测试错误处理和用户反馈增强
"""

import sys
import os
import asyncio
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
import matplotlib.pyplot as plt

# 添加项目根目录到 Python 路径
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
    """模拟 MCP 上下文"""
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
    """测试错误处理和用户反馈增强"""
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
    
    # 测试 1: 数据不存在错误
    print("\n测试 1: 数据不存在错误")
    try:
        result = await analyze_spatial_unified("non_existent_data", data_store, 
                                              SpatialAnalysisParameters(), context)
        print("错误: 应该抛出 DataNotFoundError 异常")
    except DataNotFoundError as e:
        print(f"成功捕获 DataNotFoundError: {str(e)}")
        print(f"警告消息数量: {len(context.warnings)}")
    except Exception as e:
        print(f"错误: 捕获到意外异常: {type(e).__name__}: {str(e)}")
    
    # 测试 2: 无效参数错误
    print("\n测试 2: 无效参数错误")
    try:
        params = SpatialAnalysisParameters(analysis_type="invalid_type")
        result = await analyze_spatial_unified("test_data", data_store, params, context)
        print("错误: 应该抛出 InvalidParameterError 异常")
    except InvalidParameterError as e:
        print(f"成功捕获 InvalidParameterError: {str(e)}")
        print(f"警告消息数量: {len(context.warnings)}")
    except Exception as e:
        print(f"错误: 捕获到意外异常: {type(e).__name__}: {str(e)}")
    
    # 测试 3: 可视化错误处理
    print("\n测试 3: 可视化错误处理")
    try:
        params = VisualizationParameters(plot_type="invalid_plot")
        result = await visualize_data("test_data", data_store, params, context)
        print("错误: 应该抛出 InvalidParameterError 异常")
    except InvalidParameterError as e:
        print(f"成功捕获 InvalidParameterError: {str(e)}")
        print(f"警告消息数量: {len(context.warnings)}")
    except Exception as e:
        print(f"错误: 捕获到意外异常: {type(e).__name__}: {str(e)}")
    
    # 测试 4: 处理错误但返回占位图像
    print("\n测试 4: 处理错误但返回占位图像")
    
    # 创建一个没有空间坐标的 AnnData 对象
    adata_no_spatial = ad.AnnData(X, obs=obs, var=var)
    data_store["no_spatial_data"] = {"adata": adata_no_spatial, "name": "No Spatial Data"}
    
    try:
        params = SpatialAnalysisParameters(analysis_type="neighborhood")
        result = await analyze_spatial_unified("no_spatial_data", data_store, params, context, return_type="image")
        print("错误: 应该抛出 DataCompatibilityError 异常")
    except DataCompatibilityError as e:
        print(f"成功捕获 DataCompatibilityError: {str(e)}")
        print(f"警告消息数量: {len(context.warnings)}")
    except Exception as e:
        print(f"错误: 捕获到意外异常: {type(e).__name__}: {str(e)}")
    
    # 测试 5: 测试 UMAP 计算失败但使用 PCA 作为回退
    print("\n测试 5: 测试 UMAP 计算失败但使用 PCA 作为回退")
    
    # 创建一个模拟 UMAP 计算失败的情况
    # 通过修改 adata 使其在计算 UMAP 时失败
    adata_bad = adata.copy()
    adata_bad.X = np.zeros((n_obs, n_vars))  # 全零矩阵会导致 UMAP 计算失败
    data_store["bad_data"] = {"adata": adata_bad, "name": "Bad Data"}
    
    try:
        params = VisualizationParameters(plot_type="umap")
        result = await visualize_data("bad_data", data_store, params, context)
        print(f"成功: 返回了占位图像")
        print(f"警告消息数量: {len(context.warnings)}")
        print(f"信息消息数量: {len(context.infos)}")
        
        # 保存图像到文件
        if isinstance(result, Image):
            with open("error_handling_umap_fallback.png", "wb") as f:
                f.write(result.data)
            print("图像已保存到 error_handling_umap_fallback.png")
    except Exception as e:
        print(f"错误: 捕获到意外异常: {type(e).__name__}: {str(e)}")
    
    print("\n测试完成！")


if __name__ == "__main__":
    asyncio.run(test_error_handling())
