#!/usr/bin/env python3
"""
Simple test for preprocessing improvements
"""

import asyncio
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from chatspatial.models.data import AnalysisParameters
from chatspatial.tools.preprocessing import preprocess_data


async def simple_test():
    """Simple test"""
    print("🧪 简单改进测试...")
    
    # 创建小数据集 
    np.random.seed(42)
    n_cells, n_genes = 50, 100
    X = np.random.poisson(5, (n_cells, n_genes)).astype(np.float32)
    
    adata = AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell_{i}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])
    )
    
    # 添加空间坐标
    coords = np.random.normal(0, 5, (n_cells, 2))
    adata.obsm['spatial'] = coords
    
    # 创建数据存储
    data_store = {"test": {"adata": adata}}
    
    try:
        # 测试1: 基本功能
        print("  测试1: 基本功能...")
        params = AnalysisParameters(n_hvgs=20)
        result = await preprocess_data("test", data_store, params)
        print(f"    ✅ 基本功能成功: result type = {type(result)}")
        if isinstance(result, dict):
            print(f"    Dict keys: {result.keys()}")
        else:
            print(f"    {result.n_cells}细胞, {result.clusters}聚类")
        
        # 测试2: 用户控制参数
        print("  测试2: 用户控制参数...")
        params = AnalysisParameters(n_neighbors=5, clustering_resolution=0.3, n_hvgs=20)
        result = await preprocess_data("test", data_store, params)
        print(f"    ✅ 用户参数成功: {result.clusters}聚类")
        
        print(f"\n🎉 简单测试成功!")
        return True
        
    except Exception as e:
        print(f"    ❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = asyncio.run(simple_test())
    print(f"结果: {'成功' if success else '失败'}")
    exit(0 if success else 1)