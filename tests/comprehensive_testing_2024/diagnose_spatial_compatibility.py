#!/usr/bin/env python3
"""
诊断空间分析兼容性问题
"""

import sys
import os
from pathlib import Path
import warnings
import traceback

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
warnings.filterwarnings('ignore')

import scanpy as sc
import pandas as pd
import numpy as np

# 尝试导入ChatSpatial工具
try:
    from chatspatial.tools.spatial_analysis import analyze_spatial_patterns
    from chatspatial.models.data import SpatialAnalysisParameters
    SPATIAL_TOOLS_AVAILABLE = True
    print("✓ ChatSpatial spatial tools imported successfully")
except Exception as e:
    SPATIAL_TOOLS_AVAILABLE = False
    print(f"✗ ChatSpatial spatial tools import failed: {e}")
    traceback.print_exc()

def test_spatial_requirements(adata, dataset_name):
    """测试空间分析的具体要求"""
    print(f"\n=== 测试 {dataset_name} ===")
    
    issues = []
    requirements = {
        'basic_loading': False,
        'has_spatial_coords': False,
        'spatial_coords_shape': None,
        'has_clusters': False,
        'cluster_keys': [],
        'has_spatial_neighbors': False,
        'cell_count': adata.n_obs if adata else 0,
        'gene_count': adata.n_vars if adata else 0
    }
    
    try:
        requirements['basic_loading'] = True
        print(f"  ✓ 基本加载: {adata.n_obs} cells, {adata.n_vars} genes")
        
        # 检查空间坐标
        if 'spatial' in adata.obsm:
            requirements['has_spatial_coords'] = True
            requirements['spatial_coords_shape'] = adata.obsm['spatial'].shape
            print(f"  ✓ 空间坐标: {adata.obsm['spatial'].shape}")
            
            # 检查坐标是否合理
            coords = adata.obsm['spatial']
            if np.any(np.isnan(coords)) or np.any(np.isinf(coords)):
                issues.append("空间坐标包含 NaN 或 Inf 值")
                print("  ✗ 空间坐标包含异常值")
            else:
                print("  ✓ 空间坐标数值正常")
        else:
            issues.append("缺少空间坐标 (adata.obsm['spatial'])")
            print("  ✗ 缺少空间坐标")
        
        # 检查聚类信息
        categorical_cols = list(adata.obs.select_dtypes(include=['category', 'object']).columns)
        if categorical_cols:
            requirements['has_clusters'] = True
            requirements['cluster_keys'] = categorical_cols
            print(f"  ✓ 聚类信息: {categorical_cols}")
        else:
            issues.append("缺少聚类信息")
            print("  ✗ 无聚类信息")
            
        # 检查leiden聚类(最常用)
        if 'leiden' in adata.obs.columns:
            print(f"  ✓ Leiden聚类: {adata.obs['leiden'].nunique()} clusters")
        else:
            issues.append("缺少 leiden 聚类")
            print("  ✗ 无 leiden 聚类")
        
        # 检查空间邻居
        if 'spatial_neighbors' in adata.uns:
            requirements['has_spatial_neighbors'] = True  
            print("  ✓ 空间邻居已计算")
        else:
            print("  - 空间邻居未预计算 (可以动态计算)")
            
        # 检查数据规模
        if adata.n_obs < 10:
            issues.append(f"细胞数过少: {adata.n_obs} < 10")
            print(f"  ✗ 细胞数过少: {adata.n_obs}")
        
        if adata.n_vars < 50:
            issues.append(f"基因数过少: {adata.n_vars} < 50")
            print(f"  ✗ 基因数过少: {adata.n_vars}")
            
    except Exception as e:
        issues.append(f"数据检查异常: {str(e)}")
        print(f"  ✗ 检查异常: {e}")
    
    return requirements, issues

def test_actual_spatial_analysis(adata, dataset_name):
    """测试实际的空间分析调用"""
    print(f"\n=== 测试实际空间分析: {dataset_name} ===")
    
    if not SPATIAL_TOOLS_AVAILABLE:
        print("  ✗ 空间分析工具不可用")
        return False, "空间分析工具不可用"
    
    try:
        # 创建测试数据存储
        data_store = {
            'test_data': {'adata': adata}
        }
        
        # 创建参数
        params = SpatialAnalysisParameters(
            analysis_type='moran',
            cluster_key='leiden' if 'leiden' in adata.obs.columns else adata.obs.select_dtypes(include=['category']).columns[0] if len(adata.obs.select_dtypes(include=['category']).columns) > 0 else 'default',
            n_neighbors=15
        )
        
        print(f"  尝试使用参数: cluster_key={params.cluster_key}, analysis_type={params.analysis_type}")
        
        # 这里需要使用 asyncio.run 来运行异步函数
        import asyncio
        
        async def run_analysis():
            return await analyze_spatial_patterns(
                data_id='test_data',
                data_store=data_store,
                params=params,
                context=None
            )
        
        result = asyncio.run(run_analysis())
        print("  ✓ 空间分析成功完成")
        return True, "成功"
        
    except Exception as e:
        error_msg = str(e)
        print(f"  ✗ 空间分析失败: {error_msg}")
        traceback.print_exc()
        return False, error_msg

def main():
    """主函数"""
    print("ChatSpatial 空间分析兼容性诊断")
    print("=" * 50)
    
    # 数据集目录
    datasets_dir = Path("datasets/real_datasets")
    if not datasets_dir.exists():
        datasets_dir = Path("datasets")
    
    # 寻找测试数据集
    test_files = []
    for pattern in ["*.h5ad", "core/*.h5ad", "test/*.h5ad"]:
        test_files.extend(datasets_dir.glob(pattern))
    
    test_files = list(set(test_files))[:5]  # 限制到5个数据集
    
    print(f"找到 {len(test_files)} 个测试数据集")
    
    # 汇总结果
    total_tested = 0
    spatial_analysis_success = 0
    detailed_issues = {}
    
    for dataset_path in test_files:
        if not dataset_path.exists():
            continue
            
        try:
            print(f"\n{'='*60}")
            print(f"测试数据集: {dataset_path.name}")
            print(f"{'='*60}")
            
            # 加载数据集
            adata = sc.read_h5ad(dataset_path)
            
            # 测试空间要求
            requirements, issues = test_spatial_requirements(adata, dataset_path.name)
            
            # 测试实际分析
            success, error_msg = test_actual_spatial_analysis(adata, dataset_path.name)
            
            total_tested += 1
            if success:
                spatial_analysis_success += 1
            
            # 记录详细问题
            detailed_issues[dataset_path.name] = {
                'requirements': requirements,
                'issues': issues,
                'spatial_analysis_success': success,
                'error_message': error_msg if not success else None
            }
            
        except Exception as e:
            print(f"数据集 {dataset_path.name} 测试失败: {e}")
            detailed_issues[dataset_path.name] = {
                'load_error': str(e)
            }
    
    # 生成总结报告
    print(f"\n{'='*60}")
    print("最终诊断报告")
    print(f"{'='*60}")
    
    success_rate = spatial_analysis_success / total_tested * 100 if total_tested > 0 else 0
    print(f"总测试数据集: {total_tested}")
    print(f"空间分析成功: {spatial_analysis_success}")
    print(f"成功率: {success_rate:.1f}%")
    
    # 常见问题统计
    common_issues = {}
    for dataset, info in detailed_issues.items():
        if 'issues' in info:
            for issue in info['issues']:
                if issue in common_issues:
                    common_issues[issue] += 1
                else:
                    common_issues[issue] = 1
    
    if common_issues:
        print(f"\n常见问题:")
        for issue, count in sorted(common_issues.items(), key=lambda x: x[1], reverse=True):
            print(f"  - {issue}: {count} 个数据集")
    
    # 保存详细报告
    import json
    with open(datasets_dir / "spatial_compatibility_diagnosis.json", 'w') as f:
        # 转换numpy类型为可序列化类型
        def convert_types(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, tuple):
                return list(obj)
            return obj
        
        serializable_issues = {}
        for dataset, info in detailed_issues.items():
            serializable_info = {}
            for key, value in info.items():
                if key == 'requirements' and isinstance(value, dict):
                    serializable_info[key] = {k: convert_types(v) for k, v in value.items()}
                else:
                    serializable_info[key] = convert_types(value)
            serializable_issues[dataset] = serializable_info
        
        json.dump({
            'summary': {
                'total_tested': total_tested,
                'spatial_analysis_success': spatial_analysis_success,
                'success_rate': success_rate,
                'common_issues': common_issues
            },
            'detailed_results': serializable_issues
        }, f, indent=2)
    
    print(f"\n详细诊断报告已保存到: {datasets_dir}/spatial_compatibility_diagnosis.json")

if __name__ == "__main__":
    main()