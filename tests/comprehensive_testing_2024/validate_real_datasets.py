#!/usr/bin/env python3
"""
验证真实空间转录组数据集的完整性和可用性
确保所有数据集都符合测试要求
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
import json
import traceback
from typing import Dict, List, Optional, Tuple

warnings.filterwarnings('ignore')

# 设置数据目录
REAL_DATASETS_DIR = Path(__file__).parent / 'datasets' / 'real_datasets'

def validate_single_dataset(h5ad_file: Path) -> Dict:
    """验证单个数据集"""
    import scanpy as sc
    
    validation_result = {
        'file': h5ad_file.name,
        'file_path': str(h5ad_file),
        'valid': False,
        'errors': [],
        'warnings': [],
        'metadata': {}
    }
    
    try:
        # 1. 文件读取测试
        adata = sc.read_h5ad(h5ad_file)
        validation_result['metadata']['readable'] = True
        
        # 2. 基本结构验证
        validation_result['metadata']['n_cells'] = adata.n_obs
        validation_result['metadata']['n_genes'] = adata.n_vars
        validation_result['metadata']['file_size_mb'] = round(h5ad_file.stat().st_size / 1024 / 1024, 2)
        
        # 3. 数据类型验证
        if hasattr(adata.X, 'dtype'):
            validation_result['metadata']['X_dtype'] = str(adata.X.dtype)
        
        # 4. 空间坐标验证
        has_spatial = 'spatial' in adata.obsm
        validation_result['metadata']['has_spatial'] = has_spatial
        
        if has_spatial:
            spatial_coords = adata.obsm['spatial']
            validation_result['metadata']['spatial_shape'] = list(spatial_coords.shape)
            validation_result['metadata']['spatial_dims'] = spatial_coords.shape[1]
            
            # 检查空间坐标的有效性
            if np.any(np.isnan(spatial_coords)):
                validation_result['warnings'].append('空间坐标包含NaN值')
            
            if np.all(spatial_coords == 0):
                validation_result['warnings'].append('所有空间坐标都为0')
            
            # 检查坐标范围
            coord_ranges = {
                'x_min': float(spatial_coords[:, 0].min()),
                'x_max': float(spatial_coords[:, 0].max()),
                'y_min': float(spatial_coords[:, 1].min()),
                'y_max': float(spatial_coords[:, 1].max()),
            }
            validation_result['metadata']['coordinate_ranges'] = coord_ranges
        else:
            validation_result['warnings'].append('缺少空间坐标信息')
        
        # 5. 表达数据验证
        if adata.X.shape[0] == 0 or adata.X.shape[1] == 0:
            validation_result['errors'].append('表达矩阵为空')
        
        # 检查稀疏度
        if hasattr(adata.X, 'todense'):
            # 稀疏矩阵
            sparsity = 1 - adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1])
            validation_result['metadata']['sparsity'] = round(sparsity, 4)
            validation_result['metadata']['matrix_type'] = 'sparse'
        else:
            # 密集矩阵
            sparsity = (adata.X == 0).sum() / adata.X.size
            validation_result['metadata']['sparsity'] = round(float(sparsity), 4)
            validation_result['metadata']['matrix_type'] = 'dense'
        
        # 6. 元数据验证
        obs_columns = list(adata.obs.columns)
        var_columns = list(adata.var.columns)
        validation_result['metadata']['obs_columns'] = obs_columns
        validation_result['metadata']['var_columns'] = var_columns
        validation_result['metadata']['n_obs_columns'] = len(obs_columns)
        validation_result['metadata']['n_var_columns'] = len(var_columns)
        
        # 7. obsm和varm检查
        obsm_keys = list(adata.obsm.keys())
        varm_keys = list(adata.varm.keys())
        validation_result['metadata']['obsm_keys'] = obsm_keys
        validation_result['metadata']['varm_keys'] = varm_keys
        
        # 8. uns检查
        uns_keys = list(adata.uns.keys())
        validation_result['metadata']['uns_keys'] = uns_keys
        
        # 9. 数据质量检查
        if adata.n_obs < 10:
            validation_result['warnings'].append(f'细胞数量过少: {adata.n_obs}')
        
        if adata.n_vars < 10:
            validation_result['warnings'].append(f'基因数量过少: {adata.n_vars}')
        
        if validation_result['metadata']['sparsity'] > 0.99:
            validation_result['warnings'].append(f'数据极度稀疏: {validation_result["metadata"]["sparsity"]:.3f}')
        
        # 如果没有严重错误，标记为有效
        if not validation_result['errors']:
            validation_result['valid'] = True
            
    except Exception as e:
        validation_result['errors'].append(f'读取失败: {str(e)}')
        validation_result['metadata']['exception'] = str(e)
    
    return validation_result

def validate_all_datasets() -> Dict:
    """验证所有数据集"""
    print("开始验证真实空间转录组数据集...")
    print(f"数据目录: {REAL_DATASETS_DIR}")
    
    # 查找所有h5ad文件
    h5ad_files = list(REAL_DATASETS_DIR.rglob('*.h5ad'))
    h5ad_files.sort()
    
    print(f"找到 {len(h5ad_files)} 个数据集文件")
    
    validation_results = {
        'total_files': len(h5ad_files),
        'valid_files': 0,
        'invalid_files': 0,
        'files_with_warnings': 0,
        'total_cells': 0,
        'total_genes_unique': set(),
        'total_size_mb': 0,
        'validation_timestamp': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
        'individual_results': []
    }
    
    # 验证每个文件
    for i, h5ad_file in enumerate(h5ad_files, 1):
        print(f"\n[{i:2d}/{len(h5ad_files)}] 验证 {h5ad_file.name}...")
        
        result = validate_single_dataset(h5ad_file)
        validation_results['individual_results'].append(result)
        
        # 更新统计
        if result['valid']:
            validation_results['valid_files'] += 1
            print(f"    ✓ 有效 - {result['metadata'].get('n_cells', 0)} cells, {result['metadata'].get('n_genes', 0)} genes")
            
            # 累加统计信息
            validation_results['total_cells'] += result['metadata'].get('n_cells', 0)
            validation_results['total_size_mb'] += result['metadata'].get('file_size_mb', 0)
            
        else:
            validation_results['invalid_files'] += 1
            print(f"    ❌ 无效")
            
        # 显示错误和警告
        for error in result['errors']:
            print(f"    🚫 错误: {error}")
            
        for warning in result['warnings']:
            print(f"    ⚠️  警告: {warning}")
            
        if result['warnings']:
            validation_results['files_with_warnings'] += 1
    
    return validation_results

def generate_validation_report(validation_results: Dict) -> None:
    """生成验证报告"""
    report_file = REAL_DATASETS_DIR / 'validation_report.json'
    
    # 保存详细结果
    with open(report_file, 'w', encoding='utf-8') as f:
        json.dump(validation_results, f, indent=2, ensure_ascii=False, default=str)
    
    print(f"\n=== 验证报告 ===")
    print(f"报告保存至: {report_file}")
    print(f"验证时间: {validation_results['validation_timestamp']}")
    print(f"总文件数: {validation_results['total_files']}")
    print(f"有效文件: {validation_results['valid_files']}")
    print(f"无效文件: {validation_results['invalid_files']}")
    print(f"有警告文件: {validation_results['files_with_warnings']}")
    print(f"总细胞数: {validation_results['total_cells']:,}")
    print(f"总数据大小: {validation_results['total_size_mb']:.1f} MB")
    
    # 生成简化摘要
    summary_data = []
    for result in validation_results['individual_results']:
        if result['valid']:
            summary_data.append({
                'dataset': result['file'],
                'status': '✓ 有效',
                'n_cells': result['metadata'].get('n_cells', 'N/A'),
                'n_genes': result['metadata'].get('n_genes', 'N/A'),
                'has_spatial': result['metadata'].get('has_spatial', False),
                'sparsity': result['metadata'].get('sparsity', 'N/A'),
                'file_size_mb': result['metadata'].get('file_size_mb', 'N/A'),
                'warnings': len(result['warnings']),
                'errors': len(result['errors'])
            })
        else:
            summary_data.append({
                'dataset': result['file'],
                'status': '❌ 无效',
                'n_cells': 'N/A',
                'n_genes': 'N/A',
                'has_spatial': False,
                'sparsity': 'N/A',
                'file_size_mb': result['metadata'].get('file_size_mb', 'N/A'),
                'warnings': len(result['warnings']),
                'errors': len(result['errors'])
            })
    
    summary_df = pd.DataFrame(summary_data)
    summary_csv = REAL_DATASETS_DIR / 'validation_summary.csv'
    summary_df.to_csv(summary_csv, index=False)
    print(f"验证摘要保存至: {summary_csv}")
    
    # 显示问题文件
    problem_files = [r for r in validation_results['individual_results'] 
                    if not r['valid'] or r['warnings']]
    
    if problem_files:
        print(f"\n=== 需要关注的文件 ===")
        for result in problem_files:
            print(f"\n{result['file']}:")
            for error in result['errors']:
                print(f"  🚫 {error}")
            for warning in result['warnings']:
                print(f"  ⚠️  {warning}")
    else:
        print(f"\n🎉 所有数据集都通过验证！")

def main():
    """主函数"""
    try:
        import scanpy as sc
        print(f"✓ scanpy版本: {sc.__version__}")
    except ImportError:
        print("❌ scanpy未安装，无法进行验证")
        print("请安装: pip install scanpy")
        return
    
    if not REAL_DATASETS_DIR.exists():
        print(f"❌ 数据目录不存在: {REAL_DATASETS_DIR}")
        return
    
    try:
        # 执行验证
        validation_results = validate_all_datasets()
        
        # 生成报告
        generate_validation_report(validation_results)
        
    except KeyboardInterrupt:
        print("\n❌ 用户中断验证")
    except Exception as e:
        print(f"\n❌ 验证过程中发生错误: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    main()