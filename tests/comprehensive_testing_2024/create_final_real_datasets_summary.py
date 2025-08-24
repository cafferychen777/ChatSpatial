#!/usr/bin/env python3
"""
创建真实数据集的最终综合摘要

Linus风格：简单、直接、完整
"""

import os
import pandas as pd
import scanpy as sc
from pathlib import Path
import time
import json
from typing import Dict, List

def scan_all_h5ad_files(base_dir: Path) -> List[Dict]:
    """
    扫描所有h5ad文件并提取信息
    """
    all_datasets = []
    
    print("扫描所有h5ad文件...")
    
    # 递归查找所有h5ad文件
    h5ad_files = list(base_dir.rglob("*.h5ad"))
    
    for h5ad_file in h5ad_files:
        try:
            # 获取相对路径用于分类
            rel_path = h5ad_file.relative_to(base_dir)
            category = str(rel_path.parent) if rel_path.parent != Path(".") else "root"
            
            # 读取基本信息（不加载整个文件）
            try:
                adata = sc.read_h5ad(h5ad_file)
                
                # 提取数据集信息
                dataset_info = adata.uns.get('dataset_info', {})
                
                file_size_mb = h5ad_file.stat().st_size / (1024 * 1024)
                
                # 检查是否有空间坐标
                has_spatial = 'spatial' in adata.obsm
                spatial_dims = adata.obsm['spatial'].shape[1] if has_spatial else 0
                
                # 计算稀疏度
                if hasattr(adata.X, 'toarray'):
                    # 稀疏矩阵
                    sparsity = (adata.X == 0).sum() / adata.X.size * 100
                else:
                    # 密集矩阵
                    sparsity = (adata.X == 0).sum() / adata.X.size * 100
                
                dataset_data = {
                    'file_name': h5ad_file.name,
                    'file_path': str(h5ad_file),
                    'category': category,
                    'name': dataset_info.get('name', h5ad_file.stem),
                    'source': dataset_info.get('source', 'Unknown'),
                    'technology': dataset_info.get('technology', 'Unknown'),
                    'description': dataset_info.get('description', ''),
                    'n_obs': adata.n_obs,
                    'n_vars': adata.n_vars,
                    'file_size_mb': round(file_size_mb, 1),
                    'has_spatial': has_spatial,
                    'spatial_dims': spatial_dims,
                    'sparsity_pct': round(float(sparsity), 1),
                    'obs_columns': list(adata.obs.columns),
                    'var_columns': list(adata.var.columns),
                    'obsm_keys': list(adata.obsm.keys()),
                    'uns_keys': list(adata.uns.keys()),
                    'last_modified': time.strftime('%Y-%m-%d %H:%M:%S', 
                                                  time.localtime(h5ad_file.stat().st_mtime))
                }
                
                all_datasets.append(dataset_data)
                print(f"✓ {h5ad_file.name}: {adata.n_obs} obs, {adata.n_vars} vars, {file_size_mb:.1f}MB")
                
            except Exception as e:
                print(f"✗ 无法读取 {h5ad_file.name}: {e}")
                
                # 至少记录文件基本信息
                file_size_mb = h5ad_file.stat().st_size / (1024 * 1024)
                dataset_data = {
                    'file_name': h5ad_file.name,
                    'file_path': str(h5ad_file),
                    'category': category,
                    'name': h5ad_file.stem,
                    'source': 'Unknown (read error)',
                    'technology': 'Unknown',
                    'description': f'Error reading file: {str(e)[:100]}',
                    'n_obs': 0,
                    'n_vars': 0,
                    'file_size_mb': round(file_size_mb, 1),
                    'has_spatial': False,
                    'spatial_dims': 0,
                    'sparsity_pct': 0,
                    'obs_columns': [],
                    'var_columns': [],
                    'obsm_keys': [],
                    'uns_keys': [],
                    'last_modified': time.strftime('%Y-%m-%d %H:%M:%S', 
                                                  time.localtime(h5ad_file.stat().st_mtime))
                }
                all_datasets.append(dataset_data)
                
        except Exception as e:
            print(f"✗ 处理文件 {h5ad_file} 时出错: {e}")
    
    return all_datasets

def create_comprehensive_summary(datasets: List[Dict], output_dir: Path):
    """
    创建综合摘要报告
    """
    # 转为DataFrame
    df = pd.DataFrame(datasets)
    
    if df.empty:
        print("未找到任何数据集")
        return
    
    # 保存详细CSV
    csv_file = output_dir / "comprehensive_real_datasets_summary.csv"
    df.to_csv(csv_file, index=False)
    print(f"详细摘要保存到: {csv_file}")
    
    # 保存JSON格式
    json_file = output_dir / "comprehensive_real_datasets_summary.json"
    with open(json_file, 'w', encoding='utf-8') as f:
        json.dump(datasets, f, indent=2, ensure_ascii=False)
    print(f"JSON摘要保存到: {json_file}")
    
    # 创建统计报告
    create_statistics_report(df, output_dir)
    
    # 创建使用指南
    create_usage_guide(df, output_dir)

def create_statistics_report(df: pd.DataFrame, output_dir: Path):
    """
    创建统计报告
    """
    stats_file = output_dir / "REAL_DATASETS_STATISTICS.md"
    
    total_datasets = len(df)
    total_size_mb = df['file_size_mb'].sum()
    spatial_datasets = df['has_spatial'].sum()
    
    # 按类别统计
    category_stats = df.groupby('category').agg({
        'file_name': 'count',
        'n_obs': 'sum',
        'n_vars': 'mean',
        'file_size_mb': 'sum'
    }).round(1)
    
    # 按技术统计
    tech_stats = df.groupby('technology').agg({
        'file_name': 'count',
        'n_obs': 'sum',
        'file_size_mb': 'sum'
    }).round(1)
    
    # 按来源统计
    source_stats = df.groupby('source').agg({
        'file_name': 'count',
        'file_size_mb': 'sum'
    }).round(1)
    
    with open(stats_file, 'w', encoding='utf-8') as f:
        f.write("# 真实数据集统计报告\n\n")
        f.write(f"生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## 总体统计\n\n")
        f.write(f"- **总数据集数量**: {total_datasets}\n")
        f.write(f"- **总存储大小**: {total_size_mb:.1f} MB ({total_size_mb/1024:.1f} GB)\n")
        f.write(f"- **空间数据集**: {spatial_datasets} / {total_datasets} ({spatial_datasets/total_datasets*100:.1f}%)\n")
        f.write(f"- **平均文件大小**: {total_size_mb/total_datasets:.1f} MB\n\n")
        
        f.write("## 按类别统计\n\n")
        f.write(category_stats.to_markdown())
        f.write("\n\n")
        
        f.write("## 按技术类型统计\n\n")
        f.write(tech_stats.to_markdown())
        f.write("\n\n")
        
        f.write("## 按数据来源统计\n\n")
        f.write(source_stats.to_markdown())
        f.write("\n\n")
        
        # 最大和最小数据集
        largest = df.loc[df['file_size_mb'].idxmax()]
        smallest = df.loc[df['file_size_mb'].idxmin()]
        most_obs = df.loc[df['n_obs'].idxmax()]
        most_vars = df.loc[df['n_vars'].idxmax()]
        
        f.write("## 数据集极值\n\n")
        f.write(f"- **最大文件**: {largest['file_name']} ({largest['file_size_mb']:.1f} MB)\n")
        f.write(f"- **最小文件**: {smallest['file_name']} ({smallest['file_size_mb']:.1f} MB)\n")
        f.write(f"- **最多观测**: {most_obs['file_name']} ({most_obs['n_obs']} observations)\n")
        f.write(f"- **最多变量**: {most_vars['file_name']} ({most_vars['n_vars']} variables)\n\n")
    
    print(f"统计报告保存到: {stats_file}")

def create_usage_guide(df: pd.DataFrame, output_dir: Path):
    """
    创建使用指南
    """
    guide_file = output_dir / "REAL_DATASETS_USAGE_GUIDE.md"
    
    # 按用途分类推荐数据集
    small_datasets = df[df['file_size_mb'] < 50].sort_values('file_size_mb')
    large_datasets = df[df['file_size_mb'] >= 100].sort_values('file_size_mb', ascending=False)
    spatial_datasets = df[df['has_spatial'] == True].sort_values(['technology', 'n_obs'])
    
    with open(guide_file, 'w', encoding='utf-8') as f:
        f.write("# 真实数据集使用指南\n\n")
        f.write(f"生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## 快速开始推荐\n\n")
        f.write("### 小型测试数据集 (<50MB)\n\n")
        for _, row in small_datasets.head(5).iterrows():
            f.write(f"- **{row['name']}**: {row['n_obs']} spots, {row['n_vars']} genes, {row['file_size_mb']}MB\n")
            f.write(f"  - 技术: {row['technology']}\n")
            f.write(f"  - 来源: {row['source']}\n")
            f.write(f"  - 文件: `{row['file_name']}`\n\n")
        
        f.write("### 大型基准数据集 (>100MB)\n\n")
        for _, row in large_datasets.head(3).iterrows():
            f.write(f"- **{row['name']}**: {row['n_obs']} spots, {row['n_vars']} genes, {row['file_size_mb']}MB\n")
            f.write(f"  - 技术: {row['technology']}\n")
            f.write(f"  - 来源: {row['source']}\n")
            f.write(f"  - 文件: `{row['file_name']}`\n\n")
        
        f.write("### 空间转录组数据集\n\n")
        for _, row in spatial_datasets.head(8).iterrows():
            f.write(f"- **{row['name']}**: {row['n_obs']} spots, {row['n_vars']} genes\n")
            f.write(f"  - 技术: {row['technology']}\n")
            f.write(f"  - 空间维度: {row['spatial_dims']}D\n")
            f.write(f"  - 文件: `{row['file_name']}`\n\n")
        
        f.write("## 使用示例\n\n")
        f.write("```python\n")
        f.write("import scanpy as sc\n")
        f.write("from pathlib import Path\n\n")
        f.write("# 基础目录\n")
        f.write("base_dir = Path('/Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/datasets/real_datasets')\n\n")
        f.write("# 加载小型测试数据集\n")
        if not small_datasets.empty:
            example_file = small_datasets.iloc[0]['file_name']
            f.write(f"adata = sc.read_h5ad(base_dir / '{example_file}')\n")
            f.write("print(f'Loaded: {adata.n_obs} observations, {adata.n_vars} variables')\n\n")
            f.write("# 检查空间坐标\n")
            f.write("if 'spatial' in adata.obsm:\n")
            f.write("    print(f'Spatial coordinates shape: {adata.obsm[\"spatial\"].shape}')\n\n")
        
        f.write("# 批量加载多个数据集\n")
        f.write("def load_datasets_by_category(category):\n")
        f.write("    datasets = []\n")
        f.write("    category_dir = base_dir / category\n")
        f.write("    if category_dir.exists():\n")
        f.write("        for h5ad_file in category_dir.glob('*.h5ad'):\n")
        f.write("            adata = sc.read_h5ad(h5ad_file)\n")
        f.write("            datasets.append(adata)\n")
        f.write("    return datasets\n")
        f.write("```\n\n")
        
        f.write("## 数据集类别说明\n\n")
        categories = df['category'].value_counts()
        for category, count in categories.items():
            f.write(f"- **{category}**: {count} 个数据集\n")
        
        f.write("\n## 技术类型说明\n\n")
        techs = df['technology'].value_counts()
        for tech, count in techs.items():
            f.write(f"- **{tech}**: {count} 个数据集\n")
    
    print(f"使用指南保存到: {guide_file}")

def main():
    print("真实数据集综合摘要生成器")
    print("=" * 50)
    
    base_dir = Path("/Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/datasets/real_datasets")
    
    if not base_dir.exists():
        print(f"错误: 目录不存在 {base_dir}")
        return
    
    # 扫描所有数据集
    all_datasets = scan_all_h5ad_files(base_dir)
    
    if not all_datasets:
        print("未找到任何h5ad文件")
        return
    
    print(f"\n找到 {len(all_datasets)} 个数据集")
    
    # 创建综合摘要
    create_comprehensive_summary(all_datasets, base_dir)
    
    print(f"\n综合摘要生成完成！")
    
    # 显示简要统计
    df = pd.DataFrame(all_datasets)
    total_size = df['file_size_mb'].sum()
    spatial_count = df['has_spatial'].sum()
    
    print(f"总计: {len(all_datasets)} 个数据集, {total_size:.1f}MB")
    print(f"空间数据集: {spatial_count} 个")
    print(f"平均大小: {total_size/len(all_datasets):.1f}MB")

if __name__ == "__main__":
    main()