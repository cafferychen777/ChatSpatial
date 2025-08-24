#!/usr/bin/env python3
"""
整理所有数据集到统一的data目录结构中
"""

import os
import shutil
from pathlib import Path
import json
import scanpy as sc
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

def create_data_structure():
    """创建标准的data目录结构"""
    base_dir = Path("/Users/apple/Research/SpatialTrans_MCP/chatspatial/data")
    
    # 创建标准目录结构
    directories = {
        'real_datasets': '真实的空间转录组数据集',
        'synthetic_datasets': '合成/模拟数据集',
        'benchmark_datasets': '性能基准测试数据集',
        'demo_datasets': '演示和教程数据集',
        'reference_datasets': '参考和标准数据集',
        'metadata': '数据集元信息和索引'
    }
    
    for dir_name, description in directories.items():
        dir_path = base_dir / dir_name
        dir_path.mkdir(exist_ok=True)
        
        # 创建README文件
        readme_path = dir_path / 'README.md'
        if not readme_path.exists():
            with open(readme_path, 'w') as f:
                f.write(f"# {dir_name.replace('_', ' ').title()}\n\n")
                f.write(f"{description}\n\n")
                f.write("## 数据集列表\n\n")
                f.write("更新时间: 待更新\n")
    
    return base_dir, directories

def collect_all_datasets():
    """收集所有现有的数据集文件"""
    datasets = []
    
    # 搜索路径
    search_paths = [
        Path("/Users/apple/Research/SpatialTrans_MCP/chatspatial/data"),
        Path("/Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024"),
    ]
    
    for search_path in search_paths:
        if search_path.exists():
            # 递归查找所有h5ad文件
            for h5ad_file in search_path.rglob("*.h5ad"):
                if h5ad_file.is_file():
                    datasets.append({
                        'current_path': h5ad_file,
                        'filename': h5ad_file.name,
                        'size_mb': h5ad_file.stat().st_size / (1024 * 1024),
                        'source_dir': h5ad_file.parent.name
                    })
    
    print(f"找到 {len(datasets)} 个数据集文件")
    return datasets

def categorize_dataset(dataset_info):
    """根据文件名和路径判断数据集类别"""
    filename = dataset_info['filename'].lower()
    source_dir = dataset_info['source_dir'].lower()
    
    # 分类规则
    if any(x in filename for x in ['synthetic', 'simulated', 'fake', 'artificial']):
        return 'synthetic_datasets'
    elif any(x in filename for x in ['benchmark', 'perf_test', 'stress']):
        return 'benchmark_datasets'  
    elif any(x in filename for x in ['demo', 'tutorial', 'example', 'quick']):
        return 'demo_datasets'
    elif any(x in filename for x in ['paul15', 'pancreas_tiny', 'pancreas_subset']):
        return 'reference_datasets'
    elif any(x in filename for x in ['squidpy', 'slideseq', 'visium', 'merfish', 'seqfish', 'osmfish', 'starmap']):
        return 'real_datasets'
    elif 'harmony' in source_dir:
        return 'demo_datasets'  # harmony数据主要用于演示
    elif 'test' in source_dir or 'benchmark' in source_dir:
        return 'benchmark_datasets'
    else:
        return 'real_datasets'  # 默认归类为真实数据

def analyze_dataset(file_path):
    """分析单个数据集的详细信息"""
    try:
        adata = sc.read_h5ad(file_path)
        
        info = {
            'filename': file_path.name,
            'n_cells': adata.n_obs,
            'n_genes': adata.n_vars,
            'has_spatial': 'spatial' in adata.obsm,
            'spatial_dims': adata.obsm['spatial'].shape if 'spatial' in adata.obsm else None,
            'has_clusters': len(adata.obs.select_dtypes(include=['category', 'object']).columns) > 0,
            'cluster_keys': list(adata.obs.select_dtypes(include=['category', 'object']).columns),
            'has_leiden': 'leiden' in adata.obs.columns,
            'sparsity': float((adata.X == 0).sum() / adata.X.size) if hasattr(adata.X, 'size') else 0.0,
            'file_size_mb': float(file_path.stat().st_size / (1024 * 1024)),
            'memory_usage_mb': float(adata.X.data.nbytes / (1024 * 1024)) if hasattr(adata.X, 'data') else 0.0
        }
        
        # 检测可能的技术类型
        if 'spatial' in adata.obsm:
            n_spots = adata.n_obs
            if n_spots > 20000:
                info['likely_technology'] = 'slide-seq'
            elif n_spots > 2000:
                info['likely_technology'] = 'visium_or_st'  
            elif adata.n_vars < 100:
                info['likely_technology'] = 'targeted_spatial'
            else:
                info['likely_technology'] = 'unknown_spatial'
        else:
            info['likely_technology'] = 'non_spatial'
        
        return info
        
    except Exception as e:
        return {
            'filename': file_path.name,
            'error': str(e),
            'file_size_mb': float(file_path.stat().st_size / (1024 * 1024))
        }

def move_and_organize_datasets(datasets, base_dir):
    """移动和组织所有数据集"""
    moved_datasets = {}
    errors = []
    
    for dataset in datasets:
        try:
            # 确定目标类别
            category = categorize_dataset(dataset)
            target_dir = base_dir / category
            target_path = target_dir / dataset['filename']
            
            # 检查是否已经在目标位置
            if dataset['current_path'].resolve() == target_path.resolve():
                print(f"跳过 {dataset['filename']} (已在正确位置)")
                continue
            
            # 检查目标文件是否已存在
            if target_path.exists():
                # 比较文件大小，决定是否覆盖
                current_size = dataset['current_path'].stat().st_size
                existing_size = target_path.stat().st_size
                
                if current_size == existing_size:
                    print(f"跳过 {dataset['filename']} (目标位置已有相同文件)")
                    continue
                else:
                    # 创建备份名称
                    backup_name = f"{target_path.stem}_backup_{int(existing_size/(1024*1024))}MB{target_path.suffix}"
                    backup_path = target_dir / backup_name
                    shutil.move(target_path, backup_path)
                    print(f"备份现有文件: {backup_name}")
            
            # 移动文件
            shutil.move(dataset['current_path'], target_path)
            print(f"移动 {dataset['filename']} -> {category}/")
            
            # 分析数据集
            dataset_info = analyze_dataset(target_path)
            
            if category not in moved_datasets:
                moved_datasets[category] = []
            moved_datasets[category].append(dataset_info)
            
        except Exception as e:
            error_msg = f"处理 {dataset['filename']} 时出错: {str(e)}"
            print(f"❌ {error_msg}")
            errors.append(error_msg)
    
    return moved_datasets, errors

def update_readme_files(base_dir, moved_datasets):
    """更新各目录的README文件"""
    
    for category, datasets in moved_datasets.items():
        readme_path = base_dir / category / 'README.md'
        
        # 统计信息
        total_datasets = len(datasets)
        total_cells = sum(d.get('n_cells', 0) for d in datasets if 'n_cells' in d)
        total_genes = sum(d.get('n_genes', 0) for d in datasets if 'n_genes' in d) 
        total_size_mb = sum(d.get('file_size_mb', 0) for d in datasets)
        spatial_datasets = sum(1 for d in datasets if d.get('has_spatial', False))
        
        with open(readme_path, 'w') as f:
            f.write(f"# {category.replace('_', ' ').title()}\n\n")
            
            f.write(f"## 统计信息\n\n")
            f.write(f"- **数据集数量**: {total_datasets}\n")
            f.write(f"- **总细胞数**: {total_cells:,}\n")
            f.write(f"- **总基因数**: {total_genes:,}\n") 
            f.write(f"- **总文件大小**: {total_size_mb:.1f} MB\n")
            f.write(f"- **空间数据集**: {spatial_datasets}/{total_datasets}\n\n")
            
            f.write(f"## 数据集列表\n\n")
            f.write("| 数据集 | 细胞数 | 基因数 | 大小(MB) | 空间坐标 | 技术类型 | 聚类信息 |\n")
            f.write("|--------|--------|--------|----------|----------|----------|----------|\n")
            
            for dataset in sorted(datasets, key=lambda x: x.get('n_cells', 0), reverse=True):
                if 'error' in dataset:
                    f.write(f"| {dataset['filename']} | - | - | {dataset['file_size_mb']:.1f} | ❌ 加载失败 | - | - |\n")
                else:
                    spatial_icon = "✅" if dataset.get('has_spatial', False) else "❌"
                    cluster_info = "✅" if dataset.get('has_leiden', False) else f"{len(dataset.get('cluster_keys', []))} keys" if dataset.get('has_clusters', False) else "❌"
                    
                    f.write(f"| {dataset['filename']} | {dataset.get('n_cells', 0):,} | {dataset.get('n_genes', 0):,} | {dataset.get('file_size_mb', 0):.1f} | {spatial_icon} | {dataset.get('likely_technology', 'unknown')} | {cluster_info} |\n")
            
            f.write(f"\n更新时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

def create_master_catalog(base_dir, moved_datasets):
    """创建主数据集目录"""
    catalog = {
        'metadata': {
            'created_at': pd.Timestamp.now().isoformat(),
            'total_categories': len(moved_datasets),
            'total_datasets': sum(len(datasets) for datasets in moved_datasets.values()),
            'total_size_mb': sum(d.get('file_size_mb', 0) for datasets in moved_datasets.values() for d in datasets)
        },
        'categories': moved_datasets,
        'quick_access': {
            'spatial_datasets': [],
            'large_datasets': [],
            'demo_datasets': [],
            'benchmark_datasets': []
        }
    }
    
    # 填充快速访问索引
    for category, datasets in moved_datasets.items():
        for dataset in datasets:
            if 'error' not in dataset:
                if dataset.get('has_spatial', False):
                    catalog['quick_access']['spatial_datasets'].append({
                        'filename': dataset['filename'],
                        'category': category,
                        'n_cells': dataset.get('n_cells', 0),
                        'technology': dataset.get('likely_technology', 'unknown')
                    })
                
                if dataset.get('n_cells', 0) > 10000:
                    catalog['quick_access']['large_datasets'].append({
                        'filename': dataset['filename'],
                        'category': category,
                        'n_cells': dataset.get('n_cells', 0)
                    })
    
    # 按类别添加快速访问
    if 'demo_datasets' in moved_datasets:
        catalog['quick_access']['demo_datasets'] = [
            {'filename': d['filename'], 'n_cells': d.get('n_cells', 0)} 
            for d in moved_datasets['demo_datasets'] if 'error' not in d
        ]
    
    if 'benchmark_datasets' in moved_datasets:
        catalog['quick_access']['benchmark_datasets'] = [
            {'filename': d['filename'], 'n_cells': d.get('n_cells', 0)}
            for d in moved_datasets['benchmark_datasets'] if 'error' not in d
        ]
    
    # 保存目录文件
    catalog_path = base_dir / 'metadata' / 'datasets_catalog.json'
    with open(catalog_path, 'w') as f:
        json.dump(catalog, f, indent=2)
    
    print(f"主数据集目录已保存: {catalog_path}")
    
    # 创建简化的索引文件
    index_path = base_dir / 'DATASETS_INDEX.md'
    with open(index_path, 'w') as f:
        f.write("# ChatSpatial 数据集索引\n\n")
        f.write(f"**更新时间**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"**总数据集**: {catalog['metadata']['total_datasets']}\n")
        f.write(f"**总大小**: {catalog['metadata']['total_size_mb']:.1f} MB\n\n")
        
        f.write("## 快速访问\n\n")
        
        f.write("### 🔬 空间数据集 (推荐用于空间分析)\n")
        for dataset in sorted(catalog['quick_access']['spatial_datasets'], key=lambda x: x['n_cells'], reverse=True)[:10]:
            f.write(f"- **{dataset['filename']}**: {dataset['n_cells']:,} cells ({dataset['technology']})\n")
        
        f.write("\n### 📊 大规模数据集 (推荐用于性能测试)\n")
        for dataset in sorted(catalog['quick_access']['large_datasets'], key=lambda x: x['n_cells'], reverse=True)[:5]:
            f.write(f"- **{dataset['filename']}**: {dataset['n_cells']:,} cells\n")
        
        f.write("\n### 🎯 演示数据集 (推荐用于教程)\n")
        for dataset in sorted(catalog['quick_access']['demo_datasets'], key=lambda x: x['n_cells'])[:5]:
            f.write(f"- **{dataset['filename']}**: {dataset['n_cells']:,} cells\n")
        
        f.write(f"\n## 详细分类\n\n")
        for category in moved_datasets.keys():
            f.write(f"- **[{category.replace('_', ' ').title()}]({category}/README.md)**: {len(moved_datasets[category])} datasets\n")
    
    print(f"数据集索引已创建: {index_path}")

def main():
    """主函数"""
    print("开始整理ChatSpatial数据集...")
    print("="*50)
    
    # 1. 创建目录结构
    print("1. 创建标准目录结构...")
    base_dir, directories = create_data_structure()
    
    # 2. 收集所有数据集
    print("\n2. 收集现有数据集...")
    datasets = collect_all_datasets()
    
    if not datasets:
        print("没有找到数据集文件")
        return
    
    # 3. 移动和组织数据集
    print(f"\n3. 移动和组织 {len(datasets)} 个数据集...")
    moved_datasets, errors = move_and_organize_datasets(datasets, base_dir)
    
    # 4. 更新README文件
    print("\n4. 更新文档...")
    update_readme_files(base_dir, moved_datasets)
    
    # 5. 创建主目录
    print("\n5. 创建主数据集目录...")
    create_master_catalog(base_dir, moved_datasets)
    
    # 6. 清理空目录
    print("\n6. 清理测试目录...")
    test_datasets_dir = Path("/Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/datasets")
    if test_datasets_dir.exists():
        # 检查是否还有h5ad文件
        remaining_files = list(test_datasets_dir.rglob("*.h5ad"))
        if not remaining_files:
            print("测试目录中无剩余数据集，保留目录结构")
        else:
            print(f"测试目录中还有 {len(remaining_files)} 个文件未处理")
    
    # 7. 最终报告
    print("\n" + "="*50)
    print("数据集整理完成!")
    print("="*50)
    
    total_moved = sum(len(datasets) for datasets in moved_datasets.values())
    print(f"✅ 成功整理 {total_moved} 个数据集")
    
    for category, datasets in moved_datasets.items():
        valid_datasets = [d for d in datasets if 'error' not in d]
        error_datasets = [d for d in datasets if 'error' in d]
        print(f"   📁 {category}: {len(valid_datasets)} 有效, {len(error_datasets)} 错误")
    
    if errors:
        print(f"\n⚠️  处理过程中出现 {len(errors)} 个错误:")
        for error in errors:
            print(f"   - {error}")
    
    print(f"\n📍 数据集位置: {base_dir}")
    print(f"📖 查看索引: {base_dir}/DATASETS_INDEX.md")
    print(f"🔍 详细目录: {base_dir}/metadata/datasets_catalog.json")

if __name__ == "__main__":
    main()