#!/usr/bin/env python3
"""
下载和准备全面测试所需的数据集
"""

import os
import requests
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import tarfile
from urllib.request import urlretrieve
import warnings
warnings.filterwarnings('ignore')

# 设置数据目录
BASE_DIR = Path(__file__).parent / 'datasets'
BASE_DIR.mkdir(exist_ok=True)

def create_synthetic_datasets():
    """创建合成测试数据集"""
    print("创建合成数据集...")
    
    # 小型快速测试数据 (100 cells, 500 genes)
    small_adata = sc.AnnData(
        X=np.random.negative_binomial(n=10, p=0.3, size=(100, 500)),
        obs=pd.DataFrame({
            'cell_type': np.random.choice(['TypeA', 'TypeB', 'TypeC'], 100),
            'batch': np.random.choice(['batch1', 'batch2'], 100),
        }),
        var=pd.DataFrame({
            'gene_symbol': [f'Gene_{i}' for i in range(500)],
            'highly_variable': np.random.choice([True, False], 500, p=[0.2, 0.8])
        })
    )
    
    # 添加空间坐标
    small_adata.obsm['spatial'] = np.random.rand(100, 2) * 1000
    small_adata.write(BASE_DIR / 'small_synthetic.h5ad')
    
    # 中型标准测试数据 (1000 cells, 2000 genes) 
    medium_adata = sc.AnnData(
        X=np.random.negative_binomial(n=10, p=0.3, size=(1000, 2000)),
        obs=pd.DataFrame({
            'cell_type': np.random.choice(['TypeA', 'TypeB', 'TypeC', 'TypeD'], 1000),
            'batch': np.random.choice(['batch1', 'batch2', 'batch3'], 1000),
            'condition': np.random.choice(['control', 'treatment'], 1000),
        }),
        var=pd.DataFrame({
            'gene_symbol': [f'Gene_{i}' for i in range(2000)],
            'highly_variable': np.random.choice([True, False], 2000, p=[0.3, 0.7])
        })
    )
    medium_adata.obsm['spatial'] = np.random.rand(1000, 2) * 2000
    medium_adata.write(BASE_DIR / 'medium_synthetic.h5ad')
    
    # 大型压力测试数据 (5000 cells, 3000 genes)
    large_adata = sc.AnnData(
        X=np.random.negative_binomial(n=15, p=0.2, size=(5000, 3000)),
        obs=pd.DataFrame({
            'cell_type': np.random.choice(['TypeA', 'TypeB', 'TypeC', 'TypeD', 'TypeE'], 5000),
            'batch': np.random.choice(['batch1', 'batch2', 'batch3', 'batch4'], 5000),
        }),
        var=pd.DataFrame({
            'gene_symbol': [f'Gene_{i}' for i in range(3000)],
            'highly_variable': np.random.choice([True, False], 3000, p=[0.25, 0.75])
        })
    )
    large_adata.obsm['spatial'] = np.random.rand(5000, 2) * 3000
    large_adata.write(BASE_DIR / 'large_synthetic.h5ad')
    
    print(f"合成数据集创建完成，保存在: {BASE_DIR}")

def download_10x_genomics_datasets():
    """下载10X Genomics公开数据集"""
    print("下载10X Genomics数据集...")
    
    datasets = {
        'visium_human_brain': {
            'url': 'https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_2/V1_Human_Brain_Section_2_filtered_feature_bc_matrix.h5',
            'spatial_url': 'https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_2/V1_Human_Brain_Section_2_spatial.tar.gz'
        }
    }
    
    for name, urls in datasets.items():
        try:
            print(f"  下载 {name}...")
            data_file = BASE_DIR / f"{name}_matrix.h5"
            spatial_file = BASE_DIR / f"{name}_spatial.tar.gz"
            
            if not data_file.exists():
                urlretrieve(urls['url'], data_file)
                print(f"    数据文件下载完成: {data_file}")
            
            if not spatial_file.exists():
                urlretrieve(urls['spatial_url'], spatial_file)
                print(f"    空间文件下载完成: {spatial_file}")
                
                # 解压空间数据
                with tarfile.open(spatial_file, 'r:gz') as tar:
                    tar.extractall(BASE_DIR / name)
                    
        except Exception as e:
            print(f"    下载 {name} 失败: {e}")

def download_squidpy_datasets():
    """下载squidpy提供的标准数据集"""
    print("下载squidpy标准数据集...")
    
    try:
        import squidpy as sq
        
        # 下载Visium数据
        adata = sq.datasets.visium_hne_adata()
        adata.write(BASE_DIR / 'squidpy_visium.h5ad')
        print("  Visium数据集下载完成")
        
        # 下载seqFISH数据  
        adata = sq.datasets.seqfish()
        adata.write(BASE_DIR / 'squidpy_seqfish.h5ad')
        print("  seqFISH数据集下载完成")
        
    except ImportError:
        print("  squidpy未安装，跳过")
    except Exception as e:
        print(f"  squidpy数据集下载失败: {e}")

def create_edge_case_datasets():
    """创建边界情况测试数据"""
    print("创建边界情况数据集...")
    
    # 空数据集
    empty_adata = sc.AnnData(X=np.empty((0, 100)))
    empty_adata.write(BASE_DIR / 'empty_dataset.h5ad')
    
    # 单细胞数据
    single_cell_adata = sc.AnnData(
        X=np.random.negative_binomial(n=10, p=0.3, size=(1, 1000)),
        obs=pd.DataFrame({'cell_type': ['single']})
    )
    single_cell_adata.write(BASE_DIR / 'single_cell.h5ad')
    
    # 缺失空间坐标数据
    no_spatial_adata = sc.AnnData(
        X=np.random.negative_binomial(n=10, p=0.3, size=(100, 500)),
        obs=pd.DataFrame({'cell_type': np.random.choice(['A', 'B'], 100)})
    )
    no_spatial_adata.write(BASE_DIR / 'no_spatial.h5ad')
    
    # 极高稀疏数据 (>95% zeros)
    sparse_X = np.random.negative_binomial(n=1, p=0.95, size=(500, 1000))
    sparse_adata = sc.AnnData(X=sparse_X)
    sparse_adata.obsm['spatial'] = np.random.rand(500, 2) * 1000
    sparse_adata.write(BASE_DIR / 'high_sparsity.h5ad')
    
    print("边界情况数据集创建完成")

def create_benchmark_datasets():
    """创建性能基准测试数据"""
    print("创建基准测试数据集...")
    
    sizes = [
        (500, 1000, 'benchmark_500x1k'),
        (1000, 2000, 'benchmark_1kx2k'), 
        (2000, 3000, 'benchmark_2kx3k'),
        (5000, 5000, 'benchmark_5kx5k'),
    ]
    
    for n_cells, n_genes, name in sizes:
        print(f"  创建 {name} ({n_cells} cells x {n_genes} genes)")
        adata = sc.AnnData(
            X=np.random.negative_binomial(n=10, p=0.3, size=(n_cells, n_genes)),
            obs=pd.DataFrame({
                'cell_type': np.random.choice(['A', 'B', 'C', 'D'], n_cells),
                'batch': np.random.choice(['b1', 'b2', 'b3'], n_cells),
            }),
            var=pd.DataFrame({
                'gene_symbol': [f'Gene_{i}' for i in range(n_genes)]
            })
        )
        adata.obsm['spatial'] = np.random.rand(n_cells, 2) * 1000
        adata.write(BASE_DIR / f'{name}.h5ad')

def generate_dataset_summary():
    """生成数据集摘要"""
    print("生成数据集摘要...")
    
    summary = []
    for h5ad_file in BASE_DIR.glob('*.h5ad'):
        try:
            adata = sc.read_h5ad(h5ad_file)
            summary.append({
                'dataset': h5ad_file.name,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'has_spatial': 'spatial' in adata.obsm,
                'file_size_mb': h5ad_file.stat().st_size / 1024 / 1024,
                'sparsity': (adata.X == 0).sum() / adata.X.size if hasattr(adata.X, 'size') else 'N/A'
            })
        except Exception as e:
            summary.append({
                'dataset': h5ad_file.name,
                'error': str(e)
            })
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(BASE_DIR / 'datasets_summary.csv', index=False)
    print(f"数据集摘要保存至: {BASE_DIR / 'datasets_summary.csv'}")
    print(summary_df)

if __name__ == "__main__":
    print("开始准备测试数据集...")
    
    # 创建各类数据集
    create_synthetic_datasets()
    create_edge_case_datasets() 
    create_benchmark_datasets()
    
    # 下载公开数据集
    download_10x_genomics_datasets()
    download_squidpy_datasets()
    
    # 生成摘要
    generate_dataset_summary()
    
    print(f"\n所有测试数据集准备完成！")
    print(f"数据集保存目录: {BASE_DIR}")
    print(f"共计 {len(list(BASE_DIR.glob('*.h5ad')))} 个h5ad文件")