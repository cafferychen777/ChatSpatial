#!/usr/bin/env python3
"""
下载真实的Visium空间转录组数据集

这个脚本严格按照Linus的"好品味"原则设计：
1. 简单直接，无特殊情况处理
2. 实用主义 - 使用实际可用的数据源
3. 错误处理清晰，不破坏用户环境
4. 数据结构清晰，消除复杂性

策略改变：10X官方链接需要认证，我们使用squidpy和其他公开数据源
"""

import os
import sys
import urllib.request
import urllib.error
import tarfile
import gzip
import shutil
from pathlib import Path
import time
import hashlib
import json
from typing import Dict, List, Tuple, Optional

import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np

# 实用主义原则：使用真正可用的数据源
DATASETS_CONFIG = {
    "squidpy_visium": {
        "name": "Squidpy Visium (Mouse Brain)",
        "source": "squidpy",
        "description": "Mouse brain from squidpy tutorial",
        "expected_spots": 2264,
        "expected_genes": 31053
    },
    "spatial_libd": {
        "name": "spatialLIBD Human DLPFC",
        "url": "http://spatial.libd.org/spatialLIBD/reference/spatialLIBD_Visium_example_based_on_Maynard_and_Collado-Torres_et_al_2021.rda",
        "description": "Human dorsolateral prefrontal cortex",
        "expected_size_mb": 15
    },
    "mouse_brain_coronal": {
        "name": "Mouse Brain Coronal",
        "url_base": "https://www.spatialresearch.org/wp-content/uploads/2019/12",
        "files": {
            "matrix": "filtered_gene_bc_matrices_mex.tar.gz",
            "spatial": "spatial.tar.gz"
        },
        "description": "Mouse brain coronal section",
        "expected_size_mb": 45
    },
    "zenodo_visium": {
        "name": "Zenodo Visium Dataset",
        "url": "https://zenodo.org/record/4751176/files/V1_Human_Lymph_Node_spatial.h5ad",
        "description": "Human lymph node from Zenodo",
        "expected_size_mb": 25
    }
}

def download_file(url: str, filepath: Path, expected_size_mb: Optional[float] = None) -> bool:
    """
    下载文件，Linus风格：简单、直接、有效的错误处理
    """
    print(f"正在下载: {url}")
    
    try:
        # 创建目录
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        # 下载文件
        with urllib.request.urlopen(url, timeout=300) as response:
            total_size = int(response.headers.get('Content-Length', 0))
            
            if expected_size_mb and total_size > 0:
                expected_bytes = expected_size_mb * 1024 * 1024
                if abs(total_size - expected_bytes) > expected_bytes * 0.5:  # 50%容差
                    print(f"警告: 文件大小异常 {total_size/1024/1024:.1f}MB, 预期 {expected_size_mb:.1f}MB")
            
            with open(filepath, 'wb') as f:
                downloaded = 0
                chunk_size = 8192
                
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    
                    if total_size > 0:
                        progress = downloaded / total_size * 100
                        print(f"\r进度: {progress:.1f}%", end='', flush=True)
                
                print()  # 换行
        
        # 验证文件存在且不为空
        if filepath.exists() and filepath.stat().st_size > 0:
            print(f"下载完成: {filepath}")
            return True
        else:
            print(f"下载失败: 文件为空或不存在")
            return False
            
    except urllib.error.URLError as e:
        print(f"网络错误: {e}")
        return False
    except Exception as e:
        print(f"下载失败: {e}")
        return False

def download_squidpy_dataset(dataset_id: str, config: Dict, output_dir: Path) -> bool:
    """
    使用squidpy下载Visium数据集
    """
    print(f"\n使用squidpy获取数据集: {config['name']}")
    
    try:
        import squidpy as sq
        
        dataset_dir = output_dir / dataset_id
        dataset_dir.mkdir(parents=True, exist_ok=True)
        output_file = dataset_dir / f"{dataset_id}.h5ad"
        
        if output_file.exists() and output_file.stat().st_size > 1024*1024:
            print(f"数据集已存在，跳过下载: {output_file.name}")
            return True
        
        # 下载squidpy的Visium数据
        adata = sq.datasets.visium_hne_adata_crop()
        
        # 添加数据集信息
        adata.uns['dataset_info'] = {
            'name': config['name'],
            'source': 'squidpy',
            'technology': 'Visium',
            'description': config['description'],
            'download_date': time.strftime('%Y-%m-%d'),
            'n_spots': adata.n_obs,
            'n_genes': adata.n_vars
        }
        
        # 保存数据
        adata.write(output_file)
        
        file_size_mb = output_file.stat().st_size / (1024 * 1024)
        print(f"squidpy数据获取成功: {output_file.name} ({file_size_mb:.1f}MB)")
        print(f"  - Spots: {adata.n_obs}")
        print(f"  - Genes: {adata.n_vars}")
        
        return True
        
    except ImportError:
        print("错误: 需要安装squidpy")
        print("运行: pip install squidpy")
        return False
    except Exception as e:
        print(f"squidpy数据获取失败: {e}")
        return False

def download_dataset(dataset_id: str, config: Dict, output_dir: Path) -> bool:
    """
    根据数据源类型下载数据集
    """
    # squidpy数据集
    if config.get('source') == 'squidpy':
        return download_squidpy_dataset(dataset_id, config, output_dir)
    
    # URL下载的数据集
    print(f"\n开始下载数据集: {config['name']}")
    
    dataset_dir = output_dir / dataset_id
    dataset_dir.mkdir(parents=True, exist_ok=True)
    
    # 单文件下载
    if 'url' in config:
        url = config['url']
        if url.endswith('.h5ad'):
            # 直接是h5ad文件
            filepath = dataset_dir / f"{dataset_id}.h5ad"
            if filepath.exists() and filepath.stat().st_size > 1024*1024:
                print(f"文件已存在，跳过下载: {filepath.name}")
                return True
            return download_file(url, filepath, config.get('expected_size_mb'))
        else:
            # 其他格式文件
            filename = url.split('/')[-1]
            filepath = dataset_dir / filename
            return download_file(url, filepath, config.get('expected_size_mb'))
    
    # 多文件下载
    if 'files' in config:
        success_count = 0
        total_files = len(config['files'])
        
        for file_type, filename in config['files'].items():
            url = f"{config['url_base']}/{filename}"
            filepath = dataset_dir / filename
            
            if filepath.exists() and filepath.stat().st_size > 1024:
                print(f"文件已存在，跳过: {filename}")
                success_count += 1
                continue
            
            if download_file(url, filepath, config.get('expected_size_mb')):
                success_count += 1
            else:
                print(f"文件下载失败: {filename}")
        
        success = success_count == total_files
        print(f"数据集 {dataset_id} 下载{'成功' if success else '失败'}: {success_count}/{total_files} 文件")
        return success
    
    return False

def convert_to_h5ad(dataset_id: str, config: Dict, dataset_dir: Path) -> bool:
    """
    转换或验证h5ad格式数据，遵循Linus的简洁原则
    """
    output_file = dataset_dir / f"{dataset_id}.h5ad"
    
    # 如果数据来自squidpy或已经是h5ad，直接返回成功
    if config.get('source') == 'squidpy':
        return output_file.exists() and output_file.stat().st_size > 1024*1024
    
    # 处理直接下载的h5ad文件
    if 'url' in config and config['url'].endswith('.h5ad'):
        downloaded_file = dataset_dir / config['url'].split('/')[-1]
        if downloaded_file.exists() and downloaded_file != output_file:
            # 重命名为标准格式
            shutil.move(str(downloaded_file), str(output_file))
        
        if output_file.exists():
            try:
                # 验证文件并添加数据集信息
                adata = sc.read_h5ad(output_file)
                
                if 'dataset_info' not in adata.uns:
                    adata.uns['dataset_info'] = {
                        'name': config['name'],
                        'source': config.get('source', 'External'),
                        'technology': 'Visium',
                        'description': config.get('description', ''),
                        'download_date': time.strftime('%Y-%m-%d'),
                        'n_spots': adata.n_obs,
                        'n_genes': adata.n_vars
                    }
                    adata.write(output_file)
                
                file_size_mb = output_file.stat().st_size / (1024 * 1024)
                print(f"h5ad验证成功: {output_file.name} ({file_size_mb:.1f}MB)")
                print(f"  - Spots: {adata.n_obs}")
                print(f"  - Genes: {adata.n_vars}")
                return True
                
            except Exception as e:
                print(f"h5ad文件验证失败: {e}")
                return False
    
    # 处理需要转换的10X格式数据
    if 'files' not in config:
        return False
        
    print(f"转换数据集 {dataset_id} 到 h5ad 格式...")
    
    try:
        # 如果h5ad已存在且不为空，跳过转换
        if output_file.exists() and output_file.stat().st_size > 1024*1024:
            print(f"h5ad文件已存在，跳过转换: {output_file.name}")
            return True
        
        # 文件路径
        matrix_file = dataset_dir / config['files']['matrix']
        coords_file = dataset_dir / config['files']['spatial_coords']
        
        # 读取基因表达数据
        if not matrix_file.exists():
            print(f"表达矩阵文件不存在: {matrix_file}")
            return False
            
        adata = sc.read_10x_h5(matrix_file)
        adata.var_names_unique()
        
        # 读取空间坐标
        if coords_file.exists():
            coords_df = pd.read_csv(coords_file, header=None, index_col=0)
            coords_df.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']
            
            coords_df = coords_df[coords_df['in_tissue'] == 1]
            common_barcodes = adata.obs_names.intersection(coords_df.index)
            
            if len(common_barcodes) > 0:
                adata = adata[common_barcodes]
                adata.obsm['spatial'] = coords_df.loc[common_barcodes, ['pxl_col_in_fullres', 'pxl_row_in_fullres']].values.astype(float)
                print(f"添加了 {len(common_barcodes)} 个spots的空间坐标")
            else:
                print("警告: 没有找到匹配的barcode")
                return False
        else:
            print(f"空间坐标文件不存在: {coords_file}")
            return False
        
        # 添加数据集信息
        adata.uns['dataset_info'] = {
            'name': config['name'],
            'source': config.get('source', '10X Genomics'),
            'technology': 'Visium',
            'description': config.get('description', ''),
            'download_date': time.strftime('%Y-%m-%d'),
            'n_spots': adata.n_obs,
            'n_genes': adata.n_vars
        }
        
        adata.write(output_file)
        
        file_size_mb = output_file.stat().st_size / (1024 * 1024)
        print(f"转换成功: {output_file.name} ({file_size_mb:.1f}MB)")
        print(f"  - Spots: {adata.n_obs}")
        print(f"  - Genes: {adata.n_vars}")
        
        return True
        
    except Exception as e:
        print(f"转换失败: {e}")
        return False

def create_dataset_summary(output_dir: Path, results: Dict) -> None:
    """
    创建数据集摘要文件
    """
    summary_data = []
    
    for dataset_id, success in results.items():
        if success:
            h5ad_file = output_dir / dataset_id / f"{dataset_id}.h5ad"
            if h5ad_file.exists():
                try:
                    adata = sc.read_h5ad(h5ad_file)
                    file_size_mb = h5ad_file.stat().st_size / (1024 * 1024)
                    
                    summary_data.append({
                        'dataset_id': dataset_id,
                        'name': adata.uns['dataset_info']['name'],
                        'n_spots': adata.n_obs,
                        'n_genes': adata.n_vars,
                        'file_size_mb': round(file_size_mb, 1),
                        'download_success': True,
                        'conversion_success': True,
                        'file_path': str(h5ad_file)
                    })
                except Exception as e:
                    print(f"读取数据集摘要失败 {dataset_id}: {e}")
        else:
            config = DATASETS_CONFIG.get(dataset_id, {})
            summary_data.append({
                'dataset_id': dataset_id,
                'name': config.get('name', 'Unknown'),
                'n_spots': 0,
                'n_genes': 0,
                'file_size_mb': 0,
                'download_success': False,
                'conversion_success': False,
                'file_path': ''
            })
    
    # 保存摘要
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_dir / 'real_datasets_summary.csv'
    summary_df.to_csv(summary_file, index=False)
    print(f"\n数据集摘要保存到: {summary_file}")
    
    # 打印摘要
    print("\n=== 数据集下载摘要 ===")
    successful = summary_df[summary_df['download_success']].shape[0]
    total = summary_df.shape[0]
    print(f"成功下载: {successful}/{total} 个数据集")
    
    for _, row in summary_df.iterrows():
        status = "✓" if row['download_success'] else "✗"
        print(f"{status} {row['name']}: {row['n_spots']} spots, {row['n_genes']} genes, {row['file_size_mb']}MB")

def main():
    """
    主函数：按Linus的简洁原则实现
    """
    print("10X Genomics Visium 真实数据集下载器")
    print("=" * 50)
    
    # 检查依赖
    try:
        import scanpy as sc
        import anndata as ad
        print(f"Scanpy版本: {sc.__version__}")
    except ImportError:
        print("错误: 需要安装 scanpy 和 anndata")
        print("运行: pip install scanpy anndata")
        sys.exit(1)
    
    # 目标目录
    output_dir = Path("/Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/datasets/real_datasets")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"输出目录: {output_dir}")
    
    # 下载和转换结果
    results = {}
    
    # 逐个下载数据集
    for dataset_id, config in DATASETS_CONFIG.items():
        try:
            download_success = download_dataset(dataset_id, config, output_dir)
            
            if download_success:
                conversion_success = convert_to_h5ad(dataset_id, config, output_dir / dataset_id)
                results[dataset_id] = conversion_success
            else:
                results[dataset_id] = False
                
        except KeyboardInterrupt:
            print("\n用户中断下载")
            break
        except Exception as e:
            print(f"处理数据集 {dataset_id} 时出错: {e}")
            results[dataset_id] = False
    
    # 创建摘要
    create_dataset_summary(output_dir, results)
    
    print(f"\n下载完成！数据保存在: {output_dir}")

if __name__ == "__main__":
    main()