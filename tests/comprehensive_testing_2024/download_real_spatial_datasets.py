#!/usr/bin/env python3
"""
下载真实的空间转录组数据集用于comprehensive testing
包括seqFISH、MERFISH等技术的高质量已发表数据

数据来源：
- Squidpy内置数据集
- 10X Genomics公开数据
- 已发表论文的补充数据
- GEO/SRA数据库
"""

import os
import sys
import requests
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import gzip
import tarfile
from urllib.request import urlretrieve, urlopen
from urllib.error import URLError
import warnings
import json
from typing import Dict, List, Optional, Tuple
import time

warnings.filterwarnings('ignore')

# 设置数据目录
BASE_DIR = Path(__file__).parent / 'datasets' / 'real_datasets'
BASE_DIR.mkdir(parents=True, exist_ok=True)

# 数据集元数据记录
DATASETS_METADATA = {}

def log_dataset_info(dataset_name: str, **metadata):
    """记录数据集信息"""
    DATASETS_METADATA[dataset_name] = {
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        **metadata
    }

def save_metadata():
    """保存所有数据集的元数据"""
    metadata_file = BASE_DIR / 'datasets_metadata.json'
    with open(metadata_file, 'w') as f:
        json.dump(DATASETS_METADATA, f, indent=2, ensure_ascii=False)
    print(f"数据集元数据保存至: {metadata_file}")

def download_squidpy_datasets():
    """下载squidpy内置的标准数据集"""
    print("\n=== 下载Squidpy内置数据集 ===")
    
    try:
        import squidpy as sq
        print(f"✓ Squidpy版本: {sq.__version__}")
        
        datasets_to_download = [
            ('seqfish', 'seqFISH+ Mouse Embryo', 'Lohoff et al., Nature 2022'),
            ('merfish', 'MERFISH Hypothalamus', 'Moffitt et al., Science 2018'),
            ('mibitof', 'MIBIToF Breast Cancer', 'Keren et al., Cell 2018'),
            ('slideseqv2', 'Slide-seqV2 Mouse Hippocampus', 'Stickels et al., Nature Biotechnology 2021'),
            ('visium', '10X Visium Brain Section', '10X Genomics Dataset'),
            ('imc', 'IMC Breast Cancer', 'Jackson et al., Nature 2020'),
        ]
        
        for dataset_id, description, citation in datasets_to_download:
            try:
                print(f"  下载 {dataset_id} ({description})...")
                
                # 动态获取数据集
                dataset_func = getattr(sq.datasets, dataset_id, None)
                if dataset_func is None:
                    print(f"    ❌ {dataset_id} 不可用")
                    continue
                
                adata = dataset_func()
                output_file = BASE_DIR / f'squidpy_{dataset_id}.h5ad'
                adata.write(output_file)
                
                # 记录元数据
                log_dataset_info(f'squidpy_{dataset_id}', 
                    source='squidpy',
                    technology=dataset_id.upper(),
                    description=description,
                    citation=citation,
                    n_cells=adata.n_obs,
                    n_genes=adata.n_vars,
                    has_spatial='spatial' in adata.obsm,
                    file_path=str(output_file.relative_to(BASE_DIR.parent.parent)),
                    file_size_mb=round(output_file.stat().st_size / 1024 / 1024, 2)
                )
                
                print(f"    ✓ 完成: {adata.n_obs} cells, {adata.n_vars} genes")
                
            except Exception as e:
                print(f"    ❌ 下载 {dataset_id} 失败: {e}")
        
    except ImportError:
        print("❌ squidpy未安装，跳过内置数据集下载")
    except Exception as e:
        print(f"❌ 下载squidpy数据集时出错: {e}")

def download_10x_genomics_datasets():
    """下载10X Genomics Visium公开数据集"""
    print("\n=== 下载10X Genomics公开数据集 ===")
    
    # 10X Genomics 公开数据集
    datasets = {
        'visium_human_brain': {
            'name': 'Human Brain Section (Anterior)',
            'matrix_url': 'https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_2/V1_Human_Brain_Section_2_filtered_feature_bc_matrix.h5',
            'spatial_url': 'https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Brain_Section_2/V1_Human_Brain_Section_2_spatial.tar.gz',
            'citation': '10X Genomics, Spatial Gene Expression Dataset',
        },
        'visium_mouse_brain': {
            'name': 'Mouse Brain Section (Sagittal-Anterior)',
            'matrix_url': 'https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5',
            'spatial_url': 'https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_spatial.tar.gz',
            'citation': '10X Genomics, Spatial Gene Expression Dataset',
        }
    }
    
    for dataset_id, info in datasets.items():
        try:
            print(f"  下载 {dataset_id} ({info['name']})...")
            dataset_dir = BASE_DIR / dataset_id
            dataset_dir.mkdir(exist_ok=True)
            
            # 下载表达数据
            matrix_file = dataset_dir / 'filtered_feature_bc_matrix.h5'
            if not matrix_file.exists():
                print(f"    下载表达数据...")
                urlretrieve(info['matrix_url'], matrix_file)
                print(f"    ✓ 表达数据下载完成: {matrix_file.name}")
            
            # 下载空间数据
            spatial_archive = dataset_dir / 'spatial.tar.gz'
            if not spatial_archive.exists():
                print(f"    下载空间数据...")
                urlretrieve(info['spatial_url'], spatial_archive)
                
                # 解压空间数据
                with tarfile.open(spatial_archive, 'r:gz') as tar:
                    tar.extractall(dataset_dir)
                
                print(f"    ✓ 空间数据下载并解压完成")
                
                # 删除压缩文件以节省空间
                spatial_archive.unlink()
            
            # 使用scanpy读取并保存为h5ad
            try:
                import scanpy as sc
                h5ad_file = BASE_DIR / f'{dataset_id}.h5ad'
                if not h5ad_file.exists():
                    adata = sc.read_10x_h5(matrix_file)
                    adata.var_names_unique()
                    
                    # 读取空间坐标
                    spatial_dir = dataset_dir / 'spatial'
                    if spatial_dir.exists():
                        # 读取tissue_positions_list.csv
                        positions_file = spatial_dir / 'tissue_positions_list.csv'
                        if positions_file.exists():
                            positions = pd.read_csv(positions_file, header=None, index_col=0)
                            positions.columns = ['in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
                            
                            # 匹配细胞并添加空间坐标
                            common_barcodes = adata.obs.index.intersection(positions.index)
                            adata = adata[common_barcodes].copy()
                            adata.obsm['spatial'] = positions.loc[common_barcodes, ['pxl_col_in_fullres', 'pxl_row_in_fullres']].values
                    
                    adata.write(h5ad_file)
                    
                    # 记录元数据
                    log_dataset_info(dataset_id,
                        source='10x_genomics',
                        technology='Visium',
                        description=info['name'],
                        citation=info['citation'],
                        n_cells=adata.n_obs,
                        n_genes=adata.n_vars,
                        has_spatial='spatial' in adata.obsm,
                        file_path=str(h5ad_file.relative_to(BASE_DIR.parent.parent)),
                        file_size_mb=round(h5ad_file.stat().st_size / 1024 / 1024, 2)
                    )
                    
                    print(f"    ✓ 转换为h5ad: {adata.n_obs} spots, {adata.n_vars} genes")
            
            except Exception as e:
                print(f"    ⚠️  转换h5ad失败: {e}")
                # 仍然记录原始数据
                log_dataset_info(dataset_id,
                    source='10x_genomics',
                    technology='Visium',
                    description=info['name'],
                    citation=info['citation'],
                    raw_data_dir=str(dataset_dir.relative_to(BASE_DIR.parent.parent)),
                    status='raw_data_only'
                )
                    
        except Exception as e:
            print(f"    ❌ 下载 {dataset_id} 失败: {e}")

def download_published_datasets():
    """下载已发表论文的补充数据"""
    print("\n=== 下载已发表论文数据 ===")
    
    # 已知的高质量数据集链接
    published_datasets = {
        'seqfish_plus_embryo': {
            'name': 'seqFISH+ Mouse Embryo E11.5',
            'description': '原始seqFISH+小鼠胚胎数据，空间分辨率单细胞',
            'citation': 'Lohoff et al., Nature 2022. DOI: 10.1038/s41586-021-04353-z',
            'data_urls': [
                # 这些URL需要根据实际的数据可用性进行更新
                # 通常在论文的Data Availability部分或补充材料中找到
            ],
            'geo_accession': 'GSE166692',  # 如果有GEO数据
        },
        'merfish_hypothalamus': {
            'name': 'MERFISH Mouse Hypothalamus',
            'description': '原始MERFISH小鼠下丘脑数据，高通量空间转录组',
            'citation': 'Moffitt et al., Science 2018. DOI: 10.1126/science.aau5324',
            'data_urls': [],
            'note': '原始数据可能需要从作者直接获取',
        },
        'merfish_motor_cortex': {
            'name': 'MERFISH Human Motor Cortex',
            'description': 'MERFISH人类运动皮层数据',
            'citation': 'Zhang et al., Science 2023. DOI: 10.1126/science.adf6812',
            'data_urls': [],
            'note': 'BICCN数据集，需要从BICCN门户下载',
        }
    }
    
    for dataset_id, info in published_datasets.items():
        print(f"  {dataset_id}: {info['name']}")
        print(f"    描述: {info['description']}")
        print(f"    引用: {info['citation']}")
        
        if info.get('geo_accession'):
            print(f"    GEO登录号: {info['geo_accession']}")
            print(f"    ℹ️  请手动从GEO下载: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={info['geo_accession']}")
        
        if info.get('note'):
            print(f"    注意: {info['note']}")
        
        # 记录元数据（即使没有直接下载）
        log_dataset_info(dataset_id,
            source='published_paper',
            description=info['description'],
            citation=info['citation'],
            geo_accession=info.get('geo_accession', ''),
            note=info.get('note', ''),
            status='metadata_only'
        )
        
        print()

def download_demo_datasets():
    """下载一些轻量级的演示数据集用于快速测试"""
    print("\n=== 创建演示数据集 ===")
    
    try:
        import squidpy as sq
        
        # 使用squidpy创建一些小型的演示数据集
        demo_datasets = [
            ('seqfish_demo', 'sq.datasets.seqfish', 'seqFISH+ 演示数据'),
            ('visium_demo', 'sq.datasets.visium_hne_adata', 'Visium H&E 演示数据'),
        ]
        
        for demo_name, sq_function, description in demo_datasets:
            try:
                print(f"  创建 {demo_name}...")
                
                # 动态调用squidpy函数
                parts = sq_function.split('.')
                func = sq
                for part in parts[1:]:  # 跳过'sq'
                    func = getattr(func, part)
                
                adata = func()
                
                # 如果数据太大，进行下采样
                if adata.n_obs > 1000:
                    print(f"    下采样从 {adata.n_obs} 到 1000 个细胞...")
                    sc.pp.subsample(adata, n_obs=1000, random_state=42)
                
                output_file = BASE_DIR / f'{demo_name}.h5ad'
                adata.write(output_file)
                
                log_dataset_info(demo_name,
                    source='squidpy_demo',
                    description=description,
                    n_cells=adata.n_obs,
                    n_genes=adata.n_vars,
                    has_spatial='spatial' in adata.obsm,
                    file_path=str(output_file.relative_to(BASE_DIR.parent.parent)),
                    file_size_mb=round(output_file.stat().st_size / 1024 / 1024, 2)
                )
                
                print(f"    ✓ 完成: {adata.n_obs} cells, {adata.n_vars} genes")
                
            except Exception as e:
                print(f"    ❌ 创建 {demo_name} 失败: {e}")
    
    except ImportError:
        print("❌ squidpy未安装，跳过演示数据集")

def create_datasets_readme():
    """创建数据集说明文档"""
    readme_content = f"""# 真实空间转录组数据集

本目录包含用于comprehensive testing的真实空间转录组数据集。

## 数据集来源

### 1. Squidpy内置数据集
- 经过验证的标准数据集
- 格式标准化，可直接使用
- 包含多种空间转录组技术

### 2. 10X Genomics公开数据集  
- Visium技术的官方数据
- 高质量的人类和小鼠组织数据
- 包含完整的表达和空间信息

### 3. 已发表论文数据
- 来自高影响因子期刊的原始数据
- 代表各种空间转录组技术的最新进展
- 用于验证方法在真实数据上的性能

## 数据集列表

详细的数据集信息请查看 `datasets_metadata.json` 文件。

## 使用说明

所有数据集都保存为 `.h5ad` 格式，可直接用scanpy加载：

```python
import scanpy as sc

# 加载数据集
adata = sc.read_h5ad('squidpy_seqfish.h5ad')

# 查看基本信息
print(f"细胞数: {{adata.n_obs}}")
print(f"基因数: {{adata.n_vars}}")
print(f"是否有空间坐标: {{'spatial' in adata.obsm}}")
```

## 引用信息

使用这些数据集时请引用相应的原始论文，详见metadata中的citation字段。

## 更新时间

数据集最后更新时间: {time.strftime('%Y-%m-%d %H:%M:%S')}
"""
    
    readme_file = BASE_DIR / 'README.md'
    with open(readme_file, 'w', encoding='utf-8') as f:
        f.write(readme_content)
    
    print(f"✓ 数据集说明文档保存至: {readme_file}")

def generate_datasets_summary():
    """生成数据集统计摘要"""
    print("\n=== 生成数据集摘要 ===")
    
    summary = []
    h5ad_files = list(BASE_DIR.glob('*.h5ad'))
    
    for h5ad_file in h5ad_files:
        try:
            adata = sc.read_h5ad(h5ad_file)
            summary.append({
                'dataset': h5ad_file.stem,
                'file_name': h5ad_file.name,
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'has_spatial': 'spatial' in adata.obsm,
                'spatial_dims': adata.obsm['spatial'].shape[1] if 'spatial' in adata.obsm else None,
                'file_size_mb': round(h5ad_file.stat().st_size / 1024 / 1024, 2),
                'sparsity': round((adata.X == 0).sum() / adata.X.size, 3) if hasattr(adata.X, 'size') else 'N/A',
                'cell_types': len(adata.obs.columns),
                'gene_metadata': len(adata.var.columns),
            })
        except Exception as e:
            summary.append({
                'dataset': h5ad_file.stem,
                'file_name': h5ad_file.name,
                'error': str(e),
                'file_size_mb': round(h5ad_file.stat().st_size / 1024 / 1024, 2),
            })
    
    if summary:
        summary_df = pd.DataFrame(summary)
        summary_file = BASE_DIR / 'datasets_summary.csv'
        summary_df.to_csv(summary_file, index=False)
        
        print(f"✓ 数据集摘要保存至: {summary_file}")
        print(f"\n数据集统计:")
        print(f"  总计 h5ad 文件: {len(h5ad_files)}")
        if not summary_df.empty:
            if 'error' in summary_df.columns:
                valid_datasets = summary_df[~summary_df['error'].notna()]
            else:
                valid_datasets = summary_df
            if not valid_datasets.empty:
                print(f"  有效数据集: {len(valid_datasets)}")
                print(f"  细胞总数: {valid_datasets['n_cells'].sum():,}")
                print(f"  平均基因数: {valid_datasets['n_genes'].mean():.0f}")
                print(f"  总数据大小: {valid_datasets['file_size_mb'].sum():.1f} MB")
                print(f"  有空间信息: {valid_datasets['has_spatial'].sum()}/{len(valid_datasets)}")
    else:
        print("⚠️  未找到任何h5ad文件")

def main():
    """主函数"""
    print("开始下载真实空间转录组数据集...")
    print(f"目标目录: {BASE_DIR}")
    
    # 检查依赖
    missing_deps = []
    try:
        import scanpy as sc
        print(f"✓ scanpy版本: {sc.__version__}")
    except ImportError:
        missing_deps.append('scanpy')
    
    try:
        import squidpy as sq  
        print(f"✓ squidpy版本: {sq.__version__}")
    except ImportError:
        print("⚠️  squidpy未安装，部分功能受限")
    
    if missing_deps:
        print(f"❌ 缺少依赖: {', '.join(missing_deps)}")
        print("请安装: pip install scanpy squidpy")
        return
    
    # 执行下载任务
    try:
        download_squidpy_datasets()
        download_10x_genomics_datasets() 
        download_published_datasets()
        download_demo_datasets()
        
        # 生成文档和摘要
        save_metadata()
        create_datasets_readme()
        generate_datasets_summary()
        
        print(f"\n🎉 所有数据集处理完成！")
        print(f"数据保存在: {BASE_DIR}")
        
        # 显示最终统计
        h5ad_count = len(list(BASE_DIR.glob('*.h5ad')))
        total_size = sum(f.stat().st_size for f in BASE_DIR.glob('*')) / 1024 / 1024
        
        print(f"\n最终统计:")
        print(f"  h5ad文件数: {h5ad_count}")
        print(f"  总占用空间: {total_size:.1f} MB")
        print(f"  元数据记录: {len(DATASETS_METADATA)} 个数据集")
        
    except KeyboardInterrupt:
        print("\n❌ 用户中断下载")
    except Exception as e:
        print(f"\n❌ 下载过程中发生错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()