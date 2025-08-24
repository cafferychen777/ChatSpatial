#!/usr/bin/env python3
"""
空间转录组数据格式处理演示脚本
展示如何处理Slide-seq和ST数据的不同格式
"""

import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import scipy.io
import scipy.sparse
import h5py
import warnings
warnings.filterwarnings('ignore')

# 设置scanpy
sc.settings.verbosity = 1

BASE_DIR = Path(__file__).parent / 'datasets' / 'real_datasets'

def process_slideseq_csv_format():
    """
    演示处理Slide-seq CSV格式数据
    格式: expression.csv (genes x cells) + coordinates.csv
    """
    print("\n" + "="*60)
    print("处理 Slide-seq CSV 格式数据")
    print("="*60)
    
    # 模拟创建Slide-seq CSV格式数据
    demo_dir = BASE_DIR / 'format_demos'
    demo_dir.mkdir(exist_ok=True)
    
    # 创建表达矩阵 CSV (genes x cells)
    n_genes, n_cells = 1000, 2000
    np.random.seed(42)
    
    # Slide-seq特征：高稀疏性
    expression_matrix = np.random.negative_binomial(n=2, p=0.85, size=(n_genes, n_cells))
    
    gene_names = [f'Gene_{i:04d}' for i in range(n_genes)]
    cell_barcodes = [f'Cell_{i:06d}' for i in range(n_cells)]
    
    expr_df = pd.DataFrame(expression_matrix, index=gene_names, columns=cell_barcodes)
    expr_file = demo_dir / 'slideseq_expression.csv'
    expr_df.to_csv(expr_file)
    print(f"创建表达矩阵文件: {expr_file}")
    
    # 创建坐标文件
    coords_df = pd.DataFrame({
        'barcode': cell_barcodes,
        'xcoord': np.random.uniform(-1000, 1000, n_cells),
        'ycoord': np.random.uniform(-1000, 1000, n_cells)
    })
    coords_file = demo_dir / 'slideseq_coordinates.csv'
    coords_df.to_csv(coords_file, index=False)
    print(f"创建坐标文件: {coords_file}")
    
    # 处理数据
    print("\n处理步骤:")
    
    # 1. 读取表达矩阵
    print("1. 读取表达矩阵...")
    expr_df = pd.read_csv(expr_file, index_col=0)
    print(f"   表达矩阵形状: {expr_df.shape} (genes x cells)")
    
    # 2. 转置为 cells x genes 
    print("2. 转置矩阵为 cells x genes...")
    X = expr_df.T.values
    print(f"   转置后形状: {X.shape}")
    
    # 3. 读取坐标
    print("3. 读取坐标信息...")
    coords_df = pd.read_csv(coords_file)
    print(f"   坐标形状: {coords_df.shape}")
    
    # 4. 匹配barcode
    print("4. 匹配barcode...")
    coords_df = coords_df.set_index('barcode').reindex(expr_df.columns)
    
    # 5. 创建AnnData
    print("5. 创建AnnData对象...")
    adata = sc.AnnData(
        X=X,
        obs=pd.DataFrame({
            'barcode': expr_df.columns,
            'cell_type': 'unknown'
        }, index=expr_df.columns),
        var=pd.DataFrame({
            'gene_symbol': expr_df.index
        }, index=expr_df.index)
    )
    
    # 6. 添加空间坐标
    adata.obsm['spatial'] = coords_df[['xcoord', 'ycoord']].values
    
    # 7. 添加元信息
    adata.uns['spatial'] = {
        'technology': {
            'technology': 'Slide-seq',
            'bead_diameter': 10,
            'resolution': 'subcellular'
        },
        'processing': {
            'format': 'CSV',
            'genes_x_cells': True
        }
    }
    
    # 8. 保存处理后数据
    output_file = demo_dir / 'slideseq_processed.h5ad'
    adata.write(output_file)
    print(f"6. 保存处理后数据: {output_file}")
    
    print(f"✅ Slide-seq CSV格式处理完成")
    print(f"   - Spots: {adata.n_obs:,}")
    print(f"   - Genes: {adata.n_vars:,}")
    print(f"   - 稀疏性: {(adata.X == 0).sum() / adata.X.size:.3f}")
    
    return adata

def process_st_tsv_format():
    """
    演示处理ST TSV格式数据  
    格式: expression.tsv (spots x genes) + coordinates.tsv
    """
    print("\n" + "="*60)
    print("处理 ST TSV 格式数据")
    print("="*60)
    
    demo_dir = BASE_DIR / 'format_demos'
    demo_dir.mkdir(exist_ok=True)
    
    # 创建ST格式数据
    n_spots, n_genes = 300, 2000
    np.random.seed(42)
    
    # ST特征：中等稀疏性，规则的spot位置
    expression_matrix = np.random.negative_binomial(n=5, p=0.6, size=(n_spots, n_genes))
    
    spot_names = [f'Spot_{i:04d}' for i in range(n_spots)]
    gene_names = [f'ENSG{i:08d}' for i in range(n_genes)]  # Human格式
    
    # 1. 创建表达文件 (spots x genes)
    expr_df = pd.DataFrame(expression_matrix, index=spot_names, columns=gene_names)
    expr_file = demo_dir / 'st_expression.tsv'
    expr_df.to_csv(expr_file, sep='\t')
    print(f"创建表达矩阵文件: {expr_file}")
    
    # 2. 创建坐标文件 (规则网格，模拟ST spots)
    grid_size = int(np.ceil(np.sqrt(n_spots)))
    x_coords = []
    y_coords = []
    for i in range(grid_size):
        for j in range(grid_size):
            if len(x_coords) < n_spots:
                x_coords.append(j * 100)  # 100μm间距
                y_coords.append(i * 100)
    
    # 确保坐标数量与spots数量完全匹配
    x_coords = x_coords[:n_spots]
    y_coords = y_coords[:n_spots]
    
    coords_df = pd.DataFrame({
        'x': x_coords,
        'y': y_coords
    }, index=spot_names)
    coords_file = demo_dir / 'st_coordinates.tsv'
    coords_df.to_csv(coords_file, sep='\t')
    print(f"创建坐标文件: {coords_file}")
    
    # 处理数据
    print("\n处理步骤:")
    
    # 1. 读取表达数据
    print("1. 读取表达矩阵...")
    expr_df = pd.read_csv(expr_file, sep='\t', index_col=0)
    print(f"   表达矩阵形状: {expr_df.shape} (spots x genes)")
    
    # 2. 读取坐标
    print("2. 读取坐标信息...")
    coords_df = pd.read_csv(coords_file, sep='\t', index_col=0)
    print(f"   坐标形状: {coords_df.shape}")
    
    # 3. ST数据已经是 spots x genes 格式，无需转置
    X = expr_df.values
    print(f"   数据矩阵形状: {X.shape}")
    
    # 4. 创建AnnData
    print("3. 创建AnnData对象...")
    adata = sc.AnnData(
        X=X,
        obs=pd.DataFrame({
            'spot_id': expr_df.index,
            'tissue': 'breast_cancer',
            'in_tissue': True
        }, index=expr_df.index),
        var=pd.DataFrame({
            'gene_symbol': expr_df.columns,
            'feature_type': 'Gene Expression'
        }, index=expr_df.columns)
    )
    
    # 5. 添加空间坐标
    adata.obsm['spatial'] = coords_df.values
    
    # 6. 添加元信息
    adata.uns['spatial'] = {
        'technology': {
            'technology': 'Spatial Transcriptomics',
            'spot_diameter': 100,
            'resolution': 'multi-cellular'
        },
        'processing': {
            'format': 'TSV',
            'spots_x_genes': True
        }
    }
    
    # 7. 保存
    output_file = demo_dir / 'st_processed.h5ad'
    adata.write(output_file)
    print(f"4. 保存处理后数据: {output_file}")
    
    print(f"✅ ST TSV格式处理完成")
    print(f"   - Spots: {adata.n_obs:,}")
    print(f"   - Genes: {adata.n_vars:,}")
    print(f"   - 稀疏性: {(adata.X == 0).sum() / adata.X.size:.3f}")
    
    return adata

def process_mtx_format():
    """
    演示处理MTX格式数据 (类似10X Genomics)
    格式: matrix.mtx + barcodes.tsv + features.tsv + spatial/
    """
    print("\n" + "="*60)
    print("处理 MTX 格式数据 (10X style)")
    print("="*60)
    
    demo_dir = BASE_DIR / 'format_demos' / 'mtx_demo'
    demo_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建模拟10X格式数据
    n_cells, n_genes = 1500, 2500
    np.random.seed(42)
    
    # 创建稀疏表达矩阵
    density = 0.1  # 10% 非零值
    X_sparse = scipy.sparse.random(n_genes, n_cells, density=density, format='csr')
    X_sparse.data = np.random.negative_binomial(n=3, p=0.3, size=X_sparse.nnz)
    
    # 1. 保存MTX文件
    mtx_file = demo_dir / 'matrix.mtx'
    scipy.io.mmwrite(mtx_file, X_sparse)
    print(f"创建MTX文件: {mtx_file}")
    
    # 2. 创建barcode文件
    barcodes = [f'CELL_{i:06d}-1' for i in range(n_cells)]
    barcode_file = demo_dir / 'barcodes.tsv'
    with open(barcode_file, 'w') as f:
        for bc in barcodes:
            f.write(f"{bc}\n")
    print(f"创建barcode文件: {barcode_file}")
    
    # 3. 创建features文件
    features_data = []
    for i in range(n_genes):
        gene_id = f'ENSG{i:08d}'
        gene_name = f'Gene_{i:04d}'
        features_data.append([gene_id, gene_name, 'Gene Expression'])
    
    features_file = demo_dir / 'features.tsv'
    features_df = pd.DataFrame(features_data, columns=['gene_id', 'gene_symbol', 'feature_type'])
    features_df.to_csv(features_file, sep='\t', header=False, index=False)
    print(f"创建features文件: {features_file}")
    
    # 4. 创建空间信息
    spatial_dir = demo_dir / 'spatial'
    spatial_dir.mkdir(exist_ok=True)
    
    # tissue_positions.csv
    positions = []
    for i, bc in enumerate(barcodes):
        in_tissue = np.random.choice([0, 1], p=[0.3, 0.7])  # 70% in tissue
        array_row = np.random.randint(0, 100)
        array_col = np.random.randint(0, 100) 
        pxl_row_in_fullres = array_row * 100 + np.random.randint(-50, 50)
        pxl_col_in_fullres = array_col * 100 + np.random.randint(-50, 50)
        positions.append([bc, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres])
    
    positions_file = spatial_dir / 'tissue_positions_list.csv'
    positions_df = pd.DataFrame(positions, columns=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres'])
    positions_df.to_csv(positions_file, header=False, index=False)
    print(f"创建空间位置文件: {positions_file}")
    
    # 处理数据
    print("\n处理步骤:")
    
    # 1. 读取MTX文件
    print("1. 读取MTX文件...")
    X = scipy.io.mmread(mtx_file).T.tocsr()  # 转置: genes x cells -> cells x genes
    print(f"   矩阵形状: {X.shape} (cells x genes)")
    
    # 2. 读取barcode
    print("2. 读取barcode...")
    with open(barcode_file, 'r') as f:
        barcodes = [line.strip() for line in f]
    print(f"   Barcode数量: {len(barcodes)}")
    
    # 3. 读取features  
    print("3. 读取features...")
    features_df = pd.read_csv(features_file, sep='\t', header=None, names=['gene_id', 'gene_symbol', 'feature_type'])
    print(f"   基因数量: {len(features_df)}")
    
    # 4. 读取空间位置
    print("4. 读取空间位置...")
    positions_df = pd.read_csv(positions_file, header=None, names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres'])
    positions_df = positions_df.set_index('barcode').reindex(barcodes)
    print(f"   位置信息形状: {positions_df.shape}")
    
    # 5. 创建AnnData
    print("5. 创建AnnData...")
    adata = sc.AnnData(
        X=X,
        obs=pd.DataFrame({
            'barcode': barcodes,
            'in_tissue': positions_df['in_tissue'].astype(bool),
            'array_row': positions_df['array_row'],
            'array_col': positions_df['array_col']
        }, index=barcodes),
        var=features_df.set_index('gene_id')
    )
    
    # 6. 添加空间坐标 (使用像素坐标)
    adata.obsm['spatial'] = positions_df[['pxl_col_in_fullres', 'pxl_row_in_fullres']].values
    
    # 7. 添加元信息
    adata.uns['spatial'] = {
        'technology': {
            'technology': 'Visium',
            'spot_diameter': 55,
            'resolution': 'multi-cellular'
        },
        'processing': {
            'format': 'MTX',
            'source': '10X_style'
        }
    }
    
    # 8. 保存
    output_file = demo_dir / 'mtx_processed.h5ad'
    adata.write(output_file)
    print(f"6. 保存处理后数据: {output_file}")
    
    print(f"✅ MTX格式处理完成")
    print(f"   - Spots: {adata.n_obs:,}")
    print(f"   - Genes: {adata.n_vars:,}")
    print(f"   - 稀疏性: {1.0 - (adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1])):.3f}")
    
    return adata

def demonstrate_data_validation():
    """演示数据验证过程"""
    print("\n" + "="*60)
    print("数据验证演示")
    print("="*60)
    
    # 加载现有数据集进行演示
    dataset_files = [
        'slideseq_v2_hippocampus.h5ad',
        'st_human_breast_cancer.h5ad'
    ]
    
    for dataset_file in dataset_files:
        dataset_path = BASE_DIR / dataset_file
        if not dataset_path.exists():
            print(f"⚠️ 数据集不存在: {dataset_file}")
            continue
            
        print(f"\n验证数据集: {dataset_file}")
        print("-" * 40)
        
        # 加载数据
        adata = sc.read_h5ad(dataset_path)
        
        # 1. 基本结构验证
        print("1. 基本结构验证:")
        print(f"   ✅ Spots: {adata.n_obs:,}")
        print(f"   ✅ Genes: {adata.n_vars:,}")
        print(f"   ✅ 数据类型: {type(adata.X)}")
        
        # 2. 空间坐标验证
        print("2. 空间坐标验证:")
        if 'spatial' in adata.obsm:
            coords = adata.obsm['spatial']
            print(f"   ✅ 空间坐标存在: {coords.shape}")
            print(f"   ✅ 坐标范围: X[{coords[:, 0].min():.1f}, {coords[:, 0].max():.1f}], Y[{coords[:, 1].min():.1f}, {coords[:, 1].max():.1f}]")
            
            # 检查异常值
            if np.any(np.isnan(coords)):
                print("   ⚠️ 坐标包含NaN值")
            else:
                print("   ✅ 坐标无NaN值")
        else:
            print("   ❌ 缺少空间坐标")
        
        # 3. 表达数据验证
        print("3. 表达数据验证:")
        if scipy.sparse.issparse(adata.X):
            sparsity = 1.0 - (adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]))
            has_negative = (adata.X.data < 0).any()
        else:
            sparsity = (adata.X == 0).sum() / adata.X.size  
            has_negative = (adata.X < 0).any()
        
        print(f"   ✅ 稀疏性: {sparsity:.3f}")
        if has_negative:
            print("   ⚠️ 包含负值")
        else:
            print("   ✅ 无负值")
        
        # 4. 技术信息验证
        print("4. 技术信息验证:")
        if 'spatial' in adata.uns:
            spatial_info = adata.uns['spatial']
            if 'technology' in spatial_info:
                tech = spatial_info['technology']
                print(f"   ✅ 技术类型: {tech.get('technology', 'Unknown')}")
                print(f"   ✅ 分辨率: {tech.get('resolution', 'Unknown')}")
            else:
                print("   ⚠️ 缺少技术信息")
        else:
            print("   ⚠️ 缺少空间元信息")

def create_format_comparison_report():
    """创建格式对比报告"""
    print("\n" + "="*60) 
    print("生成格式对比报告")
    print("="*60)
    
    demo_dir = BASE_DIR / 'format_demos'
    
    report_content = """# 空间转录组数据格式处理指南

## 支持的数据格式

### 1. Slide-seq CSV格式
**文件结构:**
```
slideseq_data/
├── expression.csv          # 基因表达矩阵 (genes × cells)
└── coordinates.csv         # 空间坐标 (barcode, xcoord, ycoord)
```

**特征:**
- 高分辨率 (~10μm)
- 高稀疏性 (>90%)
- 基因 × 细胞 矩阵格式

**处理步骤:**
1. 读取表达矩阵并转置为 cells × genes
2. 读取坐标信息
3. 通过barcode匹配表达和坐标数据
4. 创建AnnData对象

### 2. ST TSV格式  
**文件结构:**
```
st_data/
├── expression.tsv          # 基因表达矩阵 (spots × genes)  
└── coordinates.tsv         # 空间坐标 (spot_id, x, y)
```

**特征:**
- 中等分辨率 (~100μm)
- 中等稀疏性 (~65%)
- Spots × 基因 矩阵格式

**处理步骤:**
1. 读取表达矩阵 (已经是正确格式)
2. 读取坐标信息
3. 直接创建AnnData对象

### 3. MTX格式 (10X style)
**文件结构:**
```
10x_data/
├── matrix.mtx              # 稀疏表达矩阵
├── barcodes.tsv            # 细胞barcode
├── features.tsv            # 基因信息
└── spatial/
    └── tissue_positions_list.csv  # 空间位置
```

**特征:**  
- 稀疏矩阵存储
- 高效内存使用
- 标准化格式

**处理步骤:**
1. 读取MTX文件并转置
2. 读取barcode和features
3. 读取空间位置信息
4. 匹配所有信息创建AnnData

## 数据验证清单

### ✅ 必需验证项
- [ ] 数据矩阵非空
- [ ] 空间坐标存在且维度正确
- [ ] 无NaN或无穷值
- [ ] 基因和细胞名称唯一

### ⚠️ 推荐验证项  
- [ ] 合理的稀疏性范围
- [ ] 空间坐标分布合理
- [ ] 表达值非负
- [ ] 技术元信息完整

## 常见问题处理

### 问题1: 矩阵方向错误
**症状:** 基因数量远大于细胞数量
**解决:** 检查并转置矩阵

### 问题2: 坐标匹配失败  
**症状:** 空间坐标数量与细胞数量不匹配
**解决:** 通过barcode重新索引

### 问题3: 稀疏性异常
**症状:** 稀疏性<10% 或 >99%
**解决:** 检查数据预处理和格式转换

## 代码示例

```python
import scanpy as sc
import pandas as pd
import scipy.io

def process_slideseq_csv(expr_file, coord_file):
    # 读取并转置表达矩阵
    expr_df = pd.read_csv(expr_file, index_col=0)
    X = expr_df.T.values
    
    # 读取坐标
    coords_df = pd.read_csv(coord_file)
    coords_df = coords_df.set_index('barcode').reindex(expr_df.columns)
    
    # 创建AnnData
    adata = sc.AnnData(X=X)
    adata.obsm['spatial'] = coords_df[['xcoord', 'ycoord']].values
    
    return adata

def process_st_tsv(expr_file, coord_file):
    # 直接读取 spots × genes 矩阵
    expr_df = pd.read_csv(expr_file, sep='\\t', index_col=0)
    coords_df = pd.read_csv(coord_file, sep='\\t', index_col=0)
    
    # 创建AnnData
    adata = sc.AnnData(X=expr_df.values)
    adata.obsm['spatial'] = coords_df.values
    
    return adata

def process_mtx_10x(mtx_file, barcode_file, features_file, positions_file):
    # 读取稀疏矩阵并转置
    X = scipy.io.mmread(mtx_file).T.tocsr()
    
    # 读取barcode和features
    barcodes = pd.read_csv(barcode_file, header=None)[0].tolist()
    features = pd.read_csv(features_file, sep='\\t', header=None)
    
    # 读取空间位置
    positions = pd.read_csv(positions_file, header=None)
    positions = positions.set_index(0).reindex(barcodes)
    
    # 创建AnnData
    adata = sc.AnnData(X=X)
    adata.obsm['spatial'] = positions[[4, 5]].values
    
    return adata
```

---
*由 ChatSpatial 数据处理框架生成*
"""
    
    report_file = demo_dir / 'FORMAT_PROCESSING_GUIDE.md'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report_content)
    
    print(f"格式处理指南已保存: {report_file}")

def main():
    """主函数"""
    print("=" * 80)
    print("空间转录组数据格式处理演示")
    print("=" * 80)
    
    # 1. 处理不同格式
    slideseq_data = process_slideseq_csv_format()
    st_data = process_st_tsv_format()  
    mtx_data = process_mtx_format()
    
    # 2. 演示数据验证
    demonstrate_data_validation()
    
    # 3. 创建对比报告
    create_format_comparison_report()
    
    print("\n" + "=" * 80)
    print("格式处理演示完成！")
    print("所有演示文件保存在: datasets/real_datasets/format_demos/")
    print("=" * 80)

if __name__ == "__main__":
    main()