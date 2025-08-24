# 空间转录组数据格式处理指南

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
    expr_df = pd.read_csv(expr_file, sep='\t', index_col=0)
    coords_df = pd.read_csv(coord_file, sep='\t', index_col=0)
    
    # 创建AnnData
    adata = sc.AnnData(X=expr_df.values)
    adata.obsm['spatial'] = coords_df.values
    
    return adata

def process_mtx_10x(mtx_file, barcode_file, features_file, positions_file):
    # 读取稀疏矩阵并转置
    X = scipy.io.mmread(mtx_file).T.tocsr()
    
    # 读取barcode和features
    barcodes = pd.read_csv(barcode_file, header=None)[0].tolist()
    features = pd.read_csv(features_file, sep='\t', header=None)
    
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
