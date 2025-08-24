# 真实Visium空间转录组数据集

本目录包含用于测试ChatSpatial工具的真实格式空间转录组数据集。

## 数据集概览

**总计**: 43个数据集，5.8GB  
**生成时间**: 2025-08-24  

## 核心数据集分类

### 🔬 真实Visium数据集

#### squidpy_visium/
- **squidpy_visium.h5ad**: 684 spots × 18078 genes (89.9MB)
  - 来源: Squidpy内置数据集  
  - 技术: Visium
  - 描述: 鼠脑组织切片，经过质量验证

### 🧬 合成Visium格式数据集

#### 大型基准数据集
- **synthetic_human_brain.h5ad**: 2000 spots × 15000 genes (230.5MB)
  - 模拟人类大脑皮质，具有真实Visium特征
  - 70%稀疏度，六边形spot排列

- **synthetic_mouse_kidney.h5ad**: 1500 spots × 12000 genes (138.6MB)
  - 模拟小鼠肾脏组织，包含肾脏特异性基因模式

- **synthetic_human_heart.h5ad**: 1200 spots × 10000 genes (92.6MB)
  - 模拟人类心脏组织，心脏特异性基因表达模式

#### 测试用小型数据集
- **synthetic_small_test.h5ad**: 500 spots × 5000 genes (19.6MB)
  - 快速测试用小型数据集
  - 保留完整Visium数据结构

### 📊 参考数据集

#### scanpy_paul15/
- **scanpy_paul15.h5ad**: 2730 cells × 3451 genes (36.2MB)
  - 造血干细胞分化数据
  - 添加了合成空间坐标用于空间分析测试

#### 其他核心数据集 (core/)
- **slideseq_MOp_1217.h5ad**: 高分辨率空间转录组 (1.5GB)
- **ST_mouse_brain.h5ad**: Spatial Transcriptomics技术数据
- 多个pancreas相关数据集用于特定分析测试

## 数据格式说明

所有数据集都采用标准的AnnData h5ad格式，包含：

### 必需组件
- **X**: 基因表达矩阵 (spots × genes)
- **obs**: Spot/细胞注释信息
- **var**: 基因注释信息
- **obsm['spatial']**: 空间坐标 (2D)

### 元数据标准
每个数据集包含 `uns['dataset_info']`：
```python
{
    'name': '数据集名称',
    'source': '数据来源',
    'technology': '技术平台',
    'description': '详细描述',
    'download_date': '下载日期',
    'n_spots': '观测点数量',
    'n_genes': '基因数量'
}
```

## 使用示例

### 基础加载
```python
import scanpy as sc
from pathlib import Path

# 设置基础路径
base_dir = Path("datasets/real_datasets")

# 加载小型测试数据集
adata = sc.read_h5ad(base_dir / "synthetic_small_test/synthetic_small_test.h5ad")
print(f"加载数据集: {adata.n_obs} spots, {adata.n_vars} genes")

# 检查空间坐标
if 'spatial' in adata.obsm:
    print(f"空间坐标维度: {adata.obsm['spatial'].shape}")
```

### ChatSpatial集成测试
```python
# 使用不同大小的数据集测试性能
test_datasets = {
    'small': 'synthetic_small_test/synthetic_small_test.h5ad',
    'medium': 'squidpy_visium/squidpy_visium.h5ad', 
    'large': 'synthetic_human_brain/synthetic_human_brain.h5ad'
}

for size, path in test_datasets.items():
    adata = sc.read_h5ad(base_dir / path)
    # 运行ChatSpatial分析
    print(f"{size} dataset: {adata.shape}")
```

## 推荐使用策略

### 快速开发测试
推荐使用：
- `synthetic_small_test.h5ad` - 最小数据集，快速验证
- `squidpy_visium.h5ad` - 真实数据，中等大小

### 功能验证测试  
推荐使用：
- `synthetic_human_brain.h5ad` - 大型数据集，性能测试
- `scanpy_paul15.h5ad` - 经典参考数据集

### 生产环境测试
推荐使用：
- `core/` 目录下的真实大型数据集
- 多个数据集组合测试

## 数据质量保证

✅ **验证通过**:
- 所有数据集都能正常加载
- 空间坐标格式正确  
- 基因表达矩阵非空
- 元数据完整

✅ **兼容性测试**:
- Scanpy 1.11.0
- AnnData 0.11.4
- ChatSpatial工具链

## 故障排除

### 常见问题

1. **内存不足**
   - 使用小型数据集先测试
   - 考虑使用 `backed=True` 模式

2. **文件损坏**
   - 重新运行下载脚本
   - 检查文件大小是否正确

3. **坐标格式错误**
   - 所有数据集的spatial坐标都是float64格式
   - 坐标存储在 `adata.obsm['spatial']`

### 重新生成数据集

如需重新生成数据集：
```bash
# 重新下载和生成所有数据集
python download_real_visium_data.py
python create_additional_real_datasets.py

# 快速统计摘要
python quick_dataset_summary.py
```

## 贡献指南

添加新数据集时，请确保：

1. 遵循标准h5ad格式
2. 包含完整的`dataset_info`元数据
3. 空间坐标存储在正确位置
4. 更新本README文档

---

**维护者**: ChatSpatial开发团队  
**最后更新**: 2025-08-24  
**数据集版本**: v1.0