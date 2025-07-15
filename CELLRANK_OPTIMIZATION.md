# CellRank 性能优化指南

## 当前状态

CellRank 已经可以正常工作，但会显示以下警告：
```
WARNING: Unable to import `petsc4py` or `slepc4py`. Using `method='brandts'`
WARNING: For `method='brandts'`, dense matrix is required. Densifying
```

这意味着 CellRank 正在使用较慢的密集矩阵方法，而不是优化的稀疏矩阵方法。

## 可选优化

如果您需要处理大规模数据集（>5000个细胞），建议安装 PETSc 和 SLEPc 来提高性能。

### 方法1：使用 conda（推荐用于 macOS）

```bash
# 创建新的 conda 环境
conda create -n spatial_mcp python=3.10
conda activate spatial_mcp

# 安装 PETSc 和 SLEPc
conda install -c conda-forge petsc4py slepc4py

# 安装其他依赖
pip install -e /Users/apple/Research/SpatialTrans_MCP/chatspatial[cellrank_full]
```

### 方法2：使用 pip（需要编译）

```bash
# 安装 MPI（已完成）
pip install mpi4py

# 安装 PETSc 和 SLEPc（需要较长时间）
pip install petsc petsc4py
pip install slepc slepc4py
```

### 方法3：使用 Homebrew + pip（macOS）

```bash
# 安装系统级 PETSc
brew install petsc

# 设置环境变量
export PETSC_DIR=$(brew --prefix petsc)

# 安装 Python 绑定
pip install petsc4py slepc4py
```

## 性能比较

- **没有 PETSc**：适用于小到中等规模数据集（<5000个细胞）
  - 使用密集矩阵
  - 内存占用较高
  - 计算速度较慢

- **有 PETSc**：适用于大规模数据集（>5000个细胞）
  - 使用稀疏矩阵
  - 内存效率高
  - 计算速度快
  - 支持并行计算

## 当前建议

对于当前的测试数据集（胰腺数据集，约3700个细胞），不安装 PETSc 也完全可以接受。CellRank 会自动使用 Brandts 方法，虽然会慢一些，但结果是一样的。

如果将来需要分析更大的数据集，可以考虑安装这些优化依赖。