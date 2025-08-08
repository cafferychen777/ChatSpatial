# ChatSpatial MCP 功能完整文档

## 概述

ChatSpatial 是一个基于 Model Context Protocol (MCP) 的空间转录组学分析平台，提供了全面的空间转录组学数据分析工具集。本文档详细描述了所有可用的功能、参数和使用方法。

## 核心架构

ChatSpatial MCP 采用分层架构设计，提供了完整的空间转录组学分析能力和强大的错误处理机制。

### MCP 服务器配置
- **服务器名称**: ChatSpatial
- **传输协议**: stdio (默认) / SSE
- **Python 环境**: `/Users/apple/Research/SpatialTrans_MCP/st_mcp_env_py310/bin/python3`
- **主要模块**: 基于 FastMCP 框架的空间转录组学分析适配器

### 架构组件

#### 1. SpatialMCPAdapter (核心适配器)
**功能**: 连接 MCP 协议和空间分析功能的主要适配器
**职责**:
- 工具注册和元数据管理
- 资源和提示处理协调
- MCP 请求路由和响应格式化
- 工具装饰器集成

#### 2. 数据管理层

**SpatialDataManager (抽象接口)**:
- 数据集加载和存储管理
- 分析结果持久化
- 数据集列表和查询
- 支持自定义数据管理实现

**DefaultSpatialDataManager (默认实现)**:
- 内存数据存储 (`data_store`)
- 自动ID生成和管理
- 数据集元数据维护

#### 3. 资源管理层

**SpatialResourceManager**:
- **数据集资源**: `spatial://datasets/{data_id}` 格式
- **可视化资源**: `spatial://visualizations/{viz_id}` 格式
- **分析结果资源**: `spatial://results/{data_id}/{result_type}` 格式
- **资源缓存**: 可视化图像缓存机制
- **内容提供**: 支持文本和二进制内容

**资源类型**:
```python
MCPResource:
  - uri: 资源唯一标识符
  - name: 资源显示名称
  - mime_type: MIME 类型 (application/x-anndata, image/png, application/json)
  - description: 资源描述
  - metadata: 扩展元数据
  - content_provider: 内容提供函数
```

#### 4. 提示管理层

**SpatialPromptManager**:
- **预定义工作流**: 10+ 个分析提示模板
- **参数映射**: 提示参数到工具参数的自动转换
- **工作流执行**: 复杂分析流程的简化接口

**可用提示**:
- `analyze-spatial-expression`: 空间基因表达分析
- `find-cell-types`: 细胞类型识别
- `compare-conditions`: 条件比较分析
- `generate-visualization`: 可视化生成
- `quality-control`: 质量控制
- `batch-correction`: 批次校正
- `spatial-clustering`: 空间聚类
- `trajectory-inference`: 轨迹推断
- `spatial-deconvolution`: 空间反卷积

#### 5. 错误处理系统

**多层错误处理架构**:

**mcp_tool_error_handler** (主要错误处理器):
- 捕获所有工具执行异常
- 标准化错误格式 (ToolResult)
- **关键图像处理**: 保护 FastMCP Image 对象直接返回
- 支持同步和异步函数
- 可配置回溯信息

**mcp_pydantic_error_handler** (参数验证):
- Pydantic 模型自动验证
- 参数类型检查和转换
- 详细验证错误信息
- 与主错误处理器集成

**ProcessingError 异常体系**:
```python
SpatialMCPError (基类)
├── DataNotFoundError (数据未找到)
├── InvalidParameterError (参数无效)
├── ProcessingError (处理失败)
└── DataCompatibilityError (数据兼容性)
```

**MCP 错误格式化**:
```python
ErrorType 枚举:
- INVALID_PARAMS (-32602)
- DATASET_NOT_FOUND (-32001)
- INVALID_DATA_FORMAT (-32002)
- ANALYSIS_FAILED (-32003)
- VISUALIZATION_ERROR (-32004)
- REFERENCE_DATA_ERROR (-32005)
- INTERNAL_ERROR (-32603)
```

#### 6. 参数验证系统

**manual_parameter_validation**:
- 自定义验证函数装饰器
- 支持复杂参数验证逻辑
- 与 Pydantic 验证互补

**专用验证器**:
- `validate_analysis_params`: 分析参数验证
- `validate_visualization_params`: 可视化参数验证
- `validate_spatial_analysis_params`: 空间分析参数验证

#### 7. 工具元数据系统

**MCPToolMetadata**:
```python
- name: 工具名称
- title: 显示标题
- description: 功能描述
- read_only_hint: 只读提示
- idempotent_hint: 幂等性提示
- open_world_hint: 外部资源访问提示
```

**工具注解系统** (TOOL_ANNOTATIONS):
- 为每个工具提供行为提示
- 支持 MCP 客户端优化
- 包含资源使用和副作用信息

### 架构优势

#### 1. 模块化设计
- 清晰的职责分离
- 可扩展的组件架构
- 支持自定义实现

#### 2. 强大的错误处理
- 多层错误捕获和处理
- 标准化错误格式
- 详细的错误信息和回溯

#### 3. 资源管理
- 统一的资源访问接口
- 高效的缓存机制
- 支持多种内容类型

#### 4. 工作流简化
- 预定义分析提示
- 参数自动映射
- 复杂流程简化

#### 5. 类型安全
- Pydantic 模型验证
- 类型提示支持
- 运行时类型检查

### 服务器启动和配置

#### 启动方式
```bash
# 方式1: 直接运行模块
python -m chatspatial

# 方式2: 使用命令行参数
python -m chatspatial --transport stdio --host localhost --port 8000

# 方式3: 通过 server.py
python chatspatial/server.py --transport stdio
```

#### 配置选项
- **transport**: `stdio` (默认) 或 `sse`
- **host**: 服务器主机地址 (SSE 模式)
- **port**: 服务器端口 (SSE 模式)
- **log_level**: 日志级别

#### MCP 客户端集成
ChatSpatial MCP 可以与以下客户端集成:
- **Claude Desktop**: 主要推荐客户端
- **其他 MCP 兼容客户端**: 支持标准 MCP 协议

### 数据流架构

#### 1. 数据加载流程
```
用户请求 → load_data 工具 → SpatialDataManager.load_dataset()
→ 数据验证和预处理 → 数据存储 → 资源注册 → 返回数据ID
```

#### 2. 分析执行流程
```
分析请求 → 参数验证 (Pydantic + 自定义) → 数据检索
→ 分析工具执行 → 结果存储 → 资源创建 → 返回结果
```

#### 3. 可视化流程
```
可视化请求 → 参数验证 → 数据和结果检索 → 图表生成
→ 图像转换 → 资源缓存 → 返回 Image 对象
```

#### 4. 错误处理流程
```
异常发生 → 错误类型识别 → 错误格式化 (MCP 标准)
→ 日志记录 → 标准化错误响应 → 客户端显示
```

## 完整功能列表

### 1. 数据加载与预处理

#### 1.1 load_data - 数据加载
**功能**: 加载空间转录组学数据文件
**支持格式**: 
- 10x Visium (h5, h5ad)
- Slide-seq
- MERFISH
- seqFISH
- 通用 h5ad 格式

**参数**:
```python
file_path: str  # 数据文件路径
data_type: Literal["10x_visium", "slide_seq", "merfish", "seqfish", "other", "h5ad", "auto"]
dataset_name: Optional[str]  # 数据集名称
description: Optional[str]   # 数据集描述
```

#### 1.2 preprocess_data - 数据预处理
**功能**: 标准化和预处理空间转录组学数据
**预处理步骤**:
1. 质量控制指标计算
2. 细胞和基因过滤
3. 数据标准化
4. 高变基因识别
5. 降维分析 (PCA, UMAP)
6. 聚类分析

**参数**:
```python
# 过滤参数
min_genes: int = 200          # 每个细胞最少基因数
max_genes: int = 5000         # 每个细胞最多基因数
min_cells: int = 3            # 每个基因最少细胞数
mt_threshold: float = 20.0    # 线粒体基因比例阈值

# 标准化参数
normalization: Literal["log", "sct"] = "log"  # 标准化方法
target_sum: float = 1e4       # 标准化目标总数

# 降维参数
n_pcs: int = 50              # PCA 主成分数
n_neighbors: int = 15        # UMAP 邻居数
min_dist: float = 0.1        # UMAP 最小距离

# 聚类参数
clustering_resolution: float = 0.5  # Leiden 聚类分辨率
```

**高级预处理选项**:
- **RESOLVI**: 深度学习去噪和偏差校正
- **GLM-PCA**: 广义线性模型主成分分析
- **Pearson Residuals**: 皮尔逊残差标准化

### 2. 数据可视化

#### 2.1 visualize_data - 综合可视化
**功能**: 创建各种类型的空间转录组学可视化图表

**可视化类型**:
```python
plot_type: Literal[
    "spatial",           # 空间表达图
    "umap",             # UMAP 降维图
    "pca",              # PCA 图
    "violin",           # 小提琴图
    "heatmap",          # 热图
    "dotplot",          # 点图
    "spatial_analysis", # 空间分析结果
    "trajectory",       # 轨迹分析
    "velocity",         # RNA 速度
    "communication",    # 细胞通讯
    "deconvolution",    # 反卷积结果
    "domains",          # 空间域
    "ripley"           # Ripley 函数
]
```

**参数**:
```python
# 基本参数
feature: Optional[str]        # 要可视化的特征/基因
cluster_key: str = "leiden"   # 聚类结果键

# 图形参数
figure_size: Optional[Tuple[int, int]]  # 图形尺寸
dpi: int = 100               # 分辨率
alpha: float = 0.8           # 透明度
colormap: str = "viridis"    # 颜色映射

# 颜色参数
vmin: Optional[float]        # 颜色范围最小值
vmax: Optional[float]        # 颜色范围最大值
color_scale: Literal["linear", "log", "sqrt"] = "linear"

# 显示参数
show_legend: bool = True     # 显示图例
show_colorbar: bool = True   # 显示颜色条
show_axes: bool = True       # 显示坐标轴
```

### 3. 细胞类型注释

#### 3.1 annotate_cell_types - 细胞类型注释
**功能**: 使用多种方法进行细胞类型注释

**支持方法**:
```python
method: Literal[
    "marker_genes",    # 标记基因方法
    "correlation",     # 相关性方法
    "tangram",        # Tangram 空间映射
    "scanvi",         # scANVI 半监督学习
    "cellassign",     # CellAssign 概率分配
    "mllmcelltype",   # 大语言模型注释
    "supervised",     # 监督学习
    "popv",           # PopV 方法
    "gptcelltype",    # GPT 细胞类型
    "scrgcl"          # scRGCL 方法
]
```

**参数**:
```python
# 通用参数
marker_genes: Optional[Dict[str, List[str]]]  # 标记基因字典
reference_data: Optional[str]                # 参考数据路径
reference_data_id: Optional[str]             # 参考数据ID

# Tangram 特定参数
training_genes: Optional[List[str]]          # 训练基因
num_epochs: int = 500                        # 训练轮数
mode: Literal["cells", "clusters"] = "cells" # 映射模式

# scANVI 参数
scanvi_n_latent: int = 10                    # 潜在维度
scanvi_n_hidden: int = 128                   # 隐藏层大小
```

### 4. 空间变量基因识别

#### 4.1 find_spatial_variable_genes - 空间变量基因
**功能**: 识别具有空间表达模式的基因

**支持方法**:
```python
method: Literal["gaston", "spatialde", "spark"] = "gaston"
```

**GASTON 方法参数** (推荐):
```python
# 预处理参数
preprocessing_method: Literal["glmpca", "pearson_residuals"] = "glmpca"
n_components: int = 10                       # 降维组件数

# 分析参数
n_domains: int = 5                          # 空间域数量
num_bins: int = 70                          # 等深度分箱数
umi_threshold: int = 500                    # UMI 阈值

# 基因分类参数
continuous_quantile: float = 0.9            # 连续基因分位数
discontinuous_quantile: float = 0.9         # 不连续基因分位数
pvalue_threshold: float = 0.1               # P值阈值
```

**SpatialDE 方法参数**:
```python
spatialde_normalized: bool = True           # 数据是否已标准化
spatialde_kernel: str = "SE"                # 核函数类型
```

**SPARK 方法参数**:
```python
spark_percentage: float = 0.1               # 表达比例阈值
spark_min_total_counts: int = 10            # 最小总计数
spark_num_core: int = 1                     # 并行核心数
```

### 5. 空间分析

#### 5.1 analyze_spatial_patterns - 空间模式分析
**功能**: 执行各种空间统计分析

**分析类型**:
```python
analysis_type: Literal[
    "neighborhood",    # 邻域分析
    "co_occurrence",   # 共现分析
    "ripley",         # Ripley K 函数
    "moran",          # Moran's I 空间自相关
    "centrality",     # 中心性分析
    "getis_ord"       # Getis-Ord Gi* 统计
]
```

**参数**:
```python
cluster_key: str = "leiden"                 # 聚类键
n_neighbors: int = 15                       # 邻居数

# Getis-Ord 特定参数
getis_ord_genes: Optional[List[str]]        # 分析基因列表
getis_ord_n_genes: int = 20                 # 高变基因数量
getis_ord_correction: Literal["bonferroni", "fdr_bh", "none"] = "fdr_bh"
getis_ord_alpha: float = 0.05               # 显著性阈值
```

### 6. 空间域识别

#### 6.1 identify_spatial_domains - 空间域识别
**功能**: 识别组织中的空间功能域

**支持方法**:
```python
method: Literal["spagcn", "leiden", "louvain", "stagate", "banksy"] = "spagcn"
```

**SpaGCN 参数** (推荐):
```python
n_domains: int = 7                          # 空间域数量
spagcn_s: float = 1.0                       # 组织学权重
spagcn_b: int = 49                          # 斑点区域大小
spagcn_p: float = 0.5                       # 邻域贡献比例
spagcn_use_histology: bool = True           # 使用组织学图像
spagcn_random_seed: int = 100               # 随机种子
```

### 7. 细胞通讯分析

#### 7.1 analyze_cell_communication - 细胞通讯
**功能**: 分析细胞间配体-受体相互作用

**支持方法**:
```python
method: Literal["liana"] = "liana"          # 仅支持 LIANA+
```

**LIANA+ 参数**:
```python
species: Literal["human", "mouse", "zebrafish"] = "human"
min_cells: int = 3                          # 最少表达细胞数

# LIANA 特定参数
liana_resource: Literal["consensus", "cellchat", "cellphonedb", "connectome", "omnipath"] = "consensus"
liana_local_metric: Literal["cosine", "pearson", "spearman", "jaccard"] = "cosine"
liana_global_metric: Literal["morans", "lee"] = "morans"
liana_n_perms: int = 100                    # 置换检验次数
liana_nz_prop: float = 0.2                  # 最小表达比例
liana_bandwidth: Optional[int] = None        # 空间连接带宽
liana_cutoff: float = 0.1                   # 空间连接阈值
```

### 8. 空间反卷积

#### 8.1 perform_deconvolution - 空间反卷积
**功能**: 推断空间位点的细胞类型组成

**支持方法**:
```python
method: Literal[
    "cell2location",   # Cell2location (推荐)
    "rctd",           # RCTD
    "destvi",         # DestVI
    "stereoscope",    # Stereoscope
    "spotlight",      # SPOTlight
    "tangram",        # Tangram
    "mrvi"            # MRVI
]
```

**参数**:
```python
reference_data_id: Optional[str]            # 参考单细胞数据ID
cell_type_key: str = "cell_type"            # 细胞类型键
n_top_genes: int = 2000                     # 使用基因数
use_gpu: bool = False                       # 使用GPU
n_epochs: int = 10000                       # 训练轮数
n_cells_per_spot: Optional[int] = None      # 每个斑点细胞数

# Cell2location 特定参数
cell2location_n_cells_per_location: int = 30
cell2location_detection_alpha: float = 20.0
```

### 9. 轨迹分析

#### 9.1 analyze_trajectory - 轨迹分析
**功能**: 分析细胞状态转换和伪时间轨迹

**支持方法**:
```python
method: Literal["cellrank", "palantir", "velovi", "dpt"] = "cellrank"
```

**参数**:
```python
spatial_weight: float = 0.5                 # 空间信息权重
root_cells: Optional[List[str]] = None      # 根细胞

# CellRank 参数
cellrank_kernel_weights: Tuple[float, float] = (0.8, 0.2)
cellrank_n_states: int = 5                  # 宏状态数

# VeloVI 参数
velovi_n_hidden: int = 128                  # 隐藏层大小
velovi_n_latent: int = 10                   # 潜在维度
velovi_learning_rate: float = 1e-3          # 学习率
```

#### 9.2 analyze_rna_velocity - RNA速度分析
**功能**: 分析RNA速度和细胞动力学

**参数**:
```python
velocity_mode: Literal["deterministic", "stochastic"] = "deterministic"
n_top_genes: int = 2000                     # 速度基因数
min_shared_counts: int = 30                 # 最小共享计数
n_pcs: int = 30                            # 主成分数
```

### 10. 样本整合

#### 10.1 integrate_datasets - 多样本整合
**功能**: 整合多个空间转录组学样本

**支持方法**:
```python
method: Literal[
    "harmony",        # Harmony
    "bbknn",         # BBKNN
    "scanorama",     # Scanorama
    "mnn",           # MNN
    "scvi",          # scVI
    "multivi",       # MultiVI
    "totalvi"        # TotalVI
]
```

**参数**:
```python
batch_key: str = "batch"                    # 批次信息键
n_pcs: int = 30                            # 主成分数
align_spatial: bool = True                  # 对齐空间坐标
reference_batch: Optional[str] = None       # 参考批次

# scVI 参数
scvi_n_hidden: int = 128                   # 隐藏层大小
scvi_n_latent: int = 10                    # 潜在维度
scvi_gene_likelihood: Literal["zinb", "nb", "poisson"] = "zinb"
```

### 11. 差异表达分析

#### 11.1 perform_differential_expression - 差异表达
**功能**: 识别不同条件或细胞类型间的差异表达基因

**参数**:
```python
group_key: str                             # 分组键
reference_group: Optional[str] = None       # 参考组
test_method: Literal["wilcoxon", "t-test", "logreg"] = "wilcoxon"
n_genes: Optional[int] = None              # 返回基因数
min_fold_change: float = 1.5               # 最小倍数变化
```

### 12. 辅助功能

#### 12.1 list_datasets - 数据集列表
**功能**: 列出所有已加载的数据集

#### 12.2 register_spatial_data - 空间配准
**功能**: 配准多个空间切片

**参数**:
```python
source_id: str                             # 源数据集ID
target_id: str                             # 目标数据集ID
method: str = "paste"                      # 配准方法
landmarks: Optional[List[Dict[str, Any]]] = None  # 手动标记点
```

#### 12.3 calculate_spatial_statistics - 空间统计
**功能**: 计算特定特征的空间统计量

**参数**:
```python
feature: str                               # 分析特征
statistic: str = "morans_i"                # 统计类型
n_neighbors: int = 6                       # 邻居数
```

## 数据模型

### 输入参数模型
- `AnalysisParameters`: 预处理参数
- `VisualizationParameters`: 可视化参数
- `AnnotationParameters`: 注释参数
- `SpatialAnalysisParameters`: 空间分析参数
- `SpatialVariableGenesParameters`: 空间变量基因参数
- `CellCommunicationParameters`: 细胞通讯参数
- `DeconvolutionParameters`: 反卷积参数
- `TrajectoryParameters`: 轨迹分析参数
- `IntegrationParameters`: 整合参数
- `SpatialDomainParameters`: 空间域参数

### 输出结果模型
- `PreprocessingResult`: 预处理结果
- `AnnotationResult`: 注释结果
- `SpatialAnalysisResult`: 空间分析结果
- `SpatialVariableGenesResult`: 空间变量基因结果
- `CellCommunicationResult`: 细胞通讯结果
- `DeconvolutionResult`: 反卷积结果
- `TrajectoryResult`: 轨迹分析结果
- `IntegrationResult`: 整合结果
- `SpatialDomainResult`: 空间域结果

## 使用建议

### 推荐分析流程
1. **数据加载**: `load_data`
2. **数据预处理**: `preprocess_data`
3. **质量控制可视化**: `visualize_data`
4. **空间域识别**: `identify_spatial_domains` (SpaGCN)
5. **空间变量基因**: `find_spatial_variable_genes` (GASTON)
6. **细胞类型注释**: `annotate_cell_types` (Tangram/CellAssign)
7. **空间分析**: `analyze_spatial_patterns`
8. **细胞通讯**: `analyze_cell_communication` (LIANA+)
9. **结果可视化**: `visualize_data`

### 性能优化建议
- 大数据集使用GPU加速 (`use_gpu=True`)
- 适当调整基因数量参数 (`n_top_genes`)
- 使用高效的空间分析方法 (GASTON, SpaGCN)
- 合理设置邻居数和分辨率参数

### 常见问题
1. **超时问题**: Claude Desktop MCP 有4分钟超时限制
2. **内存问题**: 大数据集建议使用数据子采样
3. **依赖问题**: 确保所有必需的Python包已安装
4. **GPU支持**: 某些方法支持MPS加速 (Apple Silicon)

## 技术规格

### 系统要求
- **Python版本**: 3.10+ (推荐 3.11)
- **操作系统**: macOS, Linux, Windows
- **内存**: 最少 8GB RAM (推荐 16GB+)
- **存储**: 至少 5GB 可用空间

### 核心依赖
- **FastMCP**: MCP 协议实现框架
- **scanpy**: 单细胞分析核心库
- **squidpy**: 空间转录组学分析
- **scvi-tools**: 深度学习模型
- **liana**: 细胞通讯分析
- **pandas, numpy, scipy**: 数据处理
- **matplotlib, seaborn**: 可视化

### 可选依赖
- **GASTON**: 空间变量基因识别
- **SpaGCN**: 空间域识别
- **STAGATE**: 空间域识别 (备选)
- **BANKSY**: 空间域识别 (备选)
- **rpy2**: R 语言集成 (SPARK 方法)

### GPU 支持
- **CUDA**: NVIDIA GPU 加速
- **MPS**: Apple Silicon GPU 加速
- **支持的方法**: cell2location, scVI, DestVI, VeloVI

### 数据格式支持
- **输入格式**:
  - AnnData (h5ad) - 推荐
  - 10x Genomics (h5, mtx)
  - CSV/TSV 文件
  - Zarr 格式
- **输出格式**:
  - AnnData (h5ad)
  - PNG 图像
  - JSON 结果

### MCP 协议兼容性
- **MCP 版本**: 1.0+
- **传输协议**: stdio, SSE
- **消息格式**: JSON-RPC 2.0
- **资源类型**: 数据集, 可视化, 分析结果
- **提示支持**: 预定义工作流模板

## 架构图

### 整体架构
```
┌─────────────────────────────────────────────────────────────┐
│                    Claude Desktop (MCP Client)              │
└─────────────────────┬───────────────────────────────────────┘
                      │ MCP Protocol (JSON-RPC 2.0)
┌─────────────────────▼───────────────────────────────────────┐
│                  ChatSpatial MCP Server                     │
│  ┌─────────────────────────────────────────────────────┐   │
│  │              SpatialMCPAdapter                      │   │
│  │  ┌─────────────┬─────────────┬─────────────────┐   │   │
│  │  │   Tools     │  Resources  │    Prompts      │   │   │
│  │  │ Management  │ Management  │   Management    │   │   │
│  │  └─────────────┴─────────────┴─────────────────┘   │   │
│  └─────────────────────────────────────────────────────┘   │
│  ┌─────────────────────────────────────────────────────┐   │
│  │              Error Handling System                  │   │
│  │  ┌─────────────┬─────────────┬─────────────────┐   │   │
│  │  │ Tool Error  │  Pydantic   │   Parameter     │   │   │
│  │  │  Handler    │   Handler   │   Validation    │   │   │
│  │  └─────────────┴─────────────┴─────────────────┘   │   │
│  └─────────────────────────────────────────────────────┘   │
│  ┌─────────────────────────────────────────────────────┐   │
│  │               Data Management                       │   │
│  │  ┌─────────────┬─────────────┬─────────────────┐   │   │
│  │  │   Dataset   │   Results   │     Cache       │   │   │
│  │  │   Storage   │   Storage   │   Management    │   │   │
│  │  └─────────────┴─────────────┴─────────────────┘   │   │
│  └─────────────────────────────────────────────────────┘   │
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│              Spatial Analysis Tools                         │
│  ┌─────────┬─────────┬─────────┬─────────┬─────────────┐   │
│  │ Preproc │  Annot  │ Spatial │  Deconv │    Viz      │   │
│  │ essing  │  ation  │ Analysis│  olution│ ualization  │   │
│  └─────────┴─────────┴─────────┴─────────┴─────────────┘   │
│  ┌─────────┬─────────┬─────────┬─────────┬─────────────┐   │
│  │ Domains │  Genes  │  Comm   │ Traject │ Integration │   │
│  │   ID    │  Ident  │ Analysis│  ory    │             │   │
│  └─────────┴─────────┴─────────┴─────────┴─────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

### 数据流图
```
┌─────────────┐    ┌─────────────┐    ┌─────────────┐
│   Load      │───▶│ Preprocess  │───▶│  Analysis   │
│   Data      │    │    Data     │    │   Tools     │
└─────────────┘    └─────────────┘    └─────────────┘
       │                   │                   │
       ▼                   ▼                   ▼
┌─────────────┐    ┌─────────────┐    ┌─────────────┐
│   Data      │    │   QC        │    │   Results   │
│  Manager    │    │ Metrics     │    │   Storage   │
└─────────────┘    └─────────────┘    └─────────────┘
       │                   │                   │
       ▼                   ▼                   ▼
┌─────────────┐    ┌─────────────┐    ┌─────────────┐
│ Resource    │    │ Validation  │    │Visualization│
│ Creation    │    │   System    │    │   Engine    │
└─────────────┘    └─────────────┘    └─────────────┘
```

## 开发和扩展

### 添加新工具
1. 在 `tools/` 目录创建新模块
2. 实现分析函数和参数模型
3. 在 `server.py` 中注册工具
4. 添加工具元数据和注解
5. 编写测试和文档

### 自定义数据管理器
```python
class CustomDataManager(SpatialDataManager):
    async def load_dataset(self, path: str, data_type: str, name: str) -> str:
        # 自定义数据加载逻辑
        pass

    async def save_result(self, data_id: str, result_type: str, result: Any) -> None:
        # 自定义结果存储逻辑
        pass
```

### 扩展错误处理
```python
@mcp_tool_error_handler(include_traceback=True)
@mcp_pydantic_error_handler()
async def custom_analysis_tool(params: CustomParameters):
    # 工具实现
    pass
```

---

*本文档基于ChatSpatial MCP v1.0，最后更新：2025年1月*
