# ChatSpatial 10x H5 格式直接支持增强方案

## 背景分析

### 当前问题
ChatSpatial 当前对 10x Genomics H5 文件支持存在以下限制：

1. **不支持单独 H5 文件加载**：需要完整的 Visium 目录结构（包含 spatial 文件夹）
2. **错误信息不友好**：返回 "Unsupported file format for 10x_visium" 
3. **用户体验差**：需要手动组织文件结构才能加载数据

### 用户需求场景
- 用户下载了 GEO 数据集，只有 `.h5` 文件和单独的 `spatial` 文件夹
- 用户有多个 `.h5` 文件但共享同一个空间坐标信息
- 用户想快速测试基因表达分析，暂时不需要空间信息

## 技术分析

### 现有代码结构

#### 1. 数据加载流程
```
server.py (MCP endpoint)
    ↓
spatial_mcp_adapter.py (load_dataset)
    ↓
utils/data_loader.py (load_spatial_data)
    ↓
scanpy functions (sc.read_visium, sc.read_h5ad)
```

#### 2. 关键代码位置
- **主入口**: `/chatspatial/server.py:96-132` - `load_data` MCP 工具
- **数据管理**: `/chatspatial/spatial_mcp_adapter.py:547-561` - `load_dataset` 方法
- **核心加载逻辑**: `/chatspatial/utils/data_loader.py:15-169` - `load_spatial_data` 函数

#### 3. 当前支持的格式
- `10x_visium`: 完整目录结构
- `h5ad`: AnnData 格式
- `other`: 通用格式（实际调用 `sc.read_h5ad`）

### Scanpy 能力分析
Scanpy 提供了两个关键函数：
- `sc.read_10x_h5()`: 读取 10x H5 文件，返回 AnnData 对象
- `sc.read_visium()`: 读取完整 Visium 数据（包括空间信息）

## 解决方案设计

### 方案一：扩展现有 10x_visium 类型（推荐）

#### 优点
- 对现有代码改动最小
- 向后兼容性好
- 用户接口不变

#### 实现细节

1. **修改 `data_loader.py` 的 `load_spatial_data` 函数**

```python
# 在第 61-132 行的 10x_visium 处理逻辑中添加

if data_type == "10x_visium":
    try:
        # 情况1: 标准 Visium 目录结构
        if os.path.isdir(data_path):
            # ... 现有目录处理逻辑 ...
            
        # 情况2: 单独的 H5 文件（新增）
        elif os.path.isfile(data_path) and data_path.endswith('.h5'):
            # 读取 10x H5 文件
            adata = sc.read_10x_h5(data_path)
            
            # 尝试查找空间信息
            spatial_path = _find_spatial_folder(data_path)
            if spatial_path:
                adata = _add_spatial_info(adata, spatial_path)
                logger.info(f"Added spatial information from {spatial_path}")
            else:
                logger.warning("No spatial folder found. Loading expression data only.")
            
        # 情况3: H5AD 文件（保持现有逻辑）
        elif os.path.isfile(data_path) and data_path.endswith('.h5ad'):
            # ... 现有 h5ad 处理逻辑 ...
```

2. **添加辅助函数**

```python
def _find_spatial_folder(h5_path: str) -> Optional[str]:
    """
    智能查找空间信息文件夹
    
    搜索策略：
    1. 同级目录下的 'spatial' 文件夹
    2. 父目录下的 'spatial' 文件夹
    3. 同名前缀的 spatial 文件夹
    """
    base_dir = os.path.dirname(h5_path)
    base_name = os.path.splitext(os.path.basename(h5_path))[0]
    
    # 候选路径
    candidates = [
        os.path.join(base_dir, 'spatial'),
        os.path.join(base_dir, '..', 'spatial'),
        os.path.join(base_dir, f'{base_name}_spatial'),
        os.path.join(base_dir, 'spatial_data'),
    ]
    
    for candidate in candidates:
        if os.path.exists(candidate) and os.path.isdir(candidate):
            # 验证是否包含必要文件
            required_files = ['tissue_positions_list.csv', 'scalefactors_json.json']
            if all(os.path.exists(os.path.join(candidate, f)) for f in required_files):
                return candidate
    
    return None


def _add_spatial_info(adata: AnnData, spatial_path: str) -> AnnData:
    """
    将空间信息添加到 AnnData 对象
    
    参数：
        adata: 表达数据
        spatial_path: 空间信息文件夹路径
    """
    import json
    import pandas as pd
    import numpy as np
    
    # 加载 tissue positions
    positions_file = os.path.join(spatial_path, 'tissue_positions_list.csv')
    positions = pd.read_csv(positions_file, header=None)
    positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 
                        'pxl_row_in_fullres', 'pxl_col_in_fullres']
    positions.set_index('barcode', inplace=True)
    
    # 只保留数据中存在的 barcodes
    common_barcodes = adata.obs_names.intersection(positions.index)
    
    if len(common_barcodes) == 0:
        raise ValueError("No matching barcodes between expression data and spatial coordinates")
    
    # 过滤到共同的 barcodes
    adata = adata[common_barcodes, :].copy()
    positions = positions.loc[common_barcodes]
    
    # 添加空间坐标
    adata.obsm['spatial'] = positions[['pxl_col_in_fullres', 'pxl_row_in_fullres']].values
    
    # 添加组织信息
    adata.obs['in_tissue'] = positions['in_tissue'].values
    
    # 加载 scalefactors
    scalefactors_file = os.path.join(spatial_path, 'scalefactors_json.json')
    with open(scalefactors_file, 'r') as f:
        scalefactors = json.load(f)
    
    # 创建 spatial uns 结构
    adata.uns['spatial'] = {
        'library_id': 'spatial',
        'scalefactors': scalefactors
    }
    
    # 尝试加载图像
    for img_name in ['tissue_hires_image.png', 'tissue_lowres_image.png']:
        img_path = os.path.join(spatial_path, img_name)
        if os.path.exists(img_path):
            from PIL import Image
            img = np.array(Image.open(img_path))
            
            if 'images' not in adata.uns['spatial']:
                adata.uns['spatial']['images'] = {}
            
            img_key = 'hires' if 'hires' in img_name else 'lowres'
            adata.uns['spatial']['images'][img_key] = img
    
    return adata
```

3. **改进错误处理**

```python
# 在 data_loader.py 第 132 行附近
except FileNotFoundError as e:
    raise ValueError(f"File not found: {str(e)}")
except Exception as e:
    # 提供更详细的错误信息
    error_msg = f"Error loading 10x Visium data from {data_path}: {str(e)}"
    
    # 添加帮助信息
    if "No matching barcodes" in str(e):
        error_msg += "\nPlease check if the H5 file and spatial coordinates are from the same sample."
    elif ".h5" in data_path:
        error_msg += "\nTry providing the spatial folder path using the spatial_path parameter."
    
    raise ValueError(error_msg)
```

### 方案二：添加新的数据类型 "10x_h5"

#### 优点
- 更清晰的语义
- 可以有专门的参数

#### 缺点
- 需要修改更多地方
- 增加用户学习成本

#### 实现要点
1. 在 `load_spatial_data` 的 `data_type` 参数中添加 "10x_h5" 选项
2. 添加独立的处理分支
3. 更新文档和类型提示

### 方案三：扩展 load_data 参数

#### 建议添加的参数

```python
async def load_data(
    data_path: str,
    data_type: str = "auto",
    name: Optional[str] = None,
    spatial_path: Optional[str] = None,  # 新增：独立指定空间信息路径
    context: Context = None
) -> SpatialDataset:
```

这样用户可以：
```python
# 分别指定表达数据和空间数据
load_data(
    data_path="GSE198353_spleen_rep_1_filtered_feature_bc_matrix.h5",
    spatial_path="spatial/",
    data_type="10x_visium"
)
```

## 实施计划

### 第一阶段：核心功能实现
1. 修改 `data_loader.py` 支持单独 H5 文件
2. 添加智能空间信息查找
3. 改进错误处理和提示

### 第二阶段：增强功能
1. 添加 `spatial_path` 参数支持
2. 支持批量加载多个 H5 文件共享空间信息
3. 添加数据验证和自动修复功能

### 第三阶段：用户体验优化
1. 添加进度条和详细日志
2. 提供数据预览功能
3. 自动检测和建议最佳加载方式

## 测试计划

### 单元测试
```python
def test_load_10x_h5_only():
    """测试只有 H5 文件的加载"""
    adata = load_spatial_data("test.h5", data_type="10x_visium")
    assert adata is not None
    assert 'spatial' not in adata.obsm  # 没有空间信息

def test_load_10x_h5_with_spatial():
    """测试 H5 + spatial 文件夹的加载"""
    adata = load_spatial_data("test.h5", data_type="10x_visium")
    assert 'spatial' in adata.obsm
    assert 'spatial' in adata.uns

def test_barcode_mismatch():
    """测试 barcode 不匹配的错误处理"""
    with pytest.raises(ValueError, match="No matching barcodes"):
        load_spatial_data("wrong.h5", data_type="10x_visium")
```

### 集成测试
使用实际的 GEO 数据集进行测试：
- GSE198353 (脾脏数据)
- 其他常见 10x 数据集

## 向后兼容性

- ✅ 现有的目录结构加载方式不受影响
- ✅ 现有的 h5ad 加载方式不受影响
- ✅ API 接口保持不变（只是扩展功能）

## 文档更新

需要更新的文档：
1. `/docs/reference/data_formats.md` - 添加 H5 格式说明
2. `/docs/getting-started/quickstart.md` - 添加新的使用示例
3. `/docs/reference/troubleshooting/common_issues.md` - 添加常见问题解答

## 预期效果

实施后，用户可以：
```python
# 直接加载 H5 文件
mcp__chatspatial__load_data(
    data_path="GSE198353_spleen_rep_1_filtered_feature_bc_matrix.h5",
    data_type="10x_visium"
)

# 自动查找并加载空间信息
# 或手动指定空间信息路径
mcp__chatspatial__load_data(
    data_path="expression.h5",
    spatial_path="my_spatial_folder/",
    data_type="10x_visium"
)
```

## 风险评估

1. **性能影响**：最小，只是增加了文件格式判断
2. **兼容性风险**：低，使用 scanpy 标准函数
3. **维护成本**：低，代码结构清晰

## 结论

推荐采用**方案一 + 方案三的组合**：
1. 扩展现有 `10x_visium` 类型支持 H5 文件
2. 添加 `spatial_path` 可选参数
3. 实现智能空间信息查找

这样既保持了向后兼容性，又提供了灵活性，同时代码改动最小。