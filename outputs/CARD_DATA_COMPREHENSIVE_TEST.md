# CARD 数据综合功能测试报告

**测试日期**: 2025-10-16
**数据集**: CARD example (胰腺组织)
- Spatial: 428 spots × 25753 genes
- Reference: 1926 cells × 19736 genes (20 cell types)

**测试目的**:
1. 验证所有功能能否正常运行
2. 记录代码运行的 bug
3. 记录运行结果的科学性/统计性存疑问题

---

## 测试进度

- [ ] 1. Preprocessing 预处理功能
- [ ] 2. Spatial visualization 空间可视化
- [ ] 3. Cell annotation 细胞注释（多种方法）
- [ ] 4. Spatial statistics 空间统计分析
- [ ] 5. Enrichment analysis 富集分析
- [ ] 6. Spatial variable genes 空间变异基因
- [ ] 7. Cell communication 细胞通讯
- [ ] 8. Trajectory analysis 轨迹分析
- [ ] 9. Velocity analysis 速度分析
- [ ] 10. 其他功能

---

## 测试结果详情

### 1. ✅ Preprocessing 预处理功能

**状态**: 通过
**参数**: log normalization, 2000 HVGs, 15 neighbors, 30 PCs
**结果**:
- 过滤后: 428 spots × 16,056 genes (从 25,753)
- 检测到 7 个 Leiden clusters
- QC 指标: 中位数 970 genes/spot, 1569 UMIs/spot

**科学性评估**: ✅ 正常

---

### 2. ✅ Spatial Visualization 空间可视化

**状态**: 通过
**测试项**:
- Spatial clustering (leiden)
- UMAP 降维
- 基因表达空间分布 (REG3A)

**结果**: 所有可视化正常，空间模式清晰

**科学性评估**: ✅ 正常

---

### 3. Cell Annotation 细胞注释

#### 3.1 ❌ scANVI 方法

**Bug #1**: 默认参数导致 NaN 错误
- 错误信息: `Expected parameter loc...to satisfy the constraint Real(), but found invalid values: tensor([[nan, nan, ...]])`
- 原因: 数据规模较小 (428 spots)，模型过度参数化
- 解决方案: 使用小数据集参数 (n_latent=5, n_hidden=64, dropout=0.2, epochs=50, no pretrain)

**⚠️ 科学性问题 #1**: scANVI 细胞类型分布异常
- **异常高比例**:
  - Tuft_cells: 108 spots (25.2%)
  - RBCs: 81 spots (18.9%)
  - 合计 44.1% 只分配给2种细胞类型

- **与 CARD 结果矛盾**:
  | Cell Type | scANVI | CARD |
  |-----------|--------|------|
  | Tuft_cells | 25.2% | 0.75% |
  | RBCs | 18.9% | 0.03% |
  | Ductal_CRISP3 | 5.4% | 27.9% |
  | Ductal_terminal | 6.3% | 13.4% |

- **置信度偏低**: 范围 0.25-0.67 (期望 >0.8)

**结论**: scANVI 在此数据集上结果可疑，可能不适合小数据集或跨样本场景

#### 3.2 ❌ Tangram Clusters Mode 方法

**Bug #2**: 返回空的 counts 和 confidence_scores
- `counts`: {} (空字典)
- `confidence_scores`: {} (空字典)
- cell_types 列表包含20个类型名称，但未实际映射

**Bug #3**: Tangram clusters mode 完全失败
- 所有428个 spots 显示为 "NA"
- 没有任何 spot 被分配到细胞类型
- 可视化确认：所有点为灰色 (NA)

**参数**:
- mode="clusters"
- cluster_label="cellType"
- num_epochs=500

**结论**: Tangram clusters mode 在此数据集上完全失效，需要深度调试

---

### 4. Spatial Statistics 空间统计分析

#### 4.1 ✅ Neighborhood Enrichment

**状态**: 通过
**结果**: max_enrichment=21.89, min_enrichment=-12.89
**说明**: 之前修复的 bug (np.nanmax/nanmin) 生效

#### 4.2 ⚠️ Moran's I 空间自相关

**科学性问题 #2**: 分析基因数量不符
- 请求: n_top_genes=100
- 实际: n_genes_analyzed=20
- Mean Moran's I = -0.002 (接近0，空间自相关极弱)

**可能原因**: 参数未正确传递或内部限制

#### 4.3 ✅ Ripley's K

**状态**: 通过
**结果**: analysis_completed=true

---

### 5. ❌ Enrichment Analysis 富集分析

**Bug #4**: 富集分析零结果问题未完全解决

**测试结果**:
| 数据库 | n_significant | 状态 |
|--------|---------------|------|
| KEGG_Pathways | 0 | ❌ |
| GO_Biological_Process | 0 | ❌ |
| Reactome_Pathways | 0 | ❌ |

**问题**:
- 所有数据库都返回 n_significant=0
- 之前修复的 recarray 索引 bug 可能没有完全解决
- find_markers 正常返回100个基因，但富集分析无法使用

**科学性问题 #3**: 胰腺数据应该能检测到消化/代谢通路
- Marker genes 包含: FBLN1, SERPINF1, PTGDS, IGFBP4 等
- 这些基因应该富集在 ECM-receptor interaction, PI3K-Akt signaling 等通路

**可能原因**:
1. rank_genes_groups 结果未正确存储/读取
2. 基因 ID 格式不匹配（大小写、点号 vs 连字符）
3. ORA 分析的查询基因列表为空

---

### 6. ✅ Spatial Variable Genes 空间变异基因

#### 6.1 SPARK-X 方法

**状态**: 通过
**参数**: sparkx_num_core=1, n_top_genes=100
**结果**:
- n_genes_analyzed: 100
- n_significant_genes: 100 (100%)
- Top genes (p < 1e-17):
  - KRT17 (keratin 17)
  - BGN (biglycan, ECM protein)
  - SERPINF1 (serpin family F member 1)
  - IGFBP4 (insulin-like growth factor binding protein 4)
  - COL6A3, COL3A1, COL4A2 (collagens)
- Pancreatic markers detected:
  - CPA2, PRSS2, REG3A, CTRB2 (digestive enzymes)

**科学性评估**: ✅ 正常
- 检测到胰腺特异性基因
- 空间模式显著 (p值极小)
- 基因功能符合组织特征

---

## 总结

### Bug 清单

| # | 模块 | Bug 描述 | 严重性 |
|---|------|---------|--------|
| 1 | scANVI | 默认参数导致 NaN 错误（小数据集） | 中 |
| 2 | Tangram | counts 和 confidence_scores 返回空字典 | 高 |
| 3 | Tangram | Clusters mode 完全失败（所有spots=NA） | 高 |
| 4 | Enrichment | 所有数据库零结果（recarray bug未完全修复） | 高 |

### 科学性/统计性问题清单

| # | 模块 | 问题描述 | 严重性 |
|---|------|---------|--------|
| 1 | scANVI | 细胞类型分布异常（Tuft_cells 25%, RBCs 19%）| 高 |
| 2 | Moran's I | 分析基因数量不符（请求100，实际20） | 中 |
| 3 | Enrichment | 胰腺数据无法检测消化/代谢通路 | 高 |

### 功能状态汇总

| 功能模块 | 状态 | 备注 |
|---------|------|------|
| Preprocessing | ✅ 正常 | |
| Visualization | ✅ 正常 | |
| CARD Deconvolution | ✅ 正常 | 多参数测试稳定 |
| scANVI Annotation | ⚠️ 可用 | 需小数据集参数，结果可疑 |
| Tangram Annotation | ❌ 失败 | Clusters mode 完全失效 |
| Neighborhood Enrichment | ✅ 正常 | 已修复 |
| Moran's I | ⚠️ 部分正常 | 基因数量问题 |
| Ripley's K | ✅ 正常 | |
| Enrichment Analysis | ❌ 失败 | 所有数据库零结果 |
| SPARK-X | ✅ 正常 | 检测到100个显著基因 |

### 关键发现

1. **CARD deconvolution 最稳定**：多参数测试均正常，除严格过滤参数外
2. **Tangram clusters mode 需修复**：这是官方推荐的跨样本模式，但完全失效
3. **Enrichment分析核心bug**：之前的修复不完整，需要更深入调试
4. **scANVI 结果可疑**：细胞类型分布与CARD相差极大，置信度偏低

---

