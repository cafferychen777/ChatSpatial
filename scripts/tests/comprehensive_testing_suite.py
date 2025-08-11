#!/usr/bin/env python3
"""
超全面测试套件 - 深度测试可视化增强功能
包括边界情况、压力测试、真实场景、错误处理、性能测试等
"""

import asyncio
import sys
import os
import time
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Add the project root to the path
sys.path.insert(0, str(Path(__file__).parent))

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import visualize_data
from chatspatial.utils.image_utils import fig_to_image


class ComprehensiveTestSuite:
    """超全面测试套件类"""
    
    def __init__(self):
        self.test_results = {}
        self.performance_results = {}
        self.data_variants = {}
        
    def create_various_test_datasets(self):
        """创建各种类型的测试数据集"""
        print("🔄 创建多种测试数据集...")
        
        datasets = {}
        
        # 1. 小数据集 (边界情况)
        datasets['tiny'] = self._create_tiny_dataset()
        
        # 2. 中等数据集 (标准场景)
        datasets['medium'] = self._create_medium_dataset()
        
        # 3. 大数据集 (性能测试)
        datasets['large'] = self._create_large_dataset()
        
        # 4. 稀疏数据集 (特殊情况)
        datasets['sparse'] = self._create_sparse_dataset()
        
        # 5. 不平衡数据集 (真实场景)
        datasets['imbalanced'] = self._create_imbalanced_dataset()
        
        # 6. 有缺失值的数据集
        datasets['missing_data'] = self._create_missing_data_dataset()
        
        # 7. 10x Visium风格数据集
        datasets['visium'] = self._create_visium_style_dataset()
        
        # 8. 极端值数据集
        datasets['extreme_values'] = self._create_extreme_values_dataset()
        
        self.data_variants = datasets
        print(f"✅ 创建了 {len(datasets)} 个测试数据集")
        return datasets
    
    def _create_tiny_dataset(self):
        """创建极小数据集 (边界测试)"""
        n_cells, n_genes = 10, 20
        np.random.seed(42)
        
        adata = sc.AnnData(
            X=np.random.poisson(2, (n_cells, n_genes)).astype(np.float32),
            obs=pd.DataFrame(index=[f"Cell_{i}" for i in range(n_cells)]),
            var=pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])
        )
        
        # 极简空间坐标
        adata.obsm['spatial'] = np.random.normal(0, 1, (n_cells, 2))
        adata.obs['leiden'] = pd.Categorical(['A', 'B'] * 5)
        adata.obs['batch'] = pd.Categorical(['Batch1'] * 10)
        
        return adata
    
    def _create_medium_dataset(self):
        """创建中等数据集 (标准测试)"""
        n_cells, n_genes = 500, 1000
        np.random.seed(123)
        
        adata = sc.AnnData(
            X=np.random.poisson(5, (n_cells, n_genes)).astype(np.float32),
            obs=pd.DataFrame(index=[f"Cell_{i:04d}" for i in range(n_cells)]),
            var=pd.DataFrame(index=[f"Gene_{i:04d}" for i in range(n_genes)])
        )
        
        # 空间坐标 - 形成明显的cluster
        centers = np.array([[0, 0], [10, 0], [5, 8], [-5, 5]])
        cluster_labels = []
        coords = []
        
        for i in range(n_cells):
            cluster_idx = i % 4
            center = centers[cluster_idx]
            coord = center + np.random.normal(0, 2, 2)
            coords.append(coord)
            cluster_labels.append(f"Cluster_{cluster_idx}")
        
        adata.obsm['spatial'] = np.array(coords)
        adata.obs['leiden'] = pd.Categorical(cluster_labels)
        adata.obs['batch'] = pd.Categorical(np.random.choice(['Batch1', 'Batch2'], n_cells))
        adata.obs['cell_type'] = pd.Categorical(np.random.choice(['TypeA', 'TypeB', 'TypeC'], n_cells))
        
        # 添加高变基因
        sc.pp.highly_variable_genes(adata, n_top_genes=100)
        
        # 计算邻居和UMAP
        sc.pp.neighbors(adata, n_neighbors=10)
        sc.tl.umap(adata)
        
        # 添加邻域富集结果
        n_clusters = len(adata.obs['leiden'].cat.categories)
        enrichment_matrix = np.random.normal(0, 1.5, (n_clusters, n_clusters))
        enrichment_matrix = (enrichment_matrix + enrichment_matrix.T) / 2
        np.fill_diagonal(enrichment_matrix, np.abs(np.diagonal(enrichment_matrix)) + 2)
        
        adata.uns['leiden_nhood_enrichment'] = {'zscore': enrichment_matrix}
        
        return adata
    
    def _create_large_dataset(self):
        """创建大数据集 (性能测试)"""
        n_cells, n_genes = 5000, 3000
        np.random.seed(456)
        
        # 使用更高效的数据生成
        from scipy.sparse import csr_matrix
        
        # 生成稀疏矩阵以节省内存
        density = 0.1  # 10%的非零值
        X_data = np.random.poisson(3, int(n_cells * n_genes * density))
        row_indices = np.random.randint(0, n_cells, len(X_data))
        col_indices = np.random.randint(0, n_genes, len(X_data))
        X_sparse = csr_matrix((X_data, (row_indices, col_indices)), shape=(n_cells, n_genes))
        
        adata = sc.AnnData(
            X=X_sparse.astype(np.float32),
            obs=pd.DataFrame(index=[f"Cell_{i:05d}" for i in range(n_cells)]),
            var=pd.DataFrame(index=[f"Gene_{i:04d}" for i in range(n_genes)])
        )
        
        # 大规模空间坐标
        adata.obsm['spatial'] = np.random.normal(0, 20, (n_cells, 2))
        
        # 更多cluster
        n_clusters = 10
        cluster_labels = [f"Cluster_{i%n_clusters}" for i in range(n_cells)]
        adata.obs['leiden'] = pd.Categorical(cluster_labels)
        adata.obs['batch'] = pd.Categorical(np.random.choice(['B1', 'B2', 'B3', 'B4'], n_cells))
        
        return adata
    
    def _create_sparse_dataset(self):
        """创建稀疏数据集"""
        n_cells, n_genes = 200, 500
        np.random.seed(789)
        
        # 极度稀疏的数据 (90%的零值)
        X = np.random.poisson(0.5, (n_cells, n_genes)).astype(np.float32)
        X[X > 3] = 0  # 进一步增加稀疏性
        
        adata = sc.AnnData(
            X=X,
            obs=pd.DataFrame(index=[f"Cell_{i}" for i in range(n_cells)]),
            var=pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])
        )
        
        adata.obsm['spatial'] = np.random.normal(0, 5, (n_cells, 2))
        adata.obs['leiden'] = pd.Categorical(['Sparse_A', 'Sparse_B', 'Sparse_C'] * (n_cells//3 + 1))[:n_cells]
        adata.obs['batch'] = pd.Categorical(['SparseBatch'] * n_cells)
        
        return adata
    
    def _create_imbalanced_dataset(self):
        """创建不平衡数据集 (模拟真实场景)"""
        n_cells, n_genes = 800, 1200
        np.random.seed(321)
        
        adata = sc.AnnData(
            X=np.random.poisson(4, (n_cells, n_genes)).astype(np.float32),
            obs=pd.DataFrame(index=[f"Cell_{i:04d}" for i in range(n_cells)]),
            var=pd.DataFrame(index=[f"Gene_{i:04d}" for i in range(n_genes)])
        )
        
        adata.obsm['spatial'] = np.random.normal(0, 15, (n_cells, 2))
        
        # 极不平衡的cluster分布 (模拟稀有细胞类型)
        cluster_sizes = [500, 200, 50, 30, 20]  # 严重不平衡
        cluster_labels = []
        for i, size in enumerate(cluster_sizes):
            cluster_labels.extend([f"Type_{i}"] * size)
        
        adata.obs['leiden'] = pd.Categorical(cluster_labels[:n_cells])
        adata.obs['batch'] = pd.Categorical(['Imbalanced_B1', 'Imbalanced_B2'] * (n_cells//2))
        
        return adata
    
    def _create_missing_data_dataset(self):
        """创建有缺失数据的数据集"""
        n_cells, n_genes = 300, 600
        np.random.seed(654)
        
        adata = sc.AnnData(
            X=np.random.poisson(3, (n_cells, n_genes)).astype(np.float32),
            obs=pd.DataFrame(index=[f"Cell_{i}" for i in range(n_cells)]),
            var=pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])
        )
        
        adata.obsm['spatial'] = np.random.normal(0, 8, (n_cells, 2))
        
        # 部分缺失的metadata
        leiden_labels = ['Missing_A', 'Missing_B', 'Missing_C'] * (n_cells//3)
        leiden_labels = leiden_labels[:n_cells]
        # 随机设置一些为NaN
        leiden_labels = pd.Categorical(leiden_labels)
        
        adata.obs['leiden'] = leiden_labels
        adata.obs['batch'] = pd.Categorical(['MissingBatch'] * n_cells)
        
        # 故意不添加UMAP和neighbors (测试缺失数据处理)
        
        return adata
    
    def _create_visium_style_dataset(self):
        """创建10x Visium风格的数据集"""
        n_cells, n_genes = 400, 800
        np.random.seed(987)
        
        adata = sc.AnnData(
            X=np.random.negative_binomial(10, 0.3, (n_cells, n_genes)).astype(np.float32),
            obs=pd.DataFrame(index=[f"AAACAAGTATCTCCCA-1_{i}" for i in range(n_cells)]),
            var=pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])
        )
        
        # Visium式的六边形网格坐标
        coords = []
        for i in range(n_cells):
            row = i // 20
            col = i % 20
            x = col + (row % 2) * 0.5  # 六边形偏移
            y = row * 0.866  # 六边形间距
            coords.append([x * 100, y * 100])
        
        adata.obsm['spatial'] = np.array(coords)
        
        # 添加Visium风格的metadata
        adata.obs['leiden'] = pd.Categorical([f"Visium_{i%6}" for i in range(n_cells)])
        adata.obs['batch'] = pd.Categorical(['V1'] * n_cells)
        
        # 模拟组织图像信息
        adata.uns['spatial'] = {
            'V1': {
                'images': {'hires': np.random.rand(100, 100, 3)},
                'scalefactors': {'tissue_hires_scalef': 0.1}
            }
        }
        
        return adata
    
    def _create_extreme_values_dataset(self):
        """创建包含极端值的数据集"""
        n_cells, n_genes = 150, 300
        np.random.seed(111)
        
        X = np.random.poisson(2, (n_cells, n_genes)).astype(np.float32)
        
        # 添加一些极端值
        X[0, 0] = 1000000  # 极大值
        X[1, 1] = 0        # 极小值
        X[np.random.choice(n_cells, 10), np.random.choice(n_genes, 10)] = np.inf  # 无穷大
        X[np.random.choice(n_cells, 5), np.random.choice(n_genes, 5)] = -np.inf   # 负无穷大
        X[np.random.choice(n_cells, 3), np.random.choice(n_genes, 3)] = np.nan    # NaN值
        
        adata = sc.AnnData(
            X=X,
            obs=pd.DataFrame(index=[f"Extreme_{i}" for i in range(n_cells)]),
            var=pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])
        )
        
        adata.obsm['spatial'] = np.random.normal(0, 10, (n_cells, 2))
        adata.obs['leiden'] = pd.Categorical(['Extreme_A', 'Extreme_B'] * (n_cells//2 + 1))[:n_cells]
        adata.obs['batch'] = pd.Categorical(['ExtremeBatch'] * n_cells)
        
        return adata


    async def run_edge_case_tests(self):
        """运行边界情况测试"""
        print("\n🧪 运行边界情况测试...")
        
        edge_results = {}
        
        # 测试1: 空参数
        edge_results['empty_params'] = await self._test_empty_parameters()
        
        # 测试2: 无效plot_type
        edge_results['invalid_plot_type'] = await self._test_invalid_plot_type()
        
        # 测试3: 不存在的特征
        edge_results['nonexistent_feature'] = await self._test_nonexistent_feature()
        
        # 测试4: 极小数据集
        edge_results['tiny_dataset'] = await self._test_tiny_dataset()
        
        # 测试5: 缺失关键数据
        edge_results['missing_key_data'] = await self._test_missing_key_data()
        
        # 测试6: 极端参数值
        edge_results['extreme_parameters'] = await self._test_extreme_parameters()
        
        # 测试7: 内存限制
        edge_results['memory_limit'] = await self._test_memory_constraints()
        
        return edge_results
    
    async def _test_empty_parameters(self):
        """测试空参数处理"""
        try:
            data_store = {"test": {"adata": self.data_variants['medium']}}
            params = VisualizationParameters()  # 完全默认参数
            result = await visualize_data("test", data_store, params)
            return True
        except Exception as e:
            print(f"❌ 空参数测试失败: {str(e)}")
            return False
    
    async def _test_invalid_plot_type(self):
        """测试无效plot_type处理"""
        try:
            data_store = {"test": {"adata": self.data_variants['medium']}}
            params = VisualizationParameters(plot_type="invalid_type_xyz")
            result = await visualize_data("test", data_store, params)
            return False  # 应该抛出异常
        except Exception as e:
            # 预期的异常
            return "Invalid plot_type" in str(e) or "invalid_type_xyz" in str(e)
    
    async def _test_nonexistent_feature(self):
        """测试不存在的特征处理"""
        try:
            data_store = {"test": {"adata": self.data_variants['medium']}}
            params = VisualizationParameters(
                plot_type="spatial",
                feature="NonExistentGene_XYZ123"
            )
            result = await visualize_data("test", data_store, params)
            return True  # 应该优雅处理
        except Exception as e:
            print(f"❌ 不存在特征测试失败: {str(e)}")
            return False
    
    async def _test_tiny_dataset(self):
        """测试极小数据集处理"""
        try:
            data_store = {"test": {"adata": self.data_variants['tiny']}}
            
            # 测试多种plot类型
            plot_types = ["spatial", "umap", "heatmap"]
            
            for plot_type in plot_types:
                params = VisualizationParameters(plot_type=plot_type)
                result = await visualize_data("test", data_store, params)
            
            return True
        except Exception as e:
            print(f"❌ 极小数据集测试失败: {str(e)}")
            return False
    
    async def _test_missing_key_data(self):
        """测试缺失关键数据处理"""
        try:
            data_store = {"test": {"adata": self.data_variants['missing_data']}}
            
            # 测试需要UMAP但数据中没有的情况
            params = VisualizationParameters(
                plot_type="umap",
                show_velocity=True,
                show_trajectory=True
            )
            result = await visualize_data("test", data_store, params)
            return True  # 应该优雅处理缺失数据
        except Exception as e:
            print(f"❌ 缺失关键数据测试失败: {str(e)}")
            return False
    
    async def _test_extreme_parameters(self):
        """测试极端参数值"""
        try:
            data_store = {"test": {"adata": self.data_variants['medium']}}
            
            # 测试极端参数值
            extreme_params = [
                VisualizationParameters(
                    plot_type="spatial",
                    figure_size=(1, 1),  # 极小图像
                    dpi=10  # 极低DPI
                ),
                VisualizationParameters(
                    plot_type="umap",
                    figure_size=(50, 50),  # 极大图像
                    dpi=300,  # 极高DPI
                    alpha=0.0  # 完全透明
                ),
                VisualizationParameters(
                    plot_type="spatial_analysis",
                    analysis_sub_type="neighborhood",
                    network_threshold=1000.0,  # 极高阈值
                    show_network=True
                )
            ]
            
            for params in extreme_params:
                result = await visualize_data("test", data_store, params)
            
            return True
        except Exception as e:
            print(f"❌ 极端参数测试失败: {str(e)}")
            return False
    
    async def _test_memory_constraints(self):
        """测试内存约束处理"""
        try:
            data_store = {"test": {"adata": self.data_variants['large']}}
            
            # 测试大数据集的内存使用
            params = VisualizationParameters(
                plot_type="heatmap",
                obs_annotation=["leiden", "batch"]
            )
            
            # 监控内存使用
            import psutil
            process = psutil.Process()
            memory_before = process.memory_info().rss / 1024 / 1024  # MB
            
            result = await visualize_data("test", data_store, params)
            
            memory_after = process.memory_info().rss / 1024 / 1024  # MB
            memory_used = memory_after - memory_before
            
            print(f"📊 内存使用: {memory_used:.1f} MB")
            
            return memory_used < 1000  # 限制在1GB内
        except Exception as e:
            print(f"❌ 内存约束测试失败: {str(e)}")
            return False


    async def run_performance_tests(self):
        """运行性能测试"""
        print("\n⚡ 运行性能测试...")
        
        performance_results = {}
        
        # 测试不同数据集大小的性能
        for dataset_name, adata in self.data_variants.items():
            print(f"  📊 测试数据集: {dataset_name} ({adata.n_obs} cells, {adata.n_vars} genes)")
            
            data_store = {"test": {"adata": adata}}
            
            # 测试不同类型的可视化性能
            visualizations = [
                ("spatial", VisualizationParameters(plot_type="spatial", feature="Gene_0")),
                ("umap", VisualizationParameters(plot_type="umap", feature="leiden")),
                ("heatmap", VisualizationParameters(plot_type="heatmap")),
                ("spatial_interaction", VisualizationParameters(
                    plot_type="spatial_interaction",
                    lr_pairs=[("Gene_0", "Gene_1")]
                )),
                ("integration_check", VisualizationParameters(
                    plot_type="integration_check",
                    batch_key="batch"
                ))
            ]
            
            dataset_performance = {}
            
            for viz_name, params in visualizations:
                start_time = time.time()
                
                try:
                    result = await visualize_data("test", data_store, params)
                    end_time = time.time()
                    execution_time = end_time - start_time
                    
                    dataset_performance[viz_name] = {
                        'success': True,
                        'time': execution_time,
                        'time_per_cell': execution_time / adata.n_obs * 1000  # ms per cell
                    }
                    
                    print(f"    ✅ {viz_name}: {execution_time:.2f}s ({execution_time/adata.n_obs*1000:.2f}ms/cell)")
                    
                except Exception as e:
                    end_time = time.time()
                    execution_time = end_time - start_time
                    
                    dataset_performance[viz_name] = {
                        'success': False,
                        'time': execution_time,
                        'error': str(e)
                    }
                    
                    print(f"    ❌ {viz_name}: FAILED in {execution_time:.2f}s - {str(e)[:50]}...")
            
            performance_results[dataset_name] = dataset_performance
        
        self.performance_results = performance_results
        return performance_results


    async def run_stress_tests(self):
        """运行压力测试"""
        print("\n💪 运行压力测试...")
        
        stress_results = {}
        
        # 压力测试1: 快速连续调用
        stress_results['rapid_calls'] = await self._test_rapid_calls()
        
        # 压力测试2: 并发调用
        stress_results['concurrent_calls'] = await self._test_concurrent_calls()
        
        # 压力测试3: 多种参数组合
        stress_results['parameter_combinations'] = await self._test_parameter_combinations()
        
        # 压力测试4: 极端数据处理
        stress_results['extreme_data'] = await self._test_extreme_data_handling()
        
        return stress_results
    
    async def _test_rapid_calls(self):
        """测试快速连续调用"""
        try:
            data_store = {"test": {"adata": self.data_variants['medium']}}
            params = VisualizationParameters(plot_type="spatial", feature="Gene_0")
            
            start_time = time.time()
            
            # 快速连续调用50次
            for i in range(50):
                result = await visualize_data("test", data_store, params)
            
            end_time = time.time()
            total_time = end_time - start_time
            avg_time = total_time / 50
            
            print(f"    📈 50次连续调用: {total_time:.2f}s (平均 {avg_time:.3f}s/次)")
            
            return avg_time < 1.0  # 每次调用应该在1秒内
            
        except Exception as e:
            print(f"    ❌ 快速调用测试失败: {str(e)}")
            return False
    
    async def _test_concurrent_calls(self):
        """测试并发调用"""
        try:
            data_store = {"test": {"adata": self.data_variants['medium']}}
            
            # 创建不同的可视化任务
            tasks = []
            
            visualizations = [
                VisualizationParameters(plot_type="spatial", feature="Gene_0"),
                VisualizationParameters(plot_type="umap", feature="leiden"),
                VisualizationParameters(plot_type="heatmap"),
                VisualizationParameters(plot_type="spatial_interaction", lr_pairs=[("Gene_0", "Gene_1")]),
                VisualizationParameters(plot_type="integration_check", batch_key="batch")
            ]
            
            # 创建并发任务
            for i, params in enumerate(visualizations * 3):  # 每种类型重复3次
                task = asyncio.create_task(visualize_data("test", data_store, params))
                tasks.append(task)
            
            start_time = time.time()
            
            # 等待所有任务完成
            results = await asyncio.gather(*tasks, return_exceptions=True)
            
            end_time = time.time()
            total_time = end_time - start_time
            
            # 统计成功和失败
            successful = sum(1 for r in results if not isinstance(r, Exception))
            failed = len(results) - successful
            
            print(f"    🚀 15个并发任务: {total_time:.2f}s, 成功: {successful}, 失败: {failed}")
            
            return successful >= 12  # 至少80%成功
            
        except Exception as e:
            print(f"    ❌ 并发调用测试失败: {str(e)}")
            return False
    
    async def _test_parameter_combinations(self):
        """测试多种参数组合"""
        try:
            data_store = {"test": {"adata": self.data_variants['medium']}}
            
            # 生成各种参数组合
            combinations = [
                # 空间图组合
                {"plot_type": "spatial", "feature": "Gene_0", "add_outline": True, "outline_cluster_key": "leiden"},
                {"plot_type": "spatial", "feature": "leiden", "add_outline": True, "outline_color": "red"},
                
                # UMAP组合
                {"plot_type": "umap", "feature": "leiden", "size_by": "Gene_0"},
                {"plot_type": "umap", "feature": "Gene_0", "show_velocity": True, "show_trajectory": True},
                
                # 热图组合
                {"plot_type": "heatmap", "obs_annotation": ["leiden"]},
                {"plot_type": "heatmap", "obs_annotation": ["leiden", "batch", "cell_type"]},
                
                # 空间分析组合
                {"plot_type": "spatial_analysis", "analysis_sub_type": "neighborhood", "cluster_key": "leiden"},
                {"plot_type": "spatial_analysis", "analysis_sub_type": "neighborhood", "show_network": True, "network_threshold": 1.0},
                
                # 新功能组合
                {"plot_type": "spatial_interaction", "lr_pairs": [("Gene_0", "Gene_1"), ("Gene_2", "Gene_3")]},
                {"plot_type": "integration_check", "batch_key": "batch", "integration_method": "TestMethod"}
            ]
            
            successful = 0
            
            for i, combo in enumerate(combinations):
                try:
                    params = VisualizationParameters(**combo)
                    result = await visualize_data("test", data_store, params)
                    successful += 1
                    print(f"    ✅ 组合 {i+1}: {combo['plot_type']} 成功")
                except Exception as e:
                    print(f"    ❌ 组合 {i+1}: {combo['plot_type']} 失败 - {str(e)[:30]}...")
            
            success_rate = successful / len(combinations)
            print(f"    📊 参数组合成功率: {success_rate*100:.1f}% ({successful}/{len(combinations)})")
            
            return success_rate >= 0.8
            
        except Exception as e:
            print(f"    ❌ 参数组合测试失败: {str(e)}")
            return False
    
    async def _test_extreme_data_handling(self):
        """测试极端数据处理"""
        try:
            data_store = {"test": {"adata": self.data_variants['extreme_values']}}
            
            # 测试各种可视化类型对极端值的处理
            visualizations = [
                VisualizationParameters(plot_type="spatial", feature="Gene_0"),
                VisualizationParameters(plot_type="umap", feature="Gene_1"),
                VisualizationParameters(plot_type="heatmap")
            ]
            
            successful = 0
            
            for params in visualizations:
                try:
                    result = await visualize_data("test", data_store, params)
                    successful += 1
                    print(f"    ✅ 极端值处理: {params.plot_type} 成功")
                except Exception as e:
                    print(f"    ❌ 极端值处理: {params.plot_type} 失败 - {str(e)[:50]}...")
            
            return successful >= 2  # 至少2/3成功
            
        except Exception as e:
            print(f"    ❌ 极端数据测试失败: {str(e)}")
            return False


    async def run_compatibility_tests(self):
        """运行兼容性测试"""
        print("\n🔄 运行兼容性测试...")
        
        compatibility_results = {}
        
        # 测试向后兼容性
        compatibility_results['backward_compatibility'] = await self._test_backward_compatibility()
        
        # 测试不同数据格式兼容性
        compatibility_results['data_format_compatibility'] = await self._test_data_format_compatibility()
        
        # 测试库版本兼容性
        compatibility_results['library_compatibility'] = await self._test_library_compatibility()
        
        return compatibility_results
    
    async def _test_backward_compatibility(self):
        """测试向后兼容性"""
        try:
            data_store = {"test": {"adata": self.data_variants['medium']}}
            
            # 测试旧式参数调用
            old_style_params = [
                # 只使用基础参数
                VisualizationParameters(plot_type="spatial"),
                VisualizationParameters(plot_type="umap"),
                VisualizationParameters(plot_type="heatmap"),
                VisualizationParameters(plot_type="violin"),
                
                # 使用旧的参数名
                VisualizationParameters(plot_type="spatial", feature="Gene_0", colormap="viridis"),
                VisualizationParameters(plot_type="umap", feature="leiden", show_legend=True)
            ]
            
            successful = 0
            
            for params in old_style_params:
                try:
                    result = await visualize_data("test", data_store, params)
                    successful += 1
                except Exception as e:
                    print(f"    ❌ 向后兼容性失败: {params.plot_type} - {str(e)[:30]}...")
            
            compatibility_rate = successful / len(old_style_params)
            print(f"    📊 向后兼容性: {compatibility_rate*100:.1f}% ({successful}/{len(old_style_params)})")
            
            return compatibility_rate == 1.0  # 100% 向后兼容
            
        except Exception as e:
            print(f"    ❌ 向后兼容性测试失败: {str(e)}")
            return False
    
    async def _test_data_format_compatibility(self):
        """测试不同数据格式兼容性"""
        try:
            # 测试不同类型的数据集
            test_datasets = [
                ('tiny', self.data_variants['tiny']),
                ('sparse', self.data_variants['sparse']),
                ('visium', self.data_variants['visium']),
                ('imbalanced', self.data_variants['imbalanced'])
            ]
            
            successful = 0
            
            for name, adata in test_datasets:
                try:
                    data_store = {"test": {"adata": adata}}
                    params = VisualizationParameters(plot_type="spatial", feature="Gene_0")
                    result = await visualize_data("test", data_store, params)
                    successful += 1
                    print(f"    ✅ 数据格式兼容: {name}")
                except Exception as e:
                    print(f"    ❌ 数据格式不兼容: {name} - {str(e)[:40]}...")
            
            compatibility_rate = successful / len(test_datasets)
            print(f"    📊 数据格式兼容性: {compatibility_rate*100:.1f}%")
            
            return compatibility_rate >= 0.75  # 至少75%兼容
            
        except Exception as e:
            print(f"    ❌ 数据格式兼容性测试失败: {str(e)}")
            return False
    
    async def _test_library_compatibility(self):
        """测试库版本兼容性"""
        try:
            # 检查关键库的可用性
            libraries_status = {}
            
            # 测试scanpy
            try:
                import scanpy as sc
                libraries_status['scanpy'] = f"✅ {sc.__version__}"
            except ImportError:
                libraries_status['scanpy'] = "❌ 不可用"
            
            # 测试matplotlib
            try:
                import matplotlib
                libraries_status['matplotlib'] = f"✅ {matplotlib.__version__}"
            except ImportError:
                libraries_status['matplotlib'] = "❌ 不可用"
            
            # 测试numpy
            try:
                import numpy as np
                libraries_status['numpy'] = f"✅ {np.__version__}"
            except ImportError:
                libraries_status['numpy'] = "❌ 不可用"
            
            # 测试pandas
            try:
                import pandas as pd
                libraries_status['pandas'] = f"✅ {pd.__version__}"
            except ImportError:
                libraries_status['pandas'] = "❌ 不可用"
            
            # 测试可选库
            try:
                import networkx as nx
                libraries_status['networkx'] = f"✅ {nx.__version__}"
            except ImportError:
                libraries_status['networkx'] = "⚠️ 可选库不可用"
            
            try:
                import scvelo
                libraries_status['scvelo'] = f"✅ {scvelo.__version__}"
            except ImportError:
                libraries_status['scvelo'] = "⚠️ 可选库不可用"
            
            print("    📚 库兼容性检查:")
            for lib, status in libraries_status.items():
                print(f"      {lib}: {status}")
            
            # 计算必需库可用率
            required_libs = ['scanpy', 'matplotlib', 'numpy', 'pandas']
            available_required = sum(1 for lib in required_libs if "✅" in libraries_status.get(lib, ""))
            
            return available_required == len(required_libs)
            
        except Exception as e:
            print(f"    ❌ 库兼容性测试失败: {str(e)}")
            return False


    def generate_test_report(self, all_results):
        """生成详细的测试报告"""
        print("\n" + "="*80)
        print("🔬 超全面测试报告")
        print("="*80)
        
        # 总体统计
        total_tests = 0
        passed_tests = 0
        
        print("\n📊 测试结果概览:")
        
        for test_category, results in all_results.items():
            if isinstance(results, dict):
                category_total = len(results)
                category_passed = sum(1 for v in results.values() if v is True)
                
                total_tests += category_total
                passed_tests += category_passed
                
                success_rate = (category_passed / category_total * 100) if category_total > 0 else 0
                
                print(f"  {test_category}: {category_passed}/{category_total} ({success_rate:.1f}%)")
                
                # 显示失败的测试
                failed_tests = [k for k, v in results.items() if v is False]
                if failed_tests:
                    print(f"    ❌ 失败: {', '.join(failed_tests)}")
        
        overall_success_rate = (passed_tests / total_tests * 100) if total_tests > 0 else 0
        
        print(f"\n🎯 总体成功率: {overall_success_rate:.1f}% ({passed_tests}/{total_tests})")
        
        # 性能报告
        if hasattr(self, 'performance_results') and self.performance_results:
            print("\n⚡ 性能分析:")
            
            for dataset_name, dataset_perf in self.performance_results.items():
                print(f"  📊 数据集 {dataset_name}:")
                
                for viz_type, perf_data in dataset_perf.items():
                    if perf_data.get('success'):
                        time_str = f"{perf_data['time']:.2f}s"
                        if 'time_per_cell' in perf_data:
                            time_str += f" ({perf_data['time_per_cell']:.2f}ms/cell)"
                        print(f"    ✅ {viz_type}: {time_str}")
                    else:
                        print(f"    ❌ {viz_type}: FAILED")
        
        # 数据集信息
        print("\n📈 测试数据集信息:")
        for name, adata in self.data_variants.items():
            print(f"  {name}: {adata.n_obs} cells, {adata.n_vars} genes")
        
        # 建议和结论
        print("\n💡 测试结论和建议:")
        
        if overall_success_rate >= 95:
            print("  🎉 优秀！系统表现非常稳定，所有功能都工作正常。")
        elif overall_success_rate >= 85:
            print("  ✅ 良好！大部分功能正常，少数边界情况需要注意。")
        elif overall_success_rate >= 70:
            print("  ⚠️  一般。有一些问题需要解决，建议优化错误处理。")
        else:
            print("  ❌ 需要改进。存在较多问题，需要重点修复。")
        
        print("\n" + "="*80)
        
        return overall_success_rate


    async def run_all_tests(self):
        """运行所有测试"""
        print("🚀 开始超全面测试套件")
        print("="*60)
        
        # 创建测试数据
        self.create_various_test_datasets()
        
        # 运行所有测试类别
        all_results = {}
        
        # 1. 边界情况测试
        all_results['edge_cases'] = await self.run_edge_case_tests()
        
        # 2. 性能测试
        all_results['performance'] = await self.run_performance_tests()
        
        # 3. 压力测试
        all_results['stress_tests'] = await self.run_stress_tests()
        
        # 4. 兼容性测试
        all_results['compatibility'] = await self.run_compatibility_tests()
        
        # 生成报告
        final_score = self.generate_test_report(all_results)
        
        return final_score, all_results


async def main():
    """主测试函数"""
    test_suite = ComprehensiveTestSuite()
    
    start_time = time.time()
    final_score, results = await test_suite.run_all_tests()
    end_time = time.time()
    
    print(f"\n⏱️ 总测试时间: {end_time - start_time:.2f}秒")
    print(f"🏆 最终评分: {final_score:.1f}%")
    
    # 根据评分返回适当的退出码
    if final_score >= 90:
        exit_code = 0  # 优秀
    elif final_score >= 80:
        exit_code = 1  # 良好但有改进空间
    else:
        exit_code = 2  # 需要重大改进
    
    return exit_code


if __name__ == "__main__":
    exit_code = asyncio.run(main())
    sys.exit(exit_code)