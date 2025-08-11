#!/usr/bin/env python3
"""
最终特性验证测试 - 专门验证所有新增功能的工作情况
"""

import asyncio
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from chatspatial.models.data import VisualizationParameters
from chatspatial.tools.visualization import visualize_data


def create_test_data():
    """创建标准测试数据"""
    np.random.seed(42)
    n_cells, n_genes = 200, 500
    
    # 创建表达数据
    X = np.random.poisson(3, (n_cells, n_genes)).astype(np.float32)
    
    adata = sc.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"Cell_{i:03d}" for i in range(n_cells)]),
        var=pd.DataFrame(index=[f"Gene_{i:03d}" for i in range(n_genes)])
    )
    
    # 空间坐标 - 创建明显的cluster分布
    centers = np.array([[0, 0], [10, 0], [5, 8], [-5, 5]])
    coords = []
    clusters = []
    
    for i in range(n_cells):
        cluster_idx = i % 4
        center = centers[cluster_idx]
        coord = center + np.random.normal(0, 2, 2)
        coords.append(coord)
        clusters.append(f"Cluster_{cluster_idx}")
    
    adata.obsm['spatial'] = np.array(coords)
    adata.obs['leiden'] = pd.Categorical(clusters)
    adata.obs['batch'] = pd.Categorical(['Batch1', 'Batch2'] * (n_cells // 2))
    adata.obs['cell_type'] = pd.Categorical(['TypeA', 'TypeB', 'TypeC'] * (n_cells // 3 + 1))[:n_cells]
    
    # 计算基本embeddings
    sc.pp.highly_variable_genes(adata, n_top_genes=50)
    # Add basic preprocessing for PCA
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    
    # 添加邻域富集结果
    n_clusters = len(adata.obs['leiden'].cat.categories)
    enrichment_matrix = np.random.normal(0, 2, (n_clusters, n_clusters))
    enrichment_matrix = (enrichment_matrix + enrichment_matrix.T) / 2
    np.fill_diagonal(enrichment_matrix, np.abs(np.diagonal(enrichment_matrix)) + 2)
    
    adata.uns['leiden_nhood_enrichment'] = {'zscore': enrichment_matrix}
    
    return adata


async def validate_all_new_features():
    """验证所有新增功能"""
    print("🎯 最终特性验证测试")
    print("="*50)
    
    adata = create_test_data()
    data_store = {"test": {"adata": adata}}
    
    test_results = {}
    
    print("\n1️⃣ 测试新增参数模型...")
    try:
        # 验证所有新参数都可以正确设置
        enhanced_params = VisualizationParameters(
            plot_type="spatial",
            feature="Gene_001",
            # 新的空间增强参数
            add_outline=True,
            outline_color="red",
            outline_width=2.0,
            outline_cluster_key="leiden",
            # 新的UMAP增强参数
            size_by="Gene_002",
            show_velocity=True,
            show_trajectory=True,
            velocity_scale=1.5,
            # 新的热图增强参数
            obs_annotation=["leiden", "batch"],
            var_annotation=["highly_variable"],
            # 新的网络参数
            show_network=True,
            network_threshold=1.0,
            network_layout="spring",
            # 新的整合参数
            batch_key="batch",
            integration_method="TestMethod"
        )
        test_results['enhanced_parameters'] = True
        print("✅ 新增参数模型验证成功")
    except Exception as e:
        test_results['enhanced_parameters'] = False
        print(f"❌ 参数模型验证失败: {str(e)}")
    
    print("\n2️⃣ 测试空间图轮廓叠加...")
    try:
        params = VisualizationParameters(
            plot_type="spatial",
            feature="Gene_001",
            add_outline=True,
            outline_color="black",
            outline_cluster_key="leiden"
        )
        result = await visualize_data("test", data_store, params)
        test_results['spatial_outline'] = True
        print("✅ 空间图轮廓叠加功能正常")
    except Exception as e:
        test_results['spatial_outline'] = False
        print(f"❌ 空间图轮廓叠加失败: {str(e)}")
    
    print("\n3️⃣ 测试空间交互可视化...")
    try:
        params = VisualizationParameters(
            plot_type="spatial_interaction",
            lr_pairs=[("Gene_001", "Gene_002"), ("Gene_003", "Gene_004")]
        )
        result = await visualize_data("test", data_store, params)
        test_results['spatial_interaction'] = True
        print("✅ 空间交互可视化功能正常")
    except Exception as e:
        test_results['spatial_interaction'] = False
        print(f"❌ 空间交互可视化失败: {str(e)}")
    
    print("\n4️⃣ 测试UMAP双重编码...")
    try:
        params = VisualizationParameters(
            plot_type="umap",
            feature="leiden",
            size_by="Gene_001"
        )
        result = await visualize_data("test", data_store, params)
        test_results['umap_dual_encoding'] = True
        print("✅ UMAP双重编码功能正常")
    except Exception as e:
        test_results['umap_dual_encoding'] = False
        print(f"❌ UMAP双重编码失败: {str(e)}")
    
    print("\n5️⃣ 测试UMAP速度/轨迹叠加...")
    try:
        params = VisualizationParameters(
            plot_type="umap",
            feature="leiden",
            show_velocity=True,
            show_trajectory=True
        )
        result = await visualize_data("test", data_store, params)
        test_results['umap_overlays'] = True
        print("✅ UMAP速度/轨迹叠加功能正常（优雅处理缺失数据）")
    except Exception as e:
        test_results['umap_overlays'] = False
        print(f"❌ UMAP叠加功能失败: {str(e)}")
    
    print("\n6️⃣ 测试增强热图注释...")
    try:
        params = VisualizationParameters(
            plot_type="heatmap",
            obs_annotation=["leiden", "batch", "cell_type"]
        )
        result = await visualize_data("test", data_store, params)
        test_results['heatmap_annotations'] = True
        print("✅ 增强热图注释功能正常")
    except Exception as e:
        test_results['heatmap_annotations'] = False
        print(f"❌ 增强热图注释失败: {str(e)}")
    
    print("\n7️⃣ 测试整合评估可视化...")
    try:
        params = VisualizationParameters(
            plot_type="integration_check",
            batch_key="batch",
            integration_method="TestIntegration"
        )
        result = await visualize_data("test", data_store, params)
        test_results['integration_check'] = True
        print("✅ 整合评估可视化功能正常")
    except Exception as e:
        test_results['integration_check'] = False
        print(f"❌ 整合评估可视化失败: {str(e)}")
    
    print("\n8️⃣ 测试邻域富集热图模式...")
    try:
        params = VisualizationParameters(
            plot_type="spatial_analysis",
            analysis_sub_type="neighborhood",
            cluster_key="leiden"
        )
        result = await visualize_data("test", data_store, params)
        test_results['neighborhood_heatmap'] = True
        print("✅ 邻域富集热图模式正常")
    except Exception as e:
        test_results['neighborhood_heatmap'] = False
        print(f"❌ 邻域富集热图模式失败: {str(e)}")
    
    print("\n9️⃣ 测试邻域富集网络模式...")
    try:
        params = VisualizationParameters(
            plot_type="spatial_analysis",
            analysis_sub_type="neighborhood",
            cluster_key="leiden",
            show_network=True,
            network_threshold=1.5,
            network_layout="spring"
        )
        result = await visualize_data("test", data_store, params)
        test_results['neighborhood_network'] = True
        print("✅ 邻域富集网络模式功能正常")
    except Exception as e:
        test_results['neighborhood_network'] = False
        print(f"❌ 邻域富集网络模式失败: {str(e)}")
    
    print("\n🔟 测试向后兼容性...")
    try:
        # 测试旧式简单参数
        old_params = [
            VisualizationParameters(plot_type="spatial"),
            VisualizationParameters(plot_type="umap"),
            VisualizationParameters(plot_type="heatmap"),
            VisualizationParameters(plot_type="spatial", feature="Gene_001"),
        ]
        
        all_passed = True
        for params in old_params:
            result = await visualize_data("test", data_store, params)
        
        test_results['backward_compatibility'] = True
        print("✅ 向后兼容性测试通过")
    except Exception as e:
        test_results['backward_compatibility'] = False
        print(f"❌ 向后兼容性测试失败: {str(e)}")
    
    return test_results


async def test_error_handling():
    """测试错误处理能力"""
    print("\n🛡️ 错误处理能力测试")
    print("="*30)
    
    adata = create_test_data()
    data_store = {"test": {"adata": adata}}
    
    error_tests = {}
    
    print("\n📍 测试无效plot_type处理...")
    try:
        params = VisualizationParameters(plot_type="invalid_plot_type")
        result = await visualize_data("test", data_store, params)
        error_tests['invalid_plot_type'] = False  # Should have failed
        print("❌ 应该抛出错误但没有")
    except Exception as e:
        error_tests['invalid_plot_type'] = True  # Expected behavior
        print("✅ 正确处理了无效plot_type")
    
    print("\n📍 测试不存在的特征处理...")
    try:
        params = VisualizationParameters(
            plot_type="spatial",
            feature="NonExistentGene_XYZ"
        )
        result = await visualize_data("test", data_store, params)
        error_tests['nonexistent_feature'] = True  # Should handle gracefully
        print("✅ 优雅处理了不存在的特征")
    except Exception as e:
        error_tests['nonexistent_feature'] = False
        print(f"❌ 未能优雅处理不存在的特征: {str(e)}")
    
    print("\n📍 测试不存在数据集处理...")
    try:
        params = VisualizationParameters(plot_type="spatial")
        result = await visualize_data("nonexistent_dataset", data_store, params)
        error_tests['nonexistent_dataset'] = False  # Should have failed
        print("❌ 应该抛出错误但没有")
    except Exception as e:
        error_tests['nonexistent_dataset'] = True  # Expected behavior
        print("✅ 正确处理了不存在的数据集")
    
    print("\n📍 测试缺失轮廓cluster_key处理...")
    try:
        params = VisualizationParameters(
            plot_type="spatial",
            feature="Gene_001",
            add_outline=True,
            outline_cluster_key="NonExistentKey"
        )
        result = await visualize_data("test", data_store, params)
        error_tests['missing_outline_key'] = True  # Should handle gracefully
        print("✅ 优雅处理了缺失的轮廓cluster_key")
    except Exception as e:
        error_tests['missing_outline_key'] = False
        print(f"❌ 未能优雅处理缺失的轮廓cluster_key: {str(e)}")
    
    return error_tests


def generate_final_report(feature_results, error_results):
    """生成最终测试报告"""
    print("\n" + "="*60)
    print("🏆 最终特性验证报告")
    print("="*60)
    
    # 特性测试结果
    print("\n🎯 新增特性测试结果:")
    total_features = len(feature_results)
    passed_features = sum(feature_results.values())
    
    for feature, passed in feature_results.items():
        status = "✅" if passed else "❌"
        print(f"  {status} {feature}")
    
    feature_success_rate = (passed_features / total_features * 100) if total_features > 0 else 0
    print(f"\n📊 特性成功率: {feature_success_rate:.1f}% ({passed_features}/{total_features})")
    
    # 错误处理结果
    print("\n🛡️ 错误处理测试结果:")
    total_errors = len(error_results)
    passed_errors = sum(error_results.values())
    
    for error_test, passed in error_results.items():
        status = "✅" if passed else "❌"
        print(f"  {status} {error_test}")
    
    error_success_rate = (passed_errors / total_errors * 100) if total_errors > 0 else 0
    print(f"\n📊 错误处理成功率: {error_success_rate:.1f}% ({passed_errors}/{total_errors})")
    
    # 总体评估
    overall_tests = total_features + total_errors
    overall_passed = passed_features + passed_errors
    overall_success_rate = (overall_passed / overall_tests * 100) if overall_tests > 0 else 0
    
    print(f"\n🎯 总体成功率: {overall_success_rate:.1f}% ({overall_passed}/{overall_tests})")
    
    # 结论
    print("\n💡 最终结论:")
    if overall_success_rate >= 95:
        print("  🌟 卓越！所有新功能都完美工作，系统非常健壮！")
        grade = "A+"
    elif overall_success_rate >= 90:
        print("  🎉 优秀！绝大部分功能正常，系统质量很高！") 
        grade = "A"
    elif overall_success_rate >= 80:
        print("  ✅ 良好！大部分功能正常，少数需要微调。")
        grade = "B+"
    elif overall_success_rate >= 70:
        print("  ⚠️  一般。主要功能可用，但需要改进一些问题。")
        grade = "B"
    else:
        print("  ❌ 需要大量改进。许多功能存在问题。")
        grade = "C"
    
    print(f"\n🏅 系统评级: {grade}")
    print("="*60)
    
    return overall_success_rate, grade


async def main():
    """主函数"""
    print("🚀 启动最终特性验证测试")
    
    # 运行特性验证
    feature_results = await validate_all_new_features()
    
    # 运行错误处理测试
    error_results = await test_error_handling()
    
    # 生成最终报告
    success_rate, grade = generate_final_report(feature_results, error_results)
    
    return success_rate >= 80  # 80%以上算通过


if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1)