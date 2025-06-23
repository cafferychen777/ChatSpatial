#!/usr/bin/env python3
"""Demo GSEA visualization functionality with all plot types"""

import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from chatspatial.tools.visualization import create_gsea_visualization, VisualizationParameters
import asyncio

async def demo_gsea_visualization():
    # Load the data
    print("Loading ST mouse brain data...")
    adata = ad.read_h5ad("/Users/apple/Research/SpatialTrans_MCP/chatspatial/data/ST_mouse_brain.h5ad")
    print(f"Data shape: {adata.shape}")
    
    # Create comprehensive GSEA results
    print("\nCreating comprehensive GSEA results...")
    
    # 1. Detailed pathway results for enrichment plots
    pathways = {
        "GO_SYNAPTIC_SIGNALING": {
            "NES": 2.35,
            "pval": 0.001,
            "FDR": 0.01,
            "es": np.cumsum(np.random.randn(1000) * 0.01 + 0.002),
            "positions": sorted(np.random.choice(1000, 50, replace=False))
        },
        "GO_NEUROTRANSMITTER_TRANSPORT": {
            "NES": 2.12,
            "pval": 0.002,
            "FDR": 0.02
        },
        "GO_ION_CHANNEL_ACTIVITY": {
            "NES": 1.98,
            "pval": 0.003,
            "FDR": 0.03
        },
        "GO_AXON_GUIDANCE": {
            "NES": 1.87,
            "pval": 0.005,
            "FDR": 0.04
        },
        "GO_GABA_SIGNALING": {
            "NES": 1.76,
            "pval": 0.008,
            "FDR": 0.045
        },
        "GO_GLUTAMATE_RECEPTOR": {
            "NES": 1.65,
            "pval": 0.010,
            "FDR": 0.048
        },
        "GO_CALCIUM_SIGNALING": {
            "NES": 1.54,
            "pval": 0.015,
            "FDR": 0.055
        },
        "GO_NEURONAL_DEVELOPMENT": {
            "NES": 1.43,
            "pval": 0.025,
            "FDR": 0.08
        },
        "GO_CELL_CYCLE": {
            "NES": -1.85,
            "pval": 0.012,
            "FDR": 0.05,
            "es": np.cumsum(np.random.randn(1000) * 0.01 - 0.001),
            "positions": sorted(np.random.choice(range(800, 1000), 35, replace=False))
        },
        "GO_IMMUNE_RESPONSE": {
            "NES": -1.65,
            "pval": 0.025,
            "FDR": 0.08
        },
        "KEGG_METABOLISM": {
            "NES": 1.32,
            "pval": 0.035,
            "FDR": 0.10
        },
        "REACTOME_TRANSCRIPTION": {
            "NES": -1.23,
            "pval": 0.045,
            "FDR": 0.12
        }
    }
    
    # Normalize enrichment scores
    for pathway, data in pathways.items():
        if 'es' in data:
            es = data['es']
            data['es'] = es / np.abs(es).max()
    
    adata.uns['gsea_results'] = pathways
    
    # 2. Create cluster-specific results for dotplot
    clusters = ['Excitatory Neurons', 'Inhibitory Neurons', 'Astrocytes', 'Oligodendrocytes', 'Microglia']
    cluster_results = {}
    
    # Define cluster-specific enrichment patterns
    cluster_patterns = {
        'Excitatory Neurons': {'synaptic': 0.5, 'neurotrans': 0.4, 'ion': 0.3, 'immune': -0.8},
        'Inhibitory Neurons': {'synaptic': 0.3, 'gaba': 0.6, 'ion': 0.4, 'immune': -0.6},
        'Astrocytes': {'calcium': 0.5, 'glutamate': 0.4, 'metabolism': 0.3, 'synaptic': -0.3},
        'Oligodendrocytes': {'development': 0.5, 'transcription': 0.4, 'cell_cycle': 0.3, 'synaptic': -0.5},
        'Microglia': {'immune': 0.8, 'cell_cycle': 0.3, 'synaptic': -0.8, 'neurotrans': -0.6}
    }
    
    pathway_mapping = {
        'synaptic': 'GO_SYNAPTIC_SIGNALING',
        'neurotrans': 'GO_NEUROTRANSMITTER_TRANSPORT',
        'ion': 'GO_ION_CHANNEL_ACTIVITY',
        'gaba': 'GO_GABA_SIGNALING',
        'glutamate': 'GO_GLUTAMATE_RECEPTOR',
        'calcium': 'GO_CALCIUM_SIGNALING',
        'development': 'GO_NEURONAL_DEVELOPMENT',
        'cell_cycle': 'GO_CELL_CYCLE',
        'immune': 'GO_IMMUNE_RESPONSE',
        'metabolism': 'KEGG_METABOLISM',
        'transcription': 'REACTOME_TRANSCRIPTION'
    }
    
    for cluster in clusters:
        cluster_results[cluster] = {}
        pattern = cluster_patterns.get(cluster, {})
        
        for short_name, full_name in pathway_mapping.items():
            if full_name in pathways:
                base_nes = pathways[full_name]['NES']
                # Apply cluster-specific modification
                modifier = pattern.get(short_name, np.random.uniform(-0.2, 0.2))
                cluster_nes = base_nes + modifier
                
                # Calculate p-value based on NES
                if abs(cluster_nes) > 2:
                    pval = np.random.uniform(0.001, 0.01)
                elif abs(cluster_nes) > 1.5:
                    pval = np.random.uniform(0.01, 0.05)
                else:
                    pval = np.random.uniform(0.05, 0.2)
                
                cluster_results[cluster][full_name] = {
                    "NES": cluster_nes,
                    "pval": pval
                }
    
    adata.uns['gsea_cluster_results'] = cluster_results
    
    # 3. Create spatial pathway scores
    if 'spatial' in adata.obsm:
        coords = adata.obsm['spatial']
        
        # Synaptic signaling - higher in central region
        x_center = coords[:, 0].mean()
        y_center = coords[:, 1].mean()
        distances = np.sqrt((coords[:, 0] - x_center)**2 + (coords[:, 1] - y_center)**2)
        synaptic_scores = 2 * np.exp(-distances / (distances.max() * 0.4))
        synaptic_scores = (synaptic_scores - synaptic_scores.mean()) / synaptic_scores.std()
        adata.obs['GO_SYNAPTIC_SIGNALING_score'] = synaptic_scores
        
        # Immune response - higher in periphery
        immune_scores = 1.5 * (1 - np.exp(-distances / (distances.max() * 0.5)))
        immune_scores = (immune_scores - immune_scores.mean()) / immune_scores.std()
        adata.obs['GO_IMMUNE_RESPONSE_score'] = immune_scores
        
        # Cell cycle - patchy distribution
        cell_cycle_scores = np.zeros(len(adata))
        for _ in range(5):
            x_c = np.random.uniform(coords[:, 0].min(), coords[:, 0].max())
            y_c = np.random.uniform(coords[:, 1].min(), coords[:, 1].max())
            dist = np.sqrt((coords[:, 0] - x_c)**2 + (coords[:, 1] - y_c)**2)
            cell_cycle_scores += np.exp(-dist / (distances.max() * 0.2))
        cell_cycle_scores = (cell_cycle_scores - cell_cycle_scores.mean()) / cell_cycle_scores.std()
        adata.obs['GO_CELL_CYCLE_score'] = cell_cycle_scores
    
    # Create all visualization types
    print("\nCreating all GSEA visualization types...")
    
    # 1. Classic enrichment plot - positive enrichment
    print("\n1. Classic GSEA enrichment plot (positive enrichment)")
    params1 = VisualizationParameters(
        plot_type="gsea",
        gsea_plot_type="enrichment_plot",
        feature="GO_SYNAPTIC_SIGNALING",
        figure_size=(10, 8)
    )
    fig1 = await create_gsea_visualization(adata, params1)
    fig1.savefig('gsea_demo_1_enrichment_positive.png', dpi=150, bbox_inches='tight')
    plt.close(fig1)
    
    # 2. Classic enrichment plot - negative enrichment
    print("2. Classic GSEA enrichment plot (negative enrichment)")
    params2 = VisualizationParameters(
        plot_type="gsea",
        gsea_plot_type="enrichment_plot",
        feature="GO_CELL_CYCLE",
        figure_size=(10, 8)
    )
    fig2 = await create_gsea_visualization(adata, params2)
    fig2.savefig('gsea_demo_2_enrichment_negative.png', dpi=150, bbox_inches='tight')
    plt.close(fig2)
    
    # 3. Top pathways barplot
    print("3. Top enriched pathways barplot")
    params3 = VisualizationParameters(
        plot_type="gsea",
        gsea_plot_type="barplot",
        n_top_pathways=12,
        figure_size=(10, 8)
    )
    fig3 = await create_gsea_visualization(adata, params3)
    fig3.savefig('gsea_demo_3_barplot.png', dpi=150, bbox_inches='tight')
    plt.close(fig3)
    
    # 4. Multi-cluster dotplot
    print("4. Multi-cluster enrichment dotplot")
    params4 = VisualizationParameters(
        plot_type="gsea",
        gsea_plot_type="dotplot",
        gsea_results_key="gsea_cluster_results",
        figure_size=(14, 8)
    )
    fig4 = await create_gsea_visualization(adata, params4)
    fig4.savefig('gsea_demo_4_dotplot.png', dpi=150, bbox_inches='tight')
    plt.close(fig4)
    
    # 5. Spatial distribution - synaptic signaling
    print("5. Spatial distribution of synaptic signaling")
    params5 = VisualizationParameters(
        plot_type="gsea",
        gsea_plot_type="spatial",
        feature="GO_SYNAPTIC_SIGNALING",
        figure_size=(10, 8),
        colormap="RdBu_r"
    )
    fig5 = await create_gsea_visualization(adata, params5)
    fig5.savefig('gsea_demo_5_spatial_synaptic.png', dpi=150, bbox_inches='tight')
    plt.close(fig5)
    
    # 6. Create combined figure showing all spatial patterns
    print("6. Combined spatial patterns")
    fig6, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    pathways_spatial = [
        ('GO_SYNAPTIC_SIGNALING_score', 'Synaptic Signaling', 'RdBu_r'),
        ('GO_IMMUNE_RESPONSE_score', 'Immune Response', 'PuOr_r'),
        ('GO_CELL_CYCLE_score', 'Cell Cycle', 'BrBG_r')
    ]
    
    for idx, (score_col, title, cmap) in enumerate(pathways_spatial):
        ax = axes[idx]
        scores = adata.obs[score_col].values
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=scores, 
                           cmap=cmap, s=50, alpha=0.8,
                           vmin=np.percentile(scores, 5),
                           vmax=np.percentile(scores, 95))
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('Spatial X')
        ax.set_ylabel('Spatial Y')
        ax.set_aspect('equal')
        plt.colorbar(scatter, ax=ax, label='Score')
    
    plt.suptitle('Spatial Distribution of Pathway Activities', fontsize=16, fontweight='bold')
    plt.tight_layout()
    fig6.savefig('gsea_demo_6_spatial_combined.png', dpi=150, bbox_inches='tight')
    plt.close(fig6)
    
    print("\n" + "="*50)
    print("GSEA visualization demo completed!")
    print("Generated files:")
    print("  - gsea_demo_1_enrichment_positive.png")
    print("  - gsea_demo_2_enrichment_negative.png")
    print("  - gsea_demo_3_barplot.png")
    print("  - gsea_demo_4_dotplot.png")
    print("  - gsea_demo_5_spatial_synaptic.png")
    print("  - gsea_demo_6_spatial_combined.png")
    print("="*50)

if __name__ == "__main__":
    asyncio.run(demo_gsea_visualization())