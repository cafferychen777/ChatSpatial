"""
Example script demonstrating cell communication analysis functionality in ChatSpatial.

This script shows how to:
1. Load spatial transcriptomics data
2. Preprocess the data
3. Analyze cell-cell communication using COMMOT and SpatialDM
4. Visualize communication networks
5. Identify communication patterns
"""

import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Import ChatSpatial modules
from chatspatial.models.data import CellCommunicationParameters
from chatspatial.tools.cell_communication import analyze_cell_communication
from chatspatial.tools.preprocessing import preprocess_data
from chatspatial.models.data import AnalysisParameters


def create_synthetic_communication_data(n_spots=400, n_genes=150):
    """Create synthetic spatial transcriptomics data with communication patterns"""
    print("Creating synthetic spatial transcriptomics data with communication patterns...")
    
    np.random.seed(42)
    
    # Create spatial coordinates with distinct regions
    x_coords = np.random.uniform(0, 20, n_spots)
    y_coords = np.random.uniform(0, 20, n_spots)
    coords = np.column_stack([x_coords, y_coords])
    
    # Define common ligand-receptor pairs
    lr_pairs = [
        ('TGFB1', 'TGFBR1'), ('TGFB1', 'TGFBR2'),
        ('PDGFA', 'PDGFRA'), ('PDGFB', 'PDGFRB'),
        ('VEGFA', 'FLT1'), ('VEGFA', 'KDR'),
        ('TNF', 'TNFRSF1A'), ('TNF', 'TNFRSF1B'),
        ('IL1B', 'IL1R1'), ('IL6', 'IL6R'),
        ('CXCL12', 'CXCR4'), ('CCL2', 'CCR2'),
        ('WNT3A', 'FZD1'), ('WNT5A', 'FZD2'),
        ('NOTCH1', 'DLL1'), ('NOTCH2', 'JAG1')
    ]
    
    # Create base expression matrix
    base_expression = np.random.negative_binomial(8, 0.4, size=(n_spots, n_genes))
    
    # Add ligand-receptor genes to the gene list
    lr_genes = []
    for ligand, receptor in lr_pairs:
        lr_genes.extend([ligand, receptor])
    
    # Remove duplicates and ensure we have enough genes
    lr_genes = list(set(lr_genes))
    other_genes = [f"Gene_{i:03d}" for i in range(len(lr_genes), n_genes)]
    all_genes = lr_genes + other_genes
    
    # Create communication regions
    print("Creating communication hotspots...")
    
    # Region 1: High TGFB signaling (center-left)
    center1 = (5, 10)
    region1_mask = ((x_coords - center1[0])**2 + (y_coords - center1[1])**2) < 16
    if 'TGFB1' in all_genes and 'TGFBR1' in all_genes:
        tgfb1_idx = all_genes.index('TGFB1')
        tgfbr1_idx = all_genes.index('TGFBR1')
        base_expression[region1_mask, tgfb1_idx] *= 4  # High ligand expression
        base_expression[region1_mask, tgfbr1_idx] *= 3  # High receptor expression
    
    # Region 2: High VEGF signaling (center-right)
    center2 = (15, 10)
    region2_mask = ((x_coords - center2[0])**2 + (y_coords - center2[1])**2) < 16
    if 'VEGFA' in all_genes and 'FLT1' in all_genes:
        vegfa_idx = all_genes.index('VEGFA')
        flt1_idx = all_genes.index('FLT1')
        base_expression[region2_mask, vegfa_idx] *= 4
        base_expression[region2_mask, flt1_idx] *= 3
    
    # Region 3: High inflammatory signaling (top)
    region3_mask = y_coords > 15
    if 'TNF' in all_genes and 'TNFRSF1A' in all_genes:
        tnf_idx = all_genes.index('TNF')
        tnfrsf1a_idx = all_genes.index('TNFRSF1A')
        base_expression[region3_mask, tnf_idx] *= 3
        base_expression[region3_mask, tnfrsf1a_idx] *= 3
    
    # Region 4: High chemokine signaling (bottom)
    region4_mask = y_coords < 5
    if 'CXCL12' in all_genes and 'CXCR4' in all_genes:
        cxcl12_idx = all_genes.index('CXCL12')
        cxcr4_idx = all_genes.index('CXCR4')
        base_expression[region4_mask, cxcl12_idx] *= 3
        base_expression[region4_mask, cxcr4_idx] *= 3
    
    # Create AnnData object
    adata = sc.AnnData(base_expression.astype(float))
    adata.obsm['spatial'] = coords
    adata.var_names = all_genes[:n_genes]
    adata.obs_names = [f"Spot_{i:04d}" for i in range(n_spots)]
    
    # Add metadata about communication regions
    adata.obs['communication_region'] = 'background'
    adata.obs.loc[region1_mask, 'communication_region'] = 'TGFB_region'
    adata.obs.loc[region2_mask, 'communication_region'] = 'VEGF_region'
    adata.obs.loc[region3_mask, 'communication_region'] = 'inflammatory_region'
    adata.obs.loc[region4_mask, 'communication_region'] = 'chemokine_region'
    adata.obs['communication_region'] = adata.obs['communication_region'].astype('category')
    
    # Mark which genes are ligands/receptors
    adata.var['gene_type'] = 'other'
    for ligand, receptor in lr_pairs:
        if ligand in adata.var.index:
            adata.var.loc[ligand, 'gene_type'] = 'ligand'
        if receptor in adata.var.index:
            adata.var.loc[receptor, 'gene_type'] = 'receptor'
    adata.var['gene_type'] = adata.var['gene_type'].astype('category')
    
    print(f"Created synthetic data with {adata.n_obs} spots and {adata.n_vars} genes")
    print(f"  - {sum(adata.var['gene_type'] == 'ligand')} ligand genes")
    print(f"  - {sum(adata.var['gene_type'] == 'receptor')} receptor genes")
    print(f"  - 4 communication regions defined")
    
    return adata


async def demonstrate_cell_communication():
    """Demonstrate cell communication analysis"""
    print("=== ChatSpatial Cell Communication Analysis Demo ===\n")
    
    # Create synthetic data
    adata = create_synthetic_communication_data(n_spots=300, n_genes=120)
    
    # Create mock data store
    data_store = {
        "demo_data": {
            "adata": adata.copy(),
            "name": "Synthetic Communication Data",
            "type": "synthetic",
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars
        }
    }
    
    print("1. Preprocessing data...")
    
    # Mock context for logging
    class MockContext:
        async def info(self, message):
            print(f"   INFO: {message}")
        async def warning(self, message):
            print(f"   WARNING: {message}")
    
    mock_context = MockContext()
    
    # Basic preprocessing
    preprocess_params = AnalysisParameters(
        normalization="log",
        scale=False,  # Don't scale for communication analysis
        n_hvgs=None,  # Use all genes
        n_pcs=20
    )
    
    await preprocess_data("demo_data", data_store, preprocess_params, mock_context)
    
    print("\n2. Analyzing cell communication...")
    
    # Test COMMOT analysis
    print("\n--- COMMOT Analysis ---")
    commot_params = CellCommunicationParameters(
        method="commot",
        species="human",
        database="cellchat",
        commot_dis_thr=150.0,  # Distance threshold
        commot_heteromeric=True,
        min_cells=3,
        perform_global_analysis=True,
        perform_local_analysis=True,
        identify_communication_patterns=True,
        include_image=False  # Skip visualization for this demo
    )
    
    try:
        commot_result = await analyze_cell_communication(
            "demo_data", data_store, commot_params, mock_context
        )
        
        print(f"   Analyzed {commot_result.n_lr_pairs} ligand-receptor pairs")
        print(f"   Found {commot_result.n_significant_pairs} significant communication pairs")
        
        if commot_result.top_lr_pairs:
            print(f"   Top 5 communication pairs:")
            for i, pair in enumerate(commot_result.top_lr_pairs[:5]):
                print(f"     {i+1}. {pair}")
        
        if commot_result.patterns_identified:
            print(f"   Identified {commot_result.n_patterns} communication patterns")
        
        # Show some statistics
        stats = commot_result.statistics
        print(f"   Communication statistics:")
        print(f"     Distance threshold: {stats.get('distance_threshold', 'N/A')}")
        print(f"     Communicating spots: {stats.get('n_communicating_spots', 'N/A')}")
        
    except ImportError as e:
        print(f"   COMMOT not available: {e}")
        print("   Install with: pip install commot")
    except Exception as e:
        print(f"   COMMOT failed: {e}")
    
    # Test SpatialDM analysis
    print("\n--- SpatialDM Analysis ---")
    spatialdm_params = CellCommunicationParameters(
        method="spatialdm",
        species="human",
        min_cells=3,
        spatialdm_cutoff=0.1,
        spatialdm_n_permutations=100,  # Reduced for demo
        spatialdm_method="z-score",
        spatialdm_fdr=True,
        spatialdm_threshold=0.1,
        perform_global_analysis=True,
        perform_local_analysis=True,
        include_image=False
    )
    
    try:
        spatialdm_result = await analyze_cell_communication(
            "demo_data", data_store, spatialdm_params, mock_context
        )
        
        print(f"   Analyzed {spatialdm_result.n_lr_pairs} ligand-receptor pairs")
        print(f"   Found {spatialdm_result.n_significant_pairs} significant communication pairs")
        
        if spatialdm_result.top_lr_pairs:
            print(f"   Top 5 communication pairs:")
            for i, pair in enumerate(spatialdm_result.top_lr_pairs[:5]):
                print(f"     {i+1}. {pair}")
        
        if spatialdm_result.local_analysis_performed:
            print(f"   Local analysis completed")
        
        # Show some statistics
        stats = spatialdm_result.statistics
        print(f"   SpatialDM statistics:")
        print(f"     Weight cutoff: {stats.get('weight_cutoff', 'N/A')}")
        print(f"     Statistical method: {stats.get('statistical_method', 'N/A')}")
        print(f"     FDR correction: {stats.get('fdr_correction', 'N/A')}")
        
    except ImportError as e:
        print(f"   SpatialDM not available: {e}")
        print("   Install with: pip install SpatialDM")
    except Exception as e:
        print(f"   SpatialDM failed: {e}")
    
    print("\n3. Creating visualizations...")
    
    # Create visualization of the synthetic data
    final_adata = data_store["demo_data"]["adata"]
    coords = final_adata.obsm['spatial']
    
    # Plot 1: Communication regions
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Communication regions
    ax = axes[0, 0]
    regions = final_adata.obs['communication_region']
    unique_regions = regions.unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(unique_regions)))
    
    for i, region in enumerate(unique_regions):
        mask = regions == region
        ax.scatter(
            coords[mask, 0], 
            coords[mask, 1], 
            c=[colors[i]], 
            label=region,
            s=20,
            alpha=0.8
        )
    
    ax.set_title('Communication Regions')
    ax.set_xlabel('Spatial X')
    ax.set_ylabel('Spatial Y')
    ax.legend()
    ax.invert_yaxis()
    
    # Plot some key ligands and receptors
    key_genes = ['TGFB1', 'VEGFA', 'TNF', 'CXCL12']
    available_genes = [g for g in key_genes if g in final_adata.var.index]
    
    for i, gene in enumerate(available_genes[:3]):
        ax = axes[0, 1] if i == 0 else axes[1, 0] if i == 1 else axes[1, 1]
        
        gene_idx = final_adata.var.index.get_loc(gene)
        if hasattr(final_adata.X, 'toarray'):
            expr_values = final_adata.X[:, gene_idx].toarray().flatten()
        else:
            expr_values = final_adata.X[:, gene_idx]
        
        scatter = ax.scatter(
            coords[:, 0], 
            coords[:, 1], 
            c=expr_values,
            cmap='Reds',
            s=20,
            alpha=0.8
        )
        
        ax.set_title(f'{gene} Expression')
        ax.set_xlabel('Spatial X')
        ax.set_ylabel('Spatial Y')
        ax.invert_yaxis()
        plt.colorbar(scatter, ax=ax, shrink=0.8)
    
    plt.suptitle('Synthetic Cell Communication Data', fontsize=14)
    plt.tight_layout()
    
    # Save the plot
    output_path = Path("cell_communication_demo.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"   Visualization saved to: {output_path}")
    
    plt.show()
    
    print("\n4. Summary of results...")
    
    if 'commot_result' in locals():
        print(f"   COMMOT identified {commot_result.n_significant_pairs} significant communication pairs")
        if commot_result.commot_sender_key:
            print(f"   Sender signals stored in adata.obsm['{commot_result.commot_sender_key}']")
        if commot_result.commot_receiver_key:
            print(f"   Receiver signals stored in adata.obsm['{commot_result.commot_receiver_key}']")
    
    if 'spatialdm_result' in locals():
        print(f"   SpatialDM identified {spatialdm_result.n_significant_pairs} significant communication pairs")
        if spatialdm_result.global_results_key:
            print(f"   Global results stored in adata.uns['{spatialdm_result.global_results_key}']")
    
    print("\n=== Demo completed successfully! ===")
    print("\nNext steps:")
    print("1. Try with real spatial transcriptomics data")
    print("2. Install communication analysis packages:")
    print("   pip install commot SpatialDM")
    print("3. Experiment with different distance thresholds and parameters")
    print("4. Combine with spatial domain identification for region-specific communication")
    print("5. Use custom ligand-receptor databases for specific biological contexts")


if __name__ == "__main__":
    import asyncio
    asyncio.run(demonstrate_cell_communication())
