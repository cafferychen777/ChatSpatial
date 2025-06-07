"""
Example script demonstrating LIANA+ cell communication analysis in ChatSpatial.

This script shows how to:
1. Load spatial transcriptomics data
2. Preprocess the data for communication analysis
3. Analyze cell-cell communication using LIANA+
4. Visualize spatial communication patterns
5. Interpret results
"""

import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import asyncio

# Import ChatSpatial modules
from chatspatial.models.data import CellCommunicationParameters
from chatspatial.tools.cell_communication import analyze_cell_communication
from chatspatial.tools.preprocessing import preprocess_data
from chatspatial.models.data import AnalysisParameters


class MockContext:
    """Mock context for demonstration purposes"""
    async def info(self, message: str):
        print(f"INFO: {message}")
    
    async def warning(self, message: str):
        print(f"WARNING: {message}")


def create_demo_spatial_data():
    """Create synthetic spatial transcriptomics data for demonstration"""
    print("Creating synthetic spatial transcriptomics data...")
    
    # Create synthetic data with spatial structure
    n_spots = 200
    n_genes = 100
    
    # Generate spatial coordinates in a tissue-like pattern
    np.random.seed(42)
    
    # Create two main regions
    region1_spots = n_spots // 2
    region2_spots = n_spots - region1_spots
    
    # Region 1: circular pattern
    angles1 = np.random.uniform(0, 2*np.pi, region1_spots)
    radii1 = np.random.uniform(0, 3, region1_spots)
    x1 = radii1 * np.cos(angles1) + 5
    y1 = radii1 * np.sin(angles1) + 5
    
    # Region 2: circular pattern, different location
    angles2 = np.random.uniform(0, 2*np.pi, region2_spots)
    radii2 = np.random.uniform(0, 3, region2_spots)
    x2 = radii2 * np.cos(angles2) + 12
    y2 = radii2 * np.sin(angles2) + 5
    
    # Combine coordinates
    spatial_coords = np.vstack([
        np.column_stack([x1, y1]),
        np.column_stack([x2, y2])
    ])
    
    # Create gene expression data with spatial patterns
    # Include some real ligand-receptor genes
    ligand_genes = ['VEGFA', 'PDGFA', 'TGFB1', 'FGF1', 'EGF', 'IGF1', 'TNF', 'IL1A', 'CCL2', 'CXCL12']
    receptor_genes = ['KDR', 'PDGFRA', 'TGFBR1', 'FGFR1', 'EGFR', 'IGF1R', 'TNFRSF1A', 'IL1R1', 'CCR2', 'CXCR4']
    other_genes = [f'Gene_{i}' for i in range(n_genes - len(ligand_genes) - len(receptor_genes))]
    
    all_genes = ligand_genes + receptor_genes + other_genes
    
    # Generate expression data
    X = np.random.negative_binomial(5, 0.3, size=(n_spots, n_genes)).astype(float)
    
    # Add spatial patterns for ligand-receptor pairs
    for i, (ligand, receptor) in enumerate(zip(ligand_genes, receptor_genes)):
        ligand_idx = all_genes.index(ligand)
        receptor_idx = all_genes.index(receptor)
        
        # Create complementary spatial patterns
        if i % 2 == 0:
            # Ligand higher in region 1, receptor higher in region 2
            X[:region1_spots, ligand_idx] *= 3
            X[region1_spots:, receptor_idx] *= 3
        else:
            # Opposite pattern
            X[region1_spots:, ligand_idx] *= 3
            X[:region1_spots, receptor_idx] *= 3
    
    # Create AnnData object
    adata = sc.AnnData(X)
    adata.var_names = all_genes
    adata.obs_names = [f'Spot_{i}' for i in range(n_spots)]
    
    # Add spatial coordinates
    adata.obsm['spatial'] = spatial_coords
    
    # Add some metadata
    adata.obs['region'] = ['Region_1'] * region1_spots + ['Region_2'] * region2_spots
    adata.obs['region'] = adata.obs['region'].astype('category')
    
    print(f"Created synthetic data: {n_spots} spots × {n_genes} genes")
    print(f"Included {len(ligand_genes)} ligand-receptor pairs")
    
    return adata


async def run_liana_communication_analysis():
    """Run LIANA+ cell communication analysis demonstration"""
    
    print("🧬 LIANA+ Cell Communication Analysis Demo")
    print("=" * 50)
    
    # Create demo data
    adata = create_demo_spatial_data()
    
    # Create data store
    data_store = {
        "demo_data": {
            "adata": adata,
            "name": "Demo Spatial Dataset",
            "type": "synthetic"
        }
    }
    
    mock_context = MockContext()
    
    print("\n1. Preprocessing data...")
    
    # Basic preprocessing
    preprocess_params = AnalysisParameters(
        normalization="log",
        scale=False,  # Don't scale for communication analysis
        n_hvgs=None,  # Use all genes
        n_pcs=20
    )
    
    await preprocess_data("demo_data", data_store, preprocess_params, mock_context)
    
    print("\n2. Analyzing cell communication with LIANA+...")
    
    # LIANA+ spatial analysis
    print("\n--- LIANA+ Spatial Bivariate Analysis ---")
    liana_params = CellCommunicationParameters(
        method="liana",
        species="human",
        min_cells=3,
        perform_spatial_analysis=True,
        liana_local_metric="cosine",
        liana_global_metric="morans",
        liana_n_perms=50,  # Reduced for demo speed
        liana_nz_prop=0.2,
        include_image=True
    )
    
    try:
        liana_result = await analyze_cell_communication(
            "demo_data", 
            data_store, 
            liana_params, 
            mock_context
        )
        
        print(f"\n✅ LIANA+ Analysis Results:")
        print(f"   - Method: {liana_result.method}")
        print(f"   - Analysis type: {liana_result.analysis_type}")
        print(f"   - Species: {liana_result.species}")
        print(f"   - LR pairs tested: {liana_result.n_lr_pairs}")
        print(f"   - Significant pairs: {liana_result.n_significant_pairs}")
        
        if liana_result.top_lr_pairs:
            print(f"   - Top LR pairs:")
            for i, pair in enumerate(liana_result.top_lr_pairs[:5], 1):
                print(f"     {i}. {pair}")
        
        if liana_result.visualization:
            print(f"   - Visualization generated: ✅")
        
        # Display some statistics
        if liana_result.statistics:
            stats = liana_result.statistics
            print(f"\n📊 Analysis Statistics:")
            for key, value in stats.items():
                print(f"   - {key}: {value}")
        
    except Exception as e:
        print(f"❌ LIANA+ analysis failed: {e}")
        print("Make sure LIANA+ is installed: pip install liana")
        return
    
    print("\n3. Interpreting Results...")
    
    # Access the updated data
    updated_adata = data_store["demo_data"]["adata"]
    
    print(f"\n📋 Data Summary:")
    print(f"   - Original spots: {adata.n_obs}")
    print(f"   - Original genes: {adata.n_vars}")
    print(f"   - Analysis keys added to adata.uns: {list(updated_adata.uns.keys())}")
    
    if liana_result.liana_spatial_scores_key and liana_result.liana_spatial_scores_key in updated_adata.uns:
        spatial_scores = updated_adata.uns[liana_result.liana_spatial_scores_key]
        print(f"   - Spatial scores shape: {spatial_scores.shape}")
        print(f"   - Available metrics: {list(spatial_scores.var.columns)}")
    
    print("\n4. Next Steps...")
    print("   - Examine specific ligand-receptor pairs of interest")
    print("   - Correlate communication patterns with spatial domains")
    print("   - Compare different local/global metric combinations")
    print("   - Validate findings with known biological interactions")
    
    print("\n=== LIANA+ Demo completed successfully! ===")
    print("\nKey advantages of LIANA+:")
    print("✅ Fast performance (1-2 minutes vs 10-30 minutes for other methods)")
    print("✅ Spatial bivariate analysis with multiple metrics")
    print("✅ Comprehensive ligand-receptor databases")
    print("✅ Robust statistical testing")
    print("✅ Excellent visualization capabilities")


if __name__ == "__main__":
    print("🚀 Starting LIANA+ Communication Analysis Demo")
    print("=" * 50)
    
    # Check if LIANA+ is available
    try:
        import liana as li
        print(f"✅ LIANA+ is available (version: {li.__version__})")
    except ImportError:
        print("❌ LIANA+ is not installed")
        print("Please install it with: pip install liana")
        exit(1)
    
    # Run the demo
    asyncio.run(run_liana_communication_analysis())
