"""
Example script demonstrating spatial variable genes identification functionality in ChatSpatial.

This script shows how to:
1. Load spatial transcriptomics data
2. Preprocess the data
3. Identify spatial variable genes using SpatialDE
4. Perform automatic expression histology (AEH)
5. Visualize the results
"""

import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Import ChatSpatial modules
from chatspatial.models.data import SpatialVariableGenesParameters
from chatspatial.tools.spatial_variable_genes import identify_spatial_variable_genes
from chatspatial.tools.preprocessing import preprocess_data
from chatspatial.models.data import AnalysisParameters


def create_synthetic_spatial_data_with_patterns(n_spots=400, n_genes=200):
    """Create synthetic spatial transcriptomics data with known spatial patterns"""
    print("Creating synthetic spatial transcriptomics data with spatial patterns...")
    
    np.random.seed(42)
    
    # Create spatial coordinates
    x_coords = np.random.uniform(0, 20, n_spots)
    y_coords = np.random.uniform(0, 20, n_spots)
    coords = np.column_stack([x_coords, y_coords])
    
    # Create base expression matrix
    base_expression = np.random.negative_binomial(10, 0.4, size=(n_spots, n_genes))
    
    # Add spatial patterns to specific genes
    spatial_genes = []
    
    # Pattern 1: Central hotspot (genes 0-9)
    center_x, center_y = 10, 10
    for i in range(10):
        distance_from_center = np.sqrt((x_coords - center_x)**2 + (y_coords - center_y)**2)
        spatial_effect = np.exp(-distance_from_center / 3)
        base_expression[:, i] = base_expression[:, i] * (1 + 3 * spatial_effect)
        spatial_genes.append(f"Gene_{i:03d}_central")
    
    # Pattern 2: Left-right gradient (genes 10-19)
    for i in range(10, 20):
        gradient_effect = (x_coords - x_coords.min()) / (x_coords.max() - x_coords.min())
        base_expression[:, i] = base_expression[:, i] * (1 + 2 * gradient_effect)
        spatial_genes.append(f"Gene_{i:03d}_gradient")
    
    # Pattern 3: Corner clusters (genes 20-29)
    corners = [(2, 2), (18, 2), (2, 18), (18, 18)]
    for i in range(20, 30):
        corner_idx = (i - 20) % 4
        corner_x, corner_y = corners[corner_idx]
        distance_from_corner = np.sqrt((x_coords - corner_x)**2 + (y_coords - corner_y)**2)
        spatial_effect = np.exp(-distance_from_corner / 2)
        base_expression[:, i] = base_expression[:, i] * (1 + 2 * spatial_effect)
        spatial_genes.append(f"Gene_{i:03d}_corner_{corner_idx}")
    
    # Pattern 4: Ring pattern (genes 30-39)
    for i in range(30, 40):
        distance_from_center = np.sqrt((x_coords - center_x)**2 + (y_coords - center_y)**2)
        ring_effect = np.exp(-np.abs(distance_from_center - 6) / 1.5)
        base_expression[:, i] = base_expression[:, i] * (1 + 2 * ring_effect)
        spatial_genes.append(f"Gene_{i:03d}_ring")
    
    # Remaining genes are non-spatial (random noise)
    non_spatial_genes = [f"Gene_{i:03d}_random" for i in range(40, n_genes)]
    
    # Create AnnData object
    adata = sc.AnnData(base_expression.astype(float))
    adata.obsm['spatial'] = coords
    
    # Add gene names
    all_gene_names = spatial_genes + non_spatial_genes
    adata.var_names = all_gene_names[:n_genes]
    
    # Add spot names
    adata.obs_names = [f"Spot_{i:04d}" for i in range(n_spots)]
    
    # Add metadata about which genes are truly spatial
    adata.var['is_spatial'] = [True] * 40 + [False] * (n_genes - 40)
    adata.var['spatial_pattern'] = (['central'] * 10 + ['gradient'] * 10 + 
                                   ['corner'] * 10 + ['ring'] * 10 + 
                                   ['none'] * (n_genes - 40))
    
    print(f"Created synthetic data with {adata.n_obs} spots and {adata.n_vars} genes")
    print(f"  - 40 genes with known spatial patterns")
    print(f"  - {n_genes - 40} genes with random expression")
    
    return adata


async def demonstrate_spatial_variable_genes():
    """Demonstrate spatial variable genes identification"""
    print("=== ChatSpatial Spatial Variable Genes Demo ===\n")
    
    # Create synthetic data
    adata = create_synthetic_spatial_data_with_patterns(n_spots=300, n_genes=150)
    
    # Create mock data store
    data_store = {
        "demo_data": {
            "adata": adata.copy(),
            "name": "Synthetic Spatial Data with Patterns",
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
        scale=False,  # Don't scale for SpatialDE
        n_hvgs=None,  # Use all genes
        n_pcs=20
    )
    
    await preprocess_data("demo_data", data_store, preprocess_params, mock_context)
    
    print("\n2. Identifying spatial variable genes using SpatialDE...")
    
    # Test SpatialDE without AEH first
    print("\n--- SpatialDE Analysis ---")
    spatialde_params = SpatialVariableGenesParameters(
        method="spatialde",
        n_top_genes=20,
        significance_threshold=0.05,
        spatialde_normalize=True,
        spatialde_regress_out_total_counts=True,
        perform_aeh=False,
        include_image=False  # Skip visualization for this demo
    )
    
    try:
        spatialde_result = await identify_spatial_variable_genes(
            "demo_data", data_store, spatialde_params, mock_context
        )
        
        print(f"   Identified {spatialde_result.n_significant_genes} significant spatial variable genes")
        print(f"   Tested {spatialde_result.n_tested_genes} genes total")
        print(f"   Top 10 spatial variable genes:")
        for i, gene in enumerate(spatialde_result.top_genes[:10]):
            print(f"     {i+1}. {gene}")
        
        # Check accuracy against ground truth
        final_adata = data_store["demo_data"]["adata"]
        true_spatial_genes = final_adata.var[final_adata.var['is_spatial']].index.tolist()
        detected_spatial_genes = spatialde_result.top_genes[:40]  # Top 40 to match true number
        
        # Calculate precision and recall
        true_positives = len(set(true_spatial_genes) & set(detected_spatial_genes))
        precision = true_positives / len(detected_spatial_genes) if detected_spatial_genes else 0
        recall = true_positives / len(true_spatial_genes) if true_spatial_genes else 0
        
        print(f"   Accuracy metrics:")
        print(f"     True spatial genes: {len(true_spatial_genes)}")
        print(f"     Detected in top 40: {true_positives}")
        print(f"     Precision: {precision:.3f}")
        print(f"     Recall: {recall:.3f}")
        
    except ImportError as e:
        print(f"   SpatialDE not available: {e}")
        print("   Install with: pip install spatialde NaiveDE")
        return
    except Exception as e:
        print(f"   SpatialDE failed: {e}")
        return
    
    print("\n3. Performing Automatic Expression Histology (AEH)...")
    
    # Test SpatialDE with AEH
    aeh_params = SpatialVariableGenesParameters(
        method="spatialde",
        n_top_genes=20,
        significance_threshold=0.1,  # More lenient for AEH
        spatialde_normalize=True,
        spatialde_regress_out_total_counts=True,
        perform_aeh=True,
        aeh_n_patterns=4,  # We know there are 4 spatial patterns
        aeh_length_scale=None,  # Auto-determine
        include_image=False
    )
    
    try:
        aeh_result = await identify_spatial_variable_genes(
            "demo_data", data_store, aeh_params, mock_context
        )
        
        if aeh_result.aeh_performed:
            print(f"   AEH successfully identified {aeh_result.n_patterns} spatial patterns")
            
            # Analyze pattern membership
            final_adata = data_store["demo_data"]["adata"]
            membership_cols = [col for col in final_adata.var.columns 
                             if aeh_result.aeh_membership_key and aeh_result.aeh_membership_key in col]
            
            if membership_cols:
                pattern_col = [col for col in membership_cols if 'pattern' in col][0]
                patterns = final_adata.var[pattern_col].dropna()
                
                print(f"   Pattern distribution:")
                for pattern_id in sorted(patterns.unique()):
                    genes_in_pattern = patterns[patterns == pattern_id].index.tolist()
                    print(f"     Pattern {int(pattern_id)}: {len(genes_in_pattern)} genes")
                    
                    # Show some example genes
                    example_genes = genes_in_pattern[:5]
                    print(f"       Examples: {', '.join(example_genes)}")
        else:
            print("   AEH was not performed (possibly no significant genes)")
            
    except Exception as e:
        print(f"   AEH failed: {e}")
    
    print("\n4. Creating visualizations...")
    
    # Create visualization of results
    final_adata = data_store["demo_data"]["adata"]
    coords = final_adata.obsm['spatial']
    
    # Get top spatial genes for visualization
    top_genes = spatialde_result.top_genes[:6] if 'spatialde_result' in locals() else []
    
    if top_genes:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        for i, gene in enumerate(top_genes):
            ax = axes[i]
            
            if gene in final_adata.var.index:
                # Get expression values
                gene_idx = final_adata.var.index.get_loc(gene)
                if hasattr(final_adata.X, 'toarray'):
                    expr_values = final_adata.X[:, gene_idx].toarray().flatten()
                else:
                    expr_values = final_adata.X[:, gene_idx]
                
                # Create scatter plot
                scatter = ax.scatter(
                    coords[:, 0], 
                    coords[:, 1], 
                    c=expr_values,
                    cmap='viridis',
                    s=20,
                    alpha=0.8
                )
                
                # Check if this gene has a known spatial pattern
                if gene in final_adata.var.index:
                    pattern = final_adata.var.loc[gene, 'spatial_pattern']
                    is_spatial = final_adata.var.loc[gene, 'is_spatial']
                    title = f'{gene}\n({pattern}, {"spatial" if is_spatial else "random"})'
                else:
                    title = gene
                
                ax.set_title(title, fontsize=10)
                ax.set_xlabel('Spatial X')
                ax.set_ylabel('Spatial Y')
                ax.set_aspect('equal')
                
                # Add colorbar
                plt.colorbar(scatter, ax=ax, shrink=0.8)
            else:
                ax.text(0.5, 0.5, f'Gene {gene}\nnot found', 
                       ha='center', va='center', transform=ax.transAxes)
        
        plt.suptitle('Top 6 Spatial Variable Genes', fontsize=14)
        plt.tight_layout()
        
        # Save the plot
        output_path = Path("spatial_variable_genes_demo.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"   Visualization saved to: {output_path}")
        
        plt.show()
    
    print("\n5. Summary of results...")
    
    if 'spatialde_result' in locals():
        print(f"   SpatialDE identified {spatialde_result.n_significant_genes} significant spatial genes")
        print(f"   Top spatial gene: {spatialde_result.top_genes[0] if spatialde_result.top_genes else 'None'}")
        
        # Show some statistics
        stats = spatialde_result.statistics
        if stats.get('min_qval'):
            print(f"   Minimum q-value: {stats['min_qval']:.2e}")
    
    if 'aeh_result' in locals() and aeh_result.aeh_performed:
        print(f"   AEH identified {aeh_result.n_patterns} spatial expression patterns")
    
    print("\n=== Demo completed successfully! ===")
    print("\nNext steps:")
    print("1. Try with real spatial transcriptomics data")
    print("2. Install SpatialDE for full functionality:")
    print("   pip install spatialde NaiveDE")
    print("3. Experiment with different significance thresholds")
    print("4. Use AEH to discover spatial expression patterns")
    print("5. Integrate results with spatial domain identification")


if __name__ == "__main__":
    import asyncio
    asyncio.run(demonstrate_spatial_variable_genes())
