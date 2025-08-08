"""
Example script demonstrating spatial domain identification functionality in ChatSpatial.

This script shows how to:
1. Load spatial transcriptomics data
2. Preprocess the data
3. Identify spatial domains using different methods
4. Visualize the results
"""

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

# Import ChatSpatial modules
from chatspatial.models.data import SpatialDomainParameters
from chatspatial.tools.spatial_domains import identify_spatial_domains
from chatspatial.tools.preprocessing import preprocess_data
from chatspatial.models.data import AnalysisParameters


def create_synthetic_spatial_data(n_spots=500, n_genes=100):
    """Create synthetic spatial transcriptomics data for demonstration"""
    print("Creating synthetic spatial transcriptomics data...")
    
    # Create spatial coordinates with distinct regions
    np.random.seed(42)
    
    # Define 4 spatial regions
    regions = []
    region_centers = [(2, 2), (8, 2), (2, 8), (8, 8)]
    spots_per_region = n_spots // 4
    
    for i, (cx, cy) in enumerate(region_centers):
        # Generate spots around each center
        x_coords = np.random.normal(cx, 1.0, spots_per_region)
        y_coords = np.random.normal(cy, 1.0, spots_per_region)
        
        # Create expression profiles specific to each region
        base_expression = np.random.negative_binomial(5, 0.3, size=(spots_per_region, n_genes))
        
        # Add region-specific expression patterns
        if i == 0:  # Region 1: high expression of first 25 genes
            base_expression[:, :25] *= 3
        elif i == 1:  # Region 2: high expression of genes 25-50
            base_expression[:, 25:50] *= 3
        elif i == 2:  # Region 3: high expression of genes 50-75
            base_expression[:, 50:75] *= 3
        else:  # Region 4: high expression of last 25 genes
            base_expression[:, 75:] *= 3
        
        regions.append({
            'coords': np.column_stack([x_coords, y_coords]),
            'expression': base_expression,
            'true_labels': np.full(spots_per_region, i)
        })
    
    # Combine all regions
    all_coords = np.vstack([r['coords'] for r in regions])
    all_expression = np.vstack([r['expression'] for r in regions])
    true_labels = np.concatenate([r['true_labels'] for r in regions])
    
    # Create AnnData object
    adata = sc.AnnData(all_expression.astype(float))
    adata.obsm['spatial'] = all_coords
    adata.obs['true_domain'] = true_labels.astype(str)
    adata.obs['true_domain'] = adata.obs['true_domain'].astype('category')
    
    # Add gene and cell names
    adata.var_names = [f"Gene_{i:03d}" for i in range(n_genes)]
    adata.obs_names = [f"Spot_{i:04d}" for i in range(n_spots)]
    
    print(f"Created synthetic data with {adata.n_obs} spots and {adata.n_vars} genes")
    return adata


async def demonstrate_spatial_domain_identification():
    """Demonstrate spatial domain identification with different methods"""
    print("=== ChatSpatial Spatial Domain Identification Demo ===\n")
    
    # Create synthetic data
    adata = create_synthetic_spatial_data(n_spots=400, n_genes=80)
    
    # Create mock data store
    data_store = {
        "demo_data": {
            "adata": adata.copy(),
            "name": "Synthetic Spatial Data",
            "type": "synthetic",
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars
        }
    }
    
    print("1. Preprocessing data...")
    # Preprocess the data
    preprocess_params = AnalysisParameters(
        normalization="log",
        scale=True,
        n_hvgs=50,
        n_pcs=20
    )
    
    # Mock context for logging
    class MockContext:
        async def info(self, message):
            print(f"   INFO: {message}")
        async def warning(self, message):
            print(f"   WARNING: {message}")
    
    mock_context = MockContext()
    
    # Preprocess data
    await preprocess_data("demo_data", data_store, preprocess_params, mock_context)
    
    print("\n2. Testing different spatial domain identification methods...")
    
    # Method 1: Leiden clustering with spatial constraints
    print("\n--- Method 1: Leiden Clustering ---")
    leiden_params = SpatialDomainParameters(
        method="leiden",
        n_domains=4,
        resolution=0.5,
        refine_domains=True,
        cluster_n_neighbors=15,
        cluster_spatial_weight=0.3
    )
    
    leiden_result = await identify_spatial_domains(
        "demo_data", data_store, leiden_params, mock_context
    )
    
    print(f"   Identified {leiden_result.n_domains} domains")
    print(f"   Domain counts: {leiden_result.domain_counts}")
    print(f"   Domain labels stored in: {leiden_result.domain_key}")
    if leiden_result.refined_domain_key:
        print(f"   Refined labels stored in: {leiden_result.refined_domain_key}")
    
    # Method 2: Louvain clustering
    print("\n--- Method 2: Louvain Clustering ---")
    louvain_params = SpatialDomainParameters(
        method="louvain",
        n_domains=4,
        resolution=0.3,
        refine_domains=False,
        cluster_n_neighbors=20,
        cluster_spatial_weight=0.2
    )
    
    louvain_result = await identify_spatial_domains(
        "demo_data", data_store, louvain_params, mock_context
    )
    
    print(f"   Identified {louvain_result.n_domains} domains")
    print(f"   Domain counts: {louvain_result.domain_counts}")
    
    # Method 3: SpaGCN (if available)
    # Note: ChatSpatial supports multiple spatial domain methods when installed, e.g.,
    #   - method="spagcn" (SpaGCN)
    #   - method="leiden" / method="louvain" (clustering)
    #   - method="stagate" (STAGATE)
    #   - method="banksy" (BANKSY)
    # Install optional packages as needed: pip install SpaGCN STAGATE banksy-utils
    print("\n--- Method 3: SpaGCN ---")
    try:
        spagcn_params = SpatialDomainParameters(
            method="spagcn",
            n_domains=4,
            spagcn_s=1.0,
            spagcn_b=49,
            spagcn_use_histology=False  # No histology image for synthetic data
        )
        
        spagcn_result = await identify_spatial_domains(
            "demo_data", data_store, spagcn_params, mock_context
        )
        
        print(f"   Identified {spagcn_result.n_domains} domains")
        print(f"   Domain counts: {spagcn_result.domain_counts}")
        
    except ImportError as e:
        print(f"   SpaGCN not available: {e}")
        print("   Install with: pip install SpaGCN")
    except Exception as e:
        print(f"   SpaGCN failed: {e}")
    
    print("\n3. Comparing results with ground truth...")
    
    # Get the final adata with all results
    final_adata = data_store["demo_data"]["adata"]
    
    # Print available domain annotations
    domain_columns = [col for col in final_adata.obs.columns if 'domain' in col.lower()]
    print(f"   Available domain annotations: {domain_columns}")
    
    # Calculate basic accuracy metrics if we have results
    if leiden_result.domain_key in final_adata.obs.columns:
        true_labels = final_adata.obs['true_domain'].astype(int)
        pred_labels = final_adata.obs[leiden_result.domain_key].astype(int)
        
        # Simple accuracy calculation (note: this is not adjusted for label permutation)
        from sklearn.metrics import adjusted_rand_score
        ari_score = adjusted_rand_score(true_labels, pred_labels)
        print(f"   Leiden ARI score: {ari_score:.3f}")
    
    print("\n4. Creating visualization...")
    
    # Create a simple visualization
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    # Plot 1: True domains
    ax = axes[0]
    coords = final_adata.obsm['spatial']
    true_domains = final_adata.obs['true_domain']
    scatter = ax.scatter(coords[:, 0], coords[:, 1], c=true_domains.astype(int), 
                        cmap='tab10', s=20, alpha=0.7)
    ax.set_title('True Domains')
    ax.set_xlabel('Spatial X')
    ax.set_ylabel('Spatial Y')
    ax.invert_yaxis()
    
    # Plot 2: Leiden results
    if leiden_result.domain_key in final_adata.obs.columns:
        ax = axes[1]
        leiden_domains = final_adata.obs[leiden_result.domain_key]
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=leiden_domains.astype(int), 
                            cmap='tab10', s=20, alpha=0.7)
        ax.set_title('Leiden Domains')
        ax.set_xlabel('Spatial X')
        ax.set_ylabel('Spatial Y')
        ax.invert_yaxis()
    
    # Plot 3: Refined Leiden results (if available)
    if leiden_result.refined_domain_key and leiden_result.refined_domain_key in final_adata.obs.columns:
        ax = axes[2]
        refined_domains = final_adata.obs[leiden_result.refined_domain_key]
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=refined_domains.astype(int), 
                            cmap='tab10', s=20, alpha=0.7)
        ax.set_title('Refined Leiden Domains')
        ax.set_xlabel('Spatial X')
        ax.set_ylabel('Spatial Y')
        ax.invert_yaxis()
    
    # Plot 4: Louvain results
    if louvain_result.domain_key in final_adata.obs.columns:
        ax = axes[3]
        louvain_domains = final_adata.obs[louvain_result.domain_key]
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=louvain_domains.astype(int), 
                            cmap='tab10', s=20, alpha=0.7)
        ax.set_title('Louvain Domains')
        ax.set_xlabel('Spatial X')
        ax.set_ylabel('Spatial Y')
        ax.invert_yaxis()
    
    plt.tight_layout()
    
    # Save the plot
    output_path = Path("spatial_domains_comparison.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"   Visualization saved to: {output_path}")
    
    plt.show()
    
    print("\n=== Demo completed successfully! ===")
    print("\nNext steps:")
    print("1. Try with real spatial transcriptomics data")
    print("2. Install STAGATE and SpaGCN for advanced methods:")
    print("   pip install STAGATE SpaGCN")
    print("3. Experiment with different parameters")
    print("4. Use the identified domains for downstream analysis")


if __name__ == "__main__":
    import asyncio
    asyncio.run(demonstrate_spatial_domain_identification())
