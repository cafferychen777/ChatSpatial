#!/usr/bin/env python3
"""
Data preparation script for comprehensive visualization testing.

Downloads and prepares various spatial transcriptomics datasets for testing:
- 10X Visium datasets
- SlideSeq datasets  
- MERFISH datasets
- seqFISH datasets
- Stereo-seq datasets

This ensures comprehensive testing across different spatial technologies.
"""

import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
import requests
import gzip
import tarfile
import zipfile
from typing import Dict, Any, Optional, List
import warnings

# Add chatspatial to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

class SpatialDatasetDownloader:
    """Downloads and prepares spatial datasets for testing"""
    
    def __init__(self, base_dir=None):
        self.base_dir = base_dir or Path(__file__).parent.parent.parent / "data" / "test_datasets"
        self.base_dir.mkdir(parents=True, exist_ok=True)
        
    def download_file(self, url, filepath, desc=""):
        """Download file with progress indication"""
        print("Downloading {}: {}".format(desc, url))
        
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            with open(filepath, 'wb') as f:
                downloaded = 0
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total_size > 0:
                            percent = (downloaded / total_size) * 100
                            print(f"\r  Progress: {percent:.1f}%", end="", flush=True)
            
            print("\n  ✅ Download complete")
            return True
            
        except Exception as e:
            print(f"\n  ❌ Download failed: {e}")
            return False
    
    def create_enhanced_visium_dataset(self):
        """Create an enhanced Visium dataset with multiple analysis layers"""
        print("\n=== Creating Enhanced Visium Dataset ===")
        
        try:
            # Create a realistic Visium dataset with 2000 spots, 3000 genes
            n_spots = 2000
            n_genes = 3000
            
            # Generate hexagonal spatial coordinates
            coords = self._generate_hexagonal_coordinates(n_spots, spacing=1.0)
            
            # Generate realistic gene expression with spatial patterns
            X = self._generate_spatial_expression(coords, n_genes)
            
            # Create AnnData object
            adata = sc.AnnData(X=X)
            adata.obsm['spatial'] = coords
            adata.var_names = [f'Gene_{i:04d}' for i in range(n_genes)]
            adata.var_names_make_unique()
            
            # Add Visium-specific metadata
            adata.obs['in_tissue'] = np.random.choice([True, False], n_spots, p=[0.85, 0.15])
            adata.obs['array_row'] = np.random.randint(0, 78, n_spots)
            adata.obs['array_col'] = np.random.randint(0, 128, n_spots)
            
            # Add spatial metadata
            adata.uns['spatial'] = {
                'enhanced_visium': {
                    'images': {},
                    'scalefactors': {
                        'tissue_hires_scalef': 1.0,
                        'tissue_lowres_scalef': 0.5,
                        'spot_diameter_fullres': 89.43
                    }
                }
            }
            
            # Add comprehensive analysis
            self._add_comprehensive_analysis(adata)
            
            # Save dataset
            output_path = self.base_dir / "enhanced_visium.h5ad"
            adata.write(output_path)
            print(f"✅ Enhanced Visium dataset saved: {output_path}")
            print(f"   Shape: {adata.shape}, Spatial coords: {adata.obsm['spatial'].shape}")
            
            return True
            
        except Exception as e:
            print(f"❌ Failed to create enhanced Visium dataset: {e}")
            return False
    
    def create_high_resolution_merfish_dataset(self):
        """Create a high-resolution MERFISH-like dataset"""
        print("\n=== Creating High-Resolution MERFISH Dataset ===")
        
        try:
            # MERFISH typically has many cells (5000+) but fewer genes (100-500)
            n_cells = 5000
            n_genes = 300
            
            # Generate high-resolution coordinates (subcellular resolution)
            coords = np.random.uniform(0, 1000, size=(n_cells, 2))
            
            # Add spatial structure with distinct tissue regions
            coords = self._add_tissue_regions(coords, n_regions=8)
            
            # Generate sparse expression data (MERFISH is typically sparse)
            X = self._generate_merfish_expression(coords, n_genes)
            
            # Create AnnData object
            adata = sc.AnnData(X=X)
            adata.obsm['spatial'] = coords
            adata.var_names = [f'MERFISH_Gene_{i:03d}' for i in range(n_genes)]
            adata.var_names_make_unique()
            
            # Add MERFISH-specific cell type annotations
            cell_types = ['Excitatory_Neuron', 'Inhibitory_Neuron', 'Astrocyte', 
                         'Microglia', 'Oligodendrocyte', 'Endothelial', 'Pericyte']
            adata.obs['cell_type'] = pd.Categorical(np.random.choice(cell_types, n_cells))
            
            # Add subcellular compartment information (unique to MERFISH)
            compartments = ['Cytoplasm', 'Nucleus', 'Membrane']
            adata.obs['subcellular_compartment'] = pd.Categorical(
                np.random.choice(compartments, n_cells, p=[0.6, 0.3, 0.1])
            )
            
            # Add comprehensive analysis
            self._add_comprehensive_analysis(adata)
            
            # Save dataset
            output_path = self.base_dir / "high_res_merfish.h5ad"
            adata.write(output_path)
            print(f"✅ High-resolution MERFISH dataset saved: {output_path}")
            print(f"   Shape: {adata.shape}, Cell types: {len(adata.obs['cell_type'].unique())}")
            
            return True
            
        except Exception as e:
            print(f"❌ Failed to create MERFISH dataset: {e}")
            return False
    
    def create_slideseq_bead_dataset(self):
        """Create a SlideSeq bead-based dataset"""
        print("\n=== Creating SlideSeq Bead Dataset ===")
        
        try:
            # SlideSeq has many small beads
            n_beads = 8000
            n_genes = 1500
            
            # Generate random bead coordinates
            coords = np.random.uniform(0, 200, size=(n_beads, 2))
            
            # Generate highly sparse expression (SlideSeq characteristic)
            X = self._generate_slideseq_expression(coords, n_genes)
            
            # Create AnnData object
            adata = sc.AnnData(X=X)
            adata.obsm['spatial'] = coords
            adata.var_names = [f'SlideSeq_Gene_{i:04d}' for i in range(n_genes)]
            adata.var_names_make_unique()
            
            # Add SlideSeq-specific metadata
            adata.obs['bead_id'] = [f'bead_{i:06d}' for i in range(n_beads)]
            adata.obs['bead_type'] = pd.Categorical(
                np.random.choice(['Standard', 'High_Capture'], n_beads, p=[0.9, 0.1])
            )
            
            # Add UMI and gene counts
            adata.obs['n_counts'] = np.array(X.sum(axis=1)).flatten()
            adata.obs['n_genes'] = np.array((X > 0).sum(axis=1)).flatten()
            
            # Add comprehensive analysis
            self._add_comprehensive_analysis(adata)
            
            # Save dataset
            output_path = self.base_dir / "slideseq_beads.h5ad"
            adata.write(output_path)
            print(f"✅ SlideSeq bead dataset saved: {output_path}")
            print(f"   Shape: {adata.shape}, Sparsity: {(X == 0).sum() / X.size * 100:.1f}%")
            
            return True
            
        except Exception as e:
            print(f"❌ Failed to create SlideSeq dataset: {e}")
            return False
    
    def create_seqfish_dataset(self):
        """Create a seqFISH+ dataset with temporal dynamics"""
        print("\n=== Creating seqFISH+ Dataset ===")
        
        try:
            # seqFISH typically has moderate number of cells and genes
            n_cells = 1200
            n_genes = 400
            
            # Generate 3D coordinates (seqFISH can be 3D)
            coords_2d = np.random.uniform(0, 100, size=(n_cells, 2))
            z_coords = np.random.uniform(0, 20, size=(n_cells, 1))
            coords_3d = np.column_stack([coords_2d, z_coords])
            
            # Generate expression with developmental patterns
            X = self._generate_developmental_expression(coords_2d, n_genes)
            
            # Create AnnData object
            adata = sc.AnnData(X=X)
            adata.obsm['spatial'] = coords_2d  # Standard spatial for visualization
            adata.obsm['spatial_3d'] = coords_3d  # Full 3D coordinates
            adata.var_names = [f'seqFISH_Gene_{i:03d}' for i in range(n_genes)]
            adata.var_names_make_unique()
            
            # Add developmental time points
            adata.obs['timepoint'] = pd.Categorical(
                np.random.choice(['E8.5', 'E9.5', 'E10.5', 'E11.5'], n_cells)
            )
            
            # Add developmental stage
            adata.obs['dev_stage'] = pd.Categorical(
                np.random.choice(['Early', 'Mid', 'Late'], n_cells)
            )
            
            # Add z-level information
            adata.obs['z_level'] = z_coords.flatten()
            
            # Add comprehensive analysis
            self._add_comprehensive_analysis(adata)
            
            # Save dataset
            output_path = self.base_dir / "seqfish_developmental.h5ad"
            adata.write(output_path)
            print(f"✅ seqFISH+ dataset saved: {output_path}")
            print(f"   Shape: {adata.shape}, Timepoints: {len(adata.obs['timepoint'].unique())}")
            
            return True
            
        except Exception as e:
            print(f"❌ Failed to create seqFISH dataset: {e}")
            return False
    
    def create_stereo_seq_dataset(self):
        """Create a Stereo-seq dataset (high-throughput spatial)"""
        print("\n=== Creating Stereo-seq Dataset ===")
        
        try:
            # Stereo-seq can have very high resolution
            n_spots = 15000
            n_genes = 4000
            
            # Generate grid-like coordinates (Stereo-seq uses DNA nanoball arrays)
            coords = self._generate_grid_coordinates(n_spots, grid_size=200)
            
            # Generate expression with high dynamic range
            X = self._generate_high_throughput_expression(coords, n_genes)
            
            # Create AnnData object
            adata = sc.AnnData(X=X)
            adata.obsm['spatial'] = coords
            adata.var_names = [f'Stereo_Gene_{i:04d}' for i in range(n_genes)]
            adata.var_names_make_unique()
            
            # Add Stereo-seq specific metadata
            adata.obs['bin_size'] = pd.Categorical(['500nm'] * n_spots)  # Bin resolution
            adata.obs['chip_id'] = 'CHIP001'
            
            # Add tissue region annotations
            tissue_regions = ['Cortex', 'Hippocampus', 'Thalamus', 'Cerebellum', 'Brainstem']
            adata.obs['tissue_region'] = pd.Categorical(np.random.choice(tissue_regions, n_spots))
            
            # Add comprehensive analysis
            self._add_comprehensive_analysis(adata)
            
            # Save dataset
            output_path = self.base_dir / "stereo_seq_brain.h5ad"
            adata.write(output_path)
            print(f"✅ Stereo-seq dataset saved: {output_path}")
            print(f"   Shape: {adata.shape}, Resolution: 500nm")
            
            return True
            
        except Exception as e:
            print(f"❌ Failed to create Stereo-seq dataset: {e}")
            return False
    
    def _generate_hexagonal_coordinates(self, n_spots, spacing=1.0):
        """Generate hexagonal spot coordinates (Visium-like)"""
        coords = []
        rows = int(np.sqrt(n_spots * 2 / np.sqrt(3))) + 1
        
        for row in range(rows):
            cols_in_row = int(2 * n_spots / (rows * np.sqrt(3))) + 1
            for col in range(cols_in_row):
                if len(coords) >= n_spots:
                    break
                    
                x = col * spacing
                if row % 2 == 1:
                    x += spacing * 0.5
                y = row * spacing * np.sqrt(3) / 2
                
                coords.append([x, y])
        
        coords = np.array(coords[:n_spots])
        
        # Center coordinates
        coords -= coords.mean(axis=0)
        coords *= 100  # Scale to reasonable range
        
        return coords
    
    def _generate_grid_coordinates(self, n_spots, grid_size) :
        """Generate grid coordinates for Stereo-seq like data"""
        coords = []
        grid_dim = int(np.sqrt(n_spots)) + 1
        
        for i in range(grid_dim):
            for j in range(grid_dim):
                if len(coords) >= n_spots:
                    break
                coords.append([i * grid_size / grid_dim, j * grid_size / grid_dim])
        
        return np.array(coords[:n_spots])
    
    def _generate_spatial_expression(self, coords, n_genes) :
        """Generate spatially structured gene expression"""
        n_spots = coords.shape[0]
        X = np.random.negative_binomial(5, 0.3, size=(n_spots, n_genes))
        
        # Add spatial patterns to subset of genes
        center = coords.mean(axis=0)
        distances = np.linalg.norm(coords - center, axis=1)
        max_dist = distances.max()
        
        for i in range(min(50, n_genes)):  # Add patterns to first 50 genes
            if i % 5 == 0:  # Radial pattern
                pattern = np.exp(-distances / (max_dist / 3))
            elif i % 5 == 1:  # Gradient pattern
                pattern = (coords[:, 0] - coords[:, 0].min()) / (coords[:, 0].max() - coords[:, 0].min())
            elif i % 5 == 2:  # Stripe pattern  
                pattern = np.sin(coords[:, 1] / max_dist * 4 * np.pi) + 1
            elif i % 5 == 3:  # Cluster pattern
                pattern = np.exp(-((coords[:, 0] - center[0] + np.random.normal(0, max_dist/4))**2 + 
                                 (coords[:, 1] - center[1] + np.random.normal(0, max_dist/4))**2) / (max_dist/2)**2)
            else:  # Random pattern
                pattern = np.ones(n_spots)
            
            X[:, i] = (X[:, i] * (0.1 + 2 * pattern)).astype(int)
        
        return X
    
    def _generate_merfish_expression(self, coords, n_genes) :
        """Generate sparse MERFISH-like expression"""
        n_cells = coords.shape[0]
        
        # MERFISH is typically very sparse
        X = np.random.negative_binomial(1, 0.8, size=(n_cells, n_genes))
        
        # Make it more sparse
        sparsity_mask = np.random.random((n_cells, n_genes)) < 0.05  # Only 5% non-zero
        X[~sparsity_mask] = 0
        
        return X
    
    def _generate_slideseq_expression(self, coords, n_genes) :
        """Generate SlideSeq-like sparse expression"""
        n_beads = coords.shape[0]
        
        # SlideSeq is sparse but less than MERFISH
        X = np.random.negative_binomial(2, 0.6, size=(n_beads, n_genes))
        
        # Apply sparsity
        sparsity_mask = np.random.random((n_beads, n_genes)) < 0.15  # 15% non-zero
        X[~sparsity_mask] = 0
        
        return X
    
    def _generate_developmental_expression(self, coords, n_genes) :
        """Generate expression with developmental patterns"""
        n_cells = coords.shape[0]
        X = np.random.negative_binomial(3, 0.4, size=(n_cells, n_genes))
        
        # Add developmental gradients
        for i in range(min(20, n_genes)):
            if i % 4 == 0:  # Anterior-posterior gradient
                gradient = coords[:, 0] / coords[:, 0].max()
            elif i % 4 == 1:  # Dorsal-ventral gradient
                gradient = coords[:, 1] / coords[:, 1].max()
            elif i % 4 == 2:  # Central-peripheral gradient
                center = coords.mean(axis=0)
                distances = np.linalg.norm(coords - center, axis=1)
                gradient = 1 - distances / distances.max()
            else:  # Wave-like pattern
                gradient = np.sin(coords[:, 0] / coords[:, 0].max() * 2 * np.pi) + 1
            
            X[:, i] = (X[:, i] * (0.1 + 1.9 * gradient)).astype(int)
        
        return X
    
    def _generate_high_throughput_expression(self, coords, n_genes) :
        """Generate high-throughput expression with good dynamic range"""
        n_spots = coords.shape[0]
        
        # Higher counts for high-throughput methods
        X = np.random.negative_binomial(8, 0.2, size=(n_spots, n_genes))
        
        return X
    
    def _add_tissue_regions(self, coords, n_regions = 5) :
        """Add tissue region structure to coordinates"""
        from sklearn.cluster import KMeans
        
        # Cluster coordinates into regions
        kmeans = KMeans(n_clusters=n_regions, random_state=42)
        region_labels = kmeans.fit_predict(coords)
        
        # Add some structure by moving points towards region centers
        for region in range(n_regions):
            mask = region_labels == region
            if mask.sum() > 0:
                center = coords[mask].mean(axis=0)
                coords[mask] = 0.7 * coords[mask] + 0.3 * center
        
        return coords
    
    def _add_comprehensive_analysis(self, adata):
        """Add comprehensive analysis results for testing visualization functions"""
        print("  Adding comprehensive analysis layers...")
        
        try:
            # Basic QC metrics
            adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
            adata.obs['n_genes'] = np.array((adata.X > 0).sum(axis=1)).flatten()
            
            # Preprocessing
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            
            # Feature selection
            n_hvg = min(200, adata.n_vars // 3)
            sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=False)
            
            # PCA
            n_pcs = min(30, adata.n_obs - 1, adata.n_vars - 1)
            sc.tl.pca(adata, n_comps=n_pcs)
            
            # Neighbors and clustering
            sc.pp.neighbors(adata, n_neighbors=min(20, adata.n_obs - 1))
            sc.tl.leiden(adata, resolution=0.5)
            sc.tl.umap(adata)
            
            # Additional clustering at different resolutions
            sc.tl.leiden(adata, resolution=0.2, key_added='leiden_02')
            sc.tl.leiden(adata, resolution=1.0, key_added='leiden_10')
            
            # Add mock deconvolution results (multiple methods)
            cell_types = ['Neuron', 'Astrocyte', 'Microglia', 'Oligodendrocyte', 'Endothelial', 'Pericyte']
            
            for method in ['spotlight', 'cell2location', 'stereoscope']:
                # Create mock deconvolution matrix
                deconv_matrix = np.random.dirichlet(np.ones(len(cell_types)), size=adata.n_obs)
                
                # Add some structure based on spatial location if available
                if 'spatial' in adata.obsm:
                    coords = adata.obsm['spatial']
                    for i, ct in enumerate(cell_types):
                        # Create spatial patterns for each cell type
                        center_x = np.random.uniform(coords[:, 0].min(), coords[:, 0].max())
                        center_y = np.random.uniform(coords[:, 1].min(), coords[:, 1].max())
                        
                        distances = np.sqrt((coords[:, 0] - center_x)**2 + (coords[:, 1] - center_y)**2)
                        pattern = np.exp(-distances / distances.max())
                        
                        deconv_matrix[:, i] *= (0.5 + pattern)
                
                # Renormalize
                deconv_matrix = deconv_matrix / deconv_matrix.sum(axis=1, keepdims=True)
                
                # Store in obsm
                adata.obsm[f'deconvolution_{method}'] = deconv_matrix
                adata.uns[f'deconvolution_{method}_cell_types'] = cell_types
                
                # Add individual cell type columns
                for i, ct in enumerate(cell_types):
                    adata.obs[f'deconvolution_{method}_{ct}'] = deconv_matrix[:, i]
            
            # Add mock RNA velocity results
            if adata.n_obs > 100:
                velocity_genes = adata.var_names[:min(50, len(adata.var_names))]
                
                # Mock velocity vectors
                adata.layers['velocity'] = np.random.normal(0, 0.1, adata.X.shape)
                adata.var['velocity_genes'] = adata.var_names.isin(velocity_genes)
                
                # Mock velocity in UMAP space
                if 'X_umap' in adata.obsm:
                    velocity_umap = np.random.normal(0, 0.5, adata.obsm['X_umap'].shape)
                    adata.obsm['velocity_umap'] = velocity_umap
            
            # Add mock pathway enrichment results
            pathways = ['Cell_Cycle', 'Apoptosis', 'Inflammation', 'Synaptic_Transmission', 'Metabolism']
            for pathway in pathways:
                # Mock enrichment scores
                scores = np.random.normal(0, 1, adata.n_obs)
                adata.obs[f'pathway_{pathway}'] = scores
            
            # Add mock spatial analysis results
            if 'spatial' in adata.obsm and adata.n_obs > 50:
                # Mock Moran's I results
                adata.obs['morans_i'] = np.random.uniform(-0.5, 0.5, adata.n_obs)
                
                # Mock spatial domains
                n_domains = min(8, adata.n_obs // 20)
                adata.obs['spatial_domains'] = pd.Categorical(
                    [f'Domain_{i % n_domains}' for i in range(adata.n_obs)]
                )
            
            # Add mock communication results
            if adata.n_obs > 50:
                lr_pairs = ['CXCL12_CXCR4', 'TNF_TNFRSF1A', 'VEGFA_VEGFR1', 'PDGFA_PDGFRA']
                comm_matrix = np.random.exponential(0.5, size=(len(cell_types), len(cell_types), len(lr_pairs)))
                
                # Store communication results
                adata.uns['cell_communication'] = {
                    'lr_pairs': lr_pairs,
                    'cell_types': cell_types,
                    'communication_matrix': comm_matrix
                }
                
                # Add LR pair expression scores
                for lr_pair in lr_pairs:
                    adata.obs[f'lr_score_{lr_pair}'] = np.random.exponential(1, adata.n_obs)
            
            print("  ✅ Analysis layers added successfully")
            
        except Exception as e:
            print(f"  ⚠️  Some analysis layers failed: {e}")
            # Continue anyway - basic dataset should still work
    
    def run_all_preparations(self):
        """Run all dataset preparation steps"""
        print("=" * 80)
        print("SPATIAL DATASET PREPARATION")
        print("=" * 80)
        
        print(f"Output directory: {self.base_dir}")
        
        success_count = 0
        total_count = 5
        
        # Create all datasets
        if self.create_enhanced_visium_dataset():
            success_count += 1
            
        if self.create_high_resolution_merfish_dataset():
            success_count += 1
            
        if self.create_slideseq_bead_dataset():
            success_count += 1
            
        if self.create_seqfish_dataset():
            success_count += 1
            
        if self.create_stereo_seq_dataset():
            success_count += 1
        
        print("\n" + "=" * 80)
        print(f"DATASET PREPARATION SUMMARY")
        print("=" * 80)
        print(f"Successfully created: {success_count}/{total_count} datasets")
        print(f"Output directory: {self.base_dir}")
        
        # List created files
        created_files = list(self.base_dir.glob("*.h5ad"))
        print(f"\nCreated datasets:")
        for file_path in created_files:
            size_mb = file_path.stat().st_size / (1024 * 1024)
            print(f"  - {file_path.name}: {size_mb:.1f} MB")
        
        return success_count == total_count

def main():
    """Main function"""
    downloader = SpatialDatasetDownloader()
    success = downloader.run_all_preparations()
    
    if success:
        print("\n🎉 All datasets prepared successfully!")
        return 0
    else:
        print("\n⚠️  Some datasets failed to prepare")
        return 1

if __name__ == "__main__":
    sys.exit(main())