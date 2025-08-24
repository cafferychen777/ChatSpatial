"""
Graceful degradation and fallback mechanisms for ChatSpatial.

This module provides alternative implementations when optional dependencies
are missing. The philosophy is "never break userspace" - if a user can't
install a dependency, we provide a reasonable fallback.

Design Principles:
- Always provide a working alternative 
- Clearly communicate limitations to users
- Degrade gracefully, don't fail catastrophically
- Maintain API compatibility

Author: Linus-approved fallback system
"""

import warnings
from typing import Any, Dict, List, Optional, Tuple, Union, Callable
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

from .smart_import import smart_import, get_smart_importer


class FallbackWarning(UserWarning):
    """Warning for fallback method usage"""
    pass


def warn_fallback(feature: str, missing_dep: str, fallback_desc: str):
    """Issue a user-friendly fallback warning"""
    warnings.warn(
        f"Using fallback for {feature}: {missing_dep} not available. "
        f"Falling back to {fallback_desc}. "
        f"For full functionality, install: pip install {missing_dep}",
        FallbackWarning,
        stacklevel=3
    )


class SpatialFallbacks:
    """Fallback implementations for spatial analysis methods"""
    
    @staticmethod
    def spatial_neighbors_fallback(adata: sc.AnnData, 
                                 n_neighbors: int = 6,
                                 spatial_key: str = "spatial",
                                 **kwargs) -> None:
        """
        Fallback spatial neighbor computation using sklearn.
        Used when squidpy is not available.
        """
        warn_fallback("spatial neighbors", "squidpy", "sklearn NearestNeighbors")
        
        if spatial_key not in adata.obsm:
            raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")
        
        coords = adata.obsm[spatial_key]
        
        # Use sklearn NearestNeighbors
        nn = NearestNeighbors(n_neighbors=n_neighbors + 1, metric='euclidean')
        nn.fit(coords)
        distances, indices = nn.kneighbors(coords)
        
        # Remove self (first neighbor)
        distances = distances[:, 1:]
        indices = indices[:, 1:]
        
        # Create connectivity matrix
        n_obs = adata.n_obs
        connectivities = sparse.lil_matrix((n_obs, n_obs))
        
        for i in range(n_obs):
            for j, neighbor_idx in enumerate(indices[i]):
                # Simple binary connectivity (can be enhanced)
                connectivities[i, neighbor_idx] = 1.0
                
        # Store in adata
        adata.obsp['spatial_connectivities'] = connectivities.tocsr()
        adata.obsp['spatial_distances'] = sparse.lil_matrix((n_obs, n_obs))
        
        for i in range(n_obs):
            for j, (neighbor_idx, dist) in enumerate(zip(indices[i], distances[i])):
                adata.obsp['spatial_distances'][i, neighbor_idx] = dist
                
        adata.obsp['spatial_distances'] = adata.obsp['spatial_distances'].tocsr()
    
    @staticmethod
    def spatial_autocorrelation_fallback(adata: sc.AnnData, 
                                       genes: Optional[List[str]] = None,
                                       **kwargs) -> pd.DataFrame:
        """
        Fallback Moran's I calculation using scipy.
        Used when squidpy is not available.
        """
        warn_fallback("spatial autocorrelation", "squidpy", "scipy-based Moran's I")
        
        if 'spatial_connectivities' not in adata.obsp:
            raise ValueError("Spatial neighbors must be computed first")
        
        W = adata.obsp['spatial_connectivities']
        
        if genes is None:
            # Use top 100 variable genes
            if 'highly_variable' in adata.var.columns:
                genes = adata.var[adata.var.highly_variable].index.tolist()[:100]
            else:
                genes = adata.var_names[:100].tolist()
        
        results = []
        X = adata.X if not sparse.issparse(adata.X) else adata.X.toarray()
        
        for gene in genes:
            if gene not in adata.var_names:
                continue
                
            gene_idx = adata.var_names.get_loc(gene)
            x = X[:, gene_idx]
            
            # Simple Moran's I calculation
            n = len(x)
            mean_x = np.mean(x)
            
            # Numerator: sum of spatial weights * deviations
            numerator = 0
            denominator = np.sum((x - mean_x) ** 2)
            W_sum = 0
            
            W_coo = W.tocoo()
            for i, j, w in zip(W_coo.row, W_coo.col, W_coo.data):
                numerator += w * (x[i] - mean_x) * (x[j] - mean_x)
                W_sum += w
            
            if denominator > 0 and W_sum > 0:
                I = (n / W_sum) * (numerator / denominator)
                
                # Expected value under null hypothesis
                E_I = -1 / (n - 1)
                
                results.append({
                    'gene': gene,
                    'I': I,
                    'E_I': E_I,
                    'p_value': np.nan,  # P-value calculation is complex, skip for fallback
                    'z_score': np.nan
                })
        
        return pd.DataFrame(results)


class CellCommunicationFallbacks:
    """Fallback implementations for cell communication analysis"""
    
    @staticmethod
    def simple_ligand_receptor_analysis(adata: sc.AnnData,
                                      cell_type_key: str = "cell_type",
                                      **kwargs) -> Dict[str, Any]:
        """
        Simple ligand-receptor analysis fallback.
        Used when liana/cellphonedb are not available.
        """
        warn_fallback("cell communication", "liana/cellphonedb", 
                     "simple correlation-based analysis")
        
        # Hardcoded small set of common ligand-receptor pairs
        # In production, this could be loaded from a curated database
        lr_pairs = [
            ("CCL2", "CCR2"), ("CCL3", "CCR1"), ("CCL4", "CCR5"),
            ("CXCL12", "CXCR4"), ("TNF", "TNFRSF1A"), ("IL1B", "IL1R1"),
            ("VEGFA", "KDR"), ("PDGFA", "PDGFRA"), ("FGF2", "FGFR1")
        ]
        
        results = []
        
        # Get expression data
        if sparse.issparse(adata.X):
            X = adata.X.toarray()
        else:
            X = adata.X
        
        # Get cell type assignments
        if cell_type_key not in adata.obs.columns:
            raise ValueError(f"Cell type key '{cell_type_key}' not found in adata.obs")
        
        cell_types = adata.obs[cell_type_key].unique()
        
        for ligand, receptor in lr_pairs:
            if ligand not in adata.var_names or receptor not in adata.var_names:
                continue
                
            ligand_idx = adata.var_names.get_loc(ligand)
            receptor_idx = adata.var_names.get_loc(receptor)
            
            ligand_expr = X[:, ligand_idx]
            receptor_expr = X[:, receptor_idx]
            
            # Calculate expression by cell type
            for sender_type in cell_types:
                for receiver_type in cell_types:
                    sender_mask = adata.obs[cell_type_key] == sender_type
                    receiver_mask = adata.obs[cell_type_key] == receiver_type
                    
                    if not sender_mask.any() or not receiver_mask.any():
                        continue
                    
                    sender_ligand = np.mean(ligand_expr[sender_mask])
                    receiver_receptor = np.mean(receptor_expr[receiver_mask])
                    
                    # Simple communication score
                    comm_score = sender_ligand * receiver_receptor
                    
                    if comm_score > 0:  # Only include non-zero scores
                        results.append({
                            'ligand': ligand,
                            'receptor': receptor,
                            'sender': sender_type,
                            'receiver': receiver_type,
                            'ligand_expr': sender_ligand,
                            'receptor_expr': receiver_receptor,
                            'communication_score': comm_score
                        })
        
        df = pd.DataFrame(results)
        
        return {
            'communication_scores': df,
            'method': 'simple_correlation_fallback',
            'n_pairs': len(lr_pairs),
            'n_significant': len(df[df['communication_score'] > df['communication_score'].median()]) if len(df) > 0 else 0
        }


class SpatialDomainFallbacks:
    """Fallback implementations for spatial domain identification"""
    
    @staticmethod
    def kmeans_spatial_domains(adata: sc.AnnData,
                             n_clusters: int = 7,
                             spatial_key: str = "spatial",
                             use_expression: bool = True,
                             **kwargs) -> None:
        """
        K-means based spatial domain identification fallback.
        Used when SpaGCN/STAGATE are not available.
        """
        warn_fallback("spatial domains", "SpaGCN/STAGATE", "K-means clustering")
        
        features = []
        
        # Add spatial coordinates
        if spatial_key in adata.obsm:
            coords = adata.obsm[spatial_key]
            features.append(coords)
        
        # Add expression features if requested
        if use_expression:
            # Use PCA to reduce dimensionality
            if 'X_pca' in adata.obsm:
                pca_features = adata.obsm['X_pca'][:, :50]  # Top 50 PCs
            else:
                # Compute PCA on the fly
                if sparse.issparse(adata.X):
                    X = adata.X.toarray()
                else:
                    X = adata.X.copy()
                
                pca = PCA(n_components=min(50, X.shape[1], X.shape[0]))
                pca_features = pca.fit_transform(X)
            
            features.append(pca_features)
        
        if not features:
            raise ValueError("No features available for clustering")
        
        # Combine features
        X_combined = np.concatenate(features, axis=1)
        
        # Apply K-means
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        clusters = kmeans.fit_predict(X_combined)
        
        # Store results
        adata.obs['spatial_domains'] = pd.Categorical(
            [f"Domain_{i}" for i in clusters]
        )
        
        # Calculate silhouette score for quality assessment
        if len(np.unique(clusters)) > 1:
            sil_score = silhouette_score(X_combined, clusters)
            print(f"Silhouette score: {sil_score:.3f}")
    
    @staticmethod 
    def expression_based_domains(adata: sc.AnnData,
                               resolution: float = 0.5,
                               **kwargs) -> None:
        """
        Expression-only spatial domain identification using Leiden clustering.
        """
        warn_fallback("spatial domains", "SpaGCN/STAGATE", "expression-based Leiden clustering")
        
        # Use standard scanpy workflow
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.leiden(adata, resolution=resolution, key_added='spatial_domains')


class DeconvolutionFallbacks:
    """Fallback implementations for spatial deconvolution"""
    
    @staticmethod
    def nnls_deconvolution(adata_spatial: sc.AnnData,
                          adata_reference: sc.AnnData,
                          cell_type_key: str = "cell_type",
                          **kwargs) -> Dict[str, Any]:
        """
        Non-negative least squares deconvolution fallback.
        Used when cell2location is not available.
        """
        warn_fallback("deconvolution", "cell2location", "NNLS-based deconvolution")
        
        from scipy.optimize import nnls
        
        # Get common genes
        common_genes = list(set(adata_spatial.var_names) & set(adata_reference.var_names))
        if len(common_genes) < 100:
            raise ValueError(f"Too few common genes ({len(common_genes)}). Need at least 100.")
        
        # Subset to common genes
        spatial_subset = adata_spatial[:, common_genes].copy()
        ref_subset = adata_reference[:, common_genes].copy()
        
        # Get reference cell type expression profiles
        cell_types = ref_subset.obs[cell_type_key].unique()
        
        # Calculate mean expression per cell type
        ref_profiles = {}
        for ct in cell_types:
            mask = ref_subset.obs[cell_type_key] == ct
            if sparse.issparse(ref_subset.X):
                mean_expr = np.array(ref_subset.X[mask].mean(axis=0)).flatten()
            else:
                mean_expr = ref_subset.X[mask].mean(axis=0)
            ref_profiles[ct] = mean_expr
        
        # Create reference matrix
        ref_matrix = np.column_stack([ref_profiles[ct] for ct in cell_types])
        
        # Deconvolve each spatial spot
        proportions = []
        
        if sparse.issparse(spatial_subset.X):
            spatial_data = spatial_subset.X.toarray()
        else:
            spatial_data = spatial_subset.X
        
        for i in range(spatial_subset.n_obs):
            spot_expression = spatial_data[i, :]
            
            # NNLS deconvolution
            try:
                props, residual = nnls(ref_matrix, spot_expression)
                # Normalize to sum to 1
                if props.sum() > 0:
                    props = props / props.sum()
                else:
                    props = np.zeros_like(props)
            except Exception:
                # If NNLS fails, use zeros
                props = np.zeros(len(cell_types))
            
            proportions.append(props)
        
        # Store results
        prop_matrix = np.array(proportions)
        
        for i, ct in enumerate(cell_types):
            adata_spatial.obs[f"prop_{ct}"] = prop_matrix[:, i]
        
        return {
            'proportions': pd.DataFrame(prop_matrix, columns=cell_types, 
                                      index=adata_spatial.obs_names),
            'method': 'nnls_fallback',
            'n_cell_types': len(cell_types),
            'n_common_genes': len(common_genes)
        }


class VisualizationFallbacks:
    """Fallback implementations for advanced visualization"""
    
    @staticmethod
    def simple_spatial_plot(adata: sc.AnnData,
                          color: Optional[Union[str, List[str]]] = None,
                          spatial_key: str = "spatial",
                          **kwargs) -> plt.Figure:
        """
        Simple spatial plot fallback using matplotlib.
        Used when advanced plotting libraries are not available.
        """
        warn_fallback("spatial visualization", "advanced plotting", "matplotlib scatter plot")
        
        if spatial_key not in adata.obsm:
            raise ValueError(f"Spatial coordinates not found in adata.obsm['{spatial_key}']")
        
        coords = adata.obsm[spatial_key]
        
        fig, ax = plt.subplots(figsize=(8, 8))
        
        if color is None:
            # Just plot coordinates
            ax.scatter(coords[:, 0], coords[:, 1], s=20, alpha=0.7)
        elif isinstance(color, str):
            if color in adata.obs.columns:
                # Categorical or continuous color
                color_data = adata.obs[color]
                if color_data.dtype == 'category' or color_data.dtype == 'object':
                    # Categorical
                    for i, category in enumerate(color_data.cat.categories if hasattr(color_data, 'cat') else color_data.unique()):
                        mask = color_data == category
                        ax.scatter(coords[mask, 0], coords[mask, 1], 
                                 s=20, alpha=0.7, label=category)
                    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                else:
                    # Continuous
                    scatter = ax.scatter(coords[:, 0], coords[:, 1], 
                                       c=color_data, s=20, alpha=0.7, cmap='viridis')
                    plt.colorbar(scatter, ax=ax)
            elif color in adata.var_names:
                # Gene expression
                gene_idx = adata.var_names.get_loc(color)
                if sparse.issparse(adata.X):
                    expr = adata.X[:, gene_idx].toarray().flatten()
                else:
                    expr = adata.X[:, gene_idx]
                
                scatter = ax.scatter(coords[:, 0], coords[:, 1], 
                                   c=expr, s=20, alpha=0.7, cmap='viridis')
                plt.colorbar(scatter, ax=ax)
                ax.set_title(f"Expression of {color}")
        
        ax.set_xlabel("Spatial X")
        ax.set_ylabel("Spatial Y")
        ax.set_aspect('equal')
        
        plt.tight_layout()
        return fig


# Registry of fallback functions
FALLBACK_REGISTRY = {
    'spatial_neighbors': SpatialFallbacks.spatial_neighbors_fallback,
    'spatial_autocorrelation': SpatialFallbacks.spatial_autocorrelation_fallback,
    'cell_communication': CellCommunicationFallbacks.simple_ligand_receptor_analysis,
    'spatial_domains_kmeans': SpatialDomainFallbacks.kmeans_spatial_domains,
    'spatial_domains_leiden': SpatialDomainFallbacks.expression_based_domains,
    'deconvolution': DeconvolutionFallbacks.nnls_deconvolution,
    'spatial_plot': VisualizationFallbacks.simple_spatial_plot,
}


def get_fallback(method_name: str) -> Optional[Callable]:
    """Get fallback function for a method"""
    return FALLBACK_REGISTRY.get(method_name)


def list_available_fallbacks() -> List[str]:
    """List all available fallback methods"""
    return list(FALLBACK_REGISTRY.keys())


class SafeMode:
    """
    Safe mode operations using only core dependencies.
    Provides basic functionality when optional dependencies are missing.
    """
    
    def __init__(self, adata: sc.AnnData):
        self.adata = adata
        
    def basic_preprocessing(self) -> None:
        """Basic preprocessing using only scanpy"""
        print("Running safe mode preprocessing (scanpy only)...")
        
        # Basic QC
        sc.pp.calculate_qc_metrics(self.adata, inplace=True)
        
        # Filter cells and genes
        sc.pp.filter_cells(self.adata, min_genes=200)
        sc.pp.filter_genes(self.adata, min_cells=3)
        
        # Normalize
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)
        
        # Find highly variable genes
        sc.pp.highly_variable_genes(self.adata, n_top_genes=2000)
        
        # PCA
        sc.tl.pca(self.adata, n_comps=50)
        
        # Neighbors
        sc.pp.neighbors(self.adata)
        
        # UMAP
        sc.tl.umap(self.adata)
        
        # Leiden clustering
        sc.tl.leiden(self.adata)
        
        print("Safe mode preprocessing completed.")
    
    def basic_visualization(self) -> None:
        """Basic visualization using scanpy"""
        print("Creating safe mode visualizations...")
        
        # UMAP plots
        sc.pl.umap(self.adata, color='leiden', show=False)
        plt.title("Leiden Clustering (Safe Mode)")
        plt.show()
        
        # If spatial coordinates are available
        if 'spatial' in self.adata.obsm:
            fig = VisualizationFallbacks.simple_spatial_plot(
                self.adata, color='leiden'
            )
            plt.title("Spatial Plot (Safe Mode)")
            plt.show()
    
    def get_capabilities(self) -> Dict[str, bool]:
        """Get what's available in safe mode"""
        return {
            'preprocessing': True,
            'clustering': True,
            'dimensionality_reduction': True,
            'basic_visualization': True,
            'differential_expression': True,
            'spatial_neighbors': False,
            'advanced_spatial_analysis': False,
            'cell_communication': False,
            'advanced_deconvolution': False,
            'trajectory_analysis': True  # Basic with scanpy
        }


def enable_safe_mode(adata: sc.AnnData) -> SafeMode:
    """Enable safe mode for basic analysis"""
    print("🛡️  Enabling ChatSpatial Safe Mode")
    print("Using only core dependencies for basic analysis.")
    print("For full functionality, install optional dependencies.")
    
    return SafeMode(adata)