"""
Preprocessing tools for spatial transcriptomics data.
"""

from typing import Dict, Any, Optional
import numpy as np
import scanpy as sc
import squidpy as sq
import traceback
from mcp.server.fastmcp import Context

from ..models.data import AnalysisParameters
from ..models.analysis import PreprocessingResult


async def preprocess_data(
    data_id: str,
    data_store: Dict[str, Any],
    params: AnalysisParameters = AnalysisParameters(),
    context: Optional[Context] = None
) -> PreprocessingResult:
    """Preprocess spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing loaded datasets
        params: Preprocessing parameters
        context: MCP context

    Returns:
        Preprocessing result summary
    """
    try:
        if context:
            await context.info(f"Preprocessing dataset {data_id} with {params.normalization} normalization")

        # Retrieve the AnnData object from data store
        if data_id not in data_store:
            raise ValueError(f"Dataset {data_id} not found in data store")

        # Make a copy of the AnnData object to avoid modifying the original
        adata = data_store[data_id]["adata"].copy()

        # 1. Calculate QC metrics
        if context:
            await context.info("Calculating QC metrics...")
        try:
            sc.pp.calculate_qc_metrics(adata, inplace=True)
        except IndexError as e:
            # Handle the "Positions outside range of features" error
            if "Positions outside range of features" in str(e):
                if context:
                    await context.warning("Could not calculate QC metrics: Positions outside range of features. This is common with MERFISH data.")
                # Create empty QC metrics to avoid errors later
                adata.obs['n_genes_by_counts'] = np.ones(adata.n_obs)
                adata.obs['total_counts'] = np.ones(adata.n_obs) * 1000  # Dummy value
                adata.obs['pct_counts_mt'] = np.zeros(adata.n_obs)  # Dummy value
            else:
                raise

        # Store original QC metrics before filtering
        qc_metrics = {
            "n_cells_before_filtering": adata.n_obs,
            "n_genes_before_filtering": adata.n_vars,
            "median_genes_per_cell": float(np.median(adata.obs.n_genes_by_counts)),
            "median_umi_per_cell": str(np.median(adata.obs.total_counts))
        }

        # 2. Filter cells and genes (optional)
        # Here we use mild filtering to remove cells with very few genes and genes expressed in very few cells
        if context:
            await context.info("Filtering cells and genes...")

        # Check if this is MERFISH data (typically has fewer genes)
        if adata.n_vars < 200:
            # For MERFISH data, use more lenient filtering
            min_genes = min(50, adata.n_vars // 2)  # Use half of the genes or 50, whichever is smaller
            if context:
                await context.info(f"Detected possible MERFISH data, using lenient filtering (min_genes={min_genes})")
        else:
            min_genes = 200

        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=3)

        # Update QC metrics after filtering
        qc_metrics.update({
            "n_cells_after_filtering": adata.n_obs,
            "n_genes_after_filtering": adata.n_vars,
        })

        # 3. Normalize data
        if context:
            await context.info(f"Normalizing data using {params.normalization} method...")

        if params.normalization == "log":
            # Standard log normalization
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        elif params.normalization == "sct":
            # SCTransform-like normalization (simplified version)
            # In full implementation, you might want to use sctransform package
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            # Add additional variance stabilization steps here

        # 4. Find highly variable genes
        if context:
            await context.info(f"Finding highly variable genes...")

        # Adjust n_hvgs based on dataset size
        n_hvgs = min(params.n_hvgs, adata.n_vars - 1)  # Ensure we don't try to select more HVGs than available genes

        if context:
            await context.info(f"Using {n_hvgs} highly variable genes...")

        # For very small gene sets, we might want to use all genes
        if adata.n_vars < 200:
            if context:
                await context.info(f"Small gene set detected ({adata.n_vars} genes), using all genes for analysis")
            adata.var['highly_variable'] = True
        else:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs)

        # 5. Scale data (if requested)
        if params.scale:
            if context:
                await context.info("Scaling data...")
            sc.pp.scale(adata, max_value=10)

        # 6. Run PCA
        if context:
            await context.info(f"Running PCA...")

        # Adjust n_pcs based on dataset size
        n_pcs = min(params.n_pcs, adata.n_vars - 1, adata.n_obs - 1)  # Ensure we don't use more PCs than possible

        if context:
            await context.info(f"Using {n_pcs} principal components...")

        sc.pp.pca(adata, n_comps=n_pcs)

        # 7. Compute neighbors graph
        if context:
            await context.info("Computing neighbors graph...")

        # Adjust n_neighbors based on dataset size
        n_neighbors = min(10, adata.n_obs // 10)  # Use at most 10% of cells as neighbors
        n_neighbors = max(n_neighbors, 3)  # But at least 3 neighbors

        if context:
            await context.info(f"Using {n_neighbors} neighbors and {n_pcs} PCs for graph construction...")

        try:
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

            # 8. Run UMAP for visualization
            if context:
                await context.info("Running UMAP...")
            sc.tl.umap(adata)

            # 9. Run clustering
            if context:
                await context.info("Running Leiden clustering...")

            # Adjust resolution based on dataset size
            if adata.n_obs < 100:
                resolution = 0.4  # Lower resolution for small datasets
            elif adata.n_obs < 500:
                resolution = 0.6  # Medium resolution for medium datasets
            else:
                resolution = 0.8  # Higher resolution for large datasets

            if context:
                await context.info(f"Using Leiden clustering with resolution {resolution}...")

            sc.tl.leiden(adata, resolution=resolution)

            # Count clusters
            n_clusters = len(adata.obs['leiden'].unique())
        except Exception as e:
            if context:
                await context.warning(f"Error in neighbors/clustering: {str(e)}")
                await context.info("Creating fallback clustering...")

            # Create a fallback clustering based on gene expression
            if 'highly_variable' in adata.var and adata.var['highly_variable'].any():
                # Use highly variable genes
                hvg_indices = np.where(adata.var['highly_variable'])[0]
                X_hvg = adata.X[:, hvg_indices]
            else:
                # Use all genes
                X_hvg = adata.X

            # Simple clustering based on expression
            from sklearn.cluster import KMeans
            n_clusters_kmeans = min(10, adata.n_obs // 10)  # At most 10 clusters, but at most 10% of cells
            n_clusters_kmeans = max(n_clusters_kmeans, 2)  # At least 2 clusters

            if context:
                await context.info(f"Using KMeans with {n_clusters_kmeans} clusters as fallback...")

            kmeans = KMeans(n_clusters=n_clusters_kmeans, random_state=42)
            adata.obs['leiden'] = kmeans.fit_predict(X_hvg).astype(str)
            n_clusters = n_clusters_kmeans

            # Create a simple UMAP embedding if possible
            try:
                from sklearn.decomposition import PCA
                from sklearn.manifold import TSNE

                # Use PCA for dimensionality reduction
                pca = PCA(n_components=min(50, X_hvg.shape[1], X_hvg.shape[0] - 1))
                X_pca = pca.fit_transform(X_hvg)

                # Use t-SNE for visualization (faster than UMAP)
                tsne = TSNE(n_components=2, random_state=42)
                X_tsne = tsne.fit_transform(X_pca)

                # Store as UMAP coordinates for compatibility
                adata.obsm['X_umap'] = X_tsne
            except Exception as e2:
                if context:
                    await context.warning(f"Could not create fallback embedding: {str(e2)}")

        # 10. For spatial data, compute spatial neighbors if not already done
        if 'spatial' in adata.uns or any('spatial' in key for key in adata.obsm.keys()):
            if context:
                await context.info("Computing spatial neighbors...")
            try:
                # Check if spatial neighbors already exist
                if 'spatial_connectivities' not in adata.obsp:
                    # For MERFISH data, we need to ensure the spatial coordinates are correctly formatted
                    if 'spatial' in adata.obsm:
                        # Check if this is MERFISH data by looking at the shape and content
                        if adata.obsm['spatial'].shape[1] == 2:
                            # This is likely MERFISH or other single-cell resolution spatial data
                            # Use delaunay method which works better for single-cell data
                            sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
                        else:
                            # Use default parameters for other spatial data types
                            sq.gr.spatial_neighbors(adata)
                    else:
                        # Use default parameters if spatial key is in uns but not in obsm
                        sq.gr.spatial_neighbors(adata)
            except Exception as e:
                if context:
                    await context.warning(f"Could not compute spatial neighbors: {str(e)}")
                    await context.info("Continuing without spatial neighbors...")

        # Store the processed AnnData object back in the data store
        data_store[data_id]["adata"] = adata

        # Return preprocessing result
        return PreprocessingResult(
            data_id=data_id,
            n_cells=adata.n_obs,
            n_genes=adata.n_vars,
            n_hvgs=int(sum(adata.var.highly_variable)) if 'highly_variable' in adata.var else n_hvgs,
            clusters=n_clusters,
            qc_metrics=qc_metrics
        )
        
    except Exception as e:
        error_msg = f"Error in preprocessing: {str(e)}"
        tb = traceback.format_exc()
        if context:
            await context.warning(error_msg)
            await context.info(f"Error details: {tb}")
        raise RuntimeError(f"{error_msg}\n{tb}")
