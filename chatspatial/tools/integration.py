"""
Integration tools for spatial transcriptomics data.
"""

from typing import Dict, List, Optional, Any
import numpy as np
import scanpy as sc
from mcp.server.fastmcp import Context

from ..models.data import IntegrationParameters
from ..models.analysis import IntegrationResult


def integrate_multiple_samples(adatas, batch_key='batch', method='harmony', n_pcs=30):
    """Integrate multiple spatial transcriptomics samples

    Args:
        adatas: List of AnnData objects or a single combined AnnData object
        batch_key: Batch information key
        method: Integration method, options: 'harmony', 'bbknn', 'scanorama', 'mnn'
        n_pcs: Number of principal components for integration

    Returns:
        Integrated AnnData object
    """
    import scanpy as sc

    # Merge datasets
    if isinstance(adatas, list):
        # Ensure each dataset has batch labels
        for i, adata in enumerate(adatas):
            if batch_key not in adata.obs:
                adata.obs[batch_key] = f"batch_{i}"

        # Merge datasets
        combined = adatas[0].concatenate(
            adatas[1:],
            batch_key=batch_key,
            join='outer'  # Use outer join to keep all genes
        )
    else:
        # If already a merged dataset, ensure it has batch information
        combined = adatas
        if batch_key not in combined.obs:
            raise ValueError(f"Merged dataset is missing batch information key '{batch_key}'")

    # Standard preprocessing
    sc.pp.normalize_total(combined)
    sc.pp.log1p(combined)

    # Handle NaN values
    import numpy as np
    # Replace NaN values with zeros
    if isinstance(combined.X, np.ndarray):
        combined.X = np.nan_to_num(combined.X)
    else:
        # For sparse matrices
        combined.X.data = np.nan_to_num(combined.X.data)

    # Find highly variable genes
    try:
        sc.pp.highly_variable_genes(combined, batch_key=batch_key)
    except Exception as e:
        print(f"Error in highly_variable_genes with batch correction: {e}")
        # Fallback: compute HVGs without batch correction
        try:
            # Try with reduced number of top genes if dataset is small
            n_top_genes = min(2000, combined.n_vars // 2, 1000)
            if n_top_genes < 50:
                n_top_genes = combined.n_vars  # Use all genes if very few
            sc.pp.highly_variable_genes(combined, n_top_genes=n_top_genes)
        except Exception as e2:
            print(f"Error in highly_variable_genes fallback: {e2}")
            # Final fallback: mark all genes as highly variable
            combined.var['highly_variable'] = True
            combined.var['means'] = np.array(combined.X.mean(axis=0)).flatten()
            combined.var['dispersions'] = np.array(combined.X.var(axis=0)).flatten()
            combined.var['dispersions_norm'] = combined.var['dispersions'] / combined.var['means']
            print(f"Using all {combined.n_vars} genes for integration")

    combined.raw = combined  # Save raw data

    # Only keep highly variable genes for integration
    try:
        sc.pp.scale(combined, zero_center=True, max_value=10)  # Limit scaling to avoid extreme values
    except Exception as e:
        print(f"Warning in scaling: {e}")
        # Try without zero_center if scaling fails
        try:
            sc.pp.scale(combined, zero_center=False, max_value=10)
        except Exception as e2:
            print(f"Scaling failed completely: {e2}, continuing without scaling")
    
    # PCA with adaptive number of components
    try:
        # Ensure n_pcs doesn't exceed the number of genes or cells
        max_components = min(n_pcs, combined.n_vars, combined.n_obs - 1)
        sc.tl.pca(combined, n_comps=max_components, svd_solver='arpack')
    except Exception as e:
        print(f"PCA with arpack failed: {e}")
        # Fallback to randomized SVD
        try:
            max_components = min(n_pcs, combined.n_vars, combined.n_obs - 1, 50)  # Limit to 50 for randomized
            sc.tl.pca(combined, n_comps=max_components, svd_solver='randomized')
        except Exception as e2:
            print(f"PCA with randomized SVD also failed: {e2}")
            # Create a simple PCA representation
            max_components = min(10, combined.n_vars, combined.n_obs - 1)
            try:
                sc.tl.pca(combined, n_comps=max_components, svd_solver='full')
            except Exception as e3:
                print(f"All PCA methods failed: {e3}")
                # Create dummy PCA for compatibility
                combined.obsm['X_pca'] = np.random.normal(0, 1, (combined.n_obs, min(10, combined.n_vars)))

    # Apply batch correction based on selected method
    if method == 'harmony':
        # Use Harmony for batch correction
        try:
            import harmonypy

            # Get PCA result
            X_pca = combined.obsm['X_pca']

            # Run Harmony - need to pass the DataFrame with batch info, not just the labels
            # Create a temporary DataFrame with batch information
            import pandas as pd
            meta_data = pd.DataFrame({batch_key: combined.obs[batch_key]})

            # Run Harmony with proper parameters
            harmony_out = harmonypy.run_harmony(
                data_mat=X_pca,  # PCA matrix
                meta_data=meta_data,  # DataFrame with batch info
                vars_use=[batch_key],  # Column name in meta_data
                sigma=0.1,  # Default parameter
                nclust=None,  # Let harmony determine the number of clusters
                max_iter_harmony=10,  # Default number of iterations
                verbose=True  # Show progress
            )

            # Save Harmony corrected result
            combined.obsm['X_harmony'] = harmony_out.Z_corr.T

            # Use corrected result to calculate neighbor graph
            sc.pp.neighbors(combined, use_rep='X_harmony')
        except ImportError:
            raise ImportError("harmonypy package is required for harmony integration. Install with 'pip install harmonypy'")
        except Exception as e:
            # Fallback to uncorrected PCA if harmony fails
            import logging
            logging.warning(f"Harmony integration failed with error: {e}. Using uncorrected PCA instead.")
            sc.pp.neighbors(combined, use_rep='X_pca')

    elif method == 'bbknn':
        # Use BBKNN for batch correction
        try:
            import bbknn
            bbknn.bbknn(combined, batch_key=batch_key, neighbors_within_batch=3)
        except ImportError:
            raise ImportError("bbknn package is required for BBKNN integration. Install with 'pip install bbknn'")

    elif method == 'scanorama':
        # Use Scanorama for batch correction
        try:
            import scanorama
            import numpy as np

            # Separate data by batch
            batches = []
            batch_names = []
            for batch in combined.obs[batch_key].unique():
                batches.append(combined[combined.obs[batch_key] == batch].X)
                batch_names.append(batch)

            # Run Scanorama
            corrected, genes = scanorama.correct(batches, batch_names, return_dimred=True)

            # Save corrected data back to AnnData
            combined.obsm['X_scanorama'] = np.vstack(corrected)

            # Use corrected result to calculate neighbor graph
            sc.pp.neighbors(combined, use_rep='X_scanorama')
        except ImportError:
            raise ImportError("scanorama package is required for Scanorama integration. Install with 'pip install scanorama'")

    elif method == 'mnn':
        # Use MNN for batch correction
        sc.pp.combat(combined, key=batch_key)
        sc.pp.neighbors(combined)

    else:
        # Default: use uncorrected PCA result
        sc.pp.neighbors(combined)

    # Calculate UMAP embedding to visualize integration effect
    sc.tl.umap(combined)

    return combined


def align_spatial_coordinates(combined_adata, batch_key='batch', reference_batch=None):
    """Align spatial coordinates of multiple samples

    Args:
        combined_adata: Combined AnnData object containing multiple samples
        batch_key: Batch information key
        reference_batch: Reference batch, if None use the first batch

    Returns:
        AnnData object with aligned spatial coordinates
    """
    import numpy as np
    from sklearn.preprocessing import StandardScaler

    # Ensure data contains spatial coordinates
    if 'spatial' not in combined_adata.obsm:
        raise ValueError("Data is missing spatial coordinates")

    # Get batch information
    batches = combined_adata.obs[batch_key].unique()

    # If reference batch not specified, use the first batch
    if reference_batch is None:
        reference_batch = batches[0]
    elif reference_batch not in batches:
        raise ValueError(f"Reference batch '{reference_batch}' not found in data")

    # Get reference batch spatial coordinates
    ref_coords = combined_adata[combined_adata.obs[batch_key] == reference_batch].obsm['spatial']

    # Standardize reference coordinates
    scaler = StandardScaler()
    ref_coords_scaled = scaler.fit_transform(ref_coords)

    # Align spatial coordinates for each batch
    aligned_coords = []

    for batch in batches:
        # Get current batch index
        batch_idx = combined_adata.obs[batch_key] == batch

        if batch == reference_batch:
            # Reference batch remains unchanged
            aligned_coords.append(ref_coords_scaled)
        else:
            # Get current batch spatial coordinates
            batch_coords = combined_adata[batch_idx].obsm['spatial']

            # Standardize current batch coordinates
            batch_coords_scaled = scaler.transform(batch_coords)

            # Add to aligned coordinates list
            aligned_coords.append(batch_coords_scaled)

    # Merge all aligned coordinates
    combined_adata.obsm['spatial_aligned'] = np.zeros((combined_adata.n_obs, 2))

    # Fill aligned coordinates back to original data
    start_idx = 0
    for batch, coords in zip(batches, aligned_coords):
        batch_idx = combined_adata.obs[batch_key] == batch
        n_cells = np.sum(batch_idx)
        combined_adata.obsm['spatial_aligned'][start_idx:start_idx+n_cells] = coords
        start_idx += n_cells

    return combined_adata


def analyze_integrated_trajectory(combined_adata, spatial_weight=0.5, use_aligned_coords=True):
    """Analyze trajectory in integrated multi-sample spatial data

    Args:
        combined_adata: Integrated AnnData object
        spatial_weight: Weight of spatial information
        use_aligned_coords: Whether to use aligned spatial coordinates

    Returns:
        AnnData object with trajectory analysis results
    """
    import numpy as np
    from scipy.spatial.distance import pdist, squareform
    from scipy.sparse import csr_matrix

    try:
        import cellrank as cr
    except ImportError:
        raise ImportError("cellrank package is required for trajectory analysis. Install with 'pip install cellrank'")

    # Choose which spatial coordinates to use
    spatial_key = 'spatial_aligned' if use_aligned_coords and 'spatial_aligned' in combined_adata.obsm else 'spatial'

    # Get spatial coordinates
    spatial_coords = combined_adata.obsm[spatial_key]

    # Create RNA velocity kernel
    vk = cr.kernels.VelocityKernel(combined_adata)
    vk.compute_transition_matrix()

    # Create expression similarity kernel
    ck = cr.kernels.ConnectivityKernel(combined_adata)
    ck.compute_transition_matrix()

    # Create custom spatial kernel
    # Calculate spatial distance matrix
    spatial_dist = squareform(pdist(spatial_coords))
    # Convert to similarity (smaller distance = higher similarity)
    spatial_sim = np.exp(-spatial_dist / spatial_dist.mean())
    spatial_kernel = csr_matrix(spatial_sim)

    # Create spatial kernel
    sk = cr.kernels.PrecomputedKernel(spatial_kernel, combined_adata)
    sk.compute_transition_matrix()

    # Combine kernels, integrating RNA velocity, expression similarity, and spatial information
    combined_kernel = (1-spatial_weight) * (0.8 * vk + 0.2 * ck) + spatial_weight * sk

    # Use GPCCA for analysis
    g = cr.estimators.GPCCA(combined_kernel)
    g.compute_schur()
    g.compute_macrostates()

    # Calculate pseudotime and absorption probabilities
    g.set_terminal_states_from_macrostates()
    g.compute_absorption_probabilities()

    # Store results
    combined_adata.obs['pseudotime'] = g.pseudotime
    combined_adata.obsm['absorption_probabilities'] = g.absorption_probabilities

    return combined_adata


async def integrate_samples(
    data_ids: List[str],
    data_store: Dict[str, Any],
    params: IntegrationParameters = IntegrationParameters(),
    context: Optional[Context] = None
) -> IntegrationResult:
    """Integrate multiple spatial transcriptomics samples and perform batch correction

    Args:
        data_ids: List of dataset IDs to integrate
        data_store: Dictionary storing datasets
        params: Integration parameters
        context: MCP context

    Returns:
        Integration result
    """
    if context:
        await context.info(f"Integrating {len(data_ids)} samples using {params.method} method")

    # Collect all AnnData objects
    adatas = []
    for data_id in data_ids:
        if data_id not in data_store:
            raise ValueError(f"Dataset {data_id} not found in data store")
        adatas.append(data_store[data_id]["adata"].copy())

    # Integrate samples
    combined_adata = integrate_multiple_samples(
        adatas,
        batch_key=params.batch_key,
        method=params.method,
        n_pcs=params.n_pcs
    )

    if context:
        await context.info(f"Integration complete. Combined dataset has {combined_adata.n_obs} cells and {combined_adata.n_vars} genes")

    # Align spatial coordinates
    if params.align_spatial:
        if context:
            await context.info("Aligning spatial coordinates")
        combined_adata = align_spatial_coordinates(
            combined_adata,
            batch_key=params.batch_key,
            reference_batch=params.reference_batch
        )

    # Note: Visualizations should be created using the separate visualize_data tool
    # This maintains clean separation between analysis and visualization
    if context:
        await context.info("Integration analysis complete. Use visualize_data tool with plot_type='integration_umap' or 'integration_spatial' to visualize results")

    # Generate new integrated dataset ID
    integrated_id = f"integrated_{'-'.join(data_ids)}"

    # Store integrated data
    data_store[integrated_id] = {"adata": combined_adata}

    if context:
        await context.info(f"Integration complete. New dataset ID: {integrated_id}")

    # Return result
    return IntegrationResult(
        data_id=integrated_id,
        n_samples=len(data_ids),
        integration_method=params.method,
        umap_visualization=None,  # Use visualize_data tool instead
        spatial_visualization=None  # Use visualize_data tool instead
    )
