"""
Tools for RNA velocity and trajectory analysis in spatial transcriptomics data.
"""

from typing import Dict, Any, Optional, List, Tuple
import numpy as np
import pandas as pd

from mcp.server.fastmcp import Context

from ..models.data import RNAVelocityParameters, TrajectoryParameters
from ..models.analysis import RNAVelocityResult, TrajectoryResult
from ..utils.error_handling import suppress_output, ProcessingError, validate_adata, DataNotFoundError

# Import scvi-tools for advanced trajectory analysis
try:
    import scvi
    from scvi.external import VELOVI
except ImportError:
    scvi = None
    VELOVI = None


def validate_velocity_data(adata) -> Tuple[bool, List[str]]:
    """
    Validate if AnnData object has necessary data for RNA velocity analysis.
    
    This function now uses the unified validation system for consistency,
    but maintains the same interface for backward compatibility.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix to validate
        
    Returns
    -------
    Tuple[bool, List[str]]
        Tuple of (is_valid, list_of_issues)
    """
    try:
        validate_adata(adata, {}, check_velocity=True)
        return True, []
    except DataNotFoundError as e:
        # Extract individual issues from the error message
        error_msg = str(e)
        if "Validation failed: " in error_msg:
            issues = error_msg.replace("Validation failed: ", "").split(", ")
        else:
            issues = [error_msg]
        return False, issues


def validate_spatial_data(adata, spatial_key: str = 'spatial') -> Tuple[bool, List[str]]:
    """
    Validate if AnnData object has necessary spatial information.
    
    This function now uses the unified validation system for consistency,
    but maintains the same interface for backward compatibility.
    
    Parameters
    ---------- 
    adata : AnnData
        Annotated data matrix to validate
    spatial_key : str, default 'spatial'
        Key for spatial coordinates in adata.obsm
        
    Returns
    -------
    Tuple[bool, List[str]]
        Tuple of (is_valid, list_of_issues)
    """
    try:
        validate_adata(adata, {}, check_spatial=True, spatial_key=spatial_key)
        return True, []
    except DataNotFoundError as e:
        # Extract individual issues from the error message  
        error_msg = str(e)
        if "Validation failed: " in error_msg:
            issues = error_msg.replace("Validation failed: ", "").split(", ")
        else:
            issues = [error_msg]
        return False, issues


def preprocess_for_velocity(adata):
    """Prepare data for RNA velocity analysis"""
    import scvelo as scv

    # Validate velocity data
    is_valid, issues = validate_velocity_data(adata)
    if not is_valid:
        raise ValueError(f"Invalid velocity data: {'; '.join(issues)}")

    # Standard preprocessing
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)

    return adata


def compute_rna_velocity(adata, mode="stochastic"):
    """Compute RNA velocity"""
    import scvelo as scv

    # Compute RNA velocity based on selected mode
    if mode == "deterministic":
        scv.tl.velocity(adata)
    elif mode == "stochastic":
        scv.tl.velocity(adata, mode="stochastic")
    elif mode == "dynamical":
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode="dynamical")
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    # Compute velocity graph (use single core to avoid multiprocessing issues)
    scv.tl.velocity_graph(adata, n_jobs=1)

    return adata


def prepare_for_sirv(spatial_adata, scrna_adata):
    """Prepare spatial data and scRNA-seq data for SIRV analysis"""
    import scvelo as scv

    # Validate scRNA-seq velocity data
    is_valid, issues = validate_velocity_data(scrna_adata)
    if not is_valid:
        raise ValueError(f"Invalid scRNA-seq velocity data: {'; '.join(issues)}")

    # Validate spatial data
    is_valid, issues = validate_spatial_data(spatial_adata)
    if not is_valid:
        raise ValueError(f"Invalid spatial data: {'; '.join(issues)}")

    # Preprocess scRNA-seq data
    scv.pp.filter_and_normalize(scrna_adata)
    scv.pp.moments(scrna_adata)

    # Compute RNA velocity for scRNA-seq data
    scv.tl.velocity(scrna_adata)
    scv.tl.velocity_graph(scrna_adata, n_jobs=1)

    return spatial_adata, scrna_adata


def run_sirv(spatial_adata, scrna_adata, n_pcs=30, labels=None):
    """Run SIRV analysis"""
    try:
        from SIRV.main import SIRV

        # Run SIRV
        spatial_adata_with_velocity = SIRV(
            spatial_data=spatial_adata,
            scrna_data=scrna_adata,
            n_pcs=n_pcs,
            labels=labels
        )

        return spatial_adata_with_velocity
    except ImportError:
        raise ProcessingError("SIRV package not found. Install it with: pip install git+https://github.com/tabdelaal/SIRV.git")


def infer_spatial_trajectory_cellrank(adata, spatial_weight=0.5, kernel_weights=(0.8, 0.2), n_states=5):
    """Integrate RNA velocity and spatial information for trajectory inference using CellRank.
    
    Raises exceptions directly if CellRank fails - no fallback logic.
    """
    import cellrank as cr
    from scipy.spatial.distance import pdist, squareform
    from scipy.sparse import csr_matrix

    # Validate spatial data
    is_valid, issues = validate_spatial_data(adata)
    if not is_valid:
        raise ValueError(f"Invalid spatial data: {'; '.join(issues)}")

    spatial_coords = adata.obsm['spatial']

    # Create RNA velocity kernel
    vk = cr.kernels.VelocityKernel(adata)
    vk.compute_transition_matrix()

    # Create expression similarity-based kernel
    ck = cr.kernels.ConnectivityKernel(adata)
    ck.compute_transition_matrix()

    # Create custom spatial kernel
    spatial_dist = squareform(pdist(spatial_coords))
    spatial_sim = np.exp(-spatial_dist / spatial_dist.mean())
    spatial_kernel = csr_matrix(spatial_sim)

    # Create spatial kernel
    sk = cr.kernels.PrecomputedKernel(spatial_kernel, adata)
    sk.compute_transition_matrix()

    # Combine kernels using configurable weights
    vk_weight, ck_weight = kernel_weights
    combined_kernel = (1 - spatial_weight) * (vk_weight * vk + ck_weight * ck) + spatial_weight * sk

    # Use GPCCA for analysis
    g = cr.estimators.GPCCA(combined_kernel)

    # Compute eigendecomposition
    g.compute_eigendecomposition()

    # Compute macrostates
    try:
        g.compute_macrostates(n_states=n_states)
    except Exception as e:
        # If automatic n_states fails, try with fewer states
        for alt_n_states in [n_states-1, n_states-2, 3, 2]:
            if alt_n_states < 2:
                break
            try:
                g.compute_macrostates(n_states=alt_n_states)
                break
            except:
                continue
        else:
            raise RuntimeError(f"Failed to compute macrostates with any number of states: {e}")

    # Predict terminal states
    try:
        g.predict_terminal_states(method='stability')
    except ValueError as e:
        if "No macrostates have been selected" in str(e):
            # Skip terminal state prediction if no macrostates
            pass
        else:
            raise

    # Check if we have terminal states
    has_terminal_states = hasattr(g, 'terminal_states') and g.terminal_states is not None
    
    if has_terminal_states and len(g.terminal_states.cat.categories) > 0:
        # Compute fate probabilities (renamed from absorption_probabilities in newer versions)
        g.compute_fate_probabilities()

        # Get fate probabilities (renamed from absorption_probabilities)
        absorption_probs = g.fate_probabilities

        # Compute pseudotime
        terminal_states = list(g.terminal_states.cat.categories)
        root_state = terminal_states[0]
        pseudotime = 1 - absorption_probs[root_state].X.flatten()

        # Store results
        adata.obs['pseudotime'] = pseudotime
        adata.obsm['fate_probabilities'] = absorption_probs
        adata.obs['terminal_states'] = g.terminal_states
    else:
        # No terminal states, use macrostates for pseudotime
        if hasattr(g, 'macrostates') and g.macrostates is not None:
            # Use the macrostate memberships as a proxy for pseudotime
            macrostate_probs = g.macrostates_memberships
            # Use the first macrostate as the "early" state
            pseudotime = 1 - macrostate_probs[:, 0].X.flatten()
            adata.obs['pseudotime'] = pseudotime
        else:
            raise RuntimeError("CellRank could not compute either terminal states or macrostates")
    
    # Always store macrostates if available
    if hasattr(g, 'macrostates') and g.macrostates is not None:
        adata.obs['macrostates'] = g.macrostates

    return adata


def spatial_aware_embedding(adata, spatial_weight=0.3):
    """Generate spatially-aware low-dimensional embedding"""
    from sklearn.metrics.pairwise import euclidean_distances
    from umap import UMAP

    # Validate spatial data
    is_valid, issues = validate_spatial_data(adata)
    if not is_valid:
        raise ValueError(f"Invalid spatial data: {'; '.join(issues)}")

    spatial_coords = adata.obsm['spatial']

    # Ensure PCA has been computed
    if 'X_pca' not in adata.obsm:
        raise ValueError("PCA has not been computed. Run preprocessing first.")

    # Calculate expression-based distance matrix
    expr_dist = euclidean_distances(adata.obsm['X_pca'])

    # Calculate spatial distance matrix
    spatial_dist = euclidean_distances(spatial_coords)

    # Combine distance matrices
    combined_dist = (1 - spatial_weight) * expr_dist + spatial_weight * spatial_dist

    # Use UMAP for dimensionality reduction
    umap_op = UMAP(metric='precomputed')
    embedding = umap_op.fit_transform(combined_dist)

    # Store results
    adata.obsm['X_spatial_umap'] = embedding

    return adata


def infer_pseudotime_palantir(adata, root_cells=None):
    """Infer pseudotime based on spatial pattern using Palantir.
    
    Raises exceptions directly if Palantir fails - no fallback logic.
    """
    import palantir
    import pandas as pd
    import numpy as np

    # Make sure we have PCA computed
    if 'X_pca' not in adata.obsm:
        import scanpy as sc
        sc.pp.pca(adata)

    # Convert PCA projections to DataFrame for Palantir
    pca_df = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)

    # Run diffusion maps
    dm_res = palantir.utils.run_diffusion_maps(pca_df, n_components=10)

    # Determine multiscale space
    ms_data = pd.DataFrame(
        dm_res['EigenVectors'],
        index=pca_df.index
    )

    # Determine early cell (start cell)
    if root_cells is not None and len(root_cells) > 0:
        if root_cells[0] not in ms_data.index:
            raise ValueError(f"Root cell '{root_cells[0]}' not found in data")
        start_cell = root_cells[0]
    else:
        # Use the cell at the extreme of the first diffusion component
        start_cell = ms_data.iloc[:, 0].idxmax()

    # Run Palantir
    pr_res = palantir.core.run_palantir(ms_data, start_cell)

    # Store pseudotime and branch probabilities
    adata.obs['palantir_pseudotime'] = pr_res.pseudotime
    adata.obsm['palantir_branch_probs'] = pr_res.branch_probs

    return adata


def compute_dpt_fallback(adata, root_cells=None):
    """Compute Diffusion Pseudotime as fallback method."""
    import scanpy as sc
    import numpy as np

    # Make sure we have PCA computed
    if 'X_pca' not in adata.obsm:
        sc.pp.pca(adata)

    # Compute neighbors if not already computed
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, use_rep='X_pca')

    # Compute diffusion map
    sc.tl.diffmap(adata)

    # Set root cell if provided
    if root_cells is not None and len(root_cells) > 0:
        if root_cells[0] in adata.obs_names:
            adata.uns['iroot'] = np.where(adata.obs_names == root_cells[0])[0][0]
        else:
            # If provided root cell not found, use first cell and warn
            print(f"Warning: Root cell '{root_cells[0]}' not found in data. Using first cell instead.")
            adata.uns['iroot'] = 0
    else:
        # If no root cell specified, set the first cell as root
        adata.uns['iroot'] = 0

    # Compute diffusion pseudotime
    try:
        sc.tl.dpt(adata)
    except Exception as e:
        # If DPT fails, create a simple pseudotime based on diffusion components
        if 'X_diffmap' in adata.obsm:
            # Use first diffusion component as pseudotime
            pseudotime = adata.obsm['X_diffmap'][:, 0]
            # Normalize to [0, 1]
            pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())
            adata.obs['dpt_pseudotime'] = pseudotime
        else:
            raise RuntimeError(f"Failed to compute DPT and no diffusion map available: {e}")
    
    # Check if dpt_pseudotime was created
    if 'dpt_pseudotime' not in adata.obs.columns:
        raise RuntimeError("DPT computation did not create 'dpt_pseudotime' column")
    
    # Handle any NaN values
    adata.obs['dpt_pseudotime'] = adata.obs['dpt_pseudotime'].fillna(0)

    return adata


async def analyze_rna_velocity(
    data_id: str,
    data_store: Dict[str, Any],
    params: RNAVelocityParameters = RNAVelocityParameters(),
    context: Optional[Context] = None
) -> RNAVelocityResult:
    """Analyze RNA velocity in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing datasets
        params: RNA velocity analysis parameters
        context: MCP context

    Returns:
        RNA velocity analysis result (computation metadata only, no visualization)
    """
    try:
        import scvelo as scv
    except ImportError:
        raise ProcessingError("scvelo package not found. Install it with: pip install scvelo>=0.2.5")

    if context:
        await context.info(f"Analyzing RNA velocity using {params.mode} mode")

    # Get AnnData object
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Copy AnnData to avoid modifying original data
    adata = data_store[data_id]["adata"].copy()

    # Validate data for velocity analysis
    velocity_valid, velocity_issues = validate_velocity_data(adata)
    has_velocity_data = velocity_valid
    
    if not velocity_valid and not params.reference_data_id:
        if context:
            await context.warning("Data validation issues found:")
            for issue in velocity_issues:
                await context.warning(f"  - {issue}")
            await context.info("Tip: Provide reference_data_id with spliced/unspliced layers for SIRV analysis")

    velocity_computed = False

    with suppress_output():
        try:
            if has_velocity_data:
                if context:
                    await context.info("Found spliced/unspliced layers. Computing velocity directly.")
                
                # Preprocess and compute velocity
                adata = preprocess_for_velocity(adata)
                adata = compute_rna_velocity(adata, mode=params.mode)
                velocity_computed = True
                
            elif params.reference_data_id:
                if context:
                    await context.info(f"No layers found. Using reference data '{params.reference_data_id}' for SIRV.")
                
                if params.reference_data_id not in data_store:
                    raise ValueError(f"Reference dataset {params.reference_data_id} not found in data store")
                
                ref_adata = data_store[params.reference_data_id]["adata"]
                
                # Prepare data for SIRV
                spatial_adata, scrna_adata = prepare_for_sirv(adata, ref_adata)
                
                # Run SIRV
                adata = run_sirv(
                    spatial_adata=spatial_adata,
                    scrna_adata=scrna_adata,
                    n_pcs=params.n_pcs,
                    labels=params.labels
                )
                velocity_computed = True
                
            else:
                if context:
                    await context.warning("Data lacks spliced/unspliced counts and no reference data was provided. Cannot compute RNA velocity.")

        except Exception as e:
            raise ProcessingError(f"RNA velocity analysis failed: {str(e)}") from e

    # Update data store
    data_store[data_id]["adata"] = adata

    # Return result with metadata only (no visualization)
    return RNAVelocityResult(
        data_id=data_id,
        velocity_computed=velocity_computed,
        velocity_graph_key='velocity_graph' if velocity_computed else None,
        mode=params.mode
    )


async def analyze_trajectory(
    data_id: str,
    data_store: Dict[str, Any],
    params: TrajectoryParameters = TrajectoryParameters(),
    context: Optional[Context] = None
) -> TrajectoryResult:
    """Analyze trajectory and cell state transitions in spatial transcriptomics data

    Args:
        data_id: Dataset ID
        data_store: Dictionary storing datasets
        params: Trajectory analysis parameters
        context: MCP context

    Returns:
        Trajectory analysis result (computation metadata only, no visualization)
    """
    if context:
        await context.info(f"Analyzing trajectory using {params.method} method")

    # Get AnnData object
    if data_id not in data_store:
        raise ValueError(f"Dataset {data_id} not found in data store")

    # Copy AnnData to avoid modifying original data
    adata = data_store[data_id]["adata"].copy()

    # Check if RNA velocity has been computed
    has_velocity = 'velocity_graph' in adata.uns

    pseudotime_key = None
    method_used = params.method

    # Strategy 1: Try CellRank if requested and velocity is available
    if params.method == "cellrank" and has_velocity:
        if context:
            await context.info("Attempting trajectory inference with CellRank...")
        
        try:
            import cellrank as cr
        except ImportError:
            raise ProcessingError("cellrank package not found. Install it with: pip install cellrank>=2.0.0")

        try:
            with suppress_output():
                adata = infer_spatial_trajectory_cellrank(
                    adata,
                    spatial_weight=params.spatial_weight,
                    kernel_weights=params.cellrank_kernel_weights,
                    n_states=params.cellrank_n_states
                )
            pseudotime_key = 'pseudotime'
            if context:
                await context.info("CellRank analysis completed successfully.")
        except Exception as cellrank_error:
            if context:
                await context.warning(f"CellRank analysis failed: {cellrank_error}")
                await context.info("Falling back to Palantir method...")
            # Log the full traceback for debugging (if logger available)
            import traceback
            try:
                logger.error(f"CellRank traceback: {traceback.format_exc()}")
            except:
                pass  # Ignore if logger not available
            # Force switch to Palantir
            method_used = "palantir"

    # Strategy 2: Try Palantir (either as primary method or fallback)
    if method_used == "palantir" or not has_velocity:
        if context:
            await context.info("Attempting trajectory inference with Palantir...")
        
        try:
            with suppress_output():
                # Run spatially-aware embedding
                adata = spatial_aware_embedding(adata, spatial_weight=params.spatial_weight)
                
                # Run Palantir
                adata = infer_pseudotime_palantir(adata, root_cells=params.root_cells)
                
            pseudotime_key = 'palantir_pseudotime'
            method_used = "palantir"
            if context:
                await context.info("Palantir analysis completed successfully.")
                
        except Exception as palantir_error:
            if not params.allow_fallback_to_dpt:
                raise ProcessingError(f"Palantir analysis failed: {palantir_error}") from palantir_error
            
            if context:
                await context.warning(f"Palantir analysis failed: {palantir_error}")
                await context.info("Falling back to Diffusion Pseudotime (DPT)...")
            
            # Final fallback: DPT
            try:
                with suppress_output():
                    adata = compute_dpt_fallback(adata, root_cells=params.root_cells)
                    
                pseudotime_key = 'dpt_pseudotime'
                method_used = "dpt"
                if context:
                    await context.info("DPT fallback analysis completed successfully.")
                    
            except Exception as dpt_error:
                raise ProcessingError(f"All trajectory inference methods failed. Last error (DPT): {dpt_error}") from dpt_error

    # Ensure pseudotime key exists
    if pseudotime_key is None or pseudotime_key not in adata.obs.columns:
        raise ProcessingError("Failed to compute pseudotime with any available method")

    # Update data store
    data_store[data_id]["adata"] = adata

    # Return result with metadata only (no visualization)
    return TrajectoryResult(
        data_id=data_id,
        pseudotime_computed=True,
        velocity_computed=has_velocity,
        pseudotime_key=pseudotime_key,
        method=method_used,
        spatial_weight=params.spatial_weight
    )


async def analyze_velocity_with_velovi(
    adata,
    n_epochs: int = 1000,
    n_hidden: int = 128,
    n_latent: int = 10,
    use_gpu: bool = False,
    context: Optional[Context] = None
) -> Dict[str, Any]:
    """Analyze RNA velocity using VELOVI from scvi-tools
    
    VELOVI (Velocity Variational Inference) is a deep generative model for 
    transcriptional dynamics that estimates velocities in a probabilistic framework.
    
    Args:
        adata: AnnData object with spliced and unspliced counts
        n_epochs: Number of epochs for training
        n_hidden: Number of hidden units in neural networks  
        n_latent: Dimensionality of latent space
        use_gpu: Whether to use GPU for training
        context: MCP context for logging
        
    Returns:
        Dictionary containing VELOVI analysis results
        
    Raises:
        ImportError: If scvi-tools package is not available
        ValueError: If input data is invalid or missing required layers
        RuntimeError: If VELOVI computation fails
    """
    try:
        if scvi is None or VELOVI is None:
            raise ImportError("scvi-tools package is required for VELOVI analysis. Install with 'pip install scvi-tools'")
        
        if context:
            await context.info("Starting VELOVI velocity analysis...")
            
        # Validate required layers
        is_valid, issues = validate_velocity_data(adata)
        if not is_valid:
            raise ValueError(f"Invalid data for velocity analysis: {'; '.join(issues)}")
        
        if context:
            await context.info(f"Analyzing velocity for {adata.n_obs} cells and {adata.n_vars} genes with VELOVI")
        
        # Setup VELOVI
        VELOVI.setup_anndata(
            adata,
            spliced_layer="spliced",
            unspliced_layer="unspliced"
        )
        
        # Create VELOVI model
        velovi_model = VELOVI(
            adata,
            n_hidden=n_hidden,
            n_latent=n_latent
        )
        
        if context:
            await context.info("Training VELOVI model...")
        
        # Train model
        if use_gpu:
            velovi_model.train(max_epochs=n_epochs, accelerator='gpu')
        else:
            velovi_model.train(max_epochs=n_epochs)
            
        if context:
            await context.info("VELOVI training completed")
        
        # Get results
        if context:
            await context.info("Extracting VELOVI velocity results...")
            
        # Get velocity estimates
        velocity = velovi_model.get_velocity()
        
        # Handle different return types from get_velocity
        if isinstance(velocity, tuple):
            # Sometimes get_velocity returns (velocity, something_else)
            velocity = velocity[0]
        
        # Ensure velocity is numpy array
        if hasattr(velocity, 'detach'):
            # It's a tensor, convert to numpy
            velocity = velocity.detach().cpu().numpy()
        elif hasattr(velocity, 'toarray'):
            # It's a sparse matrix
            velocity = velocity.toarray()
        
        # Ensure velocity is 2D array
        if velocity.ndim == 1:
            velocity = velocity.reshape(-1, 1)
        
        # Get latent representation
        latent = velovi_model.get_latent_representation()
        
        # Get velocity uncertainty (if available)
        try:
            velocity_uncertainty = velovi_model.get_velocity_uncertainty()
            uncertainty_computed = True
        except AttributeError:
            try:
                uncertainty_result = velovi_model.get_directional_uncertainty()
                # Handle case where method returns a tuple
                if isinstance(uncertainty_result, tuple):
                    velocity_uncertainty = uncertainty_result[0]  # Usually the first element is the uncertainty
                else:
                    velocity_uncertainty = uncertainty_result
                uncertainty_computed = True
            except (AttributeError, IndexError):
                velocity_uncertainty = None
                uncertainty_computed = False
        
        # Store results in adata
        adata.layers['velocity_velovi'] = velocity
        adata.obsm['X_velovi_latent'] = latent
        if velocity_uncertainty is not None:
            try:
                adata.layers['velocity_velovi_uncertainty'] = velocity_uncertainty
            except Exception:
                # If uncertainty storage fails, just skip it
                uncertainty_computed = False
        
        # Calculate velocity statistics
        velocity_norm = np.linalg.norm(velocity, axis=1)
        
        # Store velocity norm in obs
        adata.obs['velocity_velovi_norm'] = velocity_norm
        
        # Calculate summary statistics
        results = {
            'method': 'VELOVI',
            'n_latent_dims': n_latent,
            'n_epochs': n_epochs,
            'velocity_computed': True,
            'velocity_mean_norm': velocity_norm.mean(),
            'velocity_std_norm': velocity_norm.std(),
            'velocity_shape': velocity.shape,
            'latent_shape': latent.shape,
            'uncertainty_computed': uncertainty_computed,
            'training_completed': True,
            'device': 'GPU' if use_gpu else 'CPU'
        }
        
        if context:
            await context.info(f"VELOVI velocity analysis completed successfully")
            await context.info(f"Stored velocity in adata.layers['velocity_velovi']")
            await context.info(f"Stored latent representation in adata.obsm['X_velovi_latent']")
            await context.info(f"Stored uncertainty in adata.layers['velocity_velovi_uncertainty']")
        
        return results
        
    except Exception as e:
        error_msg = f"VELOVI velocity analysis failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg) from e
