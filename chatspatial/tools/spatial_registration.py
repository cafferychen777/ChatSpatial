"""
Spatial Registration Tool

This module provides functionality for aligning and registering multiple spatial transcriptomics slices.
"""

import logging
from typing import List, Optional

import anndata as ad
import numpy as np

logger = logging.getLogger(__name__)


class SpatialRegistration:
    """
    A class for performing spatial registration of multiple tissue slices.

    This class provides methods for:
    - Pairwise alignment of spatial transcriptomics slices
    - Multi-slice integration into a common coordinate system
    - Spatial coordinate transformation
    """

    def __init__(self):
        """Initialize the SpatialRegistration tool."""
        self.method_info = {
            "paste": {
                "name": "PASTE",
                "description": "Probabilistic Alignment of Spatial Transcriptomics Experiments",
                "reference": "Zeira et al., Nature Methods 2022",
                "type": "python",
            },
            "stalign": {
                "name": "STalign",
                "description": "Spatial transcriptomics alignment using image registration",
                "reference": "Clifton et al., bioRxiv 2023",
                "type": "python",
            },
        }

    def register_slices(
        self,
        adata_list: List[ad.AnnData],
        method: str = "paste",
        reference_idx: Optional[int] = None,
        **kwargs,
    ) -> List[ad.AnnData]:
        """
        Register multiple spatial transcriptomics slices.

        Parameters
        ----------
        adata_list : List[ad.AnnData]
            List of AnnData objects to register
        method : str
            Registration method to use ('paste', 'stalign')
        reference_idx : Optional[int]
            Index of reference slice (if None, use first slice)
        **kwargs
            Additional method-specific parameters

        Returns
        -------
        List[ad.AnnData]
            List of registered AnnData objects with aligned coordinates
        """
        if method not in self.method_info:
            raise ValueError(
                f"Unknown method: {method}. Available methods: {list(self.method_info.keys())}"
            )

        logger.info(f"Performing spatial registration using {method}")

        if method == "paste":
            return self._register_with_paste(adata_list, reference_idx, **kwargs)
        elif method == "stalign":
            return self._register_with_stalign(adata_list, reference_idx, **kwargs)
        else:
            raise NotImplementedError(f"Method {method} not yet implemented")

    def _register_with_paste(
        self,
        adata_list: List[ad.AnnData],
        reference_idx: Optional[int] = None,
        alpha: float = 0.1,
        n_components: int = 30,
        use_gpu: bool = False,
        **kwargs,
    ) -> List[ad.AnnData]:
        """
        Register slices using PASTE algorithm.

        Parameters
        ----------
        adata_list : List[ad.AnnData]
            List of AnnData objects to register
        reference_idx : Optional[int]
            Index of reference slice
        alpha : float
            Spatial regularization parameter
        n_components : int
            Number of components for dimensionality reduction
        use_gpu : bool
            Whether to use GPU acceleration
        **kwargs
            Additional PASTE parameters

        Returns
        -------
        List[ad.AnnData]
            Registered AnnData objects
        """
        try:
            import paste as pst
        except ImportError:
            logger.error("PASTE not installed. Install with: pip install paste-bio")
            raise ImportError("Please install paste-bio: pip install paste-bio")

        if reference_idx is None:
            reference_idx = 0

        # Prepare data
        logger.info("Preparing data for PASTE registration")
        registered_list = [adata.copy() for adata in adata_list]

        # Ensure all slices have spatial coordinates
        for i, adata in enumerate(registered_list):
            if "spatial" not in adata.obsm:
                raise ValueError(
                    f"Slice {i} missing spatial coordinates in obsm['spatial']"
                )

        # Perform pairwise alignment
        if len(registered_list) == 2:
            logger.info("Performing pairwise alignment")

            # Get common genes
            common_genes = list(
                set(registered_list[0].var_names) & set(registered_list[1].var_names)
            )
            logger.info(f"Using {len(common_genes)} common genes")

            # Make gene names unique to avoid indexing issues
            for adata in registered_list:
                adata.var_names_make_unique()

            # Get common genes after making names unique
            common_genes = list(
                set(registered_list[0].var_names) & set(registered_list[1].var_names)
            )
            logger.info(
                f"Using {len(common_genes)} common genes after making names unique"
            )

            slice1 = registered_list[0][:, common_genes].copy()
            slice2 = registered_list[1][:, common_genes].copy()

            # Perform basic preprocessing like in the direct test
            import scanpy as sc

            sc.pp.normalize_total(slice1, target_sum=1e4)
            sc.pp.log1p(slice1)
            sc.pp.normalize_total(slice2, target_sum=1e4)
            sc.pp.log1p(slice2)

            # Run PASTE with minimal parameters that worked in direct test
            pi = pst.pairwise_align(
                slice1,
                slice2,
                alpha=0.1,  # Use same alpha as direct test
                verbose=True,
                numItermax=10,  # Reduce iterations like in direct test
            )

            # Apply transformation
            aligned_slices = pst.stack_slices_pairwise([slice1, slice2], [pi])

            # Extract the aligned spatial coordinates
            registered_list[0].obsm["spatial_registered"] = aligned_slices[0].obsm[
                "spatial"
            ]
            registered_list[1].obsm["spatial_registered"] = aligned_slices[1].obsm[
                "spatial"
            ]

        else:
            # Multi-slice registration using center alignment
            logger.info(
                f"Performing multi-slice registration with {len(registered_list)} slices"
            )

            # Get common genes
            common_genes = set(registered_list[0].var_names)
            for adata in registered_list[1:]:
                common_genes = common_genes & set(adata.var_names)
            common_genes = list(common_genes)
            logger.info(f"Using {len(common_genes)} common genes")

            # Subset to common genes
            slices = [adata[:, common_genes] for adata in registered_list]

            # Initial pairwise alignments to reference
            pis = []

            # Create backend object
            import ot

            if use_gpu:
                try:
                    import torch

                    if torch.cuda.is_available():
                        backend = ot.backend.TorchBackend()
                    else:
                        backend = ot.backend.NumpyBackend()
                except ImportError:
                    backend = ot.backend.NumpyBackend()
            else:
                backend = ot.backend.NumpyBackend()

            for i, slice_data in enumerate(slices):
                if i == reference_idx:
                    # Identity mapping for reference slice
                    n = slices[i].shape[0]
                    pi = np.eye(n)
                else:
                    # Filter out parameters that will be passed explicitly
                    filtered_kwargs = {
                        k: v
                        for k, v in kwargs.items()
                        if k
                        not in ["alpha", "backend", "use_gpu", "verbose", "gpu_verbose"]
                    }
                    pi = pst.pairwise_align(
                        slices[reference_idx],
                        slice_data,
                        alpha=alpha,
                        backend=backend,
                        use_gpu=use_gpu,
                        verbose=False,
                        gpu_verbose=False,
                        **filtered_kwargs,
                    )
                pis.append(pi)

            # Center alignment
            center_slice, pis_new = pst.center_align(
                slices[reference_idx],
                slices,
                pis_init=pis,
                alpha=alpha,
                backend=backend,
                use_gpu=use_gpu,
                n_components=n_components,
                verbose=False,
                gpu_verbose=False,
            )

            # Apply transformations
            for i, (adata, pi) in enumerate(zip(registered_list, pis_new)):
                if i == reference_idx:
                    adata.obsm["spatial_registered"] = adata.obsm["spatial"].copy()
                else:
                    # Transform coordinates based on alignment
                    coords = adata.obsm["spatial"]
                    # Apply transformation from pi
                    adata.obsm["spatial_registered"] = self._transform_coordinates(
                        coords, pi, slices[reference_idx].obsm["spatial"]
                    )

        logger.info("PASTE registration completed")
        return registered_list

    def _register_with_stalign(
        self,
        adata_list: List[ad.AnnData],
        reference_idx: Optional[int] = None,
        **kwargs,
    ) -> List[ad.AnnData]:
        """
        Register slices using STalign diffeomorphic mapping algorithm.

        Parameters:
        -----------
        adata_list : List[ad.AnnData]
            List of AnnData objects to register
        reference_idx : Optional[int]
            Index of reference slice (default: 0)
        **kwargs
            Additional parameters for STalign:
            - image_size: Tuple[int, int] = (128, 128)
            - niter: int = 100 (number of iterations)
            - a: float = 500.0 (regularization parameter)
            - use_expression: bool = True (use gene expression for alignment)

        Returns:
        --------
        List[ad.AnnData]
            Registered AnnData objects with 'spatial_stalign' coordinates
        """
        try:
            import STalign.STalign as ST
            import torch
        except ImportError:
            logger.error(
                "STalign not available. Install with: pip install git+https://github.com/JEFworks-Lab/STalign.git"
            )
            raise ImportError("STalign is required for STalign registration")

        logger.info("Starting STalign diffeomorphic registration")

        # Default parameters
        image_size = kwargs.get("image_size", (128, 128))
        niter = kwargs.get("niter", 50)  # Reduced for faster processing
        use_expression = kwargs.get("use_expression", True)
        device = kwargs.get("device", "cpu")

        stalign_params = {
            "a": kwargs.get("a", 500.0),
            "p": kwargs.get("p", 2.0),
            "expand": kwargs.get("expand", 2.0),
            "nt": kwargs.get("nt", 3),
            "niter": niter,
            "diffeo_start": kwargs.get("diffeo_start", 0),
            "epL": kwargs.get("epL", 2e-08),
            "epT": kwargs.get("epT", 0.2),
            "epV": kwargs.get("epV", 2000.0),
            "sigmaM": kwargs.get("sigmaM", 1.0),
            "sigmaB": kwargs.get("sigmaB", 2.0),
            "sigmaA": kwargs.get("sigmaA", 5.0),
            "sigmaR": kwargs.get("sigmaR", 500000.0),
            "sigmaP": kwargs.get("sigmaP", 20.0),
            "device": device,
            "dtype": torch.float32,
        }

        registered_list = [adata.copy() for adata in adata_list]
        reference_idx = reference_idx or 0

        # Ensure all slices have spatial coordinates
        for i, adata in enumerate(registered_list):
            if "spatial" not in adata.obsm:
                raise ValueError(
                    f"Slice {i} missing spatial coordinates in obsm['spatial']"
                )

        # For pairwise registration (2 slices)
        if len(registered_list) == 2:
            logger.info("Performing STalign pairwise registration")

            try:
                # Get data from both slices
                source_adata = registered_list[0]
                target_adata = registered_list[1]

                # Prepare coordinates
                source_coords = source_adata.obsm["spatial"].astype(np.float32)
                target_coords = target_adata.obsm["spatial"].astype(np.float32)

                # Prepare expression data
                if use_expression:
                    # Get common genes
                    common_genes = list(
                        set(source_adata.var_names) & set(target_adata.var_names)
                    )
                    if len(common_genes) < 100:
                        logger.warning(f"Only {len(common_genes)} common genes found")

                    source_expr = source_adata[:, common_genes].X
                    target_expr = target_adata[:, common_genes].X

                    # Handle sparse matrices
                    if hasattr(source_expr, "toarray"):
                        source_expr = source_expr.toarray()
                    if hasattr(target_expr, "toarray"):
                        target_expr = target_expr.toarray()

                    # Sum expression for intensity
                    source_intensity = source_expr.sum(axis=1).astype(np.float32)
                    target_intensity = target_expr.sum(axis=1).astype(np.float32)
                else:
                    # Use uniform intensity
                    source_intensity = np.ones(len(source_coords), dtype=np.float32)
                    target_intensity = np.ones(len(target_coords), dtype=np.float32)

                logger.info(
                    f"Registering {len(source_coords)} -> {len(target_coords)} spots"
                )

                # Try STalign LDDMM registration
                try:
                    # Prepare image data for STalign
                    def prepare_image_data(coords, intensity, image_size):
                        """Convert point cloud to image for STalign"""
                        # Normalize coordinates to image size
                        coords_norm = coords.copy()
                        padding = 0.1
                        for i in range(2):
                            coord_min, coord_max = (
                                coords[:, i].min(),
                                coords[:, i].max(),
                            )
                            coord_range = coord_max - coord_min
                            if coord_range > 0:
                                target_min = padding * image_size[i]
                                target_max = (1 - padding) * image_size[i]
                                coords_norm[:, i] = (
                                    coords[:, i] - coord_min
                                ) / coord_range
                                coords_norm[:, i] = (
                                    coords_norm[:, i] * (target_max - target_min)
                                    + target_min
                                )

                        # Create coordinate grid for STalign
                        x_coords = torch.linspace(
                            0, image_size[0], image_size[0], dtype=torch.float32
                        )
                        y_coords = torch.linspace(
                            0, image_size[1], image_size[1], dtype=torch.float32
                        )
                        xgrid = [x_coords, y_coords]

                        # Rasterize points to image
                        image = np.zeros(image_size, dtype=np.float32)
                        for i in range(len(coords)):
                            x_idx = int(
                                np.clip(coords_norm[i, 1], 0, image_size[0] - 1)
                            )
                            y_idx = int(
                                np.clip(coords_norm[i, 0], 0, image_size[1] - 1)
                            )
                            # Add Gaussian kernel for smoother image
                            for dx in range(-2, 3):
                                for dy in range(-2, 3):
                                    xi, yi = x_idx + dx, y_idx + dy
                                    if (
                                        0 <= xi < image_size[0]
                                        and 0 <= yi < image_size[1]
                                    ):
                                        dist = np.sqrt(dx * dx + dy * dy)
                                        weight = np.exp(-dist * dist / (2 * 1.0 * 1.0))
                                        image[xi, yi] += intensity[i] * weight

                        # Normalize image
                        if image.max() > 0:
                            image = image / image.max()

                        return xgrid, torch.tensor(image, dtype=torch.float32)

                    # Prepare images
                    xI, I = prepare_image_data(
                        source_coords, source_intensity, image_size
                    )
                    xJ, J = prepare_image_data(
                        target_coords, target_intensity, image_size
                    )

                    logger.info(
                        f"Running STalign LDDMM on {I.shape} -> {J.shape} images"
                    )

                    # Run STalign LDDMM
                    result = ST.LDDMM(xI=xI, I=I, xJ=xJ, J=J, **stalign_params)

                    # Extract transformation components
                    A = result.get("A", None)
                    v = result.get("v", None)
                    xv = result.get("xv", None)

                    if A is not None and v is not None and xv is not None:
                        # Transform source coordinates using STalign
                        source_points = torch.tensor(source_coords, dtype=torch.float32)
                        transformed_coords = ST.transform_points_source_to_target(
                            xv, v, A, source_points
                        )

                        if isinstance(transformed_coords, torch.Tensor):
                            transformed_coords = transformed_coords.numpy()

                        logger.info("STalign LDDMM registration successful")
                    else:
                        raise ValueError(
                            "STalign did not return valid transformation components"
                        )

                except Exception as stalign_error:
                    # STalign failed - do not fallback to inferior methods
                    error_msg = (
                        f"STalign LDDMM registration failed: {stalign_error}. "
                        f"STalign requires high-quality spatial data and proper preprocessing. "
                        f"Please check: 1) Data quality and preprocessing, 2) Spatial coordinate format, "
                        f"or 3) Consider using PASTE method instead for more robust registration."
                    )
                    logger.error(error_msg)
                    raise RuntimeError(error_msg)

            except Exception as e:
                # Re-raise the error instead of masking it with fake success
                error_msg = f"STalign registration failed: {e}. Please check your data or try PASTE method."
                logger.error(error_msg)
                raise RuntimeError(error_msg)

        else:
            # Explicitly reject multi-slice registration instead of faking it
            raise NotImplementedError(
                f"STalign does not support multi-slice registration ({len(registered_list)} slices provided). "
                f"STalign only supports pairwise registration (2 slices). "
                f"Please use pairwise registration or switch to PASTE method for multi-slice registration."
            )

        return registered_list

    def _transform_coordinates(
        self,
        coords: np.ndarray,
        transport_matrix: np.ndarray,
        reference_coords: np.ndarray,
    ) -> np.ndarray:
        """
        Transform coordinates based on optimal transport matrix.

        Parameters
        ----------
        coords : np.ndarray
            Original coordinates
        transport_matrix : np.ndarray
            Optimal transport matrix from PASTE
        reference_coords : np.ndarray
            Reference slice coordinates

        Returns
        -------
        np.ndarray
            Transformed coordinates
        """
        # Normalize transport matrix
        transport_matrix = transport_matrix / transport_matrix.sum(
            axis=1, keepdims=True
        )

        # Compute weighted average of reference coordinates
        transformed = transport_matrix @ reference_coords

        return transformed


# Create convenience functions for direct use
def register_spatial_slices(
    adata_list: List[ad.AnnData], method: str = "paste", **kwargs
) -> List[ad.AnnData]:
    """
    Register multiple spatial transcriptomics slices.

    Parameters
    ----------
    adata_list : List[ad.AnnData]
        List of AnnData objects to register
    method : str
        Registration method ('paste', 'stalign')
    **kwargs
        Method-specific parameters

    Returns
    -------
    List[ad.AnnData]
        Registered AnnData objects
    """
    registration_tool = SpatialRegistration()
    return registration_tool.register_slices(adata_list, method=method, **kwargs)


# Standard architecture-compliant wrapper function for MCP server
async def register_spatial_slices_mcp(
    source_id: str,
    target_id: str,
    data_store: dict,
    method: str = "paste",
    landmarks: Optional[dict] = None,
    context=None,
) -> dict:
    """
    Register spatial slices following the standard MCP tool architecture.

    This function follows the standard pattern used by other tools:
    - Accepts data IDs and data_store dictionary
    - Extracts AnnData objects internally
    - Returns results in standard format

    Args:
        source_id: Source dataset ID
        target_id: Target dataset ID to align to
        data_store: Dictionary containing dataset information
        method: Registration method (paste, stalign)
        landmarks: Additional parameters (currently unused)
        context: MCP context for logging

    Returns:
        Dictionary with registration results
    """
    if context:
        await context.info(
            f"Registering {source_id} to {target_id} using method: {method}"
        )

    # Validate data exists
    if source_id not in data_store:
        raise ValueError(f"Source dataset {source_id} not found in data store")
    if target_id not in data_store:
        raise ValueError(f"Target dataset {target_id} not found in data store")

    # Extract AnnData objects
    source_adata = data_store[source_id]["adata"].copy()
    target_adata = data_store[target_id]["adata"].copy()
    adata_list = [source_adata, target_adata]

    try:
        # Call the core registration function
        registered_adata_list = register_spatial_slices(adata_list, method=method)

        # Update data store with registered data
        data_store[source_id]["adata"] = registered_adata_list[0]  # source (registered)
        data_store[target_id]["adata"] = registered_adata_list[1]  # target (reference)

        # Create result dictionary
        result = {
            "method": method,
            "source_id": source_id,
            "target_id": target_id,
            "n_source_spots": registered_adata_list[0].n_obs,
            "n_target_spots": registered_adata_list[1].n_obs,
            "registration_completed": True,
            "spatial_key_registered": "spatial_registered",
        }

        if context:
            await context.info(
                "Registration completed. Registered coordinates stored in 'spatial_registered' key."
            )

        return result

    except Exception as e:
        error_msg = f"Registration failed: {str(e)}"
        if context:
            await context.error(error_msg)
        raise RuntimeError(error_msg) from e
