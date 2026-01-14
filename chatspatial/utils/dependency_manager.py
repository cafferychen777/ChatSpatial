"""
Unified dependency management for ChatSpatial MCP.

Provides a consistent API for managing optional dependencies, replacing
scattered try/except ImportError patterns with centralized handling.

Usage:
    # Require a dependency (raises if missing)
    scvi = require("scvi-tools", feature="cell type annotation")

    # Get optional dependency (returns None if missing)
    torch = get("torch")

    # Check availability
    if is_available("rpy2"):
        import rpy2
"""

import importlib
import importlib.util
import warnings
from dataclasses import dataclass
from functools import lru_cache
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext


@dataclass(frozen=True)
class DependencyInfo:
    """Metadata for an optional dependency."""

    module_name: str
    install_cmd: str
    description: str = ""


# Registry of optional dependencies with install instructions
DEPENDENCY_REGISTRY: dict[str, DependencyInfo] = {
    # Deep Learning
    "scvi-tools": DependencyInfo(
        "scvi", "pip install scvi-tools", "Single-cell variational inference tools"
    ),
    "torch": DependencyInfo(
        "torch", "pip install torch", "PyTorch deep learning framework"
    ),
    "cell2location": DependencyInfo(
        "cell2location",
        "pip install cell2location",
        "Probabilistic cell type deconvolution",
    ),
    "flashdeconv": DependencyInfo(
        "flashdeconv", "pip install flashdeconv", "Ultra-fast spatial deconvolution"
    ),
    # Spatial Analysis
    "tangram": DependencyInfo(
        "tangram",
        "pip install tangram-sc",
        "Spatial mapping of single-cell transcriptomics",
    ),
    "squidpy": DependencyInfo(
        "squidpy", "pip install squidpy", "Spatial single-cell analysis"
    ),
    "SpaGCN": DependencyInfo(
        "SpaGCN",
        "pip install SpaGCN",
        "Spatial domain identification using graph convolutional networks",
    ),
    "STAGATE": DependencyInfo(
        "STAGATE_pyG",
        "pip install STAGATE-pyG",
        "Spatial domain identification using graph attention",
    ),
    "GraphST": DependencyInfo(
        "GraphST",
        "pip install GraphST",
        "Graph self-supervised contrastive learning for spatial domains",
    ),
    "paste": DependencyInfo(
        "paste",
        "pip install paste-bio",
        "Probabilistic alignment of spatial transcriptomics",
    ),
    "stalign": DependencyInfo(
        "STalign", "pip install STalign", "Spatial transcriptomics alignment"
    ),
    # R Interface
    "rpy2": DependencyInfo(
        "rpy2", "pip install rpy2", "R-Python interface (requires R installation)"
    ),
    "anndata2ri": DependencyInfo(
        "anndata2ri",
        "pip install anndata2ri",
        "AnnData to R SingleCellExperiment conversion",
    ),
    # Cell Communication
    "liana": DependencyInfo(
        "liana", "pip install liana", "Ligand-receptor analysis framework"
    ),
    "cellphonedb": DependencyInfo(
        "cellphonedb",
        "pip install cellphonedb",
        "Statistical method for cell-cell communication",
    ),
    "ktplotspy": DependencyInfo(
        "ktplotspy", "pip install ktplotspy", "CellPhoneDB visualization toolkit"
    ),
    # RNA Velocity
    "scvelo": DependencyInfo("scvelo", "pip install scvelo", "RNA velocity analysis"),
    "velovi": DependencyInfo(
        "velovi", "pip install velovi", "Variational inference for RNA velocity"
    ),
    "cellrank": DependencyInfo(
        "cellrank", "pip install cellrank", "Trajectory inference using RNA velocity"
    ),
    "palantir": DependencyInfo(
        "palantir", "pip install palantir", "Trajectory inference for cell fate"
    ),
    # Annotation
    "singler": DependencyInfo(
        "singler",
        "pip install singler singlecellexperiment",
        "Reference-based cell type annotation",
    ),
    "mllmcelltype": DependencyInfo(
        "mllmcelltype", "pip install mllmcelltype", "LLM-based cell type annotation"
    ),
    "celldex": DependencyInfo(
        "celldex", "pip install celldex", "Cell type reference datasets for SingleR"
    ),
    # Enrichment
    "gseapy": DependencyInfo(
        "gseapy", "pip install gseapy", "Gene set enrichment analysis"
    ),
    "decoupler": DependencyInfo(
        "decoupler", "pip install decoupler", "Functional analysis of omics data"
    ),
    # Spatial Statistics
    "sparkx": DependencyInfo(
        "sparkx", "pip install SPARK-X", "SPARK-X non-parametric spatial gene detection"
    ),
    "spatialde": DependencyInfo(
        "NaiveDE",
        "pip install SpatialDE",
        "SpatialDE Gaussian process spatial gene detection",
    ),
    # CNV
    "infercnvpy": DependencyInfo(
        "infercnvpy", "pip install infercnvpy", "Copy number variation inference"
    ),
    # Visualization
    "plotly": DependencyInfo(
        "plotly", "pip install plotly", "Interactive visualization"
    ),
    "adjustText": DependencyInfo(
        "adjustText", "pip install adjustText", "Text label placement for matplotlib"
    ),
    "splot": DependencyInfo("splot", "pip install splot", "Spatial plotting for PySAL"),
    # Data handling
    "mudata": DependencyInfo(
        "mudata", "pip install mudata", "Multimodal data handling"
    ),
    # Integration
    "harmonypy": DependencyInfo(
        "harmonypy", "pip install harmonypy", "Harmony batch integration"
    ),
    "scanorama": DependencyInfo(
        "scanorama", "pip install scanorama", "Scanorama batch integration"
    ),
    "bbknn": DependencyInfo(
        "bbknn", "pip install bbknn", "Batch balanced k-nearest neighbors"
    ),
    # Spatial weights
    "esda": DependencyInfo(
        "esda", "pip install esda", "Exploratory spatial data analysis"
    ),
    "libpysal": DependencyInfo(
        "libpysal", "pip install libpysal", "Python spatial analysis library"
    ),
    # Other
    "dask": DependencyInfo("dask", "pip install dask", "Parallel computing library"),
    "ot": DependencyInfo("ot", "pip install POT", "Python Optimal Transport library"),
    "louvain": DependencyInfo(
        "louvain", "pip install louvain", "Louvain community detection algorithm"
    ),
    "pydeseq2": DependencyInfo(
        "pydeseq2", "pip install pydeseq2", "Python implementation of DESeq2"
    ),
    "enrichmap": DependencyInfo(
        "enrichmap", "pip install enrichmap", "Spatial enrichment mapping"
    ),
    "pygam": DependencyInfo(
        "pygam", "pip install pygam", "Generalized additive models"
    ),
    "skgstat": DependencyInfo(
        "skgstat", "pip install scikit-gstat", "Geostatistical analysis toolkit"
    ),
    "sklearn": DependencyInfo(
        "sklearn", "pip install scikit-learn", "Machine learning library"
    ),
    "statsmodels": DependencyInfo(
        "statsmodels", "pip install statsmodels", "Statistical models and tests"
    ),
    "scipy": DependencyInfo(
        "scipy", "pip install scipy", "Scientific computing library"
    ),
    "scanpy": DependencyInfo(
        "scanpy", "pip install scanpy", "Single-cell analysis in Python"
    ),
    "Pillow": DependencyInfo("PIL", "pip install Pillow", "Python Imaging Library"),
}


# =============================================================================
# Core Functions (using @lru_cache for thread-safe caching)
# =============================================================================


def _get_info(name: str) -> DependencyInfo:
    """Get dependency info, creating default if not in registry."""
    if name in DEPENDENCY_REGISTRY:
        return DEPENDENCY_REGISTRY[name]
    # Check by module name
    for info in DEPENDENCY_REGISTRY.values():
        if info.module_name == name:
            return info
    # Default for unknown dependencies
    return DependencyInfo(name, f"pip install {name}", f"Optional: {name}")


@lru_cache(maxsize=256)
def _try_import(module_name: str) -> Optional[Any]:
    """Import module with caching. Returns None if unavailable."""
    try:
        return importlib.import_module(module_name)
    except ImportError:
        return None


@lru_cache(maxsize=256)
def _check_spec(module_name: str) -> bool:
    """Fast availability check without importing."""
    return importlib.util.find_spec(module_name) is not None


# =============================================================================
# Public API
# =============================================================================


def is_available(name: str) -> bool:
    """Check if a dependency is available (fast, no import)."""
    return _check_spec(_get_info(name).module_name)


def get(
    name: str,
    ctx: Optional["ToolContext"] = None,
    warn_if_missing: bool = False,
) -> Optional[Any]:
    """Get optional dependency, returning None if unavailable."""
    info = _get_info(name)
    module = _try_import(info.module_name)

    if module is not None:
        return module

    if warn_if_missing:
        msg = f"{name} not available. Install: {info.install_cmd}"
        warnings.warn(msg, stacklevel=2)

    return None


def require(
    name: str,
    ctx: Optional["ToolContext"] = None,
    feature: Optional[str] = None,
) -> Any:
    """Require a dependency, raising ImportError if unavailable."""
    info = _get_info(name)
    module = _try_import(info.module_name)

    if module is not None:
        return module

    feature_msg = f" for {feature}" if feature else ""
    raise ImportError(
        f"{name} is required{feature_msg}.\n\n"
        f"Install: {info.install_cmd}\n"
        f"Description: {info.description}"
    )


# =============================================================================
# R Environment Validation
# =============================================================================


def validate_r_environment(
    ctx: Optional["ToolContext"] = None,
    required_packages: Optional[list[str]] = None,
) -> tuple[Any, ...]:
    """Validate R environment and return required modules.

    Returns:
        Tuple of (robjects, pandas2ri, numpy2ri, importr, localconverter,
                  default_converter, openrlib, anndata2ri)

    Raises:
        ImportError: If rpy2 or required R packages are not available
    """
    if not is_available("rpy2"):
        raise ImportError(
            "rpy2 is required for R-based methods. "
            "Install: pip install rpy2 (requires R installation)"
        )
    if not is_available("anndata2ri"):
        raise ImportError(
            "anndata2ri is required for R-based methods. "
            "Install: pip install anndata2ri"
        )

    try:
        import anndata2ri
        import rpy2.robjects as robjects
        from rpy2.rinterface_lib import openrlib
        from rpy2.robjects import conversion, default_converter, numpy2ri, pandas2ri
        from rpy2.robjects.conversion import localconverter
        from rpy2.robjects.packages import importr

        # Test R availability
        with openrlib.rlock:
            with conversion.localconverter(default_converter):
                robjects.r("R.version")

        # Check required R packages
        if required_packages:
            missing = []
            for pkg in required_packages:
                try:
                    with openrlib.rlock:
                        with conversion.localconverter(default_converter):
                            importr(pkg)
                except Exception:
                    missing.append(pkg)

            if missing:
                pkg_list = ", ".join(f"'{p}'" for p in missing)
                raise ImportError(
                    f"Missing R packages: {pkg_list}\n"
                    f"Install in R: install.packages(c({pkg_list}))"
                )

        return (
            robjects,
            pandas2ri,
            numpy2ri,
            importr,
            localconverter,
            default_converter,
            openrlib,
            anndata2ri,
        )

    except ImportError:
        raise
    except Exception as e:
        raise ImportError(
            f"R environment setup failed: {e}\n\n"
            "Solutions:\n"
            "  - Install R: https://www.r-project.org/\n"
            "  - Set R_HOME environment variable\n"
            "  - macOS: brew install r\n"
            "  - Ubuntu: sudo apt install r-base"
        ) from e


def validate_r_package(
    package_name: str,
    ctx: Optional["ToolContext"] = None,
    install_cmd: Optional[str] = None,
) -> bool:
    """Check if an R package is available."""
    if not is_available("rpy2"):
        raise ImportError(
            "rpy2 is required for R-based methods.\n"
            "Install: pip install rpy2 (requires R)"
        )

    try:
        from rpy2.rinterface_lib import openrlib
        from rpy2.robjects import conversion, default_converter
        from rpy2.robjects.packages import importr

        with openrlib.rlock:
            with conversion.localconverter(default_converter):
                importr(package_name)

        return True

    except Exception as e:
        install = install_cmd or f"install.packages('{package_name}')"
        raise ImportError(
            f"R package '{package_name}' not installed.\n" f"Install in R: {install}"
        ) from e


def check_r_packages(
    packages: list[str],
    ctx: Optional["ToolContext"] = None,
) -> list[str]:
    """Check availability of multiple R packages. Returns missing ones."""
    if not is_available("rpy2"):
        return packages

    missing = []
    for pkg in packages:
        try:
            validate_r_package(pkg)
        except ImportError:
            missing.append(pkg)

    return missing


def validate_scvi_tools(
    ctx: Optional["ToolContext"] = None,
    components: Optional[list[str]] = None,
) -> Any:
    """Validate scvi-tools availability and return the module."""
    scvi = require("scvi-tools", ctx, "scvi-tools methods")

    if components:
        missing = []
        for comp in components:
            try:
                if comp == "CellAssign":
                    from scvi.external import CellAssign  # noqa: F401
                elif comp == "Cell2location":
                    import cell2location  # noqa: F401
                elif comp == "SCANVI":
                    from scvi.model import SCANVI  # noqa: F401
                elif comp == "DestVI":
                    from scvi.external import DestVI  # noqa: F401
                elif comp == "Stereoscope":
                    from scvi.external import Stereoscope  # noqa: F401
                else:
                    getattr(scvi, comp, None) or getattr(
                        scvi.model, comp, None
                    ) or getattr(scvi.external, comp, None)
            except (ImportError, AttributeError):
                missing.append(comp)

        if missing:
            raise ImportError(
                f"scvi-tools components not available: {', '.join(missing)}\n"
                "Try: pip install --upgrade scvi-tools"
            )

    return scvi
