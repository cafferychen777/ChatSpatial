"""
Unified dependency management for ChatSpatial MCP.

This module provides a consistent API for managing optional dependencies across
all tool modules, replacing the scattered try/except ImportError patterns with
a centralized, type-safe, and well-documented approach.

Design principles:
1. Lazy loading: Dependencies are only imported when first requested
2. Clear error messages: Each dependency has install instructions
3. MCP context integration: Logging support for debugging
4. Thread-safe: Safe for concurrent access
5. Type-safe: Proper typing for IDE support

Usage examples:
    # For required dependencies (raises if missing)
    scvi = deps.require("scvi-tools", ctx=ctx)

    # For optional dependencies (returns None if missing)
    torch = deps.get("torch")

    # For boolean checks
    if deps.is_available("rpy2"):
        import rpy2
"""

import importlib
import importlib.util
import threading
import warnings
from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext


class DependencyCategory(Enum):
    """Categories for organizing dependencies."""

    CORE = "core"  # Essential dependencies (scanpy, anndata)
    DEEP_LEARNING = "deep_learning"  # scvi-tools, torch, tensorflow
    SPATIAL = "spatial"  # squidpy, tangram, SpaGCN
    R_INTERFACE = "r_interface"  # rpy2, anndata2ri
    COMMUNICATION = "communication"  # liana, cellphonedb
    VELOCITY = "velocity"  # scvelo, velovi
    VISUALIZATION = "visualization"  # matplotlib, seaborn


@dataclass
class DependencyInfo:
    """Information about an optional dependency.

    Attributes:
        module_name: The actual Python module name to import
        install_cmd: pip/conda install command
        description: Human-readable description
        category: Dependency category for grouping
        min_version: Minimum required version (optional)
        alternatives: Alternative packages that provide same functionality
    """

    module_name: str
    install_cmd: str
    description: str = ""
    category: DependencyCategory = DependencyCategory.CORE
    min_version: Optional[str] = None
    alternatives: List[str] = field(default_factory=list)


# Registry of all optional dependencies with metadata
# This is the single source of truth for dependency information
DEPENDENCY_REGISTRY: Dict[str, DependencyInfo] = {
    # Deep Learning Framework
    "scvi-tools": DependencyInfo(
        module_name="scvi",
        install_cmd="pip install scvi-tools",
        description="Single-cell variational inference tools for cell type annotation and deconvolution",
        category=DependencyCategory.DEEP_LEARNING,
        min_version="1.0.0",
    ),
    "torch": DependencyInfo(
        module_name="torch",
        install_cmd="pip install torch",
        description="PyTorch deep learning framework",
        category=DependencyCategory.DEEP_LEARNING,
    ),
    "cell2location": DependencyInfo(
        module_name="cell2location",
        install_cmd="pip install cell2location",
        description="Probabilistic cell type deconvolution",
        category=DependencyCategory.DEEP_LEARNING,
    ),
    "flashdeconv": DependencyInfo(
        module_name="flashdeconv",
        install_cmd="pip install flashdeconv",
        description="Ultra-fast spatial transcriptomics deconvolution using random sketching",
        category=DependencyCategory.SPATIAL,
    ),
    # Spatial Analysis
    "tangram": DependencyInfo(
        module_name="tangram",
        install_cmd="pip install tangram-sc",
        description="Spatial mapping of single-cell transcriptomics",
        category=DependencyCategory.SPATIAL,
    ),
    "squidpy": DependencyInfo(
        module_name="squidpy",
        install_cmd="pip install squidpy",
        description="Spatial single-cell analysis",
        category=DependencyCategory.SPATIAL,
    ),
    "SpaGCN": DependencyInfo(
        module_name="SpaGCN",
        install_cmd="pip install SpaGCN",
        description="Spatial domain identification using graph convolutional networks",
        category=DependencyCategory.SPATIAL,
    ),
    "STAGATE": DependencyInfo(
        module_name="STAGATE_pyG",
        install_cmd="pip install STAGATE-pyG",
        description="Spatial domain identification using graph attention",
        category=DependencyCategory.SPATIAL,
    ),
    "GraphST": DependencyInfo(
        module_name="GraphST",
        install_cmd="pip install GraphST",
        description="Graph self-supervised contrastive learning for spatial domains",
        category=DependencyCategory.SPATIAL,
    ),
    "paste": DependencyInfo(
        module_name="paste",
        install_cmd="pip install paste-bio",
        description="Probabilistic alignment of spatial transcriptomics",
        category=DependencyCategory.SPATIAL,
    ),
    "stalign": DependencyInfo(
        module_name="STalign",
        install_cmd="pip install STalign",
        description="Spatial transcriptomics alignment",
        category=DependencyCategory.SPATIAL,
    ),
    # R Interface
    "rpy2": DependencyInfo(
        module_name="rpy2",
        install_cmd="pip install rpy2",
        description="R-Python interface (requires R installation)",
        category=DependencyCategory.R_INTERFACE,
    ),
    "anndata2ri": DependencyInfo(
        module_name="anndata2ri",
        install_cmd="pip install anndata2ri",
        description="AnnData to R SingleCellExperiment conversion",
        category=DependencyCategory.R_INTERFACE,
    ),
    # Cell Communication
    "liana": DependencyInfo(
        module_name="liana",
        install_cmd="pip install liana",
        description="Ligand-receptor analysis framework",
        category=DependencyCategory.COMMUNICATION,
    ),
    "cellphonedb": DependencyInfo(
        module_name="cellphonedb",
        install_cmd="pip install cellphonedb",
        description="Statistical method for cell-cell communication",
        category=DependencyCategory.COMMUNICATION,
    ),
    "ktplotspy": DependencyInfo(
        module_name="ktplotspy",
        install_cmd="pip install ktplotspy",
        description="CellPhoneDB visualization toolkit (dotplot, chord)",
        category=DependencyCategory.COMMUNICATION,
    ),
    # RNA Velocity
    "scvelo": DependencyInfo(
        module_name="scvelo",
        install_cmd="pip install scvelo",
        description="RNA velocity analysis",
        category=DependencyCategory.VELOCITY,
    ),
    "velovi": DependencyInfo(
        module_name="velovi",
        install_cmd="pip install velovi",
        description="Variational inference for RNA velocity (requires scvi-tools)",
        category=DependencyCategory.VELOCITY,
    ),
    "cellrank": DependencyInfo(
        module_name="cellrank",
        install_cmd="pip install cellrank",
        description="Trajectory inference using RNA velocity",
        category=DependencyCategory.VELOCITY,
    ),
    "palantir": DependencyInfo(
        module_name="palantir",
        install_cmd="pip install palantir",
        description="Trajectory inference for cell fate",
        category=DependencyCategory.VELOCITY,
    ),
    # Annotation
    "singler": DependencyInfo(
        module_name="singler",
        install_cmd="pip install singler singlecellexperiment",
        description="Reference-based cell type annotation",
        category=DependencyCategory.CORE,
    ),
    "mllmcelltype": DependencyInfo(
        module_name="mllmcelltype",
        install_cmd="pip install mllmcelltype",
        description="LLM-based cell type annotation",
        category=DependencyCategory.CORE,
    ),
    "celldex": DependencyInfo(
        module_name="celldex",
        install_cmd="pip install celldex",
        description="Cell type reference datasets for SingleR",
        category=DependencyCategory.CORE,
    ),
    # Enrichment
    "gseapy": DependencyInfo(
        module_name="gseapy",
        install_cmd="pip install gseapy",
        description="Gene set enrichment analysis",
        category=DependencyCategory.CORE,
    ),
    "decoupler": DependencyInfo(
        module_name="decoupler",
        install_cmd="pip install decoupler",
        description="Functional analysis of omics data",
        category=DependencyCategory.CORE,
    ),
    # Spatial Statistics
    "sparkx": DependencyInfo(
        module_name="sparkx",
        install_cmd="pip install SPARK-X",
        description="SPARK-X non-parametric spatial gene detection",
        category=DependencyCategory.SPATIAL,
    ),
    "spatialde": DependencyInfo(
        module_name="NaiveDE",
        install_cmd="pip install SpatialDE",
        description="SpatialDE Gaussian process-based spatial gene detection",
        category=DependencyCategory.SPATIAL,
    ),
    # CNV Analysis
    "infercnvpy": DependencyInfo(
        module_name="infercnvpy",
        install_cmd="pip install infercnvpy",
        description="Copy number variation inference",
        category=DependencyCategory.CORE,
    ),
    # Visualization
    "plotly": DependencyInfo(
        module_name="plotly",
        install_cmd="pip install plotly",
        description="Interactive visualization",
        category=DependencyCategory.VISUALIZATION,
    ),
    # Data handling
    "mudata": DependencyInfo(
        module_name="mudata",
        install_cmd="pip install mudata",
        description="Multimodal data handling",
        category=DependencyCategory.CORE,
    ),
    # Harmony integration
    "harmonypy": DependencyInfo(
        module_name="harmonypy",
        install_cmd="pip install harmonypy",
        description="Harmony batch integration",
        category=DependencyCategory.CORE,
    ),
    "scanorama": DependencyInfo(
        module_name="scanorama",
        install_cmd="pip install scanorama",
        description="Scanorama batch integration",
        category=DependencyCategory.CORE,
    ),
    "bbknn": DependencyInfo(
        module_name="bbknn",
        install_cmd="pip install bbknn",
        description="Batch balanced k-nearest neighbors",
        category=DependencyCategory.CORE,
    ),
    # Parallel computing
    "dask": DependencyInfo(
        module_name="dask",
        install_cmd="pip install dask",
        description="Parallel computing library",
        category=DependencyCategory.CORE,
    ),
    # Spatial weights for statistics
    "esda": DependencyInfo(
        module_name="esda",
        install_cmd="pip install esda",
        description="Exploratory spatial data analysis",
        category=DependencyCategory.SPATIAL,
    ),
    "libpysal": DependencyInfo(
        module_name="libpysal",
        install_cmd="pip install libpysal",
        description="Python spatial analysis library",
        category=DependencyCategory.SPATIAL,
    ),
    # Optimal transport
    "ot": DependencyInfo(
        module_name="ot",
        install_cmd="pip install POT",
        description="Python Optimal Transport library",
        category=DependencyCategory.SPATIAL,
    ),
    # Clustering
    "louvain": DependencyInfo(
        module_name="louvain",
        install_cmd="pip install louvain",
        description="Louvain community detection algorithm",
        category=DependencyCategory.CORE,
    ),
    # Enrichment analysis
    "enrichmap": DependencyInfo(
        module_name="enrichmap",
        install_cmd="pip install enrichmap",
        description="Spatial enrichment mapping",
        category=DependencyCategory.CORE,
    ),
    "pygam": DependencyInfo(
        module_name="pygam",
        install_cmd="pip install pygam",
        description="Generalized additive models",
        category=DependencyCategory.CORE,
    ),
    "skgstat": DependencyInfo(
        module_name="skgstat",
        install_cmd="pip install scikit-gstat",
        description="Geostatistical analysis toolkit",
        category=DependencyCategory.SPATIAL,
    ),
    "adjustText": DependencyInfo(
        module_name="adjustText",
        install_cmd="pip install adjustText",
        description="Text label placement for matplotlib",
        category=DependencyCategory.VISUALIZATION,
    ),
    "splot": DependencyInfo(
        module_name="splot",
        install_cmd="pip install splot",
        description="Spatial plotting for PySAL",
        category=DependencyCategory.VISUALIZATION,
    ),
    # Core scientific libraries
    "sklearn": DependencyInfo(
        module_name="sklearn",
        install_cmd="pip install scikit-learn",
        description="Machine learning library",
        category=DependencyCategory.CORE,
    ),
    "statsmodels": DependencyInfo(
        module_name="statsmodels",
        install_cmd="pip install statsmodels",
        description="Statistical models and tests",
        category=DependencyCategory.CORE,
    ),
    "scipy": DependencyInfo(
        module_name="scipy",
        install_cmd="pip install scipy",
        description="Scientific computing library",
        category=DependencyCategory.CORE,
    ),
    "scanpy": DependencyInfo(
        module_name="scanpy",
        install_cmd="pip install scanpy",
        description="Single-cell analysis in Python",
        category=DependencyCategory.CORE,
    ),
    "Pillow": DependencyInfo(
        module_name="PIL",
        install_cmd="pip install Pillow",
        description="Python Imaging Library for tissue image loading",
        category=DependencyCategory.CORE,
    ),
}


class DependencyManager:
    """Centralized manager for optional dependency handling.

    Thread-safe singleton that manages lazy loading and caching of optional
    dependencies with consistent error handling and logging.

    Example:
        deps = DependencyManager()

        # Require a dependency (raises ImportError if missing)
        scvi = deps.require("scvi-tools", ctx=ctx)

        # Get optional dependency (returns None if missing)
        torch = deps.get("torch")

        # Check availability
        if deps.is_available("rpy2"):
            rpy2 = deps.get("rpy2")
    """

    _instance: Optional["DependencyManager"] = None
    _lock = threading.Lock()
    _initialized: bool = False

    def __new__(cls) -> "DependencyManager":
        """Singleton pattern for global access."""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._instance._initialized = False
        return cls._instance

    def __init__(self) -> None:
        """Initialize the dependency manager."""
        if self._initialized:
            return
        self._cache: Dict[str, Any] = {}  # module_name -> module
        self._availability: Dict[str, bool] = {}  # module_name -> available
        self._import_lock = threading.Lock()
        self._initialized = True

    def _get_info(self, name: str) -> DependencyInfo:
        """Get dependency info from registry.

        Args:
            name: Dependency name (registry key or module name)

        Returns:
            DependencyInfo for the dependency

        Raises:
            KeyError: If dependency is not in registry
        """
        # Try direct lookup first
        if name in DEPENDENCY_REGISTRY:
            return DEPENDENCY_REGISTRY[name]

        # Try lookup by module name
        for _key, info in DEPENDENCY_REGISTRY.items():
            if info.module_name == name:
                return info

        # Not found - create default info
        return DependencyInfo(
            module_name=name,
            install_cmd=f"pip install {name}",
            description=f"Optional dependency: {name}",
        )

    def _try_import(self, module_name: str) -> Tuple[bool, Optional[Any]]:
        """Try to import a module.

        Args:
            module_name: The module name to import

        Returns:
            Tuple of (success, module_or_none)
        """
        with self._import_lock:
            # Check cache first
            if module_name in self._cache:
                return True, self._cache[module_name]
            if module_name in self._availability:
                if not self._availability[module_name]:
                    return False, None

            # Try import
            try:
                module = importlib.import_module(module_name)
                self._cache[module_name] = module
                self._availability[module_name] = True
                return True, module
            except ImportError:
                self._availability[module_name] = False
                return False, None

    def is_available(self, name: str) -> bool:
        """Check if a dependency is available without importing.

        Uses importlib.util.find_spec for fast availability check.

        Args:
            name: Dependency name (registry key or module name)

        Returns:
            True if the dependency is available
        """
        info = self._get_info(name)
        module_name = info.module_name

        # Check cache first
        if module_name in self._availability:
            return self._availability[module_name]

        # Use find_spec for fast check without full import
        spec = importlib.util.find_spec(module_name)
        available = spec is not None
        self._availability[module_name] = available
        return available

    def get(
        self,
        name: str,
        ctx: Optional["ToolContext"] = None,
        warn_if_missing: bool = False,
    ) -> Optional[Any]:
        """Get an optional dependency, returning None if unavailable.

        This is the preferred method for truly optional dependencies where
        the code has a fallback path.

        Args:
            name: Dependency name (registry key or module name)
            ctx: ToolContext for logging (uses debug level for technical messages)
            warn_if_missing: Whether to warn if dependency is missing

        Returns:
            The imported module or None if unavailable
        """
        info = self._get_info(name)
        success, module = self._try_import(info.module_name)

        if success:
            if ctx:
                version = getattr(module, "__version__", "unknown")
                ctx.debug(f"Loaded {name} (version {version})")
            return module

        if warn_if_missing:
            msg = (
                f"{name} not available. Install with: {info.install_cmd}\n"
                f"Description: {info.description}"
            )
            warnings.warn(msg, stacklevel=2)
            if ctx:
                ctx.debug(msg)

        return None

    def require(
        self,
        name: str,
        ctx: Optional["ToolContext"] = None,
        feature: Optional[str] = None,
    ) -> Any:
        """Require a dependency, raising ImportError if unavailable.

        This is the preferred method for dependencies that are required
        for a specific feature to work.

        Args:
            name: Dependency name (registry key or module name)
            ctx: ToolContext for logging (uses debug level for technical messages)/error reporting
            feature: Optional feature name for better error messages

        Returns:
            The imported module

        Raises:
            ImportError: If the dependency is not available
        """
        info = self._get_info(name)
        success, module = self._try_import(info.module_name)

        if success:
            if ctx:
                version = getattr(module, "__version__", "unknown")
                ctx.debug(f"Using {name} (version {version})")
            return module

        # Build helpful error message
        feature_msg = f" for {feature}" if feature else ""
        error_msg = (
            f"{name} is required{feature_msg}.\n\n"
            f"Install with: {info.install_cmd}\n\n"
            f"Description: {info.description}"
        )

        if ctx:
            ctx.debug(f"Missing dependency: {error_msg}")

        raise ImportError(error_msg)


# Convenience functions for common operations
# Note: DependencyManager uses __new__ singleton, so DependencyManager() always
# returns the same instance. No need for module-level caching.


def require(
    name: str,
    ctx: Optional["ToolContext"] = None,
    feature: Optional[str] = None,
) -> Any:
    """Require a dependency (convenience function).

    See DependencyManager.require for details.
    """
    return DependencyManager().require(name, ctx, feature)


def get(
    name: str,
    ctx: Optional["ToolContext"] = None,
    warn_if_missing: bool = False,
) -> Optional[Any]:
    """Get an optional dependency (convenience function).

    See DependencyManager.get for details.
    """
    return DependencyManager().get(name, ctx, warn_if_missing)


def is_available(name: str) -> bool:
    """Check if a dependency is available (convenience function).

    See DependencyManager.is_available for details.
    """
    return DependencyManager().is_available(name)


# R-specific validation (common pattern across tools)
def validate_r_environment(
    ctx: Optional["ToolContext"] = None,
    required_packages: Optional[List[str]] = None,
) -> Tuple[Any, ...]:
    """Validate R environment and return required modules.

    This replaces the scattered _validate_rpy2_and_r patterns.

    Args:
        ctx: ToolContext for logging (uses debug level for technical messages)
        required_packages: Optional list of R packages to check

    Returns:
        Tuple of (robjects, pandas2ri, numpy2ri, importr, localconverter,
                  default_converter, openrlib, anndata2ri)

    Raises:
        ImportError: If rpy2 or required R packages are not available
    """
    manager = DependencyManager()

    # Check Python packages
    if not manager.is_available("rpy2"):
        raise ImportError(
            "rpy2 is required for R-based methods. "
            "Install with: pip install rpy2 (requires R installation)"
        )

    if not manager.is_available("anndata2ri"):
        raise ImportError(
            "anndata2ri is required for R-based methods. "
            "Install with: pip install anndata2ri"
        )

    # Import required modules
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

                if ctx:
                    r_version = robjects.r("R.version.string")[0]
                    ctx.debug(f"Using R: {r_version}")

        # Check required R packages if specified
        if required_packages:
            missing_r = []
            for pkg in required_packages:
                try:
                    with openrlib.rlock:
                        with conversion.localconverter(default_converter):
                            importr(pkg)
                except Exception:
                    missing_r.append(pkg)

            if missing_r:
                pkg_list = ", ".join(f"'{p}'" for p in missing_r)
                raise ImportError(
                    f"Missing R packages: {pkg_list}\n"
                    f"Install in R with: install.packages(c({pkg_list}))"
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
            f"R environment setup failed: {str(e)}\n\n"
            "Common solutions:\n"
            "  - Install R from https://www.r-project.org/\n"
            "  - Set R_HOME environment variable\n"
            "  - macOS: brew install r\n"
            "  - Ubuntu: sudo apt install r-base"
        ) from e


# scvi-tools specific validation (common pattern)
def validate_scvi_tools(
    ctx: Optional["ToolContext"] = None,
    components: Optional[List[str]] = None,
) -> Any:
    """Validate scvi-tools availability and return the module.

    This replaces the scattered _validate_scvi_tools patterns.

    Args:
        ctx: ToolContext for logging (uses debug level for technical messages)
        components: Optional list of specific components to validate
                   (e.g., ["CellAssign", "Cell2location"])

    Returns:
        The scvi module

    Raises:
        ImportError: If scvi-tools or required components are not available
    """
    manager = DependencyManager()
    scvi = manager.require("scvi-tools", ctx, "scvi-tools methods")

    if components:
        missing = []
        for comp in components:
            try:
                # Try to import the component
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
                    # Generic import attempt
                    getattr(scvi, comp, None) or getattr(
                        scvi.model, comp, None
                    ) or getattr(scvi.external, comp, None)
            except (ImportError, AttributeError):
                missing.append(comp)

        if missing:
            raise ImportError(
                f"scvi-tools components not available: {', '.join(missing)}\n"
                "This may be due to version incompatibility. "
                "Try: pip install --upgrade scvi-tools"
            )

    return scvi


def validate_r_package(
    package_name: str,
    ctx: Optional["ToolContext"] = None,
    install_cmd: Optional[str] = None,
) -> bool:
    """Check if an R package is available and can be loaded.

    This provides a simpler interface than validate_r_environment for checking
    individual R packages without returning the full rpy2 module tuple.

    Args:
        package_name: R package name (e.g., "sctransform", "numbat", "CellChat")
        ctx: ToolContext for logging (uses debug level for technical messages)
        install_cmd: Custom install command. If None, uses install.packages()

    Returns:
        True if the package is available

    Raises:
        ImportError: If rpy2 not available or R package cannot be loaded

    Example:
        # Check for sctransform
        validate_r_package("sctransform", ctx)

        # Check with custom install command
        validate_r_package(
            "numbat",
            ctx,
            install_cmd="devtools::install_github('kharchenkolab/numbat')"
        )
    """
    manager = DependencyManager()

    if not manager.is_available("rpy2"):
        raise ImportError(
            "rpy2 is required for R-based methods.\n"
            "Install with: pip install rpy2\n"
            "Note: R must also be installed on your system."
        )

    try:
        from rpy2.rinterface_lib import openrlib
        from rpy2.robjects import conversion, default_converter
        from rpy2.robjects.packages import importr

        with openrlib.rlock:
            with conversion.localconverter(default_converter):
                importr(package_name)

        if ctx:
            ctx.debug(f"R package '{package_name}' is available")
        return True

    except Exception as e:
        error_str = str(e).lower()
        if "there is no package" in error_str or package_name.lower() in error_str:
            default_install = f"install.packages('{package_name}')"
            install = install_cmd or default_install
            raise ImportError(
                f"R package '{package_name}' is not installed.\n"
                f"Install in R: {install}"
            ) from e
        else:
            raise ImportError(
                f"Failed to load R package '{package_name}': {e}\n\n"
                "Common solutions:\n"
                "  - Ensure R is properly installed\n"
                "  - Set R_HOME environment variable\n"
                "  - Check that rpy2 can find your R installation"
            ) from e


def check_r_packages(
    packages: List[str],
    ctx: Optional["ToolContext"] = None,
) -> List[str]:
    """Check availability of multiple R packages.

    Args:
        packages: List of R package names to check
        ctx: ToolContext for logging (uses debug level for technical messages)

    Returns:
        List of missing package names (empty if all available)
    """
    manager = DependencyManager()

    if not manager.is_available("rpy2"):
        return packages  # All missing if rpy2 not available

    missing = []
    for pkg in packages:
        try:
            validate_r_package(pkg)
        except ImportError:
            missing.append(pkg)

    if missing and ctx:
        ctx.debug(f"Missing R packages: {', '.join(missing)}")

    return missing
