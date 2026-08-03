"""
Unified dependency management for ChatSpatial MCP.

Provides a consistent API for managing optional dependencies, replacing
scattered try/except ImportError patterns with centralized handling.

Usage:
    # Require a dependency (raises if missing)
    scvi = require("scvi-tools", feature="cell type annotation")

    # Require a specific module provided by a dependency
    external = require_module("scvi", "scvi.external", feature="VELOVI")

    # Get optional dependency (returns None if missing)
    torch = get("torch")

    # Check availability
    if is_available("rpy2"):
        import rpy2
"""

import importlib
import importlib.util
import warnings
from contextlib import contextmanager
from dataclasses import dataclass, field
from functools import lru_cache
from typing import TYPE_CHECKING, Any, Iterator, Optional

from .exceptions import DependencyError

if TYPE_CHECKING:
    from ..spatial_mcp_adapter import ToolContext


@dataclass(frozen=True)
class DependencyInfo:
    """Metadata for an optional dependency."""

    module_name: str
    install_cmd: str
    description: str = ""


@dataclass(frozen=True)
class REnvironment:
    """Named handles for a validated rpy2 runtime."""

    robjects: Any
    pandas2ri: Any
    numpy2ri: Any
    packages: Any
    conversion: Any
    openrlib: Any
    anndata2ri: Optional[Any]
    loaded_packages: dict[str, Any] = field(default_factory=dict)

    @property
    def importr(self) -> Any:
        """Return rpy2's R-package importer."""
        return self.packages.importr

    def package(self, name: str) -> Any:
        """Return an R package loaded during environment validation."""
        try:
            return self.loaded_packages[name]
        except KeyError as exc:
            raise ValueError(
                f"R package '{name}' was not requested for validation."
            ) from exc

    def _build_converter(
        self,
        *,
        pandas: bool,
        numpy: bool,
        anndata: bool,
    ) -> Any:
        """Compose the converter required by one R execution scope."""
        converter = self.robjects.default_converter
        if anndata:
            if self.anndata2ri is None:
                raise DependencyError(
                    "anndata2ri was not loaded for this R environment."
                )
            converter += self.anndata2ri.converter
        if pandas:
            converter += self.pandas2ri.converter
        if numpy:
            converter += self.numpy2ri.converter
        return converter

    @contextmanager
    def conversion_context(
        self,
        *,
        pandas: bool = False,
        numpy: bool = False,
        anndata: bool = False,
    ) -> Iterator[None]:
        """Hold the R lock while applying the requested Python converters."""
        converter = self._build_converter(
            pandas=pandas,
            numpy=numpy,
            anndata=anndata,
        )

        with self.openrlib.rlock, self.conversion.localconverter(converter):
            yield

    @contextmanager
    def local_context(
        self,
        *,
        pandas: bool = False,
        numpy: bool = False,
        anndata: bool = False,
    ) -> Iterator[Any]:
        """Evaluate R code in an isolated child environment.

        The child inherits package visibility from the current R environment,
        but request inputs and intermediate objects never enter ``.GlobalEnv``.
        The outer lock serializes the embedded R runtime; ``use_rlock=False``
        prevents rpy2 from acquiring that same lock a second time.
        """
        converter = self._build_converter(
            pandas=pandas,
            numpy=numpy,
            anndata=anndata,
        )
        with (
            self.openrlib.rlock,
            self.conversion.localconverter(converter),
            self.robjects.local_context(use_rlock=False) as environment,
        ):
            yield environment


# Registry of optional dependencies with install instructions
DEPENDENCY_REGISTRY: dict[str, DependencyInfo] = {
    # Deep Learning
    "scvi-tools": DependencyInfo(
        "scvi",
        "pip install 'chatspatial[deep-learning]'",
        "Single-cell variational inference tools",
    ),
    "torch": DependencyInfo(
        "torch",
        "pip install 'chatspatial[deep-learning]'",
        "PyTorch deep learning framework",
    ),
    "cell2location": DependencyInfo(
        "cell2location",
        "pip install 'chatspatial[deconvolution]'",
        "Probabilistic cell type deconvolution",
    ),
    "flashdeconv": DependencyInfo(
        "flashdeconv",
        "pip install 'chatspatial[deconvolution]'",
        "Ultra-fast spatial deconvolution",
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
        "pip install 'chatspatial[spatial-domains]'",
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
    "banksy": DependencyInfo(
        "banksy",
        "pip install 'chatspatial[spatial-domains]'",
        "Spatial domain identification using neighborhood feature augmentation",
    ),
    "aestetik": DependencyInfo(
        "aestetik",
        "pip install 'chatspatial[aestetik]'",
        "Multi-modal autoencoder over expression, morphology, and spatial context",
    ),
    "paste": DependencyInfo(
        "paste",
        "pip install paste-bio",
        "Probabilistic alignment of spatial transcriptomics",
    ),
    "stalign": DependencyInfo(
        "STalign",
        'pip install --upgrade "git+https://github.com/JEFworks-Lab/STalign.git"',
        "Spatial transcriptomics alignment",
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
        "liana",
        "pip install 'chatspatial[cell-communication]'",
        "Ligand-receptor analysis framework",
    ),
    "cellphonedb": DependencyInfo(
        "cellphonedb",
        "pip install 'chatspatial[cell-communication]'",
        "Statistical method for cell-cell communication",
    ),
    "fastccc": DependencyInfo(
        "fastccc",
        "pip install 'chatspatial[cell-communication]'",
        "Fast cell-cell communication analysis",
    ),
    # RNA Velocity
    "scvelo": DependencyInfo(
        "scvelo", "pip install 'chatspatial[velocity]'", "RNA velocity analysis"
    ),
    "velovi": DependencyInfo(
        "scvi",
        "pip install 'chatspatial[deep-learning]'",
        "VELOVI model provided by scvi-tools",
    ),
    "cellrank": DependencyInfo(
        "cellrank",
        "pip install 'chatspatial[trajectory]'",
        "Trajectory inference using RNA velocity",
    ),
    "palantir": DependencyInfo(
        "palantir",
        "pip install 'chatspatial[trajectory]'",
        "Trajectory inference for cell fate",
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
        "SpatialDE",
        "pip install SpatialDE",
        "SpatialDE Gaussian process spatial gene detection",
    ),
    "naivede": DependencyInfo(
        "NaiveDE",
        "pip install SpatialDE",
        "Variance stabilization used by SpatialDE",
    ),
    "flashs": DependencyInfo(
        "flashs",
        "pip install flashs",
        "FlashS ultra-fast Python-native spatial gene detection",
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
        "mudata", "pip install 'chatspatial[deep-learning]'", "Multimodal data handling"
    ),
    # Integration
    "harmonypy": DependencyInfo(
        "harmonypy",
        "pip install 'chatspatial[integration]'",
        "Harmony batch integration",
    ),
    "scanorama": DependencyInfo(
        "scanorama",
        "pip install 'chatspatial[integration]'",
        "Scanorama batch integration",
    ),
    "bbknn": DependencyInfo(
        "bbknn",
        "pip install 'chatspatial[integration]'",
        "Batch balanced k-nearest neighbors",
    ),
    # Spatial weights
    "esda": DependencyInfo(
        "esda",
        "pip install 'chatspatial[spatial-stats]'",
        "Exploratory spatial data analysis",
    ),
    "libpysal": DependencyInfo(
        "libpysal",
        "pip install 'chatspatial[spatial-stats]'",
        "Python spatial analysis library",
    ),
    # Other
    "dask": DependencyInfo("dask", "pip install dask", "Parallel computing library"),
    "ot": DependencyInfo("ot", "pip install POT", "Python Optimal Transport library"),
    "louvain": DependencyInfo(
        "louvain",
        "pip install 'chatspatial[spatial-domains]'",
        "Louvain community detection algorithm",
    ),
    "pydeseq2": DependencyInfo(
        "pydeseq2", "pip install pydeseq2", "Python implementation of DESeq2"
    ),
    "enrichmap": DependencyInfo(
        "enrichmap", "pip install enrichmap", "Spatial enrichment mapping"
    ),
    "pygam": DependencyInfo(
        "pygam",
        "pip install 'chatspatial[trajectory]'",
        "Generalized additive models",
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
    """Import a module, returning ``None`` only when that module is absent.

    Import failures raised from inside an installed package must propagate so
    callers can distinguish a missing package from a broken dependency stack.
    """
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        missing_name = exc.name
        if missing_name and (
            module_name == missing_name or module_name.startswith(f"{missing_name}.")
        ):
            return None
        raise


@lru_cache(maxsize=256)
def _check_spec(module_name: str) -> bool:
    """Fast availability check without importing."""
    return importlib.util.find_spec(module_name) is not None


# =============================================================================
# Public API
# =============================================================================


def is_available(name: str) -> bool:
    """Check whether a dependency is installed without importing it.

    This is a fast presence check, not a guarantee that the installed package
    and all of its transitive dependencies are importable.
    """
    return _check_spec(_get_info(name).module_name)


def _load_dependency(
    name: str,
    *,
    feature: Optional[str] = None,
    module_name: Optional[str] = None,
) -> tuple[DependencyInfo, Optional[Any]]:
    """Load a dependency while preserving broken-import diagnostics."""
    info = _get_info(name)
    import_name = module_name or info.module_name
    try:
        module = _try_import(import_name)
    except Exception as exc:
        feature_msg = f" for {feature}" if feature else ""
        if import_name == info.module_name:
            subject = f"{name} is installed but could not be imported"
        else:
            subject = (
                f"{name} is installed, but required module '{import_name}' "
                "could not be imported"
            )
        raise DependencyError(
            f"{subject}{feature_msg}.\n\n"
            f"Import failure: {exc}\n"
            f"Install or repair: {info.install_cmd}\n"
            f"Description: {info.description}"
        ) from exc
    return info, module


def get(
    name: str,
    ctx: Optional["ToolContext"] = None,
    warn_if_missing: bool = False,
) -> Optional[Any]:
    """Get an optional dependency, returning ``None`` only if it is absent.

    Raises:
        DependencyError: If the package is installed but cannot be imported.
    """
    info, module = _load_dependency(name)

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
    """Return a required dependency or raise a domain dependency error.

    Raises:
        DependencyError: If the package is absent or cannot be imported.
    """
    info, module = _load_dependency(name, feature=feature)

    if module is not None:
        return module

    feature_msg = f" for {feature}" if feature else ""
    raise DependencyError(
        f"{name} is required{feature_msg}.\n\n"
        f"Install: {info.install_cmd}\n"
        f"Description: {info.description}"
    )


def require_module(
    name: str,
    module_name: str,
    ctx: Optional["ToolContext"] = None,
    feature: Optional[str] = None,
) -> Any:
    """Return a required module provided by an optional dependency.

    The dependency root is loaded first so an absent package and an installed
    package with a missing or broken submodule produce distinct diagnostics.

    Raises:
        ValueError: If ``module_name`` is outside the dependency package.
        DependencyError: If the package or requested module is unavailable.
    """
    info = _get_info(name)
    package_name = info.module_name
    if module_name != package_name and not module_name.startswith(f"{package_name}."):
        raise ValueError(
            f"Module '{module_name}' is not provided by dependency '{name}' "
            f"(expected '{package_name}' or one of its submodules)."
        )

    package = require(name, ctx, feature=feature)
    if module_name == package_name:
        return package

    _, module = _load_dependency(
        name,
        feature=feature,
        module_name=module_name,
    )
    if module is not None:
        return module

    feature_msg = f" for {feature}" if feature else ""
    raise DependencyError(
        f"{name} is installed, but required module '{module_name}' is unavailable"
        f"{feature_msg}.\n\n"
        f"Install or repair: {info.install_cmd}\n"
        f"Description: {info.description}"
    )


# =============================================================================
# R Environment Validation
# =============================================================================


def validate_r_environment(
    ctx: Optional["ToolContext"] = None,
    required_packages: Optional[list[str]] = None,
    *,
    require_anndata2ri: bool = False,
    package_install_commands: Optional[dict[str, str]] = None,
) -> REnvironment:
    """Validate R environment and return required modules.

    Returns:
        Named R environment handles and an optional AnnData converter.

    Raises:
        DependencyError: If Python/R dependencies or required R packages fail
    """
    try:
        feature = "R-based methods"
        robjects = require_module("rpy2", "rpy2.robjects", ctx, feature=feature)
        conversion = require_module(
            "rpy2", "rpy2.robjects.conversion", ctx, feature=feature
        )
        packages = require_module(
            "rpy2", "rpy2.robjects.packages", ctx, feature=feature
        )
        rinterface_lib = require_module(
            "rpy2", "rpy2.rinterface_lib", ctx, feature=feature
        )
        pandas2ri = require_module(
            "rpy2", "rpy2.robjects.pandas2ri", ctx, feature=feature
        )
        numpy2ri = require_module(
            "rpy2", "rpy2.robjects.numpy2ri", ctx, feature=feature
        )
        anndata2ri = (
            require("anndata2ri", ctx, feature=feature) if require_anndata2ri else None
        )
        environment = REnvironment(
            robjects=robjects,
            pandas2ri=pandas2ri,
            numpy2ri=numpy2ri,
            packages=packages,
            conversion=conversion,
            openrlib=rinterface_lib.openrlib,
            anndata2ri=anndata2ri,
        )

        # Test R availability
        with environment.conversion_context():
            robjects.r("R.version")

        # Check required R packages
        if required_packages:
            missing = []
            for pkg in required_packages:
                try:
                    with environment.conversion_context():
                        environment.loaded_packages[pkg] = environment.importr(pkg)
                except packages.PackageNotInstalledError:
                    missing.append(pkg)
                except packages.LibraryError as e:
                    raise DependencyError(
                        f"R package '{pkg}' is installed but could not be loaded: {e}"
                    ) from e
                except Exception as e:
                    raise DependencyError(
                        f"Failed to validate R package '{pkg}': {e}"
                    ) from e

            if missing:
                install_commands = package_install_commands or {}
                if len(missing) == 1 and missing[0] in install_commands:
                    pkg = missing[0]
                    raise DependencyError(
                        f"R package '{pkg}' not installed.\n"
                        f"Install in R: {install_commands[pkg]}"
                    )
                pkg_list = ", ".join(f"'{p}'" for p in missing)
                if install_commands:
                    install_lines = []
                    for pkg in missing:
                        command = install_commands.get(
                            pkg, f"install.packages('{pkg}')"
                        )
                        install_lines.append(f"  - {pkg}: {command}")
                    instructions = "\n".join(install_lines)
                else:
                    instructions = f"install.packages(c({pkg_list}))"
                raise DependencyError(
                    f"Missing R packages: {pkg_list}\nInstall in R:\n{instructions}"
                )

        return environment

    except DependencyError:
        raise
    except Exception as e:
        raise DependencyError(
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
    """Check whether an R package can be imported."""
    environment = validate_r_environment(ctx)

    try:
        with environment.conversion_context():
            environment.importr(package_name)

        return True

    except environment.packages.PackageNotInstalledError as e:
        install = install_cmd or f"install.packages('{package_name}')"
        raise DependencyError(
            f"R package '{package_name}' not installed.\nInstall in R: {install}"
        ) from e
    except environment.packages.LibraryError as e:
        raise DependencyError(
            f"R package '{package_name}' is installed but could not be loaded: {e}\n"
            "Check the package version and its system-library dependencies."
        ) from e
    except Exception as e:
        raise DependencyError(
            f"R environment failed while validating package '{package_name}': {e}"
        ) from e


def check_r_packages(
    packages: list[str],
    ctx: Optional["ToolContext"] = None,
) -> list[str]:
    """Check availability of multiple R packages. Returns missing ones."""
    try:
        environment = validate_r_environment(ctx)
    except DependencyError:
        return packages

    missing = []
    for pkg in packages:
        try:
            with environment.conversion_context():
                environment.importr(pkg)
        except (
            environment.packages.PackageNotInstalledError,
            environment.packages.LibraryError,
        ):
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
                    module = require_module(
                        "scvi-tools",
                        "scvi.external",
                        ctx,
                        feature="scvi-tools methods",
                    )
                    _ = module.CellAssign
                elif comp == "Cell2location":
                    require("cell2location", ctx, feature="scvi-tools methods")
                elif comp == "SCANVI":
                    module = require_module(
                        "scvi-tools",
                        "scvi.model",
                        ctx,
                        feature="scvi-tools methods",
                    )
                    _ = module.SCANVI
                elif comp in {"DestVI", "Stereoscope"}:
                    module = require_module(
                        "scvi-tools",
                        "scvi.external",
                        ctx,
                        feature="scvi-tools methods",
                    )
                    getattr(module, comp)
                else:
                    getattr(scvi, comp, None) or getattr(
                        scvi.model, comp, None
                    ) or getattr(scvi.external, comp, None)
            except AttributeError:
                missing.append(comp)

        if missing:
            raise DependencyError(
                f"scvi-tools components not available: {', '.join(missing)}\n"
                "Try: pip install --upgrade scvi-tools"
            )

    return scvi
