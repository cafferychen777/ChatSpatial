"""
Smart dependency management system for ChatSpatial.

This module provides intelligent dependency detection, graceful degradation,
and user-friendly error messages. It eliminates the chaos of scattered
try/except ImportError blocks across the codebase.

Design Philosophy (Linus-style):
- No special cases: Unified handling for all dependencies
- Good taste: Simple data structures, no complex inheritance
- Never break userspace: Existing functionality preserved
- Practical: Solves real user problems, not theoretical ones

Author: Linus-approved dependency management system
"""

import importlib
import importlib.util
import logging
import subprocess
import sys
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache, wraps
from typing import Any, Callable, Dict, List, Optional, Union

# Setup logging
logger = logging.getLogger(__name__)


class DependencyLevel(Enum):
    """Dependency importance levels"""

    CORE = "core"  # Essential - ChatSpatial won't work without these
    RECOMMENDED = "recommended"  # Most users need these
    ADVANCED = "advanced"  # Power users and specific methods
    EXPERIMENTAL = "experimental"  # Research and bleeding-edge features


@dataclass
class DependencySpec:
    """Specification for a dependency"""

    name: str  # Package name (for import)
    pip_name: str  # Package name for pip install
    level: DependencyLevel  # Importance level
    description: str  # What this enables
    install_cmd: str = ""  # Custom install command
    alternatives: List[str] = field(default_factory=list)  # Alternative packages
    min_version: Optional[str] = None  # Minimum required version
    max_version: Optional[str] = None  # Maximum compatible version
    platform_specific: Dict[str, str] = field(
        default_factory=dict
    )  # OS-specific install
    conda_name: Optional[str] = None  # Conda package name if different

    def __post_init__(self):
        """Generate install command if not provided"""
        if not self.install_cmd:
            self.install_cmd = f"pip install {self.pip_name}"


@dataclass
class ImportResult:
    """Result of an import attempt"""

    success: bool
    module: Optional[Any] = None
    error: Optional[str] = None
    suggestion: Optional[str] = None
    alternatives: List[str] = field(default_factory=list)


class SmartImporter:
    """
    Central dependency management system.

    This is the "good taste" solution - one place to handle all dependencies,
    no special cases scattered everywhere.
    """

    def __init__(self):
        self._cache: Dict[str, ImportResult] = {}
        self._dependency_registry = self._build_dependency_registry()
        self._availability_matrix: Dict[str, bool] = {}
        self._feature_matrix: Dict[str, Dict[str, bool]] = {}

    def _build_dependency_registry(self) -> Dict[str, DependencySpec]:
        """
        Build the complete dependency registry.
        This is the single source of truth for all dependencies.
        """
        registry = {}

        # CORE dependencies - ChatSpatial breaks without these
        core_deps = [
            DependencySpec(
                "numpy", "numpy", DependencyLevel.CORE, "Core numerical computing"
            ),
            DependencySpec(
                "pandas",
                "pandas",
                DependencyLevel.CORE,
                "Data manipulation and analysis",
            ),
            DependencySpec(
                "scanpy", "scanpy", DependencyLevel.CORE, "Single-cell analysis toolkit"
            ),
            DependencySpec(
                "anndata", "anndata", DependencyLevel.CORE, "Annotated data matrices"
            ),
            DependencySpec(
                "matplotlib",
                "matplotlib",
                DependencyLevel.CORE,
                "Basic plotting functionality",
            ),
            DependencySpec(
                "scipy", "scipy", DependencyLevel.CORE, "Scientific computing"
            ),
            DependencySpec(
                "sklearn",
                "scikit-learn",
                DependencyLevel.CORE,
                "Machine learning toolkit",
                pip_name="scikit-learn",
            ),
        ]

        # RECOMMENDED dependencies - Most users need these
        recommended_deps = [
            DependencySpec(
                "squidpy",
                "squidpy",
                DependencyLevel.RECOMMENDED,
                "Spatial transcriptomics analysis toolkit",
            ),
            DependencySpec(
                "umap",
                "umap-learn",
                DependencyLevel.RECOMMENDED,
                "Dimensionality reduction and visualization",
            ),
            DependencySpec(
                "igraph",
                "python-igraph",
                DependencyLevel.RECOMMENDED,
                "Graph analysis and network algorithms",
            ),
            DependencySpec(
                "leidenalg",
                "leidenalg",
                DependencyLevel.RECOMMENDED,
                "Community detection clustering",
            ),
            DependencySpec(
                "harmonypy",
                "harmonypy",
                DependencyLevel.RECOMMENDED,
                "Batch effect correction",
            ),
            DependencySpec(
                "scvelo", "scvelo", DependencyLevel.RECOMMENDED, "RNA velocity analysis"
            ),
        ]

        # ADVANCED dependencies - Specific methods
        advanced_deps = [
            # Deconvolution methods
            DependencySpec(
                "cell2location",
                "cell2location",
                DependencyLevel.ADVANCED,
                "Spatial deconvolution with cell2location",
                install_cmd="pip install cell2location",
            ),
            # Spatial domain identification
            DependencySpec(
                "SpaGCN",
                "SpaGCN",
                DependencyLevel.ADVANCED,
                "Spatial domain identification with SpaGCN",
                install_cmd="pip install SpaGCN",
            ),
            DependencySpec(
                "STAGATE",
                "STAGATE",
                DependencyLevel.ADVANCED,
                "Spatial domain identification with STAGATE",
                install_cmd="git clone https://github.com/QIFEIDKN/STAGATE.git && cd STAGATE && python setup.py install",
            ),
            DependencySpec(
                "banksy",
                "banksy",
                DependencyLevel.ADVANCED,
                "Neighborhood aggregation for spatial domains",
                install_cmd="git clone https://github.com/prabhakarlab/Banksy_py.git && cd Banksy_py && pip install -r requirements.txt",
            ),
            # Spatial gene detection
            DependencySpec(
                "spatialde",
                "spatialde",
                DependencyLevel.ADVANCED,
                "Spatial variable gene detection",
                install_cmd="pip install spatialde",
            ),
            DependencySpec(
                "NaiveDE",
                "NaiveDE",
                DependencyLevel.ADVANCED,
                "Differential expression for SpatialDE",
                install_cmd="pip install NaiveDE",
            ),
            # Cell communication
            DependencySpec(
                "liana",
                "liana",
                DependencyLevel.ADVANCED,
                "Cell-cell communication analysis with LIANA+",
                install_cmd="pip install liana",
            ),
            DependencySpec(
                "cellphonedb",
                "cellphonedb",
                DependencyLevel.ADVANCED,
                "Cell-cell communication analysis with CellPhoneDB",
                install_cmd="pip install cellphonedb",
            ),
            # Spatial registration
            DependencySpec(
                "paste",
                "paste-bio",
                DependencyLevel.ADVANCED,
                "Spatial alignment with PASTE",
                install_cmd="pip install paste-bio",
            ),
            # Deep learning frameworks
            DependencySpec(
                "torch",
                "torch",
                DependencyLevel.ADVANCED,
                "PyTorch for deep learning methods",
                install_cmd="pip install torch",
            ),
            DependencySpec(
                "tensorflow",
                "tensorflow",
                DependencyLevel.ADVANCED,
                "TensorFlow for deep learning methods",
                install_cmd="pip install 'tensorflow>=1.15.0,<2.0'",
            ),
            # R interface
            DependencySpec(
                "rpy2",
                "rpy2",
                DependencyLevel.ADVANCED,
                "R interface for sc-type and other R methods",
                install_cmd="pip install 'rpy2>=3.4.0'",
            ),
        ]

        # EXPERIMENTAL dependencies - Cutting edge research
        experimental_deps = [
            DependencySpec(
                "enrichmap",
                "enrichmap",
                DependencyLevel.EXPERIMENTAL,
                "Enhanced pathway enrichment visualization",
                install_cmd="pip install enrichmap",
            ),
            DependencySpec(
                "petsc4py",
                "petsc4py",
                DependencyLevel.EXPERIMENTAL,
                "High-performance linear algebra for CellRank",
                install_cmd="pip install petsc4py",
            ),
            DependencySpec(
                "slepc4py",
                "slepc4py",
                DependencyLevel.EXPERIMENTAL,
                "Eigenvalue solvers for CellRank",
                install_cmd="pip install slepc4py",
            ),
        ]

        # Build registry
        for dep_list in [core_deps, recommended_deps, advanced_deps, experimental_deps]:
            for dep in dep_list:
                registry[dep.name] = dep

        return registry

    @lru_cache(maxsize=128)
    def smart_import(
        self,
        package_name: str,
        fallback_name: Optional[str] = None,
        required: bool = True,
    ) -> ImportResult:
        """
        Intelligently import a package with comprehensive error handling.

        Args:
            package_name: Name of the package to import
            fallback_name: Alternative package name to try
            required: Whether this import is required for functionality

        Returns:
            ImportResult with detailed information
        """
        if package_name in self._cache:
            return self._cache[package_name]

        result = self._attempt_import(package_name)

        # Try fallback if main import failed
        if not result.success and fallback_name:
            result = self._attempt_import(fallback_name)

        # Add intelligent suggestions
        if not result.success:
            result = self._enhance_error_with_suggestions(
                package_name, result, required
            )

        # Cache the result
        self._cache[package_name] = result
        return result

    def _attempt_import(self, package_name: str) -> ImportResult:
        """Attempt to import a package"""
        try:
            module = importlib.import_module(package_name)
            return ImportResult(success=True, module=module)
        except ImportError as e:
            return ImportResult(success=False, error=str(e))
        except Exception as e:
            # Handle other import-related errors (version conflicts, etc.)
            return ImportResult(success=False, error=f"Import failed: {str(e)}")

    def _enhance_error_with_suggestions(
        self, package_name: str, result: ImportResult, required: bool
    ) -> ImportResult:
        """Add intelligent suggestions to failed imports"""
        dep_spec = self._dependency_registry.get(package_name)

        if not dep_spec:
            # Package not in our registry - generic advice
            result.suggestion = (
                f"Package '{package_name}' not found. Try: pip install {package_name}"
            )
            return result

        # Build comprehensive suggestion
        suggestion_parts = [
            f"Missing {dep_spec.level.value} dependency: {dep_spec.description}",
            f"Install with: {dep_spec.install_cmd}",
        ]

        # Add level-specific advice
        if dep_spec.level == DependencyLevel.CORE:
            suggestion_parts.append(
                "⚠️  This is required for ChatSpatial to function properly."
            )
        elif dep_spec.level == DependencyLevel.RECOMMENDED:
            suggestion_parts.append("💡 This is recommended for most users.")
        elif dep_spec.level == DependencyLevel.ADVANCED:
            suggestion_parts.append("🔬 This enables advanced functionality.")
        else:  # EXPERIMENTAL
            suggestion_parts.append("🧪 This enables experimental features.")

        # Add alternatives if available
        if dep_spec.alternatives:
            alternatives_str = ", ".join(dep_spec.alternatives)
            suggestion_parts.append(f"Alternatives: {alternatives_str}")
            result.alternatives = dep_spec.alternatives

        # Add troubleshooting for common issues
        if package_name == "torch":
            suggestion_parts.append(
                "Note: PyTorch installation may require platform-specific commands."
            )
        elif package_name == "tensorflow":
            suggestion_parts.append(
                "Note: Using TensorFlow 1.x for STAGATE compatibility."
            )
        elif package_name == "rpy2":
            suggestion_parts.append("Note: Requires R to be installed on your system.")

        result.suggestion = "\n".join(suggestion_parts)
        return result

    def get_dependency_status(self) -> Dict[str, Dict[str, Any]]:
        """Get comprehensive dependency status report"""
        status_report = {
            "core": {},
            "recommended": {},
            "advanced": {},
            "experimental": {},
        }

        for name, spec in self._dependency_registry.items():
            result = self.smart_import(name, required=False)
            level_key = spec.level.value

            status_report[level_key][name] = {
                "available": result.success,
                "description": spec.description,
                "install_cmd": spec.install_cmd,
                "error": result.error if not result.success else None,
            }

        return status_report

    def get_available_features(self) -> Dict[str, List[str]]:
        """Get matrix of available features based on installed dependencies"""
        features = {
            "preprocessing": [],
            "spatial_analysis": [],
            "cell_communication": [],
            "spatial_domains": [],
            "deconvolution": [],
            "trajectory": [],
            "visualization": [],
        }

        # Check specific feature availability
        if self.smart_import("squidpy").success:
            features["spatial_analysis"].extend(
                ["spatial_neighbors", "spatial_autocorr", "co_occurrence"]
            )

        if self.smart_import("liana").success:
            features["cell_communication"].append("liana")

        if self.smart_import("cellphonedb").success:
            features["cell_communication"].append("cellphonedb")

        if self.smart_import("SpaGCN").success:
            features["spatial_domains"].append("spagcn")

        if self.smart_import("STAGATE").success:
            features["spatial_domains"].append("stagate")

        if self.smart_import("cell2location").success:
            features["deconvolution"].append("cell2location")

        if self.smart_import("scvelo").success:
            features["trajectory"].append("rna_velocity")

        return features

    def install_dependencies(
        self, level: Union[DependencyLevel, str], dry_run: bool = False
    ) -> Dict[str, bool]:
        """
        Install dependencies by level.

        Args:
            level: Dependency level to install
            dry_run: If True, only show what would be installed

        Returns:
            Dict mapping package names to installation success
        """
        if isinstance(level, str):
            level = DependencyLevel(level)

        results = {}
        packages_to_install = [
            spec for spec in self._dependency_registry.values() if spec.level == level
        ]

        if dry_run:
            print(
                f"Would install {len(packages_to_install)} {level.value} dependencies:"
            )
            for spec in packages_to_install:
                print(f"  - {spec.name}: {spec.description}")
                print(f"    Command: {spec.install_cmd}")
            return {}

        for spec in packages_to_install:
            try:
                print(f"Installing {spec.name}...")
                subprocess.run(
                    spec.install_cmd.split(), check=True, capture_output=True, text=True
                )
                results[spec.name] = True
                print(f"✅ {spec.name} installed successfully")
            except subprocess.CalledProcessError as e:
                results[spec.name] = False
                print(f"❌ Failed to install {spec.name}: {e}")

        return results


# Global instance - singleton pattern (Linus would approve - simple and practical)
_importer_instance: Optional[SmartImporter] = None


def get_smart_importer() -> SmartImporter:
    """Get the global SmartImporter instance"""
    global _importer_instance
    if _importer_instance is None:
        _importer_instance = SmartImporter()
    return _importer_instance


def smart_import(
    package_name: str, fallback_name: Optional[str] = None, required: bool = True
) -> ImportResult:
    """
    Convenience function for smart importing.

    Usage:
        result = smart_import("liana")
        if result.success:
            liana = result.module
            # Use liana
        else:
            print(result.suggestion)  # User-friendly error message
    """
    return get_smart_importer().smart_import(package_name, fallback_name, required)


def require_dependency(package_name: str, feature_name: str = None) -> Callable:
    """
    Decorator to require a dependency for a function.

    Usage:
        @require_dependency("liana", "cell communication analysis")
        async def analyze_cell_communication(...):
            import liana as li
            # Function implementation
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        async def wrapper(*args, **kwargs):
            result = smart_import(package_name)
            if not result.success:
                feature_desc = feature_name or func.__name__
                error_msg = f"Cannot perform {feature_desc}: {result.suggestion}"
                raise ImportError(error_msg)

            # Inject the module into the function's globals if it's not there
            if package_name not in func.__globals__:
                func.__globals__[package_name] = result.module

            return await func(*args, **kwargs)

        return wrapper

    return decorator


def graceful_import(
    package_name: str, fallback_func: Optional[Callable] = None
) -> Callable:
    """
    Decorator for graceful degradation when dependencies are missing.

    Usage:
        @graceful_import("enrichmap", fallback_func=basic_enrichment)
        async def advanced_enrichment(...):
            import enrichmap as em
            # Advanced implementation
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        async def wrapper(*args, **kwargs):
            result = smart_import(package_name, required=False)

            if result.success:
                # Inject the module
                func.__globals__[package_name] = result.module
                return await func(*args, **kwargs)
            else:
                # Use fallback or raise informative error
                if fallback_func:
                    context = kwargs.get("context")
                    if context:
                        await context.warning(
                            f"Using fallback for {func.__name__}: {result.suggestion}"
                        )
                    return await fallback_func(*args, **kwargs)
                else:
                    raise ImportError(
                        f"Missing dependency for {func.__name__}: {result.suggestion}"
                    )

        return wrapper

    return decorator


def check_environment() -> Dict[str, Any]:
    """
    Comprehensive environment check.
    Returns detailed report suitable for debugging.
    """
    importer = get_smart_importer()

    # Get dependency status
    dep_status = importer.get_dependency_status()

    # Get available features
    features = importer.get_available_features()

    # Count statistics
    stats = {}
    for level, deps in dep_status.items():
        available = sum(1 for dep in deps.values() if dep["available"])
        total = len(deps)
        stats[level] = f"{available}/{total}"

    # Identify critical issues
    critical_issues = []
    for name, info in dep_status["core"].items():
        if not info["available"]:
            critical_issues.append(f"Missing core dependency: {name}")

    return {
        "dependency_status": dep_status,
        "available_features": features,
        "statistics": stats,
        "critical_issues": critical_issues,
        "python_version": sys.version,
        "total_features": sum(len(f) for f in features.values()),
    }


# Backward compatibility - maintain existing interface
def check_dependency(package_name: str) -> bool:
    """Simple dependency check (backward compatible)"""
    return smart_import(package_name, required=False).success


def get_missing_dependencies(dependencies: List[str]) -> List[str]:
    """Get list of missing dependencies (backward compatible)"""
    missing = []
    for dep in dependencies:
        if not check_dependency(dep):
            missing.append(dep)
    return missing
