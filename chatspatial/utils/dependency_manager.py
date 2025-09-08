"""
ChatSpatial Dependency Manager

This module provides intelligent dependency detection and graceful fallbacks
for optional dependencies. Following Linus's philosophy of "good taste":
- Eliminate special cases
- Clear error messages 
- Fail fast when needed, degrade gracefully when possible
"""

import sys
import warnings
import importlib
import pkg_resources
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass
from pathlib import Path


@dataclass
class DependencyInfo:
    """Information about a dependency"""
    name: str
    import_name: str
    version_constraint: str
    install_command: str
    description: str
    is_critical: bool = False
    python_version_constraint: Optional[str] = None
    alternatives: List[str] = None
    
    def __post_init__(self):
        if self.alternatives is None:
            self.alternatives = []


class DependencyError(Exception):
    """Raised when a critical dependency is missing"""
    
    def __init__(self, dependency: str, feature: str, install_command: str):
        self.dependency = dependency
        self.feature = feature
        self.install_command = install_command
        super().__init__(
            f"Missing dependency '{dependency}' for feature '{feature}'. "
            f"Install with: {install_command}"
        )


class DependencyManager:
    """
    Centralized dependency management following the UNIX philosophy:
    - Do one thing well: manage dependencies
    - Fail fast: clear error messages for critical deps
    - Degrade gracefully: optional features work without deps
    """
    
    # Core dependencies that must be available
    CRITICAL_DEPS = {
        'numpy': DependencyInfo(
            'numpy', 'numpy', '>=1.21.0,<2.0', 'pip install numpy>=1.21.0',
            'Core numerical computing'
        ),
        'pandas': DependencyInfo(
            'pandas', 'pandas', '>=1.3.0,<3.0', 'pip install pandas>=1.3.0',
            'Data manipulation and analysis'
        ),
        'scanpy': DependencyInfo(
            'scanpy', 'scanpy', '>=1.9.0,<2.0', 'pip install scanpy>=1.9.0',
            'Single-cell analysis toolkit'
        ),
        'anndata': DependencyInfo(
            'anndata', 'anndata', '>=0.8.0,<1.0', 'pip install anndata>=0.8.0',
            'Annotated data matrices'
        ),
        'scipy': DependencyInfo(
            'scipy', 'scipy', '>=1.7.0,<2.0', 'pip install scipy>=1.7.0',
            'Scientific computing'
        ),
        'scikit-learn': DependencyInfo(
            'scikit-learn', 'sklearn', '>=1.0.0,<2.0', 'pip install scikit-learn>=1.0.0',
            'Machine learning toolkit'
        ),
        'matplotlib': DependencyInfo(
            'matplotlib', 'matplotlib', '>=3.5.0,<4.0', 'pip install matplotlib>=3.5.0',
            'Plotting and visualization'
        ),
        'seaborn': DependencyInfo(
            'seaborn', 'seaborn', '>=0.11.0,<1.0', 'pip install seaborn>=0.11.0',
            'Statistical data visualization'
        ),
    }
    
    # Optional dependencies with clear version constraints
    OPTIONAL_DEPS = {
        # Deep learning frameworks - version critical
        'torch': DependencyInfo(
            'torch', 'torch', '>=2.0.0,<3.0', 
            'pip install torch>=2.0.0 --index-url https://download.pytorch.org/whl/cpu',
            'PyTorch deep learning framework',
            python_version_constraint='>=3.8,<3.13'
        ),
        'tensorflow': DependencyInfo(
            'tensorflow', 'tensorflow', '>=2.8.0,<3.0',
            'pip install tensorflow>=2.8.0',
            'TensorFlow deep learning framework (TF 1.x not supported)',
            python_version_constraint='>=3.8,<3.13'
        ),
        
        # Spatial transcriptomics tools
        'squidpy': DependencyInfo(
            'squidpy', 'squidpy', '>=1.2.0,<2.0',
            'pip install squidpy>=1.2.0',
            'Spatial transcriptomics analysis toolkit'
        ),
        'cell2location': DependencyInfo(
            'cell2location', 'cell2location', '>=0.1.3,<1.0',
            'pip install cell2location>=0.1.3',
            'Spatial deconvolution with Bayesian models',
            alternatives=['tangram', 'stereoscope']
        ),
        'scvi-tools': DependencyInfo(
            'scvi-tools', 'scvi', '>=1.0.0,<2.0',
            'pip install scvi-tools>=1.0.0',
            'Deep generative models for single-cell omics',
            python_version_constraint='>=3.8,<3.12'
        ),
        
        # Spatial domain identification
        'SpaGCN': DependencyInfo(
            'SpaGCN', 'SpaGCN', '>=1.2.5,<2.0',
            'pip install SpaGCN>=1.2.5',
            'Spatial domain identification using graph convolutional networks'
        ),
        'STAGATE': DependencyInfo(
            'STAGATE', 'STAGATE', '>=1.0.0,<2.0',
            'git clone https://github.com/QIFEIDKN/STAGATE.git && cd STAGATE && python setup.py install',
            'Spatial domain identification using attention mechanisms',
            alternatives=['SpaGCN', 'BayesSpace']
        ),
        
        # RNA velocity and trajectory
        'scvelo': DependencyInfo(
            'scvelo', 'scvelo', '>=0.2.5,<1.0',
            'pip install scvelo>=0.2.5',
            'RNA velocity analysis'
        ),
        'cellrank': DependencyInfo(
            'cellrank', 'cellrank', '>=2.0.0,<3.0',
            'pip install cellrank>=2.0.0',
            'Trajectory analysis and cellular dynamics'
        ),
        
        # Spatial genes
        'spatialde': DependencyInfo(
            'spatialde', 'spatialde', '>=1.1.3,<2.0',
            'pip install spatialde>=1.1.3',
            'Spatial variable gene identification'
        ),
        
        # Cell communication
        'liana': DependencyInfo(
            'liana', 'liana', '>=0.1.13,<1.0',
            'pip install liana>=0.1.13',
            'Ligand-receptor interaction analysis'
        ),
        'cellphonedb': DependencyInfo(
            'cellphonedb', 'cellphonedb', '>=5.0.0,<6.0',
            'pip install cellphonedb>=5.0.0',
            'Cell-cell communication inference',
            alternatives=['liana', 'squidpy']
        ),
        
        # Integration methods
        'harmonypy': DependencyInfo(
            'harmonypy', 'harmony', '>=0.0.9,<1.0',
            'pip install harmonypy>=0.0.9',
            'Batch effect correction'
        ),
        'bbknn': DependencyInfo(
            'bbknn', 'bbknn', '>=1.5.0,<2.0',
            'pip install bbknn>=1.5.0',
            'Batch balanced k-nearest neighbors'
        ),
        'scanorama': DependencyInfo(
            'scanorama', 'scanorama', '>=1.7.0,<2.0',
            'pip install scanorama>=1.7.0',
            'Panoramic stitching of single-cell data'
        ),
        
        # R interface (with clear warnings)
        'rpy2': DependencyInfo(
            'rpy2', 'rpy2', '>=3.4.0,<4.0',
            'pip install rpy2>=3.4.0 (requires R installation)',
            'R interface for advanced statistical methods',
            python_version_constraint='>=3.8,<3.12'
        ),
        
        # Visualization
        'plotly': DependencyInfo(
            'plotly', 'plotly', '>=5.0.0,<6.0',
            'pip install plotly>=5.0.0',
            'Interactive plotting'
        ),
        
        # Performance 
        'dask': DependencyInfo(
            'dask', 'dask', '>=2024.1.0,<2025.0',
            'pip install dask>=2024.1.0',
            'Parallel computing and out-of-core arrays'
        ),
        'numba': DependencyInfo(
            'numba', 'numba', '>=0.57.0,<1.0',
            'pip install numba>=0.57.0',
            'JIT compilation for numerical functions'
        ),
    }
    
    def __init__(self):
        self._cache: Dict[str, Tuple[bool, Optional[str]]] = {}
        self._python_version = sys.version_info
        self._check_python_compatibility()
    
    def _check_python_compatibility(self):
        """Check if current Python version is supported"""
        if self._python_version < (3, 8):
            raise RuntimeError(
                f"Python {self._python_version.major}.{self._python_version.minor} "
                "is not supported. ChatSpatial requires Python >= 3.8"
            )
        
        if self._python_version >= (3, 13):
            warnings.warn(
                f"Python {self._python_version.major}.{self._python_version.minor} "
                "is not fully tested. Some dependencies may not be available.",
                UserWarning
            )
    
    def check_dependency(self, dep_name: str) -> Tuple[bool, Optional[str], Optional[Any]]:
        """
        Check if a dependency is available and return module if found
        
        Returns:
            (is_available, version_or_error, module_or_none)
        """
        if dep_name in self._cache:
            available, version = self._cache[dep_name]
            if available:
                try:
                    module = importlib.import_module(self._get_import_name(dep_name))
                    return True, version, module
                except ImportError:
                    return False, "Import failed after cache hit", None
            return False, version, None
        
        dep_info = self._get_dependency_info(dep_name)
        if not dep_info:
            self._cache[dep_name] = (False, f"Unknown dependency: {dep_name}")
            return False, f"Unknown dependency: {dep_name}", None
        
        # Check Python version compatibility
        if dep_info.python_version_constraint:
            if not self._check_python_constraint(dep_info.python_version_constraint):
                error_msg = (f"{dep_name} requires Python {dep_info.python_version_constraint}, "
                           f"but you have {self._python_version.major}.{self._python_version.minor}")
                self._cache[dep_name] = (False, error_msg)
                return False, error_msg, None
        
        try:
            module = importlib.import_module(dep_info.import_name)
            
            # Check version if available
            version = None
            try:
                if hasattr(module, '__version__'):
                    version = module.__version__
                else:
                    version = pkg_resources.get_distribution(dep_info.name).version
            except Exception:
                version = "unknown"
            
            self._cache[dep_name] = (True, version)
            return True, version, module
            
        except ImportError as e:
            error_msg = f"Import failed: {str(e)}"
            self._cache[dep_name] = (False, error_msg)
            return False, error_msg, None
    
    def _get_dependency_info(self, dep_name: str) -> Optional[DependencyInfo]:
        """Get dependency info from critical or optional deps"""
        if dep_name in self.CRITICAL_DEPS:
            return self.CRITICAL_DEPS[dep_name]
        if dep_name in self.OPTIONAL_DEPS:
            return self.OPTIONAL_DEPS[dep_name]
        return None
    
    def _get_import_name(self, dep_name: str) -> str:
        """Get the import name for a dependency"""
        dep_info = self._get_dependency_info(dep_name)
        return dep_info.import_name if dep_info else dep_name
    
    def _check_python_constraint(self, constraint: str) -> bool:
        """Check if current Python version satisfies constraint"""
        # Simple constraint parsing for >=x.y,<z.w format
        if '>=' in constraint and '<' in constraint:
            parts = constraint.split(',')
            min_version = parts[0].replace('>=', '').strip()
            max_version = parts[1].replace('<', '').strip()
            
            min_tuple = tuple(map(int, min_version.split('.')))
            max_tuple = tuple(map(int, max_version.split('.')))
            
            return min_tuple <= self._python_version[:2] < max_tuple
        
        return True  # Conservative fallback
    
    def require_dependency(self, dep_name: str, feature: str = None) -> Any:
        """
        Require a dependency and raise DependencyError if not available
        
        Args:
            dep_name: Name of the dependency
            feature: Feature that requires this dependency (for error message)
            
        Returns:
            The imported module
            
        Raises:
            DependencyError: If dependency is not available
        """
        available, error_or_version, module = self.check_dependency(dep_name)
        
        if not available:
            dep_info = self._get_dependency_info(dep_name)
            install_cmd = dep_info.install_command if dep_info else f"pip install {dep_name}"
            feature = feature or dep_name
            raise DependencyError(dep_name, feature, install_cmd)
        
        return module
    
    def try_import(self, dep_name: str, feature: str = None) -> Tuple[Optional[Any], Optional[str]]:
        """
        Try to import a dependency, return None and warning if not available
        
        Args:
            dep_name: Name of the dependency
            feature: Feature that would use this dependency
            
        Returns:
            (module_or_none, warning_message_or_none)
        """
        available, error_or_version, module = self.check_dependency(dep_name)
        
        if not available:
            dep_info = self._get_dependency_info(dep_name)
            feature = feature or dep_name
            warning = (f"Optional dependency '{dep_name}' not available for {feature}. "
                      f"Install with: {dep_info.install_command if dep_info else f'pip install {dep_name}'}")
            
            # Suggest alternatives if available
            if dep_info and dep_info.alternatives:
                warning += f" Alternatives: {', '.join(dep_info.alternatives)}"
            
            return None, warning
        
        return module, None
    
    def get_available_methods(self, method_deps: Dict[str, List[str]]) -> List[str]:
        """
        Get list of methods that have all required dependencies available
        
        Args:
            method_deps: Dict mapping method names to required dependencies
            
        Returns:
            List of available method names
        """
        available = []
        for method, deps in method_deps.items():
            if all(self.check_dependency(dep)[0] for dep in deps):
                available.append(method)
        return available
    
    def get_dependency_report(self) -> Dict[str, Any]:
        """Generate a comprehensive dependency report"""
        report = {
            'python_version': f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
            'critical_dependencies': {},
            'optional_dependencies': {},
            'missing_critical': [],
            'missing_optional': [],
            'version_conflicts': [],
        }
        
        # Check critical dependencies
        for dep_name in self.CRITICAL_DEPS:
            available, version_or_error, _ = self.check_dependency(dep_name)
            report['critical_dependencies'][dep_name] = {
                'available': available,
                'version_or_error': version_or_error
            }
            if not available:
                report['missing_critical'].append(dep_name)
        
        # Check optional dependencies
        for dep_name in self.OPTIONAL_DEPS:
            available, version_or_error, _ = self.check_dependency(dep_name)
            report['optional_dependencies'][dep_name] = {
                'available': available,
                'version_or_error': version_or_error
            }
            if not available:
                report['missing_optional'].append(dep_name)
        
        return report
    
    def print_dependency_report(self):
        """Print a human-readable dependency report"""
        report = self.get_dependency_report()
        
        print("ChatSpatial Dependency Report")
        print("=" * 40)
        print(f"Python version: {report['python_version']}")
        print()
        
        print("Critical Dependencies:")
        for dep, info in report['critical_dependencies'].items():
            status = "✓" if info['available'] else "✗"
            print(f"  {status} {dep}: {info['version_or_error']}")
        
        print("\nOptional Dependencies:")
        for dep, info in report['optional_dependencies'].items():
            status = "✓" if info['available'] else "✗"
            print(f"  {status} {dep}: {info['version_or_error']}")
        
        if report['missing_critical']:
            print(f"\n⚠️  CRITICAL: Missing {len(report['missing_critical'])} critical dependencies!")
            for dep in report['missing_critical']:
                dep_info = self.CRITICAL_DEPS[dep]
                print(f"   Install: {dep_info.install_command}")
        
        if report['missing_optional']:
            print(f"\n💡 Optional: {len(report['missing_optional'])} optional dependencies missing")
            print("   Some advanced features will be unavailable")


# Global instance
dependency_manager = DependencyManager()


def require(dep_name: str, feature: str = None):
    """Convenience function to require a dependency"""
    return dependency_manager.require_dependency(dep_name, feature)


def try_import(dep_name: str, feature: str = None):
    """Convenience function to try importing a dependency"""
    return dependency_manager.try_import(dep_name, feature)


def check(dep_name: str):
    """Convenience function to check if a dependency is available"""
    return dependency_manager.check_dependency(dep_name)[0]


def get_available_methods(method_deps: Dict[str, List[str]]) -> List[str]:
    """Convenience function to get available methods"""
    return dependency_manager.get_available_methods(method_deps)


# Export commonly used functions
__all__ = [
    'DependencyManager', 'DependencyError', 'dependency_manager',
    'require', 'try_import', 'check', 'get_available_methods'
]