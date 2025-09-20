"""
User guidance and installation help system.

This module provides friendly, actionable guidance to users when they encounter
dependency issues. Instead of cryptic ImportError messages, users get clear
instructions on how to resolve problems.

Design Philosophy:
- Never leave users stranded with unhelpful errors
- Provide specific, actionable instructions
- Guide users to the right solution quickly
- Make dependency management a competitive advantage

Author: Linus-approved user experience
"""

import platform
import sys
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional

# Simple replacement for DependencyLevel
class DependencyLevel(Enum):
    CORE = "core"
    RECOMMENDED = "recommended" 
    ADVANCED = "advanced"
    EXPERIMENTAL = "experimental"


class InstallationMethod(Enum):
    """Available installation methods"""

    PIP = "pip"
    CONDA = "conda"
    PIPX = "pipx"
    SYSTEM = "system"


@dataclass
class InstallationGuide:
    """Installation guidance for a dependency"""

    package_name: str
    description: str
    methods: Dict[InstallationMethod, str]  # Method -> command mapping
    prerequisites: List[str] = None
    platform_notes: Dict[str, str] = None
    troubleshooting: List[str] = None
    estimated_time: str = "< 2 minutes"
    difficulty: str = "Easy"


class UserGuidanceSystem:
    """
    Intelligent guidance system for dependency installation and troubleshooting.

    This replaces the chaos of scattered error messages with a unified,
    helpful user experience.
    """

    def __init__(self):
        self.platform = platform.system().lower()
        self.python_version = sys.version_info
        self.guides = self._build_installation_guides()

    def _build_installation_guides(self) -> Dict[str, InstallationGuide]:
        """Build comprehensive installation guides for all dependencies"""
        guides = {}

        # Core dependencies (should rarely fail, but just in case)
        guides["numpy"] = InstallationGuide(
            package_name="numpy",
            description="Numerical computing library (core dependency)",
            methods={
                InstallationMethod.PIP: "pip install numpy",
                InstallationMethod.CONDA: "conda install numpy",
            },
            troubleshooting=[
                "If installation fails, try upgrading pip: pip install --upgrade pip",
                "On M1 Macs, ensure you're using Python with Apple Silicon support",
            ],
        )

        guides["pandas"] = InstallationGuide(
            package_name="pandas",
            description="Data manipulation library (core dependency)",
            methods={
                InstallationMethod.PIP: "pip install pandas",
                InstallationMethod.CONDA: "conda install pandas",
            },
        )

        # Spatial analysis dependencies
        guides["squidpy"] = InstallationGuide(
            package_name="squidpy",
            description="Spatial transcriptomics analysis toolkit",
            methods={
                InstallationMethod.PIP: "pip install squidpy",
                InstallationMethod.CONDA: "conda install -c conda-forge squidpy",
            },
            prerequisites=["scanpy", "anndata"],
            platform_notes={
                "windows": "May require Microsoft Visual C++ Build Tools",
                "macos": "Works best with Homebrew Python or conda",
            },
            troubleshooting=[
                "If import fails, check dask version: pip install 'dask<2025'",
                "Clear pip cache if installation stalls: pip cache purge",
            ],
            estimated_time="3-5 minutes",
        )

        # Advanced method dependencies
        guides["liana"] = InstallationGuide(
            package_name="liana",
            description="Cell-cell communication analysis with LIANA+",
            methods={
                InstallationMethod.PIP: "pip install liana",
                InstallationMethod.CONDA: "conda install -c bioconda liana",
            },
            prerequisites=["scanpy", "anndata", "pandas"],
            troubleshooting=[
                "If installation is slow, use: pip install liana --no-cache-dir",
                "For development version: pip install git+https://github.com/saezlab/liana-py.git",
            ],
            estimated_time="2-4 minutes",
            difficulty="Easy",
        )

        guides["cellphonedb"] = InstallationGuide(
            package_name="cellphonedb",
            description="Cell-cell communication analysis with CellPhoneDB",
            methods={InstallationMethod.PIP: "pip install cellphonedb"},
            prerequisites=["pandas", "scipy", "scikit-learn"],
            troubleshooting=[
                "CellPhoneDB has complex dependencies, consider using liana instead",
                "If installation fails, try: pip install cellphonedb --no-deps first",
            ],
            estimated_time="3-6 minutes",
            difficulty="Moderate",
        )

        guides["cell2location"] = InstallationGuide(
            package_name="cell2location",
            description="Spatial deconvolution using cell2location",
            methods={InstallationMethod.PIP: "pip install cell2location"},
            prerequisites=["torch", "scanpy", "anndata"],
            platform_notes={"all": "Requires PyTorch - large download (~1GB)"},
            troubleshooting=[
                "Install PyTorch first: pip install torch",
                "For GPU support, install CUDA-enabled PyTorch first",
                "If memory issues during install, close other applications",
            ],
            estimated_time="5-10 minutes",
            difficulty="Moderate",
        )

        guides["SpaGCN"] = InstallationGuide(
            package_name="SpaGCN",
            description="Spatial domain identification with SpaGCN",
            methods={InstallationMethod.PIP: "pip install SpaGCN"},
            prerequisites=["torch", "scanpy", "opencv-python"],
            platform_notes={
                "windows": "May require Visual Studio Build Tools",
                "macos": "OpenCV installation may take extra time",
            },
            troubleshooting=[
                "Install prerequisites first: pip install torch opencv-python",
                "If OpenCV fails, try: pip install opencv-python-headless",
            ],
            estimated_time="5-8 minutes",
            difficulty="Moderate",
        )

        guides["STAGATE"] = InstallationGuide(
            package_name="STAGATE",
            description="Spatial domain identification with STAGATE",
            methods={
                InstallationMethod.PIP: "git clone https://github.com/QIFEIDKN/STAGATE.git && cd STAGATE && python setup.py install"
            },
            prerequisites=["tensorflow>=1.15.0,<2.0", "scanpy"],
            platform_notes={
                "all": "Requires TensorFlow 1.x (older version for compatibility)"
            },
            troubleshooting=[
                "Install TensorFlow 1.x first: pip install 'tensorflow>=1.15.0,<2.0'",
                "TensorFlow 1.x may conflict with other packages",
                "Consider using alternative methods like clustering-based spatial domains",
            ],
            estimated_time="4-7 minutes",
            difficulty="Moderate to Hard",
        )

        guides["rpy2"] = InstallationGuide(
            package_name="rpy2",
            description="R interface for sc-type and other R methods",
            methods={InstallationMethod.PIP: "pip install rpy2"},
            prerequisites=["R installation"],
            platform_notes={
                "windows": "Install R from https://cran.r-project.org/bin/windows/base/",
                "macos": "Install R from https://cran.r-project.org/bin/macosx/ or brew install r",
                "linux": "Install R using package manager: sudo apt-get install r-base-dev",
            },
            troubleshooting=[
                "Ensure R is installed and in your PATH",
                "On Linux, install R development headers: sudo apt-get install r-base-dev",
                "Set R_HOME environment variable if needed",
            ],
            estimated_time="5-15 minutes",
            difficulty="Hard",
        )

        return guides

    def get_installation_guide(self, package_name: str) -> Optional[InstallationGuide]:
        """Get installation guide for a package"""
        return self.guides.get(package_name)

    def generate_quick_fix_command(self, package_name: str) -> str:
        """Generate the quickest installation command for a package"""
        guide = self.get_installation_guide(package_name)

        if guide and InstallationMethod.PIP in guide.methods:
            return guide.methods[InstallationMethod.PIP]
        else:
            # Generic fallback
            return f"pip install {package_name}"

    def generate_comprehensive_help(self, package_name: str) -> str:
        """Generate comprehensive help text for a missing dependency"""
        guide = self.get_installation_guide(package_name)

        if not guide:
            # Generic help for unknown packages
            return f"""
Missing dependency: {package_name}

Quick fix: pip install {package_name}

If this doesn't work:
1. Update pip: pip install --upgrade pip
2. Try with no cache: pip install {package_name} --no-cache-dir
3. Use conda: conda install {package_name}

For more help, consult the package documentation or file an issue.
"""

        # Build comprehensive help
        help_sections = []

        # Header
        help_sections.append(f"📦 Missing: {guide.package_name}")
        help_sections.append(f"   {guide.description}")
        help_sections.append(f"   Estimated install time: {guide.estimated_time}")
        help_sections.append(f"   Difficulty: {guide.difficulty}")
        help_sections.append("")

        # Quick fix
        help_sections.append("🚀 Quick Fix:")
        help_sections.append(f"   {self.generate_quick_fix_command(package_name)}")
        help_sections.append("")

        # Prerequisites
        if guide.prerequisites:
            help_sections.append("📋 Prerequisites:")
            for prereq in guide.prerequisites:
                help_sections.append(f"   • {prereq}")
            help_sections.append("")

        # Installation methods
        help_sections.append("💿 Installation Options:")
        for method, command in guide.methods.items():
            help_sections.append(f"   {method.value}: {command}")
        help_sections.append("")

        # Platform-specific notes
        if guide.platform_notes and self.platform in guide.platform_notes:
            help_sections.append(f"🖥️  {self.platform.title()} Note:")
            help_sections.append(f"   {guide.platform_notes[self.platform]}")
            help_sections.append("")

        # Troubleshooting
        if guide.troubleshooting:
            help_sections.append("🔧 Troubleshooting:")
            for i, tip in enumerate(guide.troubleshooting, 1):
                help_sections.append(f"   {i}. {tip}")
            help_sections.append("")

        return "\n".join(help_sections)

    def check_system_health(self) -> Dict[str, Any]:
        """Check system health for dependency installation"""
        health_report = {
            "python_version": f"{self.python_version.major}.{self.python_version.minor}.{self.python_version.micro}",
            "platform": self.platform,
            "issues": [],
            "recommendations": [],
        }

        # Check Python version
        if self.python_version < (3, 8):
            health_report["issues"].append(
                f"Python {self.python_version.major}.{self.python_version.minor} is too old"
            )
            health_report["recommendations"].append("Upgrade to Python 3.8+")
        elif self.python_version >= (3, 12):
            health_report["issues"].append(
                "Python 3.12+ may have compatibility issues with some packages"
            )
            health_report["recommendations"].append(
                "Consider using Python 3.10 or 3.11 for maximum compatibility"
            )

        # Check pip
        try:
            import pip

            health_report["pip_available"] = True
        except ImportError:
            health_report["pip_available"] = False
            health_report["issues"].append("pip not available")
            health_report["recommendations"].append("Install pip")

        # Check virtual environment
        if hasattr(sys, "real_prefix") or (
            hasattr(sys, "base_prefix") and sys.base_prefix != sys.prefix
        ):
            health_report["virtual_env"] = True
        else:
            health_report["virtual_env"] = False
            health_report["recommendations"].append(
                "Consider using a virtual environment"
            )

        # Platform-specific checks
        if self.platform == "windows":
            health_report["recommendations"].append(
                "Consider installing Microsoft Visual C++ Build Tools for complex packages"
            )
        elif self.platform == "darwin" and platform.machine() == "arm64":
            health_report["recommendations"].append(
                "Ensure packages support Apple Silicon"
            )

        return health_report

    def generate_installation_script(
        self,
        level: DependencyLevel,
        method: InstallationMethod = InstallationMethod.PIP,
    ) -> str:
        """Generate installation script for a dependency level"""

        # Simple package definitions without smart_import
        package_defs = {
            DependencyLevel.CORE: {
                "numpy": "Core numerical computing",
                "pandas": "Data manipulation",
                "scanpy": "Single-cell analysis",
                "anndata": "Annotated data",
            },
            DependencyLevel.ADVANCED: {
                "scvi-tools": "Deep learning",
                "liana": "Cell communication",
            },
        }
        
        packages_by_level = package_defs.get(level, {})

        script_lines = [
            "#!/bin/bash",
            f"# ChatSpatial {level.value} dependencies installation script",
            f"# Generated for {self.platform} platform",
            "",
            f"echo 'Installing ChatSpatial {level.value} dependencies...'",
            "",
        ]

        for package_name, dep_spec in packages_by_level.items():
            guide = self.get_installation_guide(package_name)

            if guide and method in guide.methods:
                cmd = guide.methods[method]
            else:
                cmd = dep_spec.install_cmd

            script_lines.extend(
                [
                    f"echo 'Installing {package_name}...'",
                    f"{cmd}",
                    "if [ $? -ne 0 ]; then",
                    f"    echo 'Failed to install {package_name}'",
                    f"    echo 'Try: {self.generate_quick_fix_command(package_name)}'",
                    "fi",
                    "",
                ]
            )

        script_lines.append("echo 'Installation complete!'")

        return "\n".join(script_lines)

    def create_requirements_file(
        self, level: DependencyLevel, output_file: Optional[str] = None
    ) -> str:
        """Create requirements.txt file for a dependency level"""

        # Simple package lists without smart_import
        package_lists = {
            DependencyLevel.CORE: [
                "numpy", "pandas", "scanpy", "anndata",
                "matplotlib", "scipy", "scikit-learn"
            ],
            DependencyLevel.ADVANCED: [
                "scvi-tools", "liana", "cellphonedb", "cellrank"
            ],
        }
        
        packages = package_lists.get(level, [])

        requirements_content = "\n".join(sorted(packages))

        if output_file:
            with open(output_file, "w") as f:
                f.write(requirements_content)
            return f"Requirements saved to {output_file}"
        else:
            return requirements_content


# Global instance
_guidance_system: Optional[UserGuidanceSystem] = None


def get_guidance_system() -> UserGuidanceSystem:
    """Get the global guidance system instance"""
    global _guidance_system
    if _guidance_system is None:
        _guidance_system = UserGuidanceSystem()
    return _guidance_system


def get_help_for_missing_package(package_name: str) -> str:
    """Get comprehensive help for a missing package"""
    return get_guidance_system().generate_comprehensive_help(package_name)


def get_quick_install_command(package_name: str) -> str:
    """Get quick installation command for a package"""
    return get_guidance_system().generate_quick_fix_command(package_name)


def check_installation_health() -> Dict[str, Any]:
    """Check system health for package installation"""
    return get_guidance_system().check_system_health()


# Backward compatibility functions
def suggest_installation_method(package_name: str) -> str:
    """Legacy function - use get_quick_install_command instead"""
    return get_quick_install_command(package_name)


def get_platform_specific_notes(package_name: str) -> Optional[str]:
    """Get platform-specific installation notes"""
    guide = get_guidance_system().get_installation_guide(package_name)
    if guide and guide.platform_notes:
        platform_name = get_guidance_system().platform
        return guide.platform_notes.get(platform_name)
    return None
