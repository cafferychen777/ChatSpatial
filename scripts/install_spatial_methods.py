#!/usr/bin/env python
"""
Installation helper for spatial transcriptomics methods in ChatSpatial.
"""

import os
import sys
import subprocess
import platform


def run_command(cmd, description):
    """Run a shell command and handle errors."""
    print(f"\n🔧 {description}...")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✅ {description} - Success")
            return True
        else:
            print(f"❌ {description} - Failed")
            print(f"Error: {result.stderr}")
            return False
    except Exception as e:
        print(f"❌ {description} - Exception: {e}")
        return False


def check_python_version():
    """Check Python version compatibility."""
    version = sys.version_info
    print(f"Python version: {version.major}.{version.minor}.{version.micro}")
    
    if version.major < 3 or (version.major == 3 and version.minor < 8):
        print("⚠️  Warning: Python 3.8+ is recommended")
        return False
    
    if version.major == 3 and version.minor >= 13:
        print("⚠️  Warning: Python 3.13+ may have compatibility issues with some packages")
    
    return True


def install_python_packages():
    """Install Python packages for spatial methods."""
    packages = [
        # Core dependencies
        ("numpy scipy pandas scikit-learn", "Core scientific packages"),
        ("scanpy squidpy", "Spatial transcriptomics core"),
        
        # Method-specific packages
        ("paste-bio POT", "PASTE - spatial registration"),
        ("SpatialDE", "SpatialDE - spatial variable genes"),
        ("pysal esda libpysal", "Spatial statistics"),
        
        # Optional but recommended
        ("leidenalg python-igraph", "Enhanced clustering"),
        ("rpy2", "R integration"),
        ("jupyterlab notebook", "Interactive analysis"),
    ]
    
    success_count = 0
    for pkg, desc in packages:
        if run_command(f"pip install {pkg}", f"Installing {desc}"):
            success_count += 1
    
    return success_count, len(packages)


def install_from_github():
    """Install packages from GitHub."""
    github_packages = [
        {
            "name": "STAGATE",
            "path": "third_party/STAGATE",
            "url": "https://github.com/QIFEIDKN/STAGATE.git",
            "install_cmd": "pip install -e .",
            "pre_install": "pip install torch scanpy pandas numpy"
        },
        {
            "name": "BANKSY",
            "path": "third_party/BANKSY_py",
            "url": "https://github.com/prabhakarlab/Banksy_py.git",
            "install_cmd": "pip install -e .",
            "pre_install": None
        }
    ]
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    success_count = 0
    
    for pkg in github_packages:
        pkg_path = os.path.join(base_dir, pkg["path"])
        
        # Clone if not exists
        if not os.path.exists(pkg_path):
            print(f"\n📦 Cloning {pkg['name']}...")
            parent_dir = os.path.dirname(pkg_path)
            os.makedirs(parent_dir, exist_ok=True)
            
            if run_command(f"git clone {pkg['url']} {pkg_path}", f"Cloning {pkg['name']}"):
                success_count += 0.5
        
        # Install
        if os.path.exists(pkg_path):
            # Pre-install dependencies if needed
            if pkg["pre_install"]:
                run_command(pkg["pre_install"], f"Installing {pkg['name']} dependencies")
            
            # Install package
            original_dir = os.getcwd()
            os.chdir(pkg_path)
            
            if run_command(pkg["install_cmd"], f"Installing {pkg['name']}"):
                success_count += 0.5
            
            os.chdir(original_dir)
    
    return success_count, len(github_packages)


def fix_compatibility_issues():
    """Fix known compatibility issues."""
    fixes_applied = 0
    
    # Fix SpatialDE scipy issue
    try:
        import SpatialDE
        spatialDE_base = os.path.join(os.path.dirname(SpatialDE.__file__), 'base.py')
        
        with open(spatialDE_base, 'r') as f:
            content = f.read()
        
        if 'from scipy.misc import derivative' in content:
            new_content = content.replace(
                'from scipy.misc import derivative',
                '''try:
    from scipy.misc import derivative
except ImportError:
    from scipy.optimize import approx_fprime
    def derivative(func, x0, dx=1e-8):
        return approx_fprime([x0], lambda x: func(x[0]), [dx])[0]'''
            )
            
            with open(spatialDE_base, 'w') as f:
                f.write(new_content)
            
            print("✅ Fixed SpatialDE scipy compatibility")
            fixes_applied += 1
    except:
        pass
    
    return fixes_applied


def install_r_packages():
    """Install R packages for methods that require them."""
    print("\n📊 R Package Installation")
    print("For R-based methods (SPOTlight, BayesSpace, SPARK), run these commands in R:")
    print("""
    # Install BiocManager if not already installed
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    # Install spatial packages
    BiocManager::install(c("SPOTlight", "SingleCellExperiment"))
    install.packages(c("SPARK", "Seurat"))
    
    # For BayesSpace (if needed)
    devtools::install_github("edward130603/BayesSpace")
    """)


def test_installations():
    """Test if packages can be imported."""
    print("\n🧪 Testing installations...")
    
    test_imports = [
        ("paste", "PASTE"),
        ("SpatialDE", "SpatialDE"),
        ("STAGATE", "STAGATE"),
        ("banksy", "BANKSY"),
        ("rpy2", "rpy2 (R integration)"),
        ("squidpy", "Squidpy"),
        ("scanpy", "Scanpy"),
    ]
    
    working = []
    failed = []
    
    for module, name in test_imports:
        try:
            __import__(module)
            working.append(name)
        except ImportError as e:
            failed.append((name, str(e)))
    
    print(f"\n✅ Working: {', '.join(working)}")
    
    if failed:
        print("\n❌ Failed imports:")
        for name, error in failed:
            print(f"  - {name}: {error}")
    
    return len(working), len(test_imports)


def main():
    """Main installation routine."""
    print("🚀 ChatSpatial Spatial Methods Installation Helper")
    print("=" * 50)
    
    # Check Python version
    check_python_version()
    
    # Install Python packages
    print("\n1️⃣ Installing Python packages...")
    py_success, py_total = install_python_packages()
    
    # Install GitHub packages
    print("\n2️⃣ Installing packages from GitHub...")
    gh_success, gh_total = install_from_github()
    
    # Fix compatibility issues
    print("\n3️⃣ Fixing compatibility issues...")
    fixes = fix_compatibility_issues()
    
    # R packages info
    install_r_packages()
    
    # Test installations
    test_success, test_total = test_installations()
    
    # Summary
    print("\n" + "=" * 50)
    print("📋 Installation Summary:")
    print(f"  Python packages: {py_success}/{py_total} successful")
    print(f"  GitHub packages: {gh_success}/{gh_total} successful")
    print(f"  Compatibility fixes: {fixes} applied")
    print(f"  Import tests: {test_success}/{test_total} working")
    
    print("\n💡 Next steps:")
    print("1. Install R packages (see above) for SPOTlight/SPARK support")
    print("2. Run 'python tests/run_new_methods_test.py' to test methods")
    print("3. Check documentation for usage examples")


if __name__ == "__main__":
    main()