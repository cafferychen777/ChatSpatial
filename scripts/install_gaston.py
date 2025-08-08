#!/usr/bin/env python3
"""
GASTON Installation Helper

This script helps install GASTON for spatial gene analysis in ChatSpatial.
It provides both production and development installation options.
"""

import subprocess
import sys
import os
import argparse
from pathlib import Path


def run_command(command, capture_output=True):
    """Run a command and return the result"""
    try:
        result = subprocess.run(
            command, shell=True, capture_output=capture_output, 
            text=True, check=True
        )
        return True, result.stdout if capture_output else ""
    except subprocess.CalledProcessError as e:
        return False, e.stderr if capture_output else str(e)


def check_gaston_installed():
    """Check if GASTON is already installed"""
    try:
        import gaston
        print("✅ GASTON is already installed")
        print(f"   Version: {getattr(gaston, '__version__', 'unknown')}")
        print(f"   Location: {gaston.__file__}")
        return True
    except ImportError:
        print("❌ GASTON is not currently installed")
        return False


def install_from_pypi():
    """Install GASTON from PyPI"""
    print("🔄 Installing GASTON from PyPI...")
    success, output = run_command("pip install gaston-spatial", capture_output=False)
    
    if success:
        print("✅ GASTON installed successfully from PyPI")
        return verify_installation()
    else:
        print("❌ Failed to install GASTON from PyPI")
        print(f"Error: {output}")
        return False


def install_from_source(dev_mode=False):
    """Install GASTON from source"""
    
    # Choose directory name based on mode
    if dev_mode:
        gaston_dir = "gaston_dev"
        print("🔄 Installing GASTON for development (with comparison capability)...")
        print("   This will clone the original repository for comparison")
    else:
        gaston_dir = "GASTON"
        print("🔄 Installing GASTON from source...")
    
    # Clone repository
    if os.path.exists(gaston_dir):
        print(f"⚠️ Directory {gaston_dir} already exists")
        response = input(f"Remove existing {gaston_dir} directory? (y/N): ")
        if response.lower() == 'y':
            success, _ = run_command(f"rm -rf {gaston_dir}")
            if not success:
                print(f"❌ Failed to remove existing {gaston_dir} directory")
                return False
        else:
            print("❌ Installation cancelled")
            return False
    
    print(f"📥 Cloning GASTON repository to {gaston_dir}...")
    success, output = run_command(
        f"git clone https://github.com/Arashz/GASTON.git {gaston_dir}",
        capture_output=False
    )
    
    if not success:
        print("❌ Failed to clone GASTON repository")
        print(f"Error: {output}")
        return False
    
    # Install package
    print("📦 Installing GASTON package...")
    success, output = run_command(
        f"cd {gaston_dir} && pip install -e .",
        capture_output=False
    )
    
    if success:
        print("✅ GASTON installed successfully from source")
        if dev_mode:
            print(f"🔧 Development setup complete:")
            print(f"   • GASTON installed from: {gaston_dir}/")
            print(f"   • Original repository available for comparison")
            print(f"   • Use 'cd {gaston_dir}' to explore the source code")
        return verify_installation()
    else:
        print("❌ Failed to install GASTON from source")
        print(f"Error: {output}")
        return False


def verify_installation():
    """Verify GASTON installation"""
    print("\n🔍 Verifying GASTON installation...")
    
    try:
        # Test basic import
        import gaston
        print("✅ Basic import successful")
        
        # Test key modules
        from gaston import neural_net, spatial_gene_classification, binning_and_plotting
        from gaston import dp_related, segmented_fit, process_NN_output
        print("✅ All required modules available")
        
        print("\n🎉 GASTON installation verified successfully!")
        print("   You can now use GASTON methods in ChatSpatial")
        print("   Example: find_spatial_genes with method='gaston'")
        
        return True
        
    except ImportError as e:
        print(f"❌ Installation verification failed: {e}")
        print("   GASTON may not be properly installed")
        return False


def create_gitignore_entry():
    """Add gaston_dev to .gitignore if needed"""
    gitignore_path = Path(".gitignore")
    gaston_entry = "gaston_dev/"
    
    if not gitignore_path.exists():
        print("📝 Creating .gitignore file...")
        with open(gitignore_path, 'w') as f:
            f.write(f"{gaston_entry}\n")
        print(f"✅ Added {gaston_entry} to .gitignore")
        return
    
    # Check if entry already exists
    with open(gitignore_path, 'r') as f:
        content = f.read()
    
    if gaston_entry not in content and "gaston_dev" not in content:
        print("📝 Updating .gitignore...")
        with open(gitignore_path, 'a') as f:
            f.write(f"\n# GASTON development directory (not committed)\n{gaston_entry}\n")
        print(f"✅ Added {gaston_entry} to .gitignore")
    else:
        print("✅ .gitignore already contains gaston_dev entry")


def main():
    parser = argparse.ArgumentParser(
        description="Install GASTON for ChatSpatial spatial gene analysis"
    )
    parser.add_argument(
        "--method", 
        choices=["pypi", "source", "dev"],
        default="pypi",
        help="Installation method (default: pypi)"
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Only check if GASTON is installed"
    )
    
    args = parser.parse_args()
    
    print("🔬 ChatSpatial GASTON Installation Helper")
    print("=" * 50)
    
    # Check current installation status
    if check_gaston_installed():
        if args.check_only:
            return
        
        response = input("\\nGASTON is already installed. Reinstall? (y/N): ")
        if response.lower() != 'y':
            print("Installation cancelled")
            return
    elif args.check_only:
        return
    
    # Choose installation method
    if args.method == "pypi":
        success = install_from_pypi()
    elif args.method == "source":
        success = install_from_source(dev_mode=False)
    elif args.method == "dev":
        success = install_from_source(dev_mode=True)
        if success:
            create_gitignore_entry()
    
    if success:
        print("\\n🚀 GASTON installation completed successfully!")
        print("   You can now use GASTON methods in ChatSpatial:")
        print("   • find_spatial_genes with method='gaston'")
        print("   • Visualizations: gaston_isodepth, gaston_domains, gaston_genes")
    else:
        print("\\n❌ GASTON installation failed")
        print("   ChatSpatial will still work with other spatial gene methods:")
        print("   • spatialde, spark (available without additional setup)")
        sys.exit(1)


if __name__ == "__main__":
    main()