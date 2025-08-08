#!/usr/bin/env python3
"""
Download 10x Genomics Harmony example datasets
"""

import os
import urllib.request
import zipfile
import tarfile
import gzip
import shutil
from pathlib import Path

def download_file(url, filename):
    """Download a file from URL"""
    print(f"Downloading {filename}...")
    try:
        urllib.request.urlretrieve(url, filename)
        print(f"✓ Downloaded {filename}")
        return True
    except Exception as e:
        print(f"✗ Failed to download {filename}: {e}")
        return False

def extract_archive(archive_path, extract_to):
    """Extract archive file"""
    print(f"Extracting {archive_path}...")
    try:
        if archive_path.endswith('.tar.gz'):
            with tarfile.open(archive_path, 'r:gz') as tar:
                tar.extractall(extract_to)
        elif archive_path.endswith('.zip'):
            with zipfile.ZipFile(archive_path, 'r') as zip_ref:
                zip_ref.extractall(extract_to)
        elif archive_path.endswith('.gz') and not archive_path.endswith('.tar.gz'):
            with gzip.open(archive_path, 'rb') as f_in:
                output_path = archive_path[:-3]  # Remove .gz extension
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        print(f"✓ Extracted {archive_path}")
        return True
    except Exception as e:
        print(f"✗ Failed to extract {archive_path}: {e}")
        return False

def main():
    """Download Harmony example datasets"""
    
    # Create data directory
    data_dir = Path("harmony_datasets")
    data_dir.mkdir(exist_ok=True)
    
    print("🔄 Downloading Harmony Example Datasets")
    print("=" * 50)
    
    # 10x Genomics dataset URLs (these are examples - actual URLs may require form submission)
    datasets = {
        "jurkat_293t_50_50": {
            "name": "Jurkat:293T 50:50 Mixture",
            "description": "50:50 mixture of Jurkat and HEK293T cells",
            "files": [
                "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/jurkat:293t_50:50/jurkat:293t_50:50_filtered_feature_bc_matrix.h5",
                "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/jurkat:293t_50:50/jurkat:293t_50:50_filtered_feature_bc_matrix.tar.gz",
            ]
        },
        "jurkat": {
            "name": "Pure Jurkat Cells",
            "description": "Pure Jurkat cell line",
            "files": [
                "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/jurkat/jurkat_filtered_feature_bc_matrix.h5",
                "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/jurkat/jurkat_filtered_feature_bc_matrix.tar.gz",
            ]
        },
        "293t": {
            "name": "Pure 293T Cells", 
            "description": "Pure HEK293T cell line",
            "files": [
                "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/293t_filtered_feature_bc_matrix.h5",
                "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/293t_filtered_feature_bc_matrix.tar.gz",
            ]
        }
    }
    
    # Alternative PBMC datasets for Harmony integration
    pbmc_datasets = {
        "pbmc_5k": {
            "name": "PBMC 5K",
            "description": "5k Peripheral Blood Mononuclear Cells",
            "files": [
                "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc5k/pbmc5k_filtered_feature_bc_matrix.h5",
            ]
        },
        "pbmc_10k": {
            "name": "PBMC 10K", 
            "description": "10k Peripheral Blood Mononuclear Cells",
            "files": [
                "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5",
            ]
        }
    }
    
    # Try to download datasets
    success_count = 0
    total_count = 0
    
    for dataset_id, info in datasets.items():
        print(f"\n📁 {info['name']}")
        print(f"   {info['description']}")
        
        dataset_dir = data_dir / dataset_id
        dataset_dir.mkdir(exist_ok=True)
        
        for file_url in info['files']:
            total_count += 1
            filename = file_url.split('/')[-1]
            filepath = dataset_dir / filename
            
            if download_file(file_url, filepath):
                success_count += 1
                
                # Extract if it's an archive
                if filename.endswith(('.tar.gz', '.zip', '.gz')):
                    extract_archive(filepath, dataset_dir)
    
    # Try PBMC datasets as alternatives
    print(f"\n📁 Alternative PBMC Datasets (commonly used with Harmony)")
    for dataset_id, info in pbmc_datasets.items():
        print(f"\n📁 {info['name']}")
        print(f"   {info['description']}")
        
        dataset_dir = data_dir / dataset_id
        dataset_dir.mkdir(exist_ok=True)
        
        for file_url in info['files']:
            total_count += 1
            filename = file_url.split('/')[-1]
            filepath = dataset_dir / filename
            
            if download_file(file_url, filepath):
                success_count += 1
    
    print("\n" + "=" * 50)
    print(f"📊 Download Summary: {success_count}/{total_count} files downloaded successfully")
    
    if success_count > 0:
        print(f"\n✅ Datasets saved to: {data_dir.absolute()}")
        print("\n📖 Usage Example:")
        print("```python")
        print("import scanpy as sc")
        print("import harmony")
        print("")
        print("# Load dataset")
        print("adata = sc.read_10x_h5('harmony_datasets/jurkat_293t_50_50/jurkat:293t_50:50_filtered_feature_bc_matrix.h5')")
        print("adata.var_names_unique()")
        print("")
        print("# Run Harmony integration")
        print("import harmonypy as hm")
        print("ho = hm.run_harmony(adata.X, adata.obs, 'batch')")
        print("```")
    else:
        print("\n⚠️  No files downloaded. You may need to:")
        print("1. Check your internet connection")
        print("2. Visit https://support.10xgenomics.com/single-cell-gene-expression/datasets")
        print("3. Fill out the form to access download links")
        print("4. Use the provided URLs manually")

if __name__ == "__main__":
    main()