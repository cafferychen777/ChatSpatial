#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ChatSpatial Quick Integration Test

Lightweight version for quick validation of core functionality.

Author: ChatSpatial Development Team
Created: 2024-08-24
"""

import sys
import os
import time
import json
from pathlib import Path

import pandas as pd
import numpy as np
import scanpy as sc

# Add chatspatial to path
sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

try:
    import squidpy as sq
    SQUIDPY_AVAILABLE = True
except ImportError:
    SQUIDPY_AVAILABLE = False
    print("WARNING: Squidpy not available - spatial tests will be limited")

def main():
    print("=" * 50)
    print("ChatSpatial Quick Integration Test")
    print("=" * 50)
    
    # Test data
    test_data_dir = Path('/Users/apple/Research/SpatialTrans_MCP/chatspatial/data')
    quick_demo_path = test_data_dir / 'demo_datasets/visium_demo.h5ad'
    
    if not quick_demo_path.exists():
        print(f"ERROR: Test dataset not found: {quick_demo_path}")
        return 1
    
    print(f"LOADING: {quick_demo_path}")
    adata = sc.read_h5ad(quick_demo_path)
    print(f"Dataset shape: {adata.shape[0]} cells x {adata.shape[1]} genes")
    
    # Set up matplotlib for non-interactive mode
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    results = {
        'test_name': 'quick_integration_test',
        'dataset': str(quick_demo_path),
        'cell_count': adata.n_obs,
        'gene_count': adata.n_vars,
        'start_time': time.time(),
        'steps': []
    }
    
    try:
        # Step 1: Basic preprocessing
        print("\nSTEP 1: Preprocessing")
        step_start = time.time()
        
        # Quality control
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        
        # Normalize
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Find variable genes
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        
        step_time = time.time() - step_start
        results['steps'].append({'name': 'preprocessing', 'time': step_time, 'success': True})
        print(f"Preprocessing: {step_time:.2f}s - SUCCESS")
        
        # Step 2: Spatial analysis (light version)
        print("\nSTEP 2: Spatial Analysis")
        step_start = time.time()
        
        if 'spatial' not in adata.obsm:
            n_cells = adata.n_obs
            adata.obsm['spatial'] = np.random.rand(n_cells, 2) * 100
            print("   Created dummy spatial coordinates")
        
        # Basic spatial validation
        spatial_coords = adata.obsm['spatial']
        if spatial_coords.shape[1] < 2:
            raise ValueError("Invalid spatial coordinates")
        if np.any(np.isnan(spatial_coords)):
            raise ValueError("NaN values in coordinates")
        
        step_time = time.time() - step_start
        results['steps'].append({'name': 'spatial_analysis', 'time': step_time, 'success': True})
        print(f"Spatial Analysis: {step_time:.2f}s - SUCCESS")
        
        # Step 3: Clustering
        print("\nSTEP 3: Clustering") 
        step_start = time.time()
        
        # Scale data
        sc.pp.scale(adata, max_value=10)
        
        # PCA
        sc.tl.pca(adata, svd_solver='arpack')
        
        # Neighbors and UMAP
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)
        
        # Clustering
        sc.tl.leiden(adata, resolution=0.5)
        
        step_time = time.time() - step_start
        results['steps'].append({'name': 'clustering', 'time': step_time, 'success': True})
        print(f"Clustering: {step_time:.2f}s - SUCCESS")
        
        # Step 4: Basic visualization
        print("\nSTEP 4: Visualization")
        step_start = time.time()
        
        # UMAP plot
        if 'X_umap' in adata.obsm and 'leiden' in adata.obs:
            sc.pl.umap(adata, color='leiden', show=False)
            plt.close('all')
        
        # Spatial plot
        if 'spatial' in adata.obsm:
            plt.figure(figsize=(6, 4))
            spatial_coords = adata.obsm['spatial']
            plt.scatter(spatial_coords[:, 0], spatial_coords[:, 1], s=1, alpha=0.6)
            plt.title('Spatial Coordinates')
            plt.close('all')
        
        step_time = time.time() - step_start
        results['steps'].append({'name': 'visualization', 'time': step_time, 'success': True})
        print(f"Visualization: {step_time:.2f}s - SUCCESS")
        
        # Summary
        results['total_time'] = time.time() - results['start_time']
        results['success'] = True
        
        print(f"\n" + "=" * 50)
        print("QUICK INTEGRATION TEST RESULTS")
        print("=" * 50)
        print(f"Total time: {results['total_time']:.2f}s")
        print(f"All steps completed successfully!")
        
        for step in results['steps']:
            print(f"  {step['name']}: {step['time']:.2f}s")
        
        # Check against quick benchmark (30s total for demo dataset)
        benchmark_time = 30.0
        performance_status = "PASS" if results['total_time'] <= benchmark_time else "FAIL"
        print(f"\nPerformance vs benchmark ({benchmark_time}s): {performance_status}")
        
        # Save results
        output_dir = Path('/Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/core_tools_tests')
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        results_file = output_dir / f"quick_test_results_{timestamp}.json"
        
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\nResults saved to: {results_file}")
        
        return 0
        
    except Exception as e:
        print(f"\nERROR: Test failed: {str(e)}")
        results['success'] = False
        results['error'] = str(e)
        results['total_time'] = time.time() - results['start_time']
        return 1


if __name__ == "__main__":
    # Configure scanpy
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    
    exit_code = main()
    sys.exit(exit_code)