#!/usr/bin/env python3
"""
Analyze existing datasets to understand their structure for creating quick demo versions
"""

import scanpy as sc
import numpy as np
import pandas as pd

def analyze_dataset(filepath):
    """Analyze a single dataset"""
    print(f"\n📊 Analyzing {filepath}")
    print("-" * 50)
    
    adata = sc.read_h5ad(filepath)
    
    print(f"Shape: {adata.shape} (cells x genes)")
    print(f"Size: {adata.X.nbytes / (1024*1024):.1f} MB")
    
    if 'batch' in adata.obs:
        batch_counts = adata.obs['batch'].value_counts()
        print(f"Batches: {dict(batch_counts)}")
    
    if 'cell_type' in adata.obs:
        cell_type_counts = adata.obs['cell_type'].value_counts()
        print(f"Cell types: {dict(cell_type_counts)}")
    
    # Memory usage breakdown
    if hasattr(adata.X, 'data'):
        sparsity = 1 - (adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]))
        print(f"Sparsity: {sparsity:.1%}")
    
    print(f"Observation keys: {list(adata.obs.columns)}")
    print(f"Variable keys: {list(adata.var.columns)}")
    
    return adata

def main():
    print("🔍 Dataset Analysis for Quick Demo Creation")
    print("=" * 60)
    
    datasets = [
        "pure_293t_simulated.h5ad",
        "pure_jurkat_simulated.h5ad", 
        "mixture_simulated.h5ad",
        "jurkat_293t_mixture_simulated.h5ad"
    ]
    
    results = {}
    
    for dataset in datasets:
        try:
            adata = analyze_dataset(dataset)
            results[dataset] = adata
        except Exception as e:
            print(f"❌ Error analyzing {dataset}: {e}")
    
    print(f"\n🎯 Recommendations for Quick Demo Dataset:")
    print("-" * 50)
    
    # Find smallest datasets for demo
    smallest_key = min(results.keys(), key=lambda k: results[k].n_obs)
    largest_key = max(results.keys(), key=lambda k: results[k].n_obs)
    
    print(f"Smallest: {smallest_key} with {results[smallest_key].n_obs} cells")
    print(f"Largest: {largest_key} with {results[largest_key].n_obs} cells")
    
    # Recommend strategy
    print(f"\n💡 Quick Demo Strategy:")
    print(f"1. Subsample large datasets to ~500 cells each")
    print(f"2. Reduce genes to ~500-1000 most variable")
    print(f"3. Target total: 1000 cells, 500 genes")
    print(f"4. Expected integration time: 2-5 minutes")

if __name__ == "__main__":
    main()