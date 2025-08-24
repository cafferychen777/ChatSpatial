#!/usr/bin/env python3
"""
Quick test to verify Harmony dataset functionality
"""

import scanpy as sc
import pandas as pd
import numpy as np

print("🔬 Loading dataset...")
try:
    adata = sc.read_h5ad("jurkat_293t_mixture_simulated.h5ad")
    print(f"✅ Dataset loaded successfully")
    print(f"📊 Shape: {adata.shape}")
    print(f"🧬 Cell types: {adata.obs['cell_type'].value_counts()}")
    print(f"🔬 Batches: {adata.obs['batch'].value_counts()}")
    
    print(f"📝 Observation columns: {list(adata.obs.columns)}")
    print(f"🧮 Variable info: {adata.var.shape}")
    
except Exception as e:
    print(f"❌ Error loading dataset: {e}")

print("\n🧪 Testing harmonypy import...")
try:
    import harmonypy as hm
    print("✅ harmonypy imported successfully")
except ImportError as e:
    print(f"❌ harmonypy import failed: {e}")

print("\n🔍 Checking other data files...")
import os
for file in ['pure_293t_simulated.h5ad', 'pure_jurkat_simulated.h5ad', 'mixture_simulated.h5ad']:
    if os.path.exists(file):
        try:
            data/test = sc.read_h5ad(file)
            print(f"✅ {file}: {data/test.shape}")
        except Exception as e:
            print(f"❌ {file}: {e}")
    else:
        print(f"❓ {file}: not found")

print("🏁 Quick test completed!")