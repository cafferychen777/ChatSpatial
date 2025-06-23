#!/usr/bin/env python
"""
Quick test script to verify new spatial transcriptomics methods.
"""

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def create_test_data():
    """Create simple test data."""
    print("Creating test data...")
    
    # Spatial data
    n_spots = 100
    n_genes = 200
    X = np.random.poisson(5, size=(n_spots, n_genes))
    coords = np.random.rand(n_spots, 2) * 50
    
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"spot_{i}" for i in range(n_spots)]),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    )
    adata.obsm['spatial'] = coords
    
    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    return adata


def test_paste():
    """Test PASTE registration."""
    print("\n=== Testing PASTE ===")
    try:
        from chatspatial.tools.spatial_registration import register_spatial_slices
        
        # Create two slices
        adata1 = create_test_data()
        adata2 = create_test_data()
        adata2.obsm['spatial'] = adata2.obsm['spatial'] + 10  # Shift
        
        registered = register_spatial_slices([adata1, adata2], method='paste')
        print("✓ PASTE registration successful")
        print(f"  - Registered {len(registered)} slices")
        return True
        
    except Exception as e:
        print(f"✗ PASTE failed: {e}")
        return False


def test_spatialDE():
    """Test SpatialDE."""
    print("\n=== Testing SpatialDE ===")
    try:
        from chatspatial.tools.spatial_statistics import find_spatial_variable_genes
        
        adata = create_test_data()
        results = find_spatial_variable_genes(adata, method='spatialDE', n_genes=10)
        print("✓ SpatialDE successful")
        print(f"  - Found {len(results)} spatial genes")
        print(f"  - Top gene: {results.iloc[0]['g']} (p={results.iloc[0]['pval']:.3e})")
        return True
        
    except Exception as e:
        print(f"✗ SpatialDE failed: {e}")
        return False


def test_stagate():
    """Test STAGATE."""
    print("\n=== Testing STAGATE ===")
    try:
        import asyncio
        from chatspatial.tools.spatial_domains import identify_spatial_domains
        from chatspatial.models.data import SpatialDomainParameters
        
        adata = create_test_data()
        sc.pp.highly_variable_genes(adata, n_top_genes=50)
        
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(
            method="stagate",
            n_domains=3,
            stagate_epochs=50  # Quick test
        )
        
        async def run():
            return await identify_spatial_domains("test", data_store, params)
        
        result = asyncio.run(run())
        print("✓ STAGATE successful")
        print(f"  - Found {result.n_domains} domains")
        return True
        
    except Exception as e:
        print(f"✗ STAGATE failed: {e}")
        return False


def test_banksy():
    """Test BANKSY."""
    print("\n=== Testing BANKSY ===")
    try:
        import asyncio
        from chatspatial.tools.spatial_domains import identify_spatial_domains
        from chatspatial.models.data import SpatialDomainParameters
        
        adata = create_test_data()
        data_store = {"test": {"adata": adata}}
        params = SpatialDomainParameters(
            method="banksy",
            n_domains=3,
            banksy_lambda=0.2
        )
        
        async def run():
            return await identify_spatial_domains("test", data_store, params)
        
        result = asyncio.run(run())
        print("✓ BANKSY successful")
        print(f"  - Found {result.n_domains} domains")
        return True
        
    except Exception as e:
        print(f"✗ BANKSY failed: {e}")
        return False


def test_spotlight():
    """Test SPOTlight."""
    print("\n=== Testing SPOTlight ===")
    try:
        from chatspatial.tools.deconvolution import is_spotlight_available, deconvolve_spotlight
        
        # Check availability first
        available, msg = is_spotlight_available()
        if not available:
            print(f"✗ SPOTlight not available: {msg}")
            return False
        
        # Create test data
        spatial = create_test_data()
        
        # Create reference
        n_cells = 200
        n_genes = 200
        ref_X = np.random.poisson(5, size=(n_cells, n_genes))
        reference = ad.AnnData(
            X=ref_X,
            obs=pd.DataFrame(
                index=[f"cell_{i}" for i in range(n_cells)],
                data={'cell_type': ['TypeA'] * 100 + ['TypeB'] * 100}
            ),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
        )
        
        proportions, stats = deconvolve_spotlight(spatial, reference)
        print("✓ SPOTlight successful")
        print(f"  - Deconvolved {proportions.shape[1]} cell types")
        return True
        
    except Exception as e:
        print(f"✗ SPOTlight failed: {e}")
        return False


def main():
    """Run all tests."""
    print("Testing newly added spatial transcriptomics methods")
    print("=" * 50)
    
    results = {
        'PASTE': test_paste(),
        'SpatialDE': test_spatialDE(),
        'STAGATE': test_stagate(),
        'BANKSY': test_banksy(),
        'SPOTlight': test_spotlight()
    }
    
    print("\n" + "=" * 50)
    print("Summary:")
    for method, success in results.items():
        status = "✓" if success else "✗"
        print(f"  {status} {method}")
    
    success_count = sum(results.values())
    total_count = len(results)
    print(f"\nPassed: {success_count}/{total_count}")
    
    if success_count == 0:
        print("\n⚠️  No methods succeeded. Please install required packages:")
        print("  - pip install paste-bio  # For PASTE")
        print("  - pip install spatialde  # For SpatialDE")
        print("  - pip install STAGATE    # For STAGATE")
        print("  - pip install rpy2       # For SPOTlight")
        print("  - R packages: BiocManager::install('SPOTlight')")


if __name__ == "__main__":
    main()