#!/usr/bin/env python
"""
Comprehensive tests for PASTE, SPOTlight, and SpatialDE.
Tests different scenarios, edge cases, and parameter variations.
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import warnings
import time
warnings.filterwarnings('ignore')

# Color codes for output
GREEN = '\033[92m'
RED = '\033[91m'
BLUE = '\033[94m'
YELLOW = '\033[93m'
RESET = '\033[0m'

def print_test_header(test_name):
    """Print formatted test header."""
    print(f"\n{BLUE}{'='*70}{RESET}")
    print(f"{BLUE}Testing: {test_name}{RESET}")
    print(f"{BLUE}{'='*70}{RESET}")

def print_result(success, message):
    """Print colored result message."""
    if success:
        print(f"{GREEN}✅ {message}{RESET}")
    else:
        print(f"{RED}❌ {message}{RESET}")

def test_paste_comprehensive():
    """Comprehensive tests for PASTE spatial registration."""
    print(f"\n{YELLOW}🔬 COMPREHENSIVE PASTE TESTS{RESET}")
    from chatspatial.tools.spatial_registration import register_spatial_slices
    
    test_results = []
    
    # Test 1: Basic two-slice alignment
    print_test_header("Test 1: Basic Two-Slice Alignment")
    try:
        n_spots = 100
        n_genes = 50
        np.random.seed(42)
        
        slice1 = ad.AnnData(X=np.random.poisson(5, size=(n_spots, n_genes)).astype(float))
        slice2 = ad.AnnData(X=np.random.poisson(5, size=(n_spots, n_genes)).astype(float))
        
        # Add gene names
        slice1.var_names = [f'gene_{i}' for i in range(n_genes)]
        slice2.var_names = [f'gene_{i}' for i in range(n_genes)]
        
        # Add spatial coordinates
        slice1.obsm['spatial'] = np.random.rand(n_spots, 2) * 10
        slice2.obsm['spatial'] = np.random.rand(n_spots, 2) * 10 + 1
        
        # Add shared spatial pattern
        for i in range(5):  # First 5 genes have spatial patterns
            pattern1 = np.sin(slice1.obsm['spatial'][:, 0] / 3 + i) * 10
            pattern2 = np.sin(slice2.obsm['spatial'][:, 0] / 3 + i) * 10
            slice1.X[:, i] += pattern1
            slice2.X[:, i] += pattern2
        
        aligned = register_spatial_slices([slice1, slice2], method='paste', alpha=0.1)
        
        success = all('spatial_registered' in s.obsm for s in aligned)
        print_result(success, f"Basic alignment: {len(aligned)} slices aligned")
        test_results.append(("Basic alignment", success))
        
    except Exception as e:
        print_result(False, f"Basic alignment failed: {e}")
        test_results.append(("Basic alignment", False))
    
    # Test 2: Multiple slice alignment (>2 slices)
    print_test_header("Test 2: Multiple Slice Alignment (4 slices)")
    try:
        n_spots = 80
        n_genes = 40
        n_slices = 4
        slices = []
        
        for i in range(n_slices):
            slice_data = ad.AnnData(X=np.random.poisson(5, size=(n_spots, n_genes)).astype(float))
            slice_data.var_names = [f'gene_{j}' for j in range(n_genes)]
            slice_data.obsm['spatial'] = np.random.rand(n_spots, 2) * 10 + i * 0.5
            
            # Add spatial gradients
            gradient = slice_data.obsm['spatial'][:, 0] / 10
            slice_data.X[:, 0] += gradient * 15
            slices.append(slice_data)
        
        aligned = register_spatial_slices(slices, method='paste', alpha=0.1)
        
        success = len(aligned) == n_slices and all('spatial_registered' in s.obsm for s in aligned)
        print_result(success, f"Multi-slice alignment: {len(aligned)} slices aligned")
        test_results.append(("Multi-slice alignment", success))
        
    except Exception as e:
        print_result(False, f"Multi-slice alignment failed: {e}")
        test_results.append(("Multi-slice alignment", False))
    
    # Test 3: Different alpha values
    print_test_header("Test 3: Different Alpha Parameters")
    alphas = [0.01, 0.1, 0.5, 1.0]
    for alpha in alphas:
        try:
            aligned = register_spatial_slices([slice1, slice2], method='paste', alpha=alpha)
            success = all('spatial_registered' in s.obsm for s in aligned)
            print_result(success, f"Alpha={alpha}: Registration successful")
            test_results.append((f"Alpha={alpha}", success))
        except Exception as e:
            print_result(False, f"Alpha={alpha} failed: {e}")
            test_results.append((f"Alpha={alpha}", False))
    
    # Test 4: Slices with different numbers of spots
    print_test_header("Test 4: Different Number of Spots")
    try:
        slice1_diff = ad.AnnData(X=np.random.poisson(5, size=(100, 50)).astype(float))
        slice2_diff = ad.AnnData(X=np.random.poisson(5, size=(120, 50)).astype(float))
        
        slice1_diff.var_names = [f'gene_{i}' for i in range(50)]
        slice2_diff.var_names = [f'gene_{i}' for i in range(50)]
        
        slice1_diff.obsm['spatial'] = np.random.rand(100, 2) * 10
        slice2_diff.obsm['spatial'] = np.random.rand(120, 2) * 10
        
        aligned = register_spatial_slices([slice1_diff, slice2_diff], method='paste')
        success = all('spatial_registered' in s.obsm for s in aligned)
        print_result(success, f"Different spot counts: {[s.shape[0] for s in aligned]} spots")
        test_results.append(("Different spot counts", success))
        
    except Exception as e:
        print_result(False, f"Different spot counts failed: {e}")
        test_results.append(("Different spot counts", False))
    
    # Test 5: Performance with larger datasets
    print_test_header("Test 5: Performance Test (Large Dataset)")
    try:
        start_time = time.time()
        n_spots_large = 500
        n_genes_large = 200
        
        slice1_large = ad.AnnData(X=np.random.poisson(5, size=(n_spots_large, n_genes_large)).astype(float))
        slice2_large = ad.AnnData(X=np.random.poisson(5, size=(n_spots_large, n_genes_large)).astype(float))
        
        slice1_large.var_names = [f'gene_{i}' for i in range(n_genes_large)]
        slice2_large.var_names = [f'gene_{i}' for i in range(n_genes_large)]
        
        slice1_large.obsm['spatial'] = np.random.rand(n_spots_large, 2) * 10
        slice2_large.obsm['spatial'] = np.random.rand(n_spots_large, 2) * 10
        
        aligned = register_spatial_slices([slice1_large, slice2_large], method='paste')
        elapsed = time.time() - start_time
        
        success = all('spatial_registered' in s.obsm for s in aligned)
        print_result(success, f"Large dataset ({n_spots_large} spots, {n_genes_large} genes): {elapsed:.2f}s")
        test_results.append(("Large dataset", success))
        
    except Exception as e:
        print_result(False, f"Large dataset failed: {e}")
        test_results.append(("Large dataset", False))
    
    return test_results


def test_spotlight_comprehensive():
    """Comprehensive tests for SPOTlight deconvolution."""
    print(f"\n{YELLOW}🔬 COMPREHENSIVE SPOTLIGHT TESTS{RESET}")
    from chatspatial.tools.deconvolution import deconvolve_spotlight
    
    test_results = []
    
    # Test 1: Basic deconvolution
    print_test_header("Test 1: Basic Deconvolution")
    try:
        n_spots = 50
        n_genes = 100
        n_cells = 200
        n_cell_types = 4
        np.random.seed(42)
        
        # Create reference
        sc_data = ad.AnnData(X=np.random.poisson(5, size=(n_cells, n_genes)))
        cell_types = np.random.choice([f'CellType_{i}' for i in range(n_cell_types)], n_cells)
        sc_data.obs['cell_type'] = cell_types
        
        # Add cell type signatures
        for i, ct in enumerate(sc_data.obs['cell_type'].unique()):
            mask = sc_data.obs['cell_type'] == ct
            sc_data.X[mask, i*20:(i+1)*20] += np.random.poisson(15, size=(mask.sum(), 20))
        
        # Create spatial data
        spatial_data = ad.AnnData(X=np.zeros((n_spots, n_genes)))
        true_props = np.random.dirichlet(np.ones(n_cell_types), n_spots)
        
        for i in range(n_spots):
            for j, ct in enumerate(sc_data.obs['cell_type'].unique()):
                ct_expr = sc_data[sc_data.obs['cell_type'] == ct].X.mean(axis=0)
                spatial_data.X[i] += true_props[i, j] * ct_expr
        
        spatial_data.X = np.round(spatial_data.X).astype(int)
        
        props, stats = deconvolve_spotlight(spatial_data, sc_data, cell_type_key='cell_type')
        
        success = props.shape == (n_spots, n_cell_types) and np.allclose(props.sum(axis=1), 1.0, atol=0.01)
        print_result(success, f"Basic deconvolution: {props.shape} proportions matrix")
        test_results.append(("Basic deconvolution", success))
        
    except Exception as e:
        print_result(False, f"Basic deconvolution failed: {e}")
        test_results.append(("Basic deconvolution", False))
    
    # Test 2: Different numbers of cell types
    print_test_header("Test 2: Various Numbers of Cell Types")
    for n_ct in [2, 5, 8, 10]:
        try:
            # Create reference with n_ct cell types
            sc_data_ct = ad.AnnData(X=np.random.poisson(5, size=(200, 100)))
            sc_data_ct.obs['cell_type'] = np.random.choice([f'CT_{i}' for i in range(n_ct)], 200)
            
            # Add signatures
            genes_per_ct = 100 // n_ct
            for i, ct in enumerate(sc_data_ct.obs['cell_type'].unique()):
                mask = sc_data_ct.obs['cell_type'] == ct
                start_gene = i * genes_per_ct
                end_gene = min((i + 1) * genes_per_ct, 100)
                sc_data_ct.X[mask, start_gene:end_gene] += 10
            
            # Create spatial data
            spatial_ct = ad.AnnData(X=np.random.poisson(10, size=(30, 100)))
            
            props, _ = deconvolve_spotlight(spatial_ct, sc_data_ct, cell_type_key='cell_type')
            success = props.shape[1] == n_ct
            print_result(success, f"{n_ct} cell types: Shape {props.shape}")
            test_results.append((f"{n_ct} cell types", success))
            
        except Exception as e:
            print_result(False, f"{n_ct} cell types failed: {e}")
            test_results.append((f"{n_ct} cell types", False))
    
    # Test 3: Different n_top_genes parameters
    print_test_header("Test 3: Different n_top_genes Parameters")
    for n_top in [500, 1000, 2000, 3000]:
        try:
            props, stats = deconvolve_spotlight(
                spatial_data, sc_data, 
                cell_type_key='cell_type',
                n_top_genes=n_top
            )
            n_used = stats.get('n_common_genes', 0)
            success = props.shape[0] == n_spots
            print_result(success, f"n_top_genes={n_top}: Used {n_used} genes")
            test_results.append((f"n_top_genes={n_top}", success))
            
        except Exception as e:
            print_result(False, f"n_top_genes={n_top} failed: {e}")
            test_results.append((f"n_top_genes={n_top}", False))
    
    # Test 4: Sparse data
    print_test_header("Test 4: Sparse Data (High Zero Content)")
    try:
        # Create very sparse data
        sparse_sc = ad.AnnData(X=np.random.poisson(0.5, size=(200, 100)))
        sparse_sc.obs['cell_type'] = np.random.choice(['Sparse_A', 'Sparse_B', 'Sparse_C'], 200)
        
        sparse_spatial = ad.AnnData(X=np.random.poisson(0.5, size=(30, 100)))
        
        props, _ = deconvolve_spotlight(sparse_spatial, sparse_sc, cell_type_key='cell_type')
        sparsity = (sparse_spatial.X == 0).sum() / sparse_spatial.X.size
        success = props.shape[0] == 30
        print_result(success, f"Sparse data ({sparsity:.1%} zeros): Completed")
        test_results.append(("Sparse data", success))
        
    except Exception as e:
        print_result(False, f"Sparse data failed: {e}")
        test_results.append(("Sparse data", False))
    
    # Test 5: Large dataset performance
    print_test_header("Test 5: Performance Test (Large Dataset)")
    try:
        start_time = time.time()
        
        # Large reference
        large_sc = ad.AnnData(X=np.random.poisson(5, size=(1000, 500)))
        large_sc.obs['cell_type'] = np.random.choice([f'CT_{i}' for i in range(8)], 1000)
        
        # Large spatial
        large_spatial = ad.AnnData(X=np.random.poisson(8, size=(100, 500)))
        
        props, _ = deconvolve_spotlight(large_spatial, large_sc, cell_type_key='cell_type')
        elapsed = time.time() - start_time
        
        success = props.shape == (100, 8)
        print_result(success, f"Large dataset (1000 cells, 500 genes): {elapsed:.2f}s")
        test_results.append(("Large dataset", success))
        
    except Exception as e:
        print_result(False, f"Large dataset failed: {e}")
        test_results.append(("Large dataset", False))
    
    return test_results


def test_spatialDE_comprehensive():
    """Comprehensive tests for SpatialDE."""
    print(f"\n{YELLOW}🔬 COMPREHENSIVE SPATIALDES TESTS{RESET}")
    from chatspatial.tools.spatial_statistics import find_spatial_variable_genes
    
    test_results = []
    
    # Test 1: Basic spatial pattern detection
    print_test_header("Test 1: Basic Spatial Pattern Detection")
    try:
        n_spots = 100
        n_genes = 30
        np.random.seed(42)
        
        coords = np.random.rand(n_spots, 2) * 10
        expression = np.random.poisson(5, size=(n_spots, n_genes))
        
        # Add different spatial patterns
        patterns = {
            'linear_x': lambda c: c[:, 0] / 10,
            'linear_y': lambda c: c[:, 1] / 10,
            'radial': lambda c: np.sqrt((c[:, 0] - 5)**2 + (c[:, 1] - 5)**2) / 10,
            'periodic': lambda c: np.sin(c[:, 0] * np.pi / 5),
            'bimodal': lambda c: (c[:, 0] > 5).astype(float)
        }
        
        for i, (pattern_name, pattern_func) in enumerate(patterns.items()):
            if i < n_genes:
                expression[:, i] = expression[:, i] + pattern_func(coords) * 20
        
        adata = ad.AnnData(X=expression.astype(float))
        adata.obsm['spatial'] = coords
        adata.var_names = [f'gene_{i}' for i in range(n_genes)]
        
        results = find_spatial_variable_genes(adata, method='spatialDE', n_genes=10, normalized=False)
        
        n_detected = sum(1 for g in results.head(5)['g'] if g in ['gene_0', 'gene_1', 'gene_2', 'gene_3', 'gene_4'])
        success = n_detected >= 3
        print_result(success, f"Pattern detection: {n_detected}/5 spatial genes in top 5")
        test_results.append(("Pattern detection", success))
        
    except Exception as e:
        print_result(False, f"Pattern detection failed: {e}")
        test_results.append(("Pattern detection", False))
    
    # Test 2: Different spatial scales
    print_test_header("Test 2: Different Spatial Scales")
    try:
        # Test with different coordinate scales
        scales = [1, 10, 100, 1000]
        for scale in scales:
            coords_scaled = coords * scale
            adata_scaled = ad.AnnData(X=expression.astype(float))
            adata_scaled.obsm['spatial'] = coords_scaled
            adata_scaled.var_names = [f'gene_{i}' for i in range(n_genes)]
            
            results = find_spatial_variable_genes(adata_scaled, method='spatialDE', n_genes=5, normalized=False)
            success = len(results) == 5
            print_result(success, f"Scale {scale}: Found {len(results)} genes")
            test_results.append((f"Scale {scale}", success))
            
    except Exception as e:
        print_result(False, f"Scale testing failed: {e}")
        test_results.append(("Scale testing", False))
    
    # Test 3: Pre-normalized vs raw data
    print_test_header("Test 3: Normalized vs Raw Data")
    try:
        # Test with pre-normalized data
        adata_norm = adata.copy()
        sc.pp.normalize_total(adata_norm, target_sum=1e4)
        sc.pp.log1p(adata_norm)
        
        results_norm = find_spatial_variable_genes(adata_norm, method='spatialDE', n_genes=10, normalized=True)
        results_raw = find_spatial_variable_genes(adata, method='spatialDE', n_genes=10, normalized=False)
        
        success = len(results_norm) == 10 and len(results_raw) == 10
        print_result(success, f"Normalized: {len(results_norm)} genes, Raw: {len(results_raw)} genes")
        test_results.append(("Normalized vs Raw", success))
        
    except Exception as e:
        print_result(False, f"Normalization test failed: {e}")
        test_results.append(("Normalized vs Raw", False))
    
    # Test 4: Different numbers of spots
    print_test_header("Test 4: Different Dataset Sizes")
    spot_counts = [30, 50, 100, 200, 500]
    for n_spots_test in spot_counts:
        try:
            coords_test = np.random.rand(n_spots_test, 2) * 10
            expr_test = np.random.poisson(5, size=(n_spots_test, 20))
            
            # Add spatial pattern
            expr_test[:, 0] += (coords_test[:, 0] / 10 * 20).astype(int)
            
            adata_test = ad.AnnData(X=expr_test.astype(float))
            adata_test.obsm['spatial'] = coords_test
            adata_test.var_names = [f'gene_{i}' for i in range(20)]
            
            results = find_spatial_variable_genes(adata_test, method='spatialDE', n_genes=5, normalized=False)
            success = 'gene_0' in results['g'].values
            print_result(success, f"{n_spots_test} spots: Spatial gene detected")
            test_results.append((f"{n_spots_test} spots", success))
            
        except Exception as e:
            print_result(False, f"{n_spots_test} spots failed: {e}")
            test_results.append((f"{n_spots_test} spots", False))
    
    # Test 5: Performance with many genes
    print_test_header("Test 5: Performance Test (Many Genes)")
    try:
        start_time = time.time()
        
        n_spots_perf = 100
        n_genes_perf = 1000
        
        coords_perf = np.random.rand(n_spots_perf, 2) * 10
        expr_perf = np.random.poisson(5, size=(n_spots_perf, n_genes_perf))
        
        # Add spatial patterns to first 50 genes
        for i in range(50):
            pattern_idx = i % len(patterns)
            pattern_func = list(patterns.values())[pattern_idx]
            expr_perf[:, i] += (pattern_func(coords_perf) * 20).astype(int)
        
        adata_perf = ad.AnnData(X=expr_perf.astype(float))
        adata_perf.obsm['spatial'] = coords_perf
        adata_perf.var_names = [f'gene_{i}' for i in range(n_genes_perf)]
        
        results = find_spatial_variable_genes(adata_perf, method='spatialDE', n_genes=100, normalized=False)
        elapsed = time.time() - start_time
        
        n_spatial_detected = sum(1 for g in results['g'] if int(g.split('_')[1]) < 50)
        success = n_spatial_detected > 30
        print_result(success, f"{n_genes_perf} genes: {n_spatial_detected}/100 spatial, {elapsed:.2f}s")
        test_results.append(("Many genes", success))
        
    except Exception as e:
        print_result(False, f"Performance test failed: {e}")
        test_results.append(("Many genes", False))
    
    return test_results


def summarize_results(all_results):
    """Print summary of all test results."""
    print(f"\n{YELLOW}{'='*70}{RESET}")
    print(f"{YELLOW}📊 COMPREHENSIVE TEST SUMMARY{RESET}")
    print(f"{YELLOW}{'='*70}{RESET}")
    
    for method, results in all_results.items():
        passed = sum(1 for _, success in results if success)
        total = len(results)
        success_rate = passed / total * 100 if total > 0 else 0
        
        if success_rate == 100:
            color = GREEN
        elif success_rate >= 80:
            color = YELLOW
        else:
            color = RED
        
        print(f"\n{BLUE}{method}:{RESET}")
        print(f"{color}  {passed}/{total} tests passed ({success_rate:.1f}%){RESET}")
        
        # Show failed tests
        failed_tests = [name for name, success in results if not success]
        if failed_tests:
            print(f"  {RED}Failed tests: {', '.join(failed_tests)}{RESET}")


def main():
    """Run all comprehensive tests."""
    print(f"{YELLOW}🧪 RUNNING COMPREHENSIVE TESTS FOR SPATIAL METHODS{RESET}")
    print(f"{YELLOW}This will test various scenarios and edge cases{RESET}")
    
    all_results = {}
    
    # Run PASTE tests
    paste_results = test_paste_comprehensive()
    all_results['PASTE'] = paste_results
    
    # Run SPOTlight tests
    spotlight_results = test_spotlight_comprehensive()
    all_results['SPOTlight'] = spotlight_results
    
    # Run SpatialDE tests
    spatialDE_results = test_spatialDE_comprehensive()
    all_results['SpatialDE'] = spatialDE_results
    
    # Summary
    summarize_results(all_results)
    
    # Overall summary
    total_passed = sum(sum(1 for _, s in results if s) for results in all_results.values())
    total_tests = sum(len(results) for results in all_results.values())
    overall_rate = total_passed / total_tests * 100 if total_tests > 0 else 0
    
    print(f"\n{YELLOW}{'='*70}{RESET}")
    if overall_rate == 100:
        print(f"{GREEN}🎉 ALL TESTS PASSED! ({total_passed}/{total_tests}){RESET}")
    elif overall_rate >= 80:
        print(f"{YELLOW}⚠️  MOSTLY PASSED: {total_passed}/{total_tests} ({overall_rate:.1f}%){RESET}")
    else:
        print(f"{RED}❌ NEEDS ATTENTION: {total_passed}/{total_tests} ({overall_rate:.1f}%){RESET}")


if __name__ == "__main__":
    main()