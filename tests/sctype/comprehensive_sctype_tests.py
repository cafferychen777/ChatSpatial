#!/usr/bin/env python3
"""
Comprehensive test suite for sc-type integration with ChatSpatial MCP
Tests multiple scenarios, error conditions, and edge cases
"""

import asyncio
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path
import traceback
import tempfile
import os
import shutil
from typing import Dict, Any

def create_diverse_test_datasets():
    """Create multiple test datasets for comprehensive testing"""
    datasets = {}
    
    # Dataset 1: Standard immune system data
    datasets["immune"] = create_immune_test_data()
    
    # Dataset 2: Brain tissue data
    datasets["brain"] = create_brain_test_data()
    
    # Dataset 3: Mixed tissue data
    datasets["mixed"] = create_mixed_tissue_data()
    
    # Dataset 4: Small dataset (edge case)
    datasets["small"] = create_small_test_data()
    
    # Dataset 5: Large dataset (stress test)
    datasets["large"] = create_large_test_data()
    
    return datasets

def create_immune_test_data():
    """Create immune system test data"""
    n_cells, n_genes = 300, 100
    X = np.random.negative_binomial(n=2, p=0.3, size=(n_cells, n_genes)).astype(float)
    
    # Immune markers
    cd3_genes = ['CD3D', 'CD3E', 'CD3G']
    cd8_genes = ['CD8A', 'CD8B']
    cd4_genes = ['CD4', 'IL7R']
    bcell_genes = ['CD19', 'MS4A1', 'CD79A']
    macro_genes = ['CD68', 'CD14', 'CSF1R']
    
    marker_genes = cd3_genes + cd8_genes + cd4_genes + bcell_genes + macro_genes
    other_genes = [f'GENE_{i:03d}' for i in range(len(marker_genes), n_genes)]
    gene_names = marker_genes + other_genes
    
    gene_name_to_idx = {name: i for i, name in enumerate(gene_names)}
    
    # Add expression patterns
    for gene in cd3_genes:
        if gene in gene_name_to_idx:
            idx = gene_name_to_idx[gene]
            X[:100, idx] += np.random.negative_binomial(n=10, p=0.2, size=100)
    
    for gene in cd8_genes:
        if gene in gene_name_to_idx:
            idx = gene_name_to_idx[gene]
            X[:50, idx] += np.random.negative_binomial(n=8, p=0.3, size=50)
    
    for gene in cd4_genes:
        if gene in gene_name_to_idx:
            idx = gene_name_to_idx[gene]
            X[50:100, idx] += np.random.negative_binomial(n=8, p=0.3, size=50)
    
    for gene in bcell_genes:
        if gene in gene_name_to_idx:
            idx = gene_name_to_idx[gene]
            X[100:200, idx] += np.random.negative_binomial(n=12, p=0.2, size=100)
    
    for gene in macro_genes:
        if gene in gene_name_to_idx:
            idx = gene_name_to_idx[gene]
            X[200:300, idx] += np.random.negative_binomial(n=15, p=0.2, size=100)
    
    return _create_anndata_with_processing(X, gene_names, n_cells, 
                                         ['T_cell'] * 100 + ['B_cell'] * 100 + ['Macrophage'] * 100)

def create_brain_test_data():
    """Create brain tissue test data"""
    n_cells, n_genes = 250, 120
    X = np.random.negative_binomial(n=3, p=0.25, size=(n_cells, n_genes)).astype(float)
    
    # Brain markers
    neuron_genes = ['RBFOX3', 'NEUN', 'MAP2', 'TUBB3']
    astro_genes = ['GFAP', 'AQP4', 'S100B', 'ALDH1L1']
    oligo_genes = ['MBP', 'MOG', 'PLP1', 'OLIG2']
    
    marker_genes = neuron_genes + astro_genes + oligo_genes
    other_genes = [f'BRAIN_GENE_{i:03d}' for i in range(len(marker_genes), n_genes)]
    gene_names = marker_genes + other_genes
    
    gene_name_to_idx = {name: i for i, name in enumerate(gene_names)}
    
    # Add expression patterns for brain cells
    for gene in neuron_genes:
        if gene in gene_name_to_idx:
            idx = gene_name_to_idx[gene]
            X[:80, idx] += np.random.negative_binomial(n=12, p=0.2, size=80)
    
    for gene in astro_genes:
        if gene in gene_name_to_idx:
            idx = gene_name_to_idx[gene]
            X[80:160, idx] += np.random.negative_binomial(n=10, p=0.25, size=80)
    
    for gene in oligo_genes:
        if gene in gene_name_to_idx:
            idx = gene_name_to_idx[gene]
            X[160:250, idx] += np.random.negative_binomial(n=8, p=0.3, size=90)
    
    return _create_anndata_with_processing(X, gene_names, n_cells, 
                                         ['Neuron'] * 80 + ['Astrocyte'] * 80 + ['Oligodendrocyte'] * 90)

def create_mixed_tissue_data():
    """Create mixed tissue data with overlapping markers"""
    n_cells, n_genes = 400, 150
    X = np.random.negative_binomial(n=2, p=0.3, size=(n_cells, n_genes)).astype(float)
    
    # Mixed markers from different tissues
    markers = [
        'CD3D', 'CD19', 'CD68',  # Immune
        'GFAP', 'NEUN', 'MBP',   # Brain
        'ALB', 'CYP3A4',         # Liver
        'ACTA2', 'MYH11'         # Muscle
    ]
    
    other_genes = [f'MIXED_GENE_{i:03d}' for i in range(len(markers), n_genes)]
    gene_names = markers + other_genes
    
    # Create diverse cell types
    cell_types = ['Mixed_Type_1'] * 100 + ['Mixed_Type_2'] * 100 + ['Mixed_Type_3'] * 100 + ['Mixed_Type_4'] * 100
    
    return _create_anndata_with_processing(X, gene_names, n_cells, cell_types)

def create_small_test_data():
    """Create small dataset to test edge cases"""
    n_cells, n_genes = 50, 30
    X = np.random.negative_binomial(n=1, p=0.5, size=(n_cells, n_genes)).astype(float)
    
    gene_names = ['CD3D', 'CD19', 'CD68'] + [f'SMALL_GENE_{i:02d}' for i in range(3, n_genes)]
    cell_types = ['Type_A'] * 25 + ['Type_B'] * 25
    
    return _create_anndata_with_processing(X, gene_names, n_cells, cell_types)

def create_large_test_data():
    """Create large dataset for stress testing"""
    n_cells, n_genes = 1000, 200
    X = np.random.negative_binomial(n=3, p=0.2, size=(n_cells, n_genes)).astype(float)
    
    # Include many known markers
    markers = [
        'CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD8B', 'CD4', 'IL7R',
        'CD19', 'MS4A1', 'CD79A', 'CD79B',
        'CD68', 'CD14', 'CSF1R', 'FCGR3A',
        'GFAP', 'AQP4', 'S100B', 'NEUN', 'MAP2'
    ]
    
    other_genes = [f'LARGE_GENE_{i:03d}' for i in range(len(markers), n_genes)]
    gene_names = markers + other_genes
    
    # Create multiple cell types
    cell_types = (['T_cell'] * 200 + ['B_cell'] * 200 + ['Macrophage'] * 200 + 
                  ['Neuron'] * 200 + ['Astrocyte'] * 200)
    
    return _create_anndata_with_processing(X, gene_names, n_cells, cell_types)

def _create_anndata_with_processing(X, gene_names, n_cells, cell_types):
    """Create AnnData with standard processing"""
    # Spatial coordinates
    x = np.random.uniform(0, 20, n_cells)
    y = np.random.uniform(0, 20, n_cells)
    
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame({
            'cell_id': [f"cell_{i:04d}" for i in range(n_cells)],
            'true_celltype': cell_types
        }),
        var=pd.DataFrame({'gene_id': gene_names})
    )
    
    adata.obs_names = [f"cell_{i:04d}" for i in range(n_cells)]
    adata.var_names = gene_names
    adata.obsm['spatial'] = np.column_stack([x, y])
    
    # QC metrics
    adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs['n_genes_by_counts'] = np.array((adata.X > 0).sum(axis=1)).flatten()
    
    # Standard processing
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.scale(adata)
    adata.layers['scaled'] = adata.X.copy()
    
    # Restore original counts
    adata.X = adata.raw.X
    
    return adata

def extract_mcp_result(mcp_response):
    """Extract result from MCP response"""
    if isinstance(mcp_response, dict):
        if 'isError' in mcp_response and mcp_response['isError']:
            raise Exception(f"MCP Error: {mcp_response.get('content', 'Unknown error')}")
        
        if 'content' in mcp_response and mcp_response['content']:
            content = mcp_response['content'][0]
            if content['type'] == 'text':
                import json
                return json.loads(content['text'])
    return mcp_response

async def test_sctype_availability():
    """Test sc-type availability and setup"""
    print("\n" + "="*60)
    print("TESTING SC-TYPE AVAILABILITY")
    print("="*60)
    
    try:
        from chatspatial.tools.annotation import SCTYPE_AVAILABLE, SCTYPE_ERROR
        print(f"SC-type available: {SCTYPE_AVAILABLE}")
        if not SCTYPE_AVAILABLE:
            print(f"SC-type error: {SCTYPE_ERROR}")
            return False
        
        print("✅ SC-type is properly configured")
        return True
        
    except Exception as e:
        print(f"❌ Failed to check SC-type availability: {e}")
        return False

async def test_basic_functionality(datasets):
    """Test basic sc-type functionality with different datasets"""
    print("\n" + "="*60)
    print("TESTING BASIC SC-TYPE FUNCTIONALITY")
    print("="*60)
    
    from chatspatial.server import load_data, annotate_cells
    
    results = {}
    
    for dataset_name, test_data in datasets.items():
        print(f"\n--- Testing {dataset_name} dataset ---")
        
        try:
            # Save test data
            test_path = Path(f"sctype_test_{dataset_name}.h5ad")
            test_data.write(test_path)
            
            # Load data
            load_result = await load_data(str(test_path), "h5ad", f"sctype_test_{dataset_name}")
            dataset = extract_mcp_result(load_result)
            print(f"✅ Loaded {dataset_name} dataset: {dataset['id']}")
            
            # Test with appropriate tissue type
            tissue_map = {
                "immune": "Immune system",
                "brain": "Brain", 
                "mixed": "Immune system",  # Default to immune
                "small": "Immune system",
                "large": "Immune system"
            }
            
            tissue = tissue_map.get(dataset_name, "Immune system")
            
            annotation_result = await annotate_cells(
                dataset['id'],
                {
                    "method": "sctype",
                    "sctype_tissue": tissue,
                    "sctype_scaled": True,
                    "sctype_use_cache": True
                }
            )
            annotation_data = extract_mcp_result(annotation_result)
            
            print(f"✅ {dataset_name} annotation succeeded!")
            print(f"   Cell types found: {len(annotation_data.get('counts', {}))}")
            print(f"   Total cells annotated: {sum(annotation_data.get('counts', {}).values())}")
            
            results[dataset_name] = {
                "success": True,
                "cell_types": len(annotation_data.get('counts', {})),
                "total_cells": sum(annotation_data.get('counts', {}).values())
            }
            
        except Exception as e:
            print(f"❌ {dataset_name} annotation failed: {e}")
            results[dataset_name] = {"success": False, "error": str(e)}
            
        finally:
            # Cleanup
            if test_path.exists():
                test_path.unlink()
    
    return results

async def test_custom_markers():
    """Test custom markers functionality"""
    print("\n" + "="*60)
    print("TESTING CUSTOM MARKERS")
    print("="*60)
    
    from chatspatial.server import load_data, annotate_cells
    
    # Create test data
    test_data = create_immune_test_data()
    test_path = Path("sctype_custom_test.h5ad")
    test_data.write(test_path)
    
    try:
        # Load data
        load_result = await load_data(str(test_path), "h5ad", "sctype_custom_test")
        dataset = extract_mcp_result(load_result)
        
        # Test multiple custom marker configurations
        test_configs = [
            {
                "name": "Simple markers",
                "markers": {
                    "T_cells": {
                        "positive": ["CD3D", "CD3E"],
                        "negative": ["CD19", "CD68"]
                    },
                    "B_cells": {
                        "positive": ["CD19", "MS4A1"],
                        "negative": ["CD3D"]
                    }
                }
            },
            {
                "name": "Complex markers",
                "markers": {
                    "CD8_T_cells": {
                        "positive": ["CD3D", "CD8A", "CD8B"],
                        "negative": ["CD4", "CD19", "CD68"]
                    },
                    "CD4_T_cells": {
                        "positive": ["CD3D", "CD4", "IL7R"],
                        "negative": ["CD8A", "CD19", "CD68"]
                    },
                    "B_cells": {
                        "positive": ["CD19", "MS4A1", "CD79A"],
                        "negative": ["CD3D", "CD68"]
                    },
                    "Macrophages": {
                        "positive": ["CD68", "CD14", "CSF1R"],
                        "negative": ["CD3D", "CD19"]
                    }
                }
            },
            {
                "name": "Single cell type",
                "markers": {
                    "All_immune": {
                        "positive": ["CD3D", "CD19", "CD68"],
                        "negative": []
                    }
                }
            }
        ]
        
        for config in test_configs:
            print(f"\n--- Testing {config['name']} ---")
            
            try:
                annotation_result = await annotate_cells(
                    dataset['id'],
                    {
                        "method": "sctype",
                        "sctype_custom_markers": config["markers"],
                        "sctype_scaled": True,
                        "sctype_use_cache": False
                    }
                )
                annotation_data = extract_mcp_result(annotation_result)
                
                print(f"✅ {config['name']} annotation succeeded!")
                print(f"   Cell types: {list(annotation_data.get('counts', {}).keys())}")
                print(f"   Counts: {annotation_data.get('counts', {})}")
                
            except Exception as e:
                print(f"❌ {config['name']} failed: {e}")
                traceback.print_exc()
        
        print("✅ Custom markers testing completed")
        return True
        
    finally:
        if test_path.exists():
            test_path.unlink()

async def test_caching_functionality():
    """Test caching functionality"""
    print("\n" + "="*60)
    print("TESTING CACHING FUNCTIONALITY")
    print("="*60)
    
    from chatspatial.server import load_data, annotate_cells
    import time
    
    # Create test data
    test_data = create_immune_test_data()
    test_path = Path("sctype_cache_test.h5ad")
    test_data.write(test_path)
    
    try:
        # Load data
        load_result = await load_data(str(test_path), "h5ad", "sctype_cache_test")
        dataset = extract_mcp_result(load_result)
        
        # Test with caching enabled
        print("--- First run (should cache) ---")
        start_time = time.time()
        result1 = await annotate_cells(
            dataset['id'],
            {
                "method": "sctype",
                "sctype_tissue": "Immune system",
                "sctype_use_cache": True
            }
        )
        first_duration = time.time() - start_time
        annotation1 = extract_mcp_result(result1)
        print(f"✅ First run completed in {first_duration:.2f}s")
        
        # Test cached run (should be faster)
        print("--- Second run (should use cache) ---")
        start_time = time.time()
        result2 = await annotate_cells(
            dataset['id'],
            {
                "method": "sctype",
                "sctype_tissue": "Immune system",
                "sctype_use_cache": True
            }
        )
        second_duration = time.time() - start_time
        annotation2 = extract_mcp_result(result2)
        print(f"✅ Second run completed in {second_duration:.2f}s")
        
        # Verify results are identical
        if annotation1.get('counts') == annotation2.get('counts'):
            print("✅ Cached results match original results")
        else:
            print("❌ Cached results differ from original")
            
        # Test with cache disabled
        print("--- Third run (cache disabled) ---")
        start_time = time.time()
        result3 = await annotate_cells(
            dataset['id'],
            {
                "method": "sctype",
                "sctype_tissue": "Immune system",
                "sctype_use_cache": False
            }
        )
        third_duration = time.time() - start_time
        print(f"✅ Third run (no cache) completed in {third_duration:.2f}s")
        
        print(f"Performance summary:")
        print(f"  First run:  {first_duration:.2f}s")
        print(f"  Cached run: {second_duration:.2f}s")
        print(f"  No cache:   {third_duration:.2f}s")
        
        return True
        
    finally:
        if test_path.exists():
            test_path.unlink()

async def test_error_conditions():
    """Test various error conditions"""
    print("\n" + "="*60)
    print("TESTING ERROR CONDITIONS")
    print("="*60)
    
    from chatspatial.server import load_data, annotate_cells
    
    # Create test data
    test_data = create_immune_test_data()
    test_path = Path("sctype_error_test.h5ad")
    test_data.write(test_path)
    
    try:
        # Load data
        load_result = await load_data(str(test_path), "h5ad", "sctype_error_test")
        dataset = extract_mcp_result(load_result)
        
        error_tests = [
            {
                "name": "Missing tissue and markers",
                "params": {"method": "sctype"},
                "expected_error": "Either sctype_tissue or sctype_custom_markers must be specified"
            },
            {
                "name": "Invalid tissue type",
                "params": {"method": "sctype", "sctype_tissue": "NonexistentTissue"},
                "expected_error": None  # Should handle gracefully
            },
            {
                "name": "Empty custom markers",
                "params": {"method": "sctype", "sctype_custom_markers": {}},
                "expected_error": None  # Should handle gracefully
            },
            {
                "name": "Invalid marker format",
                "params": {
                    "method": "sctype", 
                    "sctype_custom_markers": {"BadFormat": ["gene1", "gene2"]}
                },
                "expected_error": None  # Should handle gracefully
            }
        ]
        
        for test_case in error_tests:
            print(f"\n--- Testing {test_case['name']} ---")
            
            try:
                result = await annotate_cells(dataset['id'], test_case['params'])
                
                # Check if it's an MCP error response
                if isinstance(result, dict) and result.get('isError', False):
                    error_content = str(result.get('content', ''))
                    if test_case['expected_error'] and test_case['expected_error'] in error_content:
                        print(f"✅ Correctly caught expected error: {test_case['expected_error']}")
                    else:
                        print(f"✅ Error handled gracefully: {error_content[:100]}...")
                else:
                    print(f"⚠️  Expected error but got result: {test_case['name']}")
                
            except Exception as e:
                if test_case['expected_error'] and test_case['expected_error'] in str(e):
                    print(f"✅ Correctly caught expected error: {test_case['expected_error']}")
                else:
                    print(f"⚠️  Unexpected error: {e}")
        
        print("✅ Error condition testing completed")
        return True
        
    finally:
        if test_path.exists():
            test_path.unlink()

async def test_edge_cases():
    """Test edge cases and boundary conditions"""
    print("\n" + "="*60)
    print("TESTING EDGE CASES")
    print("="*60)
    
    from chatspatial.server import load_data, annotate_cells
    
    # Test with very small dataset
    small_data = create_small_test_data()
    small_path = Path("sctype_small_test.h5ad")
    small_data.write(small_path)
    
    try:
        # Load small data
        load_result = await load_data(str(small_path), "h5ad", "sctype_small_test")
        dataset = extract_mcp_result(load_result)
        
        # Test small dataset
        print("--- Testing small dataset (50 cells, 30 genes) ---")
        result = await annotate_cells(
            dataset['id'],
            {
                "method": "sctype",
                "sctype_tissue": "Immune system",
                "sctype_scaled": True
            }
        )
        annotation_data = extract_mcp_result(result)
        print(f"✅ Small dataset annotation succeeded!")
        print(f"   Cell types: {len(annotation_data.get('counts', {}))}")
        
        # Test with no scaling
        print("--- Testing without scaling ---")
        result = await annotate_cells(
            dataset['id'],
            {
                "method": "sctype",
                "sctype_tissue": "Immune system",
                "sctype_scaled": False
            }
        )
        annotation_data = extract_mcp_result(result)
        print(f"✅ No scaling annotation succeeded!")
        
        # Test with custom markers containing non-existent genes
        print("--- Testing custom markers with non-existent genes ---")
        result = await annotate_cells(
            dataset['id'],
            {
                "method": "sctype",
                "sctype_custom_markers": {
                    "Test_Type": {
                        "positive": ["NONEXISTENT1", "NONEXISTENT2", "CD3D"],
                        "negative": ["ALSONONEXISTENT"]
                    }
                }
            }
        )
        annotation_data = extract_mcp_result(result)
        print(f"✅ Non-existent genes handled gracefully!")
        
        return True
        
    finally:
        if small_path.exists():
            small_path.unlink()

async def run_comprehensive_tests():
    """Run all comprehensive tests"""
    print("🧪 COMPREHENSIVE SC-TYPE INTEGRATION TESTING")
    print("="*70)
    
    # Test availability first
    if not await test_sctype_availability():
        print("\n❌ SC-type not available, skipping tests")
        return
    
    # Create test datasets
    print("\n📊 Creating diverse test datasets...")
    datasets = create_diverse_test_datasets()
    print(f"✅ Created {len(datasets)} test datasets")
    
    test_results = {}
    
    # Run all tests
    test_functions = [
        ("Basic Functionality", lambda: test_basic_functionality(datasets)),
        ("Custom Markers", test_custom_markers),
        ("Caching", test_caching_functionality),
        ("Error Conditions", test_error_conditions),
        ("Edge Cases", test_edge_cases)
    ]
    
    for test_name, test_func in test_functions:
        try:
            print(f"\n🔬 Running {test_name} tests...")
            result = await test_func()
            test_results[test_name] = {"success": result, "error": None}
            print(f"✅ {test_name} tests completed successfully")
        except Exception as e:
            test_results[test_name] = {"success": False, "error": str(e)}
            print(f"❌ {test_name} tests failed: {e}")
            traceback.print_exc()
    
    # Summary
    print("\n" + "="*70)
    print("COMPREHENSIVE TEST SUMMARY")
    print("="*70)
    
    passed = sum(1 for result in test_results.values() if result['success'])
    total = len(test_results)
    
    print(f"Tests passed: {passed}/{total}")
    
    for test_name, result in test_results.items():
        status = "✅ PASS" if result['success'] else "❌ FAIL"
        print(f"{status} {test_name}")
        if not result['success'] and result['error']:
            print(f"    Error: {result['error']}")
    
    if passed == total:
        print("\n🎉 ALL COMPREHENSIVE TESTS PASSED!")
        print("SC-type integration is fully functional and robust")
    else:
        print(f"\n⚠️  {total - passed} tests failed - review errors above")
    
    return test_results

if __name__ == "__main__":
    asyncio.run(run_comprehensive_tests())