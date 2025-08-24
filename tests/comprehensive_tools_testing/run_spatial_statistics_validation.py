#!/usr/bin/env python3
"""
Advanced Spatial Statistics Validation Script

This script performs in-depth validation of spatial statistics algorithms,
including mathematical correctness verification, statistical properties testing,
and comprehensive performance benchmarking.
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import logging
from pathlib import Path
import warnings
import time
import json
from typing import Dict, List, Any, Tuple

# Add the project root to the path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from chatspatial.tools.spatial_statistics import SpatialStatistics

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

class SpatialStatisticsValidator:
    """Advanced validator for spatial statistics algorithms."""
    
    def __init__(self):
        self.spatial_stats = SpatialStatistics()
        self.validation_results = {}
        self.performance_metrics = {}
        
    def create_test_datasets(self) -> Dict[str, ad.AnnData]:
        """Create comprehensive test datasets with known spatial patterns."""
        datasets = {}
        
        # 1. Perfect spatial autocorrelation test
        datasets['perfect_autocorr'] = self._create_perfect_autocorr_data()
        
        # 2. Random (no autocorrelation) test
        datasets['random'] = self._create_random_data()
        
        # 3. Checkerboard (negative autocorrelation) test
        datasets['checkerboard'] = self._create_checkerboard_data()
        
        # 4. Hotspot pattern test
        datasets['hotspots'] = self._create_hotspot_data()
        
        # 5. Gradient pattern test
        datasets['gradient'] = self._create_gradient_data()
        
        # 6. Multi-scale pattern test
        datasets['multiscale'] = self._create_multiscale_data()
        
        return datasets
    
    def _create_perfect_autocorr_data(self, size: int = 100) -> ad.AnnData:
        """Create data with perfect positive spatial autocorrelation."""
        np.random.seed(123)
        
        # Create regular grid
        grid_dim = int(np.sqrt(size))
        x, y = np.meshgrid(np.arange(grid_dim), np.arange(grid_dim))
        coords = np.column_stack([x.ravel()[:size], y.ravel()[:size]])
        
        # Create perfectly clustered expression (smooth gradient)
        expression = np.zeros((size, 5))
        for i in range(5):
            # Smooth spatial pattern
            pattern = (coords[:, 0] + coords[:, 1] + i) / (2 * grid_dim)
            expression[:, i] = pattern * 50 + np.random.normal(0, 1, size)
        
        adata = ad.AnnData(
            X=expression,
            obs=pd.DataFrame(index=[f'spot_{i}' for i in range(size)]),
            var=pd.DataFrame(index=[f'GENE_AUTOCORR_{i}' for i in range(5)])
        )
        adata.obsm['spatial'] = coords.astype(float)
        return adata
    
    def _create_random_data(self, size: int = 100) -> ad.AnnData:
        """Create data with no spatial autocorrelation (random)."""
        np.random.seed(456)
        
        coords = np.random.uniform(0, 10, size=(size, 2))
        expression = np.random.poisson(10, size=(size, 5)).astype(float)
        
        adata = ad.AnnData(
            X=expression,
            obs=pd.DataFrame(index=[f'spot_{i}' for i in range(size)]),
            var=pd.DataFrame(index=[f'GENE_RANDOM_{i}' for i in range(5)])
        )
        adata.obsm['spatial'] = coords
        return adata
    
    def _create_checkerboard_data(self, size: int = 100) -> ad.AnnData:
        """Create checkerboard pattern (negative autocorrelation)."""
        np.random.seed(789)
        
        # Create regular grid
        grid_dim = int(np.sqrt(size))
        x, y = np.meshgrid(np.arange(grid_dim), np.arange(grid_dim))
        coords = np.column_stack([x.ravel()[:size], y.ravel()[:size]])
        
        # Create checkerboard pattern
        expression = np.zeros((size, 3))
        for i in range(3):
            checkerboard = ((coords[:, 0] + coords[:, 1]) % 2) * 20 + np.random.normal(0, 2, size)
            expression[:, i] = checkerboard
        
        adata = ad.AnnData(
            X=expression,
            obs=pd.DataFrame(index=[f'spot_{i}' for i in range(size)]),
            var=pd.DataFrame(index=[f'GENE_CHECKER_{i}' for i in range(3)])
        )
        adata.obsm['spatial'] = coords.astype(float)
        return adata
    
    def _create_hotspot_data(self, size: int = 200) -> ad.AnnData:
        """Create data with distinct hotspots."""
        np.random.seed(101)
        
        coords = np.random.uniform(0, 10, size=(size, 2))
        expression = np.zeros((size, 4))
        
        # Create hotspots at specific locations
        hotspot_centers = [(2, 2), (8, 8), (2, 8), (8, 2)]
        
        for i, center in enumerate(hotspot_centers):
            distances = np.linalg.norm(coords - np.array(center), axis=1)
            hotspot_pattern = np.exp(-distances**2 / 2) * 30
            expression[:, i] = hotspot_pattern + np.random.poisson(5, size)
        
        adata = ad.AnnData(
            X=expression,
            obs=pd.DataFrame(index=[f'spot_{i}' for i in range(size)]),
            var=pd.DataFrame(index=[f'GENE_HOTSPOT_{i}' for i in range(4)])
        )
        adata.obsm['spatial'] = coords
        return adata
    
    def _create_gradient_data(self, size: int = 150) -> ad.AnnData:
        """Create linear and non-linear gradient patterns."""
        np.random.seed(202)
        
        coords = np.random.uniform(0, 10, size=(size, 2))
        expression = np.zeros((size, 4))
        
        # Linear gradients
        expression[:, 0] = coords[:, 0] * 3 + np.random.normal(0, 2, size)  # X gradient
        expression[:, 1] = coords[:, 1] * 3 + np.random.normal(0, 2, size)  # Y gradient
        
        # Non-linear patterns  
        expression[:, 2] = (coords[:, 0]**2 + coords[:, 1]**2) + np.random.normal(0, 3, size)  # Radial
        expression[:, 3] = np.sin(coords[:, 0]) * np.cos(coords[:, 1]) * 20 + np.random.normal(0, 2, size)  # Periodic
        
        adata = ad.AnnData(
            X=expression,
            obs=pd.DataFrame(index=[f'spot_{i}' for i in range(size)]),
            var=pd.DataFrame(index=['GRAD_X', 'GRAD_Y', 'RADIAL', 'PERIODIC'])
        )
        adata.obsm['spatial'] = coords
        return adata
    
    def _create_multiscale_data(self, size: int = 300) -> ad.AnnData:
        """Create multi-scale spatial patterns."""
        np.random.seed(303)
        
        coords = np.random.uniform(0, 15, size=(size, 2))
        expression = np.zeros((size, 3))
        
        # Large scale pattern
        large_scale = np.sin(coords[:, 0] * 0.5) * np.cos(coords[:, 1] * 0.5) * 20
        expression[:, 0] = large_scale + np.random.normal(0, 2, size)
        
        # Medium scale pattern
        medium_scale = np.sin(coords[:, 0]) * np.cos(coords[:, 1]) * 15
        expression[:, 1] = medium_scale + np.random.normal(0, 2, size)
        
        # Small scale pattern
        small_scale = np.sin(coords[:, 0] * 2) * np.cos(coords[:, 1] * 2) * 10
        expression[:, 2] = small_scale + np.random.normal(0, 2, size)
        
        adata = ad.AnnData(
            X=expression,
            obs=pd.DataFrame(index=[f'spot_{i}' for i in range(size)]),
            var=pd.DataFrame(index=['LARGE_SCALE', 'MEDIUM_SCALE', 'SMALL_SCALE'])
        )
        adata.obsm['spatial'] = coords
        return adata
    
    def validate_moran_properties(self, datasets: Dict[str, ad.AnnData]) -> Dict[str, Any]:
        """Validate mathematical properties of Moran's I statistic."""
        logger.info("=== Validating Moran's I Mathematical Properties ===")
        results = {}
        
        for dataset_name, adata in datasets.items():
            logger.info(f"Testing dataset: {dataset_name}")
            dataset_results = {}
            
            genes = adata.var_names.tolist()
            
            try:
                # Test Moran's I
                moran_results = self.spatial_stats.compute_spatial_autocorrelation(
                    adata=adata,
                    genes=genes,
                    method='moran',
                    n_neighbors=8
                )
                
                if len(moran_results) > 0:
                    dataset_results['moran_statistics'] = {
                        'mean_I': float(moran_results['moran_I'].mean()),
                        'std_I': float(moran_results['moran_I'].std()),
                        'min_I': float(moran_results['moran_I'].min()),
                        'max_I': float(moran_results['moran_I'].max()),
                        'mean_expected_I': float(moran_results['expected_I'].mean()),
                        'significant_genes': int((moran_results['p_value'] < 0.05).sum()),
                        'total_genes': len(moran_results)
                    }
                    
                    # Validate expected properties
                    validations = {}
                    
                    # Expected I should be approximately -1/(n-1)
                    n_obs = adata.n_obs
                    theoretical_expected = -1.0 / (n_obs - 1)
                    actual_expected = moran_results['expected_I'].mean()
                    validations['expected_I_correct'] = abs(actual_expected - theoretical_expected) < 0.01
                    
                    # For perfectly autocorrelated data, I should be high
                    if 'autocorr' in dataset_name:
                        validations['high_autocorr_detected'] = moran_results['moran_I'].mean() > 0.3
                    
                    # For random data, I should be near expected value
                    if 'random' in dataset_name:
                        mean_deviation = abs(moran_results['moran_I'].mean() - actual_expected)
                        validations['random_near_expected'] = mean_deviation < 0.2
                    
                    # For checkerboard, I should be negative
                    if 'checker' in dataset_name:
                        validations['negative_autocorr_detected'] = moran_results['moran_I'].mean() < actual_expected
                    
                    dataset_results['validations'] = validations
                    
                    logger.info(f"  Moran's I range: [{dataset_results['moran_statistics']['min_I']:.3f}, "
                              f"{dataset_results['moran_statistics']['max_I']:.3f}]")
                    logger.info(f"  Significant genes: {dataset_results['moran_statistics']['significant_genes']}"
                              f"/{dataset_results['moran_statistics']['total_genes']}")
                
            except Exception as e:
                logger.error(f"Moran's I validation failed for {dataset_name}: {e}")
                dataset_results['error'] = str(e)
            
            results[dataset_name] = dataset_results
        
        return results
    
    def validate_geary_properties(self, datasets: Dict[str, ad.AnnData]) -> Dict[str, Any]:
        """Validate mathematical properties of Geary's C statistic."""
        logger.info("=== Validating Geary's C Mathematical Properties ===")
        results = {}
        
        for dataset_name, adata in datasets.items():
            logger.info(f"Testing dataset: {dataset_name}")
            dataset_results = {}
            
            genes = adata.var_names[:3].tolist()  # Test subset for speed
            
            try:
                # Test Geary's C
                geary_results = self.spatial_stats.compute_spatial_autocorrelation(
                    adata=adata,
                    genes=genes,
                    method='geary',
                    n_neighbors=8
                )
                
                if len(geary_results) > 0:
                    dataset_results['geary_statistics'] = {
                        'mean_C': float(geary_results['geary_C'].mean()),
                        'std_C': float(geary_results['geary_C'].std()),
                        'min_C': float(geary_results['geary_C'].min()),
                        'max_C': float(geary_results['geary_C'].max()),
                        'mean_expected_C': float(geary_results['expected_C'].mean()),
                        'significant_genes': int((geary_results['p_value'] < 0.05).sum())
                    }
                    
                    # Validate expected properties
                    validations = {}
                    
                    # Expected C should be approximately 1
                    validations['expected_C_correct'] = abs(geary_results['expected_C'].mean() - 1.0) < 0.1
                    
                    # For positively autocorrelated data, C < 1
                    if 'autocorr' in dataset_name or 'gradient' in dataset_name or 'hotspot' in dataset_name:
                        validations['positive_autocorr_detected'] = geary_results['geary_C'].mean() < 1.0
                    
                    # For negatively autocorrelated data, C > 1
                    if 'checker' in dataset_name:
                        validations['negative_autocorr_detected'] = geary_results['geary_C'].mean() > 1.0
                    
                    dataset_results['validations'] = validations
                    
                    logger.info(f"  Geary's C range: [{dataset_results['geary_statistics']['min_C']:.3f}, "
                              f"{dataset_results['geary_statistics']['max_C']:.3f}]")
                
            except Exception as e:
                logger.error(f"Geary's C validation failed for {dataset_name}: {e}")
                dataset_results['error'] = str(e)
            
            results[dataset_name] = dataset_results
        
        return results
    
    def validate_local_statistics(self, datasets: Dict[str, ad.AnnData]) -> Dict[str, Any]:
        """Validate Local Indicators of Spatial Association (LISA)."""
        logger.info("=== Validating Local Spatial Statistics ===")
        results = {}
        
        for dataset_name, adata in datasets.items():
            logger.info(f"Testing dataset: {dataset_name}")
            dataset_results = {}
            
            genes = adata.var_names[:2].tolist()
            
            try:
                # Test Local Moran's I
                local_moran = self.spatial_stats.local_spatial_statistics(
                    adata=adata,
                    genes=genes,
                    method='local_moran',
                    n_neighbors=8
                )
                
                if local_moran:
                    moran_stats = {}
                    for gene in genes:
                        if gene in local_moran:
                            local_Is = local_moran[gene]['Is']
                            p_vals = local_moran[gene]['p_values']
                            clusters = local_moran[gene]['clusters']
                            
                            moran_stats[gene] = {
                                'mean_local_I': float(np.mean(local_Is)),
                                'significant_spots': int(np.sum(p_vals < 0.05)),
                                'total_spots': len(local_Is),
                                'cluster_types': list(np.unique(clusters))
                            }
                    
                    dataset_results['local_moran'] = moran_stats
                
                # Test Getis-Ord G*
                getis_ord = self.spatial_stats.local_spatial_statistics(
                    adata=adata,
                    genes=genes,
                    method='getis_ord',
                    n_neighbors=8
                )
                
                if getis_ord:
                    getis_stats = {}
                    for gene in genes:
                        if gene in getis_ord:
                            gi_values = getis_ord[gene]['Gi']
                            p_vals = getis_ord[gene]['p_values']
                            z_scores = getis_ord[gene]['z_scores']
                            
                            getis_stats[gene] = {
                                'mean_Gi': float(np.mean(gi_values)),
                                'mean_z_score': float(np.mean(z_scores)),
                                'significant_spots': int(np.sum(p_vals < 0.05)),
                                'hotspots': int(np.sum((p_vals < 0.05) & (z_scores > 1.96))),
                                'coldspots': int(np.sum((p_vals < 0.05) & (z_scores < -1.96)))
                            }
                    
                    dataset_results['getis_ord'] = getis_stats
                
                # Validate that results are written to adata.obs
                validations = {}
                for gene in genes:
                    local_moran_cols = [f'{gene}_local_moran_I', f'{gene}_local_moran_pval', f'{gene}_local_moran_cluster']
                    getis_ord_cols = [f'{gene}_getis_ord_Gi', f'{gene}_getis_ord_pval', f'{gene}_getis_ord_zscore']
                    
                    validations[f'{gene}_local_moran_in_obs'] = all(col in adata.obs.columns for col in local_moran_cols)
                    validations[f'{gene}_getis_ord_in_obs'] = all(col in adata.obs.columns for col in getis_ord_cols)
                
                dataset_results['validations'] = validations
                
            except Exception as e:
                logger.error(f"Local statistics validation failed for {dataset_name}: {e}")
                dataset_results['error'] = str(e)
            
            results[dataset_name] = dataset_results
        
        return results
    
    def validate_spatial_weights_matrix(self, test_coords: np.ndarray, k_values: List[int]) -> Dict[str, Any]:
        """Validate spatial weights matrix properties."""
        logger.info("=== Validating Spatial Weights Matrix ===")
        results = {}
        
        try:
            from libpysal.weights import KNN
            
            for k in k_values:
                logger.info(f"Testing k={k} neighbors")
                k_results = {}
                
                try:
                    # Create weights matrix
                    w = KNN.from_array(test_coords, k=k)
                    w.transform = 'R'  # Row standardization
                    
                    # Test basic properties
                    k_results['n_observations'] = w.n
                    k_results['total_neighbors'] = sum(len(neighbors) for neighbors in w.neighbors)
                    k_results['average_neighbors'] = k_results['total_neighbors'] / w.n
                    
                    # Test row standardization
                    row_sums = []
                    for i in range(w.n):
                        if len(w.weights[i]) > 0:
                            row_sums.append(sum(w.weights[i]))
                    
                    k_results['row_sum_mean'] = float(np.mean(row_sums)) if row_sums else 0
                    k_results['row_sum_std'] = float(np.std(row_sums)) if row_sums else 0
                    k_results['properly_standardized'] = abs(k_results['row_sum_mean'] - 1.0) < 1e-10
                    
                    # Test connectivity
                    k_results['min_neighbors'] = min(len(neighbors) for neighbors in w.neighbors)
                    k_results['max_neighbors'] = max(len(neighbors) for neighbors in w.neighbors)
                    k_results['fully_connected'] = k_results['min_neighbors'] > 0
                    
                    logger.info(f"  k={k}: {k_results['average_neighbors']:.1f} avg neighbors, "
                              f"standardized={k_results['properly_standardized']}")
                    
                except Exception as e:
                    logger.error(f"Weights matrix test failed for k={k}: {e}")
                    k_results['error'] = str(e)
                
                results[f'k_{k}'] = k_results
            
        except ImportError:
            logger.warning("libpysal not available for weights matrix testing")
            results['error'] = "libpysal not available"
        
        return results
    
    def performance_benchmark(self, datasets: Dict[str, ad.AnnData]) -> Dict[str, Any]:
        """Comprehensive performance benchmarking."""
        logger.info("=== Performance Benchmarking ===")
        results = {}
        
        for dataset_name, adata in datasets.items():
            logger.info(f"Benchmarking dataset: {dataset_name}")
            dataset_results = {}
            
            # Benchmark different methods
            methods_to_test = [
                ('moran', 'compute_spatial_autocorrelation'),
                ('geary', 'compute_spatial_autocorrelation'),
                ('local_moran', 'local_spatial_statistics'),
                ('getis_ord', 'local_spatial_statistics')
            ]
            
            test_genes = adata.var_names[:5].tolist()
            
            for method_name, function_name in methods_to_test:
                try:
                    start_time = time.time()
                    
                    if function_name == 'compute_spatial_autocorrelation':
                        results_data = self.spatial_stats.compute_spatial_autocorrelation(
                            adata=adata,
                            genes=test_genes,
                            method=method_name,
                            n_neighbors=15
                        )
                    else:  # local_spatial_statistics
                        results_data = self.spatial_stats.local_spatial_statistics(
                            adata=adata,
                            genes=test_genes[:3],  # Use fewer genes for local methods
                            method=method_name,
                            n_neighbors=15
                        )
                    
                    elapsed_time = time.time() - start_time
                    
                    dataset_results[method_name] = {
                        'execution_time': elapsed_time,
                        'time_per_gene': elapsed_time / len(test_genes) if test_genes else 0,
                        'results_generated': len(results_data) if hasattr(results_data, '__len__') else bool(results_data)
                    }
                    
                    logger.info(f"  {method_name}: {elapsed_time:.3f}s for {len(test_genes)} genes")
                    
                except Exception as e:
                    logger.warning(f"Performance test failed for {method_name} on {dataset_name}: {e}")
                    dataset_results[method_name] = {'error': str(e)}
            
            results[dataset_name] = dataset_results
        
        return results
    
    def run_comprehensive_validation(self) -> Dict[str, Any]:
        """Run all validation tests and compile results."""
        logger.info("Starting Comprehensive Spatial Statistics Validation")
        
        # Create test datasets
        logger.info("Creating test datasets...")
        datasets = self.create_test_datasets()
        logger.info(f"Created {len(datasets)} test datasets")
        
        validation_results = {}
        
        # 1. Validate Moran's I properties
        validation_results['moran_validation'] = self.validate_moran_properties(datasets)
        
        # 2. Validate Geary's C properties
        validation_results['geary_validation'] = self.validate_geary_properties(datasets)
        
        # 3. Validate local statistics
        validation_results['local_statistics_validation'] = self.validate_local_statistics(datasets)
        
        # 4. Validate spatial weights matrix
        test_coords = np.random.uniform(0, 10, size=(100, 2))
        validation_results['weights_matrix_validation'] = self.validate_spatial_weights_matrix(
            test_coords, [5, 10, 15, 20, 30]
        )
        
        # 5. Performance benchmarking
        validation_results['performance_benchmark'] = self.performance_benchmark(datasets)
        
        # 6. Dataset information
        validation_results['datasets_info'] = {}
        for name, adata in datasets.items():
            validation_results['datasets_info'][name] = {
                'n_obs': adata.n_obs,
                'n_vars': adata.n_vars,
                'spatial_coords_shape': adata.obsm['spatial'].shape,
                'description': self._get_dataset_description(name)
            }
        
        return validation_results
    
    def _get_dataset_description(self, dataset_name: str) -> str:
        """Get description for dataset."""
        descriptions = {
            'perfect_autocorr': 'Data with perfect positive spatial autocorrelation (smooth gradients)',
            'random': 'Randomly distributed data with no spatial structure',
            'checkerboard': 'Checkerboard pattern with negative spatial autocorrelation',
            'hotspots': 'Data with distinct spatial hotspots',
            'gradient': 'Linear and non-linear gradient patterns',
            'multiscale': 'Multi-scale spatial patterns (large, medium, small)'
        }
        return descriptions.get(dataset_name, 'Unknown pattern')
    
    def generate_validation_report(self, results: Dict[str, Any]) -> str:
        """Generate comprehensive validation report."""
        report_lines = []
        
        report_lines.append("# Comprehensive Spatial Statistics Validation Report")
        report_lines.append(f"Generated on: {pd.Timestamp.now()}")
        report_lines.append("")
        
        # Summary
        report_lines.append("## Executive Summary")
        total_tests = 0
        passed_tests = 0
        
        # Count validations in Moran results
        for dataset_name, result in results['moran_validation'].items():
            if 'validations' in result:
                for validation_name, passed in result['validations'].items():
                    total_tests += 1
                    if passed:
                        passed_tests += 1
        
        # Count validations in Geary results  
        for dataset_name, result in results['geary_validation'].items():
            if 'validations' in result:
                for validation_name, passed in result['validations'].items():
                    total_tests += 1
                    if passed:
                        passed_tests += 1
        
        success_rate = (passed_tests / total_tests * 100) if total_tests > 0 else 0
        report_lines.append(f"- **Overall Success Rate**: {success_rate:.1f}% ({passed_tests}/{total_tests} tests passed)")
        report_lines.append("")
        
        # Dataset Information
        report_lines.append("## Test Datasets")
        for name, info in results['datasets_info'].items():
            report_lines.append(f"- **{name}**: {info['n_obs']} spots, {info['n_vars']} genes - {info['description']}")
        report_lines.append("")
        
        # Moran's I Validation
        report_lines.append("## Moran's I Validation Results")
        for dataset_name, result in results['moran_validation'].items():
            if 'error' not in result:
                stats = result['moran_statistics']
                report_lines.append(f"### {dataset_name}")
                report_lines.append(f"- Mean Moran's I: {stats['mean_I']:.3f}")
                report_lines.append(f"- Range: [{stats['min_I']:.3f}, {stats['max_I']:.3f}]")
                report_lines.append(f"- Significant genes: {stats['significant_genes']}/{stats['total_genes']}")
                
                if 'validations' in result:
                    report_lines.append("- Validations:")
                    for validation_name, passed in result['validations'].items():
                        status = "✓ PASS" if passed else "✗ FAIL"
                        report_lines.append(f"  - {validation_name}: {status}")
                report_lines.append("")
        
        # Geary's C Validation
        report_lines.append("## Geary's C Validation Results")
        for dataset_name, result in results['geary_validation'].items():
            if 'error' not in result:
                stats = result['geary_statistics']
                report_lines.append(f"### {dataset_name}")
                report_lines.append(f"- Mean Geary's C: {stats['mean_C']:.3f}")
                report_lines.append(f"- Range: [{stats['min_C']:.3f}, {stats['max_C']:.3f}]")
                report_lines.append(f"- Significant genes: {stats['significant_genes']}")
                report_lines.append("")
        
        # Performance Results
        report_lines.append("## Performance Benchmark Results")
        for dataset_name, result in results['performance_benchmark'].items():
            report_lines.append(f"### {dataset_name}")
            for method_name, perf_data in result.items():
                if 'error' not in perf_data:
                    report_lines.append(f"- {method_name}: {perf_data['execution_time']:.3f}s "
                                      f"({perf_data['time_per_gene']:.3f}s per gene)")
            report_lines.append("")
        
        # Spatial Weights Matrix Validation
        report_lines.append("## Spatial Weights Matrix Validation")
        weights_results = results['weights_matrix_validation']
        if 'error' not in weights_results:
            for k_name, k_result in weights_results.items():
                if 'error' not in k_result:
                    report_lines.append(f"- {k_name}: {k_result['average_neighbors']:.1f} avg neighbors, "
                                      f"properly standardized: {k_result['properly_standardized']}")
        report_lines.append("")
        
        return "\n".join(report_lines)

def main():
    """Run the comprehensive validation."""
    validator = SpatialStatisticsValidator()
    
    print("Starting comprehensive spatial statistics validation...")
    print("This will test mathematical correctness, statistical properties, and performance.")
    print()
    
    # Run validation
    results = validator.run_comprehensive_validation()
    
    # Generate and save report
    report_text = validator.generate_validation_report(results)
    
    # Save results
    output_dir = Path(__file__).parent
    
    # Save JSON results
    results_file = output_dir / 'spatial_statistics_validation_results.json'
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    # Save text report
    report_file = output_dir / 'spatial_statistics_validation_report.md'
    with open(report_file, 'w') as f:
        f.write(report_text)
    
    print("\n" + "="*60)
    print("COMPREHENSIVE VALIDATION COMPLETED")
    print("="*60)
    print(report_text)
    print("\n" + "="*60)
    print(f"Detailed results saved to: {results_file}")
    print(f"Human-readable report saved to: {report_file}")

if __name__ == "__main__":
    main()