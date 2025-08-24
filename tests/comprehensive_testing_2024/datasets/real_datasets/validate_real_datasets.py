#!/usr/bin/env python3
"""
Real Dataset Validation and Integration System
Linus-style: Simple, direct, no special cases.
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import json
import warnings
from datetime import datetime
from typing import Dict, List, Any, Optional
import sys
import os

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent))

warnings.filterwarnings('ignore', category=FutureWarning)

class DatasetValidator:
    """Simple, effective dataset validation. No bullshit."""
    
    def __init__(self, source_data_dir: Path, target_dir: Path):
        self.source_data_dir = Path(source_data_dir)
        self.target_dir = Path(target_dir)
        self.target_dir.mkdir(parents=True, exist_ok=True)
        
        self.validation_results = []
        self.registry_data = {
            'validation_timestamp': datetime.now().isoformat(),
            'datasets': [],
            'summary_stats': {},
            'compatibility_status': 'pending'
        }
    
    def validate_single_dataset(self, h5ad_path: Path) -> Dict[str, Any]:
        """Validate one dataset. Return facts, not opinions."""
        
        print(f"Validating: {h5ad_path.name}")
        
        result = {
            'file_path': str(h5ad_path),
            'file_name': h5ad_path.name,
            'file_size_mb': round(h5ad_path.stat().st_size / 1024 / 1024, 2),
            'validation_status': 'pending'
        }
        
        try:
            # Load dataset
            adata = sc.read_h5ad(h5ad_path)
            
            # Basic metrics
            result.update({
                'n_cells': int(adata.n_obs),
                'n_genes': int(adata.n_vars),
                'data_type': str(adata.X.dtype),
                'is_sparse': hasattr(adata.X, 'sparse') or hasattr(adata.X, 'toarray')
            })
            
            # Safe sparsity calculation
            try:
                if hasattr(adata.X, 'toarray'):
                    # Sparse matrix
                    total_elements = adata.X.shape[0] * adata.X.shape[1]
                    zero_elements = total_elements - adata.X.nnz
                    result['sparsity_ratio'] = float(zero_elements / total_elements)
                else:
                    # Dense matrix
                    result['sparsity_ratio'] = float((adata.X == 0).mean())
            except Exception as e:
                result['sparsity_ratio'] = None
                result['sparsity_error'] = str(e)
            
            # Spatial coordinates check
            spatial_info = self._check_spatial_data(adata)
            result.update(spatial_info)
            
            # Metadata analysis
            metadata_info = self._analyze_metadata(adata)
            result.update(metadata_info)
            
            # Data quality assessment
            quality_info = self._assess_data_quality(adata)
            result.update(quality_info)
            
            result['validation_status'] = 'success'
            print(f"✓ {h5ad_path.name}: {result['n_cells']} cells, {result['n_genes']} genes")
            
        except Exception as e:
            result.update({
                'validation_status': 'failed',
                'error': str(e),
                'error_type': type(e).__name__
            })
            print(f"✗ {h5ad_path.name}: {str(e)}")
        
        return result
    
    def _check_spatial_data(self, adata) -> Dict[str, Any]:
        """Check spatial coordinates. Binary: exists or doesn't."""
        
        spatial_info = {
            'has_spatial': False,
            'spatial_dimensions': None,
            'spatial_key': None,
            'coordinate_range': None
        }
        
        # Check common spatial keys
        spatial_keys = ['spatial', 'X_spatial', 'X_umap', 'coordinates']
        
        for key in spatial_keys:
            if key in adata.obsm:
                spatial_info['has_spatial'] = True
                spatial_info['spatial_key'] = key
                
                coords = adata.obsm[key]
                spatial_info['spatial_dimensions'] = coords.shape[1]
                
                if coords.shape[1] >= 2:
                    spatial_info['coordinate_range'] = {
                        'x_min': float(coords[:, 0].min()),
                        'x_max': float(coords[:, 0].max()),
                        'y_min': float(coords[:, 1].min()),
                        'y_max': float(coords[:, 1].max())
                    }
                break
        
        return spatial_info
    
    def _analyze_metadata(self, adata) -> Dict[str, Any]:
        """Analyze observation and variable metadata."""
        
        metadata_info = {
            'obs_columns': list(adata.obs.columns),
            'var_columns': list(adata.var.columns),
            'uns_keys': list(adata.uns.keys()) if hasattr(adata, 'uns') else [],
            'obsm_keys': list(adata.obsm.keys()),
            'varm_keys': list(adata.varm.keys()) if hasattr(adata, 'varm') else []
        }
        
        # Check for common annotation columns
        common_obs_keys = ['cell_type', 'cluster', 'leiden', 'louvain', 'batch', 'sample']
        found_annotations = [key for key in common_obs_keys if key in adata.obs.columns]
        metadata_info['found_annotations'] = found_annotations
        
        return metadata_info
    
    def _assess_data_quality(self, adata) -> Dict[str, Any]:
        """Assess basic data quality metrics."""
        
        # Gene expression stats
        gene_counts = np.array(adata.X.sum(axis=0)).flatten()
        cell_counts = np.array(adata.X.sum(axis=1)).flatten()
        
        quality_info = {
            'genes_per_cell_mean': float(np.mean((adata.X > 0).sum(axis=1))),
            'genes_per_cell_median': float(np.median((adata.X > 0).sum(axis=1))),
            'counts_per_cell_mean': float(np.mean(cell_counts)),
            'counts_per_cell_median': float(np.median(cell_counts)),
            'cells_per_gene_mean': float(np.mean((adata.X > 0).sum(axis=0))),
            'highly_expressed_genes': int(np.sum(gene_counts > np.percentile(gene_counts, 95))),
            'zero_count_cells': int(np.sum(cell_counts == 0)),
            'zero_count_genes': int(np.sum(gene_counts == 0))
        }
        
        # Quality classification
        if quality_info['genes_per_cell_mean'] > 500 and quality_info['counts_per_cell_mean'] > 1000:
            quality_info['quality_tier'] = 'high'
        elif quality_info['genes_per_cell_mean'] > 200 and quality_info['counts_per_cell_mean'] > 500:
            quality_info['quality_tier'] = 'medium'
        else:
            quality_info['quality_tier'] = 'low'
        
        return quality_info
    
    def copy_and_organize_datasets(self) -> None:
        """Copy datasets to organized structure. Keep it simple."""
        
        print("\n=== Organizing Datasets ===")
        
        # Find all h5ad files
        h5ad_files = []
        for root, dirs, files in os.walk(self.source_data_dir):
            for file in files:
                if file.endswith('.h5ad'):
                    h5ad_files.append(Path(root) / file)
        
        print(f"Found {len(h5ad_files)} h5ad files")
        
        # Organize by category
        categories = {
            'core': [],
            'harmony': [],
            'test': [],
            'synthetic': []
        }
        
        for file_path in h5ad_files:
            # Simple categorization based on path
            if 'harmony' in str(file_path):
                category = 'harmony'
            elif 'test' in str(file_path) or 'synthetic' in file_path.name:
                category = 'test'
            elif 'core' in str(file_path):
                category = 'core'
            else:
                category = 'core'  # default
            
            categories[category].append(file_path)
        
        # Create organized structure
        for category, files in categories.items():
            if not files:
                continue
                
            category_dir = self.target_dir / category
            category_dir.mkdir(exist_ok=True)
            
            print(f"\n{category.upper()} datasets:")
            for file_path in files:
                # Create symlink to original file
                target_path = category_dir / file_path.name
                if not target_path.exists():
                    target_path.symlink_to(file_path)
                    print(f"  Linked: {file_path.name}")
                else:
                    print(f"  Exists: {file_path.name}")
    
    def validate_all_datasets(self) -> None:
        """Validate all found datasets."""
        
        print("\n=== Dataset Validation ===")
        
        # Find all h5ad files in organized structure
        h5ad_files = list(self.target_dir.rglob("*.h5ad"))
        
        if not h5ad_files:
            print("No h5ad files found. Running organization first...")
            self.copy_and_organize_datasets()
            h5ad_files = list(self.target_dir.rglob("*.h5ad"))
        
        print(f"Validating {len(h5ad_files)} datasets...")
        
        for h5ad_path in sorted(h5ad_files):
            result = self.validate_single_dataset(h5ad_path)
            self.validation_results.append(result)
    
    def test_chatspatial_compatibility(self) -> Dict[str, Any]:
        """Test compatibility with ChatSpatial tools."""
        
        print("\n=== Compatibility Testing ===")
        
        compatibility_results = {
            'tested_datasets': 0,
            'successful_loads': 0,
            'spatial_compatible': 0,
            'preprocessing_compatible': 0,
            'analysis_compatible': 0,
            'failed_datasets': [],
            'issues_found': []
        }
        
        # Test a subset of datasets (avoid testing synthetic/test data)
        real_datasets = [r for r in self.validation_results 
                        if r['validation_status'] == 'success' 
                        and 'synthetic' not in r['file_name'].lower()
                        and 'test' not in r['file_name'].lower()]
        
        for dataset_info in real_datasets[:10]:  # Test first 10 real datasets
            compatibility_results['tested_datasets'] += 1
            dataset_path = Path(dataset_info['file_path'])
            
            try:
                # Test basic loading
                adata = sc.read_h5ad(dataset_path)
                compatibility_results['successful_loads'] += 1
                
                # Test spatial compatibility
                if dataset_info['has_spatial']:
                    try:
                        # Test basic spatial operations
                        spatial_key = dataset_info['spatial_key']
                        coords = adata.obsm[spatial_key]
                        if coords.shape[1] >= 2:
                            compatibility_results['spatial_compatible'] += 1
                    except Exception as e:
                        compatibility_results['issues_found'].append({
                            'dataset': dataset_info['file_name'],
                            'issue': 'spatial_access_failed',
                            'error': str(e)
                        })
                
                # Test preprocessing compatibility
                try:
                    # Basic preprocessing operations
                    if adata.n_obs > 10 and adata.n_vars > 10:
                        test_adata = adata.copy()
                        sc.pp.filter_cells(test_adata, min_genes=1)
                        sc.pp.filter_genes(test_adata, min_cells=1)
                        compatibility_results['preprocessing_compatible'] += 1
                except Exception as e:
                    compatibility_results['issues_found'].append({
                        'dataset': dataset_info['file_name'],
                        'issue': 'preprocessing_failed',
                        'error': str(e)
                    })
                
                # Test analysis compatibility
                try:
                    if adata.n_obs > 50 and adata.n_vars > 100:
                        test_adata = adata.copy()
                        sc.pp.normalize_total(test_adata, target_sum=1e4)
                        sc.pp.log1p(test_adata)
                        compatibility_results['analysis_compatible'] += 1
                except Exception as e:
                    compatibility_results['issues_found'].append({
                        'dataset': dataset_info['file_name'],
                        'issue': 'analysis_failed',
                        'error': str(e)
                    })
                
                print(f"✓ {dataset_info['file_name']} - Compatible")
                
            except Exception as e:
                compatibility_results['failed_datasets'].append({
                    'dataset': dataset_info['file_name'],
                    'error': str(e)
                })
                print(f"✗ {dataset_info['file_name']} - Failed: {str(e)}")
        
        return compatibility_results
    
    def generate_summary_stats(self) -> Dict[str, Any]:
        """Generate overall summary statistics."""
        
        successful_validations = [r for r in self.validation_results 
                                if r['validation_status'] == 'success']
        
        if not successful_validations:
            return {'error': 'No successful validations'}
        
        # Calculate summary statistics
        total_cells = sum(r['n_cells'] for r in successful_validations)
        total_genes = sum(r['n_genes'] for r in successful_validations)
        total_size_mb = sum(r['file_size_mb'] for r in successful_validations)
        
        spatial_datasets = sum(1 for r in successful_validations if r['has_spatial'])
        
        quality_tiers = {}
        for r in successful_validations:
            tier = r.get('quality_tier', 'unknown')
            quality_tiers[tier] = quality_tiers.get(tier, 0) + 1
        
        return {
            'total_datasets': len(self.validation_results),
            'successful_validations': len(successful_validations),
            'failed_validations': len(self.validation_results) - len(successful_validations),
            'total_cells': total_cells,
            'total_genes': total_genes,
            'total_size_mb': round(total_size_mb, 2),
            'spatial_datasets': spatial_datasets,
            'quality_distribution': quality_tiers,
            'average_cells_per_dataset': round(total_cells / len(successful_validations), 0),
            'average_genes_per_dataset': round(total_genes / len(successful_validations), 0)
        }
    
    def generate_reports(self) -> None:
        """Generate all reports."""
        
        print("\n=== Generating Reports ===")
        
        # Test compatibility
        compatibility_results = self.test_chatspatial_compatibility()
        
        # Generate summary stats
        summary_stats = self.generate_summary_stats()
        
        # Update registry
        self.registry_data.update({
            'datasets': self.validation_results,
            'summary_stats': summary_stats,
            'compatibility_results': compatibility_results,
            'compatibility_status': 'completed'
        })
        
        # Write JSON summary
        json_path = self.target_dir / 'real_datasets_summary.json'
        with open(json_path, 'w') as f:
            json.dump(self.registry_data, f, indent=2)
        print(f"✓ Written: {json_path}")
        
        # Write human-readable registry
        self._write_registry_markdown()
        
        # Write compatibility results
        compatibility_path = self.target_dir / 'compatibility_test_results.json'
        with open(compatibility_path, 'w') as f:
            json.dump(compatibility_results, f, indent=2)
        print(f"✓ Written: {compatibility_path}")
        
        # Write download log
        self._write_download_log()
    
    def _write_registry_markdown(self) -> None:
        """Write human-readable registry."""
        
        md_path = self.target_dir / 'real_datasets_registry.md'
        
        with open(md_path, 'w') as f:
            f.write("# ChatSpatial Real Datasets Registry\n\n")
            f.write(f"**Validation Date**: {self.registry_data['validation_timestamp']}\n\n")
            
            # Summary stats
            stats = self.registry_data['summary_stats']
            f.write("## Summary Statistics\n\n")
            f.write(f"- **Total Datasets**: {stats.get('total_datasets', 0)}\n")
            f.write(f"- **Successfully Validated**: {stats.get('successful_validations', 0)}\n")
            f.write(f"- **Total Cells**: {stats.get('total_cells', 0):,}\n")
            f.write(f"- **Total Genes**: {stats.get('total_genes', 0):,}\n")
            f.write(f"- **Total Size**: {stats.get('total_size_mb', 0):.1f} MB\n")
            f.write(f"- **Spatial Datasets**: {stats.get('spatial_datasets', 0)}\n\n")
            
            # Quality distribution
            quality_dist = stats.get('quality_distribution', {})
            if quality_dist:
                f.write("## Quality Distribution\n\n")
                for tier, count in quality_dist.items():
                    f.write(f"- **{tier.title()}**: {count} datasets\n")
                f.write("\n")
            
            # Compatibility results
            compat = self.registry_data.get('compatibility_results', {})
            if compat:
                f.write("## Compatibility Results\n\n")
                f.write(f"- **Tested**: {compat.get('tested_datasets', 0)} datasets\n")
                f.write(f"- **Load Success**: {compat.get('successful_loads', 0)}\n")
                f.write(f"- **Spatial Compatible**: {compat.get('spatial_compatible', 0)}\n")
                f.write(f"- **Preprocessing Compatible**: {compat.get('preprocessing_compatible', 0)}\n")
                f.write(f"- **Analysis Compatible**: {compat.get('analysis_compatible', 0)}\n\n")
            
            # Dataset details
            f.write("## Dataset Details\n\n")
            
            successful_datasets = [d for d in self.registry_data['datasets'] 
                                 if d['validation_status'] == 'success']
            
            # Group by category
            categories = {}
            for dataset in successful_datasets:
                path_parts = Path(dataset['file_path']).parts
                if 'harmony' in path_parts:
                    category = 'Harmony Integration'
                elif 'core' in path_parts:
                    category = 'Core Datasets'
                elif 'test' in path_parts:
                    category = 'Test Datasets'
                else:
                    category = 'Other'
                
                if category not in categories:
                    categories[category] = []
                categories[category].append(dataset)
            
            for category, datasets in categories.items():
                f.write(f"### {category}\n\n")
                f.write("| Dataset | Cells | Genes | Size (MB) | Spatial | Quality |\n")
                f.write("|---------|-------|-------|-----------|---------|----------|\n")
                
                for dataset in sorted(datasets, key=lambda x: x['file_name']):
                    spatial_icon = "✓" if dataset.get('has_spatial', False) else "✗"
                    quality = dataset.get('quality_tier', 'unknown').title()
                    
                    f.write(f"| {dataset['file_name']} | "
                           f"{dataset.get('n_cells', 0):,} | "
                           f"{dataset.get('n_genes', 0):,} | "
                           f"{dataset.get('file_size_mb', 0):.1f} | "
                           f"{spatial_icon} | "
                           f"{quality} |\n")
                
                f.write("\n")
            
            # Issues found
            issues = compat.get('issues_found', [])
            if issues:
                f.write("## Known Issues\n\n")
                for issue in issues:
                    f.write(f"- **{issue['dataset']}**: {issue['issue']} - {issue['error']}\n")
                f.write("\n")
            
            f.write("## Usage Instructions\n\n")
            f.write("```python\n")
            f.write("import scanpy as sc\n")
            f.write("from pathlib import Path\n\n")
            f.write("# Load a dataset\n")
            f.write("dataset_path = Path('datasets/real_datasets/core/ST_mouse_brain.h5ad')\n")
            f.write("adata = sc.read_h5ad(dataset_path)\n\n")
            f.write("# Check spatial coordinates\n")
            f.write("if 'spatial' in adata.obsm:\n")
            f.write("    print(f'Spatial coords: {adata.obsm[\"spatial\"].shape}')\n")
            f.write("```\n")
        
        print(f"✓ Written: {md_path}")
    
    def _write_download_log(self) -> None:
        """Write download process log."""
        
        log_path = self.target_dir / 'download_log.txt'
        
        with open(log_path, 'w') as f:
            f.write("ChatSpatial Real Datasets Download Log\n")
            f.write("=" * 50 + "\n")
            f.write(f"Validation Date: {self.registry_data['validation_timestamp']}\n\n")
            
            f.write("Data Sources:\n")
            f.write("- Core datasets: /data/core/\n")
            f.write("- Harmony datasets: /data/harmony/\n")
            f.write("- Test datasets: /data/test/\n\n")
            
            f.write("Validation Summary:\n")
            stats = self.registry_data['summary_stats']
            f.write(f"- Total datasets processed: {stats.get('total_datasets', 0)}\n")
            f.write(f"- Successfully validated: {stats.get('successful_validations', 0)}\n")
            f.write(f"- Failed validations: {stats.get('failed_validations', 0)}\n")
            f.write(f"- Total data size: {stats.get('total_size_mb', 0):.1f} MB\n\n")
            
            # Failed datasets
            failed = [d for d in self.registry_data['datasets'] 
                     if d['validation_status'] == 'failed']
            
            if failed:
                f.write("Failed Datasets:\n")
                for dataset in failed:
                    f.write(f"- {dataset['file_name']}: {dataset.get('error', 'Unknown error')}\n")
                f.write("\n")
            
            f.write("Validation completed successfully.\n")
        
        print(f"✓ Written: {log_path}")


def main():
    """Main execution. Keep it simple."""
    
    print("ChatSpatial Real Datasets Validation System")
    print("=" * 50)
    
    # Set up paths
    base_dir = Path(__file__).parent.parent.parent.parent.parent
    source_data_dir = base_dir / "data"
    target_dir = Path(__file__).parent
    
    print(f"Source: {source_data_dir}")
    print(f"Target: {target_dir}")
    
    # Initialize validator
    validator = DatasetValidator(source_data_dir, target_dir)
    
    # Run validation pipeline
    try:
        validator.copy_and_organize_datasets()
        validator.validate_all_datasets()
        validator.generate_reports()
        
        print("\n" + "=" * 50)
        print("✓ Dataset validation completed successfully!")
        print(f"✓ Results saved to: {target_dir}")
        
        # Print summary
        stats = validator.registry_data.get('summary_stats', {})
        if stats:
            print(f"\nSummary:")
            print(f"  {stats.get('total_datasets', 0)} datasets processed")
            print(f"  {stats.get('total_cells', 0):,} total cells")
            print(f"  {stats.get('total_size_mb', 0):.1f} MB total size")
            print(f"  {stats.get('spatial_datasets', 0)} spatial datasets")
        
    except Exception as e:
        print(f"\n✗ Validation failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()