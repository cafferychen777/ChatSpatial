#!/usr/bin/env python3
"""
Demo script for Linus's "Good Taste" Data Validation System

This demonstrates the comprehensive data validation without dependencies on scanpy
"""

import numpy as np
import pandas as pd
from scipy import sparse
from typing import Dict, Any, Optional, List


class MockContext:
    """Mock context for demonstration"""
    def __init__(self):
        self.messages = []
    
    async def info(self, msg):
        self.messages.append(f"INFO: {msg}")
        print(f"ℹ️  {msg}")
    
    async def warning(self, msg):
        self.messages.append(f"WARNING: {msg}")
        print(f"⚠️  {msg}")


class MockAnnData:
    """Mock AnnData class for demonstration"""
    def __init__(self, X, obs, var):
        self.X = X
        self.obs = obs
        self.var = var
        self.obsm = {}
        self.uns = {}
        self.n_obs = X.shape[0]
        self.n_vars = X.shape[1]
        self.obsp = {}


# Inline the validation functions for demo
async def validate_basic_structure(adata: Any, context: Optional = None) -> Dict[str, Any]:
    """Validate basic AnnData structure"""
    result = {"passed": True, "errors": [], "warnings": [], "suggestions": []}
    
    try:
        # Check if it's actually an AnnData object
        if not hasattr(adata, 'X') or not hasattr(adata, 'obs') or not hasattr(adata, 'var'):
            result["errors"].append("Invalid AnnData object: missing X, obs, or var attributes")
            result["suggestions"].append("Ensure you're passing a valid AnnData object")
            result["passed"] = False
            return result
        
        # Check dimensions consistency
        if adata.X.shape[0] != len(adata.obs):
            result["errors"].append(f"Dimension mismatch: X has {adata.X.shape[0]} cells but obs has {len(adata.obs)} entries")
            result["suggestions"].append("Check data loading process - cells and observations must match")
            result["passed"] = False
        
        if adata.X.shape[1] != len(adata.var):
            result["errors"].append(f"Dimension mismatch: X has {adata.X.shape[1]} genes but var has {len(adata.var)} entries")
            result["suggestions"].append("Check data loading process - genes and variables must match")
            result["passed"] = False
        
        # Check for empty data
        if adata.n_obs == 0:
            result["errors"].append("Dataset is empty: no cells found")
            result["suggestions"].append("Load data with actual cell measurements")
            result["passed"] = False
        
        if adata.n_vars == 0:
            result["errors"].append("Dataset is empty: no genes found")
            result["suggestions"].append("Load data with actual gene measurements")
            result["passed"] = False
        
        # Data size warnings
        if adata.n_obs < 50:
            result["warnings"].append(f"Very few cells ({adata.n_obs}) - cell communication analysis may be unreliable")
        
        if adata.n_vars < 1000:
            result["warnings"].append(f"Few genes ({adata.n_vars}) - consider using more comprehensive gene set")
        
    except Exception as e:
        result["errors"].append(f"Structure validation failed: {str(e)}")
        result["suggestions"].append("Check if the input is a valid AnnData object")
        result["passed"] = False
    
    return result


async def validate_expression_matrix(adata: Any, context: Optional = None) -> Dict[str, Any]:
    """Validate expression matrix data quality"""
    result = {"passed": True, "errors": [], "warnings": [], "suggestions": []}
    
    try:
        X = adata.X
        
        # Convert sparse to dense for NaN/Inf checking (sample if too large)
        if sparse.issparse(X):
            # For large matrices, sample to check quality
            if X.shape[0] > 5000 or X.shape[1] > 2000:
                sample_cells = min(1000, X.shape[0])
                sample_genes = min(500, X.shape[1])
                cell_idx = np.random.choice(X.shape[0], sample_cells, replace=False)
                gene_idx = np.random.choice(X.shape[1], sample_genes, replace=False)
                X_sample = X[cell_idx, :][:, gene_idx].toarray()
            else:
                X_sample = X.toarray()
        else:
            X_sample = X if X.size <= 10_000_000 else X[:1000, :500]  # Sample large dense matrices
        
        # Check for NaN values
        if np.isnan(X_sample).any():
            nan_count = np.isnan(X_sample).sum()
            nan_fraction = nan_count / X_sample.size
            if nan_fraction > 0.1:  # More than 10% NaN
                result["errors"].append(f"High proportion of NaN values ({nan_fraction:.2%}) in expression matrix")
                result["suggestions"].append("Remove or impute NaN values before analysis")
                result["passed"] = False
            else:
                result["warnings"].append(f"Found {nan_count} NaN values ({nan_fraction:.2%}) in expression matrix")
        
        # Check for infinite values
        if np.isinf(X_sample).any():
            inf_count = np.isinf(X_sample).sum()
            result["errors"].append(f"Found {inf_count} infinite values in expression matrix")
            result["suggestions"].append("Replace infinite values with finite numbers or remove affected cells/genes")
            result["passed"] = False
        
        # Check for negative values (suspicious in count data)
        if (X_sample < 0).any():
            neg_count = (X_sample < 0).sum()
            neg_fraction = neg_count / X_sample.size
            if neg_fraction > 0.01:  # More than 1% negative
                result["warnings"].append(f"Found {neg_count} negative values ({neg_fraction:.2%}) - unusual for count data")
        
        # Check for extremely large values
        max_val = X_sample.max()
        if max_val > 1e6:
            result["warnings"].append(f"Very large expression values detected (max: {max_val:.2e}) - consider normalization")
        
        # Check data type
        if not np.issubdtype(X.dtype, np.number):
            result["errors"].append(f"Expression matrix has non-numeric data type: {X.dtype}")
            result["suggestions"].append("Convert expression data to numeric format")
            result["passed"] = False
        
    except Exception as e:
        result["errors"].append(f"Expression matrix validation failed: {str(e)}")
        result["suggestions"].append("Check expression matrix format and values")
        result["passed"] = False
    
    return result


async def validate_spatial_coordinates(adata: Any, context: Optional = None) -> Dict[str, Any]:
    """Validate spatial coordinates"""
    result = {"passed": True, "errors": [], "warnings": [], "suggestions": []}
    
    try:
        # Check for spatial coordinates existence
        spatial_key = None
        for key in adata.obsm.keys():
            if 'spatial' in key.lower():
                spatial_key = key
                break
        
        if spatial_key is None:
            result["errors"].append("No spatial coordinates found in adata.obsm")
            result["suggestions"].append("Add spatial coordinates to adata.obsm['spatial'] or similar key")
            result["passed"] = False
            return result
        
        spatial_coords = adata.obsm[spatial_key]
        
        # Check for NaN/Inf in spatial coordinates
        if np.isnan(spatial_coords).any():
            nan_count = np.isnan(spatial_coords).sum()
            result["errors"].append(f"Found {nan_count} NaN values in spatial coordinates")
            result["suggestions"].append("Remove cells with missing spatial coordinates")
            result["passed"] = False
        
        if np.isinf(spatial_coords).any():
            inf_count = np.isinf(spatial_coords).sum()
            result["errors"].append(f"Found {inf_count} infinite values in spatial coordinates")
            result["suggestions"].append("Replace infinite spatial coordinates with valid values")
            result["passed"] = False
        
    except Exception as e:
        result["errors"].append(f"Spatial coordinates validation failed: {str(e)}")
        result["suggestions"].append("Check spatial coordinates format and values")
        result["passed"] = False
    
    return result


async def comprehensive_data_validation(adata: Any, context: Optional = None) -> Dict[str, Any]:
    """Comprehensive data quality validation for AnnData objects"""
    
    validation_result = {
        "passed": True,
        "error_message": "",
        "warnings": [],
        "suggestions": "",
        "validation_details": {}
    }
    
    errors = []
    warnings_list = []
    suggestions = []
    
    try:
        # 1. Basic structure validation
        structure_check = await validate_basic_structure(adata, context)
        validation_result["validation_details"]["structure"] = structure_check
        if not structure_check["passed"]:
            errors.extend(structure_check["errors"])
            suggestions.extend(structure_check["suggestions"])
        warnings_list.extend(structure_check["warnings"])
        
        # 2. Expression matrix validation
        expression_check = await validate_expression_matrix(adata, context)
        validation_result["validation_details"]["expression"] = expression_check
        if not expression_check["passed"]:
            errors.extend(expression_check["errors"])
            suggestions.extend(expression_check["suggestions"])
        warnings_list.extend(expression_check["warnings"])
        
        # 3. Spatial coordinates validation
        spatial_check = await validate_spatial_coordinates(adata, context)
        validation_result["validation_details"]["spatial"] = spatial_check
        if not spatial_check["passed"]:
            errors.extend(spatial_check["errors"])
            suggestions.extend(spatial_check["suggestions"])
        warnings_list.extend(spatial_check["warnings"])
        
        # Compile final results
        if errors:
            validation_result["passed"] = False
            validation_result["error_message"] = "\\n".join([f"• {error}" for error in errors])
        
        validation_result["warnings"] = warnings_list
        if suggestions:
            validation_result["suggestions"] = "\\n".join([f"• {suggestion}" for suggestion in suggestions])
        
        # Log summary
        if context:
            if validation_result["passed"]:
                await context.info(f"✅ Data validation passed with {len(warnings_list)} warnings")
            else:
                await context.warning(f"❌ Data validation failed: {len(errors)} critical issues found")
        
        return validation_result
        
    except Exception as e:
        validation_result["passed"] = False
        validation_result["error_message"] = f"Validation process failed: {str(e)}"
        validation_result["suggestions"] = "Please check your data format and try again"
        return validation_result


def create_good_data():
    """Create well-formed test data"""
    n_obs, n_vars = 200, 1500
    
    X = sparse.random(n_obs, n_vars, density=0.1, format='csr')
    X.data = np.abs(X.data) * 10
    
    obs = pd.DataFrame({
        'cell_type': np.random.choice(['T_cell', 'B_cell', 'NK_cell'], n_obs),
    }, index=[f"CELL_{i:06d}" for i in range(n_obs)])
    
    var = pd.DataFrame({
        'gene_name': [f"GENE_{i:06d}" for i in range(n_vars)],
    }, index=[f"GENE_{i:06d}" for i in range(n_vars)])
    
    adata = MockAnnData(X=X, obs=obs, var=var)
    adata.obsm['spatial'] = np.random.uniform(0, 100, (n_obs, 2))
    
    return adata


def create_bad_data():
    """Create problematic test data"""
    n_obs, n_vars = 30, 200  # Too small
    
    X = sparse.random(n_obs, n_vars, density=0.05, format='csr')
    X_dense = X.toarray()
    X_dense[10:15, 20:25] = np.nan  # NaN values
    X_dense[20:25, 50:55] = np.inf  # Inf values
    X = sparse.csr_matrix(X_dense)
    
    # Duplicate cell IDs
    cell_ids = [f"CELL_{i:06d}" for i in range(n_obs)]
    cell_ids[10] = cell_ids[9]
    
    obs = pd.DataFrame({
        'single_type': ['same_type'] * n_obs,  # Only one type
    }, index=cell_ids)
    
    var = pd.DataFrame(index=[f"GENE_{i}" for i in range(n_vars)])
    
    adata = MockAnnData(X=X, obs=obs, var=var)
    
    # Bad spatial coordinates
    spatial_coords = np.random.uniform(0, 10, (n_obs, 2))
    spatial_coords[15:20, :] = np.nan  # NaN coordinates
    adata.obsm['spatial'] = spatial_coords
    
    return adata


async def demo_validation():
    """Demo the validation system"""
    print("🔬 Linus's 'Good Taste' Data Validation System")
    print("=" * 55)
    print("Based on Linus Torvalds' engineering principles:")
    print("• Check real problems that occur in production")
    print("• Fail early with clear error messages")
    print("• Provide actionable suggestions for fixes")
    print("• Use unified validation framework")
    
    context = MockContext()
    
    # Test 1: Good data
    print("\n🟢 TEST 1: Well-formed data")
    print("-" * 30)
    good_data = create_good_data()
    result = await comprehensive_data_validation(good_data, context)
    
    print(f"Result: {'✅ PASSED' if result['passed'] else '❌ FAILED'}")
    if result['warnings']:
        print(f"Warnings ({len(result['warnings'])}):")
        for warning in result['warnings'][:2]:
            print(f"  ⚠️  {warning}")
    
    # Test 2: Bad data
    print("\n🔴 TEST 2: Problematic data")
    print("-" * 30)
    bad_data = create_bad_data()
    result = await comprehensive_data_validation(bad_data, context)
    
    print(f"Result: {'✅ PASSED' if result['passed'] else '❌ FAILED'}")
    if result['error_message']:
        print("Critical Issues:")
        for line in result['error_message'].split('\\n')[:3]:
            print(f"  🚫 {line}")
    
    if result['suggestions']:
        print("Suggestions:")
        for line in result['suggestions'].split('\\n')[:2]:
            print(f"  💡 {line}")
    
    # Detailed breakdown
    print("\n📊 VALIDATION BREAKDOWN")
    print("-" * 30)
    
    test_data = [("Good Data", good_data), ("Bad Data", bad_data)]
    for name, data in test_data:
        result = await comprehensive_data_validation(data, None)
        print(f"\n{name}:")
        
        for check_type, details in result['validation_details'].items():
            status = "✅" if details['passed'] else "❌"
            errors = len(details['errors'])
            warnings = len(details['warnings'])
            print(f"  {status} {check_type.capitalize()}: {errors} errors, {warnings} warnings")
    
    print("\n✨ This validation system implements Linus's 'good taste':")
    print("• Eliminates special cases by using unified framework")
    print("• Catches real production issues (NaN, Inf, dimension mismatch)")
    print("• Provides clear, actionable error messages")
    print("• Fails fast to prevent downstream mysterious crashes")


if __name__ == "__main__":
    import asyncio
    print("🚀 Starting Linus-style data validation demo...")
    asyncio.run(demo_validation())
    print("\n🎉 Demo completed!")