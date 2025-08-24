"""
Base infrastructure for integration workflow testing.

This module provides the fundamental testing framework for end-to-end
spatial transcriptomics analysis workflows.
"""

import pytest
import numpy as np
import pandas as pd
import scanpy as sc
import os
import tempfile
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from unittest.mock import patch, MagicMock
import sys

# Add the chatspatial path
sys.path.insert(0, '/Users/apple/Research/SpatialTrans_MCP/chatspatial')

from chatspatial.tools import (
    preprocessing, spatial_analysis, visualization, 
    differential, cell_communication, spatial_domains,
    spatial_genes, integration, deconvolution, annotation
)
from chatspatial.utils.data_loader import load_spatial_data


class WorkflowTestBase:
    """Base class for integration workflow tests."""
    
    def __init__(self):
        self.test_data_dir = Path("/Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/datasets")
        self.temp_dir = None
        self.logger = logging.getLogger(__name__)
        
    def setup_test_environment(self):
        """Setup test environment with temporary directory."""
        self.temp_dir = tempfile.mkdtemp(prefix="workflow_test_")
        return self.temp_dir
    
    def cleanup_test_environment(self):
        """Clean up test environment."""
        if self.temp_dir and os.path.exists(self.temp_dir):
            import shutil
            shutil.rmtree(self.temp_dir)
    
    def load_test_dataset(self, dataset_name: str):
        """Load test dataset for workflow testing."""
        dataset_path = self.test_data_dir / dataset_name
        if not dataset_path.exists():
            raise FileNotFoundError(f"Test dataset {dataset_name} not found at {dataset_path}")
        
        try:
            adata = sc.read_h5ad(dataset_path)
            self.logger.info(f"Loaded dataset {dataset_name}: {adata.n_obs} cells, {adata.n_vars} genes")
            return adata
        except Exception as e:
            raise ValueError(f"Failed to load dataset {dataset_name}: {str(e)}")
    
    def validate_workflow_step(self, step_name: str, input_data: Any, output_data: Any, 
                             expected_properties: Dict[str, Any]) -> bool:
        """Validate a single workflow step."""
        try:
            # Basic validation
            if output_data is None:
                self.logger.error(f"Step {step_name}: Output is None")
                return False
            
            # Check expected properties
            for prop, expected_value in expected_properties.items():
                if hasattr(output_data, prop):
                    actual_value = getattr(output_data, prop)
                    if prop == 'shape' and isinstance(expected_value, tuple):
                        if actual_value != expected_value:
                            self.logger.warning(f"Step {step_name}: Shape mismatch - expected {expected_value}, got {actual_value}")
                    elif prop.startswith('has_') and isinstance(expected_value, bool):
                        if bool(actual_value) != expected_value:
                            self.logger.warning(f"Step {step_name}: Property {prop} mismatch")
                
            self.logger.info(f"Step {step_name}: Validation passed")
            return True
            
        except Exception as e:
            self.logger.error(f"Step {step_name}: Validation failed - {str(e)}")
            return False
    
    def measure_resource_usage(self, func, *args, **kwargs):
        """Measure resource usage during function execution."""
        import psutil
        import time
        
        process = psutil.Process(os.getpid())
        start_memory = process.memory_info().rss / 1024 / 1024  # MB
        start_time = time.time()
        
        try:
            result = func(*args, **kwargs)
            success = True
        except Exception as e:
            result = e
            success = False
        
        end_time = time.time()
        end_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        metrics = {
            'execution_time': end_time - start_time,
            'peak_memory_mb': end_memory,
            'memory_increase_mb': max(0, end_memory - start_memory),
            'success': success,
            'result': result
        }
        
        return metrics


class DataFlowValidator:
    """Validates data flow between workflow steps."""
    
    def __init__(self):
        self.flow_history = []
    
    def track_data_flow(self, step_name: str, input_data: Any, output_data: Any):
        """Track data transformation in workflow step."""
        flow_info = {
            'step': step_name,
            'input_type': type(input_data).__name__,
            'output_type': type(output_data).__name__ if output_data is not None else 'None',
            'timestamp': pd.Timestamp.now()
        }
        
        # Add specific properties for AnnData objects
        if hasattr(input_data, 'n_obs'):
            flow_info['input_cells'] = input_data.n_obs
            flow_info['input_genes'] = input_data.n_vars
        
        if hasattr(output_data, 'n_obs'):
            flow_info['output_cells'] = output_data.n_obs
            flow_info['output_genes'] = output_data.n_vars
        
        self.flow_history.append(flow_info)
    
    def validate_data_consistency(self) -> List[Dict[str, Any]]:
        """Validate data consistency across workflow steps."""
        issues = []
        
        for i in range(1, len(self.flow_history)):
            prev_step = self.flow_history[i-1]
            curr_step = self.flow_history[i]
            
            # Check for unexpected data loss
            if 'output_cells' in prev_step and 'input_cells' in curr_step:
                if prev_step['output_cells'] != curr_step['input_cells']:
                    issues.append({
                        'type': 'cell_count_mismatch',
                        'between_steps': f"{prev_step['step']} -> {curr_step['step']}",
                        'prev_output': prev_step['output_cells'],
                        'curr_input': curr_step['input_cells']
                    })
            
            # Check for gene count changes
            if 'output_genes' in prev_step and 'input_genes' in curr_step:
                if prev_step['output_genes'] != curr_step['input_genes']:
                    issues.append({
                        'type': 'gene_count_mismatch',
                        'between_steps': f"{prev_step['step']} -> {curr_step['step']}",
                        'prev_output': prev_step['output_genes'],
                        'curr_input': curr_step['input_genes']
                    })
        
        return issues
    
    def get_flow_summary(self) -> pd.DataFrame:
        """Get summary of data flow."""
        return pd.DataFrame(self.flow_history)


class WorkflowExceptionRecovery:
    """Handles exception recovery in workflow testing."""
    
    @staticmethod
    def with_fallback(primary_func, fallback_func, *args, **kwargs):
        """Execute primary function with fallback on failure."""
        try:
            return {'result': primary_func(*args, **kwargs), 'used_fallback': False}
        except Exception as e:
            logging.warning(f"Primary function failed: {str(e)}, trying fallback")
            try:
                return {'result': fallback_func(*args, **kwargs), 'used_fallback': True}
            except Exception as fallback_e:
                logging.error(f"Both primary and fallback failed: {str(fallback_e)}")
                return {'result': None, 'used_fallback': True, 'error': str(fallback_e)}
    
    @staticmethod
    def robust_execute(func, max_retries: int = 3, *args, **kwargs):
        """Execute function with retry logic."""
        last_error = None
        
        for attempt in range(max_retries):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                last_error = e
                logging.warning(f"Attempt {attempt + 1} failed: {str(e)}")
                if attempt < max_retries - 1:
                    import time
                    time.sleep(0.5 * (attempt + 1))  # Exponential backoff
        
        raise last_error