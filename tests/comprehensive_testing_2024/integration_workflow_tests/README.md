# ChatSpatial Integration Workflow Tests

This directory contains comprehensive end-to-end integration tests for ChatSpatial's spatial transcriptomics analysis workflows. These tests validate the complete system integration from data loading through analysis to visualization.

## Test Structure

### Core Test Modules

#### 1. `test_complete_analysis_workflow.py`
**Complete Spatial Transcriptomics Analysis Pipeline Testing**

Tests the full end-to-end analysis pipeline including:
- Basic Visium workflow (preprocessing → spatial analysis → visualization)
- seqFISH workflow with fallback mechanisms  
- Synthetic data complete pipeline with performance metrics
- Error recovery and diagnostic capabilities

**Key Validations:**
- End-to-end execution success rate
- Data shape consistency through pipeline
- Memory usage and execution time monitoring
- Graceful error handling

#### 2. `test_multi_tool_chaining.py`
**Multi-Tool Chain Validation**

Tests complex workflows chaining multiple ChatSpatial tools:
- Preprocessing → Spatial Analysis → Visualization chains
- Differential Analysis → Enrichment → Visualization chains  
- Complex spatial domain analysis chains
- Data flow validation between tools

**Key Validations:**
- Seamless tool integration
- Data consistency across tool boundaries
- Chain success rates and performance
- Fallback mechanism effectiveness

#### 3. `test_batch_processing_workflow.py`
**Batch Processing and Scalability Testing**

Tests batch processing capabilities across multiple datasets:
- Sequential batch processing workflows
- Parallel processing with multiple workers
- Spatial-specific batch workflows
- Error handling and recovery in batch mode

**Key Validations:**
- Batch processing consistency
- Parallel processing efficiency
- Resource usage optimization
- Error isolation and recovery

#### 4. `test_real_user_scenarios.py`
**Real-World Usage Pattern Simulation**

Simulates realistic user scenarios and workflows:
- New user first-time analysis experience
- Experienced user comprehensive analysis patterns
- Comparative analysis across datasets
- Troubleshooting and recovery scenarios
- Mixed user session patterns

**Key Validations:**
- User experience quality
- Workflow success rates for different user types
- Recovery from common user errors
- Session continuity and data preservation

#### 5. `test_data_flow_validation.py`
**Data Integrity and Flow Validation**

Tests data integrity and consistency across complex workflows:
- Preprocessing pipeline data integrity
- Spatial analysis data flow consistency
- Cross-tool data compatibility
- Data preservation through complex workflows

**Key Validations:**
- Data checksums and integrity verification
- Metadata preservation across steps
- Shape and content consistency
- Essential data preservation

### Support Infrastructure

#### `workflow_test_base.py`
**Core Testing Infrastructure**

Provides fundamental testing framework components:
- `WorkflowTestBase`: Base class for workflow testing with data loading and environment management
- `DataFlowValidator`: Tracks and validates data flow between workflow steps  
- `WorkflowExceptionRecovery`: Handles exception recovery with fallback mechanisms
- Resource usage measurement and performance monitoring

#### `run_all_workflow_tests.py`
**Comprehensive Test Runner**

Automated test execution and reporting:
- Runs all workflow test modules sequentially
- Generates comprehensive HTML reports
- Provides JSON and CSV result exports
- Calculates success rates and performance metrics
- Automated report generation with timestamps

## Usage

### Running Individual Tests

```bash
# Run a specific test module
pytest test_complete_analysis_workflow.py -v

# Run with detailed output
pytest test_multi_tool_chaining.py -v -s

# Run specific test within module
pytest test_batch_processing_workflow.py::TestBatchProcessingWorkflow::test_basic_batch_processing -v
```

### Running Full Test Suite

```bash
# Run all workflow tests with comprehensive reporting
python run_all_workflow_tests.py

# Alternative: Run via pytest
pytest . -v
```

### Test Data Requirements

The tests use datasets from:
- `/Users/apple/Research/SpatialTrans_MCP/chatspatial/tests/comprehensive_testing_2024/datasets/`

Required datasets:
- `medium_synthetic.h5ad` - Main testing dataset
- `benchmark_1kx2k.h5ad` - Large-scale testing
- `squidpy_visium.h5ad` - Visium-format testing
- `squidpy_seqfish.h5ad` - seqFISH-format testing
- `small_synthetic.h5ad` - Quick testing
- `high_sparsity.h5ad` - Edge case testing
- `empty_dataset.h5ad` - Error handling testing

## Test Categories and Coverage

### 1. End-to-End Execution
- **Coverage**: Complete analysis pipelines from raw data to final visualization
- **Validation**: Success rates, execution times, resource usage
- **Scenarios**: Visium, seqFISH, synthetic data, real datasets

### 2. Tool Integration
- **Coverage**: All major ChatSpatial tool combinations
- **Validation**: Data consistency, interface compatibility, error propagation
- **Scenarios**: Preprocessing chains, analysis chains, visualization chains

### 3. Batch Processing
- **Coverage**: Multiple dataset processing, parallel execution
- **Validation**: Consistency across datasets, performance scaling, error isolation
- **Scenarios**: Sequential processing, parallel processing, mixed dataset types

### 4. User Experience
- **Coverage**: Realistic usage patterns, error recovery, session continuity
- **Validation**: User journey success, error handling quality, recovery mechanisms
- **Scenarios**: New users, experienced users, troubleshooting sessions

### 5. Data Integrity  
- **Coverage**: Data preservation, transformation validation, consistency checking
- **Validation**: Checksums, metadata preservation, shape consistency
- **Scenarios**: Complex workflows, tool chaining, data preservation

## Performance Benchmarks

### Success Rate Targets
- **Complete Analysis Workflows**: ≥ 80% success rate
- **Multi-Tool Chaining**: ≥ 70% success rate  
- **Batch Processing**: ≥ 75% success rate
- **User Scenarios**: ≥ 60% success rate
- **Data Flow Validation**: ≥ 85% success rate

### Performance Targets
- **Single Dataset Analysis**: < 300 seconds
- **Batch Processing**: < 600 seconds total
- **Memory Usage**: < 4GB peak for medium datasets
- **Parallel Speedup**: > 1.2x for 2+ datasets

### Error Recovery Targets
- **Fallback Success Rate**: ≥ 70%
- **Error Isolation**: 100% (errors shouldn't crash entire workflow)
- **Recovery Time**: < 30 seconds additional overhead

## Reporting and Analysis

### Automated Reports
- **HTML Report**: Comprehensive visual report with test results, timing, and success rates
- **JSON Export**: Machine-readable results for CI/CD integration
- **CSV Summary**: Tabular results for analysis and tracking
- **JUnit XML**: Integration with testing frameworks and CI systems

### Report Contents
- Overall success rates and performance metrics
- Individual test module results and timing
- Error analysis and common failure patterns  
- Performance benchmarks and resource usage
- Data flow validation results
- Trend analysis (when run repeatedly)

### Report Locations
- `../reports/workflow_test_report_[timestamp].html` - Timestamped reports
- `../reports/latest_workflow_report.html` - Latest results (always current)
- `../reports/workflow_test_results_[timestamp].json` - Full results data
- `../reports/latest_workflow_results.json` - Latest results data

## Maintenance and Development

### Adding New Tests
1. Follow the existing test patterns in workflow test modules
2. Use `WorkflowTestBase` for consistent setup/teardown
3. Implement proper data flow validation using `DataFlowValidator`
4. Add comprehensive assertions with meaningful error messages
5. Update this README with new test descriptions

### Test Data Management
- Test datasets are managed centrally in the `datasets/` directory
- Use `load_test_dataset()` method for consistent data loading
- Add new datasets to `datasets_summary.csv` for tracking
- Ensure test datasets represent realistic use cases

### Performance Monitoring
- All tests include timing and resource usage measurement
- Performance regressions are flagged in reports  
- Memory usage is monitored to prevent resource leaks
- Execution times are tracked for performance optimization

### Error Handling Standards
- All tests implement graceful error handling
- Fallback mechanisms are tested explicitly
- Error messages provide actionable information
- Recovery scenarios are validated comprehensively

## Integration with CI/CD

The test suite is designed for integration with continuous integration systems:
- Exit codes indicate overall success/failure
- JUnit XML output for test result parsing
- JSON results for automated analysis
- Performance metrics for regression detection
- Comprehensive logging for debugging

Example CI integration:
```bash
# Run tests and capture results
python run_all_workflow_tests.py

# Check exit code for pass/fail
if [ $? -eq 0 ]; then
    echo "All workflow tests passed"
else
    echo "Workflow tests failed"
    exit 1
fi
```