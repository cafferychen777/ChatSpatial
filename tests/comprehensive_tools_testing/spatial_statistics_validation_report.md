# Comprehensive Spatial Statistics Validation Report
Generated on: 2025-08-24 04:14:24.847573

## Executive Summary
- **Overall Success Rate**: 89.5% (17/19 tests passed)

## Test Datasets
- **perfect_autocorr**: 100 spots, 5 genes - Data with perfect positive spatial autocorrelation (smooth gradients)
- **random**: 100 spots, 5 genes - Randomly distributed data with no spatial structure
- **checkerboard**: 100 spots, 3 genes - Checkerboard pattern with negative spatial autocorrelation
- **hotspots**: 200 spots, 4 genes - Data with distinct spatial hotspots
- **gradient**: 150 spots, 4 genes - Linear and non-linear gradient patterns
- **multiscale**: 300 spots, 3 genes - Multi-scale spatial patterns (large, medium, small)

## Moran's I Validation Results
### perfect_autocorr
- Mean Moran's I: 0.903
- Range: [0.898, 0.908]
- Significant genes: 5/5
- Validations:
  - expected_I_correct: ✓ PASS
  - high_autocorr_detected: ✓ PASS

### random
- Mean Moran's I: -0.026
- Range: [-0.079, 0.034]
- Significant genes: 0/5
- Validations:
  - expected_I_correct: ✓ PASS
  - random_near_expected: ✓ PASS

### checkerboard
- Mean Moran's I: 0.059
- Range: [0.054, 0.065]
- Significant genes: 0/3
- Validations:
  - expected_I_correct: ✓ PASS
  - negative_autocorr_detected: ✗ FAIL

### hotspots
- Mean Moran's I: 0.676
- Range: [0.536, 0.769]
- Significant genes: 4/4
- Validations:
  - expected_I_correct: ✓ PASS

### gradient
- Mean Moran's I: 0.804
- Range: [0.562, 0.919]
- Significant genes: 4/4
- Validations:
  - expected_I_correct: ✓ PASS

### multiscale
- Mean Moran's I: 0.460
- Range: [0.049, 0.845]
- Significant genes: 2/3
- Validations:
  - expected_I_correct: ✓ PASS

## Geary's C Validation Results
### perfect_autocorr
- Mean Geary's C: 0.067
- Range: [0.065, 0.068]
- Significant genes: 3

### random
- Mean Geary's C: 1.013
- Range: [0.939, 1.060]
- Significant genes: 0

### checkerboard
- Mean Geary's C: 0.929
- Range: [0.922, 0.933]
- Significant genes: 0

### hotspots
- Mean Geary's C: 0.355
- Range: [0.293, 0.450]
- Significant genes: 3

### gradient
- Mean Geary's C: 0.071
- Range: [0.041, 0.100]
- Significant genes: 3

### multiscale
- Mean Geary's C: 0.541
- Range: [0.165, 0.945]
- Significant genes: 3

## Performance Benchmark Results
### perfect_autocorr
- moran: 0.041s (0.008s per gene)
- geary: 0.071s (0.014s per gene)
- local_moran: 0.027s (0.005s per gene)
- getis_ord: 0.200s (0.040s per gene)

### random
- moran: 0.068s (0.014s per gene)
- geary: 0.092s (0.018s per gene)
- local_moran: 0.036s (0.007s per gene)
- getis_ord: 0.156s (0.031s per gene)

### checkerboard
- moran: 0.030s (0.010s per gene)
- geary: 0.080s (0.027s per gene)
- local_moran: 0.021s (0.007s per gene)
- getis_ord: 0.203s (0.068s per gene)

### hotspots
- moran: 0.071s (0.018s per gene)
- geary: 0.072s (0.018s per gene)
- local_moran: 0.029s (0.007s per gene)
- getis_ord: 0.217s (0.054s per gene)

### gradient
- moran: 0.055s (0.014s per gene)
- geary: 0.081s (0.020s per gene)
- local_moran: 0.029s (0.007s per gene)
- getis_ord: 0.275s (0.069s per gene)

### multiscale
- moran: 0.066s (0.022s per gene)
- geary: 0.324s (0.108s per gene)
- local_moran: 0.163s (0.054s per gene)
- getis_ord: 0.346s (0.115s per gene)

## Spatial Weights Matrix Validation
