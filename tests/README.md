# ChatSpatial Test Suite

This directory contains all test scripts for ChatSpatial, organized by functionality.

## Directory Structure

```text
tests/
├── README.md                          # This document
├── visualization_tests/               # Visualization tests
│   ├── test_visualization_comprehensive.py
│   ├── test_visualization_scenarios.py
│   ├── test_direct_visualization.py
│   ├── test_final_visualization.py
│   ├── test_mcp_server_visualization.py
│   └── test_image_extraction.py
├── claude_tests/                      # Claude frontend interaction tests
│   ├── test_claude_conversation_simulation.py    # Original conversation simulation
│   ├── test_claude_conversation_improved.py      # Improved version (some failures)
│   ├── test_claude_conversation_fixed.py         # Fixed version (100% success)
│   ├── test_claude_frontend_complete.py
│   ├── test_claude_frontend_simulation.py
│   └── test_claude_functionality_simple.py
├── stress_tests/                      # Stress tests & real-world scenarios
│   ├── test_comprehensive_stress.py
│   ├── test_real_analysis_session.py
│   └── test_real_world_scenarios.py
├── test_all_features_claude.py        # Full feature tests
└── test_fixes_verification.py         # Fix verification tests
```

## Run Tests

### 1. Run all feature tests (recommended)

```bash
python tests/test_all_features_claude.py
```

### 2. Run Claude conversation simulation tests (most complete)

```bash
python tests/claude_tests/test_claude_conversation_fixed.py
```

### 3. Run visualization tests

```bash
python tests/visualization_tests/test_visualization_scenarios.py
```

### 4. Run stress tests

```bash
python tests/stress_tests/test_comprehensive_stress.py
```

## Test Results

Latest results:

- Success rate: 100% (all tests passed)
- Test report: `docs/test_reports/FINAL_TEST_REPORT.md`

## Dependencies

Before running tests, ensure required dependencies are installed:

```bash
# Base dependencies
pip install -r requirements.txt

# Full dependencies (including advanced features)
pip install -r requirements-full.txt
```

## Notes

1. Some tests generate temporary data and visualization files
2. Test runtime is about 2–5 minutes (depending on the scope)
3. If you encounter dependency issues, see `requirements-full.txt`