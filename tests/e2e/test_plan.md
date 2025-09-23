# ChatSpatial MCP Server - End-to-End Test Plan

## Test Objectives
Verify the functional completeness, user experience, and client integration of ChatSpatial MCP server in real environments.

## Test Environment
- **MCP Inspector**: For direct protocol testing
- **Claude Desktop**: For real user scenario testing
- **Test Data**: Standardized spatial transcriptomics datasets

---

## Test Case Checklist

### E2E-001: Basic Workflow Testing
**Objective**: Verify the completeness of standard analysis workflow

| Step | Operation | Expected Result | Status |
|------|-----------|-----------------|--------|
| 1 | Launch MCP Inspector | Server connects successfully, tool list displays | ⬜ |
| 2 | Call `load_data` to load visium data | Returns success message and data ID | ⬜ |
| 3 | Call `preprocess_data` for preprocessing | Returns preprocessing completion message | ⬜ |
| 4 | Call `visualize_data` to create spatial plot | Returns valid PNG image | ⬜ |
| 5 | Verify image quality | Image is clear, contains spatial coordinates and gene expression | ⬜ |

**Success Criteria**: All steps execute successfully, final image quality meets expectations

---

### E2E-002: Claude Desktop Integration Testing
**Objective**: Verify seamless integration with Claude Desktop

| Step | Operation | Expected Result | Status |
|------|-----------|-----------------|--------|
| 1 | Configure Claude Desktop MCP settings | Tool icons display in interface | ⬜ |
| 2 | Restart Claude Desktop | Successfully loads ChatSpatial tools | ⬜ |
| 3 | Send message: "Load visium data and show spatial expression plot for Gad1 gene" | Claude understands request and calls appropriate tools | ⬜ |
| 4 | Verify tool call sequence | Correctly calls load_data → visualize_data | ⬜ |
| 5 | Check final result | Chat window displays Gad1 gene spatial expression plot | ⬜ |

**Success Criteria**: Claude can understand natural language requests and correctly call tool chains

---

### E2E-003: Progress Reporting Testing
**Objective**: Verify progress reporting functionality for long-running tasks

| Step | Operation | Expected Result | Status |
|------|-----------|-----------------|--------|
| 1 | Launch MCP Inspector | Server connects successfully | ⬜ |
| 2 | Call time-consuming tool (e.g., `identify_spatial_domains`) with progressToken | Task starts execution | ⬜ |
| 3 | Observe Notifications panel | Real-time progress updates display | ⬜ |
| 4 | Wait for task completion | Receive completion notification and final result | ⬜ |
| 5 | Verify result correctness | Spatial domain identification results meet expectations | ⬜ |

**Success Criteria**: Progress reporting is accurate, user experience is good

---

### E2E-004: Error Handling Testing
**Objective**: Verify handling of various error scenarios

| Test Scenario | Operation | Expected Result | Status |
|---------------|-----------|-----------------|--------|
| Invalid Data ID | Request non-existent dataset | Return clear error message | ⬜ |
| Invalid Gene Name | Request visualization of non-existent gene | Graceful handling, provide alternatives or warnings | ⬜ |
| Missing Data | Request UMAP on unprocessed data | Auto-preprocess or provide clear guidance | ⬜ |
| Parameter Error | Provide invalid plot_type | Return parameter validation error | ⬜ |
| Memory Insufficient | Process oversized dataset | Graceful degradation or memory management | ⬜ |

**Success Criteria**: All error scenarios have appropriate error messages and recovery suggestions

---

### E2E-005: Enhanced Features Testing
**Objective**: Verify newly added enhanced features

| Feature | Test Operation | Expected Result | Status |
|---------|----------------|-----------------|--------|
| Spatial Contour Overlay | Request spatial plot with contours | Image displays cluster boundary contours | ⬜ |
| UMAP Dual Encoding | Request color+size encoded UMAP | Points have meaningful color and size | ⬜ |
| Spatial Interaction Visualization | Request ligand-receptor pair visualization | Display spatial distribution of LR pairs | ⬜ |
| Enhanced Heatmap | Request annotated heatmap | Display row/column annotation information | ⬜ |
| Integration Assessment | Request batch effect evaluation | Multi-panel display of integration quality | ⬜ |
| Network Visualization | Request neighborhood enrichment network plot | Display network layout | ⬜ |

**Success Criteria**: All new features work properly, image quality meets professional standards

---

### E2E-006: Performance Testing
**Objective**: Verify system performance under different loads

| Test Scenario | Data Scale | Expected Response Time | Status |
|---------------|------------|------------------------|--------|
| Small Dataset | 100 cells | < 5 seconds | ⬜ |
| Standard Dataset | 1000 cells | < 15 seconds | ⬜ |
| Large Dataset | 5000 cells | < 60 seconds | ⬜ |
| Concurrent Requests | 5 simultaneous requests | Normal processing, no crashes | ⬜ |
| Memory Usage | 10 consecutive visualizations | Memory growth controlled within reasonable range | ⬜ |

**Success Criteria**: Response times meet user expectations, system runs stably

---

### E2E-007: Compatibility Testing
**Objective**: Verify compatibility across different environments

| Environment | Test Content | Expected Result | Status |
|-------------|--------------|-----------------|--------|
| macOS | Complete functionality test | All features work normally | ⬜ |
| Linux | Complete functionality test | All features work normally | ⬜ |
| Windows | Complete functionality test | All features work normally | ⬜ |
| Python 3.10 | Basic functionality test | Runs normally | ⬜ |
| Python 3.11 | Basic functionality test | Runs normally | ⬜ |
| Python 3.12 | Basic functionality test | Runs normally | ⬜ |

**Success Criteria**: Stable operation on mainstream operating systems and Python versions

---

## Execution Guidelines

### Pre-test Preparation
1. **Environment Setup**
   ```bash
   # Install test dependencies
   pip install pytest pytest-asyncio pytest-cov

   # Install MCP Inspector
   npm install -g @modelcontextprotocol/inspector

   # Prepare test data
   python scripts/prepare_test_data.py
   ```

2. **Server Startup**
   ```bash
   # Launch MCP Inspector
   npx @modelcontextprotocol/inspector python -m chatspatial

   # Or start server for Claude Desktop
   python -m chatspatial
   ```

### Test Execution
1. **Manual Testing**: Execute test cases one by one according to the checklist
2. **Automated Testing**: Run end-to-end automation scripts
3. **Record Results**: Mark ✅ (pass) or ❌ (fail) in status column
4. **Issue Recording**: Detailed recording of discovered issues and solutions

### Test Reporting
After each test completion, generate a test report containing:
- Test execution date and environment information
- Detailed results of each test case
- Discovered issues and resolution status
- Performance data and trend analysis
- Improvement recommendations

### Test Frequency
- **Complete E2E Testing**: Before each major version release
- **Core Functionality Testing**: After each feature update
- **Regression Testing**: Weekly regular execution
- **Performance Testing**: Monthly execution

---

## Success Criteria

### Overall Success Criteria
- ✅ 95%+ test cases pass
- ✅ 100% success for critical workflows
- ✅ Error handling covers all known scenarios
- ✅ Performance metrics meet user expectations
- ✅ User experience is smooth and natural

### Quality Gates
- 🚫 Serious functional defects found → Block release
- ⚠️  Minor issues found → Record and schedule fixes
- 📊 Performance below standard → Optimize and retest
- 🐛 New bugs discovered → Add corresponding test cases

Through this comprehensive E2E test plan, ensure that ChatSpatial MCP server provides high-quality, reliable service experience across various real-world scenarios.