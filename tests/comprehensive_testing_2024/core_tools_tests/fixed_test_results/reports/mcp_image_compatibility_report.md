# MCP Image Transmission Compatibility Report

Generated: 2025-08-24 04:38:52

## Executive Summary

- **Images Tested**: 15
- **Total Tests**: 60
- **Passed Tests**: 30
- **Failed Tests**: 30
- **Success Rate**: 50.0%

🔴 **Overall Status**: POOR - Significant compatibility issues detected

## Detailed Results by Dataset

### Dataset: quick_test

#### Image: quick_test_umap_clusters.png

**base64_conversion**: ✅ SUCCESS
- Original size: 69037 bytes
- Base64 size: 92052 bytes
- Compression ratio: 1.33
- Duration: 0.001s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.234s

**file_size_limits**: ✅ SUCCESS
- File size: 0.07 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: quick_test_qc_violins.png

**base64_conversion**: ✅ SUCCESS
- Original size: 133699 bytes
- Base64 size: 178268 bytes
- Compression ratio: 1.33
- Duration: 0.001s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.061s

**file_size_limits**: ✅ SUCCESS
- File size: 0.13 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: quick_test_expression_heatmap.png

**base64_conversion**: ✅ SUCCESS
- Original size: 58908 bytes
- Base64 size: 78544 bytes
- Compression ratio: 1.33
- Duration: 0.001s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.079s

**file_size_limits**: ✅ SUCCESS
- File size: 0.06 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: quick_test_spatial_clusters.png

**base64_conversion**: ✅ SUCCESS
- Original size: 68166 bytes
- Base64 size: 90888 bytes
- Compression ratio: 1.33
- Duration: 0.001s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.059s

**file_size_limits**: ✅ SUCCESS
- File size: 0.07 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: quick_test_marker_dotplot.png

**base64_conversion**: ✅ SUCCESS
- Original size: 108051 bytes
- Base64 size: 144068 bytes
- Compression ratio: 1.33
- Duration: 0.001s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.052s

**file_size_limits**: ✅ SUCCESS
- File size: 0.10 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

### Dataset: standard

#### Image: standard_qc_violins.png

**base64_conversion**: ✅ SUCCESS
- Original size: 100502 bytes
- Base64 size: 134004 bytes
- Compression ratio: 1.33
- Duration: 0.001s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.057s

**file_size_limits**: ✅ SUCCESS
- File size: 0.10 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: standard_expression_heatmap.png

**base64_conversion**: ✅ SUCCESS
- Original size: 44532 bytes
- Base64 size: 59376 bytes
- Compression ratio: 1.33
- Duration: 0.001s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.067s

**file_size_limits**: ✅ SUCCESS
- File size: 0.04 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: standard_umap_clusters.png

**base64_conversion**: ✅ SUCCESS
- Original size: 382218 bytes
- Base64 size: 509624 bytes
- Compression ratio: 1.33
- Duration: 0.002s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.057s

**file_size_limits**: ✅ SUCCESS
- File size: 0.36 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: standard_spatial_clusters.png

**base64_conversion**: ✅ SUCCESS
- Original size: 714230 bytes
- Base64 size: 952308 bytes
- Compression ratio: 1.33
- Duration: 0.004s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.051s

**file_size_limits**: ✅ SUCCESS
- File size: 0.68 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: standard_marker_dotplot.png

**base64_conversion**: ✅ SUCCESS
- Original size: 107156 bytes
- Base64 size: 142876 bytes
- Compression ratio: 1.33
- Duration: 0.001s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.053s

**file_size_limits**: ✅ SUCCESS
- File size: 0.10 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

### Dataset: performance

#### Image: performance_spatial_clusters.png

**base64_conversion**: ✅ SUCCESS
- Original size: 258618 bytes
- Base64 size: 344824 bytes
- Compression ratio: 1.33
- Duration: 0.002s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.051s

**file_size_limits**: ✅ SUCCESS
- File size: 0.25 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: performance_qc_violins.png

**base64_conversion**: ✅ SUCCESS
- Original size: 57507 bytes
- Base64 size: 76676 bytes
- Compression ratio: 1.33
- Duration: 0.001s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.050s

**file_size_limits**: ✅ SUCCESS
- File size: 0.05 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: performance_umap_clusters.png

**base64_conversion**: ✅ SUCCESS
- Original size: 279309 bytes
- Base64 size: 372412 bytes
- Compression ratio: 1.33
- Duration: 0.002s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.053s

**file_size_limits**: ✅ SUCCESS
- File size: 0.27 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: performance_marker_dotplot.png

**base64_conversion**: ✅ SUCCESS
- Original size: 19487 bytes
- Base64 size: 25984 bytes
- Compression ratio: 1.33
- Duration: 0.000s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.053s

**file_size_limits**: ✅ SUCCESS
- File size: 0.02 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

#### Image: performance_expression_heatmap.png

**base64_conversion**: ✅ SUCCESS
- Original size: 53562 bytes
- Base64 size: 71416 bytes
- Compression ratio: 1.33
- Duration: 0.001s

**mcp_image_creation**: ❌ FAILED
- Error: Image.__init__() got an unexpected keyword argument 'media_type'
- Duration: 0.000s

**image_utils_compatibility**: ❌ FAILED
- Error: fig_to_image output validation failed
- Duration: 0.049s

**file_size_limits**: ✅ SUCCESS
- File size: 0.05 MB
- Category: small
- Recommendation: Optimal size for MCP transmission
- Duration: 0.000s

## Recommendations

⚠️ **30 test(s) failed.** Review the following:


## Technical Details

### MCP Image Protocol Requirements
- **Format**: PNG/JPEG images encoded as base64
- **Size limit**: Recommended < 5MB for optimal performance
- **Media type**: Proper MIME type specification required
- **Encoding**: UTF-8 base64 encoding

---
*Report generated by ChatSpatial MCP Image Compatibility Tester*
