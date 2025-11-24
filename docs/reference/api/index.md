# API Reference

Complete technical documentation for all ChatSpatial MCP tools and data models.

## API Documentation

- [Data Models](data_models.md) - Parameter schemas, data structures, and type definitions
- [Error Handling](error_handling.md) - Error codes and exception handling
- [Tool Reference](../../reference/quick-reference/all-tools.md) - Complete list of available tools

## Understanding ChatSpatial's API

ChatSpatial exposes its functionality through the Model Context Protocol (MCP), which allows LLMs like Claude to call analysis tools with structured parameters.

### Tool Call Structure

Every ChatSpatial tool follows this pattern:

```json
{
  "tool_name": "analyze_spatial_data",
  "parameters": {
    "data_id": "dataset_001",
    "method": "spagcn",
    "params": {
      "n_domains": 7,
      "resolution": 0.5
    }
  }
}
```

### Parameter Validation

All parameters are validated using Pydantic models to ensure:

- **Type safety** - Correct data types for all inputs
- **Range validation** - Parameters within acceptable bounds
- **Required fields** - All mandatory parameters provided
- **Default values** - Sensible defaults when not specified

### Return Values

Tools return structured results containing:

- **Analysis results** - Numerical data, statistics, classifications
- **Metadata** - Method used, parameters, execution time
- **Visualizations** - Optional images in base64 format
- **Messages** - Informative feedback about the analysis

## Data Model Categories

### Input Parameters

Models defining parameters accepted by tools:

- `PreprocessingParameters` - Data preprocessing options
- `AnnotationParameters` - Cell type annotation settings
- `SpatialDomainParameters` - Spatial domain identification
- `DeconvolutionParameters` - Cell type deconvolution
- `VisualizationParameters` - Plotting and visualization

### Analysis Results

Models defining analysis outputs:

- `AnalysisResult` - Generic analysis result container
- `AnnotationResult` - Cell type annotation results
- `SpatialStatisticsResult` - Spatial statistics outputs
- `CellCommunicationResult` - Cell-cell interaction results

### Data Types

Common data structures:

- `DataID` - Dataset identifier
- `GeneList` - List of gene names
- `CellTypeMapping` - Cell type annotations
- `SpatialCoordinates` - Spatial location data

## Error Handling

ChatSpatial provides detailed error messages for:

- **Invalid parameters** - Wrong types, out-of-range values
- **Missing dependencies** - Optional packages not installed
- **Data incompatibility** - Unsuitable data for chosen method
- **Runtime errors** - Analysis execution failures

See [Error Handling](error_handling.md) for details.

## Next Steps

- Browse [Data Models](data_models.md) for detailed schemas
- Check [All Tools](../../reference/quick-reference/all-tools.md) for available methods
- Try [Tutorials](../../tutorials/index.md) for hands-on learning
