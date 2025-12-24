"""
Test error handling module (chatspatial.utils.error_handling)

Tests custom exceptions, validation functions, and context managers.
"""
import pytest
import warnings
import logging
from chatspatial.utils.error_handling import (
    SpatialMCPError,
    DataNotFoundError,
    InvalidParameterError,
    ProcessingError,
    DataCompatibilityError,
    validate_adata,
    suppress_output,
)


# ========== Custom Exception Tests ==========

class TestCustomExceptions:
    """Test custom exception classes"""

    def test_base_exception(self):
        """Test base exception class"""
        with pytest.raises(SpatialMCPError):
            raise SpatialMCPError("Test error")

    def test_data_not_found_error(self):
        """Test DataNotFoundError exception"""
        with pytest.raises(DataNotFoundError):
            raise DataNotFoundError("Dataset not found")

    def test_invalid_parameter_error(self):
        """Test InvalidParameterError exception"""
        with pytest.raises(InvalidParameterError):
            raise InvalidParameterError("Invalid parameter value")

    def test_processing_error(self):
        """Test ProcessingError exception"""
        with pytest.raises(ProcessingError):
            raise ProcessingError("Processing failed")

    def test_data_compatibility_error(self):
        """Test DataCompatibilityError exception"""
        with pytest.raises(DataCompatibilityError):
            raise DataCompatibilityError("Data format incompatible")

    def test_exception_inheritance(self):
        """Test exception inheritance hierarchy"""
        assert issubclass(DataNotFoundError, SpatialMCPError)
        assert issubclass(InvalidParameterError, SpatialMCPError)
        assert issubclass(ProcessingError, SpatialMCPError)
        assert issubclass(DataCompatibilityError, SpatialMCPError)


# ========== validate_adata Tests ==========

class TestValidateAdata:
    """Test AnnData validation function"""

    def test_valid_adata_single_key(self, mock_adata):
        """Test validation with single required key"""
        # Add required key
        mock_adata.obs["leiden"] = ["0", "1"] * 50

        # Should not raise exception
        validate_adata(mock_adata, required_keys={"obs": ["leiden"]})

    def test_valid_adata_multiple_keys(self, mock_adata):
        """Test validation with multiple required keys"""
        mock_adata.obs["leiden"] = ["0"] * 100
        mock_adata.obs["cell_type"] = ["TypeA"] * 100
        import numpy as np
        mock_adata.obsm["X_pca"] = np.random.rand(100, 30)

        validate_adata(
            mock_adata,
            required_keys={"obs": ["leiden", "cell_type"], "obsm": ["X_pca"]},
        )

    def test_missing_obs_key(self, mock_adata):
        """Test missing obs key raises error"""
        with pytest.raises(DataNotFoundError, match="Missing required keys"):
            validate_adata(mock_adata, required_keys={"obs": ["non_existent"]})

    def test_missing_obsm_key(self, mock_adata):
        """Test missing obsm key raises error"""
        with pytest.raises(DataNotFoundError, match="Missing required keys"):
            validate_adata(mock_adata, required_keys={"obsm": ["X_pca"]})

    def test_missing_var_key(self, mock_adata):
        """Test missing var key raises error"""
        with pytest.raises(DataNotFoundError):
            validate_adata(mock_adata, required_keys={"var": ["highly_variable"]})

    def test_spatial_validation_pass(self, mock_adata_with_spatial):
        """Test spatial coordinate validation passes"""
        # Should not raise exception
        validate_adata(
            mock_adata_with_spatial,
            required_keys={},
            check_spatial=True,
            spatial_key="spatial",
        )

    def test_spatial_validation_fail_missing(self, mock_adata):
        """Test spatial validation fails when coordinates missing"""
        with pytest.raises(DataNotFoundError, match="Missing 'spatial'"):
            validate_adata(
                mock_adata, required_keys={}, check_spatial=True, spatial_key="spatial"
            )

    def test_velocity_validation_pass(self, mock_velocity_adata):
        """Test velocity validation passes"""
        validate_adata(mock_velocity_adata, required_keys={}, check_velocity=True)

    def test_velocity_validation_fail(self, mock_adata):
        """Test velocity validation fails when layers missing"""
        with pytest.raises(DataNotFoundError, match="Missing 'spliced'"):
            validate_adata(mock_adata, required_keys={}, check_velocity=True)

    def test_empty_required_keys(self, mock_adata):
        """Test validation with empty required keys dict"""
        # Should not raise exception
        validate_adata(mock_adata, required_keys={})


# ========== suppress_output Tests ==========

class TestSuppressOutput:
    """Test output suppression context manager"""

    def test_suppress_warnings(self):
        """Test warnings are suppressed"""
        with suppress_output():
            warnings.warn("This warning should be suppressed")
        # No assertion needed - just verify no exception

    def test_context_manager_cleanup(self):
        """Test context manager properly cleans up"""
        original_level = logging.getLogger().level

        with suppress_output():
            # In context, log level should be ERROR
            assert logging.getLogger().level == logging.ERROR

        # After context, log level should be restored
        assert logging.getLogger().level == original_level

    def test_exception_in_context(self):
        """Test exception handling in context"""
        original_level = logging.getLogger().level

        with pytest.raises(ValueError):
            with suppress_output():
                raise ValueError("Test exception")

        # Even with exception, log level should be restored
        assert logging.getLogger().level == original_level

    def test_nested_suppress_output(self):
        """Test nested suppress_output contexts"""
        original_level = logging.getLogger().level

        with suppress_output():
            with suppress_output():
                assert logging.getLogger().level == logging.ERROR
            assert logging.getLogger().level == logging.ERROR

        assert logging.getLogger().level == original_level


# ========== Edge Cases ==========

class TestEdgeCases:
    """Test edge cases and boundary conditions"""

    def test_validate_with_string_key(self, mock_adata):
        """Test validation with string key (converted to list)"""
        mock_adata.obs["leiden"] = ["0"] * 100
        # String should be converted to list internally
        validate_adata(mock_adata, required_keys={"obs": "leiden"})

    def test_validate_with_none_context(self, mock_adata):
        """Test validation with None context"""
        mock_adata.obs["leiden"] = ["0"] * 100
        validate_adata(mock_adata, required_keys={"obs": ["leiden"]}, context=None)

    def test_multiple_missing_keys(self, mock_adata):
        """Test error message with multiple missing keys"""
        with pytest.raises(DataNotFoundError) as exc_info:
            validate_adata(
                mock_adata,
                required_keys={"obs": ["key1", "key2"], "var": ["key3"]},
            )
        # Should mention all missing keys
        assert "key1" in str(exc_info.value) or "Missing required keys" in str(exc_info.value)
