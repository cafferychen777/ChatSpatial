"""
Test Pydantic data models (chatspatial.models.data)

These tests verify that parameter model validation works correctly.
"""
import pytest
from pydantic import ValidationError

from chatspatial.models.data import (
    PreprocessingParameters,
    AnnotationParameters,
    DeconvolutionParameters,
    SpatialDomainParameters,
    SpatialStatisticsParameters,
    CellCommunicationParameters,
    ColumnInfo,
    SpatialDataset,
)


# ========== ColumnInfo Tests ==========

class TestColumnInfo:
    """Test ColumnInfo model"""

    def test_categorical_column_valid(self):
        """Test valid categorical column"""
        col = ColumnInfo(
            name="cell_type",
            dtype="categorical",
            n_unique=5,
            sample_values=["TypeA", "TypeB", "TypeC"],
        )
        assert col.name == "cell_type"
        assert col.dtype == "categorical"
        assert col.n_unique == 5

    def test_numerical_column_valid(self):
        """Test valid numerical column"""
        col = ColumnInfo(
            name="total_counts", dtype="numerical", n_unique=100, range=(0.0, 1000.0)
        )
        assert col.name == "total_counts"
        assert col.dtype == "numerical"
        assert col.range == (0.0, 1000.0)

    def test_invalid_dtype(self):
        """Test invalid data type"""
        with pytest.raises(ValidationError):
            ColumnInfo(name="test", dtype="invalid_type", n_unique=10)


# ========== SpatialDataset Tests ==========

class TestSpatialDataset:
    """Test SpatialDataset model"""

    def test_basic_dataset(self):
        """Test basic dataset"""
        dataset = SpatialDataset(
            id="test_001",
            name="Test Dataset",
            data_type="10x_visium",
            n_cells=1000,
            n_genes=500,
        )
        assert dataset.id == "test_001"
        assert dataset.n_cells == 1000

    def test_invalid_data_type(self):
        """Test invalid data type"""
        with pytest.raises(ValidationError):
            SpatialDataset(
                id="test",
                name="Test",
                data_type="invalid_type",
                n_cells=100,
                n_genes=50,
            )


# ========== PreprocessingParameters Tests ==========

class TestPreprocessingParameters:
    """Test preprocessing parameters model"""

    def test_default_parameters(self):
        """Test default parameters"""
        params = PreprocessingParameters()
        assert params.normalization == "log"
        assert params.n_hvgs == 2000
        assert params.n_pcs == 30
        assert params.filter_genes_min_cells == 3

    def test_valid_normalization_methods(self):
        """Test all valid normalization methods"""
        valid_methods = ["log", "sct", "pearson_residuals", "none", "scvi"]
        for method in valid_methods:
            params = PreprocessingParameters(normalization=method)
            assert params.normalization == method

    def test_invalid_normalization_method(self):
        """Test invalid normalization method"""
        with pytest.raises(ValidationError):
            PreprocessingParameters(normalization="invalid_method")

    def test_n_hvgs_out_of_range(self):
        """Test n_hvgs out of range"""
        with pytest.raises(ValidationError):
            PreprocessingParameters(n_hvgs=10000)  # Exceeds limit of 5000

        with pytest.raises(ValidationError):
            PreprocessingParameters(n_hvgs=0)  # Must be > 0

    def test_n_pcs_out_of_range(self):
        """Test n_pcs out of range"""
        with pytest.raises(ValidationError):
            PreprocessingParameters(n_pcs=150)  # Exceeds limit of 100

        with pytest.raises(ValidationError):
            PreprocessingParameters(n_pcs=-1)  # Must be > 0

    def test_custom_parameters(self):
        """Test custom parameters"""
        params = PreprocessingParameters(
            normalization="sct",
            n_hvgs=1000,
            n_pcs=20,
            filter_genes_min_cells=5,
            filter_cells_min_genes=50,
        )
        assert params.normalization == "sct"
        assert params.n_hvgs == 1000
        assert params.filter_genes_min_cells == 5


# ========== AnnotationParameters Tests ==========

class TestAnnotationParameters:
    """Test cell annotation parameters model"""

    def test_valid_methods(self):
        """Test all valid annotation methods"""
        valid_methods = [
            "tangram",
            "scanvi",
            "cellassign",
            "mllmcelltype",  # Note: no underscore
            "sctype",
            "singler",
        ]
        for method in valid_methods:
            params = AnnotationParameters(method=method, reference_id="ref_001")
            assert params.method == method

    def test_invalid_method(self):
        """Test invalid method"""
        with pytest.raises(ValidationError):
            AnnotationParameters(method="invalid_method", reference_id="ref_001")


# ========== DeconvolutionParameters Tests ==========

class TestDeconvolutionParameters:
    """Test deconvolution parameters model"""

    def test_valid_methods(self):
        """Test all valid deconvolution methods"""
        valid_methods = [
            "cell2location",
            "destvi",
            "rctd",
            "stereoscope",
            "tangram",
            "spotlight",
        ]
        for method in valid_methods:
            params = DeconvolutionParameters(
                method=method,
                reference_data_id="ref_001",
                cell_type_key="cell_type"  # Required field
            )
            assert params.method == method


# ========== SpatialDomainParameters Tests ==========

class TestSpatialDomainParameters:
    """Test spatial domain identification parameters model"""

    def test_valid_methods(self):
        """Test all valid methods"""
        valid_methods = ["spagcn", "stagate", "graphst", "leiden"]
        for method in valid_methods:
            params = SpatialDomainParameters(method=method)
            assert params.method == method

    def test_n_domains_validation(self):
        """Test domain count validation"""
        # Valid range
        params = SpatialDomainParameters(method="spagcn", n_domains=5)
        assert params.n_domains == 5

        # Out of range
        with pytest.raises(ValidationError):
            SpatialDomainParameters(method="spagcn", n_domains=0)

        with pytest.raises(ValidationError):
            SpatialDomainParameters(method="spagcn", n_domains=51)


# ========== SpatialStatisticsParameters Tests ==========

class TestSpatialStatisticsParameters:
    """Test spatial statistics parameters model"""

    def test_valid_analysis_types(self):
        """Test all valid analysis types"""
        valid_types = [
            "moran",
            "local_moran",
            "geary",
            "getis_ord",
            "neighborhood",
            "co_occurrence",
            "ripley",
            "bivariate_moran",
            "join_count",
            "local_join_count",
            "network_properties",
            "centrality",
        ]
        for analysis_type in valid_types:
            params = SpatialStatisticsParameters(analysis_type=analysis_type)
            assert params.analysis_type == analysis_type

    def test_invalid_analysis_type(self):
        """Test invalid analysis type"""
        with pytest.raises(ValidationError):
            SpatialStatisticsParameters(analysis_type="invalid_stat_method")


# ========== CellCommunicationParameters Tests ==========

class TestCellCommunicationParameters:
    """Test cell communication parameters model"""

    def test_valid_methods(self):
        """Test all valid methods"""
        valid_methods = ["liana", "cellphonedb", "cellchat_r"]
        for method in valid_methods:
            params = CellCommunicationParameters(
                method=method,
                species="human",  # Required field
                cell_type_key="cell_type",  # Required field
                groupby="cell_type"
            )
            assert params.method == method

    def test_species_required(self):
        """Test species is required"""
        params = CellCommunicationParameters(
            method="liana",
            species="human",
            cell_type_key="cell_type",
            groupby="leiden"
        )
        assert params.species == "human"


# ========== Parametrized Tests ==========

@pytest.mark.parametrize(
    "normalization,expected",
    [
        ("log", "log"),
        ("sct", "sct"),
        ("pearson_residuals", "pearson_residuals"),
        ("none", "none"),
        ("scvi", "scvi"),
    ],
)
def test_all_normalization_methods(normalization, expected):
    """Parametrized test for all normalization methods"""
    params = PreprocessingParameters(normalization=normalization)
    assert params.normalization == expected


@pytest.mark.parametrize(
    "n_hvgs,should_pass",
    [
        (100, True),
        (2000, True),
        (5000, True),
        (0, False),  # Must be > 0
        (10000, False),  # Exceeds limit
    ],
)
def test_n_hvgs_validation(n_hvgs, should_pass):
    """Parametrized test for n_hvgs validation"""
    if should_pass:
        params = PreprocessingParameters(n_hvgs=n_hvgs)
        assert params.n_hvgs == n_hvgs
    else:
        with pytest.raises(ValidationError):
            PreprocessingParameters(n_hvgs=n_hvgs)
