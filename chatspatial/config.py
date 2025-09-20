"""
Configuration module for ChatSpatial
"""

import os
import warnings

# Configure Dask to use new DataFrame implementation
os.environ.setdefault("DASK_DATAFRAME__QUERY_PLANNING", "True")

# Try to import and configure dask
try:
    import dask

    dask.config.set({"dataframe.query-planning": True})
except ImportError:
    pass  # Dask not installed


# Suppress specific warnings
def configure_warnings():
    """Configure warning filters for known issues"""
    # Suppress dask legacy dataframe warning
    warnings.filterwarnings(
        "ignore",
        message="The legacy Dask DataFrame implementation is deprecated",
        category=FutureWarning,
    )

    # Suppress spatialdata functools.partial warnings
    warnings.filterwarnings(
        "ignore",
        message="functools.partial will be a method descriptor",
        category=FutureWarning,
    )

    # Suppress numba nopython warning
    warnings.filterwarnings(
        "ignore",
        message="nopython is set for njit and is ignored",
        category=RuntimeWarning,
    )

    # Suppress anndata read_text warning
    warnings.filterwarnings(
        "ignore",
        message="Importing read_text from `anndata` is deprecated",
        category=FutureWarning,
    )


# Configure on import
configure_warnings()
