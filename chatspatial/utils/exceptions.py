"""
Exception classes for ChatSpatial.

All exceptions in one place. No duplication.
"""


class ChatSpatialError(Exception):
    """Base exception for all ChatSpatial errors."""

    pass


class DataError(ChatSpatialError):
    """Data-related errors (missing, invalid format, etc.)."""

    pass


class DataNotFoundError(DataError):
    """Required data not found."""

    pass


class DataCompatibilityError(DataError):
    """Data format or compatibility issues."""

    pass


class ParameterError(ChatSpatialError):
    """Invalid parameter errors."""

    pass


class ProcessingError(ChatSpatialError):
    """Errors during analysis processing."""

    pass


class DependencyError(ChatSpatialError):
    """Missing or incompatible dependency."""

    pass
