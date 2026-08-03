"""Shared model foundations for public ChatSpatial parameters."""

from pydantic import BaseModel, ConfigDict


class FiniteModel(BaseModel):
    """Base class for exact, standards-compliant public JSON models."""

    model_config = ConfigDict(allow_inf_nan=False, extra="forbid")


class StrictParameters(FiniteModel):
    """Base class that rejects unsupported public analysis parameters."""

    model_config = ConfigDict(extra="forbid")
