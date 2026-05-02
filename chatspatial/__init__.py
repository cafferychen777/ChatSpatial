"""
ChatSpatial

Agentic workflow orchestration platform for spatial transcriptomics analysis.
Integrates 60 methods from Python and R ecosystems via Model Context Protocol.
"""

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import tomllib


try:
    __version__ = version("chatspatial")
except PackageNotFoundError:
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    with pyproject_path.open("rb") as f:
        __version__ = tomllib.load(f)["project"]["version"]


# Import configuration to set up environment
from . import config as config  # noqa: F401
