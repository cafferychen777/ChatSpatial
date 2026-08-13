#!/usr/bin/env python3
"""Exercise the default FlashS backend from an installed ChatSpatial wheel."""

from __future__ import annotations

import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np
from flashs import FlashS

import chatspatial


def main() -> None:
    """Run a deterministic, minimal FlashS spatial gene calculation."""
    package_path = Path(chatspatial.__file__).resolve()
    environment_path = Path(sys.prefix).resolve()
    if not package_path.is_relative_to(environment_path):
        raise RuntimeError(
            f"ChatSpatial was not imported from the test environment: {package_path}"
        )

    coordinates = np.asarray(
        [(row, column) for row in range(6) for column in range(6)], dtype=float
    )
    rng = np.random.default_rng(0)
    expression = rng.poisson(2.0, size=(36, 4)).astype(float) + 1.0
    gene_names = [f"gene_{index}" for index in range(expression.shape[1])]

    result = FlashS(
        n_features=16,
        n_scales=2,
        min_expressed=3,
        random_state=0,
    ).fit_test(coordinates, expression, gene_names=gene_names)

    if result.n_tested != expression.shape[1]:
        raise RuntimeError(
            f"FlashS tested {result.n_tested} genes; expected {expression.shape[1]}"
        )
    if len(result.pvalues) != expression.shape[1]:
        raise RuntimeError("FlashS returned an unexpected number of p-values")

    print(
        "[release] Installed wheel smoke test passed: "
        f"chatspatial {version('chatspatial')}, flashs {version('flashs')}, "
        f"{result.n_tested} genes tested"
    )


if __name__ == "__main__":
    main()
