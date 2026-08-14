#!/usr/bin/env bash
# The one definition of "lint" for this repository.
#
# It used to live in two places that disagreed: the CI lint job checked
# chatspatial and reproducibility/scripts, while the release gate checked
# chatspatial, tests and scripts. Neither covered everything, so a formatting
# error in tests/ passed the lint job and one in reproducibility/scripts/
# passed the gate. Both now call this.
set -euo pipefail

cd "$(dirname "$0")/../.."

TARGETS=(chatspatial tests scripts reproducibility/scripts)

echo "[lint] Checking formatting"
black --check "${TARGETS[@]}"
isort --check-only "${TARGETS[@]}"

echo "[lint] Running Ruff"
ruff check "${TARGETS[@]}"
# Complexity is enforced on the package only: tests and scripts are allowed to
# be long and linear where that reads better than a decomposition.
ruff check chatspatial --select C901

echo "[lint] PASS"
