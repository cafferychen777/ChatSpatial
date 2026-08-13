#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

distribution_dir="${1:-dist}"
if [[ -e "$distribution_dir" ]]; then
  echo "[release] Refusing to reuse existing output: $distribution_dir" >&2
  exit 1
fi

build_constraints="$PWD/constraints/release-build.txt"
source_date_epoch="${SOURCE_DATE_EPOCH:-$(git show -s --format=%ct HEAD)}"
export PIP_CONSTRAINT="$build_constraints"
export SOURCE_DATE_EPOCH="$source_date_epoch"

scripts/quality/check_test_gates.sh

echo "[release] Building distributions"
python -m build --outdir "$distribution_dir"

sdists=("$distribution_dir"/*.tar.gz)
if [[ "${#sdists[@]}" -ne 1 || ! -f "${sdists[0]}" ]]; then
  echo "[release] Expected exactly one source distribution" >&2
  exit 1
fi
python scripts/quality/normalize_sdist.py "${sdists[0]}" "$source_date_epoch"

echo "[release] Checking package metadata"
twine check "$distribution_dir"/*

echo "[release] Auditing package contents"
python scripts/quality/audit_distribution.py "$distribution_dir"

wheel_paths=("$distribution_dir"/*.whl)
if [[ "${#wheel_paths[@]}" -ne 1 || ! -f "${wheel_paths[0]}" ]]; then
  echo "[release] Expected exactly one wheel" >&2
  exit 1
fi

echo "[release] Testing the wheel in an isolated environment"
wheel_smoke_dir=$(mktemp -d "${TMPDIR:-/tmp}/chatspatial-wheel-smoke.XXXXXX")
trap 'rm -rf "$wheel_smoke_dir"' EXIT
uv venv --python 3.11 "$wheel_smoke_dir"
uv pip install --python "$wheel_smoke_dir/bin/python" "${wheel_paths[0]}"
"$wheel_smoke_dir/bin/python" scripts/quality/smoke_test_installed_wheel.py
rm -rf "$wheel_smoke_dir"
trap - EXIT
