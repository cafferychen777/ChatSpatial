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
