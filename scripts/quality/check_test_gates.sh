#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/../.."

scripts/quality/check_lint.sh

echo "[gate] Running mypy"
mypy chatspatial

echo "[gate] Collecting tests"
collect_output=$(pytest --collect-only -q 2>&1 || true)
echo "$collect_output"

collected=$(echo "$collect_output" | sed -n 's/.*collected \([0-9][0-9]*\) item.*/\1/p' | tail -n1)
if [[ -z "${collected:-}" || "$collected" -eq 0 ]]; then
  echo "[gate] FAIL: pytest collected 0 tests"
  exit 1
fi

echo "[gate] Running default fast suite"
pytest -m "not slow"

echo "[gate] Verifying at least one e2e workflow passes"
pytest tests/e2e/test_core_workflows.py -q

echo "[gate] PASS: collected=$collected, e2e workflow validated"
