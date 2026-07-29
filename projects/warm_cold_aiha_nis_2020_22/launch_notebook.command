#!/bin/bash
set -e
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$PROJECT_DIR/../.." && pwd)"
cd "$REPO_ROOT"
exec .venv/bin/jupyter lab projects/warm_cold_aiha_nis_2020_22/Warm_Cold_AIHA_NIS_2020_2022.ipynb
