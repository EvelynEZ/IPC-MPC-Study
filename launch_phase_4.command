#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

if [ ! -x ".venv/bin/jupyter" ]; then
  echo "The notebook environment is not installed."
  read -r -p "Press Enter to close..."
  exit 1
fi

exec .venv/bin/jupyter lab notebooks/phase_4_complications.ipynb
