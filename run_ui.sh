#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
exec conda run --no-capture-output -n ecotype python "$ROOT/ui/ui.py"
