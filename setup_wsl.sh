#!/bin/bash
# One-shot setup for Ubuntu on WSL2 (also works on plain Ubuntu/Debian).
# Installs build tools + Miniforge (if conda is missing), creates the `ecotype` env,
# builds the EcoSim helpers, and downloads the light Bakta database to ~/bakta_db.
# Usage: bash setup_wsl.sh [light|full]   (default: light)
set -e
cd "$(dirname "$0")"
db_type="${1:-light}"

sudo apt-get update
sudo apt-get install -y build-essential gfortran git curl

if ! command -v conda >/dev/null; then
    curl -fsSLo /tmp/miniforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
    bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
    "$HOME/miniforge3/bin/conda" init bash
fi
source "$(conda info --base 2>/dev/null || echo "$HOME/miniforge3")/etc/profile.d/conda.sh"

conda env list | grep -q '^ecotype ' || conda env create -f environment.yml
bash tools/build_ecosim_linux.sh
conda run -n ecotype bakta_db download --output "$HOME/bakta_db" --type "$db_type"

echo
echo "Done. Open a new terminal, then:  cd $(pwd) && conda activate ecotype && python ui/ui.py"
echo "Bakta database: ~/bakta_db/db-light (or ~/bakta_db/db for full)"
