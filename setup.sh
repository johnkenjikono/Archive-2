#!/bin/bash
# Setup script for full pipeline dependencies

set -e
cd "$(dirname "$0")"  # paths below (venv, tools/) are repo-relative

echo "================================================================"
echo "Ecotype Pipeline Setup Tool"
echo "================================================================"
echo ""

# 0. Arch Linux system packages (compilers for EcoSim's tools, Java). Only prompts for sudo when something is missing.
if command -v pacman &> /dev/null; then
    if ! command -v gfortran &> /dev/null || ! command -v java &> /dev/null || ! command -v git &> /dev/null; then
        echo "--- Installing Arch Linux packages ---"
        sudo pacman -S --needed --noconfirm base-devel gcc-fortran git jre-openjdk
    fi
fi

# 1. Check for conda/mamba
echo "--- Checking for Conda/Mamba ---"
if command -v mamba &> /dev/null; then
    CONDA_EXE="mamba"
    echo "✓ Mamba found at: $(which mamba)"
elif command -v conda &> /dev/null; then
    CONDA_EXE="conda"
    echo "✓ Conda found at: $(which conda)"
else
    echo "❌ Neither conda nor mamba found. Please install Miniconda or Miniforge first:"
    echo "   https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi
# conda-forge + bioconda only: skips the "defaults" channel (and its Terms of Service prompt).
CHANNELS="--override-channels -c conda-forge -c bioconda --strict-channel-priority"

# 2. Setup Bakta Environment
echo ""
echo "--- Setting up Bakta Environment (bakta_env) ---"
if $CONDA_EXE env list | grep -q "^bakta_env "; then
    echo "✓ bakta_env already exists"
else
    echo "Creating bakta_env..."
    $CONDA_EXE create -n bakta_env $CHANNELS bakta -y
fi

# 3. Setup Panaroo Environment
echo ""
echo "--- Setting up Panaroo Environment (panaroo_env) ---"
if $CONDA_EXE env list | grep -q "^panaroo_env "; then
    echo "✓ panaroo_env already exists"
else
    echo "Creating panaroo_env..."
    $CONDA_EXE create -n panaroo_env $CHANNELS panaroo mafft -y
fi

# 4. Check Python version for tree pipeline
echo ""
echo "--- Setting up Python Virtual Environment ---"
python_version=$(python3 --version 2>&1 | awk '{print $2}')
echo "✓ Python version: $python_version"

# Activate virtual environment
if [ -d "venv" ]; then
    echo "✓ Virtual environment found"
else
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Install requirements
echo "Installing python dependencies..."
source venv/bin/activate
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
else
    echo "⚠ requirements.txt not found. Please ensure biopython and numpy are installed."
fi
deactivate

# 5. NCBI datasets CLI + VeryFastTree. They live in their own env (keeps conda base clean)
# and get linked into venv/bin, so `source venv/bin/activate` puts them on PATH.
echo ""
echo "--- Setting up CLI tools (ecotools env: datasets, veryfasttree) ---"
if $CONDA_EXE env list | grep -q "^ecotools "; then
    echo "✓ ecotools already exists"
else
    $CONDA_EXE create -n ecotools $CHANNELS ncbi-datasets-cli veryfasttree -y
fi
ECOTOOLS_BIN="$($CONDA_EXE env list | awk '$1=="ecotools" {print $NF}')/bin"
# bioconda names the binary VeryFastTree; link it as veryfasttree too.
for tool in datasets dataformat VeryFastTree; do
    if [ -x "$ECOTOOLS_BIN/$tool" ]; then
        ln -sf "$ECOTOOLS_BIN/$tool" "venv/bin/$tool"
        [ "$tool" = VeryFastTree ] && ln -sf "$ECOTOOLS_BIN/$tool" venv/bin/veryfasttree
        echo "✓ Linked $tool into venv/bin"
    else
        echo "⚠ $tool not found in $ECOTOOLS_BIN"
    fi
done
export PATH="$PWD/venv/bin:$PATH"  # so the checks below see the linked tools

# 5b. Check Java (required for EcoSim)
echo ""
echo "--- Checking Java (required for EcoSim, Step 8) ---"
if command -v java &> /dev/null; then
    echo "✓ Java found: $(java -version 2>&1 | head -1)"
else
    echo "⚠ Java not found. EcoSim (Step 8) will not run without it."
    echo "  Arch: sudo pacman -S jre-openjdk   (or https://adoptium.net)"
fi

# 5c. EcoSim's native tools. tools/bin ships macOS arm64 builds; Linux needs its own,
# built from the EcoSim 2.1.7 source that matches tools/ecosim.jar.
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo ""
    echo "--- Building EcoSim native tools for Linux (tools/linux/bin) ---"
    if [ -x tools/linux/bin/hillclimb ] && [ -x tools/linux/bin/fasttree ]; then
        echo "✓ Already built"
    else
        bash tools/build_ecosim_linux.sh tools/linux
    fi
fi

# 6. Check VeryFastTree / FastTree
echo ""
echo "--- Checking Tree Builders ---"
if command -v veryfasttree &> /dev/null; then
    echo "✓ VeryFastTree found at: $(which veryfasttree)"
elif command -v fasttree &> /dev/null; then
    echo "✓ FastTree found at: $(which fasttree)"
else
    echo "⚠ VeryFastTree/FastTree not found"
    
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo ""
        echo "Installing VeryFastTree via Homebrew..."
        if command -v brew &> /dev/null; then
            brew install veryfasttree
            echo "✓ VeryFastTree installed!"
        else
            echo "❌ Homebrew not found. Install Homebrew first: https://brew.sh"
            echo "   Then run: brew install veryfasttree"
        fi
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo ""
        echo "Installing VeryFastTree via apt..."
        if command -v apt-get &> /dev/null; then
            sudo apt-get update && sudo apt-get install -y veryfasttree
            echo "✓ VeryFastTree installed!"
        else
            echo "ℹ No veryfasttree package here; the pipeline falls back to tools/linux/bin/fasttree."
            echo "  For faster trees: conda install -c conda-forge -c bioconda veryfasttree"
        fi
    fi
fi

# Verify installation
if command -v veryfasttree &> /dev/null || command -v fasttree &> /dev/null || [ -x tools/linux/bin/fasttree ]; then
    echo "✓ Tree builder is ready!"
else
    echo "⚠ VeryFastTree/FastTree still not found. You may need to install manually."
fi

# Check for optional alignment tools
echo ""
echo "Checking optional alignment tools (for manual tree checks)..."
if command -v mafft &> /dev/null; then
    echo "✓ System MAFFT found at: $(which mafft)"
else
    echo "ℹ System MAFFT not installed (optional, Panaroo environment has its own)"
fi

echo ""
echo "================================================================"
echo "Setup Complete!"
echo "================================================================"
echo ""
echo "Next steps:"
echo "1. Activate your python virtual environment:"
echo "   source venv/bin/activate"
echo "2. Download the Bakta database if you haven't already:"
echo "   conda run -n bakta_env bakta_db download --output bakta_db --type light   (or --type full, ~70 GB)"
echo "3. EcoSim needs no config: it uses tools/ecosim.jar and the Linux helpers in tools/linux/bin."
echo "4. Run the full pipeline (including genome download):"
echo "   python pipeline.py --species 'Treponema pallidum' --sample-size 10 --outgroup 'Treponema paraluiscuniculi' --db bakta_db/db-light"
echo ""
