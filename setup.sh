#!/bin/bash
# Setup script for full pipeline dependencies

set -e
cd "$(dirname "$0")"  # paths below (venv, tools/) are repo-relative

echo "================================================================"
echo "Ecotype Pipeline Setup Tool"
echo "================================================================"
echo ""

# 0. Linux system packages: compilers for EcoSim's tools, Java, and Python's venv/Tk support.
# Only prompts for sudo when something is missing.
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    missing=0
    for cmd in gfortran gcc make git java; do
        command -v "$cmd" &> /dev/null || missing=1
    done
    python3 -c "import ensurepip, venv, tkinter" &> /dev/null || missing=1
    if [ "$missing" = 1 ]; then
        if command -v pacman &> /dev/null; then
            echo "--- Installing Arch Linux packages ---"
            sudo pacman -S --needed --noconfirm base-devel gcc-fortran git jre-openjdk tk
        elif command -v apt-get &> /dev/null; then
            echo "--- Installing Debian/Ubuntu packages ---"
            sudo apt-get update
            sudo apt-get install -y build-essential gfortran git default-jre python3-venv python3-tk
        elif command -v dnf &> /dev/null; then
            echo "--- Installing Fedora packages ---"
            sudo dnf install -y gcc gcc-gfortran make git java-21-openjdk python3-tkinter
        else
            echo "⚠ Missing some of: gfortran gcc make git java python3-venv."
            echo "  Install them with your package manager, then rerun setup.sh."
        fi
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
# requirements.txt pins matplotlib 3.11, which needs Python 3.11+.
if ! python3 -c 'import sys; sys.exit(sys.version_info < (3, 11))'; then
    echo "❌ Python $python_version found; the pipeline needs Python 3.11 or newer."
    echo "   Install a newer python3, or run setup.sh from an activated conda env with Python 3.11+"
    echo "   (e.g. conda create -n py312 -c conda-forge python=3.12 && conda activate py312)."
    exit 1
fi
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
    elif command -v gfortran &> /dev/null && command -v gcc &> /dev/null && command -v make &> /dev/null && command -v git &> /dev/null; then
        bash tools/build_ecosim_linux.sh tools/linux
    else
        echo "⚠ Skipped: needs gfortran, gcc, make and git. Install them and rerun setup.sh,"
        echo "  or EcoSim (Step 8) will fail."
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
        echo "ℹ The pipeline falls back to tools/linux/bin/fasttree."
        echo "  For faster trees, recreate the ecotools env: conda env remove -n ecotools && bash setup.sh"
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
echo "   conda run -n bakta_env bakta_db download --output ~/bakta_db --type light   (or drop --type light for the full ~70 GB database)"
echo "3. EcoSim needs no config: it uses tools/ecosim.jar and the Linux helpers in tools/linux/bin."
echo "4. Run the full pipeline (including genome download):"
echo "   python pipeline.py --species 'Treponema pallidum' --sample-size 10 --outgroup 'Treponema paraluiscuniculi' --db ~/bakta_db/db-light"
echo ""
