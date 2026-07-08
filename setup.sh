#!/bin/bash
# Setup script for full pipeline dependencies

set -e

echo "================================================================"
echo "Ecotype Pipeline Setup Tool"
echo "================================================================"
echo ""

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

# 2. Setup Bakta Environment
echo ""
echo "--- Setting up Bakta Environment (bakta_env) ---"
if $CONDA_EXE env list | grep -q "bakta_env"; then
    echo "✓ bakta_env already exists"
else
    echo "Creating bakta_env..."
    $CONDA_EXE create -n bakta_env -c conda-forge -c bioconda bakta -y
fi

# 3. Setup Panaroo Environment
echo ""
echo "--- Setting up Panaroo Environment (panaroo_env) ---"
if $CONDA_EXE env list | grep -q "panaroo_env"; then
    echo "✓ panaroo_env already exists"
else
    echo "Creating panaroo_env..."
    $CONDA_EXE create -n panaroo_env -c conda-forge -c bioconda panaroo mafft -y
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

# 5. Check NCBI Datasets CLI
echo ""
echo "--- Checking NCBI Datasets CLI ---"
if command -v datasets &> /dev/null; then
    echo "✓ NCBI datasets CLI found at: $(which datasets)"
else
    echo "⚠ NCBI datasets CLI not found (needed for downloading genomes in Step 1)."
    echo "  Installing via conda..."
    $CONDA_EXE install -c conda-forge ncbi-datasets-cli -y || \
        echo "  ❌ Auto-install failed. Install manually: conda install -c conda-forge ncbi-datasets-cli"
fi

# 5b. Check Java (required for EcoSim)
echo ""
echo "--- Checking Java (required for EcoSim, Step 8) ---"
if command -v java &> /dev/null; then
    echo "✓ Java found: $(java -version 2>&1 | head -1)"
else
    echo "⚠ Java not found. EcoSim (Step 8) will not run without it."
    echo "  Install a JDK: https://adoptium.net  or  conda install -c conda-forge openjdk"
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
            echo "❌ apt-get not found. Install with your package manager: veryfasttree"
        fi
    fi
fi

# Verify installation
if command -v veryfasttree &> /dev/null || command -v fasttree &> /dev/null; then
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
    echo "ℹ System MAFFT not installed (optional, Roary environment has its own)"
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
echo "   conda activate bakta_env"
echo "   bakta_db download --output /path/to/db"
echo "   conda deactivate"
echo "3. Configure EcoSim paths (only needed if running step 8):"
echo "   export ECOSIM_JAR=/path/to/ecosim.jar"
echo "   export ECOSIM_DIR=/path/to/ecosim"
echo "4. Run the full pipeline (including genome download):"
echo "   python pipeline.py --species 'Treponema paraluiscuniculi' --sample-size 5 --db /path/to/bakta/db"
echo ""
