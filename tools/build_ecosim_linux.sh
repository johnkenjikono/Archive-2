#!/bin/bash
# Build EcoSim's native helpers for Linux into <outdir>/bin (default: tools/linux).
# tools/bin ships macOS arm64 builds; this uses the EcoSim 2.1.7 source that matches tools/ecosim.jar.
# Needs git, make, gcc, gfortran. Used by setup.sh and the Modal image (modal_app.py).
set -e
out="${1:-$(dirname "$0")/linux}/bin"
src=$(mktemp -d)
git clone -q https://github.com/sandain/ecosim.git "$src"
git -C "$src" checkout -q fc1e01973e4ee882fba8f01bd08bbe5064145c64
# FC must be explicit: make's built-in FC=f77 overrides the Makefile's "FC ?= gfortran".
make -C "$src" FC=gfortran CC=gcc build/c/fasttree \
    build/fortran/hillclimb build/fortran/npopCI build/fortran/omegaCI \
    build/fortran/sigmaCI build/fortran/demarcation
mkdir -p "$out"
cp "$src"/build/c/fasttree "$src"/build/fortran/{hillclimb,npopCI,omegaCI,sigmaCI,demarcation} "$out"/
rm -rf "$src"
echo "✓ Built EcoSim helpers into $out"
