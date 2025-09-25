#!/bin/bash
#
# Build script for SHG Python bindings
# 
# This script builds the Python extension module for the SHG library
#

set -e  # Exit on any error

echo "=== Building SHG Python Bindings ==="

# Check if we're in the right directory
if [ ! -f "setup.py" ]; then
    echo "Error: setup.py not found. Please run this script from the project root directory."
    exit 1
fi

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is required but not found."
    exit 1
fi

# Check if pip is available
if ! command -v pip3 &> /dev/null && ! command -v pip &> /dev/null; then
    echo "Error: pip is required but not found."
    exit 1
fi

# Set pip command
PIP_CMD="pip3"
if ! command -v pip3 &> /dev/null; then
    PIP_CMD="pip"
fi

echo "Installing Python dependencies..."
$PIP_CMD install -r requirements.txt

echo "Building Python extension..."
python3 setup.py build_ext --inplace

echo "Testing import..."
if python3 -c "import pyshg; print('Import successful!')"; then
    echo "Build completed successfully!"
    echo ""
    echo "You can now use the library:"
    echo "  python3 python/example_usage.py"
    echo ""
    echo "Or install it system-wide:"
    echo "  pip3 install ."
else
    echo "Build completed but import failed. Check for errors above."
    exit 1
fi