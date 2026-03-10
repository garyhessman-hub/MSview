#!/bin/bash
echo "============================================"
echo " MSview - First Time Setup"
echo "============================================"
echo

# Check Python 3 is available
if ! command -v python3 &>/dev/null; then
    echo "ERROR: Python 3 not found."
    echo
    echo "Please install Python 3.10 or newer from https://www.python.org"
    echo "Or via Homebrew: brew install python"
    echo
    read -p "Press Enter to exit."
    exit 1
fi

PYTHON=$(command -v python3)
echo "Using: $PYTHON ($($PYTHON --version))"
echo

# Create a virtual environment next to this script if it doesn't exist
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV="$SCRIPT_DIR/.venv"

if [ ! -d "$VENV" ]; then
    echo "Creating virtual environment..."
    $PYTHON -m venv "$VENV"
    echo
    echo "Installing required packages (one-time, may take a minute)..."
    "$VENV/bin/pip" install --upgrade pip --quiet
    "$VENV/bin/pip" install PyQt6 pyqtgraph numpy --quiet
    echo "Done."
else
    echo "Virtual environment found, skipping install."
fi

echo
echo "Launching MSview..."
echo
"$VENV/bin/python" "$SCRIPT_DIR/msview.py"
