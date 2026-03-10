@echo off
echo ============================================
echo  MSview - First Time Setup
echo ============================================
echo.

REM Check Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python not found.
    echo.
    echo Please install Python 3.10 or newer from https://www.python.org
    echo Make sure to tick "Add Python to PATH" during installation.
    echo.
    pause
    exit /b 1
)

echo Python found. Installing required packages...
echo This only happens once and may take a minute.
echo.

python -m pip install --upgrade pip
python -m pip install PyQt6 pyqtgraph numpy

echo.
echo Setup complete. Launching MSview...
echo.
python msview.py
pause
