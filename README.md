# MSview — Spectrum Viewer

A desktop application for viewing and annotating mass spectra.

---

## What you need first

### Python 3.10 or newer
Download from **https://www.python.org/downloads/**

- **Windows:** Run the installer. On the first screen, **tick the box that says "Add Python to PATH"** before clicking Install.
- **Mac:** Download and run the .pkg installer. If you already have Homebrew, you can also run `brew install python` in Terminal.

---

## Running the app

### Windows
1. Double-click **`run_windows.bat`**
2. The first time only, it will install the required packages (takes ~1 minute)
3. The app will open automatically

### Mac
1. Open **Terminal** (press Cmd+Space, type "Terminal", press Enter)
2. Type the following and press Enter:
   ```
   cd ~/Downloads/msview
   ```
   (If you saved the folder somewhere else, replace `Downloads` with `Desktop` or wherever it is)
3. Then type the following and press Enter:
   ```
   chmod +x run_mac.sh && ./run_mac.sh
   ```
4. The first time only, it will install required packages (~1 minute)
5. The app will open automatically

**Next time:** just repeat steps 1–3. The install step is skipped after the first run.

---

## Loading a spectrum

- Click **Open file…** or drag a file directly onto the app window
- Supported formats: `.xy`, `.ascii`, `.txt` — any two-column text file with m/z in the first column and intensity in the second

---

## Features

| Feature | How to use |
|---|---|
| Zoom | Scroll wheel, or click and drag to box-zoom |
| Pan | Right-click and drag |
| Reset view | Click **Reset view** in the toolbar |
| Isotope overlay | Enter formula → click **Calculate & Overlay** |
| Adjust overlay intensity | Use the **Overlay intensity** slider |
| Annotate a peak | Click directly on a peak in the spectrum |
| Custom annotation | Enter m/z and label text → click **+** |
| Auto-label top peaks | Click **Top 5** or **Top 10** |
| Remove annotation | Double-click it in the list |
| Export image | Click **Export PNG** or **Export SVG** |

---

## Files in this folder

| File | Purpose |
|---|---|
| `msview.py` | Main application |
| `isotopes.py` | Isotope distribution calculator |
| `run_windows.bat` | Windows launcher |
| `run_mac.sh` | Mac launcher |
| `README.md` | This file |

---

## Troubleshooting

**"Python not found" on Windows**
Reinstall Python from python.org and make sure to tick "Add Python to PATH".

**App opens but looks odd on a high-resolution screen**
This is a known Qt scaling issue on some Windows machines. It will be addressed in a future version.

**Mac: "run_mac.sh cannot be opened because it is from an unidentified developer"**
Right-click the file → Open → Open anyway. This only happens once.

**Any other issue**
Note down what happened and report it for the next iteration.
