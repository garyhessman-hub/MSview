# MSview — Spectrum Viewer

A desktop application for viewing and annotating mass spectra.

---

## Installation

### macOS
1. Download **MSview.dmg** from the latest build under the [Actions tab](../../actions)
2. Open the .dmg file
3. Drag **MSview** into your Applications folder
4. Double-click MSview to launch it

> First launch: if macOS says the app is from an unidentified developer, right-click the app → Open → Open anyway. This only happens once.

### Windows
1. Download **MSview_Setup.exe** from the latest build under the [Actions tab](../../actions)
2. Run the installer and follow the prompts
3. Launch MSview from the Start menu or desktop shortcut

---

## Loading a spectrum

- Click **Open file…** or drag a file directly onto the app window
- Supported formats: `.mzML`, `.mzXML`, `.xml` (Bruker DataAnalysis), `.ascii` (Bruker DataAnalysis), `.mgf`, `.xy`, `.txt` — and any two-column text file with m/z in the first column and intensity in the second

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
| Export data | Click **Export data…** → choose CSV, TXT, or .dat |

---

## Files in this repository

| File | Purpose |
|---|---|
| `msview.py` | Main application |
| `isotopes.py` | Isotope distribution calculator |
| `MSview_QuickStart.pdf` | Quick start guide for testers |
| `icon.png` / `msview.icns` / `msview.ico` | App icons |
| `.github/workflows/` | GitHub Actions build pipeline |

---

## Troubleshooting

**Mac: "MSview cannot be opened because it is from an unidentified developer"**
Right-click the app → Open → Open anyway. This only happens once.

**App opens but looks odd on a high-resolution screen**
This is a known Qt scaling issue on some Windows machines. It will be addressed in a future version.

**Any other issue**
Please report bugs or suggestions by email to [gary.hessman@tcd.ie](mailto:gary.hessman@tcd.ie) with subject line `[MSview Bug]` or `[MSview Suggestion]`.
