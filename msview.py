"""
MSview — Mass Spectrum Viewer  v4
PyQt6 + pyqtgraph

New in v4:
  - Multi-spectrum overlay (load several files, each with its own colour/visibility)
  - Δm/z ruler tool (click two peaks, see mass difference + possible formula losses)
  - XML / mzML / mzXML / MGF file support
  - All previous features preserved
"""

import sys, os, re, math, tempfile
import xml.etree.ElementTree as ET
import struct, base64, zlib
import numpy as np

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QSplitter, QLabel, QLineEdit, QPushButton, QComboBox, QSlider,
    QFileDialog, QMessageBox, QScrollArea, QFrame,
    QDoubleSpinBox, QSpinBox, QCheckBox, QListWidget, QListWidgetItem,
    QStatusBar, QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
    QStackedWidget, QButtonGroup, QColorDialog,
)
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QFont, QColor, QAction, QGuiApplication, QImage

import pyqtgraph as pg
import pyqtgraph.exporters  # required: pg.exporters is not auto-loaded by the base import

from isotopes import parse_formula, format_formula, isotope_distribution, gaussian_profile

# ── Palettes ──────────────────────────────────────────────────────────────────
# First colour is the default for the primary spectrum
SPECTRUM_PALETTE = [
    "#2563eb",  # steel blue
    "#dc2626",  # crimson
    "#059669",  # emerald
    "#7c3aed",  # violet
    "#ea580c",  # orange
    "#0891b2",  # cyan
    "#be185d",  # pink
    "#65a30d",  # lime
]
SPECTRUM_COLORS = {
    "Steel blue": "#2563eb", "Black": "#111111", "Slate grey": "#475569",
    "Emerald": "#059669", "Crimson": "#dc2626", "Violet": "#7c3aed",
    "Orange": "#ea580c", "Cyan": "#0891b2",
}
OVERLAY_COLORS = {
    "Red": "#dc2626", "Orange": "#ea580c", "Green": "#059669",
    "Purple": "#7c3aed", "Blue": "#2563eb",
}
ANN_COLOR  = "#2563eb"
RULER_COLOR = "#f59e0b"   # amber

# Common neutral losses for Δm/z hints
NEUTRAL_LOSSES = {
    1.0079:  "H",    17.0027: "OH",   18.0106: "H₂O",  27.9949: "CO",
    28.0101: "CO / C₂H₄", 34.0055: "H₂S", 35.9898: "HCl",
    44.9977: "CO₂ / C₂H₅N", 46.0055: "C₂H₆O", 56.0262: "C₃H₄O",
    57.0215: "Gly", 64.0637: "C₄H₈O", 78.9590: "HBr", 79.9663: "SO₃",
    97.9769: "H₂SO₄–H₂O", 162.0528: "Hex", 176.0320: "GlcUA",
    203.0794: "HexNAc",
}


# ── File parsers ──────────────────────────────────────────────────────────────
def _parse_text(path):
    """Two-column ASCII: m/z  intensity"""
    mz, intensity = [], []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line[0].isalpha():
                continue
            parts = line.replace(",", " ").replace(";", " ").split()
            if len(parts) >= 2:
                try:
                    x, y = float(parts[0]), float(parts[1])
                    if x > 0:
                        mz.append(x); intensity.append(abs(y))
                except ValueError:
                    pass
    return np.array(mz, np.float64), np.array(intensity, np.float64)


def _decode_mzml_array(encoded, compression, dtype_str):
    """Decode a base64+zlib mzML binary array."""
    raw = base64.b64decode(encoded)
    if compression in ("zlib compression", "zlib"):
        raw = zlib.decompress(raw)
    dt = np.float64 if "64" in dtype_str else np.float32
    return np.frombuffer(raw, dtype=dt)


def _parse_mzml(path):
    """
    Parse first spectrum from mzML.
    Each <binaryDataArray> is parsed as a self-contained unit containing
    cvParams (array type, compression, precision) and a <binary> element.
    """
    tree = ET.parse(path)
    root = tree.getroot()

    def strip_ns(tag):
        return re.sub(r'[^}]+}', '', tag)

    def decode_array(binary_el, cv_names):
        text = binary_el.text
        if not text or not text.strip():
            return np.array([])
        raw = base64.b64decode(text.strip())
        if any("zlib" in n for n in cv_names):
            raw = zlib.decompress(raw)
        elif not any("no compression" in n for n in cv_names):
            try:
                raw = zlib.decompress(raw)
            except Exception:
                pass
        dt = np.float64 if any("64" in n for n in cv_names) else np.float32
        return np.frombuffer(raw, dtype=dt)

    for spectrum_el in root.iter():
        if strip_ns(spectrum_el.tag) != "spectrum":
            continue
        mz_arr = None
        int_arr = None
        for bda in spectrum_el.iter():
            if strip_ns(bda.tag) != "binaryDataArray":
                continue
            cv_names = []
            array_type = None
            binary_el = None
            for child in bda:
                ctag = strip_ns(child.tag)
                if ctag == "cvParam":
                    n = child.attrib.get("name", "")
                    cv_names.append(n)
                    if "m/z array" in n:
                        array_type = "mz"
                    elif "intensity array" in n:
                        array_type = "intensity"
                elif ctag == "binary":
                    binary_el = child
            if binary_el is None or array_type is None:
                continue
            arr = decode_array(binary_el, cv_names)
            if array_type == "mz":
                mz_arr = arr
            elif array_type == "intensity":
                int_arr = arr
        if mz_arr is not None and int_arr is not None and len(mz_arr) and len(int_arr):
            return mz_arr.astype(np.float64), int_arr.astype(np.float64)

    raise ValueError("No spectrum data found in mzML file.")


def _parse_mzxml(path):
    """
    Parse first scan from mzXML.
    Tolerates namespace variants and files where root contains peaks directly.
    """
    tree = ET.parse(path)
    root = tree.getroot()

    def strip_ns(tag):
        return re.sub(r'[^}]+}', '', tag)

    for el in root.iter():
        for child in el:
            if strip_ns(child.tag) == "peaks" and child.text and child.text.strip():
                try:
                    precision = int(child.attrib.get("precision", "32"))
                    compression = child.attrib.get("compressionType", "none").lower()
                    raw = base64.b64decode(child.text.strip())
                    if "zlib" in compression:
                        raw = zlib.decompress(raw)
                    dt = ">f8" if precision == 64 else ">f4"
                    data = np.frombuffer(raw, dtype=dt).astype(np.float64)
                    if len(data) >= 2:
                        mz = data[0::2]
                        intensity = data[1::2]
                        if len(mz):
                            return mz, intensity
                except Exception:
                    continue
    raise ValueError("No scan data found in mzXML file.")

def _parse_bruker_xml(path):
    """
    Parse Bruker DataAnalysis XML export.
    Spectrum data is stored as <pk mz="..." i="..."/> elements
    inside an <ms_peaks> element within <ms_spectrum>.
    """
    tree = ET.parse(path)
    root = tree.getroot()

    def strip_ns(tag):
        return re.sub(r'\{[^}]+\}', '', tag)

    for el in root.iter():
        if strip_ns(el.tag) == "ms_peaks":
            mz_list, int_list = [], []
            for child in el:
                if strip_ns(child.tag) == "pk":
                    try:
                        mz_list.append(float(child.attrib["mz"]))
                        int_list.append(float(child.attrib["i"]))
                    except (KeyError, ValueError):
                        pass
            if mz_list:
                return np.array(mz_list, np.float64), np.array(int_list, np.float64), "centroid"

    raise ValueError("No peak data found in Bruker XML file.")


def _parse_mgf(path):
    """Parse first spectrum from MGF (Mascot Generic Format)."""
    mz, intensity = [], []
    inside = False
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line == "BEGIN IONS":
                inside = True
                mz, intensity = [], []
            elif line == "END IONS":
                if mz:
                    return np.array(mz, np.float64), np.array(intensity, np.float64)
                inside = False
            elif inside and line and not "=" in line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz.append(float(parts[0]))
                        intensity.append(float(parts[1]))
                    except ValueError:
                        pass
    if mz:
        return np.array(mz, np.float64), np.array(intensity, np.float64)
    raise ValueError("No spectrum data found in MGF file.")


def load_spectrum(path):
    """
    Dispatch to the right parser based on extension.
    Returns (mz, intensity, mode) where mode is "centroid" or "profile".
    """
    ext = os.path.splitext(path)[1].lower()

    if ext in (".mzml", ".xml"):
        errors = []
        for parser, label in [(_parse_mzml,       "mzML"),
                               (_parse_mzxml,      "mzXML"),
                               (_parse_bruker_xml, "Bruker XML")]:
            try:
                result = parser(path)
                # parsers may return 2 or 3 values
                if len(result) == 3:
                    return result
                mz, intensity = result
                # mzML/mzXML — check metadata for centroid flag
                mode = _detect_mode_xml(path)
                return mz, intensity, mode
            except Exception as exc:
                errors.append(f"{label}: {exc}")
        raise ValueError("Could not parse XML file.\n" + "\n".join(errors))

    elif ext == ".mzxml":
        mz, intensity = _parse_mzxml(path)
        return mz, intensity, _detect_mode_xml(path)

    elif ext == ".mgf":
        mz, intensity = _parse_mgf(path)
        return mz, intensity, "centroid"   # MGF is always centroid

    else:
        mz, intensity = _parse_text(path)
        return mz, intensity, _detect_mode_text(mz)


def _detect_mode_xml(path):
    """
    Sniff the XML file for centroid/profile keywords in cvParams or attributes.
    Falls back to "profile" if ambiguous.
    """
    try:
        # Read first 8 KB — the spectrum metadata is always near the top
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            head = f.read(8192).lower()
        if "centroid" in head:
            return "centroid"
        if "profile" in head:
            return "profile"
        # Bruker XML uses scanmode attribute: fullscan = profile-like, but
        # <pk> elements mean it is already centroided
        if "<pk " in head or "<pk	" in head:
            return "centroid"
    except Exception:
        pass
    return "profile"


def _detect_mode_text(mz):
    """
    Heuristic for plain-text files: if peaks are evenly spaced it is profile,
    if spacings are irregular (centroid peak list) it is centroid.
    Uses coefficient of variation of m/z differences.
    """
    if len(mz) < 10:
        return "centroid"
    diffs = np.diff(mz)
    diffs = diffs[diffs > 0]
    if len(diffs) == 0:
        return "centroid"
    cv = diffs.std() / diffs.mean()
    # Profile data has very consistent spacing (cv < 0.1)
    return "profile" if cv < 0.1 else "centroid"


# ── Helpers ───────────────────────────────────────────────────────────────────
def make_label(text, bold=False, small=False):
    lbl = QLabel(text)
    f = lbl.font()
    if bold: f.setBold(True)
    if small: f.setPointSize(f.pointSize() - 1)
    lbl.setFont(f)
    return lbl


def make_section(title):
    """Plain widget with bold label — no QGroupBox margin/padding issues."""
    w = QWidget()
    lay = QVBoxLayout(w)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(6)
    if title:
        lbl = QLabel(title)
        lbl.setStyleSheet("font-weight:bold; font-size:12px; color:#0f172a; padding:0;")
        lay.addWidget(lbl)
    # Callers do QVBoxLayout(gb) which attaches to this widget — we store it
    w._inner_lay = lay
    return w


def delta_hint(delta):
    """Return a string hint if delta matches a known neutral loss."""
    best, best_name = None, None
    for mass, name in NEUTRAL_LOSSES.items():
        if abs(delta - mass) < 0.02:
            if best is None or abs(delta - mass) < abs(delta - best):
                best, best_name = mass, name
    return f"  ≈ {best_name}" if best_name else ""


# ── Custom ViewBox ─────────────────────────────────────────────────────────────
class SpectrumViewBox(pg.ViewBox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._rb_origin = None
        self._rb_rect   = None
        self._zoom_x    = True
        self._zoom_y    = False
        self.setMouseMode(pg.ViewBox.PanMode)
        self.setMenuEnabled(False)

    def wheelEvent(self, ev, axis=None):
        if not len(self.addedItems): return
        factor = 1.15 if ev.delta() < 0 else 1 / 1.15
        pos = self.mapToView(ev.pos())
        if self._zoom_x:
            cx = pos.x(); xr = self.viewRange()[0]
            lo = max(cx + (xr[0] - cx) * factor, 0.0)
            hi = max(cx + (xr[1] - cx) * factor, lo + 1e-9)
            self.setXRange(lo, hi, padding=0)
            self._clamp_y_to_visible()
        elif self._zoom_y:
            cy = pos.y(); yr = self.viewRange()[1]
            lo = max(cy + (yr[0] - cy) * factor, 0.0)
            hi = max(cy + (yr[1] - cy) * factor, lo + 1e-9)
            self.setYRange(lo, hi, padding=0)
        ev.accept()

    def mousePressEvent(self, ev):
        if ev.button() == Qt.MouseButton.RightButton:
            self._rb_origin = ev.pos()
            if self._rb_rect:
                self.scene().removeItem(self._rb_rect)
            from PyQt6.QtWidgets import QGraphicsRectItem
            from PyQt6.QtGui import QPen, QBrush
            self._rb_rect = QGraphicsRectItem(0, 0, 0, 0)
            self._rb_rect.setPen(QPen(QColor("#2563eb"), 1, Qt.PenStyle.DashLine))
            self._rb_rect.setBrush(QBrush(QColor(37, 99, 235, 30)))
            self._rb_rect.setZValue(1e9)
            self.scene().addItem(self._rb_rect)
            ev.accept()
        else:
            super().mousePressEvent(ev)

    def mouseMoveEvent(self, ev):
        if self._rb_origin is not None:
            from PyQt6.QtCore import QRectF
            sr = QRectF(self.mapToScene(self._rb_origin), self.mapToScene(ev.pos())).normalized()
            self._rb_rect.setRect(sr)
            ev.accept()
        else:
            super().mouseMoveEvent(ev)

    def mouseReleaseEvent(self, ev):
        if ev.button() == Qt.MouseButton.RightButton and self._rb_origin is not None:
            x0 = self.mapToView(self._rb_origin).x()
            x1 = self.mapToView(ev.pos()).x()
            if self._rb_rect:
                self.scene().removeItem(self._rb_rect); self._rb_rect = None
            self._rb_origin = None
            xr = self.viewRange()[0]
            if abs(x1 - x0) > (xr[1] - xr[0]) * 0.005:
                lo, hi = sorted([x0, x1])
                self.setXRange(max(lo, 0), hi, padding=0.02)
                self._clamp_y_to_visible()
            ev.accept()
        else:
            super().mouseReleaseEvent(ev)

    def mouseDragEvent(self, ev, axis=None):
        super().mouseDragEvent(ev, axis)
        if ev.isFinish(): self._clamp_to_data()

    def setRange(self, *args, **kwargs):
        super().setRange(*args, **kwargs)
        self._clamp_to_data()

    def _clamp_to_data(self):
        xr = self.viewRange()[0]
        yr = self.viewRange()[1]
        if xr[0] < 0:
            self.setXRange(0, xr[1] - xr[0], padding=0, update=False)
        if yr[0] < 0:
            self.setYRange(0, yr[1], padding=0, update=False)

    def _clamp_y_to_visible(self):
        xr = self.viewRange()[0]
        best = 0.0
        for item in self.addedItems:
            if hasattr(item, 'xData') and item.xData is not None:
                try:
                    mask = (item.xData >= xr[0]) & (item.xData <= xr[1])
                    vis  = item.yData[mask]
                    vis  = vis[~np.isnan(vis)]
                    if len(vis): best = max(best, float(vis.max()))
                except Exception: pass
        if best > 0:
            self.setYRange(0, best * 1.08, padding=0)


# ── Main Window ───────────────────────────────────────────────────────────────
class MSView(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MSview — Spectrum Viewer")
        self.setMinimumSize(1050, 640)
        self.resize(1320, 740)

        # App icon — use platform-appropriate format
        import sys as _sys
        _base = os.path.dirname(__file__)
        if _sys.platform == 'darwin':
            _icon_path = os.path.join(_base, 'msview.icns')
        elif _sys.platform == 'win32':
            _icon_path = os.path.join(_base, 'msview.ico')
        else:
            _icon_path = os.path.join(_base, 'icon.png')
        if os.path.exists(_icon_path):
            from PyQt6.QtGui import QIcon
            self.setWindowIcon(QIcon(_icon_path))

        # ── Multi-spectrum state ──────────────────────────────────────────────
        # Each entry: {name, mz, intensity, color, visible, plot_item}
        self.spectra = []

        # ── Single-spectrum helpers (always points at spectra[0] if present) ─
        @property
        def spec_mz(self):
            return self.spectra[0]["mz"] if self.spectra else np.array([])
        @property
        def spec_int(self):
            return self.spectra[0]["intensity"] if self.spectra else np.array([])

        # ── Other state ───────────────────────────────────────────────────────
        self.annotations   = []
        self.free_labels   = []   # list of {x, y, text, angle}
        self.iso_peaks     = None
        self.iso_params    = None
        self.overlay_scale = 1.0

        # ── Free-label placement mode ─────────────────────────────────────────
        self._fl_active     = False   # click-to-place mode on/off
        self._fl_pending_text = ""    # text waiting to be placed

        # ── Ruler state ──────────────────────────────────────────────────────
        self._ruler_active  = False
        self._ruler_pts     = []   # list of (mz, intensity) clicked
        self._ruler_items   = []   # plot items for ruler visuals

        # ── Plot items ───────────────────────────────────────────────────────
        self._overlay_item = None
        self._ann_items    = []
        self._fl_items     = []   # free-label TextItems

        self._rescale_timer = QTimer()
        self._rescale_timer.setSingleShot(True)
        self._rescale_timer.setInterval(150)
        self._rescale_timer.timeout.connect(self._rescale_overlay)

        self._build_ui()
        self._build_menu()
        self._apply_stylesheet()

    # ── Convenience accessors ─────────────────────────────────────────────────
    def _primary_mz(self):
        return self.spectra[0]["mz"] if self.spectra else np.array([])

    def _primary_int(self):
        return self.spectra[0]["intensity"] if self.spectra else np.array([])

    def _primary_disp(self):
        """Display intensity of primary spectrum (normalised if checkbox ticked)."""
        raw = self._primary_int()
        if self.chk_normalise.isChecked() and len(raw):
            mx = raw.max()
            return raw / mx * 100.0 if mx > 0 else raw
        return raw

    def _spec_disp(self, s):
        """Display intensity of any given spectrum (normalised if checkbox ticked).
        Generalisation of _primary_disp for use when iterating over multiple spectra,
        e.g. data export. Accepts a spectrum dict from self.spectra."""
        raw = s["intensity"]
        if self.chk_normalise.isChecked() and len(raw):
            mx = raw.max()
            return raw / mx * 100.0 if mx > 0 else raw
        return raw

    # ── UI construction ───────────────────────────────────────────────────────
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._build_sidebar())
        splitter.addWidget(self._build_plot_panel())
        splitter.setSizes([290, 1030])
        splitter.setHandleWidth(1)
        root.addWidget(splitter)

        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)
        self.statusbar.showMessage("No spectrum loaded")

        # Permanent credit label — right-aligned, always visible
        credit = QLabel("Built by Dr Gary Hessman with Claude")
        credit.setStyleSheet("color:#94a3b8; font-size:10px; padding-right:8px;")
        self.statusbar.addPermanentWidget(credit)

    def _build_sidebar(self):
        """
        Full-width sidebar. Top strip has two rows of word-tab buttons.
        Below that a QStackedWidget shows only the active panel.
        No side rail — panel content gets the full 290px width.
        """
        sidebar = QWidget()
        sidebar.setFixedWidth(290)
        sidebar.setObjectName("sidebar")

        outer = QVBoxLayout(sidebar)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # ── Tab strip (two rows of word buttons) ──────────────────────────────
        strip = QWidget()
        strip.setObjectName("tabStrip")
        strip_lay = QVBoxLayout(strip)
        strip_lay.setContentsMargins(6, 6, 6, 4)
        strip_lay.setSpacing(4)

        self._tab_group = QButtonGroup(self)
        self._tab_group.setExclusive(True)
        self._tab_btns = []

        # Row 1: File  Display  Ruler
        # Row 2: Pattern  Annotations
        rows = [["File", "Display", "Ruler"],
                ["Pattern", "Annotations"]]

        idx = 0
        for row_labels in rows:
            row_w = QWidget()
            row_lay = QHBoxLayout(row_w)
            row_lay.setContentsMargins(0, 0, 0, 0)
            row_lay.setSpacing(4)
            for label in row_labels:
                btn = QPushButton(label)
                btn.setObjectName("tabBtn")
                btn.setCheckable(True)
                btn.clicked.connect(lambda _, i=idx: self._switch_panel(i))
                self._tab_group.addButton(btn, idx)
                self._tab_btns.append(btn)
                row_lay.addWidget(btn)
                idx += 1
            strip_lay.addWidget(row_w)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setObjectName("tabLine")
        strip_lay.addWidget(sep)

        outer.addWidget(strip)

        # ── Stacked panels ────────────────────────────────────────────────────
        self._stack = QStackedWidget()
        self._stack.setObjectName("panelStack")

        builders = [self._build_file_section,
                    self._build_display_section,
                    self._build_ruler_section,
                    self._build_isotope_section,
                    self._build_annotation_section]

        for builder in builders:
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
            scroll.setObjectName("panelScroll")
            scroll.setFrameShape(QFrame.Shape.NoFrame)
            w = QWidget()
            w.setMinimumWidth(260)
            lay = QVBoxLayout(w)
            lay.setContentsMargins(10, 8, 10, 10)
            lay.setSpacing(8)
            lay.addWidget(builder())
            lay.addStretch()
            scroll.setWidget(w)
            self._stack.addWidget(scroll)

        outer.addWidget(self._stack)

        # Default: File tab
        self._tab_btns[0].setChecked(True)
        self._stack.setCurrentIndex(0)

        return sidebar

    def _switch_panel(self, idx):
        self._stack.setCurrentIndex(idx)

    def _build_file_section(self):
        gb = make_section("Spectra")
        lay = gb._inner_lay
        lay.setSpacing(5)

        row = QHBoxLayout()
        btn_add = QPushButton("Add spectrum…")
        btn_add.setObjectName("primaryBtn")
        btn_add.clicked.connect(self.open_file)
        row.addWidget(btn_add)
        btn_clr = QPushButton("Clear all")
        btn_clr.setObjectName("dangerBtn")
        btn_clr.clicked.connect(self.clear_all)
        row.addWidget(btn_clr)
        lay.addLayout(row)

        # Spectrum list
        self.spec_list = QListWidget()
        self.spec_list.setMaximumHeight(120)
        self.spec_list.setFont(QFont("Segoe UI", 10))
        self.spec_list.currentRowChanged.connect(self._on_spec_selected)
        self.spec_list.itemDoubleClicked.connect(self._toggle_spectrum_visibility)
        lay.addWidget(self.spec_list)
        lay.addWidget(make_label("Double-click to show/hide", small=True))

        # Per-spectrum controls row (colour + visibility + remove)
        colour_row = QHBoxLayout()
        colour_row.setSpacing(4)
        self.btn_colour = QPushButton("Change colour")
        self.btn_colour.setObjectName("smallBtn")
        self.btn_colour.clicked.connect(self._change_spectrum_colour)
        colour_row.addWidget(self.btn_colour)
        self._swatch = QLabel("  ")
        self._swatch.setFixedWidth(22)
        self._swatch.setFixedHeight(22)
        self._swatch.setObjectName("colourSwatch")
        self._swatch.setStyleSheet("background:#2563eb; border:1px solid #94a3b8; border-radius:3px;")
        colour_row.addWidget(self._swatch)
        colour_row.addStretch()
        lay.addLayout(colour_row)

        btn_rem = QPushButton("Remove selected")
        btn_rem.setObjectName("dangerBtn")
        btn_rem.clicked.connect(self._remove_selected_spectrum)
        lay.addWidget(btn_rem)

        row2 = QHBoxLayout()
        row2.addWidget(make_label("Ion mode:"))
        self.combo_mode = QComboBox()
        self.combo_mode.addItems(["Positive", "Negative"])
        row2.addWidget(self.combo_mode)
        lay.addLayout(row2)
        return gb

    def _build_display_section(self):
        gb = make_section("Display")
        lay = gb._inner_lay
        lay.setSpacing(5)

        lay.addWidget(make_label("Plot style:"))
        self.combo_style = QComboBox()
        self.combo_style.addItems(["Line (profile)", "Stick (centroid)"])
        self.combo_style.currentTextChanged.connect(self._replot_all)
        lay.addWidget(self.combo_style)

        self.chk_normalise = QCheckBox("Normalise each spectrum to 100%")
        self.chk_normalise.stateChanged.connect(self._replot_all)
        lay.addWidget(self.chk_normalise)
        return gb

    def _build_ruler_section(self):
        gb = make_section("Δm/z Ruler")
        lay = gb._inner_lay
        lay.setSpacing(5)

        lay.addWidget(make_label("Click two peaks to measure the mass difference:", small=True))

        self.btn_ruler = QPushButton("Activate ruler")
        self.btn_ruler.setCheckable(True)
        self.btn_ruler.setObjectName("rulerBtn")
        self.btn_ruler.clicked.connect(self._toggle_ruler)
        lay.addWidget(self.btn_ruler)

        self.lbl_ruler = QLabel("—")
        self.lbl_ruler.setObjectName("rulerLabel")
        self.lbl_ruler.setWordWrap(True)
        self.lbl_ruler.setFont(QFont("Courier", 11))
        lay.addWidget(self.lbl_ruler)

        btn_clr = QPushButton("Clear ruler")
        btn_clr.setObjectName("smallBtn")
        btn_clr.clicked.connect(self._clear_ruler)
        lay.addWidget(btn_clr)
        return gb

    def _build_isotope_section(self):
        gb = make_section("Isotope Overlay")
        lay = gb._inner_lay
        lay.setSpacing(5)

        lay.addWidget(make_label("Molecular formula:"))
        self.edit_formula = QLineEdit()
        self.edit_formula.setPlaceholderText("e.g. C26H25Cl2 (enter ion formula)")
        self.edit_formula.setFont(QFont("Courier", 11))
        lay.addWidget(self.edit_formula)

        r1 = QHBoxLayout()
        r1.addWidget(make_label("Charge (z):"))
        self.spin_charge = QSpinBox(); self.spin_charge.setRange(1, 20); self.spin_charge.setValue(1)
        r1.addWidget(self.spin_charge)
        lay.addLayout(r1)

        r2 = QHBoxLayout()
        r2.addWidget(make_label("FWHM (Da):"))
        self.spin_fwhm = QDoubleSpinBox()
        self.spin_fwhm.setRange(0.001, 10.0); self.spin_fwhm.setSingleStep(0.01)
        self.spin_fwhm.setDecimals(3); self.spin_fwhm.setValue(0.05)
        r2.addWidget(self.spin_fwhm)
        lay.addLayout(r2)

        lay.addWidget(make_label("Overlay colour:"))
        self.combo_iso_color = QComboBox()
        for name in OVERLAY_COLORS: self.combo_iso_color.addItem(name)
        lay.addWidget(self.combo_iso_color)

        self.lbl_formula_err = QLabel("")
        self.lbl_formula_err.setObjectName("errLabel")
        self.lbl_formula_err.setWordWrap(True)
        lay.addWidget(self.lbl_formula_err)

        btn_calc = QPushButton("Calculate && Overlay")
        btn_calc.setObjectName("primaryBtn")
        btn_calc.clicked.connect(self.calculate_overlay)
        lay.addWidget(btn_calc)

        self.scale_group = QWidget()
        sl = QVBoxLayout(self.scale_group)
        sl.setContentsMargins(0, 6, 0, 0); sl.setSpacing(3)
        sr = QHBoxLayout()
        sr.addWidget(make_label("Overlay intensity:"))
        self.lbl_scale = QLabel("100%"); self.lbl_scale.setObjectName("scaleLabel")
        sr.addWidget(self.lbl_scale, alignment=Qt.AlignmentFlag.AlignRight)
        sl.addLayout(sr)
        self.slider_scale = QSlider(Qt.Orientation.Horizontal)
        self.slider_scale.setRange(1, 300); self.slider_scale.setValue(100)
        self.slider_scale.valueChanged.connect(self._on_scale_slider)
        sl.addWidget(self.slider_scale)
        btn_rst = QPushButton("Reset to 100%"); btn_rst.setObjectName("smallBtn")
        btn_rst.clicked.connect(lambda: self.slider_scale.setValue(100))
        sl.addWidget(btn_rst)
        lay.addWidget(self.scale_group)
        self.scale_group.setVisible(False)

        self.iso_result_group = QWidget()
        irl = QVBoxLayout(self.iso_result_group); irl.setContentsMargins(0, 4, 0, 0)
        self.lbl_iso_peaks = QLabel("")
        self.lbl_iso_peaks.setFont(QFont("Courier", 10))
        self.lbl_iso_peaks.setObjectName("isoPeakLabel")
        self.lbl_iso_peaks.setWordWrap(True)
        irl.addWidget(self.lbl_iso_peaks)
        lay.addWidget(self.iso_result_group)
        self.iso_result_group.setVisible(False)

        self.btn_clear_overlay = QPushButton("Clear overlay")
        self.btn_clear_overlay.setObjectName("dangerBtn")
        self.btn_clear_overlay.clicked.connect(self.clear_overlay)
        self.btn_clear_overlay.setVisible(False)
        lay.addWidget(self.btn_clear_overlay)
        return gb

    def _build_annotation_section(self):
        gb = make_section("Annotations")
        lay = gb._inner_lay
        lay.setSpacing(5)

        lay.addWidget(make_label("Click a peak to label it, or add manually:", small=True))
        row = QHBoxLayout()
        self.edit_ann_mz = QLineEdit()
        self.edit_ann_mz.setPlaceholderText("m/z")
        self.edit_ann_mz.setMaximumWidth(88)
        self.edit_ann_mz.setFont(QFont("Courier", 11))
        row.addWidget(self.edit_ann_mz)
        self.edit_ann_label = QLineEdit()
        self.edit_ann_label.setPlaceholderText("label (optional)")
        row.addWidget(self.edit_ann_label)
        btn_add = QPushButton("+"); btn_add.setFixedWidth(28)
        btn_add.setObjectName("primaryBtn")
        btn_add.clicked.connect(self._add_manual_annotation)
        row.addWidget(btn_add)
        lay.addLayout(row)

        brow = QHBoxLayout()
        for lbl, n in [("Top 5", 5), ("Top 10", 10)]:
            b = QPushButton(lbl); b.setObjectName("smallBtn")
            b.clicked.connect(lambda _, n=n: self._auto_label(n))
            brow.addWidget(b)
        clr = QPushButton("Clear all"); clr.setObjectName("dangerBtn")
        clr.clicked.connect(self._clear_annotations)
        brow.addWidget(clr)
        lay.addLayout(brow)

        self.ann_table = QTableWidget(0, 4)
        self.ann_table.setHorizontalHeaderLabels(["m/z", "Intensity", "Rel %", "Label"])
        self.ann_table.setFont(QFont("Courier", 10))
        self.ann_table.setMaximumHeight(190)
        self.ann_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.ann_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.ann_table.setAlternatingRowColors(True)
        hdr = self.ann_table.horizontalHeader()
        for i, m in enumerate([QHeaderView.ResizeMode.ResizeToContents]*3 +
                               [QHeaderView.ResizeMode.Stretch]):
            hdr.setSectionResizeMode(i, m)
        self.ann_table.verticalHeader().setVisible(False)
        lay.addWidget(self.ann_table)

        btn_del = QPushButton("Remove selected row")
        btn_del.setObjectName("dangerBtn")
        btn_del.clicked.connect(self._remove_selected_annotation)
        lay.addWidget(btn_del)
        lay.addWidget(make_label("Click to highlight  ·  select + Remove to delete", small=True))

        # ── Free-floating label section ───────────────────────────────────────
        sep = QFrame(); sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("color:#d0d7de; margin-top:4px;")
        lay.addWidget(sep)

        lay.addWidget(make_label("Free text labels:", bold=True))
        lay.addWidget(make_label(
            "Type a label, activate placement, then click anywhere on the plot.",
            small=True))

        fl_row = QHBoxLayout(); fl_row.setSpacing(4)
        self.edit_fl_text = QLineEdit()
        self.edit_fl_text.setPlaceholderText("Label text…")
        fl_row.addWidget(self.edit_fl_text)

        self.combo_fl_angle = QComboBox()
        self.combo_fl_angle.addItems(["Horizontal", "45°", "Vertical"])
        self.combo_fl_angle.setFixedWidth(90)
        fl_row.addWidget(self.combo_fl_angle)
        lay.addLayout(fl_row)

        self.btn_fl_place = QPushButton("Place label on plot")
        self.btn_fl_place.setCheckable(True)
        self.btn_fl_place.setObjectName("rulerBtn")
        self.btn_fl_place.clicked.connect(self._toggle_fl_mode)
        lay.addWidget(self.btn_fl_place)

        # Free-label list
        self.fl_list = QListWidget()
        self.fl_list.setMaximumHeight(90)
        self.fl_list.setFont(QFont("Courier", 10))
        lay.addWidget(self.fl_list)

        btn_fl_del = QPushButton("Remove selected label")
        btn_fl_del.setObjectName("dangerBtn")
        btn_fl_del.clicked.connect(self._remove_selected_fl)
        lay.addWidget(btn_fl_del)

        btn_fl_clr = QPushButton("Clear all free labels")
        btn_fl_clr.setObjectName("smallBtn")
        btn_fl_clr.clicked.connect(self._clear_free_labels)
        lay.addWidget(btn_fl_clr)

        return gb

    def _build_plot_panel(self):
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)

        # Toolbar
        tb = QWidget(); tb.setObjectName("toolbar")
        tbl = QHBoxLayout(tb)
        tbl.setContentsMargins(8, 4, 8, 4); tbl.setSpacing(6)

        for label, slot in [("⟲  Reset view",       self.reset_view),
                             ("⎘  Copy to clipboard", self.copy_to_clipboard),
                             ("↓  Save PNG…",         lambda: self.export_image()),
                             ("⇟  Export data…",      lambda: self.export_data())]:
            b = QPushButton(label); b.setObjectName("toolBtn")
            b.clicked.connect(slot); tbl.addWidget(b)
            if label == "⟲  Reset view":
                sep = QFrame(); sep.setFrameShape(QFrame.Shape.VLine)
                sep.setObjectName("tbSep"); tbl.addWidget(sep)

        sep2 = QFrame(); sep2.setFrameShape(QFrame.Shape.VLine)
        sep2.setObjectName("tbSep"); tbl.addWidget(sep2)

        tbl.addWidget(make_label("Scroll zoom:", small=True))
        for axis in ["X", "Y"]:
            b = QPushButton(axis); b.setObjectName("zoomBtn")
            b.setCheckable(True); b.setChecked(axis == "X")
            b.setFixedWidth(30)
            b.clicked.connect(self._on_zoom_axis_changed)
            tbl.addWidget(b)
            if axis == "X":
                self.btn_zoom_x = b
            else:
                self.btn_zoom_y = b

        sep3 = QFrame(); sep3.setFrameShape(QFrame.Shape.VLine)
        sep3.setObjectName("tbSep"); tbl.addWidget(sep3)
        clr = QPushButton("✕  Clear all"); clr.setObjectName("toolBtnDanger")
        clr.clicked.connect(self.clear_all); tbl.addWidget(clr)
        tbl.addStretch()
        lay.addWidget(tb)

        # Plot
        pg.setConfigOptions(antialias=False, background="w", foreground="#334155")
        self._vb = SpectrumViewBox()
        self._vb._zoom_x = True; self._vb._zoom_y = False
        self.plot_widget = pg.PlotWidget(viewBox=self._vb)
        self.plot_widget.showGrid(x=True, y=True, alpha=0.25)
        self.plot_widget.setLabel("bottom", "m/z")
        self.plot_widget.setLabel("left", "Intensity")
        self.plot_widget.getAxis("bottom").setStyle(tickFont=QFont("Courier", 9))
        self.plot_widget.getAxis("left").setStyle(tickFont=QFont("Courier", 9))

        self._crosshair = pg.InfiniteLine(
            angle=90, movable=False,
            pen=pg.mkPen("#94a3b8", width=1, style=Qt.PenStyle.DashLine)
        )
        self.plot_widget.addItem(self._crosshair, ignoreBounds=True)

        self._mouse_proxy = pg.SignalProxy(
            self.plot_widget.scene().sigMouseMoved, rateLimit=30, slot=self._on_mouse_move
        )
        self.plot_widget.scene().sigMouseClicked.connect(self._on_plot_click)
        self.plot_widget.sigRangeChanged.connect(self._on_range_changed)
        lay.addWidget(self.plot_widget)
        return panel

    def _build_menu(self):
        mb = self.menuBar()
        fm = mb.addMenu("File")
        for label, shortcut, slot in [
            ("Add spectrum…",     "Ctrl+O",       self.open_file),
            (None, None, None),
            ("Copy to clipboard", "Ctrl+Shift+C", self.copy_to_clipboard),
            ("Save PNG…",         "Ctrl+Shift+P", self.export_image),
            ("Export data…",     "Ctrl+Shift+E", self.export_data),
            (None, None, None),
            ("Quit",              "Ctrl+Q",       self.close),
        ]:
            if label is None: fm.addSeparator()
            else:
                a = QAction(label, self)
                if shortcut: a.setShortcut(shortcut)
                a.triggered.connect(slot)
                fm.addAction(a)

    def _apply_stylesheet(self):
        self.setStyleSheet("""
            QMainWindow, QWidget { background:#fff; color:#0f172a;
                font-family:"Inter","Segoe UI","Helvetica Neue",sans-serif; font-size:12px; }
            #sidebar { background:#f4f6f8; border-right:1px solid #d0d7de; }
            #sidebarScroll { background:#f4f6f8; border:none; }
            #toolbar { background:#f4f6f8; border-bottom:1px solid #d0d7de; min-height:36px; }
            #tbSep { color:#d0d7de; max-height:18px; }
            QPushButton { background:#fff; border:1px solid #d0d7de; border-radius:4px;
                padding:5px 12px; color:#334155; }
            QPushButton:hover { background:#f0f7ff; border-color:#2563eb; color:#2563eb; }
            QPushButton:pressed { background:#dbeafe; }
            #primaryBtn { background:#2563eb; color:#fff; border:none; font-weight:bold; }
            #primaryBtn:hover { background:#1d4ed8; color:#fff; }
            #dangerBtn { color:#dc2626; border-color:#fca5a5; }
            #dangerBtn:hover { background:#fef2f2; border-color:#dc2626; }
            #smallBtn { padding:3px 8px; font-size:11px; }
            #toolBtn { background:transparent; border:1px solid #d0d7de;
                padding:3px 10px; font-size:11px; }
            #toolBtn:hover { background:#f0f7ff; border-color:#2563eb; color:#2563eb; }
            #toolBtnDanger { background:transparent; border:1px solid #d0d7de;
                padding:3px 10px; font-size:11px; color:#dc2626; }
            #toolBtnDanger:hover { background:#fef2f2; border-color:#dc2626; }
            #zoomBtn { padding:2px 4px; font-size:11px; font-weight:bold;
                border:1px solid #d0d7de; border-radius:3px; background:#fff; color:#475569; }
            #zoomBtn:checked { background:#2563eb; color:#fff; border-color:#2563eb; }
            #zoomBtn:hover:!checked { background:#f0f7ff; border-color:#2563eb; color:#2563eb; }
            #rulerBtn { background:#fff; border:1px solid #f59e0b; color:#92400e;
                border-radius:4px; padding:4px 10px; }
            #rulerBtn:checked { background:#fef3c7; border-color:#d97706; color:#78350f;
                font-weight:bold; }
            #rulerBtn:hover { background:#fef9c3; }
            #rulerLabel { color:#92400e; background:#fef3c7; border:1px solid #fde68a;
                border-radius:4px; padding:5px 8px; min-height:30px; }
            QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox { background:#fff;
                border:1px solid #d0d7de; border-radius:4px; padding:4px 7px; color:#0f172a; }
            QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus
                { border-color:#2563eb; }
            QComboBox::drop-down { border:none; }
            QComboBox::down-arrow { image:none; width:12px; }
            QSlider::groove:horizontal { height:4px; background:#e2e8f0; border-radius:2px; }
            QSlider::handle:horizontal { background:#2563eb; width:14px; height:14px;
                margin:-5px 0; border-radius:7px; }
            QSlider::sub-page:horizontal { background:#93c5fd; border-radius:2px; }
            #errLabel { color:#dc2626; font-size:11px; }
            #fileLabel { color:#475569; font-size:11px; }
            #scaleLabel { color:#2563eb; font-weight:bold; font-size:12px; }
            #isoPeakLabel { color:#334155; line-height:1.6; }
            QListWidget { background:#fff; border:1px solid #e2e8f0; border-radius:4px; }
            QListWidget::item { padding:4px 6px; }
            QListWidget::item:hover { background:#f0f7ff; }
            QListWidget::item:selected { background:#dbeafe; color:#1e3a5f; }
            QScrollBar:vertical { width:6px; background:transparent; }
            QScrollBar::handle:vertical { background:#d0d7de; border-radius:3px; min-height:20px; }
            QStatusBar { color:#64748b; font-size:11px; }
            QTableWidget { border:1px solid #e2e8f0; border-radius:4px;
                gridline-color:#f1f5f9; font-size:11px; }
            QTableWidget::item { padding:2px 4px; }
            QTableWidget::item:selected { background:#dbeafe; color:#1e3a5f; }
            QTableWidget::item:alternate { background:#f8fafc; }
            QHeaderView::section { background:#f4f6f8; border:none;
                border-bottom:1px solid #d0d7de; padding:3px 4px;
                font-size:10px; font-weight:bold; color:#475569; }
            #tabStrip { background:#f0f2f5; border-bottom:1px solid #d0d7de; }
            #tabLine { color:#d0d7de; margin:0; }
            #tabBtn { background:#e2e8f0; border:1px solid #cbd5e1;
                border-radius:5px; font-size:11px; font-weight:600;
                color:#475569; padding:5px 6px; min-width:0; }
            #tabBtn:hover { background:#dde4ec; color:#1e293b; border-color:#94a3b8; }
            #tabBtn:checked { background:#2563eb; color:#ffffff; border-color:#1d4ed8; }
            #panelStack { background:#f8fafc; }
            #panelScroll { background:#f8fafc; border:none; }
        """)

    # ── File loading ──────────────────────────────────────────────────────────
    def open_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Spectrum File", "",
            "Spectrum files (*.xy *.ascii *.txt *.csv *.mzml *.mzxml *.xml *.mgf)"
            ";;All files (*)"
        )
        if path: self._load_file(path)

    def _load_file(self, path):
        try:
            mz, intensity, mode = load_spectrum(path)
            if not len(mz):
                QMessageBox.warning(self, "Empty file", "No spectrum data found in this file.")
                return
            # Auto-switch plot style on first spectrum only
            if len(self.spectra) == 0:
                target = "Stick (centroid)" if mode == "centroid" else "Line (profile)"
                if self.combo_style.currentText() != target:
                    self.combo_style.blockSignals(True)
                    self.combo_style.setCurrentText(target)
                    self.combo_style.blockSignals(False)
            # Assign next colour from palette
            color = SPECTRUM_PALETTE[len(self.spectra) % len(SPECTRUM_PALETTE)]
            spec = {
                "name":      os.path.basename(path),
                "mz":        mz,
                "intensity": intensity,
                "color":     color,
                "visible":   True,
                "plot_item": None,
            }
            self.spectra.append(spec)
            self._add_spectrum_to_plot(spec)
            self._update_spec_list()
            self._update_status()
            # Reset annotations when first spectrum loaded
            if len(self.spectra) == 1:
                self.annotations.clear()
                self._refresh_annotation_list()
        except Exception as e:
            QMessageBox.critical(self, "Error loading file", str(e))

    def dragEnterEvent(self, ev):
        if ev.mimeData().hasUrls(): ev.acceptProposedAction()

    def dropEvent(self, ev):
        for url in ev.mimeData().urls():
            path = url.toLocalFile()
            if path: self._load_file(path)

    # ── Spectrum list management ──────────────────────────────────────────────
    def _update_spec_list(self):
        self.spec_list.clear()
        for i, spec in enumerate(self.spectra):
            txt = ("● " if spec["visible"] else "○ ") + spec["name"]
            item = QListWidgetItem(txt)
            item.setForeground(QColor(spec["color"]))
            self.spec_list.addItem(item)
        # Refresh swatch to currently selected row
        idx = self.spec_list.currentRow()
        if 0 <= idx < len(self.spectra):
            c = self.spectra[idx]["color"]
            self._swatch.setStyleSheet(
                f"background:{c}; border:1px solid #94a3b8; border-radius:3px;")

    def _on_spec_selected(self, idx):
        """Update the colour swatch to show the selected spectrum's colour."""
        if 0 <= idx < len(self.spectra):
            c = self.spectra[idx]["color"]
            self._swatch.setStyleSheet(
                f"background:{c}; border:1px solid #94a3b8; border-radius:3px;")

    def _change_spectrum_colour(self):
        """Open a colour dialog and update the selected spectrum."""
        idx = self.spec_list.currentRow()
        if not (0 <= idx < len(self.spectra)):
            return
        current = QColor(self.spectra[idx]["color"])
        colour = QColorDialog.getColor(current, self, "Choose spectrum colour")
        if not colour.isValid():
            return
        hex_col = colour.name()
        self.spectra[idx]["color"] = hex_col
        self._swatch.setStyleSheet(
            f"background:{hex_col}; border:1px solid #94a3b8; border-radius:3px;")
        # Rebuild only this spectrum's plot item
        self._add_spectrum_to_plot(self.spectra[idx])
        self._update_spec_list()

    def _toggle_spectrum_visibility(self, item):
        idx = self.spec_list.row(item)
        if 0 <= idx < len(self.spectra):
            self.spectra[idx]["visible"] = not self.spectra[idx]["visible"]
            pi = self.spectra[idx]["plot_item"]
            if pi:
                pi.setVisible(self.spectra[idx]["visible"])
            self._update_spec_list()

    def _remove_selected_spectrum(self):
        rows = sorted(set(
            self.spec_list.currentRow()
            for _ in self.spec_list.selectedItems()
        ), reverse=True)
        if not rows and self.spec_list.currentRow() >= 0:
            rows = [self.spec_list.currentRow()]
        for idx in rows:
            if 0 <= idx < len(self.spectra):
                pi = self.spectra[idx]["plot_item"]
                if pi: self.plot_widget.removeItem(pi)
                self.spectra.pop(idx)
        self._update_spec_list()
        self._update_status()

    # ── Plotting ──────────────────────────────────────────────────────────────
    def _add_spectrum_to_plot(self, spec):
        """Create and add a PlotDataItem for this spectrum entry."""
        if spec["plot_item"] is not None:
            self.plot_widget.removeItem(spec["plot_item"])

        mz, intensity = spec["mz"], spec["intensity"]
        if self.chk_normalise.isChecked():
            mx = intensity.max()
            intensity = intensity / mx * 100.0 if mx > 0 else intensity

        color = spec["color"]
        if self.combo_style.currentText().startswith("Stick"):
            n = len(mz)
            x = np.empty(n * 3); y = np.empty(n * 3)
            x[0::3] = mz; x[1::3] = mz; x[2::3] = mz
            y[0::3] = 0;  y[1::3] = intensity; y[2::3] = np.nan
            pi = pg.PlotDataItem(x=x, y=y,
                pen=pg.mkPen(color=color, width=1.2),
                connect="finite", skipFiniteCheck=False)
        else:
            pi = pg.PlotDataItem(x=mz, y=intensity,
                pen=pg.mkPen(color=color, width=1.5),
                fillLevel=0,
                brush=pg.mkBrush(QColor(color).lighter(185)))

        pi.setVisible(spec["visible"])
        self.plot_widget.addItem(pi)
        spec["plot_item"] = pi

        ylabel = "Relative Intensity (%)" if self.chk_normalise.isChecked() else "Intensity"
        self.plot_widget.setLabel("left", ylabel)

    def _replot_all(self):
        """Rebuild all spectrum items (style/normalise changed)."""
        for spec in self.spectra:
            self._add_spectrum_to_plot(spec)
        self._refresh_overlay()
        self._redraw_annotations()
        self._redraw_free_labels()

    # ── Ruler ─────────────────────────────────────────────────────────────────
    def _toggle_ruler(self):
        self._ruler_active = self.btn_ruler.isChecked()
        if self._ruler_active:
            self._ruler_pts.clear()
            self._clear_ruler_items()
            self.lbl_ruler.setText("Click first peak…")
            self.plot_widget.setCursor(Qt.CursorShape.CrossCursor)
        else:
            self._clear_ruler()
            self.plot_widget.setCursor(Qt.CursorShape.ArrowCursor)

    def _clear_ruler(self):
        self._ruler_active = False
        self.btn_ruler.setChecked(False)
        self._ruler_pts.clear()
        self._clear_ruler_items()
        self.lbl_ruler.setText("—")
        self.plot_widget.setCursor(Qt.CursorShape.ArrowCursor)

    def _clear_ruler_items(self):
        for item in self._ruler_items:
            self.plot_widget.removeItem(item)
        self._ruler_items.clear()

    def _ruler_click(self, clicked_mz):
        """Called from _on_plot_click when ruler is active."""
        if not self.spectra: return

        # Snap to local peak maximum (same logic as annotation)
        mz = self._primary_mz()
        disp = self._primary_disp()
        xr = self.plot_widget.viewRange()[0]
        search_r = (xr[1] - xr[0]) * 0.03
        mask = (mz >= clicked_mz - search_r) & (mz <= clicked_mz + search_r)
        if not mask.any(): return
        local_mz = mz[mask]; local_int = disp[mask]
        snapped_mz = float(local_mz[int(np.argmax(local_int))])
        snapped_int = float(local_int[int(np.argmax(local_int))])

        # Draw a vertical marker
        vline = pg.InfiniteLine(pos=snapped_mz, angle=90, movable=False,
            pen=pg.mkPen(RULER_COLOR, width=2, style=Qt.PenStyle.DashLine))
        self.plot_widget.addItem(vline)
        self._ruler_items.append(vline)
        self._ruler_pts.append((snapped_mz, snapped_int))

        if len(self._ruler_pts) == 1:
            self.lbl_ruler.setText(f"A: {snapped_mz:.4f}\nClick second peak…")
        elif len(self._ruler_pts) >= 2:
            mz_a, int_a = self._ruler_pts[-2]
            mz_b, int_b = self._ruler_pts[-1]
            delta = abs(mz_b - mz_a)
            hint  = delta_hint(delta)

            # Draw horizontal connector line
            y_mid = max(int_a, int_b) * 0.5
            hline = pg.PlotDataItem(
                x=[mz_a, mz_b], y=[y_mid, y_mid],
                pen=pg.mkPen(RULER_COLOR, width=1.5))
            self.plot_widget.addItem(hline)
            self._ruler_items.append(hline)

            # Label
            txt = pg.TextItem(f"Δ {delta:.4f}{hint}", color=RULER_COLOR, anchor=(0.5, 1))
            txt.setFont(QFont("Courier", 10))
            txt.setPos((mz_a + mz_b) / 2, y_mid)
            self.plot_widget.addItem(txt)
            self._ruler_items.append(txt)

            self.lbl_ruler.setText(
                f"A: {mz_a:.4f}\nB: {mz_b:.4f}\nΔ: {delta:.4f} Da{hint}"
            )
            # Keep ruler active for more measurements; clear pts for next pair
            self._ruler_pts.clear()

    # ── Isotope overlay ───────────────────────────────────────────────────────
    def calculate_overlay(self):
        formula_str = self.edit_formula.text().strip()
        self.lbl_formula_err.setText("")
        if not formula_str:
            self.lbl_formula_err.setText("Enter a molecular formula."); return
        try: atoms = parse_formula(formula_str)
        except ValueError as e: self.lbl_formula_err.setText(str(e)); return

        z     = self.spin_charge.value()
        fwhm  = self.spin_fwhm.value()
        color = OVERLAY_COLORS.get(self.combo_iso_color.currentText(), "#dc2626")

        self.iso_peaks  = isotope_distribution(atoms, z)
        self.iso_params = {"fwhm": fwhm, "color": color,
                           "label": f"Calc. {format_formula(atoms)}  z={z}",
                           "formula": format_formula(atoms), "z": z}

        self.slider_scale.blockSignals(True)
        self.slider_scale.setValue(100)
        self.slider_scale.blockSignals(False)
        self.lbl_scale.setText("100%"); self.overlay_scale = 1.0

        self.scale_group.setVisible(True)
        self.btn_clear_overlay.setVisible(True)

        txt  = f"<b>{format_formula(atoms)}</b>  z={z}<br>"
        txt += f"Monoisotopic m/z: <b>{self.iso_peaks[0][0]:.4f}</b><br><br>"
        txt += "<b>m/z&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Rel.%</b><br>"
        for mz, rel in self.iso_peaks:
            txt += f"{mz:.4f}&nbsp;&nbsp;{rel*100:.1f}%<br>"
        self.lbl_iso_peaks.setText(txt)
        self.iso_result_group.setVisible(True)
        self._refresh_overlay()

    def _visible_max(self):
        if not self.spectra: return 1.0
        xr = self.plot_widget.viewRange()[0]
        best = 0.0
        for spec in self.spectra:
            if not spec["visible"]: continue
            mz, intensity = spec["mz"], spec["intensity"]
            if self.chk_normalise.isChecked():
                mx = intensity.max()
                intensity = intensity / mx * 100.0 if mx > 0 else intensity
            mask = (mz >= xr[0]) & (mz <= xr[1])
            vis  = intensity[mask]
            if len(vis): best = max(best, float(vis.max()))
        return best if best > 0 else 1.0

    def _refresh_overlay(self):
        if self._overlay_item is not None:
            self.plot_widget.removeItem(self._overlay_item)
            self._overlay_item = None
        if not self.iso_peaks or not self.iso_params: return
        mz_arr, int_arr = gaussian_profile(
            self.iso_peaks, self.iso_params["fwhm"],
            self._visible_max() * self.overlay_scale)
        self._overlay_item = pg.PlotDataItem(
            x=np.array(mz_arr), y=np.array(int_arr),
            pen=pg.mkPen(color=self.iso_params["color"], width=2),
            name=self.iso_params["label"])
        self.plot_widget.addItem(self._overlay_item)

    def _rescale_overlay(self):
        if not self.iso_peaks or not self.iso_params: return
        mz_arr, int_arr = gaussian_profile(
            self.iso_peaks, self.iso_params["fwhm"],
            self._visible_max() * self.overlay_scale)
        if self._overlay_item is not None:
            self._overlay_item.setData(x=np.array(mz_arr), y=np.array(int_arr))
        else:
            self._refresh_overlay()

    def _on_scale_slider(self, value):
        self.overlay_scale = value / 100.0
        self.lbl_scale.setText(f"{value}%")
        if self.iso_peaks: self._rescale_overlay()

    def clear_overlay(self):
        self.iso_peaks = None; self.iso_params = None
        if self._overlay_item is not None:
            self.plot_widget.removeItem(self._overlay_item)
            self._overlay_item = None
        self.scale_group.setVisible(False)
        self.iso_result_group.setVisible(False)
        self.btn_clear_overlay.setVisible(False)

    def _on_range_changed(self):
        if self.iso_peaks: self._rescale_timer.start()

    # ── Annotations ───────────────────────────────────────────────────────────
    def _on_plot_click(self, event):
        if event.button() != Qt.MouseButton.LeftButton: return
        pos = event.scenePos()
        vb  = self.plot_widget.getViewBox()
        pt  = vb.mapSceneToView(pos)
        clicked_mz = pt.x()
        clicked_y  = pt.y()

        # Free-label placement takes priority over everything
        if self._fl_active:
            self._place_free_label(clicked_mz, clicked_y)
            return

        if not self.spectra: return

        if self._ruler_active:
            self._ruler_click(clicked_mz)
            return

        mz   = self._primary_mz()
        disp = self._primary_disp()
        xr   = self.plot_widget.viewRange()[0]
        search_r = (xr[1] - xr[0]) * 0.03
        mask = (mz >= clicked_mz - search_r) & (mz <= clicked_mz + search_r)
        if not mask.any(): return

        local_mz  = mz[mask]; local_int = disp[mask]
        snapped   = float(local_mz[int(np.argmax(local_int))])
        if any(abs(a["mz"] - snapped) < 0.0001 for a in self.annotations): return
        self._add_annotation(snapped, f"{snapped:.4f}")

    def _add_manual_annotation(self):
        try: mz = float(self.edit_ann_mz.text())
        except ValueError:
            QMessageBox.warning(self, "Invalid m/z", "Enter a numeric m/z value."); return
        label = self.edit_ann_label.text().strip() or f"{mz:.4f}"
        self._add_annotation(mz, label)
        self.edit_ann_mz.clear(); self.edit_ann_label.clear()

    def _auto_label(self, n):
        if not self.spectra: return
        disp = self._primary_disp()
        mz   = self._primary_mz()
        top  = np.argsort(disp)[-n:][::-1]
        self.annotations.clear()
        for idx in top:
            self.annotations.append({"mz": float(mz[idx]), "label": f"{mz[idx]:.4f}"})
        self._refresh_annotation_list()
        self._redraw_annotations()

    def _add_annotation(self, mz, label):
        self.annotations.append({"mz": mz, "label": label})
        self._refresh_annotation_list()
        self._redraw_annotations()

    def _remove_selected_annotation(self):
        rows = sorted(set(i.row() for i in self.ann_table.selectedItems()), reverse=True)
        sorted_anns = sorted(self.annotations, key=lambda x: x["mz"])
        for row in rows:
            if 0 <= row < len(sorted_anns):
                try: self.annotations.remove(sorted_anns[row])
                except ValueError: pass
        self._refresh_annotation_list(); self._redraw_annotations()

    def _clear_annotations(self):
        self.annotations.clear()
        self._refresh_annotation_list(); self._redraw_annotations()

    def _refresh_annotation_list(self):
        self.ann_table.setRowCount(0)
        if not self.annotations or not self.spectra: return
        disp = self._primary_disp()
        mz   = self._primary_mz()
        rows = []
        for a in sorted(self.annotations, key=lambda x: x["mz"]):
            idx = int(np.argmin(np.abs(mz - a["mz"])))
            rows.append((a["mz"], float(disp[idx]), a["label"]))
        cluster_max = max(r[1] for r in rows) or 1.0
        for pmz, intensity, label in rows:
            rel = intensity / cluster_max * 100.0
            r   = self.ann_table.rowCount()
            self.ann_table.insertRow(r)
            for col, val in enumerate([f"{pmz:.4f}", f"{intensity:.0f}",
                                        f"{rel:.1f}%", label]):
                item = QTableWidgetItem(val)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.ann_table.setItem(r, col, item)

    # ── Free-floating labels ─────────────────────────────────────────────────

    def _toggle_fl_mode(self):
        self._fl_active = self.btn_fl_place.isChecked()
        if self._fl_active:
            txt = self.edit_fl_text.text().strip()
            if not txt:
                QMessageBox.warning(self, "No label text",
                    "Enter label text before activating placement.")
                self.btn_fl_place.setChecked(False)
                self._fl_active = False
                return
            self._fl_pending_text = txt
            self.plot_widget.setCursor(Qt.CursorShape.CrossCursor)
            self.statusbar.showMessage(
                f'Click anywhere on the plot to place "{txt}"', 0)
        else:
            self.plot_widget.setCursor(Qt.CursorShape.ArrowCursor)
            self.statusbar.showMessage("")

    def _place_free_label(self, x, y):
        angle_map = {"Horizontal": 0, "45°": 45, "Vertical": 90}
        angle = angle_map.get(self.combo_fl_angle.currentText(), 0)
        self.free_labels.append({
            "x": x, "y": y,
            "text": self._fl_pending_text,
            "angle": angle,
        })
        self._fl_active = False
        self.btn_fl_place.setChecked(False)
        self.plot_widget.setCursor(Qt.CursorShape.ArrowCursor)
        self.statusbar.showMessage(
            f'Label "{self._fl_pending_text}" placed at m/z {x:.4f}', 3000)
        self.edit_fl_text.clear()
        self._refresh_fl_list()
        self._redraw_free_labels()

    def _refresh_fl_list(self):
        self.fl_list.clear()
        for fl in self.free_labels:
            angle_str = f"{fl['angle']}°" if fl['angle'] else "H"
            self.fl_list.addItem(f"{fl['x']:.4f}  [{angle_str}]  {fl['text']}")

    def _remove_selected_fl(self):
        idx = self.fl_list.currentRow()
        if 0 <= idx < len(self.free_labels):
            self.free_labels.pop(idx)
            self._refresh_fl_list()
            self._redraw_free_labels()

    def _clear_free_labels(self):
        self.free_labels.clear()
        self._refresh_fl_list()
        self._redraw_free_labels()

    def _redraw_free_labels(self):
        for item in self._fl_items:
            self.plot_widget.removeItem(item)
        self._fl_items.clear()
        for fl in self.free_labels:
            ti = pg.TextItem(
                text=fl["text"],
                color="#0f172a",
                anchor=(0, 1),
                angle=fl["angle"],
            )
            ti.setFont(QFont("Segoe UI", 10))
            ti.setPos(fl["x"], fl["y"])
            self.plot_widget.addItem(ti)
            self._fl_items.append(ti)

    def _redraw_annotations(self):
        for line, text in self._ann_items:
            self.plot_widget.removeItem(line)
            self.plot_widget.removeItem(text)
        self._ann_items.clear()
        if not self.spectra: return
        disp = self._primary_disp()
        mz   = self._primary_mz()
        for ann in self.annotations:
            idx = int(np.argmin(np.abs(mz - ann["mz"])))
            y   = float(disp[idx])
            line = pg.InfiniteLine(pos=ann["mz"], angle=90, movable=False,
                pen=pg.mkPen(ANN_COLOR, width=1, style=Qt.PenStyle.DashLine))
            text = pg.TextItem(text=ann["label"], color=ANN_COLOR, anchor=(0, 1), angle=90)
            text.setFont(QFont("Courier", 9))
            text.setPos(ann["mz"], y * 1.02)
            self.plot_widget.addItem(line)
            self.plot_widget.addItem(text)
            self._ann_items.append((line, text))

    # ── Mouse ─────────────────────────────────────────────────────────────────
    def _on_mouse_move(self, evt):
        pos = evt[0]
        if not self.plot_widget.sceneBoundingRect().contains(pos):
            self._crosshair.setVisible(False); return
        mz = self.plot_widget.getViewBox().mapSceneToView(pos).x()
        self._crosshair.setPos(mz); self._crosshair.setVisible(True)
        if self.spectra:
            pmz  = self._primary_mz()
            disp = self._primary_disp()
            idx  = int(np.argmin(np.abs(pmz - mz)))
            suffix = "%" if self.chk_normalise.isChecked() else ""
            ruler_hint = "  [RULER ACTIVE — click a peak]" if self._ruler_active else ""
            self.statusbar.showMessage(
                f"m/z: {mz:.4f}    Nearest: {pmz[idx]:.4f}  "
                f"Int: {disp[idx]:.0f}{suffix}{ruler_hint}"
            )

    # ── Zoom buttons ──────────────────────────────────────────────────────────
    def _on_zoom_axis_changed(self):
        sender = self.sender()
        if sender == self.btn_zoom_x and self.btn_zoom_x.isChecked():
            self.btn_zoom_y.setChecked(False)
        elif sender == self.btn_zoom_y and self.btn_zoom_y.isChecked():
            self.btn_zoom_x.setChecked(False)
        else:
            sender.setChecked(True)
        self._vb._zoom_x = self.btn_zoom_x.isChecked()
        self._vb._zoom_y = self.btn_zoom_y.isChecked()

    # ── View controls ─────────────────────────────────────────────────────────
    def reset_view(self):
        self.plot_widget.autoRange()

    def clear_all(self):
        for spec in self.spectra:
            if spec["plot_item"]: self.plot_widget.removeItem(spec["plot_item"])
        self.spectra.clear()
        self.annotations.clear()
        self.iso_peaks = None; self.iso_params = None
        if self._overlay_item:
            self.plot_widget.removeItem(self._overlay_item)
            self._overlay_item = None
        for line, text in self._ann_items:
            self.plot_widget.removeItem(line); self.plot_widget.removeItem(text)
        self._ann_items.clear()
        self._clear_ruler_items(); self._ruler_pts.clear()
        for item in self._fl_items:
            self.plot_widget.removeItem(item)
        self._fl_items.clear()
        self.free_labels.clear()
        self.scale_group.setVisible(False)
        self.iso_result_group.setVisible(False)
        self.btn_clear_overlay.setVisible(False)
        self._update_spec_list()
        self._refresh_annotation_list()
        self.statusbar.showMessage("No spectrum loaded")

    # ── Export ────────────────────────────────────────────────────────────────
    def _render_png_to_file(self, path):
        exp = pg.exporters.ImageExporter(self.plot_widget.plotItem)
        exp.parameters()["width"]  = 1600
        exp.parameters()["height"] = 800
        exp.export(path)

    def copy_to_clipboard(self):
        if not self.spectra:
            QMessageBox.warning(self, "No data", "Load a spectrum first."); return
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f: tmp = f.name
        try:
            self._render_png_to_file(tmp)
            QGuiApplication.clipboard().setImage(QImage(tmp))
            self.statusbar.showMessage("Spectrum copied to clipboard  ✓", 3000)
        except Exception as e:
            QMessageBox.critical(self, "Copy failed", str(e))
        finally:
            try: os.unlink(tmp)
            except: pass

    def export_image(self):
        if not self.spectra:
            QMessageBox.warning(self, "No data", "Load a spectrum first."); return
        path, _ = QFileDialog.getSaveFileName(self, "Save PNG", "spectrum.png", "PNG (*.png)")
        if not path: return
        try:
            self._render_png_to_file(path)
            self.statusbar.showMessage(f"Saved: {path}", 4000)
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e))

    def export_data(self):
        """Export peak data for all visible spectra as CSV or tab-delimited TXT/DAT.
        If an isotopic pattern overlay is active, the user is offered the option
        to append the theoretical peak list as a labelled block in the same file."""
        if not self.spectra:
            QMessageBox.warning(self, "No data", "Load a spectrum first."); return

        path, _ = QFileDialog.getSaveFileName(
            self, "Export Spectrum Data", "spectrum_data",
            "CSV — comma separated (*.csv);;"
            "Text — tab separated (*.txt);;"
            "Origin ASCII (*.dat)"
        )
        if not path: return

        # If a pattern is active, ask whether to include it
        include_pattern = False
        if self.iso_peaks and self.iso_params:
            # Build dialog manually so we can use object-identity comparison on the
            # clicked button — avoids any fragile StandardButton enum comparisons.
            msg = QMessageBox(self)
            msg.setWindowTitle("Isotopic Pattern")
            msg.setText(
                f"Include isotopic pattern data?\n"
                f"({self.iso_params['formula']}  z={self.iso_params['z']})"
            )
            msg.setIcon(QMessageBox.Icon.Question)
            yes_btn = msg.addButton("Include pattern", QMessageBox.ButtonRole.YesRole)
            no_btn  = msg.addButton("Spectra only",    QMessageBox.ButtonRole.NoRole)
            msg.setDefaultButton(yes_btn)
            msg.exec()
            include_pattern = (msg.clickedButton() is yes_btn)

        try:
            import csv as _csv

            # Determine delimiter from chosen format
            ext = os.path.splitext(path)[1].lower()
            delim = "\t" if ext in (".txt", ".dat") else ","

            visible = [s for s in self.spectra if s["visible"]]
            if not visible:
                QMessageBox.warning(self, "No data", "No visible spectra to export."); return

            # Pre-compute display intensities for each visible spectrum
            disps = [self._spec_disp(s) for s in visible]

            # Total row count must span the longest spectrum and (if requested)
            # the isotopic pattern, so shorter columns are padded with blanks.
            spec_max = max(len(s["mz"]) for s in visible)
            max_len  = max(spec_max, len(self.iso_peaks)) if include_pattern else spec_max

            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = _csv.writer(f, delimiter=delim)

                # ── Header row ────────────────────────────────────────────────
                header = []
                if len(visible) == 1:
                    header += ["m/z", "Intensity"]
                else:
                    for s in visible:
                        header += [f"m/z_{s['name']}", f"Int_{s['name']}"]
                if include_pattern:
                    header += ["m/z_pattern", "Rel.Int._pattern"]
                writer.writerow(header)

                # ── Data rows ─────────────────────────────────────────────────
                for i in range(max_len):
                    row = []
                    # Spectrum column(s)
                    for s, d in zip(visible, disps):
                        if i < len(s["mz"]):
                            row += [f"{s['mz'][i]:.6f}", f"{d[i]:.4f}"]
                        else:
                            row += ["", ""]
                    # Optional pattern columns
                    if include_pattern:
                        if i < len(self.iso_peaks):
                            mz, rel = self.iso_peaks[i]
                            row += [f"{mz:.5f}", f"{rel * 100:.4f}"]
                        else:
                            row += ["", ""]
                    writer.writerow(row)

            self.statusbar.showMessage(
                f"Exported: {path}" + ("  (pattern included)" if include_pattern else ""),
                4000
            )
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))

    def _update_status(self):
        if not self.spectra:
            self.statusbar.showMessage("No spectrum loaded"); return
        n = sum(len(s["mz"]) for s in self.spectra)
        mz_all = np.concatenate([s["mz"] for s in self.spectra])
        self.statusbar.showMessage(
            f"{len(self.spectra)} spectrum{'a' if len(self.spectra)>1 else ''}  |  "
            f"{n:,} points  |  "
            f"m/z {mz_all.min():.2f} – {mz_all.max():.2f}"
        )


# ── Entry point ───────────────────────────────────────────────────────────────
def main():
    app = QApplication(sys.argv)
    app.setApplicationName("MSview")
    window = MSView()
    window.setAcceptDrops(True)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
