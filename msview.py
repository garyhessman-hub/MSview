# ─────────────────────────────────────────────────────────────────────────────
# MSview patch — Bruker DataAnalysis ASCII support
# ─────────────────────────────────────────────────────────────────────────────
# Two changes to msview.py:
#
#   1. PASTE the _parse_bruker_ascii function below alongside the other
#      file parsers (next to _parse_bruker_xml is the natural home).
#
#   2. REPLACE the final 'else' branch of load_spectrum so that .ascii files
#      try the Bruker ASCII parser first, then fall back to _parse_text.
#
# That's it. No other code paths change.
# ─────────────────────────────────────────────────────────────────────────────


# ══════════════════════════════════════════════════════════════════════════════
# CHANGE 1 — Add this function next to _parse_bruker_xml
# ══════════════════════════════════════════════════════════════════════════════

def _parse_bruker_ascii(path):
    """
    Parse a Bruker DataAnalysis ASCII export.

    The file is logically one comma-separated line with the structure

        RT, polarity, ionization, MS-level, ?, mode, mass_range, count,
        "m/z intensity", "m/z intensity", ...

    e.g.  0.0864333,+,ESI,ms1,-,profile,44.9381-2005.3472,292352,44.9381 0,...

    The exporter inserts hard newlines at arbitrary character positions
    (sometimes mid-value), so we join the file into one logical line before
    parsing. Raises ValueError on any non-Bruker content so the caller can
    fall back to the generic two-column reader.

    Returns (mz, intensity, mode).
    """
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()

    # Strip the exporter's arbitrary hard wraps — the file is one logical line.
    text = re.sub(r"\s*\n\s*", "", text).strip()

    fields = text.split(",")
    if len(fields) < 9:
        raise ValueError("Bruker ASCII: too few fields.")

    # Header signature — strict enough to avoid false positives on plain
    # two-column ASCII, lenient enough to accept format variations.
    try:
        float(fields[0])   # retention time
    except ValueError:
        raise ValueError("Bruker ASCII: field 0 is not a retention time.")

    pol  = fields[1].strip()
    ms   = fields[3].strip().lower()
    mode = fields[5].strip().lower()
    mzr  = fields[6].strip()

    if pol not in ("+", "-"):
        raise ValueError("Bruker ASCII: polarity field not '+' or '-'.")
    if not re.match(r"ms\d*$", ms):
        raise ValueError("Bruker ASCII: MS-level field not 'msN'.")
    if mode not in ("profile", "centroid", "line", "stick"):
        raise ValueError("Bruker ASCII: mode field not profile/centroid.")
    if not re.match(r"[\d.]+\s*-\s*[\d.]+$", mzr):
        raise ValueError("Bruker ASCII: no mass range 'low-high' in field 6.")

    detected_mode = "profile" if "profile" in mode else "centroid"

    # The remaining fields are space-separated 'm/z intensity' pairs.
    mz_list, int_list = [], []
    for field in fields[8:]:
        parts = field.split()
        if len(parts) >= 2:
            try:
                m = float(parts[0])
                i = float(parts[1])
                if m > 0:
                    mz_list.append(m)
                    int_list.append(abs(i))
            except ValueError:
                pass

    if not mz_list:
        raise ValueError("Bruker ASCII: no peak data after header.")

    return (np.asarray(mz_list,  np.float64),
            np.asarray(int_list, np.float64),
            detected_mode)


# ══════════════════════════════════════════════════════════════════════════════
# CHANGE 2 — Update the final 'else' branch of load_spectrum
# ══════════════════════════════════════════════════════════════════════════════
#
# FIND this block at the bottom of load_spectrum:
#
#     else:
#         mz, intensity = _parse_text(path)
#         return mz, intensity, _detect_mode_text(mz)
#
# REPLACE it with:

    else:
        # For .ascii files, try the Bruker DataAnalysis exporter format first.
        # It looks like a text file but has a metadata header that confuses
        # the generic two-column reader. If the signature doesn't match,
        # fall back to plain text.
        if ext == ".ascii":
            try:
                return _parse_bruker_ascii(path)
            except ValueError:
                pass
        mz, intensity = _parse_text(path)
        return mz, intensity, _detect_mode_text(mz)
