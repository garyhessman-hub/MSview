"""
Isotope distribution calculator.
Pure Python — no external dependencies.
"""

ELEMENTS = {
    "H":  {"mono": 1.00782503207,  "avg": 1.00794,   "iso_m": [1.00782503207, 2.01410178],               "iso_p": [0.999885, 0.000115]},
    "C":  {"mono": 12.0,           "avg": 12.0107,   "iso_m": [12.0, 13.00335484],                       "iso_p": [0.9893, 0.0107]},
    "N":  {"mono": 14.0030740048,  "avg": 14.0067,   "iso_m": [14.0030740048, 15.0001089],               "iso_p": [0.9963, 0.0037]},
    "O":  {"mono": 15.99491461956, "avg": 15.9994,   "iso_m": [15.99491461956, 16.99913170, 17.9991610], "iso_p": [0.99757, 0.00038, 0.00205]},
    "P":  {"mono": 30.97376163,    "avg": 30.97376,  "iso_m": [30.97376163],                             "iso_p": [1.0]},
    "S":  {"mono": 31.97207100,    "avg": 32.065,    "iso_m": [31.97207100, 32.97145876, 33.96786690, 35.96708076], "iso_p": [0.9493, 0.0076, 0.0429, 0.0002]},
    "F":  {"mono": 18.99840322,    "avg": 18.9984,   "iso_m": [18.99840322],                             "iso_p": [1.0]},
    "Cl": {"mono": 34.96885268,    "avg": 35.453,    "iso_m": [34.96885268, 36.96590260],                "iso_p": [0.7576, 0.2424]},
    "Br": {"mono": 78.9183371,     "avg": 79.904,    "iso_m": [78.9183371, 80.9162897],                  "iso_p": [0.5069, 0.4931]},
    "I":  {"mono": 126.904468,     "avg": 126.90447, "iso_m": [126.904468],                              "iso_p": [1.0]},
    "Si": {"mono": 27.9769265325,  "avg": 28.0855,   "iso_m": [27.9769265325, 28.9764947, 29.9737702],  "iso_p": [0.92297, 0.04685, 0.03018]},
    "Na": {"mono": 22.9897692809,  "avg": 22.98977,  "iso_m": [22.9897692809],                          "iso_p": [1.0]},
    "K":  {"mono": 38.96370668,    "avg": 39.0983,   "iso_m": [38.96370668, 39.96399848, 40.96182576],  "iso_p": [0.93258, 0.00012, 0.06730]},
    "Li": {"mono": 7.01600455,     "avg": 6.941,     "iso_m": [6.01512279, 7.01600455],                  "iso_p": [0.0759, 0.9241]},
    "B":  {"mono": 11.0093054,     "avg": 10.811,    "iso_m": [10.0129370, 11.0093054],                  "iso_p": [0.199, 0.801]},
    "Se": {"mono": 79.9165218,     "avg": 78.96,     "iso_m": [73.9224764, 75.9192136, 76.9199140, 77.9173091, 79.9165218, 81.9166994], "iso_p": [0.0089, 0.0937, 0.0763, 0.2377, 0.4961, 0.0873]},
    "Fe": {"mono": 55.9349393,     "avg": 55.845,    "iso_m": [53.9396105, 55.9349393, 56.9353940, 57.9332756], "iso_p": [0.05845, 0.91754, 0.02119, 0.00282]},
    "Cu": {"mono": 62.9295975,     "avg": 63.546,    "iso_m": [62.9295975, 64.9277895],                  "iso_p": [0.6917, 0.3083]},
    "Zn": {"mono": 63.9291422,     "avg": 65.38,     "iso_m": [63.9291422, 65.9260334, 66.9271273, 67.9248442, 69.9253193], "iso_p": [0.4917, 0.2773, 0.0404, 0.1845, 0.0061]},
    "Ag": {"mono": 106.905097,     "avg": 107.868,   "iso_m": [106.905097, 108.904752],                  "iso_p": [0.51839, 0.48161]},
    "Au": {"mono": 196.966569,     "avg": 196.967,   "iso_m": [196.966569],                              "iso_p": [1.0]},
    "Pt": {"mono": 194.964791,     "avg": 195.084,   "iso_m": [191.961035, 193.962664, 194.964791, 195.964951, 197.967869], "iso_p": [0.00782, 0.32967, 0.33832, 0.25242, 0.07163]},
    "Sn": {"mono": 119.902195,     "avg": 118.710,   "iso_m": [111.904818, 113.902779, 114.903342, 115.901741, 116.902952, 117.901603, 118.903308, 119.902195, 121.903440, 123.905274], "iso_p": [0.0097, 0.0066, 0.0034, 0.1454, 0.0768, 0.2422, 0.0859, 0.3258, 0.0463, 0.0579]},
    "Hg": {"mono": 201.970643,     "avg": 200.59,    "iso_m": [195.965833, 197.966769, 198.968279, 199.968326, 200.970302, 201.970643, 203.973494], "iso_p": [0.0015, 0.0997, 0.1687, 0.2310, 0.1318, 0.2986, 0.0687]},
    "Mo": {"mono": 97.9054078,     "avg": 95.96,     "iso_m": [91.906810, 93.905088, 94.905841, 95.904679, 96.906021, 97.905407, 99.907477], "iso_p": [0.1453, 0.0915, 0.1584, 0.1667, 0.0960, 0.2439, 0.0982]},
    "Pd": {"mono": 105.903486,     "avg": 106.42,    "iso_m": [101.905609, 103.904036, 104.905085, 105.903486, 107.903892, 109.905153], "iso_p": [0.0102, 0.1114, 0.2233, 0.2733, 0.2646, 0.1172]},
    "Ru": {"mono": 101.905634,     "avg": 101.07,    "iso_m": [95.907598, 97.905287, 98.905939, 99.904219, 100.905582, 101.904350, 103.905430], "iso_p": [0.0554, 0.0187, 0.1276, 0.1260, 0.1706, 0.3155, 0.1862]},
    "Ni": {"mono": 57.9353429,     "avg": 58.693,    "iso_m": [57.9353429, 59.9307884, 60.9310560, 61.9283451, 63.9279660], "iso_p": [0.68077, 0.26223, 0.01140, 0.03634, 0.00926]},
}

PROTON = 1.00727646677


def parse_formula(formula: str) -> dict:
    """Parse a molecular formula string into {element: count} dict.
    Raises ValueError on invalid input."""
    formula = formula.strip()
    if not formula:
        raise ValueError("Empty formula")
    result = {}
    i = 0
    while i < len(formula):
        if not formula[i].isupper():
            raise ValueError(f"Expected element symbol at position {i}: '{formula[i]}'")
        el = formula[i]
        i += 1
        if i < len(formula) and formula[i].islower():
            el += formula[i]
            i += 1
        if el not in ELEMENTS:
            raise ValueError(f"Unknown element: '{el}'")
        num = ""
        while i < len(formula) and formula[i].isdigit():
            num += formula[i]
            i += 1
        result[el] = result.get(el, 0) + (int(num) if num else 1)
    return result


def format_formula(atoms: dict) -> str:
    """Format atom dict as Hill-order formula string."""
    order = ["C", "H", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
    out = ""
    for el in order:
        if el in atoms:
            out += el + (str(atoms[el]) if atoms[el] > 1 else "")
    for el, n in atoms.items():
        if el not in order:
            out += el + (str(n) if n > 1 else "")
    return out


def monoisotopic_mass(atoms: dict) -> float:
    return sum(ELEMENTS[el]["mono"] * n for el, n in atoms.items())


def average_mass(atoms: dict) -> float:
    return sum(ELEMENTS[el]["avg"] * n for el, n in atoms.items())


def _merge_dist(dist: list) -> list:
    """Merge isotope peaks within 0.003 Da of each other."""
    dist.sort(key=lambda d: d[0])
    out = []
    for mass, prob in dist:
        if out and abs(out[-1][0] - mass) < 0.003:
            total = out[-1][1] + prob
            merged_mass = (out[-1][0] * out[-1][1] + mass * prob) / total
            out[-1] = (merged_mass, total)
        else:
            out.append((mass, prob))
    return out


def isotope_distribution(atoms: dict, charge: int = 1) -> list:
    """
    Calculate isotope distribution for a formula at a given charge state.
    Returns list of (mz, relative_intensity) sorted by mz,
    where the most abundant peak = 1.0.
    """
    dist = [(0.0, 1.0)]

    for el, count in atoms.items():
        iso_m = ELEMENTS[el]["iso_m"]
        iso_p = ELEMENTS[el]["iso_p"]
        for _ in range(count):
            new_dist = []
            for mass, prob in dist:
                for im, ip in zip(iso_m, iso_p):
                    p = prob * ip
                    if p > 1e-9:
                        new_dist.append((mass + im, p))
            dist = _merge_dist(new_dist)

    # Normalise
    max_p = max(p for _, p in dist)
    dist = [(m, p / max_p) for m, p in dist if p / max_p > 0.001]

    # Apply charge — convert to m/z
    return sorted(
        [((mass + charge * PROTON) / charge, rel) for mass, rel in dist],
        key=lambda x: x[0]
    )


def gaussian_profile(peaks: list, fwhm: float, scale: float) -> tuple:
    """
    Build a Gaussian-profiled envelope from stick peaks.
    Returns (mz_array, intensity_array) as numpy-compatible lists.
    peaks: list of (mz, relative_intensity)
    fwhm: full width at half maximum in Da
    scale: absolute intensity for a peak with relative_intensity = 1.0
    """
    import math
    if not peaks:
        return [], []

    sigma = fwhm / 2.3548
    step = fwhm / 25.0
    pad = sigma * 5.0
    mz_min = peaks[0][0] - pad
    mz_max = peaks[-1][0] + pad

    mz_out, int_out = [], []
    m = mz_min
    while m <= mz_max:
        y = sum(
            rel * scale * math.exp(-0.5 * ((m - peak_mz) / sigma) ** 2)
            for peak_mz, rel in peaks
        )
        mz_out.append(round(m, 6))
        int_out.append(y)
        m += step

    return mz_out, int_out
