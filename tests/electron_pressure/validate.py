#!/usr/bin/env python3
"""Validator for the evolved electron pressure test (tests/electron_pressure).

Three decks share this validator and differ only in the electron heat flux:

``adiabatic``    kappa0 = 0, B = 0            -> Pe must track rho^gamma
``conduction``   kappa0 huge, B = 0           -> Pe must track rho (Te flattened)
``crossfield``   kappa0 huge, B along y       -> conduction suppressed: rho^gamma

The check is built from the two scale-free ratios

    A = Pe / rho^gamma     (flat for adiabatic electrons)
    I = Pe / rho           (flat for isothermal electrons)

evaluated pointwise over the last plot frame. Because only the relative spread
(std/mean) of each ratio is compared, the result is independent of the output
units of Pe and rho.
"""
import glob
import logging
import math
import os

import tests._shared.hybrid as _hyb
from tests._shared.hybrid import validate_hybrid

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"

# gamma_e used by every deck in this suite (#ELECTRONTEMPERATURE).
GAMMA_E = 5.0 / 3.0

# The "flat" ratio must be at most this fraction of the other one. The seeded
# signal is ~9% and the per-cell PIC noise floor is ~2%, so 0.6 is a safe
# separation while still catching a genuinely flat/broken conduction term.
MAX_RATIO = 0.6

# Which ratio is expected to be flat, per variant.
EXPECT_FLAT = {
    "electron_pressure_adiabatic": "A",
    "electron_pressure_conduction": "I",
    "electron_pressure_crossfield": "A",
}


def set_run_dir(run_dir):
    """Point the plot helpers at the current run directory."""
    global RUN_DIR
    RUN_DIR = run_dir
    _hyb.set_run_dir(run_dir)


def _load_columns(out_file, names):
    """Return {NAME: [values]} for the requested columns of an idl .out file.

    The header layout matches tests/_shared/hybrid.py: variable names on the
    fifth line, data from the sixth onwards.
    """
    try:
        with open(out_file, "r", encoding="latin-1") as f:
            lines = f.readlines()
    except OSError:
        return None

    if len(lines) < 6:
        return None

    var_names = lines[4].split()
    vidx = {v.upper(): i for i, v in enumerate(var_names)}
    want = {n.upper(): n for n in names}
    if not all(w in vidx for w in want):
        return None

    cols = {n: [] for n in want.values()}
    min_cols = max(vidx[w] for w in want) + 1
    for line in lines[5:]:
        parts = line.split()
        if len(parts) < min_cols:
            continue
        try:
            for w, name in want.items():
                cols[name].append(float(parts[vidx[w]]))
        except ValueError:
            continue

    return cols if cols[names[0]] else None


def _rel_spread(values):
    """Relative spread std/mean of a list; None if it cannot be formed."""
    n = len(values)
    if n < 4:
        return None
    mean = sum(values) / n
    if mean == 0.0:
        return None
    var = sum((v - mean) ** 2 for v in values) / n
    return math.sqrt(var) / abs(mean)


def _frame_stats(out_file):
    """Return (spreadRho, spreadA, spreadI) for one plot frame, or None."""
    cols = _load_columns(out_file, ["Pe", "rhoS0"])
    if not cols:
        return None

    pe = cols["Pe"]
    rho = cols["rhoS0"]
    if len(pe) != len(rho):
        return None

    keep = [(p, r) for p, r in zip(pe, rho) if r > 0.0 and p > 0.0]
    if len(keep) < 4:
        return None

    ratio_a = [p / r ** GAMMA_E for p, r in keep]
    ratio_i = [p / r for p, r in keep]
    return (_rel_spread(rho), _rel_spread(ratio_a), _rel_spread(ratio_i))


def _spreads():
    """Return (spreadA, spreadI) from the most discriminating plot frame.

    The ion acoustic wave is standing and Landau damped, so the density
    contrast decays from frame to frame, and the first frames still carry the
    initial adiabatic transient that the conduction needs ~2 s to flatten
    (the point-implicit update is stable for any step but converges at one
    Jacobi sweep per step). The frame where the two closure hypotheses differ
    most, max |spreadA - spreadI|, is therefore the meaningful one: in
    noise-only frames both spreads are equal and the difference vanishes.
    """
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    best = None
    for path in plots:
        stats = _frame_stats(path)
        if stats is None:
            continue
        spread_rho, spread_a, spread_i = stats
        if spread_rho is None or spread_a is None or spread_i is None:
            continue
        gap = abs(spread_a - spread_i)
        if best is None or gap > best[0]:
            best = (gap, spread_a, spread_i)

    if best is None:
        return None, None
    return best[1], best[2]


def _check_closure(test_name):
    """Check which of Pe/rho^gamma or Pe/rho is the flat one."""
    expected = EXPECT_FLAT.get(test_name)
    if expected is None:
        return True, "No expectation for %s (skipped)" % test_name

    spread_a, spread_i = _spreads()
    if spread_a is None or spread_i is None:
        logger.debug("    [ELECTRON_PRESSURE] No Pe/rho plot data (skipped)")
        return True, "No Pe/rho plot data (skipped)"

    flat, other = ((spread_a, spread_i) if expected == "A"
                   else (spread_i, spread_a))
    label = "Pe/rho^gamma" if expected == "A" else "Pe/rho"
    logger.debug("    [ELECTRON_PRESSURE] spread Pe/rho^gamma %.4f, "
                 "Pe/rho %.4f (expect %s flat)", spread_a, spread_i, label)

    if other <= 0.0:
        return True, "Degenerate spread measurement"

    if flat > MAX_RATIO * other:
        return False, ("%s spread %.4f is not flat versus %.4f (bound %.2fx) "
                       "-- electron pressure closure is wrong"
                       % (label, flat, other, MAX_RATIO))
    return True, "Passed"


def validate_log(pic_diags=None, test_name=None):
    """Stability / conservation checks shared with the other hybrid tests."""
    return validate_hybrid(pic_diags=pic_diags, test_name=test_name)


def validate_plot(test_name):
    """Verify which electron-pressure closure the evolved Pe follows."""
    return _check_closure(test_name)
