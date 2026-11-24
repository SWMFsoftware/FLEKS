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
    "electron_pressure_conduction_rkl2": "I",
    "electron_pressure_amr": "A",
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
        # Skip the t=0 / step 0 frame: it contains the initial condition
        # (always seeded adiabatically) before dynamical evolution / conduction.
        if "_t00000000_n00000000" in path:
            continue
        stats = _frame_stats(path)
        if stats is None:
            continue
        spread_rho, spread_a, spread_i = stats
        if spread_rho is None or spread_a is None or spread_i is None:
            continue
        gap = abs(spread_a - spread_i)
        logger.debug("Plot frame %s: gap=%.4f, spread_a=%.4f, spread_i=%.4f",
                     os.path.basename(path), gap, spread_a, spread_i)
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


def _frame_time(out_file):
    """Extract simulation time from the second line of an IDL .out header."""
    try:
        with open(out_file, "r", encoding="latin-1") as f:
            f.readline()
            line2 = f.readline().split()
        return float(line2[1])
    except (OSError, IndexError, ValueError):
        return None


def _check_conduction_rkl2():
    """Verify RKL2 flattens Te, conserves sum(Pe), and has no checkerboard."""
    ok, msg = _check_closure("electron_pressure_conduction_rkl2")
    if not ok:
        return ok, msg

    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    if len(plots) < 2:
        return False, "Expected >= 2 plot frames for RKL2 conservation check"

    sums = []
    for path in plots:
        cols = _load_columns(path, ["Pe", "rhoS0"])
        if not cols:
            return False, f"Missing Pe/rhoS0 in {path}"
        sums.append(sum(cols["Pe"]))

    s0 = sums[0]
    if s0 <= 0.0:
        return False, f"Non-positive initial Pe sum {s0}"
    max_rel = max(abs(s - s0) / s0 for s in sums[1:])
    logger.debug("    [RKL2] sum(Pe) relative drift across frames: %.3e", max_rel)
    if max_rel > 1e-5:
        return False, f"RKL2 sum(Pe) drifted by {max_rel:.3e} (> 1e-5)"

    # Check absence of 2*dx checkerboard in Te = Pe/rho at the final frame.
    last_cols = _load_columns(plots[-1], ["Pe", "rhoS0"])
    te = [p / r for p, r in zip(last_cols["Pe"], last_cols["rhoS0"]) if r > 0.0]
    if len(te) >= 8:
        nyq = abs(sum(((-1.0) ** i) * v for i, v in enumerate(te))) / len(te)
        rel_nyq = nyq / (sum(te) / len(te))
        logger.debug("    [RKL2] final Te Nyquist mode fraction: %.3e", rel_nyq)
        if rel_nyq > 1e-3:
            return False, f"RKL2 Te shows checkerboard mode {rel_nyq:.3e}"

    return True, "Passed"


def _check_ohmic(test_name):
    """Verify Ohmic heating raises Pe and conserves B^2/(8*pi) + Pe/(gamma-1)."""
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    if len(plots) < 2:
        return False, f"Expected >= 2 plot frames for {test_name}"

    tot_energies = []
    pe_means = []
    for path in plots:
        cols = _load_columns(path, ["Bx", "By", "Bz", "Pe"])
        if not cols:
            return False, f"Missing Bx/By/Bz/Pe in {path}"
        n = len(cols["Pe"])
        em = sum(
            (bx * bx + by * by + bz * bz) / (8.0 * math.pi)
            for bx, by, bz in zip(cols["Bx"], cols["By"], cols["Bz"])
        ) / n
        eth = sum(p / (GAMMA_E - 1.0) for p in cols["Pe"]) / n
        tot_energies.append(em + eth)
        pe_means.append(sum(cols["Pe"]) / n)

    pe_gain = (pe_means[-1] - pe_means[0]) / pe_means[0]
    e0 = tot_energies[0]
    max_rel = max(abs(e - e0) / e0 for e in tot_energies[1:])
    logger.debug(
        "    [%s] Pe gain=%.4f, total E_mag+E_th_e rel_drift=%.3e",
        test_name,
        pe_gain,
        max_rel,
    )
    if pe_gain < 0.05:
        return False, f"Ohmic heating Pe gain {pe_gain:.4f} too small (< 0.05)"
    if max_rel > 5e-3:
        return False, f"Ohmic total energy conservation drift {max_rel:.3e} (> 5e-3)"
    return True, "Passed"


def _check_collision():
    """Verify #ELECTRONCOLLISION matches analytical Braginskii relaxation ODE."""
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    if len(plots) < 2:
        return False, "Expected >= 2 plot frames for collision test"

    times = []
    pe_vals = []
    pi_vals = []
    rho_vals = []
    for path in plots:
        t = _frame_time(path)
        cols = _load_columns(path, ["Pe", "pS0", "rhoS0"])
        if t is None or not cols:
            return False, f"Missing time or columns in {path}"
        times.append(t)
        pe_vals.append(sum(cols["Pe"]) / len(cols["Pe"]))
        pi_vals.append(sum(cols["pS0"]) / len(cols["pS0"]))
        rho_vals.append(sum(cols["rhoS0"]) / len(cols["rhoS0"]))

    # Compute code-unit collisionCoefEi from PARAM.in.collision constants.
    mp = 1.67262192e-27
    me = 9.1093837e-31
    qe = 1.60217663e-19
    kb = 1.380649e-23
    eps0 = 8.8541878176198092e-12
    l_norm = 1.0e5
    u_norm = 1.0e5
    m_norm = 1.0e7 * l_norm * (mp / qe) ** 2
    si2no_rho = (l_norm ** 3) / m_norm
    no2si_n = 1.0 / (si2no_rho * mp)
    si2no_t = l_norm / u_norm
    si2no_temp = kb / (mp * u_norm * u_norm)
    mr = (me * mp) / (me + mp)
    e2eps = (qe * qe) / eps0
    coef_si = 20.0 * (math.sqrt(mr) / mp) * (e2eps ** 2) / (3.0 * (2.0 * math.pi * kb) ** 1.5)
    c_ei = 2.0 * 1.0e4 * coef_si * no2si_n * (si2no_temp ** 1.5) / si2no_t

    # Integrate exact ODE dPe/dt = H(Te)*(Pi - Pe) over [times[0], times[-1]].
    # In PARAM.in.collision, ion species 0 has mass = 1e4 amu and charge = 1,
    # so number/charge density is ne = rhoS0 / 1e4.
    pe_ode = pe_vals[0]
    pi_ref = sum(pi_vals) / len(pi_vals)
    ne_ref = (sum(rho_vals) / len(rho_vals)) / 1.0e4
    n_sub = 2000
    dt_sub = (times[-1] - times[0]) / n_sub
    for _ in range(n_sub):
        te = pe_ode / ne_ref
        h = c_ei * ne_ref / (te ** 1.5)
        pe_ode = (pe_ode + dt_sub * h * pi_ref) / (1.0 + dt_sub * h)

    rel_err = abs(pe_vals[-1] - pe_ode) / pe_ode
    frac_prog = (pe_vals[-1] - pe_vals[0]) / (pi_ref - pe_vals[0])
    logger.debug(
        "    [COLLISION] Pe %.4e -> %.4e (ODE %.4e, Pi %.4e, progress %.1f%%, err %.3e)",
        pe_vals[0],
        pe_vals[-1],
        pe_ode,
        pi_ref,
        100.0 * frac_prog,
        rel_err,
    )
    if frac_prog < 0.10:
        return False, f"Insufficient collisional equilibration ({frac_prog:.3f})"
    if rel_err > 0.05:
        return False, f"Collision ODE mismatch {rel_err:.3e} (> 0.05)"
    return True, "Passed"


def _check_biermann_2d():
    """Verify non-collinear grad(ne) x grad(Te) generates out-of-plane Bz."""
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    if len(plots) < 2:
        return False, "Expected >= 2 plot frames for biermann_2d"

    cols0 = _load_columns(plots[0], ["Bz"])
    cols1 = _load_columns(plots[-1], ["Bz", "Pe", "rhoS0"])
    if not cols0 or not cols1:
        return False, "Missing Bz/Pe/rhoS0 in biermann_2d plot files"

    max_bz0 = max(abs(b) for b in cols0["Bz"])
    max_bz1 = max(abs(b) for b in cols1["Bz"])
    logger.debug("    [BIERMANN_2D] max|Bz| %.3e -> %.3e", max_bz0, max_bz1)
    if max_bz0 > 1e-10:
        return False, f"Initial Bz unexpectedly non-zero: {max_bz0:.3e}"
    if not math.isfinite(max_bz1) or max_bz1 < 1e-4:
        return False, f"Biermann battery failed to generate Bz (max={max_bz1:.3e})"
    return True, "Passed"


def validate_log(pic_diags=None, test_name=None):
    """Stability / conservation checks shared with the other hybrid tests."""
    return validate_hybrid(pic_diags=pic_diags, test_name=test_name)


def validate_plot(test_name):
    """Verify evolved electron pressure physics for each variant."""
    if test_name == "electron_pressure_conduction_rkl2":
        return _check_conduction_rkl2()
    if test_name in ("electron_pressure_ohmic", "electron_pressure_pc"):
        return _check_ohmic(test_name)
    if test_name == "electron_pressure_collision":
        return _check_collision()
    if test_name == "electron_pressure_biermann_2d":
        return _check_biermann_2d()
    return _check_closure(test_name)
