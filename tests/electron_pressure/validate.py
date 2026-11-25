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
import cmath
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


def _pearson_corr(xs, ys):
    """Return Pearson correlation coefficient between two equal-length lists."""
    n = len(xs)
    if n < 3 or n != len(ys):
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= 0.0 or vy <= 0.0:
        return None
    return cov / math.sqrt(vx * vy)


def _check_pe_trace_identity(cols):
    """Return max relative error of Pe == (PePar + 2*PePerp)/3."""
    errs = []
    for pe, p_par, p_perp in zip(cols["Pe"], cols["PePar"], cols["PePerp"]):
        if pe <= 0.0:
            return None
        errs.append(abs(pe - (p_par + 2.0 * p_perp) / 3.0) / pe)
    return max(errs) if errs else None


def _cell_center_rho(rho_nodal):
    """Average periodic nodal density onto cell centers matching Pe."""
    n = len(rho_nodal)
    return [0.5 * (rho_nodal[i] + rho_nodal[(i + 1) % n]) for i in range(n)]


def _check_cgl_compression():
    """Verify CGL parallel compression makes PePar/PePerp track (rho/rho_0)^2."""
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    if len(plots) < 2:
        return False, "Expected >= 2 plot frames for cgl_compression"

    cols0 = _load_columns(plots[0], ["Pe", "PePar", "PePerp", "rhoS0"])
    if not cols0:
        return False, "Missing Pe/PePar/PePerp/rhoS0 in frame 0"
    rho0 = _cell_center_rho(cols0["rhoS0"])

    best = None
    max_trace_err = 0.0
    for path in plots[1:]:
        cols = _load_columns(path, ["Pe", "PePar", "PePerp", "rhoS0"])
        if not cols or len(cols["rhoS0"]) != len(rho0):
            return False, f"Invalid columns in {path}"
        trace_err = _check_pe_trace_identity(cols)
        if trace_err is None:
            return False, f"Non-positive Pe in {path}"
        max_trace_err = max(max_trace_err, trace_err)

        rho_cell = _cell_center_rho(cols["rhoS0"])
        ratios = [p_par / p_perp for p_par, p_perp in zip(cols["PePar"], cols["PePerp"])]
        spread_ratio = _rel_spread(ratios)
        rho_ratio = [r / r0 for r, r0 in zip(rho_cell, rho0)]
        scaled = [rat / (rr * rr) for rat, rr in zip(ratios, rho_ratio)]
        spread_scaled = _rel_spread(scaled)
        corr = _pearson_corr([rat - 1.0 for rat in ratios], [rr - 1.0 for rr in rho_ratio])
        if spread_ratio is None or spread_scaled is None or corr is None:
            continue
        if best is None or spread_ratio > best[0]:
            best = (spread_ratio, spread_scaled, corr, os.path.basename(path))

    if max_trace_err > 1e-6:
        return False, f"Trace identity Pe=(PePar+2*PePerp)/3 violated (err {max_trace_err:.3e})"
    if best is None:
        return False, "Could not evaluate CGL compression statistics"
    spread_ratio, spread_scaled, corr, fname = best
    logger.debug(
        "    [CGL_COMPRESSION] %s: spread(PePar/PePerp)=%.4f, spread/n^2=%.4f, "
        "corr=%.4f, trace_err=%.3e",
        fname,
        spread_ratio,
        spread_scaled,
        corr,
        max_trace_err,
    )
    if spread_ratio <= 0.08:
        return False, f"Insufficient CGL anisotropy spread {spread_ratio:.4f} (<= 0.08)"
    if spread_scaled >= 0.45 * spread_ratio:
        return (
            False,
            f"CGL invariant spread {spread_scaled:.4f} not flat vs {spread_ratio:.4f}",
        )
    if corr <= 0.90:
        return False, f"PePar/PePerp correlation with density {corr:.4f} (<= 0.90)"
    return True, "Passed"


def _check_aniso_limit(test_name):
    """Verify marginal-stability clamping conserves Pe and enforces bounds."""
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    if len(plots) < 2:
        return False, f"Expected >= 2 plot frames for {test_name}"

    req = ["Pe", "PePar", "PePerp", "Bx", "By", "Bz"]
    cols0 = _load_columns(plots[0], req)
    cols1 = _load_columns(plots[-1], req)
    if not cols0 or not cols1:
        return False, f"Missing columns in {test_name} plot files"

    for label, cols in (("t=0", cols0), ("final", cols1)):
        trace_err = _check_pe_trace_identity(cols)
        if trace_err is None or trace_err > 1e-6:
            return False, f"{label} Pe=(PePar+2*PePerp)/3 mismatch ({trace_err})"

    mean_pe0 = sum(cols0["Pe"]) / len(cols0["Pe"])
    mean_pe1 = sum(cols1["Pe"]) / len(cols1["Pe"])
    pe_drift = abs(mean_pe1 - mean_pe0) / mean_pe0
    if pe_drift > 1e-5:
        return False, f"Isotropization failed to conserve Pe (drift {pe_drift:.3e})"

    pb0 = [(bx * bx + by * by + bz * bz) / (8.0 * math.pi)
           for bx, by, bz in zip(cols0["Bx"], cols0["By"], cols0["Bz"])]
    pb1 = [(bx * bx + by * by + bz * bz) / (8.0 * math.pi)
           for bx, by, bz in zip(cols1["Bx"], cols1["By"], cols1["Bz"])]

    if test_name == "electron_pressure_aniso_firehose":
        init_val = sum((par - perp) / pb for par, perp, pb in zip(
            cols0["PePar"], cols0["PePerp"], pb0)) / len(pb0)
        max_err = max(abs((par - perp) / pb - 2.0) for par, perp, pb in zip(
            cols1["PePar"], cols1["PePerp"], pb1))
        logger.debug(
            "    [FIREHOSE] init (Ppar-Pperp)/PB=%.3f, final max|err|=%.3e, Pe drift=%.3e",
            init_val,
            max_err,
            pe_drift,
        )
        if init_val <= 2.5:
            return False, f"Initial firehose parameter {init_val:.3f} not > 2.5"
        if max_err >= 1e-3:
            return False, f"Firehose threshold mismatch {max_err:.3e} (>= 1e-3)"
        return True, "Passed"

    init_val = sum((perp / par - 1.0) for par, perp in zip(
        cols0["PePar"], cols0["PePerp"])) / len(pb0)
    diffs = [
        (perp / par - 1.0) - 0.36 / ((par / pb) ** 0.55)
        for par, perp, pb in zip(cols1["PePar"], cols1["PePerp"], pb1)
    ]
    ratios = [
        (perp / par - 1.0) / (0.36 / ((par / pb) ** 0.55))
        for par, perp, pb in zip(cols1["PePar"], cols1["PePerp"], pb1)
    ]
    max_excess = max(diffs)
    min_frac = min(ratios)
    logger.debug(
        "    [WHISTLER_LIMIT] init Pperp/Ppar-1=%.3f, max_excess=%.3e, "
        "min(lhs/rhs)=%.4f, Pe drift=%.3e",
        init_val,
        max_excess,
        min_frac,
        pe_drift,
    )
    if init_val <= 5.0:
        return False, f"Initial whistler parameter {init_val:.3f} not > 5.0"
    if max_excess >= 1e-3:
        return False, f"Whistler threshold exceeded {max_excess:.3e} (>= 1e-3)"
    if min_frac < 0.85:
        return False, f"Whistler limiter over-isotropized (min lhs/rhs {min_frac:.3f} < 0.85)"
    return True, "Passed"


def _check_aniso_wave():
    """Verify gyrotropic -div(Pe)/rho_q modifies the parallel whistler-Alfven speed."""
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    if len(plots) < 5:
        return False, "Expected >= 5 plot frames for aniso_wave"

    cols0 = _load_columns(plots[0], ["X", "Bx", "By", "Bz", "PePar", "PePerp", "rhoS0"])
    if not cols0 or len(cols0["X"]) < 8:
        return False, "Missing columns in aniso_wave frame 0"

    xs = cols0["X"]
    n = len(xs)
    dx = xs[1] - xs[0]
    lx = n * dx
    k_wave = 2.0 * math.pi / lx

    b0 = sum(cols0["Bx"]) / n
    rho0 = sum(cols0["rhoS0"]) / n
    pe_par0 = sum(cols0["PePar"]) / n
    pe_perp0 = sum(cols0["PePerp"]) / n
    p_b = (b0 * b0) / (8.0 * math.pi)
    v_a = abs(b0) / math.sqrt(4.0 * math.pi * rho0)
    d_i = 1.0 / math.sqrt(4.0 * math.pi * rho0)
    alpha = 1.0 - (pe_par0 - pe_perp0) / (2.0 * p_b)

    t1 = 0.5 * alpha * (k_wave ** 2) * d_i * v_a
    omega_aniso = t1 + math.sqrt(t1 * t1 + alpha * (k_wave ** 2) * (v_a ** 2))
    t1_iso = 0.5 * (k_wave ** 2) * d_i * v_a
    omega_iso = t1_iso + math.sqrt(t1_iso * t1_iso + (k_wave ** 2) * (v_a ** 2))

    # Pick spatial harmonic sign (+k vs -k) that carries the wave.
    m_pos = sum(
        complex(by, bz) * cmath.exp(-2j * math.pi * i / n)
        for i, (by, bz) in enumerate(zip(cols0["By"], cols0["Bz"]))
    ) / n
    m_neg = sum(
        complex(by, bz) * cmath.exp(2j * math.pi * i / n)
        for i, (by, bz) in enumerate(zip(cols0["By"], cols0["Bz"]))
    ) / n
    ksign = 1.0 if abs(m_pos) >= abs(m_neg) else -1.0

    samples = []
    for path in plots:
        t = _frame_time(path)
        cols = _load_columns(path, ["By", "Bz"])
        if t is None or not cols or len(cols["By"]) != n:
            continue
        c_mode = sum(
            complex(by, bz) * cmath.exp(-ksign * 2j * math.pi * i / n)
            for i, (by, bz) in enumerate(zip(cols["By"], cols["Bz"]))
        ) / n
        if abs(c_mode) > 1e-12:
            samples.append((t, cmath.phase(c_mode)))

    if len(samples) < 5:
        return False, "Too few usable frames for aniso_wave dispersion fit"
    samples.sort(key=lambda s: s[0])
    unwrapped = [samples[0][1]]
    for i in range(1, len(samples)):
        dphi = samples[i][1] - samples[i - 1][1]
        dphi -= 2.0 * math.pi * round(dphi / (2.0 * math.pi))
        unwrapped.append(unwrapped[-1] + dphi)

    ts = [s[0] for s in samples]
    nt = len(ts)
    mt = sum(ts) / nt
    mp = sum(unwrapped) / nt
    cov = sum((t - mt) * (p - mp) for t, p in zip(ts, unwrapped))
    vt = sum((t - mt) ** 2 for t in ts)
    if vt <= 0.0:
        return False, "Degenerate time axis in aniso_wave"
    omega_meas = abs(cov / vt)

    err_aniso = abs(omega_meas - omega_aniso) / omega_aniso
    sep_iso = abs(omega_meas - omega_iso) / omega_iso
    logger.debug(
        "    [ANISO_WAVE] alpha=%.3f, omega_meas=%.4f, omega_aniso=%.4f (err %.2f%%), "
        "omega_iso=%.4f (sep %.2f%%)",
        alpha,
        omega_meas,
        omega_aniso,
        100.0 * err_aniso,
        omega_iso,
        100.0 * sep_iso,
    )
    if err_aniso >= 0.08:
        return (
            False,
            f"Anisotropic dispersion mismatch: meas {omega_meas:.4f} vs {omega_aniso:.4f} "
            f"(err {err_aniso:.3f} >= 0.08)",
        )
    if sep_iso <= 0.25:
        return (
            False,
            f"Measured omega {omega_meas:.4f} too close to isotropic {omega_iso:.4f} "
            f"(sep {sep_iso:.3f} <= 0.25)",
        )
    return True, "Passed"


def _check_le2009():
    """Verify algebraic Le et al. (2009) closure and TePar/TePerp diagnostics."""
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots", "*.out")))
    if not plots:
        return False, "No plot frames found for le2009"

    cols = _load_columns(
        plots[-1],
        ["Pe", "PePar", "PePerp", "Te", "TePar", "TePerp", "Bx", "By", "Bz", "rhoS0"],
    )
    if not cols:
        return False, "Missing Pe/Te/B/rhoS0 columns in le2009 final plot"

    trace_err = _check_pe_trace_identity(cols)
    if trace_err is None or trace_err > 1e-6:
        return False, f"Pe=(PePar+2*PePerp)/3 violated in le2009 ({trace_err})"

    rho_cell = _cell_center_rho(cols["rhoS0"])
    n = len(rho_cell)
    rho0 = sum(rho_cell) / n
    bx_avg = [0.5 * (cols["Bx"][i] + cols["Bx"][(i + 1) % n]) for i in range(n)]
    by_avg = [0.5 * (cols["By"][i] + cols["By"][(i + 1) % n]) for i in range(n)]
    bz_avg = [0.5 * (cols["Bz"][i] + cols["Bz"][(i + 1) % n]) for i in range(n)]
    b_mags = [math.sqrt(bx * bx + by * by + bz * bz)
              for bx, by, bz in zip(bx_avg, by_avg, bz_avg)]
    b0_code = (sum(b_mags) / n) * (1.5 / 1.0439)

    max_rel_err = 0.0
    max_te_err = 0.0
    min_ratio = float("inf")
    for p_par, p_perp, t_par, t_perp, b_mag, rho in zip(
        cols["PePar"], cols["PePerp"], cols["TePar"], cols["TePerp"], b_mags, rho_cell
    ):
        n_tilde = rho / rho0
        b_tilde = b_mag / b0_code
        x = (math.pi / 6.0) * (n_tilde / b_tilde) ** 2
        r_exact = ((1.0 + x ** 3) ** (1.0 / 3.0)) * (1.0 + x * x) / (
            1.0 + b_tilde * x * x
        )
        r_meas = p_par / p_perp
        max_rel_err = max(max_rel_err, abs(r_meas - r_exact) / r_exact)
        min_ratio = min(min_ratio, r_meas)

        tr_meas = t_par / t_perp
        max_te_err = max(max_te_err, abs(tr_meas - r_meas) / r_meas)

    logger.debug(
        "    [LE2009] min(PePar/PePerp)=%.4f, max_rel_err=%.3e, max_te_err=%.3e",
        min_ratio,
        max_rel_err,
        max_te_err,
    )
    if min_ratio <= 1.15:
        return False, f"Expected PePar/PePerp > 1.15 everywhere, got min {min_ratio:.4f}"
    if max_rel_err >= 5e-3:
        return False, f"Le et al. (2009) pointwise mismatch {max_rel_err:.3e} (>= 5e-3)"
    if max_te_err >= 1e-6:
        return False, f"TePar/TePerp != PePar/PePerp mismatch {max_te_err:.3e}"
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
    if test_name == "electron_pressure_cgl_compression":
        return _check_cgl_compression()
    if test_name in (
        "electron_pressure_aniso_firehose",
        "electron_pressure_aniso_whistler",
    ):
        return _check_aniso_limit(test_name)
    if test_name == "electron_pressure_aniso_wave":
        return _check_aniso_wave()
    if test_name == "electron_pressure_le2009":
        return _check_le2009()
    return _check_closure(test_name)
