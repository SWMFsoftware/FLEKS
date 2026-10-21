#!/usr/bin/env python3
"""Validation for the hybrid-PIC upwind (#HYBRIDUPWIND) dissipation test.

Variants (selected by base_name):
  * hybrid_upwind       -- first-order upwind with a fixed signal speed in a
                           frozen plasma: the seeded mode decays at the
                           analytic discrete rate.
  * hybrid_upwind_mc    -- MC-limited reconstruction: the rate lies between
                           the unlimited second-order and the first-order rate.
  * hybrid_upwind_divb  -- noisy 2D plasma, Hall term, whistler signal speed:
                           the discrete div(B) does not change.
  * hybrid_upwind_whistler -- seeded n = 4 whistler (k d_i = 3.9, 16 cells per
                           wavelength) with kinetic ions and the default signal
                           speed: the frequency stays on the Hall branch and the
                           damping per radian stays in a narrow band.
  * hybrid_upwind_whistler_sensor -- the same with upwindSensorFloor = 0.1: the
                           damping drops by an order of magnitude.
"""
import cmath
import logging
import math
import os

import numpy as np

from tests._shared import run_dir as _run_dir
from tests._shared.validators import validate_hybrid_energy_bounded
from tests.hyper_resistivity.validate import (_block_numbers, _fit_decay,
                                              _named_number)

logger = logging.getLogger(__name__)

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(os.path.dirname(TEST_DIR))
RUN_DIR = "run_test"  # default; overridden by set_run_dir() if called


def set_run_dir(run_dir):
    """Mirror the runner's RUN_DIR into this module."""
    global RUN_DIR
    RUN_DIR = run_dir
    _run_dir.set_run_dir(run_dir)


TOL = 0.005           # Relative tolerance of the first-order damping rate
MAX_RESIDUAL = 1e-6   # Fit residual bound (first order is exactly linear)
MC_TOL = 0.02         # Slack below the unlimited (Fromm) rate for MC
MIN_SAMPLES = 20      # Minimum log samples required
# |div(B)_last - div(B)_first| relative to max|B|/dx. The update changes B by
# a curl of the nodal E, whose discrete divergence is zero up to round-off.
DIVB_TOL = 1e-10
# whistler variant: damping per radian gamma/omega of the seeded whistler.
# Measured 0.0022 (no dissipation: 0; hyper-resistivity with etaHyperCh = 0.01:
# 0.0007; upwindCoefWhistler = 0.25: 0.020, = 1: 0.070). The lower bound shows
# the term acts; the upper one catches a much more diffusive scheme.
WHISTLER_DAMP_MIN = 0.0005
WHISTLER_DAMP_MAX = 0.005
# whistler_sensor: measured 0.00015 with upwindSensorFloor = 0.1 (0.00007 at
# 0.05, 0.00036 at 0.2); the band ends where the unscaled one starts.
SENSOR_DAMP_MIN = 0.00005
SENSOR_DAMP_MAX = 0.0005
# Frequency vs the Hall whistler branch at k d_i = k (the decks set d_i ~ 1);
# every whistler deck measures +2.4 % from that normalization.
WHISTLER_OMEGA_TOL = 0.04


def _variant(test_name):
    name = test_name or ""
    if name.endswith("_divb"):
        return "divb"
    if name.endswith("_whistler"):
        return "whistler"
    if name.endswith("_whistler_sensor"):
        return "whistler_sensor"
    if name.endswith("_mc"):
        return "mc"
    return "upwind1"


def _param_text():
    """The PARAM.in of the run (the runner writes the active one to RUN_DIR)."""
    path = os.path.join(REPO_DIR, RUN_DIR, "PARAM.in")
    if not os.path.isfile(path):
        path = os.path.join(TEST_DIR, "PARAM.in")
    with open(path, "r") as f:
        return f.read()


def _rk4_rate(lam, dt):
    """Damping rate per unit time of RK4 applied to dB/dt = lam*B."""
    z = lam * dt
    amp = abs(1.0 + z + z * z / 2.0 + z ** 3 / 6.0 + z ** 4 / 24.0)
    return -math.log(amp) / dt if 0.0 < amp < 1.0 else -lam.real


def _analytic_rates(text):
    """Return (gamma_first_order, gamma_unlimited_second_order) of the mode.

    The upwind term adds (a/2) * jump to the nodal E, and the node-to-center
    curl turns it into dB_j/dt = (a/2dx) (jump_{j+1/2} - jump_{j-1/2}) for
    each transverse component. With B_j ~ exp(i j theta):
      first order:  jump = B_j - B_{j-1}
      unlimited:    jump = (B_j - s_j/2) - (B_{j-1} + s_{j-1}/2),
                    s_j = (B_{j+1} - B_{j-1})/2  (central slope).
    """
    l_norm, u_norm = _block_numbers(text, "#NORMALIZATION", 2)
    x_lo, x_hi = _block_numbers(text, "#GEOMETRY", 2)
    n_cell = _block_numbers(text, "#NCELL", 1)[0]
    mode = _named_number(text, "#WAVEIC", "waveMode")
    speed = _named_number(text, "#HYBRIDUPWIND", "upwindFixedSpeed")
    dt = _named_number(text, "#TIMESTEPPING", "dt")
    if None in (mode, speed, dt) or n_cell <= 0 or speed <= 0:
        return None, None

    dx = abs(x_hi - x_lo) / n_cell
    theta = 2.0 * math.pi * mode / n_cell
    e = cmath.exp(-1j * theta)  # B_{j-1} / B_j

    jump1 = 1.0 - e
    s_j = (cmath.exp(1j * theta) - e) / 2.0
    s_jm1 = (1.0 - e * e) / 2.0
    jump2 = (1.0 - 0.5 * s_j) - (e + 0.5 * s_jm1)

    pref = speed / (2.0 * dx) * (cmath.exp(1j * theta) - 1.0)
    return _rk4_rate(pref * jump1, dt), _rk4_rate(pref * jump2, dt)


def validate_log(pic_diags=None, test_name=None):
    """Decay-rate checks (frozen-plasma variants) or bounded energy (divb)."""
    variant = _variant(test_name)
    if not pic_diags:
        return False, "No energy log"
    if variant in ("divb", "whistler", "whistler_sensor"):
        return validate_hybrid_energy_bounded(
            pic_diags, growth_max=10.0, test_label=f"Upwind {variant}",
            logger=logger)

    if len(pic_diags) < MIN_SAMPLES:
        return False, "Too few energy-log samples"
    t = [row["time"] for row in pic_diags]
    eb = [row["Eb"] for row in pic_diags]
    if any(not math.isfinite(v) for v in eb):
        return False, "Eb is not finite"
    if any(eb[i + 1] > eb[i] for i in range(len(eb) - 1)):
        return False, "Eb is not monotonically decreasing"

    g1, g2 = _analytic_rates(_param_text())
    if g1 is None:
        return False, "Could not derive the analytic rates from PARAM.in"

    gamma, _, resid = _fit_decay(t, eb)
    if gamma is None or gamma <= 0.0:
        return False, "Seeded mode did not decay"
    logger.debug("    [UPWIND] %s: gamma_fit = %.6g, first order = %.6g, "
                 "unlimited second order = %.6g, residual = %s",
                 variant, gamma, g1, g2,
                 "n/a" if resid is None else "%.2e" % resid)

    if variant == "upwind1":
        if abs(gamma / g1 - 1.0) > TOL:
            return False, ("measured decay rate %.6g differs from the analytic "
                           "first-order rate %.6g by %.2f%% (tolerance %.1f%%)"
                           % (gamma, g1, 100.0 * (gamma / g1 - 1.0),
                              100.0 * TOL))
        if resid is not None and resid > MAX_RESIDUAL:
            return False, ("decay is not a single exponential: residual %.2e "
                           "(tolerance %.1e)" % (resid, MAX_RESIDUAL))
        return True, "Passed"

    if not (1.0 - MC_TOL) * g2 <= gamma <= g1:
        return False, ("MC decay rate %.6g outside [%.6g (unlimited), %.6g "
                       "(first order)]" % (gamma, g2, g1))
    return True, "Passed"


def _whistler_fit(text):
    """Frequency and decay rate of the seeded whistler from the plot frames.

    c(t) = <(By + i hand Bz) exp(-i k x)> holds the seeded whistler and, more
    weakly, the slow branch of the same helicity. The frequency is the phase
    slope of c; the decay rate is fitted to the envelope of c demodulated by
    that frequency and averaged over one period, which removes the slow branch.
    Returns (omega, gamma, n_frames) in code units, or None.
    """
    x_lo, x_hi = _block_numbers(text, "#GEOMETRY", 2)
    mode = _named_number(text, "#WAVEIC", "waveMode")
    dt = _named_number(text, "#TIMESTEPPING", "dt")
    if None in (mode, dt):
        return None
    lx = abs(x_hi - x_lo)
    k = 2.0 * math.pi * mode / lx
    ts, cp, cm = [], [], []
    for f in _run_dir.plot_files():
        vidx, rows = _run_dir._read_out_file(f)
        a = np.array(rows)
        x = a[:, vidx["X"]]
        keep = x < x_lo + lx - 1e-6        # drop a periodic end node
        e = np.exp(-1j * k * x[keep])
        by, bz = a[keep, vidx["BY"]], a[keep, vidx["BZ"]]
        cp.append(np.mean((by + 1j * bz) * e))
        cm.append(np.mean((by - 1j * bz) * e))
        ts.append(int(os.path.basename(f).split("_n")[-1][:8]) * dt)
    if len(ts) < MIN_SAMPLES:
        return None
    t = np.array(ts)
    c = np.array(cp if abs(cp[0]) > abs(cm[0]) else cm)   # seeded helicity
    omega = np.polyfit(t, np.unwrap(np.angle(c)), 1)[0]
    m = max(1, int(round(2.0 * math.pi / abs(omega) / dt)))
    env = np.abs(np.convolve(c * np.exp(-1j * omega * t), np.ones(m) / m,
                             mode="valid"))
    gamma = -np.polyfit(t[m // 2:m // 2 + len(env)], np.log(env), 1)[0]
    return abs(omega), gamma, len(t)


def _validate_whistler(text, damp_min, damp_max):
    fit = _whistler_fit(text)
    if fit is None:
        return False, "Could not fit the seeded whistler from the plot frames"
    omega, gamma, n = fit
    x_lo, x_hi = _block_numbers(text, "#GEOMETRY", 2)
    kappa = 2.0 * math.pi * _named_number(text, "#WAVEIC", "waveMode") / abs(
        x_hi - x_lo)
    omega_hall = 0.5 * (kappa ** 2 + math.sqrt(kappa ** 4 + 4.0 * kappa ** 2))
    damp = gamma / omega
    logger.debug("    [UPWIND] whistler k d_i = %.3f: omega = %.4f (Hall %.4f, "
                 "%+.2f%%), gamma = %.4g, gamma/omega = %.4g over %d frames",
                 kappa, omega, omega_hall, 100.0 * (omega / omega_hall - 1.0),
                 gamma, damp, n)
    if abs(omega / omega_hall - 1.0) > WHISTLER_OMEGA_TOL:
        return False, ("whistler frequency %.4f is off the Hall branch %.4f by "
                       "%.1f%% (tolerance %.0f%%)"
                       % (omega, omega_hall, 100.0 * (omega / omega_hall - 1.0),
                          100.0 * WHISTLER_OMEGA_TOL))
    if not damp_min <= damp <= damp_max:
        return False, ("whistler damping per radian %.4g outside [%g, %g]"
                       % (damp, damp_min, damp_max))
    return True, "Passed"


def validate_plot(test_name):
    """divb: div(B) of the last plot equals that of the first.
    whistler: frequency and damping of the seeded whistler."""
    if _variant(test_name) == "whistler":
        return _validate_whistler(_param_text(), WHISTLER_DAMP_MIN,
                                  WHISTLER_DAMP_MAX)
    if _variant(test_name) == "whistler_sensor":
        return _validate_whistler(_param_text(), SENSOR_DAMP_MIN,
                                  SENSOR_DAMP_MAX)
    if _variant(test_name) != "divb":
        return True, "Passed (no plot-file check)"

    vidx0, rows0 = _run_dir.load_first_out()
    vidx1, rows1 = _run_dir.load_last_out()
    div0 = _run_dir.col(vidx0, rows0, "DIVB")
    div1 = _run_dir.col(vidx1, rows1, "DIVB")
    if div0 is None or div1 is None or len(div0) != len(div1):
        return False, "divB missing from the first/last plot"
    if len(_run_dir.plot_files()) < 2:
        return False, "need the first and the last plot"

    comps = [_run_dir.col(vidx1, rows1, c) for c in ("BX", "BY", "BZ")]
    if any(c is None for c in comps):
        return False, "B missing from the last plot"
    b_max = max(math.sqrt(bx * bx + by * by + bz * bz)
                for bx, by, bz in zip(*comps))
    by0 = _run_dir.col(vidx0, rows0, "BY")
    db_max = max(abs(a - b) for a, b in zip(comps[1], by0))

    text = _param_text()
    x_lo, x_hi = _block_numbers(text, "#GEOMETRY", 2)
    dx = abs(x_hi - x_lo) / _block_numbers(text, "#NCELL", 1)[0]

    d_div = max(abs(a - b) for a, b in zip(div1, div0))
    rel = d_div * dx / b_max if b_max > 0 else float("inf")
    logger.debug("    [UPWIND] divb: max|divB_end - divB_0| = %.3e "
                 "(%.2e of max|B|/dx), max|dBy| = %.3e", d_div, rel, db_max)
    if db_max <= 0.0:
        return False, "B did not evolve"
    if not rel <= DIVB_TOL:
        return False, ("div(B) changed by %.2e of max|B|/dx (tolerance %.0e)"
                       % (rel, DIVB_TOL))
    return True, "Passed"
