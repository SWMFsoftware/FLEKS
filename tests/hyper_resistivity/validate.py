#!/usr/bin/env python3
"""Validator for the hybrid-PIC hyper-resistivity tests.

The runner executes one variant per PARAM file and passes the variant name
(PARAM.in.<suffix> -> "<test>_<suffix>") as ``test_name``:

  hyper_resistivity
      Full-solver integration run (grid-mode hyper-resistivity on top of
      convection + Hall + eta*J + grad Pe).  Generic completion check only;
      the physics is exercised but the late-time state is noise driven, so
      there is nothing quantitative to assert.

  hyper_resistivity_damping
      Frozen plasma (no macroparticles -> rho = 0 -> every other Ohm term is
      inert) + "si" mode hyper-resistivity + one seeded transverse mode.  The
      measured exponential decay rate of that mode is compared against the
      analytic discrete symbol of the term,

          gamma = (eta_h/4*pi) * 4 sin^2(theta) sin^2(theta/2) / dx^4,
          theta = k*dx,

      which follows from the two operators the term is built from: the
      collocated 2*dx curl (symbol i*sin(theta)/dx) and the 3-point Laplacian
      (symbol -4 sin^2(theta/2)/dx^2).  This is the check that catches a
      silently disabled or mis-converted eta_h, a wrong sign, and any change
      to the curl / Laplacian stencils or to the RK trial-state ghosts.

  hyper_resistivity_damping_nyquist
      Same, seeded at the Nyquist mode.  The 2*dx curl annihilates it
      (sin(pi) = 0), so the term cannot damp it by construction and the
      magnetic energy must stay constant.
"""
import logging
import math
import os

logger = logging.getLogger(__name__)

# Relative tolerance of the measured damping rate against the analytic value.
TOL = 0.02
# The Nyquist mode must be untouched down to round-off.
FROZEN_TOL = 1e-10
# Minimum number of energy-log samples for a fit.
MIN_SAMPLES = 20


def _param_path(variant):
    """Path of the PARAM.in file driving the given variant."""
    base = os.path.join("tests", "hyper_resistivity", "PARAM.in")
    return base if not variant else base + "." + variant


def _block_numbers(text, command, n):
    """First n numeric values of the PARAM block `command`."""
    vals, capture = [], False
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        tok = s.split()[0]
        if tok.startswith("#"):
            capture = (tok == command)
            continue
        if capture:
            try:
                vals.append(float(tok))
            except ValueError:
                pass
    return vals[:n]


def _named_number(text, command, name):
    """Value of the `value   name` pair inside the PARAM block `command`."""
    capture = False
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        tok = s.split()[0]
        if tok.startswith("#"):
            capture = (tok == command)
            continue
        if capture:
            cols = s.split()
            if len(cols) >= 2 and cols[1] == name:
                try:
                    return float(cols[0])
                except ValueError:
                    return None
    return None


def _analytic_gamma(variant):
    """Return (gamma, gamma_as_measured_by_RK4) for the seeded mode.

    gamma is the exact discrete rate of the term; the field integrator is
    classical RK4, whose amplification on a real eigenvalue -gamma is
    R(-z) = 1 - z + z^2/2 - z^3/6 + z^4/24 with z = gamma*dt, so a fit of
    exp(-2*gamma_fit*t) measures -ln(R)/dt rather than gamma itself.
    """
    text = open(_param_path(variant)).read()
    l_norm, u_norm = _block_numbers(text, "#NORMALIZATION", 2)
    x_lo, x_hi = _block_numbers(text, "#GEOMETRY", 2)
    n_cell = _block_numbers(text, "#NCELL", 1)[0]
    mode = _named_number(text, "#WAVEIC", "waveMode")
    eta_si = _named_number(text, "#HYPERRESISTIVITY", "etaHyperSI")
    dt = _named_number(text, "#TIMESTEPPING", "dt")
    if None in (mode, eta_si, dt) or n_cell <= 0:
        return None, None

    dx = abs(x_hi - x_lo) / n_cell
    # eta_h_code = 4*pi*eta_h_SI*Si2NoV*Si2NoL^3 with
    # Si2NoV = 100/Unorm = 1/uNormSI and Si2NoL = 100/Lnorm = 1/lNormSI.
    eta_code = 4.0 * math.pi * eta_si / (u_norm * l_norm ** 3)

    theta = 2.0 * math.pi * mode / n_cell
    shape = 4.0 * math.sin(theta) ** 2 * math.sin(theta / 2.0) ** 2
    gamma = (eta_code / (4.0 * math.pi)) * shape / dx ** 4

    z = gamma * dt
    amp = 1.0 - z + z * z / 2.0 - z ** 3 / 6.0 + z ** 4 / 24.0
    return gamma, (-math.log(amp) / dt if 0.0 < amp < 1.0 else gamma)


def _fit_gamma(t, eb):
    """Fit Eb(t) = C + A*exp(-2*gamma*t) via dEb/dt = -2*gamma*(Eb - C).

    Regressing the derivative against Eb avoids having to know the constant
    (guide-field) part of the magnetic energy.
    """
    n = len(t) - 1
    if n < 2:
        return None
    x = [0.5 * (eb[i] + eb[i + 1]) for i in range(n)]
    y = [(eb[i + 1] - eb[i]) / (t[i + 1] - t[i]) for i in range(n)]
    sx, sy = sum(x), sum(y)
    sxx = sum(v * v for v in x)
    sxy = sum(x[i] * y[i] for i in range(n))
    den = n * sxx - sx * sx
    if den == 0.0:
        return None
    return -(n * sxy - sx * sy) / den / 2.0


def _check_damping(pic_diags, variant):
    if not pic_diags or len(pic_diags) < MIN_SAMPLES:
        return True, "Too few energy-log samples (skipped)"
    t = [row["time"] for row in pic_diags]
    eb = [row["Eb"] for row in pic_diags]
    if any(not math.isfinite(v) for v in eb):
        return False, "Eb is not finite (hyper-resistive term unstable?)"
    if eb[-1] > eb[0]:
        return False, "Eb grew: the hyper-resistive term is not dissipative"

    g_exact, g_rk4 = _analytic_gamma(variant)
    if g_exact is None:
        return False, "Could not derive the analytic decay rate from PARAM.in"
    if g_exact <= 0.0:
        return False, ("PARAM.in gives a zero analytic decay rate (etaHyperSI "
                       "= 0, or the seeded mode is in the curl null space): "
                       "this variant must damp the seeded mode")
    gamma = _fit_gamma(t, eb)
    if gamma is None or gamma <= 0:
        return False, "The seeded mode did not decay (hyper-resistivity inert?)"

    ratio = gamma / g_rk4
    logger.debug("    [HYPER] gamma_fit = %.6g, gamma_analytic = %.6g "
                 "(RK4: %.6g), ratio = %.4f", gamma, g_exact, g_rk4, ratio)
    if abs(ratio - 1.0) > TOL:
        return False, ("measured decay rate %.6g differs from the analytic "
                       "%.6g by %.1f%% (tolerance %.1f%%)"
                       % (gamma, g_rk4, 100.0 * (ratio - 1.0), 100.0 * TOL))
    return True, "Passed"


def _check_frozen(pic_diags, variant):
    if not pic_diags or len(pic_diags) < 2:
        return True, "Too few energy-log samples (skipped)"
    eb0 = pic_diags[0]["Eb"]
    eb1 = pic_diags[-1]["Eb"]
    if eb0 <= 0 or not (math.isfinite(eb0) and math.isfinite(eb1)):
        return False, "Eb is not usable (hyper-resistive term unstable?)"
    drift = abs(eb1 - eb0) / eb0
    logger.debug("    [HYPER] Nyquist mode: Eb %.6e -> %.6e (rel. change %.2e)",
                 eb0, eb1, drift)
    if drift > FROZEN_TOL:
        return False, ("Nyquist-mode energy changed by %.2e (tolerance %.1e): "
                       "the 2*dx curl should annihilate this mode"
                       % (drift, FROZEN_TOL))
    return True, "Passed"


def validate_log(pic_diags=None, test_name=None):
    """Energy-log checks for the hyper-resistivity variants."""
    if test_name == "hyper_resistivity_damping":
        return _check_damping(pic_diags, "damping")
    if test_name == "hyper_resistivity_damping_nyquist":
        return _check_frozen(pic_diags, "damping_nyquist")
    return True, "Passed (generic completion check)"
