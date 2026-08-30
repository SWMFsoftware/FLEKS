#!/usr/bin/env python3
"""Validator for the hybrid-PIC hyper-resistivity test (tests/hyper_resistivity).

The test is a frozen-plasma setup: no macroparticles are loaded, so rho = 0
everywhere, every other term of the generalized Ohm's law (convection, Hall,
resistive eta*J, grad Pe) is short-circuited by the rho > 0 guard, and the
hyper-resistive term is the only thing driving B.  A single circularly
polarized transverse mode is seeded on the uniform guide field, so the
magnetic energy must decay as

    Eb(t) - Eb_guide = A * exp(-2*gamma*t)

with the discrete rate of the term (see README.md)

    gamma = (eta_h/4*pi) * 4 sin^2(theta) sin^2(theta/2) / dx^4,  theta = k*dx

which follows from the two operators the term is built from: the collocated
2*dx curl (symbol i*sin(theta)/dx) and the 3-point Laplacian (symbol
-4 sin^2(theta/2)/dx^2).

Three assertions, all derived from the PARAM.in values:

1. dissipative: Eb decreases and stays finite;
2. the decay is a clean single exponential -- a corrupted stencil or stale
   ghost cells show up as a non-exponential decay long before they show up as
   a wrong rate;
3. the fitted rate matches the analytic symbol within TOL (the RK4
   amplification -ln(R(-gamma*dt))/dt is folded into the expectation).

This is what catches a silently disabled or mis-converted eta_h, a wrong sign,
and any regression in the stencils or in the RK trial-state ghosts.
"""
import logging
import math
import os

logger = logging.getLogger(__name__)

PARAM_PATH = os.path.join("tests", "hyper_resistivity", "PARAM.in")

# Relative tolerance of the measured damping rate against the analytic value.
# A correct run agrees to ~1e-10.
TOL = 0.005
# Bound on the non-exponential part of the decay, relative to its amplitude.
# A correct run sits at ~1e-12; stale RK-stage ghosts gave 2e-2.
MAX_RESIDUAL = 1e-6
# Minimum number of energy-log samples for a fit.
MIN_SAMPLES = 20


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


def _analytic_gamma():
    """Return (gamma, gamma_as_measured) for the seeded mode.

    gamma is the exact discrete rate of the term.  The field integrator is
    classical RK4, whose amplification on a real eigenvalue -gamma is
    R(-z) = 1 - z + z^2/2 - z^3/6 + z^4/24 with z = gamma*dt, so a fit of
    exp(-2*gamma_fit*t) measures -ln(R)/dt rather than gamma itself.

    eta_h_code = 4*pi * etaHyperSI * Si2NoV * Si2NoL^3, with
    Si2NoV = 100/Unorm = 1/uNormSI and Si2NoL = 100/Lnorm = 1/lNormSI.
    """
    text = open(PARAM_PATH).read()
    l_norm, u_norm = _block_numbers(text, "#NORMALIZATION", 2)
    x_lo, x_hi = _block_numbers(text, "#GEOMETRY", 2)
    n_cell = _block_numbers(text, "#NCELL", 1)[0]
    mode = _named_number(text, "#WAVEIC", "waveMode")
    eta_si = _named_number(text, "#HYPERRESISTIVITY", "etaHyperSI")
    dt = _named_number(text, "#TIMESTEPPING", "dt")
    if None in (mode, eta_si, dt) or n_cell <= 0:
        return None, None

    dx = abs(x_hi - x_lo) / n_cell
    eta_code = 4.0 * math.pi * eta_si / (u_norm * l_norm ** 3)

    theta = 2.0 * math.pi * mode / n_cell
    shape = 4.0 * math.sin(theta) ** 2 * math.sin(theta / 2.0) ** 2
    gamma = (eta_code / (4.0 * math.pi)) * shape / dx ** 4

    z = gamma * dt
    amp = 1.0 - z + z * z / 2.0 - z ** 3 / 6.0 + z ** 4 / 24.0
    return gamma, (-math.log(amp) / dt if 0.0 < amp < 1.0 else gamma)


def _median(vals):
    vals = sorted(vals)
    return vals[len(vals) // 2]


def _fit_decay(t, eb):
    """Fit Eb(t) = C + A*exp(-2*gamma*t); return (gamma, A, residual).

    The log is sampled every step and a single decaying mode is exactly
    geometric, so the ratio rho = exp(-2*gamma*dt) follows from three
    consecutive samples and the constant (guide-field) part C then follows
    from two:

        rho   = (eb[i+2] - eb[i+1]) / (eb[i+1] - eb[i])
        C     = (eb[i+1] - rho*eb[i]) / (1 - rho)

    Using the median over the samples that are well above the round-off floor
    makes both robust.  This matters: estimating C by regression instead
    biases gamma at the 1e-4 level and leaves a ~1e-5 residual, because the
    run only covers about three and a half decades of decay.  Falls back to
    that regression if the sample spacing is not uniform.
    """
    n = len(t)
    if n < 10:
        return None, None, None

    dts = [t[i + 1] - t[i] for i in range(n - 1)]
    dt = dts[0]
    uniform = dt > 0.0 and all(abs(d - dt) <= 1e-6 * dt for d in dts)
    ref = abs(eb[1] - eb[0])
    # Samples whose step-to-step change is still a sizeable fraction of the
    # first one, i.e. safely above the round-off floor.
    strong = [i for i in range(n - 2)
              if ref > 0.0 and abs(eb[i + 1] - eb[i]) > 1e-3 * ref]

    if uniform and len(strong) >= 5:
        rho = _median([(eb[i + 2] - eb[i + 1]) / (eb[i + 1] - eb[i])
                       for i in strong])
        if 0.0 < rho < 1.0:
            gamma = -math.log(rho) / (2.0 * dt)
            c = _median([(eb[i + 1] - rho * eb[i]) / (1.0 - rho)
                         for i in strong])
            amp = eb[0] - c
            if amp != 0.0:
                resid = max(
                    abs(eb[i] - c - amp * math.exp(-2.0 * gamma * (t[i] - t[0])))
                    for i in range(n)
                    if (eb[i] - c) > 1e-4 * abs(amp)) / abs(amp)
                return gamma, amp, resid

    # Fallback: dEb/dt = -2*gamma*(Eb - C) regression (non-uniform sampling).
    x = [0.5 * (eb[i] + eb[i + 1]) for i in range(n - 1)]
    y = [(eb[i + 1] - eb[i]) / dts[i] for i in range(n - 1)]
    m = len(x)
    sx, sy = sum(x), sum(y)
    sxx = sum(v * v for v in x)
    sxy = sum(x[i] * y[i] for i in range(m))
    den = m * sxx - sx * sx
    if den == 0.0:
        return None, None, None
    slope = (m * sxy - sx * sy) / den
    intercept = (sy - slope * sx) / m
    gamma = -slope / 2.0
    if gamma <= 0.0:
        return gamma, None, None
    c = intercept / (2.0 * gamma)
    amp = eb[0] - c
    if amp == 0.0:
        return gamma, amp, None
    resid = max(abs(eb[i] - c - amp * math.exp(-2.0 * gamma * (t[i] - t[0])))
                for i in range(n) if (eb[i] - c) > 1e-4 * abs(amp)) / abs(amp)
    return gamma, amp, resid


def validate_log(pic_diags=None, test_name=None):
    """Energy-log check: rate and purity of the hyper-resistive damping."""
    if not pic_diags or len(pic_diags) < MIN_SAMPLES:
        return True, "Too few energy-log samples (skipped)"

    t = [row["time"] for row in pic_diags]
    eb = [row["Eb"] for row in pic_diags]
    if any(not math.isfinite(v) for v in eb):
        return False, "Eb is not finite (hyper-resistive term unstable?)"
    if eb[-1] > eb[0]:
        return False, "Eb grew: the hyper-resistive term is not dissipative"

    g_exact, g_expect = _analytic_gamma()
    if g_exact is None:
        return False, "Could not derive the analytic decay rate from PARAM.in"
    if g_exact <= 0.0:
        return False, ("PARAM.in gives a zero analytic decay rate (etaHyperSI"
                       " = 0, or the seeded mode sits in the curl null space)")

    gamma, amp, resid = _fit_decay(t, eb)
    if gamma is None or gamma <= 0.0:
        return False, "The seeded mode did not decay (hyper-resistivity inert?)"

    ratio = gamma / g_expect
    logger.debug("    [HYPER] gamma_fit = %.6g, gamma_analytic = %.6g "
                 "(RK4: %.6g), ratio = %.4f, residual = %.2e",
                 gamma, g_exact, g_expect, ratio, resid)

    if abs(ratio - 1.0) > TOL:
        return False, ("measured decay rate %.6g differs from the analytic "
                       "%.6g by %.1f%% (tolerance %.1f%%)"
                       % (gamma, g_expect, 100.0 * (ratio - 1.0), 100.0 * TOL))
    if resid is not None and resid > MAX_RESIDUAL:
        return False, ("decay is not a single exponential: residual %.2e "
                       "(tolerance %.1e) -- check the stencils and the RK "
                       "trial-state ghosts" % (resid, MAX_RESIDUAL))
    return True, "Passed"
