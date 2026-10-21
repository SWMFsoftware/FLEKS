#!/usr/bin/env python3
"""Validator for the ionospheric emission shell (tests/ionosphereemission).

H+ (1e16 /s) and O+ (2e15 /s) are emitted at r = 2 above an absorbing body
(R = 1.5) into vacuum, with E = B = 0. In code units one proton has the mass
m1 = e^2 / (1e7 lNormSI m_p) (FluidInterface: Mnorm = 1e3 mNormSI [g]), so

* the injected mass grows at (1e16 * 1 + 2e15 * 16) m1 per second (1e-6);
* the emission is outward, so the body absorbs nothing (< 1e-3);
* within 1 s no ion reaches the walls: the mass of each species in the box
  grows at its own rate times its mass (1e-3; the vacuum seed is negligible).

validate_plot: step 1 holds the particles of one injection only, so the
ratio of the ion mass on the nodes with x > 0 to those with x < 0 is the ratio
of the shell weights, sum(area * a) per hemisphere, less half of the cells
next to the plane x = 0 (deposited onto it). The weights are computed
here as in Pic::build_emission_weights (deterministic stratified sampling):
a = night + (noon - night) cos(SZA) on the dayside and night on the nightside
(5.5 for a smooth sphere). Tolerance 3%.
"""
import glob
import logging
import math
import os

import numpy as np

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"


def set_run_dir(run_dir):
    global RUN_DIR
    RUN_DIR = run_dir


E_SI = 1.602176634e-19
MP_SI = 1.67262192369e-27
L_NORM = 1.0e5
M1 = E_SI**2 / (1.0e7 * L_NORM * MP_SI)
RATES = [1.0e16, 2.0e15]
MASSES = [1.0, 16.0]


def _non_finite(diags):
    """Reason string if a log value is NaN/inf (comparisons with NaN are false,
    so a tolerance check alone would pass a broken run), else None."""
    keys = sorted({k for e in (diags or []) for k, v in e.items()
                   if isinstance(v, float) and not math.isfinite(v)})
    if keys:
        return f"non-finite log values in {', '.join(keys[:6])}"
    return None


def validate_log(pic_diags=None, test_name=None):
    bad = _non_finite(pic_diags)
    if bad:
        return False, bad
    d = pic_diags
    if not d or len(d) < 10 or "mBodyInject" not in d[0]:
        return False, "missing log or mBodyInject column"

    reasons = []

    def check(cond, msg):
        logger.debug("    %s: %s", "OK  " if cond else "FAIL", msg)
        if not cond:
            reasons.append(msg)

    a, b = d[1], d[-1]
    dt = b["time"] - a["time"]
    expected = sum(r * m for r, m in zip(RATES, MASSES)) * M1
    rate = (b["mBodyInject"] - a["mBodyInject"]) / dt
    check(abs(rate / expected - 1) < 1e-6,
          f"injected mass rate = {rate / expected:.8f} x expected (tol 1e-6)")

    absorbed = b["mBodyAbsorb"] / max(b["mBodyInject"], 1e-300)
    check(absorbed < 1e-3,
          f"absorbed / injected mass = {absorbed:.2e} (outward emission, "
          f"tol 1e-3)")

    for i, (r, m) in enumerate(zip(RATES, MASSES)):
        grow = (b[f"Mpart{i}"] - a[f"Mpart{i}"]) / dt / (r * m * M1)
        check(abs(grow - 1) < 1e-3,
              f"species {i}: box mass grows at {grow:.6f} x its emission rate "
              f"(tol 1e-3)")

    if reasons:
        return False, "; ".join(reasons)
    return True, "Passed"


# ---------------------------------------------------------------------------
# Shell weights, as in Pic::build_emission_weights
# ---------------------------------------------------------------------------
NCELL, LO, RSHELL, RBODY = 24, -4.0, 2.0, 1.5   # code units (lNorm = 1e5 m)
NOON, NIGHT = 1.0, 0.1


def _hemisphere_weights():
    dx = 8.0 / NCELL
    c = LO + (np.arange(NCELL) + 0.5) * dx
    nSub, h = 24, 0.2 * dx
    sub = ((np.arange(nSub) + 0.5) / nSub - 0.5) * dx
    sx, sy, sz = np.meshgrid(sub, sub, sub, indexing="ij")
    dV = dx**3 / nSub**3
    # Node sums off the plane x = 0: the linear deposit of particles placed
    # uniformly in a cell gives half the mass of the cells next to the plane
    # to the nodes on it.
    day = night = 0.0
    for xc in c:
        for yc in c:
            for zc in c:
                rc = math.sqrt(xc * xc + yc * yc + zc * zc)
                if rc < RBODY or abs(rc - RSHELL) > 0.9 * dx + h:
                    continue
                r = np.sqrt((xc + sx)**2 + (yc + sy)**2 + (zc + sz)**2)
                nHit = np.count_nonzero(np.abs(r - RSHELL) < 0.5 * h)
                if nHit == 0:
                    continue
                cosSza = xc / rc
                a = NIGHT + (NOON - NIGHT) * cosSza if cosSza > 0 else NIGHT
                w = nHit * dV / h * a
                if abs(xc) < dx:
                    w *= 0.5
                if xc > 0:
                    day += w
                else:
                    night += w
    return day, night


def _read_out(path):
    with open(path, "r", encoding="latin-1") as f:
        lines = f.readlines()
    names = lines[4].split()
    data = np.array([[float(v) for v in l.split()] for l in lines[5:]
                     if l.split()])
    return names, data


def validate_plot(test_name=None):
    plots = sorted(glob.glob(os.path.join(RUN_DIR, "PC", "plots",
                                          "3d*_n00000001*.out")))
    if not plots:
        return False, "no 3d plot of step 1"
    names, d = _read_out(plots[0])
    x = d[:, 0]
    ratios = []
    for iS in (0, 1):
        rho = d[:, names.index(f"rhoS{iS}")]
        eps = 1e-9 * (x.max() - x.min())
        ratios.append(rho[x > eps].sum() / rho[x < -eps].sum())
    wDay, wNight = _hemisphere_weights()
    expected = wDay / wNight
    reasons = []
    for iS, ratio in zip((0, 1), ratios):
        logger.debug("    species %d: day/night emission %.4f, expected %.4f",
                     iS, ratio, expected)
        if not math.isfinite(ratio) or abs(ratio / expected - 1) > 0.03:
            reasons.append(f"species {iS}: day/night = {ratio:.3f}, "
                           f"expected {expected:.3f} (tol 3%)")
    if reasons:
        return False, "; ".join(reasons)
    return True, f"Passed (day/night {ratios[0]:.3f} vs {expected:.3f})"
