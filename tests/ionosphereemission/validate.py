#!/usr/bin/env python3
"""Validator for the ionospheric emission shell (tests/ionosphereemission).

H+ (1e16 /s) and O+ (2e15 /s) are emitted at r = 2 above an absorbing body
(R = 1.5) into vacuum, with E = B = 0. In code units one proton has the mass
m1 = e^2 / (1e7 lNormSI m_p) (FluidInterface: Mnorm = 1e3 mNormSI [g]), so

* the injected mass grows at (1e16 * 1 + 2e15 * 16) m1 per second (1e-6);
* the emission is outward, so the body absorbs nothing (< 1e-3);
* within 2 s no ion reaches the walls: the mass of each species in the box
  grows at its own rate times its mass (1e-3; the vacuum seed is negligible).

validate_plot: step 1 holds the particles of one injection only, so the
ratio of the ion mass on the nodes with x > 0 to those with x < 0 is the ratio
of the shell weights, sum(area * a) per hemisphere, less half of the cells
next to the plane x = 0 (deposited onto it). The weights are computed
here as in Pic::build_emission_weights (deterministic stratified sampling):
a = night + (noon - night) cos(SZA) for the default cosSza profile (5.5 for a
smooth sphere), and a = night + (noon - night) (Phi/Phi_subsolar)^0.5 with the
photoionization column Phi of the photochem variant, an optically thick layer
for which Phi ~ cos(SZA) (7.0 for a smooth sphere). Tolerance 3%.
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


def _ions(test_name):
    """Index of the first ion species (photochem has the electron slot)."""
    return 1 if test_name and test_name.endswith("photochem") else 0

E_SI = 1.602176634e-19
MP_SI = 1.67262192369e-27
L_NORM = 1.0e5
M1 = E_SI**2 / (1.0e7 * L_NORM * MP_SI)
RATES = [1.0e16, 2.0e15]
MASSES = [1.0, 16.0]


def validate_log(pic_diags=None, test_name=None):
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

    i0 = _ions(test_name)
    for i, (r, m) in enumerate(zip(RATES, MASSES)):
        grow = (b[f"Mpart{i + i0}"] - a[f"Mpart{i + i0}"]) / dt / (r * m * M1)
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
# photochem layer (SI)
N0, H0, SIGMA, NU = 1.0e12, 1.0e3, 1.0e-13, 1.0e-3
RP, R1 = 1.5e5, 2.0e5
MINPROD, TAUFLOOR, TAUCUT, MUFLOOR = 1.0e-10, 6.0e-3, 11.5, 5.0e-4


def _attenuation(x, r):
    """ExoSource::photo_attenuation (non-Chapman) for the solar dir +x."""
    proj = x
    if proj <= 0:
        p = math.sqrt(max(r * r - proj * proj, 0.0))
        if p <= RP:
            return MINPROD
        n = N0 * math.exp(-(p - RP) / H0)
        tau = max(SIGMA * n * H0 * math.sqrt(2 * math.pi * p / H0), TAUFLOOR)
    else:
        n = N0 * math.exp(-(r - RP) / H0)
        tau = max(SIGMA * n * H0, TAUFLOOR) / max(proj / r, MUFLOOR)
    if tau >= TAUCUT:
        return MINPROD
    return max(math.exp(-tau), MINPROD)


def _column(d):
    """ExoSource::column_photoproduction_si along the unit direction d."""
    col, r = 0.0, RP
    dr0 = min((R1 - RP) / 200, H0 / 8)
    while r < R1:
        dr = min(dr0, R1 - r)
        rm = r + 0.5 * dr
        col += NU * N0 * math.exp(-(rm - RP) / H0) * \
            _attenuation(rm * d[0], rm) * dr
        r += dr
    return col


def _hemisphere_weights(photochem):
    dx = 8.0 / NCELL
    c = LO + (np.arange(NCELL) + 0.5) * dx
    nSub, h = 24, 0.2 * dx
    sub = ((np.arange(nSub) + 0.5) / nSub - 0.5) * dx
    sx, sy, sz = np.meshgrid(sub, sub, sub, indexing="ij")
    dV = dx**3 / nSub**3
    phiSS = _column((1.0, 0.0, 0.0)) if photochem else 1.0
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
                if photochem:
                    phi = _column((xc / rc, yc / rc, zc / rc)) / phiSS
                    a = NIGHT + (NOON - NIGHT) * math.sqrt(max(phi, 0.0))
                else:
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
    photochem = bool(test_name and test_name.endswith("photochem"))
    i0 = _ions(test_name)
    ratios = []
    for iS in (i0, i0 + 1):
        rho = d[:, names.index(f"rhoS{iS}")]
        eps = 1e-9 * (x.max() - x.min())
        ratios.append(rho[x > eps].sum() / rho[x < -eps].sum())
    wDay, wNight = _hemisphere_weights(photochem)
    expected = wDay / wNight
    reasons = []
    for iS, ratio in zip((i0, i0 + 1), ratios):
        logger.debug("    species %d: day/night emission %.4f, expected %.4f",
                     iS, ratio, expected)
        if abs(ratio / expected - 1) > 0.03:
            reasons.append(f"species {iS}: day/night = {ratio:.3f}, "
                           f"expected {expected:.3f} (tol 3%)")
    if reasons:
        return False, "; ".join(reasons)
    return True, f"Passed (day/night {ratios[0]:.3f} vs {expected:.3f})"
