#!/usr/bin/env python3
"""Validator for the Langmuir test (tests/langmuir).

A 1D uniform thermal electron-ion plasma (m_i/m_e = 100, periodic, no drive).
Its thermal noise excites Langmuir oscillations at every wavenumber, so the
test checks the particle-field coupling of the full-PIC solvers in two ways:

1. Energy (log): nothing drives the plasma, so the total energy stays
   constant. Bound on max |Etot/Etot0 - 1| over all logged steps:
     explicit, explicit_sub, explicit_mag   1%
     explicit_pc                            3% (dissipative: -1.2% in 16 s)
     implicit                               no heating (< +1%), cooling < 20%
                                            (theta = 0.51 damps: -11% in 16 s)
2. Langmuir branch (plot): the omega-k spectrum of E_x(x, t) must peak at the
   Bohm-Gross frequency omega^2 = omega_pe^2 (1 + 3 k^2 lambda_De^2) for the
   12 lowest modes (k lambda_De <= 0.11, so the thermal term reaches 1.7%),
   with each peak at least 100 times the median of its spectrum. The
   frequency of a mode is the power-weighted centroid of its spectrum within
   8% of omega_pe around the peak, over the 16 s record (one noise
   realization scatters single modes by about 1%; an 8 s record or the bare
   peak bin scatter by up to 4%). Leapfrog variants: each mode within 4%, the
   mean within 1%. The predictor-corrector
   and implicit schemes run 1-2% low, more at high k (numerical dispersion):
   each mode within 6%, the mean within 3%. Skipped for explicit_mag, where
   the oblique field splits the electrostatic branch.
3. Gauss's law (plot, explicit variants): the explicit solver deposits a
   charge-conserving (Esirkepov) current and runs no div(E) correction, so
   the residual r = div_node_to_center(E) - 4 pi rho_c, with rho_c the node
   charge sum_s (q_s/m_s) rhoS_s averaged to the cell centers, must keep its
   initial value: max |r(t) - r(t_0)| < 1e-5 max |div E| over the record.
   The ascii output (8 digits) sets a floor near 1e-6; particle resampling,
   which would break the conservation, is off in these decks.

A current with the wrong sign or magnitude moves the branch, and one that is
centered at the wrong time makes it grow; FLEKS issue 170 (explicit solver
using the mass flux, centered at t^n) heated this plasma by a factor of ~560.
"""
import logging
import math

import numpy as np

from tests._shared import run_dir as _run_dir

logger = logging.getLogger(__name__)

set_run_dir = _run_dir.set_run_dir

# Deck parameters in code units (c = 1, n = 1, q = 1, m_e = 0.01,
# T_e = 1.25e-3): omega_pe = sqrt(4 pi n q^2 / m_e), lambda_De = v_te/omega_pe.
OMEGA_PE = math.sqrt(4.0 * math.pi / 0.01)
LAMBDA_DE = math.sqrt(1.25e-3 / 0.01) / OMEGA_PE
BOX_LENGTH = 2.0 * 3.544908
MODES = range(1, 13)  # k lambda_De <= 0.11
MIN_PEAK_RATIO = 100.0
GAUSS_TOL = 1e-5  # max |r(t) - r(t_0)| / max |div E|
CENTROID_WIDTH = 0.08  # half-width of the centroid window / omega_pe
# (per-mode, mean) tolerance on omega / omega_Bohm-Gross - 1.
FREQ_TOL = {"explicit": (0.04, 0.01), "explicit_sub": (0.04, 0.01),
            "explicit_pc": (0.06, 0.03), "implicit": (0.06, 0.03)}


def _variant(test_name):
    name = (test_name or "").lower()
    for v in ("explicit_sub", "explicit_mag", "explicit_pc", "implicit"):
        if name.endswith(v):
            return v
    return "explicit"


def validate_log(pic_diags=None, test_name=None):
    variant = _variant(test_name)
    if not pic_diags or len(pic_diags) < 10:
        n = len(pic_diags) if pic_diags else 0
        return False, f"Expected >= 10 energy log entries, got {n}"

    e0 = pic_diags[0].get("Etot", 0.0)
    if not (math.isfinite(e0) and e0 > 0.0):
        return False, f"Initial Etot = {e0} is not finite and positive"

    ratios = []
    for d in pic_diags:
        e = d.get("Etot", float("nan"))
        if not math.isfinite(e):
            return False, f"Non-finite Etot at cycle {d.get('cycle')}"
        ratios.append(e / e0)
    rmax, rmin = max(ratios), min(ratios)
    logger.debug("    Etot/Etot0: final %.5f, max %.5f, min %.5f", ratios[-1],
                 rmax, rmin)

    if variant == "implicit":
        if rmax > 1.01 or rmin < 0.80:
            return False, (f"Etot/Etot0 left [0.80, 1.01]: max {rmax:.4f}, "
                           f"min {rmin:.4f}")
        return True, f"Passed (Etot/Etot0 in [{rmin:.4f}, {rmax:.4f}])"

    tol = 0.03 if variant == "explicit_pc" else 0.01
    worst = max(rmax - 1.0, 1.0 - rmin)
    if worst > tol:
        return False, (f"Energy not conserved: max |Etot/Etot0 - 1| = "
                       f"{worst:.4f} > {tol} (final ratio {ratios[-1]:.4f})")
    return True, f"Passed (max |Etot/Etot0 - 1| = {worst:.4f} < {tol})"


def _read_frames():
    """Return the ascii z=0 'var' plot frames sorted in time as a list of
    (t, dx, {var: values on the x-sorted nodes}, {param: value})."""
    frames = []
    for f in _run_dir.plot_files():
        with open(f, "r", encoding="latin-1") as fh:
            lines = fh.readlines()
        if len(lines) < 6:
            continue
        names = lines[4].split()
        params = [float(v) for v in lines[3].split()]
        rows = np.array([[float(c) for c in l.split()] for l in lines[5:]
                         if l.strip()])
        nVar = rows.shape[1]
        if names[0] != "x" or "Ex" not in names[:nVar]:
            continue
        rows = rows[np.argsort(rows[:, 0])]
        data = {n: rows[:, i] for i, n in enumerate(names[:nVar])}
        frames.append((float(lines[1].split()[1]), rows[1, 0] - rows[0, 0],
                       data, dict(zip(names[nVar:], params))))
    frames.sort(key=lambda fr: fr[0])
    return frames


def _check_gauss(frames):
    """Drift of the Gauss's law residual (see check 3 in the module doc)."""
    if "rhoS1" not in frames[0][2]:
        return False, "Expected rhoS0 and rhoS1 in the plots (Gauss check)"
    r0, worst, scale = None, 0.0, 0.0
    for t, dx, data, par in frames:
        ex = data["Ex"]
        rho = sum(data[f"rhoS{s}"] * par[f"qS{s}"] / par[f"mS{s}"]
                  for s in (0, 1))
        # 1 cell along y (periodic): the y terms of div_node_to_center vanish
        # and the node charge averages over the two x nodes of each cell.
        divE = (np.roll(ex, -1) - ex) / dx
        r = divE - 4.0 * math.pi * 0.5 * (rho + np.roll(rho, -1))
        if r0 is None:
            r0 = r
        worst = max(worst, float(np.max(np.abs(r - r0))))
        scale = max(scale, float(np.max(np.abs(divE))))
    rel = worst / scale
    logger.debug("    Gauss residual drift: max |r - r0| = %.3e (%.2e of "
                 "max |div E|)", worst, rel)
    if rel > GAUSS_TOL:
        return False, (f"Gauss's law residual drifted by {rel:.2e} of max "
                       f"|div E| (> {GAUSS_TOL:.0e})")
    return True, f"Gauss residual drift {rel:.1e} of max |div E|"


def _langmuir_frequencies(t, ex):
    """Centroid |omega| of each mode's spectrum within CENTROID_WIDTH *
    omega_pe of its peak in the band [0.5, 2] omega_pe, the Bohm-Gross
    frequency, and the peak-to-median ratio."""
    amp = np.fft.rfft(ex, axis=1)  # a_m(t)
    window = np.hanning(len(t))[:, None]
    power = np.abs(np.fft.fft(amp * window, axis=0)) ** 2
    omega = 2.0 * np.pi * np.fft.fftfreq(len(t), t[1] - t[0])
    band = (np.abs(omega) > 0.5 * OMEGA_PE) & (np.abs(omega) < 2.0 * OMEGA_PE)
    results = []
    absOmega = np.abs(omega)
    for m in MODES:
        p = power[:, m]
        i = int(np.argmax(np.where(band, p, 0.0)))
        near = band & (np.abs(absOmega - absOmega[i]) < CENTROID_WIDTH * OMEGA_PE)
        measured = float(np.sum(absOmega[near] * p[near]) / np.sum(p[near]))
        k = 2.0 * math.pi * m / BOX_LENGTH
        theory = OMEGA_PE * math.sqrt(1.0 + 3.0 * (k * LAMBDA_DE) ** 2)
        results.append((m, measured, theory, p[i] / np.median(p[band])))
    return results


def validate_plot(test_name):
    variant = _variant(test_name)
    frames = _read_frames()
    if len(frames) < 50:
        return False, "Expected >= 50 E_x frames (z=0 var ascii plots)"
    # Periodic output may repeat the x = xMax node; the deck has 256 cells.
    if len(frames[0][2]["x"]) == 257:
        frames = [(t, dx, {k: v[:-1] for k, v in d.items()}, p)
                  for t, dx, d, p in frames]

    gauss = ""
    if variant != "implicit":
        ok, gauss = _check_gauss(frames)
        if not ok:
            return False, gauss
        gauss = "; " + gauss

    if variant == "explicit_mag":
        return True, ("No Langmuir-branch check in an oblique field" + gauss)
    mode_tol, mean_tol = FREQ_TOL[variant]

    t = np.array([fr[0] for fr in frames])
    ex = np.array([fr[2]["Ex"] for fr in frames])

    devs = []
    for m, measured, theory, peak in _langmuir_frequencies(t, ex):
        dev = measured / theory - 1.0
        devs.append(dev)
        logger.debug("    mode %2d: omega %.3f, Bohm-Gross %.3f (%+.2f%%), "
                     "peak/median %.0f", m, measured, theory, 100 * dev, peak)
        if peak < MIN_PEAK_RATIO:
            return False, (f"No clear Langmuir peak for mode {m} "
                           f"(peak/median {peak:.1f} < {MIN_PEAK_RATIO:.0f})")
        if abs(dev) > mode_tol:
            return False, (f"Mode {m}: omega = {measured:.3f}, Bohm-Gross "
                           f"{theory:.3f} ({100 * dev:+.1f}% > "
                           f"{100 * mode_tol:.0f}%)")
    mean_dev = sum(devs) / len(devs)
    if abs(mean_dev) > mean_tol:
        return False, (f"Mean Langmuir frequency off by {100 * mean_dev:+.2f}% "
                       f"(> {100 * mean_tol:.0f}%)")
    return True, (f"Passed (Langmuir branch: mean {100 * mean_dev:+.2f}%, "
                  f"max {100 * max(abs(d) for d in devs):.2f}% from Bohm-Gross"
                  f"{gauss})")
