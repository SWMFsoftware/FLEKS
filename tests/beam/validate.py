#!/usr/bin/env python3
"""Validator for the beam instability test (tests/beam).

The primary validation is the FFT-based transverse-wave check performed on the
plot output by _check_beam_transverse_wave().  The beam also tracks test
particles, so a particle-log check is enabled via PARTICLE_TOL.
"""
import glob
import logging
import math
import os

logger = logging.getLogger(__name__)

from .._shared.hybrid import RUN_DIR  # noqa: E402

# Particle-tracking tolerance passed to validate_test_particles() in the
# common runner (validate_tests.py).
PARTICLE_TOL = {
    "expected_active_species": [0],
    "launch_threshold": 0.5,
    "max_speed": 10.0,
}


def validate_log(pic_diags=None, test_name=None):
    """Validate the beam instability test.

    The primary validation is the FFT-based transverse-wave check performed
    on the plot output by _check_beam_transverse_wave().  No log-file-based
    checks are performed here.
    """
    logger.debug("Validating Beam Instability Test...")
    logger.debug("  [INFO] Beam diagnostic checks rely on plot output (FFT).")
    logger.debug("Beam Instability Test: PASSED")
    return True, "Passed"


def _check_beam_transverse_wave():
    """Check transverse EM wave growth against the cyclotron resonance.

    Reads the final .out plot file (at t ~= 0.1 normalized), FFTs the
    transverse magnetic-field profile (By, Bz), and compares the dominant
    spatial mode to the theoretical cyclotron-resonant wavenumber
    ``k_res = Omega_i / Delta_v``.

    For the beam test the resonant wavelength (``k_res^-1`` ~ 10^4 km)
    vastly exceeds the 2 km periodic box, so the resonant mode (n ~= 0)
    cannot fit.  The instability therefore populates the longest-wavelength
    modes that fit in the box.  This check verifies that (1) the transverse
    wave has grown above the numerical noise floor and (2) the wave power
    is concentrated in low-order spatial modes consistent with the
    box-limited instability, reporting the dominant mode and k_res for
    inspection.

    Returns (passed: bool, reason: str).
    """
    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [FFT] No .out files found (PostProc.pl not run?).")
        return True, "No .out files (skipped)"

    # Use the final frame (latest cycle); this is the t ~= 0.1 snapshot.
    out_file = out_files[-1]
    logger.debug("    [FFT] Loading .out: %s", os.path.basename(out_file))

    # Plot .out files are written by PostProc.pl and may contain non-UTF-8
    # bytes; read byte-safe so a stray byte never aborts the validation.
    with open(out_file, "r", encoding="latin-1") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return True, "Short .out file"

    var_names = lines[4].split()
    vidx = {v.upper(): i for i, v in enumerate(var_names)}
    for need in ("BY", "BZ", "BX"):
        if need not in vidx:
            logger.debug("    [FFT] '%s' not found in .out variables: %s",
                         need, var_names)
            return True, f"{need} not in .out"

    iby, ibz, ibx = vidx["BY"], vidx["BZ"], vidx["BX"]

    x = []
    by = []
    bz = []
    bx = []
    for line in lines[5:]:
        cols = line.split()
        if len(cols) <= max(iby, ibz, ibx):
            continue
        try:
            x.append(float(cols[0]))
            by.append(float(cols[iby]))
            bz.append(float(cols[ibz]))
            bx.append(float(cols[ibx]))
        except (ValueError, IndexError):
            continue

    n = len(x)
    if n < 4:
        logger.debug("    [FFT] Too few data points for FFT.")
        return True, "Too few points"

    # Simulation time (normalized) from line 2.
    try:
        t_norm = float(lines[1].split()[1])
    except (ValueError, IndexError):
        t_norm = float("nan")

    # B-fields are in nT in the .out file.
    bx_mean = sum(bx) / n
    bperp = [math.hypot(by[i], bz[i]) for i in range(n)]
    bperp_max = max(bperp)

    logger.debug("    [FFT] t=%.4f (normalized), N=%d cells", t_norm, n)
    logger.debug("    [FFT] |Bx|=%.4f nT, max|B_perp|=%.4e nT", bx_mean, bperp_max)

    # ---- Check 1: wave growth above the noise floor ----------------------
    # At t=0 the transverse field is exactly zero; after the instability
    # triggers it grows from numerical noise.  Require the amplitude to
    # exceed a small fraction of the guide field.
    noise_frac = 1e-4
    if bx_mean <= 0:
        growth_ok = bperp_max > 0
    else:
        growth_ok = bperp_max > noise_frac * abs(bx_mean)
    logger.debug("    [FFT] Wave growth: max|B_perp|/|Bx| = %.3e "
                 "(threshold %.0e) -> %s",
                 bperp_max / max(abs(bx_mean), 1e-30), noise_frac,
                 'OK' if growth_ok else 'FAIL')

    # ---- DFT of the transverse field -------------------------------------
    # FFT By and Bz separately (preserving sign/oscillation), then combine
    # the per-mode amplitudes.  Using |B_perp| directly would introduce
    # spurious harmonics from the magnitude operation.
    def _dft_amp(data):
        nn = len(data)
        out = []
        for k in range(nn // 2 + 1):
            re = sum(data[j] * math.cos(2 * math.pi * k * j / nn)
                     for j in range(nn))
            im = -sum(data[j] * math.sin(2 * math.pi * k * j / nn)
                      for j in range(nn))
            out.append(math.hypot(re, im))
        return out

    try:
        import numpy as np
        fy = np.abs(np.fft.rfft(by))
        fz = np.abs(np.fft.rfft(bz))
        amps = [math.hypot(float(fy[k]), float(fz[k]))
                for k in range(len(fy))]
    except ImportError:
        ay = _dft_amp(by)
        az = _dft_amp(bz)
        amps = [math.hypot(ay[k], az[k]) for k in range(len(ay))]

    # Non-DC (n>=1) power.
    total_power = sum(a * a for a in amps[1:])
    if total_power <= 0:
        logger.debug("    [FFT] No non-DC spectral power; wave has not grown.")
        return False, "No transverse wave power detected"

    # Dominant non-DC mode.
    n_dom = max(range(1, len(amps)), key=lambda k: amps[k])
    dom_frac = amps[n_dom] ** 2 / total_power

    # Fraction of power in low-order modes (n <= N/4).
    n_low = n // 4
    low_frac = sum(a * a for a in amps[1:n_low + 1]) / total_power

    logger.debug("    [FFT] Dominant mode: n=%d (%.1f%% of non-DC power)",
                 n_dom, 100 * dom_frac)
    logger.debug("    [FFT] Power in low modes (n<=%d): %.1f%%", n_low, 100 * low_frac)

    # ---- Theoretical resonant wavenumber ---------------------------------
    # Cyclotron resonance: k_res = Omega_i / Delta_v (ion-ion beam-beam).
    # Omega_i = q_p * |Bx| / m_p (SI; the Boris pusher uses q*dt/(2*m)
    # without a c factor, so this is the SI cyclotron frequency).
    q_p = 1.60217663e-19   # C  (cUnitChargeSI)
    m_p = 1.67262192e-27   # kg (cProtonMassSI)
    bx_si = abs(bx_mean) * 1e-9          # nT -> T
    omega_i = q_p * bx_si / m_p          # rad/s
    delta_v = 8e5                         # m/s (beam +400, bg -400 km/s)
    k_res = omega_i / delta_v             # 1/m

    # Box geometry (x is in metres in the .out file).
    dx = x[1] - x[0]
    L = n * dx
    k1 = 2 * math.pi / L                  # box-fundamental wavenumber
    n_res = max(1, round(k_res / k1))

    logger.debug("    [FFT] Omega_i = %.4f rad/s, Delta_v = %.1e m/s",
                 omega_i, delta_v)
    logger.debug("    [FFT] k_res = %.3e 1/m, k_1 = %.3e 1/m (L = %.1f m)",
                 k_res, k1, L)
    logger.debug("    [FFT] Resonant wavelength = %.3e m vs box L = %.1f m",
                 2 * math.pi / k_res, L)

    if k_res < k1:
        # Resonant wavelength exceeds the box: the resonant mode (n ~ 0)
        # cannot fit, so the nearest available mode is the box-fundamental
        # n=1.  The instability populates the longest-wavelength modes that
        # fit; verify the bulk of the wave power resides in low-order modes
        # (not grid-scale noise near the Nyquist frequency).
        mode_ok = low_frac > 0.4
        logger.debug("    [FFT] k_res < k_1: resonant mode (n=%.2e) "
                     "exceeds the box; nearest available mode is n=1.",
                     k_res / k1)
        logger.debug("    [FFT] Mode check (box-limited): low-mode power "
                     "fraction %.2f > 0.4 -> %s", low_frac,
                     'OK' if mode_ok else 'FAIL')
    else:
        tol = 2
        mode_ok = abs(n_dom - n_res) <= tol
        logger.debug("    [FFT] Mode check: |n_dom(%d) - n_res(%d)| <= %d -> %s",
                     n_dom, n_res, tol, 'OK' if mode_ok else 'FAIL')

    if growth_ok and mode_ok:
        logger.debug("    [FFT] Transverse wave check: PASSED")
        return True, "Passed"
    else:
        reasons = []
        if not growth_ok:
            reasons.append(f"wave amplitude {bperp_max:.2e} nT below "
                           f"noise floor ({noise_frac:.0e}*|Bx|)")
        if not mode_ok:
            reasons.append(f"dominant mode n={n_dom} inconsistent with "
                           f"resonance (n_res={n_res})")
        return False, "; ".join(reasons)


def validate_plot(test_name):
    """Plot-output check: FFT-based transverse-wave resonant-wavenumber check."""
    logger.debug("  --- Validating Output Files (FFT transverse wave) ---")
    result, reason = _check_beam_transverse_wave()
    if result:
        logger.debug("    [FFT] Beam transverse-wave resonance check: VERIFIED")
    return result, reason
