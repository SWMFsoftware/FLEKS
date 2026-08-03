#!/usr/bin/env python3
"""Convection-wave plot validation helpers for tests/hybrid_convection_wave.

Kept in _shared because they are closely related to the other hybrid tests and
share the _frame_time() helper defined in _shared/hybrid.
"""
import logging
import math
import os
import glob

logger = logging.getLogger(__name__)

from .hybrid import RUN_DIR, _frame_time  # noqa: E402


def _dft_mode(c, k):
    """Complex DFT coefficient C_k = (1/N) sum_i c_i exp(-2*pi*i*k*i/N)."""
    import cmath
    n = len(c)
    s = 0.0 + 0.0j
    for i in range(n):
        s += c[i] * cmath.exp(-2j * math.pi * k * i / n)
    return s / n


def _conv_read_params(test_name):
    """Read bulk ux [km/s], uNormSI/lNormSI [m/s, m], Lx [code], TimeMax [s]."""
    p = os.path.join("tests", test_name, "PARAM.in")

    def numeric_after(command):
        toks = []
        capture = False
        with open(p) as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                tok = s.split()[0]
                if tok.startswith("#"):
                    capture = (tok == command)
                    continue
                if capture:
                    try:
                        toks.append(float(tok))
                    except ValueError:
                        pass
        return toks

    uxs = numeric_after("#UNIFORMSTATE")
    norms = numeric_after("#NORMALIZATION")
    geoms = numeric_after("#GEOMETRY")
    stops = numeric_after("#STOP")
    ux = uxs[1] if len(uxs) > 1 else 0.0
    lNormSI = norms[0] if len(norms) > 0 else 1.0e5
    uNormSI = norms[1] if len(norms) > 1 else 5.0e4
    Lx = (geoms[1] - geoms[0]) if len(geoms) >= 2 else 6.4
    # TimeMax lives in #STOP (MaxIter, TimeMax), not #TIMESTEPPING.
    TimeMax = stops[1] if len(stops) > 1 else 10.0
    return ux, uNormSI, lNormSI, Lx, TimeMax


def _load_convection_profile(out_file):
    """Load the x-sorted transverse field (By, Bz) profile from a hybrid .out plot.

    The 2D ascii plot writes one row per (ix, iy) cell with a header line of
    variable names (which may include X and Y coordinate columns).  We group the
    cells by their x coordinate and average By, Bz over y (the wave is uniform in
    y) to recover the 1D x-profile needed for the advection-phase check.
    Returns (by_list, bz_list) sorted by x, or None if the frame is degenerate.
    """
    try:
        with open(out_file, "r", encoding="latin-1") as f:
            lines = f.readlines()
    except OSError:
        return None
    vidx = None
    data_start = None
    for li, line in enumerate(lines):
        toks = line.split()
        up = [t.upper() for t in toks]
        if "BY" in up and "BZ" in up:
            vidx = {t: i for i, t in enumerate(up)}
            data_start = li + 1
            break
    if vidx is None:
        return None
    iby, ibz = vidx["BY"], vidx["BZ"]
    ix = vidx.get("X")
    prof = {}
    for line in lines[data_start:]:
        c = line.split()
        if len(c) <= max(iby, ibz):
            continue
        try:
            by = float(c[iby])
            bz = float(c[ibz])
            xv = float(c[ix]) if ix is not None else None
        except (ValueError, IndexError, TypeError):
            continue
        if xv is None:
            xv = float(len(prof))
        prof.setdefault(round(xv, 8), []).append((by, bz))
    if len(prof) < 4:
        return None
    xs = sorted(prof.keys())
    by = [sum(p[0] for p in prof[x]) / len(prof[x]) for x in xs]
    bz = [sum(p[1] for p in prof[x]) / len(prof[x]) for x in xs]
    return by, bz


def _convection_phase_shift(by_e, bz_e, by_l, bz_l):
    """Phase shift (rad) of the k=1 transverse DFT mode between early/late."""
    import cmath
    def c1(by, bz):
        n = len(by)
        c = [complex(by[i], bz[i]) for i in range(n)]
        return sum(c[i] * cmath.exp(-2j * math.pi * 1 * i / n)
                   for i in range(n)) / n
    C1_e = c1(by_e, bz_e)
    C1_l = c1(by_l, bz_l)
    if abs(C1_e) < 1e-9 or abs(C1_l) < 1e-9:
        return None, None
    dphi = cmath.phase(C1_l) - cmath.phase(C1_e)
    dphi = dphi - 2.0 * math.pi * round(dphi / (2.0 * math.pi))
    return dphi, abs(C1_l) / abs(C1_e)


def _check_convection_advection(test_name):
    """Verify the transverse wave is advected rigidly at the bulk-flow speed.

    With the Hall term OFF the only E-field source is the convection term
    E = -U_i x B, so the induction equation reduces to dB/dt = -U . grad B and a
    transverse wave must translate rigidly at speed U.  We seed the wave with
    #TESTCASE ConvectionWave (fill_hybrid_wave(0.2), no velocity perturbation) and
    impose a uniform bulk flow ux in #UNIFORMSTATE on a true 1D grid (32x1x1).
    The run lasts exactly ONE advection period (TimeMax = Lx*lNormSI/ux) with
    plots at quarter periods, and two complementary checks are applied:

    1. RATE: the spatial Fourier k=1 phase of the (By, Bz) profile advances by
       -kx * Ux_code * dt_code between the two consecutive plot frames with the
       SHORTEST time gap (frame times are SI seconds; dt_code = dt_SI / tNorm
       with tNorm = lNormSI/uNormSI).  With quarter-period frames this is a
       well-resolved -pi/2 that cannot wrap.
    2. RETURN: between the FIRST and LAST frames the wave has translated one
       full wavelength, so the wrapped phase shift must be ~0 and the pattern
       must coincide with the initial one.  This closes the loop the rate check
       alone leaves open (any per-frame rate error accumulates 4x here).

    Both checks require the transverse amplitude to be conserved (no damping
    or growth).
    """
    out_dir = os.path.join(RUN_DIR, "PC", "plots")
    if not os.path.isdir(out_dir):
        return False, "no plots dir (%s)" % out_dir
    out_files = sorted(glob.glob(os.path.join(out_dir, "*.out")))
    valid = []
    for f in out_files:
        prof = _load_convection_profile(f)
        if prof is None:
            continue
        t = _frame_time(f)
        if t is None:
            continue
        valid.append((t, prof[0], prof[1], f))
    if len(valid) < 2:
        return False, "need >=2 valid profile frames, found %d" % len(valid)
    valid.sort(key=lambda v: v[0])
    # Measure the phase shift over the SHORTEST consecutive-frame gap. This is
    # robust to 2*pi phase wrapping over long runs (e.g. a full-period run where
    # the first and last frames would wrap back to ~0) and keeps the shift
    # well-resolved. The shift is scaled to code time via tNorm below.
    best = None  # (gap, index)
    for i in range(1, len(valid)):
        gap = valid[i][0] - valid[i - 1][0]
        if gap > 0 and (best is None or gap < best[0]):
            best = (gap, i)
    if best is None:
        return False, "need >=2 distinct-time frames, found %d" % len(valid)
    gap, i = best
    t_e, by_e, bz_e, fe = valid[i - 1]
    t_l, by_l, bz_l, fl = valid[i]
    if len(by_e) != len(by_l):
        return False, "profile length mismatch between frames"
    dphi, amp_ratio = _convection_phase_shift(by_e, bz_e, by_l, bz_l)
    if dphi is None:
        return False, "no coherent transverse wave (|C1|~0)"
    ux_kms, uNormSI, lNormSI, Lx, _ = _conv_read_params(test_name)
    # Plot-header times are in SI seconds; kx*ux_code is per CODE time unit.
    # Convert the frame separation to code units via tNorm = lNormSI/uNormSI.
    tNorm = lNormSI / uNormSI
    dt = (t_l - t_e) / tNorm
    ux_code = ux_kms * 1000.0 / uNormSI
    kx = 2.0 * math.pi / Lx
    expected = -kx * ux_code * dt
    tol = max(0.15, 0.1 * abs(expected))
    ok_phase = abs(dphi - expected) <= tol
    ok_amp = 0.7 <= amp_ratio <= 1.3
    msg = ("phase shift dphi=%.4f rad (expected %.4f over dt_code=%.2f, tol %.3f); "
           "amp ratio=%.3f; ux=%.3f km/s -> ux_code=%.4f, Lx=%.3f"
           % (dphi, expected, dt, tol, amp_ratio, ux_kms, ux_code, Lx))
    if not ok_phase:
        msg += " [PHASE MISMATCH]"
    if not ok_amp:
        msg += " [AMPLITUDE OUT OF RANGE]"

    # RETURN check: over the whole run (first vs last frame) the WRAPPED phase
    # shift must match the wrapped expectation.  For the standard one-period
    # run (TimeMax = Lx*lNormSI/ux) the expected total is -2*pi, which wraps to
    # 0: the wave must have returned to its initial pattern.
    t_0, by_0, bz_0, f0 = valid[0]
    t_n, by_n, bz_n, fn = valid[-1]
    ok_ret = True
    if t_n > t_0 and len(by_0) == len(by_n):
        dphi_full, amp_full = _convection_phase_shift(by_0, bz_0, by_n, bz_n)
        if dphi_full is None:
            return False, msg + "; RETURN: no coherent wave in first/last frame"
        exp_full = -kx * ux_code * (t_n - t_0) / tNorm
        exp_wrap = exp_full - 2.0 * math.pi * round(exp_full / (2.0 * math.pi))
        derr = dphi_full - exp_wrap
        derr = derr - 2.0 * math.pi * round(derr / (2.0 * math.pi))
        tol_ret = max(0.15, 0.1 * abs(exp_full))
        ok_ret = abs(derr) <= tol_ret and 0.7 <= amp_full <= 1.3
        msg += ("; RETURN over [%.1f, %.1f]s: wrapped dphi=%.4f rad "
                "(expected %.4f, total %.4f, tol %.3f), amp ratio=%.3f"
                % (t_0, t_n, dphi_full, exp_wrap, exp_full, tol_ret, amp_full))
        if not ok_ret:
            msg += " [RETURN-TO-START FAILED]"

    ok = ok_phase and ok_amp and ok_ret
    return ok, msg
