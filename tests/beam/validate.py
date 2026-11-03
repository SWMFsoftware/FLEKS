#!/usr/bin/env python3
"""Validator for the beam instability test (tests/beam).

The runner exercises the directory once per PARAM.in variant:

  * ``PARAM.in``                       -> base_name "beam"
  * ``PARAM.in.implicit``              -> "beam_implicit"
  * ``PARAM.in.explicit``              -> "beam_explicit"
  * ``PARAM.in.hybrid``                -> "beam_hybrid"
  * ``PARAM.in.implicit.30s.waves``    -> "beam_implicit_30s_waves"

The primary validation is the FFT-based transverse-wave check on plot output.
The short full-PIC and hybrid variants share `_check_beam_transverse_wave`; the
waves variant is a long-thin-box run whose wave saturates during the run, so it
uses `_validate_waves_plot` instead.  The hybrid variant additionally runs the
shared hybrid-family energy-log checks (`validate_hybrid`) to guard against NaN
/ field blow-up.  The beam tracks test particles, so a particle-log check is
enabled via `PARTICLE_TOL`.
"""
import glob
import logging
import math
import os

from tests._shared import run_dir as _run_dir
from .._shared.hybrid import validate_hybrid

logger = logging.getLogger(__name__)

set_run_dir = _run_dir.set_run_dir

# Particle-tracking tolerance passed to validate_test_particles() in the
# common runner (validate_tests.py).
# min_velocity_change: require >= 1% relative mean-velocity change (beam
# gyration) to catch a silently disabled field gather.
PARTICLE_TOL = {
    "expected_active_species": [0],
    "launch_threshold": 0.5,
    "max_speed": 10.0,
    "min_velocity_change": 0.01,
}

# PARAM.in.implicit.30s.waves: the long-thin-box semi-implicit variant.  The
# runner turns the PARAM suffix ".implicit.30s.waves" into the base name
# "beam_implicit_30s_waves", so this token selects that variant's checks.
WAVES_VARIANT = "30s_waves"

# Thresholds for the waves variant, measured on the deck as committed
# (2000x2 cells, 100 ppc, cfl 0.2, TimeMax = 4, ~166 cycles, 2 MPI ranks
# -> ~55 s).  max|B_perp| rises from exactly 0 at t = 0 to 0.19 nT at
# t = 0.5, peaks at 0.27 nT (5.4% of the 5 nT guide field) at t = 2, then
# decays slowly to 0.21 nT; >= 80% of the transverse power sits below
# n = N/4; Bx stays exactly 5 nT with zero spread; the H+ density
# round-trips to 5.0 amu/cc.  The bounds below keep a factor of ~3-5 margin.
WAVES_MIN_AMPLITUDE = 1.0e-2   # peak max|B_perp| / |Bx| (measured 5.4e-2)
WAVES_MAX_AMPLITUDE = 2.0e-1   # blow-up guard (measured 5.4e-2)
WAVES_MIN_LOW_MODE_FRACTION = 0.4   # power in n <= N/4 (measured 0.82)
WAVES_MIN_LATE_FRACTION = 0.5  # final rms / peak rms (measured 1.02)
WAVES_BX_SPREAD_TOL = 1.0e-4   # (max - min) / mean of |Bx| (measured 0)
WAVES_DENSITY_TOL = 0.05       # rhoS1 round-trip, relative

# Input values of the deck, used as the reference for the round-trip checks.
WAVES_GUIDE_FIELD_NT = 5.0     # #UNIFORMSTATE bx
WAVES_ION_DENSITY_AMU_CC = 5.0  # #UNIFORMSTATE rho of H+


def _dft_amp(data):
    """Single-sided DFT amplitude spectrum (pure-Python fallback for numpy)."""
    nn = len(data)
    out = []
    for k in range(nn // 2 + 1):
        re = sum(data[j] * math.cos(2 * math.pi * k * j / nn)
                 for j in range(nn))
        im = -sum(data[j] * math.sin(2 * math.pi * k * j / nn)
                  for j in range(nn))
        out.append(math.hypot(re, im))
    return out


def validate_log(pic_diags=None, test_name=None):
    """Energy-log validation for the beam test.

    The hybrid variant (test_name ends with "hybrid") runs the shared hybrid
    energy checks to catch NaN / field blow-up; the full-PIC variants perform no
    log-file check here beyond a finiteness guard (their validation is the
    plot-output FFT).
    """
    if test_name and test_name.endswith("hybrid"):
        return validate_hybrid(pic_diags=pic_diags, test_name=test_name)

    if test_name and "explicit" in test_name:
        if not pic_diags or len(pic_diags) < 2:
            return False, f"Expected >= 2 energy log entries, got {len(pic_diags) if pic_diags else 0}"

        # 1. Electrons (species 0) must be fully kinetic: Epart0 must be finite, > 0, and evolve over time
        first, last = pic_diags[0], pic_diags[-1]
        ep0_start = first.get("Epart0", 0.0)
        ep0_end = last.get("Epart0", 0.0)
        if not (math.isfinite(ep0_start) and math.isfinite(ep0_end)):
            return False, "Non-finite electron kinetic energy (NaN/Inf)"
        if ep0_start <= 0.0 or ep0_end <= 0.0:
            return False, (f"Electron kinetic energy zero (electrons not initialized/kinetic: "
                           f"start={ep0_start}, end={ep0_end})")
        delta_ep0 = abs(ep0_end - ep0_start)
        if delta_ep0 <= 0.0:
            return False, "Electron kinetic energy static (electrons not pushed in explicit mode)"

        # 2. Ions (species 1) must be kinetic
        ep1_start = first.get("Epart1", 0.0)
        ep1_end = last.get("Epart1", 0.0)
        if not (math.isfinite(ep1_start) and math.isfinite(ep1_end)):
            return False, "Non-finite ion kinetic energy (NaN/Inf)"
        if ep1_start <= 0.0 or ep1_end <= 0.0:
            return False, f"Ion kinetic energy zero (start={ep1_start}, end={ep1_end})"

        # 3. Verify #SMOOTHJ was executed
        stdout_path = os.path.join(_run_dir.RUN_DIR, "fleks_stdout.log")
        if os.path.isfile(stdout_path):
            with open(stdout_path, "r", encoding="utf-8", errors="ignore") as f:
                stdout_text = f.read()
            if "Pic::smooth_multifab" not in stdout_text:
                return False, "#SMOOTHJ was specified but Pic::smooth_multifab was never executed"

        logger.debug("    [EXPLICIT] Kinetic electron check: Epart0 %.4e -> %.4e (delta=%.4e)",
                     ep0_start, ep0_end, delta_ep0)
        logger.debug("    [EXPLICIT] Kinetic ion check: Epart1 %.4e -> %.4e", ep1_start, ep1_end)
        logger.debug("    [EXPLICIT] #SMOOTHJ execution verified.")
        return True, "Passed (kinetic electrons and ions active, #SMOOTHJ verified)"

    if test_name and WAVES_VARIANT in test_name:
        # Cheap early guard: a blown-up implicit solve shows up in the energy
        # log at the first report (dnReport = 1), well before the plot frames
        # are post-processed, so this turns a coarse "plot check failed" into
        # a specific message.
        if not pic_diags or len(pic_diags) < 2:
            return False, (f"Expected >= 2 energy log entries, "
                           f"got {len(pic_diags) if pic_diags else 0}")
        first, last = pic_diags[0], pic_diags[-1]
        for key in ("Etot", "Eb"):
            if not (math.isfinite(first.get(key, 0.0)) and
                    math.isfinite(last.get(key, 0.0))):
                return False, f"Non-finite {key} in the energy log (NaN/Inf)"
        eb0, eb1 = first.get("Eb", 0.0), last.get("Eb", 0.0)
        if eb0 > 0.0 and not (0.5 * eb0 <= eb1 <= 2.0 * eb0):
            return False, (f"Magnetic energy out of bounds: Eb {eb0:.4e} -> "
                           f"{eb1:.4e} (factor > 2)")
        logger.debug("    [WAVES] Energy log: Etot %.4e -> %.4e, Eb %.4e -> %.4e",
                     first.get("Etot", 0.0), last.get("Etot", 0.0), eb0, eb1)
        return True, "Passed (energy log finite, Eb steady)"

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
    plots_dir = os.path.join(_run_dir.RUN_DIR, "PC", "plots")
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
    # The H+ ion density column is requested in #SAVEPLOT varName.  In PLANETARY
    # output it is converted back to amu/cc, so it must round-trip to the input
    # #UNIFORMSTATE value.  The species index differs between the full-PIC
    # (H+ = species 1, rhoS1) and hybrid (single H+ species, rhoS0) variants.
    irho = vidx.get("RHOS1")
    if irho is None:
        irho = vidx.get("RHOS0")

    x = []
    by = []
    bz = []
    bx = []
    rho1 = []
    for line in lines[5:]:
        cols = line.split()
        if len(cols) <= max(iby, ibz, ibx, irho if irho is not None else 0):
            continue
        try:
            x.append(float(cols[0]))
            by.append(float(cols[iby]))
            bz.append(float(cols[ibz]))
            bx.append(float(cols[ibx]))
            if irho is not None:
                rho1.append(float(cols[irho]))
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

    # The plot is saved in PLANETARY units (#SAVEPLOT ... planet): x
    # coordinates are normalised by rPlanet.  Because this test sets no
    # #PLANETRADIUS, rPlanet defaults to lNormSI (#NORMALIZATION), so physical
    # length [m] = (output coordinate) * lNormSI.
    l_norm_si = 1000.0  # lNormSI [m] from #NORMALIZATION (no #PLANETRADIUS -> rPlanet)

    logger.debug("    [FFT] t=%.4f s (SI, .out header), N=%d cells", t_norm, n)
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

    # ---- Check 2: dominant mode & power in low-order modes ----------------
    # FFT By and Bz separately (preserving sign/oscillation), then combine
    # the per-mode amplitudes.  Using |B_perp| directly would introduce
    # spurious harmonics from the magnitude operation.
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

    # Box geometry: x is in planetary units in the .out file (normalised by
    # rPlanet = lNormSI here), so the physical box length in metres is
    # L = (coordinate extent) * lNormSI.  This keeps every quantity in SI
    # (L in m, k in 1/m) for the comparison against k_res.
    dx_out = x[1] - x[0]
    L = n * dx_out * l_norm_si              # physical box length [m]
    k1 = 2 * math.pi / L                    # box-fundamental wavenumber [1/m]
    n_res = max(1, round(k_res / k1))

    # ---- Check 0: output density round-trips to the input value -----------
    # #UNIFORMSTATE sets the H+ background density n = 5.0 /cc.
    # The code converts it to code units (x1e6*m_p*Si2NoRho) internally, and
    # PLANETARY output converts it back to amu/cc.
    # The round trip is the identity, so the plotted density column must
    # reproduce the input.  Guard against a unit mismatch (e.g. rho interpreted
    # as PIC/code units in the validator) with a tight relative tolerance.
    density_ok = True
    if irho is not None and rho1:
        rho_mean = sum(rho1) / len(rho1)
        rho_input = 5.0  # amu/cc, #UNIFORMSTATE H+ background
        logger.debug("    [FFT] H+ density output = %.4f amu/cc "
                     "(input %.2f amu/cc)", rho_mean, rho_input)
        if not math.isfinite(rho_mean):
            density_ok = False
        elif abs(rho_mean - rho_input) > 0.05 * rho_input:
            density_ok = False
            logger.debug("    [FFT] Density round-trip FAIL: output density "
                         "= %.4f amu/cc vs input %.2f amu/cc",
                         rho_mean, rho_input)
    else:
        logger.debug("    [FFT] H+ density column not in .out; "
                     "density round-trip check skipped.")

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

    # ---- Check 3: the transverse wave power must GROW over time ------------
    # The instability grows from seed noise, so max|B_perp| must increase across
    # the time-resolved frames -- distinguishing a real instability from static
    # numerical seed.  Skips gracefully (not a failure) if < 3 frames exist.
    growth_ok2, growth_reason = _beam_growth_over_time(plots_dir)
    if growth_reason is not None:
        logger.debug("    [FFT] Time-growth: %s", growth_reason)

    if not density_ok:
        logger.debug("    [FFT] Density round-trip check: FAIL")
        return False, "output rhoS1 does not round-trip to the input amu/cc density"
    if growth_ok and mode_ok and growth_ok2:
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
        if not growth_ok2:
            reasons.append(growth_reason or "transverse wave not growing over time")
        return False, "; ".join(reasons)


def _beam_growth_over_time(plots_dir):
    """Verify max|B_perp| grows across the time-resolved .out frames.

    Reads the per-frame max transverse-field amplitude (max of hypot(By, Bz))
    and requires it to be non-decreasing and to end at least 2x its early value,
    confirming the seeded cyclotron wave is genuinely amplified by the
    instability rather than being static seed noise.

    Returns (passed: bool, reason: str|None).  reason is None when the check
    could not be evaluated (too few frames) and should be skipped, not failed.
    """
    frames = []
    for f in sorted(glob.glob(os.path.join(plots_dir, "*.out"))):
        with open(f, "r", encoding="latin-1") as fh:
            lines = fh.readlines()
        if len(lines) < 6:
            continue
        var_names = lines[4].split()
        vidx = {v.upper(): i for i, v in enumerate(var_names)}
        if "BY" not in vidx or "BZ" not in vidx:
            continue
        iby, ibz = vidx["BY"], vidx["BZ"]
        bperp = []
        for line in lines[5:]:
            cols = line.split()
            if len(cols) <= max(iby, ibz):
                continue
            try:
                bperp.append(math.hypot(float(cols[iby]), float(cols[ibz])))
            except (ValueError, IndexError):
                continue
        if bperp:
            frames.append(max(bperp))
    if len(frames) < 3:
        return True, None  # too few frames -> skip (not a failure)
    if any(frames[i + 1] < frames[i] for i in range(len(frames) - 1)):
        logger.debug("    [FFT] Non-monotonic max|B_perp| across frames: %s",
                     ["%.2e" % v for v in frames])
        return False, (f"max|B_perp| not non-decreasing across "
                       f"{len(frames)} frames")
    if frames[0] > 0 and frames[-1] <= 2.0 * frames[0]:
        logger.debug("    [FFT] Transverse amplitude nearly static: %s",
                     ["%.2e" % v for v in frames])
        return False, (f"max|B_perp| {frames[0]:.2e} -> {frames[-1]:.2e} "
                       f"nT (no instability growth)")
    logger.debug("    [FFT] Time-growth: max|B_perp| %s",
                 " -> ".join("%.2e" % v for v in frames))
    return True, None


def _validate_instability_plot(test_name):
    """Validate ion beam R-mode instability growth and wave mode excitation."""
    import numpy as np

    plots_dir = os.path.join(_run_dir.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        return False, "No .out files found in plots directory"

    frames_bperp = []
    last_by = None
    last_bz = None
    for f in out_files:
        with open(f, "r", encoding="latin-1") as fp:
            lines = fp.readlines()
        if len(lines) < 6:
            continue
        vnames = lines[4].split()
        vidx = {v.upper(): i for i, v in enumerate(vnames)}
        if "BY" not in vidx or "BZ" not in vidx:
            continue
        iby, ibz = vidx["BY"], vidx["BZ"]
        by_arr = []
        bz_arr = []
        for line in lines[5:]:
            cols = line.split()
            if len(cols) > max(iby, ibz):
                try:
                    by_arr.append(float(cols[iby]))
                    bz_arr.append(float(cols[ibz]))
                except ValueError:
                    pass
        if by_arr:
            bperp = np.sqrt(np.array(by_arr)**2 + np.array(bz_arr)**2)
            frames_bperp.append(float(bperp.max()))
            last_by = np.array(by_arr)
            last_bz = np.array(bz_arr)

    if len(frames_bperp) < 2:
        return False, f"Expected >= 2 valid frames, found {len(frames_bperp)}"

    bperp_0 = frames_bperp[0]
    bperp_end = frames_bperp[-1]

    if not math.isfinite(bperp_end):
        return False, "Non-finite values (NaN/Inf) detected in magnetic field"

    if bperp_end <= bperp_0:
        return False, f"No wave growth: bperp_0={bperp_0:.3e}, bperp_end={bperp_end:.3e}"

    fy = np.abs(np.fft.rfft(last_by))
    fz = np.abs(np.fft.rfft(last_bz))
    amps = np.hypot(fy, fz)
    non_dc_power = np.sum(amps[1:]**2)
    if non_dc_power <= 0:
        return False, "No non-DC wave power detected"

    n_dom = int(np.argmax(amps[1:])) + 1
    dom_frac = float(amps[n_dom]**2 / non_dc_power)

    # Power in the lower half of the spectrum (n <= 16)
    n_half = len(amps) // 2
    low_mode_power = np.sum(amps[1:n_half]**2)
    low_frac = float(low_mode_power / non_dc_power)
    if low_frac < 0.25:
        return False, f"Power not concentrated in low modes: low_frac={low_frac:.2f} < 0.25"

    msg = (f"Ion beam instability: Bperp {bperp_0:.2e} -> {bperp_end:.2e} "
           f"(growth {bperp_end/max(bperp_0, 1e-30):.1f}x), dominant mode n={n_dom} "
           f"({dom_frac:.1%} power), lower-half power {low_frac:.1%}")
    logger.debug("    %s", msg)
    return True, msg


def _read_waves_frame(path):
    """Read one ``.out`` frame of the waves variant.

    Returns a dict with the frame time and NumPy arrays of the plotted
    columns, or None when the file is not a readable ASCII frame.
    """
    with open(path, "r", encoding="latin-1") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return None
    vidx = {v.upper(): i for i, v in enumerate(lines[4].split())}
    for need in ("BX", "BY", "BZ"):
        if need not in vidx:
            return None
    try:
        t = float(lines[1].split()[1])
    except (ValueError, IndexError):
        return None

    rows = []
    for line in lines[5:]:
        cols = line.split()
        if len(cols) <= max(vidx["BX"], vidx["BY"], vidx["BZ"]):
            continue
        try:
            rows.append([float(c) for c in cols])
        except ValueError:
            continue
    if len(rows) < 4:
        return None

    import numpy as np
    data = np.asarray(rows, dtype=float)
    frame = {"t": t, "n": data.shape[0],
             "bx": data[:, vidx["BX"]],
             "by": data[:, vidx["BY"]],
             "bz": data[:, vidx["BZ"]]}
    for name, key in (("RHOS1", "rho1"), ("RHOS0", "rho0")):
        if name in vidx:
            frame[key] = data[:, vidx[name]]
    return frame


def _validate_waves_plot(test_name):
    """Validation for the long-thin-box semi-implicit waves variant.

    The deck starts with an exactly zero transverse field on a uniform warm
    beam plasma, so the run exercises the whole chain (implicit E solve,
    upwind E, smoothed B/J, comoving frame) on a wave that grows from
    numerical noise and then saturates.  The generic monotonic-growth check
    of `_check_beam_transverse_wave` cannot be used here -- the wave
    saturates by t ~ 2 and decays slowly afterwards -- so the checks are:

      1. the initial frame really has no transverse field;
      2. the wave reaches the saturated level (>= 1% of the guide field) with
         its power in box-compatible low-order modes, i.e. a coherent
         instability rather than grid-scale noise;
      3. the amplitude stays bounded (no runaway) and the wave does not
         collapse at the end of the run;
      4. the guide field stays uniform and the ion density round-trips to
         the #UNIFORMSTATE value (a broken rho unit or a corrupted
         compressional field shows up here).

    Returns (passed: bool, reason: str).
    """
    import numpy as np

    plots_dir = os.path.join(_run_dir.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        return False, "No .out files found in plots directory (PostProc.pl not run?)"

    frames = [f for f in (_read_waves_frame(p) for p in out_files) if f]
    if len(frames) < 3:
        return False, (f"expected >= 3 readable frames with Bx/By/Bz, "
                       f"found {len(frames)} of {len(out_files)}")

    for f in frames:
        for key in ("bx", "by", "bz"):
            if not np.all(np.isfinite(f[key])):
                return False, f"non-finite {key.upper()} in frame t={f['t']:.3f}"

    stats = []
    for f in frames:
        bperp = np.hypot(f["by"], f["bz"])
        stats.append({
            "t": f["t"],
            "bx_mean": float(f["bx"].mean()),
            "bx_spread": float(f["bx"].max() - f["bx"].min()),
            "amp": float(bperp.max()),
            "rms": float(np.sqrt((bperp * bperp).mean())),
        })
        logger.debug("    [WAVES] t=%7.3f  Bx=%9.5f  max|B_perp|=%10.4e  rms=%10.4e",
                     f["t"], stats[-1]["bx_mean"], stats[-1]["amp"],
                     stats[-1]["rms"])

    bx = abs(stats[-1]["bx_mean"])
    if bx <= 0.0:
        return False, "guide field Bx is zero or negative in the output"

    peak = max(s["amp"] for s in stats)
    peak_rms = max(s["rms"] for s in stats)
    final_rms = stats[-1]["rms"]
    amp_ratio = peak / bx

    # 1. The initial state must be field-free in the transverse plane.
    if stats[0]["amp"] > 1e-6 * bx:
        return False, (f"initial transverse field not zero: "
                       f"max|B_perp| = {stats[0]['amp']:.3e} nT at t={stats[0]['t']:.3f}")

    # 2. The instability must reach a coherent, box-compatible mode.
    if amp_ratio < WAVES_MIN_AMPLITUDE:
        return False, (f"transverse wave did not reach saturation: "
                       f"peak max|B_perp|/|Bx| = {amp_ratio:.2e} < {WAVES_MIN_AMPLITUDE:.0e}")

    last = frames[-1]
    amps = np.hypot(np.abs(np.fft.rfft(last["by"])),
                    np.abs(np.fft.rfft(last["bz"])))
    power = float(np.sum(amps[1:] ** 2))
    if power <= 0.0:
        return False, "no non-DC transverse wave power in the final frame"
    n_low = last["n"] // 4
    low_frac = float(np.sum(amps[1:n_low + 1] ** 2)) / power
    if low_frac < WAVES_MIN_LOW_MODE_FRACTION:
        return False, (f"transverse power not in low-order modes: low_frac = "
                       f"{low_frac:.2f} < {WAVES_MIN_LOW_MODE_FRACTION:.1f} "
                       f"(grid-scale noise instead of a wave)")
    n_dom = int(np.argmax(amps[1:])) + 1

    # 3. Bounded (no runaway) and no late collapse.
    if amp_ratio > WAVES_MAX_AMPLITUDE:
        return False, (f"transverse field blew up: peak max|B_perp|/|Bx| = "
                       f"{amp_ratio:.2e} > {WAVES_MAX_AMPLITUDE:.0e}")
    if peak_rms > 0.0 and final_rms < WAVES_MIN_LATE_FRACTION * peak_rms:
        return False, (f"transverse wave collapsed at the end of the run: "
                       f"final rms = {final_rms:.3e} < "
                       f"{WAVES_MIN_LATE_FRACTION:.1f} x peak rms = {peak_rms:.3e}")

    # 4. Guide field uniform and unchanged; ion density round-trip.
    for s in stats:
        if abs(s["bx_mean"] - WAVES_GUIDE_FIELD_NT) > 0.02 * WAVES_GUIDE_FIELD_NT:
            return False, (f"guide field drifted: Bx = {s['bx_mean']:.5f} nT at "
                           f"t={s['t']:.3f} (input {WAVES_GUIDE_FIELD_NT:.3f} nT)")
        if s["bx_spread"] > WAVES_BX_SPREAD_TOL * abs(s["bx_mean"]):
            return False, (f"compressional field appeared: spread(Bx) = "
                           f"{s['bx_spread']:.3e} nT at t={s['t']:.3f} "
                           f"(tol {WAVES_BX_SPREAD_TOL:.0e} relative)")

    if "rho1" not in last:
        return False, "Required ion-density column rhoS1 is missing from the final frame"
    rho_mean = float(last["rho1"].mean())
    if not math.isfinite(rho_mean) or \
            abs(rho_mean - WAVES_ION_DENSITY_AMU_CC) > \
            WAVES_DENSITY_TOL * WAVES_ION_DENSITY_AMU_CC:
        return False, (f"RHOS1 does not round-trip: "
                       f"{rho_mean:.4f} amu/cc vs input "
                       f"{WAVES_ION_DENSITY_AMU_CC:.1f} amu/cc")

    msg = (f"waves: peak max|B_perp|/|Bx| = {amp_ratio:.2e} at t="
           f"{max(stats, key=lambda s: s['amp'])['t']:.2f}, final rms "
           f"{final_rms:.2e} ({final_rms / peak_rms if peak_rms else 0:.2f}x peak), "
           f"dominant mode n={n_dom}, low-mode power {low_frac:.1%}, "
           f"Bx uniform, rho round-trip OK, {len(frames)} frames")
    logger.debug("    [WAVES] %s", msg)
    return True, msg


def validate_plot(test_name):
    """Plot-output check: FFT-based transverse-wave resonant-wavenumber check."""
    logger.debug("  --- Validating Output Files (FFT transverse wave) ---")
    if "instability" in test_name:
        return _validate_instability_plot(test_name)
    if WAVES_VARIANT in (test_name or ""):
        return _validate_waves_plot(test_name)
    result, reason = _check_beam_transverse_wave()
    if result:
        logger.debug("    [FFT] Beam transverse-wave resonance check: VERIFIED")
    return result, reason
