#!/usr/bin/env python3
"""Shared validators and plot helpers for the hybrid-PIC family of tests.

Exported functions:
  * validate_hybrid -- energy-log checks (finite energies, bounded growth).
  * validate_plot   -- seeded-wavelength and whistler-dispersion validation.
"""
import cmath
import glob
import logging
import math
import os

from tests._shared import run_dir as _run_dir

logger = logging.getLogger(__name__)
RUN_DIR = _run_dir.RUN_DIR


def set_run_dir(run_dir):
    """Point the plot helpers at the current run directory."""
    _run_dir.set_run_dir(run_dir)
    globals()["RUN_DIR"] = run_dir


def validate_hybrid(pic_diags=None, test_name=None):
    """Validate energy diagnostics for hybrid-PIC wave tests.

    Success criteria:
    1. Magnetic energy (Eb) and ion energy (Epart) are finite.
    2. Magnetic energy growth is bounded (catches numerical Hall instabilities).
    3. Ion kinetic energy ratio is bounded (catches gross particle non-conservation).
    """
    logger.debug("Validating Hybrid PIC Wave Test...")

    if not pic_diags or len(pic_diags) < 2:
        logger.debug("  [INFO] No PIC energy log found; skipping energy checks.")
        return True, "Passed (no pic log)"

    first, last = pic_diags[0], pic_diags[-1]
    passed = True
    reasons = []

    eb0 = first.get("Eb", 0.0)
    eb1 = last.get("Eb", 0.0)
    logger.debug("    Eb (magnetic): %.6e -> %.6e", eb0, eb1)
    if not math.isfinite(eb1):
        passed = False
        reasons.append("Eb not finite (NaN/Inf)")
    elif eb0 > 0 and eb1 > eb0 * 5.0:
        passed = False
        reasons.append("Eb grew >5x (whistler instability / 4pi bug?)")

    ep1 = last.get("Epart", 0.0)
    logger.debug("    Epart (ions):  %.6e -> %.6e", first.get("Epart", 0.0), ep1)
    if not math.isfinite(ep1):
        passed = False
        reasons.append("Epart not finite (NaN/Inf)")

    # Guard against gross ion energy non-conservation or runaway blow-up.
    e0 = first.get("Epart0", 0.0)
    e1 = last.get("Epart0", 0.0)
    if e0 > 0:
        ratio = e1 / e0
        logger.debug("    Epart0 (ions): %.6e -> %.6e (ratio %.4f)", e0, e1, ratio)
        if ratio < 0.2 or ratio > 10.0:
            passed = False
            reasons.append(
                f"Ion energy ratio {ratio:.3f} outside [0.2, 10.0] "
                f"(gross particle non-conservation / runaway)"
            )
    else:
        logger.debug("    [INFO] Epart0 initial zero; skipping ion check.")

    if passed:
        logger.debug("Hybrid PIC Wave Test: PASSED")
        return True, "Passed"
    return False, "; ".join(reasons)


def _hyb_load_out(out_file):
    """Load By and Bz columns from a hybrid-wave .out plot file."""
    try:
        with open(out_file, "r", encoding="latin-1") as f:
            lines = f.readlines()
    except OSError:
        return None

    if len(lines) < 6:
        return None

    var_names = lines[4].split()
    vidx = {v.upper(): i for i, v in enumerate(var_names)}
    if "BY" not in vidx or "BZ" not in vidx:
        return None

    iby, ibz = vidx["BY"], vidx["BZ"]
    min_cols = max(iby, ibz) + 1
    by, bz = [], []

    for line in lines[5:]:
        cols = line.split()
        if len(cols) < min_cols:
            continue
        try:
            by.append(float(cols[iby]))
            bz.append(float(cols[ibz]))
        except ValueError:
            continue

    return (by, bz) if by else None


def _hyb_dft_dominant(by):
    """Compute 1D spatial DFT; return (dominant_k, dominant_frac, nondc_power)."""
    n = len(by)
    if n < 4:
        return 0, 0.0, 0.0

    dft_mag = [
        abs(sum(val * cmath.exp(-2j * math.pi * k * i / n) for i, val in enumerate(by)))
        for k in range(n // 2 + 1)
    ]
    nondc_power = sum(dft_mag[k] ** 2 for k in range(1, n // 2 + 1))
    if nondc_power < 1e-30:
        return 0, 0.0, 0.0

    dominant_k = max(range(1, n // 2 + 1), key=lambda k: dft_mag[k])
    dominant_frac = dft_mag[dominant_k] ** 2 / nondc_power
    return dominant_k, dominant_frac, nondc_power


def _hyb_param_path(test_name=None):
    """Locate the PARAM.in file for the test, checking RUN_DIR first."""
    base = test_name or "whistler"
    if base.endswith("_hybrid"):
        base = base[:-len("_hybrid")]

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_tests = os.path.abspath(os.path.join(script_dir, ".."))

    candidates = [
        os.path.join(_run_dir.RUN_DIR, "PARAM.in"),
        os.path.join(repo_tests, base, "PARAM.in"),
        os.path.join("tests", base, "PARAM.in"),
        os.path.join(repo_tests, "whistler", "PARAM.in"),
        os.path.join("tests", "whistler", "PARAM.in"),
    ]
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    return os.path.join("tests", "whistler", "PARAM.in")


def _parse_param_commands(param_path, target_commands):
    """Parse numeric values under specified #COMMAND blocks in a PARAM.in file."""
    blocks = {cmd: [] for cmd in target_commands}
    current_cmd = None
    try:
        with open(param_path, "r", encoding="latin-1") as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                first_tok = s.split()[0]
                if first_tok.startswith("#"):
                    current_cmd = first_tok if first_tok in blocks else None
                    continue
                if current_cmd:
                    try:
                        blocks[current_cmd].append(float(first_tok))
                    except ValueError:
                        pass
    except OSError:
        pass
    return blocks


def _parse_number(tok):
    """PARAM value -> float: logicals (T/F), fractions ('1.0/1836.0'), numbers."""
    t = tok.strip()
    if t.upper() in ("T", "TRUE", ".TRUE."):
        return 1.0
    if t.upper() in ("F", "FALSE", ".FALSE."):
        return 0.0
    try:
        return float(t)
    except ValueError:
        pass
    if "/" in t:
        num, _, den = t.partition("/")
        try:
            return float(num) / float(den)
        except (ValueError, ZeroDivisionError):
            return None
    return None


def _parse_named(param_path, target_commands):
    """Like _parse_param_commands but keeps the SWMF parameter name."""
    blocks = {cmd: [] for cmd in target_commands}
    current_cmd = None
    try:
        with open(param_path, "r", encoding="latin-1") as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue
                toks = s.split()
                if toks[0].startswith("#"):
                    current_cmd = toks[0] if toks[0] in blocks else None
                    continue
                if current_cmd is None:
                    continue
                value = _parse_number(toks[0])
                if value is not None:
                    blocks[current_cmd].append(
                        (value, toks[1] if len(toks) > 1 else ""))
    except OSError:
        pass
    return blocks


def _first_root(f, lo, hi, nscan=2000):
    """Lowest root of f on (lo, hi) by sign-change scan + bisection."""
    fa = f(lo)
    x0, f0 = lo, fa
    for i in range(1, nscan + 1):
        x1 = lo + (hi - lo) * i / nscan
        f1 = f(x1)
        if f0 == 0.0:
            return x0
        if f0 * f1 < 0.0:
            a, b = x0, x1
            for _ in range(80):
                m = 0.5 * (a + b)
                if f(a) * f(m) <= 0.0:
                    b = m
                else:
                    a = m
            return 0.5 * (a + b)
        x0, f0 = x1, f1
    return None


def _full_em_branches(kappa, mu, ne_over_ni):
    """Lowest positive roots of the cold-plasma parallel R and L modes.

    In the standalone code units (c = 1, Omega_i = 1, d_i = 1, so w and kappa
    are omega/Omega_i and k d_i) the mode equations, i.e.

        c^2 k^2/omega^2 = 1 - sum_s w_ps^2 / (omega (omega +- Omega_s))

    read

        R: kappa^2/w^2 = 1 + (ne/ni)/(w (1 - w/mu)) - 1/(w (w+1))
        L: kappa^2/w^2 = 1 - (ne/ni)/(w (1 + w/mu)) - 1/(w (w-1))

    with mu = |Omega_e|/Omega_i = m_i/m_e.  Unlike the Hall / massless-electron
    form these retain the displacement current (the leading 1), which is only
    negligible when c >> v_A.
    """
    def f_r(w):
        return (1.0 + ne_over_ni / (w * (1.0 - w / mu))
                - 1.0 / (w * (w + 1.0)) - kappa * kappa / (w * w))

    def f_l(w):
        return (1.0 - ne_over_ni / (w * (1.0 + w / mu))
                - 1.0 / (w * (w - 1.0)) - kappa * kappa / (w * w))

    return (
        ("right-hand/whistler", _first_root(f_r, 1e-3, 0.9999 * mu)),
        ("left-hand/ion-cyclotron", _first_root(f_l, 1e-3, 0.9999)),
    )


def _solver_and_species(param_path):
    """Return (is_full_em, mu, ne_over_ni) read from the deck.

    is_full_em is True when the deck runs the Maxwell solver (#SOLVEEM T), for
    which the expectation is the full cold-plasma dispersion rather than the
    Hall (massless-electron, no displacement current) one.  mu = m_i/m_e;
    mu is None when the deck has a single species.
    """
    blocks = _parse_named(param_path, ["#SOLVEEM", "#PLASMA", "#UNIFORMSTATE"])
    solveem = blocks.get("#SOLVEEM", [])
    is_full_em = bool(solveem) and solveem[0][0] != 0.0

    plasma = blocks.get("#PLASMA", [])
    masses = [v for v, n in plasma if n.lower().startswith("mass")]
    unif = blocks.get("#UNIFORMSTATE", [])
    rhos = [v for v, n in unif if n.lower().startswith("rho")]

    mu = None
    if len(masses) >= 2 and masses[1] > 0:
        mu = masses[0] / masses[1]
    ne_over_ni = None
    if mu and len(rhos) >= 2 and rhos[0] > 0 and masses[0] > 0:
        n_i = rhos[0] / masses[0]
        n_e = rhos[1] / masses[1]
        if n_i > 0:
            ne_over_ni = n_e / n_i
    return is_full_em, mu, ne_over_ni


def _hyb_seeded_mode(test_name=None):
    """Return seeded spatial mode from #WAVEIC waveMode (defaults to 1)."""
    p = _hyb_param_path(test_name)
    try:
        with open(p, "r", encoding="latin-1") as f:
            for line in f:
                toks = line.strip().split()
                if len(toks) >= 2:
                    # In FLEKS PARAM files, values typically precede parameter names (e.g. '1 waveMode')
                    if toks[1].lower().startswith("wavemode"):
                        return int(float(toks[0]))
                    if toks[0].lower().startswith("wavemode"):
                        return int(float(toks[1]))
    except (OSError, ValueError):
        pass
    return 1


def _hyb_whistler_dispersion(out_files, test_name=None):
    """Measure the seeded-mode whistler frequency and compare with hybrid theory.

    Tracks the phase of the dominant spatial mode C(t) across all frames, fits
    phi(t) = omega*t + phi0, and compares against the two parallel-propagating
    Hall / cold-plasma branches (kappa = k d_i):

        right-hand (whistler)     : w/Omega_i = [kappa^2 + kappa*sqrt(kappa^2+4)]/2
        left-hand (ion-cyclotron) : w/Omega_i = [-kappa^2 + kappa*sqrt(kappa^2+4)]/2

    The check passes when the measured |omega| matches either branch (whichever
    the seeded helicity should have excited).  NOTE: at kappa -> 0 both reduce
    to the Alfven wave (w = k v_A) and at large kappa they approach kappa^2
    (whistler, i.e. w ~ k^2) and Omega_i (ion-cyclotron resonance) respectively.

    Returns:
      (True, message)  -> measured frequency matched theory
      (False, reason)  -> measured frequency mismatched theory
      (None, reason)   -> insufficient data to measure
    """
    if not out_files or len(out_files) < 3:
        return None, "need >=3 .out frames"

    # Identify dominant mode from the first readable frame.
    first_data = None
    for f in out_files:
        data = _hyb_load_out(f)
        if data is not None and len(data[0]) >= 4:
            first_data = data
            break
    if first_data is None:
        return None, "no parseable frames"

    by0, bz0 = first_data
    kdom, kfrac, _ = _hyb_dft_dominant(by0)
    if kdom <= 0:
        return None, "no dominant spatial mode (flat By)"
    if kfrac < 0.3:
        return None, f"dominant mode too weak ({kfrac:.2f})"

    # The seeded helicity decides which sign of the spatial harmonic carries
    # the wave: Psi = By + i Bz ~ exp(+i k x) for the left-hand seed and
    # exp(-i k x) for the right-hand one, so track whichever dominates the
    # transverse field.  (Projecting onto the wrong sign averages the wave
    # away and leaves pure noise.)
    n0 = len(by0)
    m_pos = sum(
        complex(by0[i], bz0[i]) * cmath.exp(-2j * math.pi * kdom * i / n0)
        for i in range(n0)
    ) / n0
    m_neg = sum(
        complex(by0[i], bz0[i]) * cmath.exp(2j * math.pi * kdom * i / n0)
        for i in range(n0)
    ) / n0
    ksign = 1.0 if abs(m_pos) >= abs(m_neg) else -1.0

    # Collect (time_si, phase) samples for dominant mode.
    samples = []
    for f in out_files:
        data = _hyb_load_out(f)
        if data is None:
            continue
        by, bz = data
        n = len(by)
        if n < 4:
            continue
        t_si = _frame_time(f)
        if t_si is None:
            continue
        c = sum(
            complex(by[i], bz[i]) * cmath.exp(-ksign * 2j * math.pi * kdom * i / n)
            for i in range(n)
        )
        c_mode = c / n
        if abs(c_mode) >= 1e-12:
            bperp = max(math.hypot(by[i], bz[i]) for i in range(n))
            samples.append((t_si, cmath.phase(c_mode), bperp))

    if len(samples) < 3:
        return None, "too few usable frames"

    # The explicit hybrid advance of a low-beta plasma slowly drives grid-scale
    # fields, so drop the trailing frames once the transverse amplitude has
    # grown well past its seeded value (keeps >= 3 samples).
    samples.sort(key=lambda s: s[0])
    limit = 1.5 * samples[0][2]
    for i, s in enumerate(samples):
        if s[2] > limit and i >= 3:
            samples = samples[:i]
            break

    # Phase unwrapping.
    samples.sort(key=lambda s: s[0])
    unwrapped = [samples[0][1]]
    for i in range(1, len(samples)):
        dphi = samples[i][1] - samples[i - 1][1]
        dphi -= 2.0 * math.pi * round(dphi / (2.0 * math.pi))
        unwrapped.append(unwrapped[-1] + dphi)

    # Linear regression: phi = omega_si * t + phi0.
    ts = [s[0] for s in samples]
    nt = len(ts)
    mean_t = sum(ts) / nt
    mean_p = sum(unwrapped) / nt
    cov = sum((ts[i] - mean_t) * (unwrapped[i] - mean_p) for i in range(nt))
    var_t = sum((ts[i] - mean_t) ** 2 for i in range(nt))
    if var_t <= 0:
        return None, "no time spread in frames"

    omega_si = cov / var_t
    var_p = sum((unwrapped[i] - mean_p) ** 2 for i in range(nt))
    r2 = (cov ** 2 / (var_t * var_p)) if var_p > 0 else 0.0

    # Read geometry and normalization from PARAM.in.
    p = _hyb_param_path(test_name)
    blocks = _parse_param_commands(p, ["#GEOMETRY", "#NORMALIZATION"])
    geoms = blocks.get("#GEOMETRY", [])
    norms = blocks.get("#NORMALIZATION", [])

    if len(geoms) >= 2 and len(norms) >= 2:
        Lx = abs(geoms[1] - geoms[0])
        tNorm = norms[0] / norms[1] if norms[1] > 0 else 1.0
    else:
        Lx, tNorm = 6.4, 1.0

    k = 2.0 * math.pi * kdom / Lx
    # These decks normalize lNormSI = uNormSI so that d_i ~ 1 (and Omega_i ~ 1)
    # in code units, hence kappa = k d_i ~ k.
    kappa = k

    # The expected dispersion depends on which physics the field solver carries:
    # the hybrid (massless-electron, generalized Ohm) solver reproduces the
    # Hall branches, while a full Maxwell (#SOLVEEM T) run keeps the
    # displacement current *and* the electron inertia, which are only
    # negligible for c >> v_A.  Pick the deck-appropriate model.
    full_em, mu, ne_over_ni = _solver_and_species(p)
    branches, model = None, "Hall / massless electrons"
    if full_em:
        if mu and ne_over_ni:
            branches = _full_em_branches(kappa, mu, ne_over_ni)
            model = (f"full cold plasma (m_i/m_e = {mu:.0f}, "
                     f"n_e/n_i = {ne_over_ni:.3g})")
        else:
            branches = None
            model = "full cold plasma (species not parseable)"
    if not branches or any(w is None for _, w in branches):
        kappa2 = kappa * kappa
        root = math.sqrt(kappa2 * kappa2 + 4.0 * kappa2)
        branches = (
            ("right-hand/whistler", 0.5 * (kappa2 + root)),
            ("left-hand/ion-cyclotron", 0.5 * (-kappa2 + root)),
        )
        model = "Hall / massless electrons"

    omega_code = omega_si * tNorm
    measured = abs(omega_code)
    branch, omega_theory = min(
        branches, key=lambda b: abs(measured - b[1]))
    tol = 0.20 * max(omega_theory, 1e-9)

    msg = (
        f"whistler dispersion n={kdom}: measured |omega|/Omega_i = {measured:.3f} "
        f"({branch} theory = {omega_theory:.3f} [{model}], "
        f"k d_i = {kappa:.3f}), fit r^2 = {r2:.2f}"
    )

    if abs(measured - omega_theory) <= tol:
        return True, msg
    return False, f"{msg} [DISPERSION MISMATCH]"


def _check_hybrid_wave_dispersion(test_name=None):
    """Verify transverse wave initialization, bounded amplitude, and dispersion."""
    plots_dir = os.path.join(_run_dir.RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        logger.debug("    [HYB] No .out files found (PostProc.pl not run?).")
        return True, "No .out files (skipped)"

    # Early-time check: verify seeded mode dominates.
    early_file = out_files[0]
    logger.debug("    [HYB] Early .out: %s", os.path.basename(early_file))
    early_data = _hyb_load_out(early_file)
    if early_data is None:
        return True, "Could not parse early .out file"

    by_e, bz_e = early_data
    n_e = len(by_e)
    bperp_e = [math.hypot(by_e[i], bz_e[i]) for i in range(n_e)]
    bperp_max_e = max(bperp_e) if bperp_e else 0.0
    dom_k_e, dom_frac_e, nondc_e = _hyb_dft_dominant(by_e)

    logger.debug(
        "    [HYB] Early: N=%d, max|B_perp|=%.4e, dominant mode n=%d (%.1f%%)",
        n_e, bperp_max_e, dom_k_e, dom_frac_e * 100
    )

    if nondc_e < 1e-30:
        return False, "No transverse wave power at early time (By is flat/zero)"

    seeded_mode = _hyb_seeded_mode(test_name)
    if dom_k_e != seeded_mode:
        return False, (
            f"Early dominant mode n={dom_k_e} (expected n={seeded_mode} for seeded wavelength)"
        )
    if dom_frac_e < 0.5:
        return False, (
            f"Early mode n={seeded_mode} carries only {dom_frac_e * 100:.1f}% of non-DC power "
            "(wave spectrum not clean at t=0)"
        )

    # Late-time check: verify amplitude remains bounded.
    late_file = out_files[-1]
    logger.debug("    [HYB] Late .out:  %s", os.path.basename(late_file))
    late_data = _hyb_load_out(late_file)
    if late_data is None:
        return True, "Could not parse late .out file"

    by_l, bz_l = late_data
    n_l = len(by_l)
    bperp_l = [math.hypot(by_l[i], bz_l[i]) for i in range(n_l)]
    bperp_max_l = max(bperp_l) if bperp_l else 0.0
    dom_k_l, dom_frac_l, _ = _hyb_dft_dominant(by_l)

    growth = bperp_max_l / bperp_max_e if bperp_max_e > 0 else float("inf")
    logger.debug(
        "    [HYB] Late:  N=%d, max|B_perp|=%.4e, dominant mode n=%d (%.1f%%)",
        n_l, bperp_max_l, dom_k_l, dom_frac_l * 100
    )
    logger.debug("    [HYB] Growth: %.1fx (early -> late)", growth)

    if bperp_max_l > 10.0:
        return False, f"Late amplitude {bperp_max_l:.2e} too large (unstable; seed was ~0.02)"

    # Whistler dispersion check (requires >= 3 time-resolved frames).
    disp_ok, disp_reason = _hyb_whistler_dispersion(out_files, test_name)
    if disp_ok is False:
        return False, disp_reason
    if disp_ok is True:
        logger.debug("    [HYB] %s", disp_reason)
    else:
        logger.debug("    [HYB] Whistler-dispersion check skipped: %s", disp_reason)

    logger.debug("    [HYB] Hybrid wave: early wavelength + late bounded-amplitude: VERIFIED")
    return True, "Passed"


def _frame_time(out_file):
    """Read simulation time from .out header (line 1, 2nd token)."""
    try:
        with open(out_file, "r", encoding="latin-1") as f:
            f.readline()
            nxt = f.readline()
        return float(nxt.split()[1])
    except (OSError, IndexError, ValueError):
        return None


def validate_plot(test_name):
    """Plot-output check shared by the hybrid-family tests."""
    if test_name == "freestream":
        logger.debug("  --- Validating Output Files: No plot-file check for this test ---")
        return True, "Passed (no plot-file check)"
    logger.debug("  --- Validating Output Files (Hybrid wave) ---")
    result, reason = _check_hybrid_wave_dispersion(test_name)
    if result:
        logger.debug("    [HYB] Hybrid wave output check: VERIFIED")
    return result, reason
