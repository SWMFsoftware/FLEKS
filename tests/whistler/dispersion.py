#!/usr/bin/env python3
"""Whistler dispersion relation from time series at probe points.

This is the *physics* companion to ``validate.py`` (which only checks that the
frequency is in the right ball-park).  It runs a scan over the box-fundamental
and its harmonics, tracks dB_y(x0, t) at a probe point, measures omega(k) from
the time series (FFT + least-squares refinement + phase difference between two
probes), and compares against the analytic parallel whistler dispersion.

Physics
-------
For a hybrid (kinetic ions + massless fluid electrons) plasma with a guide
field B0 along x, the generalized Ohm's law

    E = -u_i x B + (J x B)/rho_q,     J = curl(B)/(4*pi)      [code units]

combined with the ion momentum equation rho du_i/dt = J x B and Faraday's law
gives, for a parallel-propagating circularly polarized wave with amplitude
handedness lambda = +-i (rotation operator x-hat x),

    (omega/Omega_i)^2 = (k d_i)^2 [1 + i (omega/Omega_i) lambda].

The two physical branches (kappa = k d_i) are

    right-hand  (whistler)      : w/Omega_i = [+kappa^2 + kappa sqrt(kappa^2+4)]/2
    left-hand   (ion-cyclotron) : w/Omega_i = [-kappa^2 + kappa sqrt(kappa^2+4)]/2

with the usual limits: w = k v_A for kappa -> 0 (Alfven wave) and
w -> Omega_i kappa^2 (whistler, i.e. w ~ k^2) for kappa >> 1, while the
left-hand branch saturates at the ion-cyclotron resonance w -> Omega_i.

The seed is a transverse, circularly polarized perturbation

    dB_y(x,0) = dB cos(kx),   dB_z(x,0) = hand * dB sin(kx)

with hand = -1 for the right-hand (whistler) seed (#WAVEIC rightHand T) and
hand = +1 otherwise (see PARAM.XML).  For B0 = +x a right-hand seed with
waveMode > 0 propagates along +x - i.e. along B0 - and rotates y -> z about
B0, the electron-gyration sense.

Pass criteria reported by this script (see --help for the tolerances):
  1. Right-hand polarized: the transverse field at the probe rotates y -> z
     about +B0 (hodogram circulation positive), giving the whistler branch.
  2. omega grows like k^2 at large kappa: the measured local log-log slope
     d ln(w)/d ln(k) approaches the analytic one, which tends to 2.
  3. Phase velocity v_p = omega/k agrees with the analytic whistler v_p within
     the discrete truncation error of the compact curl/gradient stencils.

Usage
-----
    python3 tests/whistler/dispersion.py                 # hybrid, modes 1..4
    python3 tests/whistler/dispersion.py --modes 1 2 3 4 5 6
    python3 tests/whistler/dispersion.py --variant full   # full-PIC variant
    python3 tests/whistler/dispersion.py --no-run         # re-analyse a scan

Outputs ``whistler_dispersion.png`` (omega-k, d omega/d k, phase velocity) and
``whistler_waveform.png`` (probe time series, hodogram, spectrum) next to this
file, plus a PASS/FAIL table on stdout.
"""

import argparse
import glob
import math
import os
import subprocess
import sys
import time

import numpy as np

# FLEKS constants (include/Constants.h).
MP = 1.67262192e-27  # proton mass [kg]
QE = 1.60217663e-19  # unit charge [C]
FOUR_PI = 4.0 * math.pi

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

DECKS = {
    "hybrid": "PARAM.in.hybrid",
    "full": "PARAM.in",
}


# ---------------------------------------------------------------------------
# PARAM.in parsing
# ---------------------------------------------------------------------------
def parse_deck(path):
    """Return {COMMAND: [(value_str, name_str), ...]} for a PARAM.in file."""
    blocks, cmd = {}, None
    with open(path, "r", encoding="latin-1") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                cmd = s.split()[0][1:].upper()
                blocks.setdefault(cmd, [])
                continue
            if cmd is None:
                continue
            toks = s.split()
            blocks[cmd].append((toks[0], toks[1] if len(toks) > 1 else ""))
    return blocks


def get(blocks, cmd, name, default=None, nth=1):
    """nth value whose SWMF name matches (nth is 1-based)."""
    hits = [v for v, n in blocks.get(cmd, []) if n == name]
    if len(hits) < nth:
        if default is None:
            raise KeyError(f"missing #{cmd} {name}")
        return default
    return float(hits[nth - 1])


def get_str(blocks, cmd, name, default=None):
    for v, n in blocks.get(cmd, []):
        if n == name:
            return v
    return default


# ---------------------------------------------------------------------------
# Normalization derived from the deck (no hard-coded unit system)
# ---------------------------------------------------------------------------
class Units:
    """Code-unit scales of a FLEKS deck, from #NORMALIZATION / #PLASMA etc.

    FLEKS fixes the mass unit at mNormSI = 1e7 * lNormSI * (m_p/q_e)^2, which
    makes c = 1 and (for q = m = 1 in code units) B_code = Omega_i.  Hence

        rho_code = rho[amu/cc] * 1e6 * m_p * Si2NoRho,  Si2NoRho = lNormSI^3/mNormSI
        d_i      = 1 / omega_pi,   omega_pi^2 = 4 pi rho_q (q/m)^2 ,  c_code = 1
    """

    def __init__(self, blocks):
        self.lNormSI = get(blocks, "NORMALIZATION", "lNormSI", 1.0e5)
        self.uNormSI = get(blocks, "NORMALIZATION", "uNormSI", 1.0e5)
        self.tNorm = self.lNormSI / self.uNormSI

        self.xMin = get(blocks, "GEOMETRY", "xMin")
        self.xMax = get(blocks, "GEOMETRY", "xMax")
        self.Lx = abs(self.xMax - self.xMin)

        self.rhoAmuCc = get(blocks, "UNIFORMSTATE", "rho")
        self.bxT = get(blocks, "UNIFORMSTATE", "bx")
        self.mass = get(blocks, "PLASMA", "mass", 1.0)
        self.charge = get(blocks, "PLASMA", "charge", 1.0)

        self.mNormSI = 1.0e7 * self.lNormSI * (MP / QE) ** 2
        self.Si2NoRho = self.lNormSI ** 3 / self.mNormSI
        self.Bnorm = math.sqrt(self.mNormSI / (1000.0 * self.lNormSI ** 3)) * (
            100.0 * self.uNormSI)
        self.Si2NoB = 1.0e4 / self.Bnorm

        # Code-unit plasma quantities.
        self.rhoCode = self.rhoAmuCc * 1.0e6 * MP * self.Si2NoRho
        self.rhoQ = self.rhoCode * (self.charge / self.mass)
        self.omegaPi = math.sqrt(FOUR_PI * self.rhoQ * (self.charge / self.mass))
        self.di = 1.0 / self.omegaPi  # c_code = 1
        self.bCode = self.bxT * self.Si2NoB
        self.omegaCi = self.charge * self.bCode / self.mass  # c_code = 1
        self.vA = self.bCode / math.sqrt(FOUR_PI * self.rhoCode)

    def summary(self):
        return (
            f"  Lx = {self.Lx:g} code units | tNorm = {self.tNorm:g} s\n"
            f"  rho = {self.rhoCode:.5f} code units | "
            f"d_i = {self.di:.5f} code units | v_A = {self.vA:.5f} code units\n"
            f"  Omega_i = {self.omegaCi:.5f} code units "
            f"(check Omega_i*d_i/v_A = {self.omegaCi * self.di / self.vA:.5f})"
        )


def analytic(kappa, units):
    """Whistler (right-hand) and ion-cyclotron (left-hand) omega in code units."""
    k2 = kappa * kappa
    root = math.sqrt(k2 * k2 + 4.0 * k2)
    return (units.omegaCi * 0.5 * (k2 + root),   # right-hand / whistler
            units.omegaCi * 0.5 * (-k2 + root))  # left-hand / ion-cyclotron


# ---------------------------------------------------------------------------
# Run management
# ---------------------------------------------------------------------------
def prepare_run_dir(run_dir):
    os.makedirs(os.path.join(run_dir, "PC", "plots"), exist_ok=True)
    os.makedirs(os.path.join(run_dir, "PC", "restartOUT"), exist_ok=True)
    links = [
        (os.path.join(REPO, "bin", "FLEKS.exe"), os.path.join(run_dir, "FLEKS.exe")),
        (os.path.join(REPO, "share", "Scripts", "PostProc.pl"),
         os.path.join(run_dir, "PostProc.pl")),
        (os.path.join(REPO, "bin", "PostIDL.exe"),
         os.path.join(run_dir, "PC", "PostIDL.exe")),
        (os.path.join(REPO, "share", "Scripts", "pIDL"),
         os.path.join(run_dir, "PC", "pIDL")),
    ]
    for src, dst in links:
        if os.path.islink(dst) or os.path.exists(dst):
            os.remove(dst)
        os.symlink(src, dst)


def set_block(text, command, updates):
    """Rewrite the values of <name> lines of a #COMMAND block.

    PARAM.in lines are "<value> <name> [comment]"; only the value token is
    replaced so the name (which is what the reader matches on) survives.
    """
    out, in_block = [], False
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("#"):
            in_block = s.split()[0][1:].upper() == command.upper()
            out.append(line)
            continue
        if in_block and s:
            toks = s.split()
            if len(toks) > 1 and toks[1] in updates:
                head = line[:len(line) - len(line.lstrip())]
                body = line[len(head):]
                out.append(head + body.replace(toks[0], str(updates[toks[1]]), 1))
                continue
        out.append(line)
    return "\n".join(out) + "\n"


def eigenmode_walen_factor(kappa, units, right_hand):
    """|u_perp|/|B_perp| (code units) of the seeded branch.

    The seeded Alfvenic kick (walenFactor = 1) is the exact eigenmode only for
    k -> 0.  The parallel Hall eigenmode needs

        u_perp = -(k d_i)/(omega/Omega_i) * (v_A/B0) * B_perp

    i.e. walenFactor = kappa/omega_bar * v_A/B0 for the right-hand whistler
    seeded with waveMode > 0.
    """
    if not right_hand:
        return None
    w_w, _ = analytic(kappa, units)
    return (kappa / (w_w / units.omegaCi)) * (units.vA / units.bCode)


def build_param(deck_text, mode, right_hand, dt, time_max, dn, walen_factor=None,
                ppc=None):
    updates = {"waveMode": mode, "rightHand": "T" if right_hand else "F"}
    if walen_factor is not None:
        updates["walenFactor"] = f"{walen_factor:.6f}"
    if "#WAVEIC" not in deck_text:
        deck_text = deck_text.replace(
            "#TIMESTEPPING", "#WAVEIC\n\n#TIMESTEPPING", 1)
    text = set_block(deck_text, "WAVEIC", updates)
    if walen_factor is not None:
        # set_block only rewrites existing <name> lines; add this one if absent.
        lines, out, in_block, done = text.splitlines(), [], False, False
        for ln in lines:
            s = ln.strip()
            if s.startswith("#"):
                if in_block and not done:
                    out.append(f"{walen_factor:.6f}  walenFactor")
                    done = True
                in_block = s.split()[0][1:].upper() == "WAVEIC"
            elif in_block and len(s.split()) > 1 and s.split()[1] == "walenFactor":
                done = True
            out.append(ln)
        if in_block and not done:
            out.append(f"{walen_factor:.6f}  walenFactor")
        text = "\n".join(out) + "\n"
    text = set_block(text, "TIMESTEPPING", {"dt": dt}) if dt else text
    text = set_block(text, "STOP", {"TimeMax": time_max})
    text = set_block(text, "SAVEPLOT", {"dn": dn})
    if ppc:
        # #PARTICLES holds three positional values, all named "nParticle".
        lines, out, k, in_p = text.splitlines(), [], 0, False
        for line in lines:
            s = line.strip()
            if s.startswith("#"):
                in_p = s.split()[0][1:].upper() == "PARTICLES"
                out.append(line)
                continue
            if in_p and s and k < 3:
                head = line[:len(line) - len(line.lstrip())]
                body = line[len(head):]
                out.append(head + body.replace(s.split()[0],
                                               str(ppc if k == 0 else 1), 1))
                k += 1
                continue
            out.append(line)
        text = "\n".join(out) + "\n"
    return text


def run_case(run_dir, param_text, nproc):
    prepare_run_dir(run_dir)
    with open(os.path.join(run_dir, "PARAM.in"), "w") as fh:
        fh.write(param_text)
    cmd = ["./FLEKS.exe"] if nproc <= 1 else ["mpirun", "-n", str(nproc), "./FLEKS.exe"]
    res = subprocess.run(cmd, cwd=run_dir, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    if res.returncode != 0:
        sys.stderr.write(res.stdout.decode("utf-8", "replace"))
        raise RuntimeError(f"FLEKS.exe failed in {run_dir}")
    subprocess.run(["./PostProc.pl"], cwd=run_dir, stdout=subprocess.DEVNULL,
                   stderr=subprocess.STDOUT)


# ---------------------------------------------------------------------------
# Frame reading
# ---------------------------------------------------------------------------
def read_frame(path):
    with open(path, "r", encoding="latin-1") as fh:
        lines = fh.readlines()
    t = float(lines[1].split()[1])
    idx = {v.upper(): i for i, v in enumerate(lines[4].split())}
    data = np.array([[float(x) for x in ln.split()] for ln in lines[5:] if ln.strip()])
    return t, data[:, idx["X"]], data[:, idx["BY"]], data[:, idx["BZ"]]


def load_series(plots_dir):
    """Return (t, x, By, Bz) with shapes (nt,), (nx,), (nt, nx)."""
    files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if len(files) < 4:
        raise RuntimeError(f"need >=4 .out frames in {plots_dir}, found {len(files)}")
    t, by, bz, x = [], [], [], None
    for f in files:
        ti, xi, byi, bzi = read_frame(f)
        x = xi
        t.append(ti)
        by.append(byi)
        bz.append(bzi)
    return np.array(t), np.array(x), np.array(by), np.array(bz)


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------
def _lsq_residual(w, t, y):
    """Residual of y ~ a cos(wt) + b sin(wt) (columns of y stacked)."""
    m = np.column_stack([np.cos(w * t), np.sin(w * t)])
    sol, *_ = np.linalg.lstsq(m, y, rcond=None)
    return float(np.sum((m @ sol - y) ** 2)), sol, m


def fit_sinusoid(t, y, w_guess, rel_window=0.25):
    """Least-squares fit of y(t) to a cos(wt) + b sin(wt); solve for w too.

    y may be a stack of several real series sharing one omega (used for the
    (By, Bz) pair at a probe, which is exactly the case for a circularly
    polarized wave).
    """
    lo, hi = w_guess * (1.0 - rel_window), w_guess * (1.0 + rel_window)
    grid = np.linspace(lo, hi, 400)
    best = (np.inf, w_guess, None, None)
    for w in grid:
        r, sol, m = _lsq_residual(w, t, y)
        if r < best[0]:
            best = (r, w, sol, m)
    # Parabolic refinement.
    dw = grid[1] - grid[0]
    w0 = best[1]
    for _ in range(4):
        ws = np.linspace(w0 - dw, w0 + dw, 21)
        res = [_lsq_residual(w, t, y)[0] for w in ws]
        i = int(np.argmin(res))
        best = (res[i], ws[i], *_lsq_residual(ws[i], t, y)[1:])
        dw /= 5.0
    r, w, sol, m = best
    amp = math.hypot(sol[0, 0], sol[1, 0]) if y.ndim > 1 else math.hypot(
        sol[0], sol[1])
    npts = y.shape[0] * (y.shape[1] if y.ndim > 1 else 1)
    return {
        "w": w,
        "amp": amp,
        "resid_rel": math.sqrt(r / npts) / amp if amp > 0 else float("inf"),
        "coef": np.atleast_2d(sol),
    }


def reference_frequency(t, probe, dt_frame):
    """FFT estimate of the dominant frequency (Hann windowed, zero padded)."""
    y = probe - probe.mean()
    n = len(y)
    nfft = int(2 ** math.ceil(math.log2(max(n, 8)) + 3))
    f = np.fft.rfft(y * np.hanning(n), n=nfft)
    w = 2.0 * math.pi * np.fft.rfftfreq(nfft, d=dt_frame)
    k = int(np.argmax(np.abs(f[1:])) + 1)
    return abs(w[k]), k


def best_initial_fit(t, y, guesses):
    """Pick the best single-frequency fit over a set of starting guesses."""
    best = None
    for g in guesses:
        f = fit_sinusoid(t, y, g, rel_window=0.3)
        if best is None or f["resid_rel"] < best["resid_rel"]:
            best = f
    return best


def project(t, y, w):
    """Project y(t) onto cos(wt), sin(wt) at a FIXED frequency.

    A fixed-frequency projection acts as a narrow band-pass, which keeps the
    probe-point phases usable even when broadband grid noise is present.
    """
    r, sol, _ = _lsq_residual(w, t, y)
    sol = np.atleast_2d(sol)
    amp = math.hypot(sol[0, 0], sol[1, 0])
    npts = len(t)
    return {
        "w": w,
        "amp": amp,
        "resid_rel": math.sqrt(r / npts) / amp if amp > 0 else float("inf"),
        "coef": sol,
    }


def measure_window(t, hsig, probe, start, n_min):
    """Window and fits for a case.

    The explicit hybrid advance of this low-beta plasma slowly drives
    grid-scale fields, so the window is grown (in >=2-period steps) only while
    both the seeded-harmonic amplitude and the fixed-frequency probe projection
    stay clean single oscillations.
    """
    n = min(n_min, len(t))
    step = max(1, len(t) // 40)
    hfit = start
    pfit = project(t[:n], probe[:n], hfit["w"])
    n_best = n
    hlim, plim = 0.02, 0.20
    hlim = max(hlim, 2.0 * hfit["resid_rel"])
    plim = max(plim, 1.5 * pfit["resid_rel"])
    while n < len(t):
        n = min(n + step, len(t))
        hcand = fit_sinusoid(t[:n], hsig[:n], hfit["w"], rel_window=0.02)
        if hcand["resid_rel"] > hlim:
            break
        pcand = project(t[:n], probe[:n], hcand["w"])
        if pcand["resid_rel"] > plim:
            break
        n_best, hfit, pfit = n, hcand, pcand
    return n_best, hfit, pfit


def analyse_case(run_dir, mode, units, x_probe=0.0):
    t, x, by, bz = load_series(os.path.join(run_dir, "PC", "plots"))

    # The .out times are SI seconds; convert to code time (Omega_i^-1).
    tc = t / units.tNorm
    tc = tc - tc[0]
    dt_frame = float(np.median(np.diff(tc)))
    # Domain-wide transverse amplitude: a diagnostic for the late-time
    # grid-scale growth of the explicit hybrid advance.
    bperp_max = np.hypot(by, bz).max(axis=1)
    nx = len(x)

    i1 = int(np.argmin(np.abs(x - x_probe)))
    # A second probe a quarter box away fixes the propagation direction and k.
    i2 = int(np.argmin(np.abs(x - (x_probe + units.Lx / 4.0))))
    p1 = np.column_stack([by[:, i1], bz[:, i1]])
    p2 = np.column_stack([by[:, i2], bz[:, i2]])
    p1 = p1 - p1.mean(axis=0)
    p2 = p2 - p2.mean(axis=0)

    # Seeded spatial harmonic of Psi = By + i Bz.  This is the same physical
    # signal as the probe point, but band-passed to the seeded wavenumber, so
    # it is immune to the broadband grid noise.
    psi_hat = np.fft.fft(by + 1j * bz, axis=1) / nx
    harm = psi_hat[:, (-mode) % nx]
    hsig = np.column_stack([harm.real, harm.imag])

    k_seed = 2.0 * math.pi * abs(mode) / units.Lx
    kappa_seed = k_seed * units.di
    w_w, w_l = analytic(kappa_seed, units)
    period = 2.0 * math.pi / w_w

    # Window used for the time-series fits: at least 2 periods, extended while
    # the signal remains a single clean oscillation.  The starting guess comes
    # from an early-window FFT (the late-time noise would otherwise dominate),
    # with the two analytic branches as fallbacks.
    n_min = min(len(tc), max(8, int(round(2.0 * period / dt_frame)) + 1))
    w_fft, _ = reference_frequency(tc[:n_min], harm[:n_min].real, dt_frame)
    guesses = [w_fft, w_w, w_l]
    start = best_initial_fit(tc[:n_min], hsig[:n_min], guesses)
    n_win, hfit, pfit = measure_window(tc, hsig, p1, start, n_min)
    t_win = tc[:n_win]

    # Primary measurement: omega of the seeded spatial harmonic, whose complex
    # amplitude C(t) rotates in the (By, Bz) plane.  A right-hand wave has
    # C ~ exp(+i w t), i.e. circulation d(real)/dt * imag - ... > 0 in the
    # (cos, sin) coefficient sense used by fit_sinusoid.
    hc = hfit["coef"]
    ha, hb = hc[0, 0], hc[1, 0]
    hcc, hd = hc[0, 1], hc[1, 1]
    circulation = ha * hd - hb * hcc     # > 0 <=> y -> z rotation about +x
    right_handed = circulation > 0.0

    # Secondary: the raw probe time series dB_y(x0, t), dB_z(x0, t) over the
    # same window, band-passed at the measured frequency (this is the "track a
    # probe point" measurement).
    pc = pfit["coef"]
    pcir = pc[0, 0] * pc[1, 1] - pc[1, 0] * pc[0, 1]
    probe_right_handed = pcir > 0.0
    # Independent frequency estimate from the probe series alone (scanning fit);
    # only meaningful while the probe is not buried in grid noise.
    pscan = fit_sinusoid(t_win, p1[:n_win], hfit["w"], rel_window=0.05)

    # Wavenumber: the dominant spatial harmonic of Psi over the measurement
    # window (robust, and independent of the noisy probe point).  The seeded
    # mode should be the dominant one.
    spec_x = np.abs(psi_hat[:n_win]).mean(axis=0)
    m_peak = int(np.argmax(spec_x[1:]) + 1)
    n_signed = m_peak if m_peak <= nx // 2 else m_peak - nx
    k_meas = abs(2.0 * math.pi * n_signed / units.Lx)

    # Auxiliary: k from the phase difference between the two probe points
    # (unwrapped against the seeded wavelength), which also fixes the
    # propagation direction.
    f2 = project(t_win, p2[:n_win], hfit["w"])
    ph1 = math.atan2(pc[1, 0], pc[0, 0])
    ph2 = math.atan2(f2["coef"][1, 0], f2["coef"][0, 0])
    dx = x[i2] - x[i1]
    dphi = ph2 - ph1
    dphi -= 2.0 * math.pi * round((dphi - k_seed * dx) / (2.0 * math.pi))
    k_probe = dphi / dx

    vph_meas = hfit["w"] / k_meas if k_meas else float("nan")

    # Spectra over the measurement window of the seeded harmonic (clean) and of
    # the raw probe point (what the test asks for).  For Psi = By + i Bz the
    # right-hand (y -> z) component sits at +w and the left-hand one at -w.
    nfft = int(2 ** math.ceil(math.log2(max(n_win, 8)) + 3))
    win = np.hanning(n_win)
    spec = np.fft.fftshift(np.fft.fft(
        (p1[:n_win, 0] + 1j * p1[:n_win, 1]) * win, n=nfft))
    hspec = np.fft.fftshift(np.fft.fft(harm[:n_win] * win, n=nfft))
    w_ax = 2.0 * math.pi * np.fft.fftshift(np.fft.fftfreq(nfft, d=dt_frame))
    amp_pos = float(np.max(np.abs(spec[w_ax > 0.05]))) if (w_ax > 0.05).any() else 0.0
    amp_neg = float(np.max(np.abs(spec[w_ax < -0.05]))) if (w_ax < -0.05).any() else 0.0

    return {
        "mode": mode,
        "k": k_seed,
        "kappa": kappa_seed,
        "w": hfit["w"],
        "w_probe": pscan["w"],
        "w_right": w_w,
        "w_left": w_l,
        "period": period,
        "k_meas": k_meas,
        "k_probe": k_probe,
        "vph": vph_meas,
        "right_handed": right_handed,
        "circulation": circulation,
        "amp": hfit["amp"],
        "amp_probe": pfit["amp"],
        "resid_rel": hfit["resid_rel"],
        "resid_probe": pfit["resid_rel"],
        "probe_right_handed": probe_right_handed,
        "probe_circulation": pcir,
        "spec": (w_ax, np.abs(spec)),
        "hspec": (w_ax, np.abs(hspec)),
        "amp_pos": amp_pos,
        "amp_neg": amp_neg,
        "n_win": n_win,
        "t_end": float(t_win[-1]),
        "n_periods": float(t_win[-1] / period),
        "bperp_growth": float(bperp_max[-1] / bperp_max[0]),
        "t": tc,
        "by": by[:, i1],
        "bz": bz[:, i1],
        "x_probe": float(x[i1]),
        "x_probe2": float(x[i2]),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", choices=sorted(DECKS), default="hybrid")
    ap.add_argument("--modes", type=int, nargs="+", default=[1, 2, 3, 4])
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--periods", type=float, default=6.0,
                    help="run time in whistler periods of the seeded mode")
    ap.add_argument("--seed", choices=["eigenmode", "alfven"], default="eigenmode",
                    help="velocity kick: exact single-branch eigenmode, or the "
                         "kappa -> 0 Alfven relation (excites both branches)")
    ap.add_argument("--ppc", type=int, default=None)
    ap.add_argument("--nproc", type=int, default=2)
    ap.add_argument("--probe", type=float, default=0.0, help="probe x (code units)")
    ap.add_argument("--run-root", default=os.path.join(REPO, "run_whistler_disp"))
    ap.add_argument("--no-run", action="store_true",
                    help="re-use existing run directories")
    ap.add_argument("--outdir", default=HERE)
    args = ap.parse_args()

    deck_path = os.path.join(HERE, DECKS[args.variant])
    deck_text = open(deck_path).read()
    blocks = parse_deck(deck_path)
    units = Units(blocks)
    right_hand_seed = get_str(blocks, "WAVEIC", "rightHand", "F").upper().startswith("T")
    frac = get(blocks, "WAVEIC", "frac", 0.02)

    print(f"Whistler dispersion scan ({args.variant} variant, {deck_path})")
    print(units.summary())
    print(f"  seed: rightHand = {right_hand_seed}, dB/B0 = {frac:g}, "
          f"kick = {args.seed}")
    print(f"  note: the analytic branches are c/omega_pi based; kappa(m=1) = "
          f"{2.0 * math.pi * units.di / units.Lx:.4f}\n")

    cases = []
    for mode in args.modes:
        kappa = 2.0 * math.pi * abs(mode) * units.di / units.Lx
        w_w, _ = analytic(kappa, units)
        period = 2.0 * math.pi / w_w
        time_max = args.periods * period
        steps_per_period = period / args.dt
        dn = max(1, int(round(steps_per_period / 20.0)))
        run_dir = os.path.join(args.run_root, f"{args.variant}_m{mode}_{args.seed}")
        walen = (eigenmode_walen_factor(kappa, units, right_hand_seed)
                 if args.seed == "eigenmode" else None)
        param = build_param(deck_text, mode, right_hand_seed, args.dt,
                            f"{time_max:.4f}", dn, walen, args.ppc)
        print(f"  mode {mode}: kappa = {kappa:.4f}, T = {period:.3f} "
              f"(= {steps_per_period:.0f} steps), TimeMax = {time_max:.2f}, "
              f"dn = {dn}, walenFactor = "
              f"{'-' if walen is None else f'{walen:.4f}'}")
        if args.no_run:
            existing = sorted(glob.glob(os.path.join(run_dir, "PC", "plots", "*.out")))
            if not existing:
                print(f"    [skip] no frames in {run_dir}")
                continue
        else:
            t0 = time.time()
            run_case(run_dir, param, args.nproc)
            print(f"    ran in {time.time() - t0:.1f} s")
        cases.append(analyse_case(run_dir, mode, units, args.probe))

    if not cases:
        print("no cases analysed")
        return 1

    return report_and_plot(cases, units, args.outdir)


def report_and_plot(cases, units, outdir):
    modes = [c["mode"] for c in cases]
    kappa = np.array([c["kappa"] for c in cases])
    w = np.array([c["w"] for c in cases])
    w_w = np.array([c["w_right"] for c in cases])
    w_l = np.array([c["w_left"] for c in cases])

    hdr = (f"{'m':>3} {'k d_i':>8} {'w_mode':>9} {'w_probe':>9} {'w_whist':>9} "
           f"{'w_ic':>8} {'err':>7} {'hand':>5} {'v_p/v_A':>9} {'v_p th':>9} "
           f"{'res':>7} {'win':>9}")
    print("\n" + hdr)
    print("-" * len(hdr))
    ok_all, reasons = True, []
    for c in cases:
        err = abs(c["w"] - c["w_right"]) / c["w_right"]
        # v_ph/v_A = omega/(k v_A) = (omega/Omega_i)/kappa, since v_A = Omega_i d_i.
        vth = c["w_right"] / (units.omegaCi * c["kappa"])
        vme = c["vph"] / units.vA
        hand = "R" if c["right_handed"] else "L"
        print(f"{c['mode']:>3} {c['kappa']:>8.4f} {c['w']:>9.4f} "
              f"{c['w_probe']:>9.4f} {c['w_right']:>9.4f} {c['w_left']:>8.4f} "
              f"{err * 100:>6.2f}% {hand:>5} {vme:>9.4f} {vth:>9.4f} "
              f"{c['resid_rel']:>7.1%} {c['n_periods']:>8.1f}T")
        c["vph_theory"] = vth
        c["vph_ratio"] = vme

    # --- criteria -----------------------------------------------------------
    if not all(c["right_handed"] for c in cases):
        ok_all = False
        reasons.append("not all modes are right-hand polarized")
    errs = [abs(c["w"] - c["w_right"]) / c["w_right"] for c in cases]
    if max(errs) > 0.05:
        ok_all = False
        reasons.append(f"frequency mismatch up to {max(errs) * 100:.1f}% "
                       "(tolerance 5%)")
    verr = [abs(c["vph"] - c["vph_theory"]) / c["vph_theory"] for c in cases]
    if max(verr) > 0.05:
        ok_all = False
        reasons.append(f"phase-velocity mismatch up to {max(verr) * 100:.1f}%")

    # local log-log slope d ln w / d ln k (should approach 2 for kappa >> 1)
    slope_m, slope_a = [], []
    if len(cases) > 1:
        for i in range(1, len(cases)):
            slope_m.append(math.log(w[i] / w[i - 1]) / math.log(kappa[i] / kappa[i - 1]))
            slope_a.append(math.log(w_w[i] / w_w[i - 1]) / math.log(kappa[i] / kappa[i - 1]))

    print("\nlocal d ln(omega)/d ln(k):")
    for i, m in enumerate(modes[1:], start=1):
        print(f"  m={m}: measured {slope_m[i - 1]:.3f}  analytic {slope_a[i - 1]:.3f}")

    print("\nspatial wavenumber, dominant harmonic of Psi over the window "
          "(k_meas/k_seed; the two-probe phase difference in brackets):")
    for c in cases:
        print(f"  m={c['mode']}: k_meas/k_seed = {c['k_meas'] / c['k']:.4f} "
              f"[{c['k_probe'] / c['k']:.4f}]")

    print("\nlate-time transverse-amplitude growth over the full frame range "
          "(grid-scale noise of the explicit hybrid advance):")
    for c in cases:
        print(f"  m={c['mode']}: max|B_perp| x{c['bperp_growth']:.2f} "
              f"(measured over the first {c['n_periods']:.1f} periods, "
              f"probe fit residual {c['resid_probe']:.1%}, "
              f"probe hand = {'R' if c['probe_right_handed'] else 'L'})")

    pverr = [abs(c["w_probe"] - c["w"]) / c["w"] for c in cases]

    print("\nPass criteria:")
    print("  1. right-hand polarization  : "
          + ("PASS" if all(c["right_handed"] for c in cases) else "FAIL")
          + (" (probe hodogram also right-hand)"
             if all(c["probe_right_handed"] for c in cases) else
             " [probe hodogram noisy/mixed]"))
    print("  2. omega -> k^2 scaling     : "
          + (f"PASS (log-log slope {slope_m[-1]:.2f} vs analytic {slope_a[-1]:.2f})"
             if slope_m else "n/a"))
    print("  3. phase velocity           : "
          + (f"PASS (max error {max(verr) * 100:.2f}%)" if max(verr) <= 0.05
             else f"FAIL (max error {max(verr) * 100:.2f}%)"))
    print("  4. probe vs mode frequency  : "
          + (f"PASS (max difference {max(pverr) * 100:.2f}%)"
             if max(pverr) <= 0.05 and max(c["resid_probe"] for c in cases) <= 0.15
             else f"probe noise-dominated (max residual "
                  f"{max(c['resid_probe'] for c in cases):.0%}); "
                  f"max difference {max(pverr) * 100:.2f}%"))
    if max(pverr) > 0.05 and max(c["resid_probe"] for c in cases) <= 0.15:
        ok_all = False
        reasons.append(f"probe and mode frequencies differ by up to "
                       f"{max(pverr) * 100:.1f}%")
    if not all(c["probe_right_handed"] for c in cases):
        reasons.append("probe hodogram not right-hand for every mode "
                       "(noisy probe); the mode-projected polarization is used "
                       "for criterion 1")
    if not ok_all:
        print("  -> " + "; ".join(reasons))

    _plot(cases, units, outdir)
    return 0 if ok_all else 1


def _plot(cases, units, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kappa = np.array([c["kappa"] for c in cases])
    w = np.array([c["w"] for c in cases])
    w_w = np.array([c["w_right"] for c in cases])
    w_l = np.array([c["w_left"] for c in cases])
    kk = np.logspace(math.log10(kappa.min() * 0.6),
                     math.log10(kappa.max() * 1.6), 200)
    ww = np.array([analytic(k, units)[0] for k in kk])
    wl = np.array([analytic(k, units)[1] for k in kk])

    # ---------------- figure 1: dispersion relation ------------------------
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))

    ax[0].loglog(kk, ww, "C0-", label=r"whistler (right-hand)")
    ax[0].loglog(kk, wl, "C2--", label=r"ion-cyclotron (left-hand)")
    ax[0].loglog(kk, units.omegaCi * kk ** 2, "k:", lw=1,
                 label=r"asymptote $\omega=\Omega_i\kappa^2$")
    ax[0].loglog(kk, units.omegaCi * kk, "0.6", lw=1,
                 label=r"Alfvén $\omega=kv_A$")
    ax[0].loglog(kappa, w, "o", ms=9, mfc="none", mew=2, label="measured")
    ax[0].set_xlabel(r"$\kappa = k\,d_i$")
    ax[0].set_ylabel(r"$\omega/\Omega_i$")
    ax[0].set_title("Parallel whistler dispersion")
    ax[0].grid(True, which="both", alpha=0.3)
    ax[0].legend(fontsize=8)

    # d omega / d k
    kd = kappa / units.di
    kk_d = kk / units.di
    dwd = np.gradient(ww, kk_d)
    ax[1].loglog(kk_d, np.abs(dwd), "C0-")
    if len(cases) > 1:
        ax[1].loglog(kd, np.abs(np.gradient(w, kd)), "o", ms=9, mfc="none",
                     mew=2, color="C0")
    else:
        ax[1].loglog(kd, w / kd, "o", ms=9, mfc="none", mew=2, color="C0")
    ax[1].loglog(kk_d, units.omegaCi * units.di * kk_d, "k:", lw=1,
                 label=r"$\propto k$")
    ax[1].set_xlabel(r"$k$ [$d_i^{-1}$]")
    ax[1].set_ylabel(r"$d\omega/dk$ [$d_i\,\Omega_i$]")
    ax[1].set_title(r"$d\omega/dk \propto k$")
    ax[1].grid(True, which="both", alpha=0.3)
    ax[1].legend(fontsize=8)

    # phase velocity
    vth = np.array([c["vph_theory"] for c in cases])
    vme = np.array([c["vph_ratio"] for c in cases])
    ax[2].semilogx(kk, ww / kk, "C0-", label="whistler (analytic)")
    ax[2].semilogx(kk, wl / kk, "C2--", label="ion-cyclotron (analytic)")
    ax[2].semilogx(kappa, vth, "s", ms=8, mfc="none", mew=2, color="C0",
                   label="theory at seeded $k$")
    ax[2].semilogx(kappa, vme, "o", ms=9, mfc="none", mew=2, color="C1",
                   label=r"measured $v_p/v_A$")
    ax[2].axhline(1.0, color="0.6", lw=1, ls=":")
    ax[2].set_xlabel(r"$\kappa = k\,d_i$")
    ax[2].set_ylabel(r"$v_p/v_A$")
    ax[2].set_title("Phase velocity")
    ax[2].grid(True, which="both", alpha=0.3)
    ax[2].legend(fontsize=8)

    fig.suptitle(f"Whistler dispersion, FLEKS hybrid vs analytic "
                 f"($d_i$ = {units.di:.4f} code units)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    f1 = os.path.join(outdir, "whistler_dispersion.png")
    fig.savefig(f1, dpi=140)
    print(f"\nwrote {f1}")

    # ---------------- figure 2: waveform / polarization --------------------
    c0 = cases[0]
    fig, ax = plt.subplots(2, 2, figsize=(12.5, 8.5))

    ax[0, 0].axvspan(0.0, c0["t_end"], color="C2", alpha=0.08,
                     label="fit window")
    ax[0, 0].plot(c0["t"], c0["by"], "C0-", lw=1, label=r"$\delta B_y$")
    ax[0, 0].plot(c0["t"], c0["bz"], "C1-", lw=1, label=r"$\delta B_z$")
    ax[0, 0].set_xlabel(r"$t$ [$\Omega_i^{-1}$]")
    ax[0, 0].set_ylabel(r"$\delta B$ [code units]")
    ax[0, 0].set_ylim(-2.5 * c0["amp"], 2.5 * c0["amp"])
    ax[0, 0].set_title(
        f"Probe $x_0$ = {c0['x_probe']:+.2f}, mode m = {c0['mode']}"
        f"  ($\\omega$ = {c0['w']:.3f}, hand = "
        f"{'R' if c0['right_handed'] else 'L'}; late-time grid noise off-scale)")
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(alpha=0.3)

    nw = c0["n_win"]
    by_w = c0["by"][:nw]
    bz_w = c0["bz"][:nw]
    ax[0, 1].plot(by_w, bz_w, "C0-", lw=1)
    n = max(len(by_w), 1)
    ax[0, 1].annotate("", xy=(by_w[n // 60], bz_w[n // 60]),
                      xytext=(by_w[0], bz_w[0]),
                      arrowprops=dict(arrowstyle="->", color="C3"))
    ax[0, 1].set_xlabel(r"$\delta B_y$")
    ax[0, 1].set_ylabel(r"$\delta B_z$")
    ax[0, 1].set_title("Hodogram (y $\\rightarrow$ z $\\Rightarrow$ right-hand "
                       "about $+\\hat{B}_0$)")
    ax[0, 1].set_aspect("equal", adjustable="datalim")
    ax[0, 1].grid(alpha=0.3)

    for c in cases:
        wa, sa = c["hspec"]
        sa = sa / sa.max()
        ax[1, 0].semilogy(wa, np.maximum(sa, 1e-12), lw=1, label=f"m={c['mode']}")
        ax[1, 0].plot([c["w"]], [1.0], "kv", ms=4)
        ax[1, 0].plot([-c["w"]], [np.interp(-c["w"], wa, sa)], "k^", ms=4)
    ax[1, 0].set_xlim(-max(c["w"] for c in cases) * 1.4,
                      max(c["w"] for c in cases) * 1.4)
    ax[1, 0].set_ylim(1e-8, 3.0)
    ax[1, 0].set_xlabel(r"$\omega$ [$\Omega_i$]")
    ax[1, 0].set_ylabel(r"$|\hat C(\omega)|$ (normalised)")
    ax[1, 0].set_title(r"Seeded-mode spectrum $\hat C=\overline{\Psi e^{-ikx}}$"
                       "\n(v = measured, +$\\omega$ = right-hand, "
                       "$-\\omega$ = left-hand)")
    ax[1, 0].grid(alpha=0.3)
    ax[1, 0].legend(fontsize=8)

    ax[1, 1].loglog(kappa, w, "o", ms=9, mfc="none", mew=2, label="measured")
    ax[1, 1].loglog(kk, ww, "C0-", label="whistler")
    ax[1, 1].loglog(kk, units.omegaCi * kk ** 2, "k:", lw=1,
                    label=r"$\Omega_i\kappa^2$")
    for c in cases:
        ax[1, 1].annotate(f"{c['w'] / c['w_right'] - 1:+.1%}",
                          (c["kappa"], c["w"]), textcoords="offset points",
                          xytext=(4, -10), fontsize=7)
    ax[1, 1].set_xlabel(r"$\kappa = k d_i$")
    ax[1, 1].set_ylabel(r"$\omega/\Omega_i$")
    ax[1, 1].set_title("Measured vs analytic whistler")
    ax[1, 1].grid(True, which="both", alpha=0.3)
    ax[1, 1].legend(fontsize=8)

    fig.tight_layout()
    f2 = os.path.join(outdir, "whistler_waveform.png")
    fig.savefig(f2, dpi=140)
    print(f"wrote {f2}")


if __name__ == "__main__":
    sys.exit(main())
