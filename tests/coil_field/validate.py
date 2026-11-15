#!/usr/bin/env python3
"""Independent validator for the generic external-coil and probe test suite.

Covers all five variants under ``tests/coil_field/``:
  - helmholtz_toroidal : 3D multi-turn solenoids + regularized 1/R toroidal field
  - segment_loop       : 3D square loop built from 4 Biot-Savart segments
  - filament2d         : fake-2D dual wires with linear_ramp & table waveforms
  - waveform_static    : 3D cosine_ramp + sinusoid loops, curl(E0)=-dB0/dt,
                         plus subprocess tests of tools/extract_probes.py
  - waveform_hybrid    : Hybrid PIC with crowbar_rlc waveform & Faraday response
"""

from __future__ import annotations

import csv
import glob
import logging
import math
import os
import subprocess
import sys

import numpy as np

import tests._shared.run_dir as _run_dir


def _carlson_rf(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Vectorized Carlson symmetric elliptic integral R_F(x, y, z)."""
    xn = np.asarray(x, dtype=np.float64).copy()
    yn = np.asarray(y, dtype=np.float64).copy()
    zn = np.asarray(z, dtype=np.float64).copy()
    errtol = 1.0e-4

    while True:
        a = (xn + yn + zn) / 3.0
        dx = np.abs(a - xn)
        dy = np.abs(a - yn)
        dz = np.abs(a - zn)
        if np.max(np.maximum(np.maximum(dx, dy), dz) / a) <= errtol:
            break
        lam = np.sqrt(xn * yn) + np.sqrt(yn * zn) + np.sqrt(zn * xn)
        xn = 0.25 * (xn + lam)
        yn = 0.25 * (yn + lam)
        zn = 0.25 * (zn + lam)

    X = (a - xn) / a
    Y = (a - yn) / a
    Z = -(X + Y)
    E2 = X * Y - Z * Z
    E3 = X * Y * Z
    return (
        1.0
        - E2 / 10.0
        + E3 / 14.0
        + (E2 * E2) / 24.0
        - (3.0 / 44.0) * E2 * E3
    ) / np.sqrt(a)


def _carlson_rd(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Vectorized Carlson symmetric elliptic integral R_D(x, y, z)."""
    xn = np.asarray(x, dtype=np.float64).copy()
    yn = np.asarray(y, dtype=np.float64).copy()
    zn = np.asarray(z, dtype=np.float64).copy()
    sum_acc = np.zeros_like(xn)
    fac = 1.0
    errtol = 1.0e-4

    while True:
        a = (xn + yn + 3.0 * zn) / 5.0
        dx = np.abs(a - xn)
        dy = np.abs(a - yn)
        dz = np.abs(a - zn)
        if np.max(np.maximum(np.maximum(dx, dy), dz) / a) <= errtol:
            break
        lam = np.sqrt(xn * yn) + np.sqrt(yn * zn) + np.sqrt(zn * xn)
        sum_acc += fac / (np.sqrt(zn) * (zn + lam))
        fac *= 0.25
        xn = 0.25 * (xn + lam)
        yn = 0.25 * (yn + lam)
        zn = 0.25 * (zn + lam)

    X = (a - xn) / a
    Y = (a - yn) / a
    Z = -(X + Y) / 3.0
    E2 = X * Y - 6.0 * Z * Z
    E3 = (3.0 * X * Y - 8.0 * Z * Z) * Z
    E4 = 3.0 * (X * Y - Z * Z) * Z * Z
    E5 = X * Y * Z * Z * Z
    poly = (
        1.0
        - (3.0 / 14.0) * E2
        + E3 / 6.0
        + (9.0 / 88.0) * E2 * E2
        - (3.0 / 22.0) * E4
        - (9.0 / 52.0) * E2 * E3
        + (3.0 / 26.0) * E5
    )
    return 3.0 * sum_acc + fac * poly / (a * np.sqrt(a))


def _ellipke(m: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    """Pure-NumPy complete elliptic integrals K(m), E(m) for 0 <= m < 1."""
    m_arr = np.asarray(m, dtype=np.float64)
    m_safe = np.clip(m_arr, 0.0, 1.0 - 1.0e-15)
    zeros = np.zeros_like(m_safe)
    ones = np.ones_like(m_safe)
    one_minus_m = 1.0 - m_safe
    K = _carlson_rf(zeros, one_minus_m, ones)
    E = K - (m_safe / 3.0) * _carlson_rd(zeros, one_minus_m, ones)
    singular = m_arr >= 1.0
    if np.any(singular):
        K = np.where(singular, np.inf, K)
        E = np.where(singular, 1.0, E)
    return K, E


def ellipk(m: np.ndarray | float) -> np.ndarray:
    """Complete elliptic integral of the first kind K(m)."""
    return _ellipke(m)[0]


def ellipe(m: np.ndarray | float) -> np.ndarray:
    """Complete elliptic integral of the second kind E(m)."""
    return _ellipke(m)[1]

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from tools.extract_probes import read_frame  # noqa: E402

MU0_OVER_4PI = 1.0e-7


def set_run_dir(run_dir):
    """Synchronize module RUN_DIR with tests/validate_tests.py."""
    global RUN_DIR
    RUN_DIR = run_dir
    _run_dir.set_run_dir(run_dir)


# ---------------------------------------------------------------------------
# Independent SI analytical evaluators (using SciPy elliptic integrals)
# ---------------------------------------------------------------------------
def analytic_loop_ba(
    coords: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
    radius: float,
    current: float,
    core_radius: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute SI (B_ref, A_ref) of a unit-scale circular loop at (N, 3) points."""
    axis_u = np.asarray(axis, dtype=np.float64)
    axis_u = axis_u / np.linalg.norm(axis_u)
    center = np.asarray(center, dtype=np.float64)
    k_i = MU0_OVER_4PI * float(current)

    d = coords - center[np.newaxis, :]
    zeta = d @ axis_u
    r_perp = d - zeta[:, np.newaxis] * axis_u[np.newaxis, :]
    rho2 = np.sum(r_perp * r_perp, axis=1)
    rho = np.sqrt(rho2)

    r2 = radius * radius
    z_eff2 = zeta * zeta + core_radius * core_radius
    eps2 = 0.0 if core_radius > 0.0 else (1.0e-24 * r2)

    alpha2 = (radius - rho) ** 2 + z_eff2 + eps2
    beta2 = (radius + rho) ** 2 + z_eff2 + eps2
    beta = np.sqrt(beta2)

    m = np.minimum((4.0 * radius * rho) / beta2, 1.0 - 1.0e-15)
    K, E = _ellipke(m)

    b_zeta = (2.0 * k_i / beta) * (K + ((r2 - rho2 - z_eff2) / alpha2) * E)

    # Smooth on-axis handling (m -> 0 gives exact finite B_rho/rho and A_phi/rho).
    small = m < 1.0e-5
    b_rho_over_rho = np.empty_like(m)
    a_phi_over_rho = np.empty_like(m)

    if np.any(~small):
        ms = m[~small]
        b_rho_over_rho[~small] = (
            2.0 * k_i * zeta[~small] / (rho2[~small] * beta[~small])
        ) * (-K[~small] + ((r2 + rho2[~small] + z_eff2[~small]) / alpha2[~small]) * E[~small])
        a_phi_over_rho[~small] = (k_i * beta[~small] / rho2[~small]) * (
            (2.0 - ms) * K[~small] - 2.0 * E[~small]
        )

    if np.any(small):
        ms = m[small]
        b5 = beta2[small] ** 2 * beta[small]
        b3 = beta2[small] * beta[small]
        b_rho_over_rho[small] = (
            8.0 * math.pi * k_i * r2 * zeta[small] / b5
        ) * (0.375 + ms * (0.46875 + ms * (525.0 / 1024.0 + ms * (2205.0 / 4096.0))))
        a_phi_over_rho[small] = (8.0 * math.pi * k_i * r2 / b3) * (
            0.125 + ms * (0.09375 + ms * (75.0 / 1024.0 + ms * (245.0 / 4096.0)))
        )

    b = b_zeta[:, np.newaxis] * axis_u[np.newaxis, :] + b_rho_over_rho[:, np.newaxis] * r_perp
    phi_vec = np.cross(axis_u[np.newaxis, :], r_perp)
    a = a_phi_over_rho[:, np.newaxis] * phi_vec
    return b, a


def analytic_solenoid_ba(
    coords: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
    radius: float,
    length: float,
    n_turns: int,
    current: float,
    core_radius: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute SI (B_ref, A_ref) of a stack of *n_turns* coaxial loops."""
    axis_u = np.asarray(axis, dtype=np.float64)
    axis_u = axis_u / np.linalg.norm(axis_u)
    center = np.asarray(center, dtype=np.float64)
    b_tot = np.zeros_like(coords)
    a_tot = np.zeros_like(coords)
    for k in range(n_turns):
        offset = (
            0.0
            if n_turns == 1
            else (-0.5 * length + (k + 0.5) * (length / n_turns))
        )
        b_k, a_k = analytic_loop_ba(
            coords, center + offset * axis_u, axis_u, radius, current, core_radius
        )
        b_tot += b_k
        a_tot += a_k
    return b_tot, a_tot


def analytic_segment_b(
    coords: np.ndarray,
    x1: np.ndarray,
    x2: np.ndarray,
    current: float,
    core_radius: float = 0.0,
) -> np.ndarray:
    """Compute regularized Biot-Savart B [T] of a finite wire segment x1 -> x2."""
    x1 = np.asarray(x1, dtype=np.float64)
    x2 = np.asarray(x2, dtype=np.float64)
    delta = x2 - x1
    seg_len = float(np.linalg.norm(delta))
    d_hat = delta / seg_len
    rx1 = coords - x1[np.newaxis, :]
    s1 = rx1 @ d_hat
    s2 = s1 - seg_len
    perp = rx1 - s1[:, np.newaxis] * d_hat[np.newaxis, :]
    rho_eff2 = (
        np.sum(perp * perp, axis=1)
        + core_radius * core_radius
        + 1.0e-24 * seg_len * seg_len
    )
    r1 = np.sqrt(s1 * s1 + rho_eff2)
    r2 = np.sqrt(s2 * s2 + rho_eff2)
    fac = (MU0_OVER_4PI * current / rho_eff2) * (s1 / r1 - s2 / r2)
    return fac[:, np.newaxis] * np.cross(d_hat[np.newaxis, :], perp)


def analytic_toroidal_b(
    coords: np.ndarray,
    center: np.ndarray,
    axis: np.ndarray,
    b_ref: float,
    r_ref: float,
    r_core: float,
) -> np.ndarray:
    """Compute regularized 1/R toroidal field [T]."""
    axis_u = np.asarray(axis, dtype=np.float64)
    axis_u = axis_u / np.linalg.norm(axis_u)
    d = coords - np.asarray(center, dtype=np.float64)[np.newaxis, :]
    zeta = d @ axis_u
    perp = d - zeta[:, np.newaxis] * axis_u[np.newaxis, :]
    rho_eff2 = np.sum(perp * perp, axis=1) + r_core * r_core + 1.0e-24 * r_ref * r_ref
    fac = (b_ref * r_ref) / rho_eff2
    return fac[:, np.newaxis] * np.cross(axis_u[np.newaxis, :], perp)


def analytic_filament2d_ba(
    xy: np.ndarray, x0: float, y0: float, current: float, core_radius: float
) -> tuple[np.ndarray, np.ndarray]:
    """Compute SI (B_ref, A_ref) of an infinite wire along +z at (x0, y0)."""
    dx = xy[:, 0] - x0
    dy = xy[:, 1] - y0
    r_eff2 = dx * dx + dy * dy + core_radius * core_radius + 1.0e-24
    k_i = MU0_OVER_4PI * current
    b = np.zeros((xy.shape[0], 3), dtype=np.float64)
    a = np.zeros((xy.shape[0], 3), dtype=np.float64)
    fac_b = 2.0 * k_i / r_eff2
    b[:, 0] = -fac_b * dy
    b[:, 1] = fac_b * dx
    a[:, 2] = -k_i * np.log(r_eff2)
    return b, a


# ---------------------------------------------------------------------------
# Waveform evaluators
# ---------------------------------------------------------------------------
def wf_linear_ramp(t: float, t0: float, t1: float, f0: float, f1: float) -> tuple[float, float]:
    if t < t0:
        return f0, 0.0
    if t >= t1:
        return f1, 0.0
    slope = (f1 - f0) / (t1 - t0)
    return f0 + slope * (t - t0), slope


def wf_cosine_ramp(t: float, t0: float, t1: float, f0: float, f1: float) -> tuple[float, float]:
    if t < t0:
        return f0, 0.0
    if t >= t1:
        return f1, 0.0
    u = (t - t0) / (t1 - t0)
    df = f1 - f0
    return f0 + 0.5 * df * (1.0 - math.cos(math.pi * u)), (0.5 * math.pi / (t1 - t0)) * df * math.sin(math.pi * u)


def wf_sinusoid(
    t: float, t0: float, f_off: float, f_amp: float, freq_hz: float, phase_deg: float = 0.0
) -> tuple[float, float]:
    ph0 = math.radians(phase_deg)
    if t < t0:
        return f_off + f_amp * math.sin(ph0), 0.0
    omega = 2.0 * math.pi * freq_hz
    ang = omega * (t - t0) + ph0
    return f_off + f_amp * math.sin(ang), omega * f_amp * math.cos(ang)


def wf_crowbar_rlc(
    t: float, t0: float, t_rise: float, tau_decay: float, do_crowbar: bool, tau_crowbar: float
) -> tuple[float, float]:
    tau = t - t0
    if tau <= 0.0:
        return 0.0, 0.0
    omega = 0.5 * math.pi / t_rise
    gamma = (1.0 / tau_decay) if tau_decay > 0.0 else 0.0
    t_peak = (math.atan(omega / gamma) / omega) if gamma > 0.0 else t_rise
    a_norm = math.exp(gamma * t_peak) / math.sin(omega * t_peak)
    if not do_crowbar or tau <= t_peak:
        env = a_norm * math.exp(-gamma * tau)
        return env * math.sin(omega * tau), env * (omega * math.cos(omega * tau) - gamma * math.sin(omega * tau))
    if tau_crowbar > 0.0:
        val = math.exp(-(tau - t_peak) / tau_crowbar)
        return val, -val / tau_crowbar
    return 1.0, 0.0


# ---------------------------------------------------------------------------
# Frame helpers
# ---------------------------------------------------------------------------
def _load_frames():
    files = sorted(glob.glob(os.path.join(RUN_DIR, "PC/plots", "*.out")))
    return [read_frame(p) for p in files]


def _vec(frame: dict, prefix: str) -> np.ndarray:
    return np.column_stack(
        [frame["fields"][prefix + axis] for axis in ("x", "y", "z")]
    )


def _coords3(frame: dict) -> np.ndarray:
    coords = frame["coords"]
    if coords.shape[1] == 3:
        return coords
    out = np.zeros((coords.shape[0], 3), dtype=np.float64)
    out[:, : coords.shape[1]] = coords
    return out


def _rel_max_err(got: np.ndarray, ref: np.ndarray) -> float:
    ref_scale = max(float(np.max(np.linalg.norm(ref, axis=-1))), 1.0e-12)
    return float(np.max(np.linalg.norm(got - ref, axis=-1)) / ref_scale)


# ---------------------------------------------------------------------------
# Log validator
# ---------------------------------------------------------------------------
def validate_log(pic_diags=None, test_name=None):
    """Check finite PIC diagnostics and static Eb conservation where applicable."""
    if not pic_diags or len(pic_diags) < 2:
        return False, f"Insufficient PIC log rows for {test_name}"

    for entry in pic_diags:
        for key in ("Etot", "Ee", "Eb", "Epart"):
            if not math.isfinite(entry.get(key, 0.0)):
                return False, f"Non-finite {key} in PIC log"

    if test_name in ("coil_field_helmholtz_toroidal", "coil_field_segment_loop"):
        eb0, eb1 = pic_diags[0]["Eb"], pic_diags[-1]["Eb"]
        if abs(eb0) > 0 and abs(eb1 - eb0) > 1.0e-10 * abs(eb0):
            return False, f"Static Eb drifted ({eb0:.6e} -> {eb1:.6e})"
        return True, "Passed (finite logs, static Eb constant)"

    if pic_diags[-1]["Eb"] <= 0.0 and pic_diags[0]["Eb"] <= 0.0:
        return False, "Expected positive magnetic energy in time-dependent coil run"

    return True, "Passed (finite logs, magnetic energy active)"


# ---------------------------------------------------------------------------
# Plot + probe validators
# ---------------------------------------------------------------------------
def _check_helmholtz_toroidal(frames):
    if len(frames) < 2:
        return False, "Expected >= 2 frames for helmholtz_toroidal"
    f0, f1 = frames[0], frames[-1]
    coords = _coords3(f1)

    b_s1, _ = analytic_solenoid_ba(
        coords, [0.0, 0.0, -0.125], [0.0, 0.0, 1.0], 0.25, 0.04, 4, 500.0, 0.0
    )
    b_s2, _ = analytic_solenoid_ba(
        coords, [0.0, 0.0, 0.125], [0.0, 0.0, 1.0], 0.25, 0.04, 4, 500.0, 0.0
    )
    b_tor = analytic_toroidal_b(
        coords, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.02, 0.30, 0.06
    )
    b_ref = b_s1 + b_s2 + b_tor

    b0 = _vec(f1, "B0")
    btot = _vec(f1, "B")
    eind = _vec(f1, "Eind")
    if not (np.all(np.isfinite(b0)) and np.all(np.isfinite(btot))):
        return False, "Non-finite field values"

    err_b0 = _rel_max_err(b0, b_ref)
    err_b = _rel_max_err(btot, b_ref)
    err_static = _rel_max_err(_vec(f0, "B0"), b0)
    max_eind = float(np.max(np.abs(eind)))
    if err_b0 > 1.0e-10 or err_b > 1.0e-10:
        return False, f"B0/B analytic mismatch (err_B0={err_b0:.3e}, err_B={err_b:.3e})"
    if err_static > 1.0e-12:
        return False, f"Static coil B0 drifted between frames ({err_static:.3e})"
    if max_eind > 1.0e-14:
        return False, f"Expected zero Eind for constant waveform, got {max_eind:.3e}"
    return True, f"Passed (helmholtz_toroidal max rel err = {err_b0:.2e})"


def _check_segment_loop(frames):
    f1 = frames[-1]
    coords = _coords3(f1)
    verts = [
        ([-0.2, -0.2, 0.0], [0.2, -0.2, 0.0]),
        ([0.2, -0.2, 0.0], [0.2, 0.2, 0.0]),
        ([0.2, 0.2, 0.0], [-0.2, 0.2, 0.0]),
        ([-0.2, 0.2, 0.0], [-0.2, -0.2, 0.0]),
    ]
    b_ref = np.zeros_like(coords)
    for x1, x2 in verts:
        b_ref += analytic_segment_b(coords, x1, x2, 600.0, 0.02)

    b0 = _vec(f1, "B0")
    btot = _vec(f1, "B")
    err_b0 = _rel_max_err(b0, b_ref)
    err_b = _rel_max_err(btot, b_ref)
    # ASCII format writes 7 significant digits (~1e-6).
    if err_b0 > 1.0e-5 or err_b > 1.0e-5:
        return False, f"Segment Biot-Savart mismatch (err_B0={err_b0:.3e}, err_B={err_b:.3e})"
    return True, f"Passed (segment_loop ASCII max rel err = {err_b0:.2e})"


def _check_filament2d(frames):
    if len(frames) < 4:
        return False, f"Expected >= 4 frames for filament2d, got {len(frames)}"

    # Constant reference spatial templates on 2D slice.
    coords0 = _coords3(frames[0])
    b_pos_ref, a_pos_ref = analytic_filament2d_ba(coords0, -0.25, 0.0, 1000.0, 0.05)
    b_neg_ref, a_neg_ref = analytic_filament2d_ba(coords0, 0.25, 0.0, -1000.0, 0.05)

    for fr in frames:
        t = fr["time"]
        f_lin, df_lin = wf_linear_ramp(t, 0.0, 1.0e-5, 0.0, 1.0)
        # waveform_table.txt encodes the identical linear ramp over [0, 1e-5] s.
        b_ref = f_lin * (b_pos_ref + b_neg_ref)
        e_ref = -df_lin * (a_pos_ref + a_neg_ref)

        err_b0 = _rel_max_err(_vec(fr, "B0"), b_ref)
        err_b = _rel_max_err(_vec(fr, "B"), b_ref)
        err_eind = _rel_max_err(_vec(fr, "Eind"), e_ref)
        err_etot = _rel_max_err(_vec(fr, "E"), e_ref)
        err_e0 = _rel_max_err(_vec(fr, "E0"), e_ref)
        if max(err_b0, err_b) > 1.0e-10:
            return False, f"filament2d B0/B mismatch at t={t:.3e}: {max(err_b0, err_b):.3e}"
        if t > 0.0 and t < 0.99e-5 and max(err_eind, err_etot, err_e0) > 1.0e-10:
            return False, f"filament2d E mismatch at t={t:.3e}: {max(err_eind, err_etot, err_e0):.3e}"

    # Verify discrete curl(E_z) = -dB/dt on the interior 2D grid at midpoint.
    f_mid = frames[2]
    xy = f_mid["coords"]
    xs = np.unique(xy[:, 0])
    ys = np.unique(xy[:, 1])
    nx, ny = len(xs), len(ys)
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    order = np.lexsort((xy[:, 1], xy[:, 0]))
    ez = f_mid["fields"]["Eindz"][order].reshape(nx, ny)
    bx_t0 = frames[1]["fields"]["B0x"][order].reshape(nx, ny)
    bx_t2 = frames[3]["fields"]["B0x"][order].reshape(nx, ny)
    by_t0 = frames[1]["fields"]["B0y"][order].reshape(nx, ny)
    by_t2 = frames[3]["fields"]["B0y"][order].reshape(nx, ny)
    dt2 = frames[3]["time"] - frames[1]["time"]

    curl_ex = (ez[1:-1, 2:] - ez[1:-1, :-2]) / (2.0 * dy)
    curl_ey = -(ez[2:, 1:-1] - ez[:-2, 1:-1]) / (2.0 * dx)
    neg_dbxdt = -(bx_t2[1:-1, 1:-1] - bx_t0[1:-1, 1:-1]) / dt2
    neg_dbydt = -(by_t2[1:-1, 1:-1] - by_t0[1:-1, 1:-1]) / dt2

    scale = max(float(np.max(np.abs(neg_dbxdt))), float(np.max(np.abs(neg_dbydt))), 1.0e-12)
    faraday_rms = float(
        np.sqrt(np.mean((curl_ex - neg_dbxdt) ** 2 + (curl_ey - neg_dbydt) ** 2)) / scale
    )
    if faraday_rms > 0.05:
        return False, f"filament2d discrete Faraday mismatch (rms={faraday_rms:.3e})"

    return True, f"Passed (filament2d linear+table match, Faraday rms={faraday_rms:.2e})"


def _check_waveform_static(frames):
    if len(frames) < 5:
        return False, f"Expected >= 5 frames for waveform_static, got {len(frames)}"

    coords = _coords3(frames[0])
    b_r_ref, a_r_ref = analytic_loop_ba(
        coords, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.25, 1200.0, 0.02
    )
    b_a_ref, a_a_ref = analytic_loop_ba(
        coords, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.20, 400.0, 0.02
    )

    for fr in frames:
        t = fr["time"]
        framp, dframp = wf_cosine_ramp(t, 0.0, 1.0e-5, 0.2, 1.2)
        fac, dfac = wf_sinusoid(t, 0.0, 0.5, 0.5, 5.0e4, 0.0)
        b_ref = framp * b_r_ref + fac * b_a_ref
        e_ref = -dframp * a_r_ref - dfac * a_a_ref

        err_b0 = _rel_max_err(_vec(fr, "B0"), b_ref)
        err_b = _rel_max_err(_vec(fr, "B"), b_ref)
        err_eind = _rel_max_err(_vec(fr, "Eind"), e_ref)
        err_e = _rel_max_err(_vec(fr, "E"), e_ref)
        err_e0 = _rel_max_err(_vec(fr, "E0"), e_ref)
        if max(err_b0, err_b) > 1.0e-10:
            return False, f"waveform_static B mismatch at t={t:.3e}: {max(err_b0, err_b):.3e}"
        if np.max(np.abs(e_ref)) > 1.0e-12 and max(err_eind, err_e, err_e0) > 1.0e-10:
            return False, f"waveform_static E mismatch at t={t:.3e}: {max(err_eind, err_e, err_e0):.3e}"

    # Check discrete curl(Eind) == -dB0/dt at midpoint frame index 2 (t = 2e-6 s).
    is_3d = frames[0]["coords"].shape[1] == 3 and len(np.unique(coords[:, 2])) > 1
    xs = np.unique(coords[:, 0])
    ys = np.unique(coords[:, 1])
    nx, ny = len(xs), len(ys)
    dx, dy = xs[1] - xs[0], ys[1] - ys[0]

    if is_3d:
        zs = np.unique(coords[:, 2])
        nz = len(zs)
        dz = zs[1] - zs[0]
        order = np.lexsort((coords[:, 2], coords[:, 1], coords[:, 0]))

        def grid3(fr, name):
            return fr["fields"][name][order].reshape(nx, ny, nz)

        ex = grid3(frames[2], "Eindx")
        ey = grid3(frames[2], "Eindy")
        ez = grid3(frames[2], "Eindz")
        curl_x = (ez[1:-1, 2:, 1:-1] - ez[1:-1, :-2, 1:-1]) / (2.0 * dy) - (
            ey[1:-1, 1:-1, 2:] - ey[1:-1, 1:-1, :-2]
        ) / (2.0 * dz)
        curl_y = (ex[1:-1, 1:-1, 2:] - ex[1:-1, 1:-1, :-2]) / (2.0 * dz) - (
            ez[2:, 1:-1, 1:-1] - ez[:-2, 1:-1, 1:-1]
        ) / (2.0 * dx)
        curl_z = (ey[2:, 1:-1, 1:-1] - ey[:-2, 1:-1, 1:-1]) / (2.0 * dx) - (
            ex[1:-1, 2:, 1:-1] - ex[1:-1, :-2, 1:-1]
        ) / (2.0 * dy)

        dt2 = frames[3]["time"] - frames[1]["time"]
        dbx = -(grid3(frames[3], "B0x")[1:-1, 1:-1, 1:-1] - grid3(frames[1], "B0x")[1:-1, 1:-1, 1:-1]) / dt2
        dby = -(grid3(frames[3], "B0y")[1:-1, 1:-1, 1:-1] - grid3(frames[1], "B0y")[1:-1, 1:-1, 1:-1]) / dt2
        dbz = -(grid3(frames[3], "B0z")[1:-1, 1:-1, 1:-1] - grid3(frames[1], "B0z")[1:-1, 1:-1, 1:-1]) / dt2

        scale = max(float(np.max(np.sqrt(dbx**2 + dby**2 + dbz**2))), 1.0e-12)
        faraday_rms = float(
            np.sqrt(np.mean((curl_x - dbx) ** 2 + (curl_y - dby) ** 2 + (curl_z - dbz) ** 2))
            / scale
        )
        if faraday_rms > 0.08:
            return False, f"3D discrete Faraday mismatch (rms={faraday_rms:.3e})"
        z_line = 0.05
        z_plane = 0.025
    else:
        order = np.lexsort((coords[:, 1], coords[:, 0]))

        def grid2(fr, name):
            return fr["fields"][name][order].reshape(nx, ny)

        ex = grid2(frames[2], "Eindx")
        ey = grid2(frames[2], "Eindy")
        curl_z = (ey[2:, 1:-1] - ey[:-2, 1:-1]) / (2.0 * dx) - (
            ex[1:-1, 2:] - ex[1:-1, :-2]
        ) / (2.0 * dy)
        dt2 = frames[3]["time"] - frames[1]["time"]
        dbz = -(grid2(frames[3], "B0z")[1:-1, 1:-1] - grid2(frames[1], "B0z")[1:-1, 1:-1]) / dt2
        scale = max(float(np.max(np.abs(dbz))), 1.0e-12)
        faraday_rms = float(np.sqrt(np.mean((curl_z - dbz) ** 2)) / scale)
        if faraday_rms > 0.08:
            return False, f"2D discrete Faraday mismatch (rms={faraday_rms:.3e})"
        z_line = 0.0
        z_plane = 0.0

    # Subprocess CLI validation for tools/extract_probes.py on real8 .out frames.
    probe_script = os.path.join(REPO_ROOT, "tools", "extract_probes.py")
    plot_glob = os.path.join(RUN_DIR, "PC", "plots", "*.out")
    csv_out = os.path.join(RUN_DIR, "probe_line.csv")
    npz_out = os.path.join(RUN_DIR, "probe_plane.npz")

    subprocess.run(
        [
            sys.executable,
            probe_script,
            plot_glob,
            "--line",
            f"-0.25,0.0,{z_line}",
            f"0.25,0.0,{z_line}",
            "11",
            "--vars",
            "Bx,By,Bz,Eindx,Eindy,Eindz",
            "-o",
            csv_out,
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            probe_script,
            plot_glob,
            "--plane",
            f"-0.08,-0.08,{z_plane}",
            "0.16,0,0",
            "0,0.16,0",
            "6",
            "6",
            "--vars",
            "Bx,Bz,Eindy",
            "-o",
            npz_out,
        ],
        check=True,
    )

    with open(csv_out, "r", encoding="utf-8") as cf:
        csv_rows = list(csv.DictReader(cf))
    if len(csv_rows) != len(frames) * 11:
        return False, f"Expected {len(frames) * 11} CSV probe rows, got {len(csv_rows)}"

    # Final-frame probe points sit on z=0.05 (exact grid node plane) and y=0 (exact node);
    # compare interpolated CSV values against analytic values.
    last_rows = [r for r in csv_rows if int(r["frame"]) == len(frames) - 1]
    line_pts = np.array([[float(r["x"]), float(r["y"]), float(r["z"])] for r in last_rows])
    line_b = np.array([[float(r["Bx"]), float(r["By"]), float(r["Bz"])] for r in last_rows])
    t_last = frames[-1]["time"]
    framp, _ = wf_cosine_ramp(t_last, 0.0, 1.0e-5, 0.2, 1.2)
    fac, _ = wf_sinusoid(t_last, 0.0, 0.5, 0.5, 5.0e4, 0.0)
    b_l_r, _ = analytic_loop_ba(line_pts, [0, 0, 0], [0, 1, 0], 0.25, 1200.0, 0.02)
    b_l_a, _ = analytic_loop_ba(line_pts, [0, 0, 0], [0, 0, 1], 0.20, 400.0, 0.02)
    b_line_ref = framp * b_l_r + fac * b_l_a
    if not np.all(np.isfinite(line_b)):
        return False, "Non-finite values in extract_probes CSV"
    # Points include on-node x steps of 0.05; allow float formatting tolerance.
    line_err = _rel_max_err(line_b, b_line_ref)
    if line_err > 1.0e-8:
        return False, f"CLI line probe mismatch against exact nodes: {line_err:.3e}"

    npz_data = np.load(npz_out)
    if tuple(npz_data["plane_shape"]) != (6, 6):
        return False, f"Unexpected plane_shape in NPZ: {npz_data['plane_shape']}"
    if npz_data["Bx"].shape != (len(frames), 36):
        return False, f"Unexpected NPZ array shape: {npz_data['Bx'].shape}"
    if not np.all(np.isfinite(npz_data["Bz"])):
        return False, "Non-finite values in extract_probes NPZ"
    plane_pts = npz_data["coords"]
    b_p_r, _ = analytic_loop_ba(plane_pts, [0, 0, 0], [0, 1, 0], 0.25, 1200.0, 0.02)
    b_p_a, _ = analytic_loop_ba(plane_pts, [0, 0, 0], [0, 0, 1], 0.20, 400.0, 0.02)
    bz_plane_ref = (framp * b_p_r + fac * b_p_a)[:, 2]
    bz_plane_err = float(
        np.max(np.abs(npz_data["Bz"][-1] - bz_plane_ref))
        / max(float(np.max(np.abs(bz_plane_ref))), 1.0e-12)
    )
    if bz_plane_err > 0.05:
        return False, f"CLI plane probe Bz mismatch: {bz_plane_err:.3e}"

    return True, f"Passed (waveform_static exact + Faraday rms={faraday_rms:.2e} + CLI probes)"


def _check_waveform_hybrid(frames):
    if len(frames) < 5:
        return False, f"Expected >= 5 frames for waveform_hybrid, got {len(frames)}"

    coords0 = _coords3(frames[0])
    b_pos_ref, a_pos_ref = analytic_filament2d_ba(coords0, -0.20, 0.0, 400.0, 0.06)
    b_neg_ref, a_neg_ref = analytic_filament2d_ba(coords0, 0.20, 0.0, -400.0, 0.06)
    b_base = b_pos_ref + b_neg_ref
    a_base = a_pos_ref + a_neg_ref

    saw_eddy = False
    for fr in frames:
        for vname in ("B0x", "B0y", "B0z", "Bx", "By", "Bz", "Eindx", "Eindy", "Eindz", "E0x", "E0y", "E0z", "rhoS0"):
            if not np.all(np.isfinite(fr["fields"][vname])):
                return False, f"Non-finite {vname} in hybrid frame at t={fr['time']:.3e}"

        t = fr["time"]
        f_rlc, df_rlc = wf_crowbar_rlc(t, 0.0, 0.4e-5, 2.0e-5, True, 1.0e-5)
        b_ref = f_rlc * b_base
        e_ref = -df_rlc * a_base

        err_b0 = _rel_max_err(_vec(fr, "B0"), b_ref)
        err_eind = _rel_max_err(_vec(fr, "Eind"), e_ref)
        err_e0 = _rel_max_err(_vec(fr, "E0"), e_ref)
        if err_b0 > 1.0e-10:
            return False, f"Hybrid B0 crowbar_rlc mismatch at t={t:.3e}: {err_b0:.3e}"
        if abs(df_rlc) > 0.0 and max(err_eind, err_e0) > 1.0e-10:
            return False, f"Hybrid Eind/E0 mismatch at t={t:.3e}: {max(err_eind, err_e0):.3e}"

        b1 = _vec(fr, "B") - _vec(fr, "B0")
        b0 = _vec(fr, "B0")
        if t > 0.0 and float(np.max(np.linalg.norm(b1, axis=1))) > 0.0:
            # Faraday splitting dB1/dt = -curl(E_Ohm - E0) drives an eddy B1 opposing rising B0.
            if float(np.sum(b1 * b0)) < 0.0:
                saw_eddy = True

    if not saw_eddy:
        return False, "Expected plasma eddy B1 opposing rising coil B0 in hybrid run"

    return True, "Passed (hybrid crowbar_rlc + Eind/E0 + Faraday shielding response)"


def validate_plot(test_name):
    """Dispatch plot validation by coil_field variant name."""
    frames = _load_frames()
    if not frames:
        return False, "No .out frames produced"

    if test_name == "coil_field_helmholtz_toroidal":
        return _check_helmholtz_toroidal(frames)
    if test_name == "coil_field_segment_loop":
        return _check_segment_loop(frames)
    if test_name == "coil_field_filament2d":
        return _check_filament2d(frames)
    if test_name == "coil_field_waveform_static":
        return _check_waveform_static(frames)
    if test_name == "coil_field_waveform_hybrid":
        return _check_waveform_hybrid(frames)
    return False, f"Unknown coil_field variant '{test_name}'"
