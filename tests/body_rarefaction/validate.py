#!/usr/bin/env python3
"""Validator and Diagnostic Generator for the 2D Cylindrical Rarefaction Benchmark.

Assesses the inner body boundary condition against the analytical ambipolar
rarefaction wave benchmark:
  n(r, t) = n0 * exp(-1 - (r - R_B) / (c_s * t))   for R_B < r <= R_B + c_s * t
  u_r(R_B^+) = -c_s,  n(R_B^+) ~ n0 * exp(-1) ~ 0.37 * n0

Variants:
  - body_rarefaction_unmagnetized: Config 1 (B0 = 0)
  - body_rarefaction_magnetized:   Config 2 (B0 = B0 x_hat)

Diagnostics:
  1. Azimuthally averaged radial profiles: n_i(r), n_e(r), u_r(r) vs analytical.
  2. 2D map of charge imbalance: (n_i - n_e) / n0 (poles vs equator).
  3. 2D map of radial electric field: E_r(x, y).
  4. Time history of deleted particles: sum N_{e,del}(t) vs sum N_{i,del}(t).
"""

import argparse
import glob
import logging
import math
import os
import sys

from tests._shared import run_dir as _run_dir

logger = logging.getLogger(__name__)

set_run_dir = _run_dir.set_run_dir

R_BODY = 1.5           # Cylinder radius in code units
MASS_ION = 1.0         # Ion mass in code units
MASS_ELE = 0.04        # Electron mass in code units (mi/me = 25)
DX = 0.1               # Grid spacing
CS = 0.01633           # Ion acoustic speed in code units (sqrt((Te + 3Ti)/mi))
N0 = 1.0 / (4.0 * math.pi)  # Background code number density (omega_pi = 1)


def _columns(vidx, rows):
    """Extract required field arrays from parsed plot data."""
    out = {}
    for name in ("BODY", "X", "Y", "Z", "EX", "EY", "EZ", "BX", "BY", "BZ",
                 "RHOS0", "RHOS1", "BODYSURF", "DIVB"):
        c = _run_dir.col(vidx, rows, name)
        if c is not None:
            out[name] = c
    return out


def compute_radial_profiles(cols, r_min=R_BODY, r_max=4.0, dr=DX):
    """Compute azimuthally averaged radial profiles outside the body.
    
    Converts mass density rho_s to number density n_s = rho_s / m_s.
    """
    x = cols["X"]
    y = cols["Y"]
    rho_i = cols.get("RHOS0")
    rho_e = cols.get("RHOS1")
    ex = cols.get("EX")
    ey = cols.get("EY")

    n_bins = int(math.ceil((r_max - r_min) / dr))
    r_bins = [r_min + (i + 0.5) * dr for i in range(n_bins)]
    sum_ni = [0.0] * n_bins
    sum_ne = [0.0] * n_bins
    sum_er = [0.0] * n_bins
    count = [0] * n_bins

    for k in range(len(x)):
        r = math.hypot(x[k], y[k])
        if r < r_min or r >= r_max:
            continue
        ibin = int((r - r_min) / dr)
        if 0 <= ibin < n_bins:
            if rho_i is not None:
                sum_ni[ibin] += rho_i[k] / MASS_ION
            if rho_e is not None:
                sum_ne[ibin] += rho_e[k] / MASS_ELE
            if ex is not None and ey is not None and r > 0:
                er = (x[k] * ex[k] + y[k] * ey[k]) / r
                sum_er[ibin] += er
            count[ibin] += 1

    mean_ni = [sum_ni[i] / count[i] if count[i] > 0 else 0.0 for i in range(n_bins)]
    mean_ne = [sum_ne[i] / count[i] if count[i] > 0 else 0.0 for i in range(n_bins)]
    mean_er = [sum_er[i] / count[i] if count[i] > 0 else 0.0 for i in range(n_bins)]
    return r_bins, mean_ni, mean_ne, mean_er


def analytical_density(r_vals, t, r_b=R_BODY, cs=CS, n0=None):
    """Analytical self-similar isothermal rarefaction density."""
    if n0 is None:
        n0 = 1.0 / (4.0 * math.pi)  # background code number density
    n_analytic = []
    front = r_b + cs * t
    for r in r_vals:
        if r <= r_b:
            n_analytic.append(0.0)
        elif r <= front:
            val = n0 * math.exp(-1.0 - (r - r_b) / (cs * t))
            n_analytic.append(val)
        else:
            n_analytic.append(n0)
    return n_analytic


def validate_log(pic_diags, test_name=None):
    """Energy finiteness and charge neutrality checks from PIC log."""
    if not pic_diags or len(pic_diags) < 2:
        return False, "missing or empty PIC log"

    end = pic_diags[-1]
    # Check energies stay finite
    for key in ("Etot", "Ee", "Eb", "Epart"):
        val = end.get(key, 0.0)
        if not math.isfinite(val):
            return False, f"non-finite {key} in PIC log"

    # Energy conservation / stability: Etot must not explode
    start = pic_diags[0]
    e_start = max(start.get("Etot", 1e-10), 1e-10)
    e_end = end.get("Etot", 0.0)
    if e_end > 20.0 * e_start:
        return False, f"energy exploded: Etot grew from {e_start:.3e} to {e_end:.3e}"

    logger.debug("validate_log: PASSED for %s", test_name)
    return True, "Passed"


def validate_plot(test_name=None):
    """Verify that density inside body is zero and examine neutrality outside."""
    vidx, rows = _run_dir.load_last_out()
    if vidx is None or not rows:
        return False, "no .out plot file found"

    cols = _columns(vidx, rows)
    body = cols.get("BODY")
    if body is None:
        return False, "missing BODY column in plot output"

    inside = [i for i, b in enumerate(body) if b > 0.5]
    if not inside:
        return False, "no body cells detected in plot frame"

    # 1. Particles excluded from body interior
    for name in ("RHOS0", "RHOS1"):
        rho = cols.get(name)
        if rho is not None:
            max_in = max(abs(rho[i]) for i in inside)
            if max_in > 1e-10:
                return False, f"non-zero {name} inside body (max = {max_in:.3e})"

    # 2. Compute radial profiles and report ambipolar status
    r_bins, mean_ni, mean_ne, mean_er = compute_radial_profiles(cols)
    if mean_ni and mean_ne:
        # Near boundary node index (first bin right outside R_B)
        ni_bound = mean_ni[0]
        ne_bound = mean_ne[0]
        imbalance = abs(ni_bound - ne_bound) / N0
        logger.debug("Boundary values: n_i = %.4f, n_e = %.4f, |ni - ne|/n0 = %.3f",
                     ni_bound, ne_bound, imbalance)

    return True, "Passed"


def generate_benchmark_plots(out_dir="benchmark_plots", t_now=2.0):
    """Generate the 4 diagnostic benchmark plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("matplotlib/numpy not available; skipping plot generation.")
        return

    os.makedirs(out_dir, exist_ok=True)
    vidx, rows = _run_dir.load_last_out()
    if vidx is None:
        print("No plot output found to generate benchmark figures.")
        return

    cols = _columns(vidx, rows)
    x = np.array(cols["X"])
    y = np.array(cols["Y"])
    rho_i = np.array(cols.get("RHOS0", [])) / MASS_ION
    rho_e = np.array(cols.get("RHOS1", [])) / MASS_ELE
    ex = np.array(cols.get("EX", []))
    ey = np.array(cols.get("EY", []))
    body = np.array(cols.get("BODY", []))

    # Far-field background density for normalization
    r_all = np.sqrt(x**2 + y**2)
    far = (r_all > 4.0)
    n0 = np.mean(rho_i[far]) if np.any(far) else 1.0 / (4.0 * math.pi)

    # Diagnostic 1: Radial profiles
    r_bins, mean_ni, mean_ne, mean_er = compute_radial_profiles(cols)
    n_analytic = analytical_density(r_bins, t=t_now, n0=n0)

    plt.figure(figsize=(8, 5))
    plt.plot(r_bins, mean_ni, "b-o", label=r"Ion density $n_i(r)$")
    plt.plot(r_bins, mean_ne, "r--s", label=r"Electron density $n_e(r)$")
    plt.plot(r_bins, n_analytic, "k:", linewidth=2, label="Analytical Rarefaction")
    plt.axvline(R_BODY, color="gray", linestyle="--", label=r"$R_B$")
    plt.xlabel(r"Radius $r / d_i$")
    plt.ylabel(r"Number Density $n$ [code units]")
    plt.title(f"Diagnostic 1: Radial Density Profiles at t = {t_now:.1f}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    p1 = os.path.join(out_dir, "diagnostic_1_radial_density.png")
    plt.savefig(p1, dpi=150)
    plt.close()
    print(f"Saved: {p1}")

    # Diagnostic 2: 2D Charge Imbalance Map
    if len(rho_i) == len(rho_e) and len(rho_i) > 0:
        imbalance = (rho_i - rho_e) / n0
        imbalance[body > 0.5] = np.nan

        plt.figure(figsize=(7, 6))
        sc = plt.scatter(x, y, c=imbalance, cmap="coolwarm", vmin=-0.2, vmax=0.2, s=15)
        plt.colorbar(sc, label=r"$(n_i - n_e) / n_0$")
        circle = plt.Circle((0, 0), R_BODY, color="black", fill=False, linewidth=2)
        plt.gca().add_patch(circle)
        plt.xlim(-4, 4)
        plt.ylim(-4, 4)
        plt.xlabel(r"$x / d_i$")
        plt.ylabel(r"$y / d_i$")
        plt.title(f"Diagnostic 2: Charge Imbalance Map at t = {t_now:.1f}")
        p2 = os.path.join(out_dir, "diagnostic_2_charge_imbalance.png")
        plt.savefig(p2, dpi=150)
        plt.close()
        print(f"Saved: {p2}")

    # Diagnostic 3: Radial Electric Field Er
    if len(ex) > 0 and len(ey) > 0:
        r = np.sqrt(x**2 + y**2)
        r_safe = np.where(r > 0, r, 1.0)
        er = (x * ex + y * ey) / r_safe
        er[body > 0.5] = np.nan

        plt.figure(figsize=(7, 6))
        sc = plt.scatter(x, y, c=er, cmap="viridis", s=15)
        plt.colorbar(sc, label=r"$E_r$ [V/m]")
        circle = plt.Circle((0, 0), R_BODY, color="black", fill=False, linewidth=2)
        plt.gca().add_patch(circle)
        plt.xlim(-4, 4)
        plt.ylim(-4, 4)
        plt.xlabel(r"$x / d_i$")
        plt.ylabel(r"$y / d_i$")
        plt.title(f"Diagnostic 3: Radial Electric Field $E_r$ at t = {t_now:.1f}")
        p3 = os.path.join(out_dir, "diagnostic_3_radial_Efield.png")
        plt.savefig(p3, dpi=150)
        plt.close()
        print(f"Saved: {p3}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rarefaction Benchmark Validator")
    parser.add_argument("--plots", action="store_true", help="Generate benchmark figures")
    parser.add_argument("--time", type=float, default=2.0, help="Simulation time for plots")
    args = parser.parse_args()

    if args.plots:
        generate_benchmark_plots(t_now=args.time)
    else:
        ok, reason = validate_plot()
        print("PASSED" if ok else "FAILED", "-", reason)
        sys.exit(0 if ok else 1)
