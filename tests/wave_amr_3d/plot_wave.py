#!/usr/bin/env python3
"""Plot 3D AMR traveling wave simulation results.

Generates:
  1. A 2D contour plot of By on the cut plane (z=0) with the AMR box overlay.
  2. A 1D slice along y=0 showing the transverse magnetic field components (By, Bz).
  3. Optionally assembles an MP4 animation across all time steps if ffmpeg is available.

Usage:
    python3 plot_wave.py [run_dir_or_plots_dir] [--animate]
"""

import argparse
import glob
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

NHEADER = 5


def load_frame(path):
    """Load a fluid cut-plane .out file."""
    try:
        with open(path, "r", encoding="latin-1") as f:
            lines = f.read().splitlines()
    except Exception as e:
        print(f"Error opening {path}: {e}")
        return None

    if len(lines) <= NHEADER:
        return None

    time = 0.0
    try:
        time = float(lines[1].split()[1])
    except Exception:
        pass

    all_names = lines[NHEADER - 1].split()
    rows = []
    for ln in lines[NHEADER:]:
        c = ln.split()
        if len(c) >= 5:
            rows.append(c)

    if not rows:
        return None

    ncol = min(len(rows[0]), len(all_names))
    data = np.array([r[:ncol] for r in rows], dtype=float)
    result = {all_names[i]: data[:, i] for i in range(ncol)}
    result["time"] = time
    return result


def plot_frame(frame, output_path, amr_box=(-2.0, 2.0, -1.0, 1.0)):
    """Render a 2-panel plot (2D cut plane + 1D slice) for a given frame."""
    x = frame["x"]
    y = frame["y"]
    By = frame.get("By")
    Bz = frame.get("Bz")
    time = frame.get("time", 0.0)

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8, 6), gridspec_kw={"height_ratios": [1.2, 1]}
    )

    # 1. 2D contour of By on cut plane z=0
    levels = np.linspace(-0.25, 0.25, 25)
    cnt = ax1.tricontourf(x, y, By, levels=levels, cmap="RdBu_r", extend="both")
    cbar = fig.colorbar(cnt, ax=ax1, pad=0.02)
    cbar.set_label("$B_y$ [code units]")

    # Overlay AMR Level 1 box
    xmin, xmax, ymin, ymax = amr_box
    rect = Rectangle(
        (xmin, ymin),
        xmax - xmin,
        ymax - ymin,
        fill=False,
        edgecolor="black",
        linestyle="--",
        linewidth=1.5,
        label="AMR Level 1 Box",
    )
    ax1.add_patch(rect)
    ax1.set_xlim(-4.0, 4.0)
    ax1.set_ylim(-2.0, 2.0)
    ax1.set_xlabel("$x / d_i$")
    ax1.set_ylabel("$y / d_i$")
    ax1.set_title(f"3D AMR Traveling Wave ($z=0$) at $t = {time:.2f}$ s")
    ax1.legend(loc="upper right", framealpha=0.85)

    # 2. 1D slice along y ~ 0
    mask = np.abs(y) < 0.15
    xm = x[mask]
    Bym = By[mask]
    order = np.argsort(xm)

    ax2.plot(xm[order], Bym[order], "b.-", lw=1.5, ms=4, label="$B_y(x, y=0)$")
    if Bz is not None:
        Bzm = Bz[mask]
        ax2.plot(xm[order], Bzm[order], "g.--", lw=1.2, ms=3, label="$B_z(x, y=0)$")

    ax2.axvspan(xmin, xmax, color="gray", alpha=0.15, label="AMR Region ($2\\times$)")
    ax2.set_xlim(-4.0, 4.0)
    ax2.set_ylim(-0.26, 0.26)
    ax2.set_xlabel("$x / d_i$")
    ax2.set_ylabel("$B_\\perp$")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right", framealpha=0.85)

    fig.tight_layout()
    plt.savefig(output_path, dpi=130)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "target",
        nargs="?",
        default=".",
        help="Run directory containing PC/plots, or direct plots directory",
    )
    parser.add_argument(
        "--animate", action="store_true", help="Generate an MP4 animation across all frames"
    )
    args = parser.parse_args()

    plots_dir = (
        os.path.join(args.target, "PC", "plots")
        if os.path.isdir(os.path.join(args.target, "PC", "plots"))
        else args.target
    )

    out_files = sorted(glob.glob(os.path.join(plots_dir, "*_fluid*.out")))
    if not out_files:
        print(f"No fluid .out files found in {plots_dir}")
        sys.exit(1)

    print(f"Found {len(out_files)} fluid files in {plots_dir}")

    # Plot final frame as static figure
    last_frame = load_frame(out_files[-1])
    if last_frame:
        out_png = "wave_amr_3d.png"
        plot_frame(last_frame, out_png)
        print(f"Saved snapshot to {out_png}")

    if args.animate or len(out_files) > 1:
        tmp_dir = "frames_tmp"
        os.makedirs(tmp_dir, exist_ok=True)
        print(f"Rendering {len(out_files)} frames...")
        for i, fpath in enumerate(out_files):
            fr = load_frame(fpath)
            if fr:
                plot_frame(fr, os.path.join(tmp_dir, f"frame_{i:04d}.png"))

        video_path = "wave_amr_3d.mp4"
        cmd = f"ffmpeg -y -framerate 4 -i {tmp_dir}/frame_%04d.png -c:v libx264 -pix_fmt yuv420p {video_path} > /dev/null 2>&1"
        ret = os.system(cmd)
        if ret == 0:
            print(f"Successfully generated animation: {video_path}")
            # Clean up temp frames
            for f in glob.glob(os.path.join(tmp_dir, "*.png")):
                os.remove(f)
            os.rmdir(tmp_dir)
        else:
            print(f"ffmpeg failed or not found. Individual frames kept in {tmp_dir}/")


if __name__ == "__main__":
    main()
