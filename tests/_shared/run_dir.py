#!/usr/bin/env python3
"""Shared run-directory and plot helpers for the validation scripts."""

import glob
import os

RUN_DIR = "run_test"


def set_run_dir(run_dir):
    """Keep the shared run directory in sync with the active test runner."""
    global RUN_DIR
    RUN_DIR = run_dir


def plot_files(run_dir=None, pattern="*.out", plots_subdir="PC/plots"):
    """Return all matching plot frames for the active run directory."""
    base_dir = run_dir or RUN_DIR
    plots_dir = os.path.join(base_dir, plots_subdir)
    return sorted(glob.glob(os.path.join(plots_dir, pattern)))


def load_last_out(run_dir=None, pattern="*.out", plots_subdir="PC/plots"):
    """Return the final .out frame in the run directory as (vidx, rows)."""
    out_files = plot_files(run_dir=run_dir, pattern=pattern, plots_subdir=plots_subdir)
    if not out_files:
        return None, None
    return _read_out_file(out_files[-1])


def _read_out_file(out_file):
    """Parse a standard FLEKS PostProc .out file into column map and rows."""
    with open(out_file, "r", encoding="latin-1") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return None, None

    vidx = {v.upper(): i for i, v in enumerate(lines[4].split())}
    rows = []
    for line in lines[5:]:
        cols = line.split()
        if not cols:
            continue
        try:
            rows.append([float(c) for c in cols])
        except ValueError:
            continue
    return vidx, rows
