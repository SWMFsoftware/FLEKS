#!/usr/bin/env python3
"""Validator for tests/restart (implicit, explicit, hybrid PIC restarts).

Flow per variant:
  1. Verify the uninterrupted two-session run reached cycle 10 with valid logs.
  2. Verify cycle-5 checkpoint files exist in PC/restartOUT/.
  3. Snapshot the uninterrupted cycle-10 real8 plot fields from PC/plots/.
  4. Copy PC/restartOUT/ -> PC/restartIN/, generate a standalone restart
     PARAM.in using `#INCLUDE PC/restartIN/FLEKS1_restart.H`, run FLEKS.exe
     (matching the MPI rank count in fleks_stdout.log) + PostProc.pl -v in a
     clean PC/plots/ directory.
  5. Compare overlapping energy-log cycles (5..10) and cycle-10 plot fields
     between the restarted run and uninterrupted reference (< 1e-10 rel diff).
"""
import glob
import logging
import math
import os
import re
import shutil
import struct
import subprocess

import numpy as np

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"
_PLOT_RESULT = None

RTOL = 1e-10
ATOL = 1e-14

_EXCLUDED_RESTART_COMMANDS = {
    "#RESTART",
    "#NSTEP",
    "#TIMESIMULATION",
    "#TIMESTEP",
    "#DT",
    "#GEOMETRY",
    "#NCELL",
    "#MAXBLOCKSIZE",
    "#ELECTRON",
    "#PARTICLES",
    "#SUPID",
    "#TIMESTEPPING",
    "#TESTCASE",
    "#WAVEPERTURBATION",
    "#WAVEIC",
    "#GEMIC",
    "#BEAMIC",
    "#TOPHAT",
    "#FADEEV",
    "#FORCEFREE",
    "#SAVERESTART",
    "#STOP",
    "#RUN",
}


def set_run_dir(run_dir):
    """Point validator at the current run directory and reset cached state."""
    global RUN_DIR, _PLOT_RESULT
    RUN_DIR = run_dir
    _PLOT_RESULT = None


def _rel_close(a, b, rtol=RTOL, atol=ATOL):
    """Return True if arrays/scalars a and b are finite and match within rtol/atol."""
    a_arr = np.asarray(a, dtype=np.float64)
    b_arr = np.asarray(b, dtype=np.float64)
    if not (np.all(np.isfinite(a_arr)) and np.all(np.isfinite(b_arr))):
        return False
    diff = np.abs(a_arr - b_arr)
    scale = np.maximum(np.abs(a_arr), np.abs(b_arr))
    return bool(np.all(diff <= atol + rtol * scale))


def _max_rel_err(a, b):
    """Return max relative difference with ATOL floor, or inf if non-finite."""
    a_arr = np.asarray(a, dtype=np.float64)
    b_arr = np.asarray(b, dtype=np.float64)
    if not (np.all(np.isfinite(a_arr)) and np.all(np.isfinite(b_arr))):
        return float("inf")
    diff = np.abs(a_arr - b_arr)
    denom = np.maximum(np.maximum(np.abs(a_arr), np.abs(b_arr)), ATOL)
    return float(np.max(diff / denom))


def _read_rec(f):
    """Read one Fortran unformatted sequential record."""
    hdr = f.read(4)
    if len(hdr) < 4:
        raise ValueError("Unexpected EOF reading Fortran record header")
    (nbytes,) = struct.unpack("<i", hdr)
    payload = f.read(nbytes)
    if len(payload) < nbytes:
        raise ValueError(
            f"Truncated Fortran record payload: expected {nbytes}, got {len(payload)}"
        )
    ftr = f.read(4)
    if len(ftr) < 4:
        raise ValueError("Unexpected EOF reading Fortran record trailer")
    (nbytes2,) = struct.unpack("<i", ftr)
    if nbytes2 != nbytes:
        raise ValueError(f"Record marker mismatch: {nbytes} != {nbytes2}")
    return payload


def _read_real8_out(path):
    """Parse a real8 Fortran unformatted IDL .out plot file."""
    with open(path, "rb") as f:
        _header_str = _read_rec(f)
        n_step, sim_time, n_dim_out, n_param, n_var = struct.unpack(
            "<idiii", _read_rec(f)
        )
        n_dim = abs(n_dim_out)
        n_sizes = struct.unpack(f"<{n_dim}i", _read_rec(f))
        n_point = math.prod(n_sizes)
        params = (
            np.frombuffer(_read_rec(f), dtype="<f8").copy()
            if n_param > 0
            else np.array([], dtype=np.float64)
        )
        names = _read_rec(f).decode("ascii", errors="replace").split()
        coords = (
            np.frombuffer(_read_rec(f), dtype="<f8")
            .reshape(n_dim, n_point)
            .copy()
        )
        fields = [
            np.frombuffer(_read_rec(f), dtype="<f8").copy()
            for _ in range(n_var)
        ]

    return {
        "nStep": n_step,
        "time": sim_time,
        "dims": n_sizes,
        "params": params,
        "names": names[: n_dim + n_var],
        "coords": coords,
        "fields": fields,
    }


def _find_cycle10_out(plots_dir):
    """Locate the cycle-10 .out plot file inside plots_dir."""
    candidates = sorted(glob.glob(os.path.join(plots_dir, "*n*10.out")))
    if not candidates:
        candidates = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    for path in reversed(candidates):
        data = _read_real8_out(path)
        if data["nStep"] == 10:
            return path, data
    raise FileNotFoundError(
        f"No cycle-10 real8 .out plot file found in {plots_dir} "
        f"(candidates: {[os.path.basename(p) for p in candidates]})"
    )


def _read_pic_log_file(log_path):
    """Parse a FLEKS PC/plots/log_*.log diagnostic table into row dicts."""
    with open(log_path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if len(lines) < 2:
        raise ValueError(f"Diagnostic log {log_path} has too few lines")

    def _is_number(tok):
        try:
            float(tok)
            return True
        except ValueError:
            return False

    # Find the column header row (last non-numeric row before data).
    hdr_idx = 0
    for idx, ln in enumerate(lines):
        toks = ln.split()
        if toks and _is_number(toks[0]):
            hdr_idx = max(0, idx - 1)
            break

    cols = lines[hdr_idx].split()
    rows = []
    for ln in lines[hdr_idx + 1 :]:
        parts = ln.split()
        if len(parts) < min(len(cols), 6) or not _is_number(parts[0]):
            continue
        row = {}
        for idx, (name, tok) in enumerate(zip(cols, parts)):
            if idx == 1 or name.lower() in ("it", "cycle", "step", "nstep", "istep"):
                row["cycle"] = int(float(tok))
            else:
                row[name] = float(tok)
        if "cycle" not in row and len(parts) >= 2:
            row["cycle"] = int(float(parts[1]))
        # Standard energy columns (indices 2..5: Etot, Ee, Eb, Epart).
        if len(parts) >= 6:
            row.setdefault("Etot", float(parts[2]))
            row.setdefault("Ee", float(parts[3]))
            row.setdefault("Eb", float(parts[4]))
            row.setdefault("Epart", float(parts[5]))
        rows.append(row)
    return rows


def _check_restart_out(restart_out_dir, test_name):
    """Verify required checkpoint header and MultiFab/particle files exist."""
    if not os.path.isdir(restart_out_dir):
        return f"Missing restart output directory {restart_out_dir}"

    required_files = ["FLEKS1_restart.H", "FLEKS1_amrex_restart.H"]
    for fname in required_files:
        if not os.path.isfile(os.path.join(restart_out_dir, fname)):
            return f"Missing {fname} in {restart_out_dir}"

    mf_prefixes = ["FLEKS1_nodeE", "FLEKS1_nodeB", "FLEKS1_centerB"]
    if test_name and "hybrid" in test_name.lower():
        mf_prefixes.extend(["FLEKS1_nodePlasmaPrev", "FLEKS1_centerPeState"])

    for prefix in mf_prefixes:
        matches = glob.glob(os.path.join(restart_out_dir, f"{prefix}*_H"))
        if not matches:
            return f"Missing MultiFab header matching {prefix}*_H in {restart_out_dir}"

    part_hdr = os.path.join(restart_out_dir, "FLEKS1_particles0", "Header")
    if not os.path.isfile(part_hdr):
        return f"Missing particle checkpoint header {part_hdr}"

    return None


def _detect_nprocs(run_dir):
    """Detect MPI process count used by the initial run from fleks_stdout.log."""
    stdout_log = os.path.join(run_dir, "fleks_stdout.log")
    if os.path.isfile(stdout_log):
        with open(stdout_log, "r", errors="replace") as f:
            text = f.read()
        m = re.search(r"MPI initialized with\s+(\d+)\s+MPI\s+process", text)
        if m:
            return max(1, int(m.group(1)))
    return 1


def _build_restart_param(orig_param_text):
    """Construct standalone restart PARAM.in text from session 1 of orig_param_text."""
    lines = orig_param_text.splitlines()
    session1_lines = []
    for ln in lines:
        toks = ln.split()
        if toks and toks[0].upper() == "#RUN":
            break
        session1_lines.append(ln)

    blocks = []
    cur_cmd = None
    cur_lines = []
    for ln in session1_lines:
        stripped = ln.strip()
        if stripped.startswith("#"):
            if cur_cmd is not None or cur_lines:
                blocks.append((cur_cmd, cur_lines))
            cur_cmd = stripped.split()[0].upper()
            cur_lines = [ln]
        else:
            cur_lines.append(ln)
    if cur_cmd is not None or cur_lines:
        blocks.append((cur_cmd, cur_lines))

    kept = []
    for cmd, blines in blocks:
        if cmd in _EXCLUDED_RESTART_COMMANDS:
            continue
        kept.append("\n".join(blines))

    parts = [
        "#INCLUDE\nPC/restartIN/FLEKS1_restart.H",
        *kept,
        "#SAVERESTART\nF                       doSaveRestart",
        "#STOP\n10                      MaxIter\n-1.0                    TimeMax",
    ]
    return "\n\n".join(parts) + "\n"


def validate_log(pic_diags=None, test_name=None):
    """Run checkpoint-restart comparison and validate energy logs (cycles 5..10)."""
    global _PLOT_RESULT
    _PLOT_RESULT = None

    if not pic_diags or len(pic_diags) < 11:
        got = len(pic_diags) if pic_diags else 0
        return False, f"Uninterrupted run log incomplete (expected >= 11 entries, got {got})"

    if pic_diags[-1].get("cycle") != 10:
        return False, f"Expected final uninterrupted cycle 10, got {pic_diags[-1].get('cycle')}"

    for row in pic_diags:
        for key in ("Etot", "Eb", "Ee", "Epart"):
            if key in row and not math.isfinite(row[key]):
                return False, f"Non-finite {key} in uninterrupted log at cycle {row.get('cycle')}"

    restart_out = os.path.join(RUN_DIR, "PC", "restartOUT")
    err = _check_restart_out(restart_out, test_name)
    if err is not None:
        return False, err

    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    try:
        _, ref_plot = _find_cycle10_out(plots_dir)
    except Exception as exc:
        return False, f"Failed to read uninterrupted cycle-10 plot: {exc}"

    param_path = os.path.join(RUN_DIR, "PARAM.in")
    with open(param_path, "r") as f:
        orig_param = f.read()

    restart_param = _build_restart_param(orig_param)
    restart_in = os.path.join(RUN_DIR, "PC", "restartIN")
    plots_backup = os.path.join(RUN_DIR, "PC", "plots_uninterrupted")
    nprocs = _detect_nprocs(RUN_DIR)

    try:
        if os.path.exists(restart_in):
            shutil.rmtree(restart_in)
        shutil.copytree(restart_out, restart_in)

        if os.path.exists(plots_backup):
            shutil.rmtree(plots_backup)
        shutil.move(plots_dir, plots_backup)
        os.makedirs(plots_dir, exist_ok=True)

        with open(param_path, "w") as f:
            f.write(restart_param)

        cmd = (
            ["mpirun", "-n", str(nprocs), "./FLEKS.exe"]
            if nprocs > 1
            else ["./FLEKS.exe"]
        )
        logger.debug("Running restart phase in %s: %s", RUN_DIR, " ".join(cmd))
        res = subprocess.run(
            cmd,
            cwd=RUN_DIR,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if res.returncode != 0:
            tail = (res.stdout + "\n" + res.stderr)[-1200:]
            return False, f"Restarted FLEKS.exe failed (rc={res.returncode}):\n{tail}"

        pp_cmd = ["./PostProc.pl", "-v"]
        if not os.path.exists(os.path.join(RUN_DIR, "PostProc.pl")):
            pp_cmd = ["perl", "PostProc.pl", "-v"]
        pp_res = subprocess.run(
            pp_cmd,
            cwd=RUN_DIR,
            capture_output=True,
            text=True,
            timeout=60,
        )
        if pp_res.returncode != 0:
            tail = (pp_res.stdout + "\n" + pp_res.stderr)[-1200:]
            return False, f"PostProc.pl failed after restart (rc={pp_res.returncode}):\n{tail}"

        log_files = sorted(glob.glob(os.path.join(plots_dir, "log_*.log")))
        if not log_files:
            return False, "No diagnostic log_*.log written by restarted run"

        rst_rows = _read_pic_log_file(log_files[-1])
        if not rst_rows or rst_rows[-1]["cycle"] != 10:
            last_c = rst_rows[-1]["cycle"] if rst_rows else None
            return False, f"Restarted run ended at cycle {last_c} (expected 10)"

        ref_by_cycle = {r["cycle"]: r for r in pic_diags}
        expected_cycles = list(range(5, 11))
        overlap_cycles = [r["cycle"] for r in rst_rows]
        if overlap_cycles != expected_cycles:
            return (
                False,
                f"Expected restarted log cycles {expected_cycles}, got {overlap_cycles}",
            )

        max_log_err = 0.0
        for rrow in rst_rows:
            cyc = rrow["cycle"]
            if cyc not in ref_by_cycle:
                return False, f"Cycle {cyc} missing from uninterrupted reference log"
            urow = ref_by_cycle[cyc]
            for key in ("Etot", "Eb", "Ee", "Epart"):
                if key not in rrow or key not in urow:
                    return (
                        False,
                        f"Missing required energy column '{key}' at cycle {cyc}",
                    )
                err_val = _max_rel_err(rrow[key], urow[key])
                max_log_err = max(max_log_err, err_val)
                if not _rel_close(rrow[key], urow[key]):
                    return (
                        False,
                        f"Log mismatch at cycle {cyc} {key}: "
                        f"restart={rrow[key]:.15e}, uninterrupted={urow[key]:.15e} "
                        f"(rel_err={err_val:.2e})",
                    )

        try:
            _, rst_plot = _find_cycle10_out(plots_dir)
        except Exception as exc:
            _PLOT_RESULT = (False, f"Failed to read restarted cycle-10 plot: {exc}")
        else:
            if not _rel_close(rst_plot["coords"], ref_plot["coords"]):
                _PLOT_RESULT = (False, "Cycle-10 plot coordinates differ after restart")
            elif len(rst_plot["fields"]) != len(ref_plot["fields"]):
                _PLOT_RESULT = (
                    False,
                    f"Cycle-10 variable count mismatch: {len(rst_plot['fields'])} vs {len(ref_plot['fields'])}",
                )
            else:
                max_field_err = 0.0
                worst_var = ""
                for idx, (ra, ua) in enumerate(zip(rst_plot["fields"], ref_plot["fields"])):
                    verr = _max_rel_err(ra, ua)
                    if verr > max_field_err:
                        max_field_err = verr
                        n_dim = rst_plot["coords"].shape[0]
                        worst_var = (
                            rst_plot["names"][n_dim + idx]
                            if n_dim + idx < len(rst_plot["names"])
                            else f"var{idx}"
                        )
                    if not _rel_close(ra, ua):
                        _PLOT_RESULT = (
                            False,
                            f"Cycle-10 field mismatch in {worst_var}: "
                            f"max_rel_err={verr:.2e} > {RTOL:.1e}",
                        )
                        break
                if _PLOT_RESULT is None:
                    _PLOT_RESULT = (
                        True,
                        f"Passed (cycle-10 {len(rst_plot['fields'])} real8 fields match; "
                        f"max rel diff {max_field_err:.2e} in {worst_var})",
                    )

        return (
            True,
            f"Passed (cycles {overlap_cycles[0]}..{overlap_cycles[-1]} log energies match "
            f"uninterrupted run, nprocs={nprocs}, max rel diff={max_log_err:.2e})",
        )
    finally:
        with open(param_path, "w") as f:
            f.write(orig_param)
        if os.path.exists(restart_in):
            shutil.rmtree(restart_in)
        if os.path.exists(plots_backup):
            shutil.rmtree(plots_backup)


def validate_plot(test_name=None):
    """Return cached cycle-10 real8 plot comparison from validate_log()."""
    if _PLOT_RESULT is not None:
        return _PLOT_RESULT
    return False, "Plot comparison was not completed during validate_log()"
