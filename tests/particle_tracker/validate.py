#!/usr/bin/env python3
"""Validator for the pure test-particle tracker tests (tests/particle_tracker/).

Four variants are discovered from this directory:
  - PARAM.in      -> "particle_tracker"      Relativistic Boris E x B benchmark
  - PARAM.in.vay  -> "particle_tracker_vay"  Relativistic Vay (2008) benchmark
  - PARAM.in.hc   -> "particle_tracker_hc"   Relativistic Higuera-Cary (2017) + multi-flush TPOUTPUT
  - PARAM.in.amr  -> "particle_tracker_amr"  2-level AMR leaf seeding, reflecting walls, body, dipole
"""
import glob
import logging
import math
import os
import struct

logger = logging.getLogger(__name__)

RUN_DIR = "run_test"

C_LIGHT = 2.99792458e8
M_PROTON = 1.67262192369e-27
BETA_E = 0.8
V_E = BETA_E * C_LIGHT
GAMMA_E = 1.0 / math.sqrt(1.0 - BETA_E * BETA_E)  # 5/3

# Do not invoke the runner's generic PIC-oriented launchThreshold log checker,
# since these runs intentionally set launchThreshold = 0.0 and check full
# double-precision moments directly inside validate_log().
PARTICLE_TOL = None


def set_run_dir(run_dir):
    global RUN_DIR
    RUN_DIR = run_dir


def _read_pt_log_rows(run_dir):
    files = sorted(glob.glob(os.path.join(run_dir, "PC", "plots", "log_pt_n*.log")))
    if not files:
        return []
    rows = []
    with open(files[-1], "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if len(lines) < 2:
        return []
    header = lines[0].split()
    n_species = sum(1 for col in header if col.startswith("mass_"))
    for ln in lines[1:]:
        vals = ln.split()
        if len(vals) < 2 + 5 * n_species:
            continue
        row = {"time": float(vals[0]), "cycle": int(vals[1])}
        idx = 2
        for s in range(n_species):
            row[f"mass{s}"] = float(vals[idx])
            row[f"moment_x{s}"] = float(vals[idx + 1])
            row[f"moment_y{s}"] = float(vals[idx + 2])
            row[f"moment_z{s}"] = float(vals[idx + 3])
            row[f"energy{s}"] = float(vals[idx + 4])
            idx += 5
        rows.append(row)
    return rows


def _read_tp_header(tp_dir):
    header_path = os.path.join(tp_dir, "Header")
    if not os.path.isfile(header_path):
        raise FileNotFoundError(f"Missing Header in {tp_dir}")
    with open(header_path, "r") as f:
        lines = [ln.strip() for ln in f.readlines()]
    if len(lines) < 2:
        raise ValueError(f"Expected >= 2 lines in {header_path}, got {len(lines)}")
    n_real = int(lines[0])
    cols = lines[1].split()
    if len(cols) != n_real:
        raise ValueError(
            f"Header column count mismatch: line 1 has {n_real}, line 2 has {len(cols)} ({cols})"
        )
    return n_real, cols


def _read_species_trajectories(tp_dir, i_species, n_real):
    """Read and validate binary trajectory & index list files for one species.

    Returns (data_files, combined_trajectories_by_pid), where each PID maps to a
    sorted list of records ordered across flush windows.
    """
    data_files = sorted(
        glob.glob(
            os.path.join(tp_dir, f"FLEKS1_particle_species_{i_species}_n*")
        )
    )
    if not data_files:
        raise AssertionError(f"No trajectory data files for species {i_species} in {tp_dir}")

    combined = {}
    list_struct = struct.Struct("iiQ")
    hdr_struct = struct.Struct("iiif")

    for df in data_files:
        lf = df.replace("_particle_species_", "_particle_list_species_")
        if not os.path.isfile(lf):
            raise AssertionError(f"Missing list file {lf} matching {df}")

        list_bytes = os.path.getsize(lf)
        if list_bytes % list_struct.size != 0:
            raise AssertionError(f"Corrupt list file size {list_bytes} in {lf}")
        n_list = list_bytes // list_struct.size
        list_map = {}
        with open(lf, "rb") as f:
            for _ in range(n_list):
                cpu, pid, offset = list_struct.unpack(f.read(list_struct.size))
                key = (cpu, pid)
                if key in list_map:
                    raise AssertionError(f"Duplicate PID {key} in list file {lf}")
                list_map[key] = offset

        seen_in_file = set()
        file_size = os.path.getsize(df)
        with open(df, "rb") as f:
            while f.tell() < file_size:
                rec_offset = f.tell()
                chunk = f.read(hdr_struct.size)
                if len(chunk) != hdr_struct.size:
                    raise AssertionError(f"Truncated record header at {rec_offset} in {df}")
                cpu, pid, n_rec, weight = hdr_struct.unpack(chunk)
                key = (cpu, pid)
                if key in seen_in_file:
                    raise AssertionError(f"Duplicate PID {key} in data file {df}")
                seen_in_file.add(key)

                if key not in list_map:
                    raise AssertionError(f"PID {key} missing from list file {lf}")
                if list_map[key] != rec_offset:
                    raise AssertionError(
                        f"List offset mismatch for {key} in {lf}: "
                        f"{list_map[key]} != actual {rec_offset}"
                    )

                payload_size = 4 * n_real * n_rec
                raw = f.read(payload_size)
                if len(raw) != payload_size:
                    raise AssertionError(f"Truncated trajectory payload for {key} in {df}")
                floats = struct.unpack(f"{n_real * n_rec}f", raw)
                recs = [
                    floats[i * n_real : (i + 1) * n_real] for i in range(n_rec)
                ]
                combined.setdefault(key, []).extend(recs)

        if seen_in_file != set(list_map.keys()):
            raise AssertionError(f"List/data key set mismatch in {df}")

    for key, recs in combined.items():
        if not recs:
            raise AssertionError(f"Empty trajectory for {key}")
        for i in range(len(recs) - 1):
            if recs[i + 1][0] <= recs[i][0]:
                raise AssertionError(
                    f"Non-increasing timestamps for {key}: {recs[i][0]} -> {recs[i+1][0]}"
                )

    return data_files, combined


def validate_log(pic_diags=None, test_name="particle_tracker"):
    """Verify PIC log suppression (#PIC F) and double-precision PT moments."""
    if pic_diags:
        return False, f"Expected no log_pic_n*.log when PIC=F, found {len(pic_diags)} rows"

    pic_logs = glob.glob(os.path.join(RUN_DIR, "PC", "plots", "log_pic_n*.log"))
    if pic_logs:
        return False, f"Unexpected log_pic file created with PIC=F: {pic_logs}"

    pt_rows = _read_pt_log_rows(RUN_DIR)
    if len(pt_rows) < 18:
        return False, f"Too few log_pt rows: {len(pt_rows)}"

    for i in range(len(pt_rows) - 1):
        if pt_rows[i + 1]["time"] < pt_rows[i]["time"]:
            return False, "PT log time went backwards"
        if pt_rows[i + 1]["cycle"] <= pt_rows[i]["cycle"]:
            return False, "PT log cycle did not strictly increase"

    if test_name in ("particle_tracker", "particle_tracker_vay", "particle_tracker_hc"):
        vx_vals = []
        for r in pt_rows:
            m0 = r["mass0"]
            if not math.isfinite(m0) or m0 <= 0.0:
                return False, f"Non-positive mass0={m0} in PT log"
            vx_vals.append(r["moment_x0"] / m0)

        max_err = max(abs(vx - BETA_E) / BETA_E for vx in vx_vals)
        if test_name == "particle_tracker":
            if max_err <= 0.02:
                return (
                    False,
                    f"Expected relativistic Boris drift error > 2%, got {max_err:.3e}",
                )
        else:
            if max_err >= 1e-10:
                return (
                    False,
                    f"{test_name} E x B drift error {max_err:.3e} >= 1e-10",
                )

    return True, "Passed"


def _check_flekspy_load(tp_dir, i_species, expected_n_real):
    try:
        from flekspy import FLEKSTP
    except ImportError:
        logger.debug("  flekspy not installed; skipping FLEKSTP read check.")
        return

    tp = FLEKSTP(tp_dir, iDomain=1, iSpecies=i_species, unit="SI")
    if len(tp) == 0:
        raise AssertionError("flekspy.FLEKSTP loaded 0 particles")
    traj0 = tp[0].collect().to_numpy()
    if traj0.ndim != 2 or traj0.shape[1] != expected_n_real or traj0.shape[0] < 2:
        raise AssertionError(
            f"Unexpected flekspy trajectory shape {traj0.shape}, expected (*, {expected_n_real})"
        )


def validate_plot(test_name="particle_tracker"):
    """Validate binary particle trajectory files, Header columns, and pusher physics."""
    tp_dir = os.path.join(RUN_DIR, "PC", "plots", "test_particles")
    if not os.path.isdir(tp_dir):
        return False, f"Missing test_particles output directory: {tp_dir}"

    try:
        n_real, cols = _read_tp_header(tp_dir)
    except Exception as exc:
        return False, str(exc)

    expected_cols_13 = [
        "time", "x", "y", "z", "vx", "vy", "vz", "bx", "by", "bz", "ex", "ey", "ez"
    ]
    expected_cols_18 = expected_cols_13 + [
        "v_exb_x", "v_exb_y", "v_exb_z", "energy", "gamma"
    ]
    expected_cols_26 = expected_cols_13 + [
        "v_gradb_x", "v_gradb_y", "v_gradb_z",
        "v_curv_x", "v_curv_y", "v_curv_z",
        "v_exb_x", "v_exb_y", "v_exb_z",
        "mu", "pitch", "energy", "gamma",
    ]

    try:
        if test_name in ("particle_tracker", "particle_tracker_vay"):
            if n_real != 13 or cols != expected_cols_13:
                return False, f"Unexpected Header for {test_name}: nReal={n_real}, cols={cols}"

            files0, tr0 = _read_species_trajectories(tp_dir, 0, n_real)
            files1, tr1 = _read_species_trajectories(tp_dir, 1, n_real)
            if len(files0) != 2 or len(files1) != 2:
                return False, f"Expected 2 flush files per species, got {len(files0)}, {len(files1)}"

            # Species 0: force-free E x B drift
            max_vx_err = 0.0
            max_vy_err = 0.0
            for recs in tr0.values():
                if len(recs) < 20:
                    return False, f"Expected >= 20 records per particle, got {len(recs)}"
                for rec in recs:
                    vx, vy, vz = rec[4], rec[5], rec[6]
                    speed = math.sqrt(vx * vx + vy * vy + vz * vz)
                    if speed >= C_LIGHT:
                        return False, f"Superluminal speed {speed} >= {C_LIGHT}"
                    max_vx_err = max(max_vx_err, abs(vx - V_E) / V_E)
                    max_vy_err = max(max_vy_err, abs(vy) / C_LIGHT)

            if test_name == "particle_tracker":
                if max_vx_err <= 0.02:
                    return False, f"Expected Boris vx error > 2%, got {max_vx_err:.3e}"
            else:
                if max_vx_err >= 1e-5 or max_vy_err >= 1e-5:
                    return (
                        False,
                        f"Vay E x B drift error too large: vx_err={max_vx_err:.3e}, vy_err={max_vy_err:.3e}",
                    )

            # Species 1: relativistic cycloid from rest
            max_gamma = 0.0
            max_inv_dev = 0.0
            for recs in tr1.values():
                for rec in recs:
                    vx, vy, vz = rec[4], rec[5], rec[6]
                    v2 = vx * vx + vy * vy + vz * vz
                    if v2 >= C_LIGHT * C_LIGHT:
                        return False, f"Species 1 superluminal velocity v={math.sqrt(v2)}"
                    gamma = 1.0 / math.sqrt(1.0 - v2 / (C_LIGHT * C_LIGHT))
                    max_gamma = max(max_gamma, gamma)
                    inv = gamma * (1.0 - BETA_E * (vx / C_LIGHT))
                    max_inv_dev = max(max_inv_dev, abs(inv - 1.0))

            if test_name == "particle_tracker_vay":
                if max_gamma <= 3.5:
                    return False, f"Expected Vay peak gamma > 3.5, got {max_gamma:.3f}"
                if max_inv_dev >= 0.35:
                    return False, f"Vay boosted-frame invariant drift {max_inv_dev:.3e} >= 0.35"

            _check_flekspy_load(tp_dir, 0, 13)
            return True, (
                f"Passed (nPart={len(tr0)}+{len(tr1)}, max_vx_err={max_vx_err:.2e}, "
                f"max_gamma={max_gamma:.2f})"
            )

        if test_name == "particle_tracker_hc":
            if n_real != 18 or cols != expected_cols_18:
                return False, f"Unexpected Header for HC: nReal={n_real}, cols={cols}"

            files0, tr0 = _read_species_trajectories(tp_dir, 0, n_real)
            files1, tr1 = _read_species_trajectories(tp_dir, 1, n_real)
            # n00000000 + n00000008 + n00000016 + n00000024 = 4 files
            if len(files0) != 4 or len(files1) != 4:
                return False, f"Expected 4 files per species (init + 3 flushes), got {len(files0)}, {len(files1)}"

            max_vx_err = 0.0
            max_vy_err = 0.0
            max_vexb_err = 0.0
            max_gam_err = 0.0
            for recs in tr0.values():
                if len(recs) != 25:
                    return False, f"Expected 25 records (t=0 + 24 steps), got {len(recs)}"
                for rec in recs:
                    vx, vy, vz = rec[4], rec[5], rec[6]
                    vexbx = rec[13]
                    gamma_out = rec[17]
                    max_vx_err = max(max_vx_err, abs(vx - V_E) / V_E)
                    max_vy_err = max(max_vy_err, abs(vy) / C_LIGHT)
                    max_vexb_err = max(max_vexb_err, abs(vexbx - V_E) / V_E)
                    max_gam_err = max(max_gam_err, abs(gamma_out - GAMMA_E))

            if max_vx_err >= 1e-5 or max_vy_err >= 1e-5:
                return False, f"HC E x B error too large: vx={max_vx_err:.3e}, vy={max_vy_err:.3e}"
            if max_vexb_err >= 1e-5 or max_gam_err >= 1e-5:
                return False, f"HC output diagnostic error: vexb={max_vexb_err:.3e}, gamma={max_gam_err:.3e}"

            max_energy_rel_err = 0.0
            mc2 = M_PROTON * C_LIGHT * C_LIGHT
            for recs in tr1.values():
                for rec in recs:
                    vx, vy, vz = rec[4], rec[5], rec[6]
                    v2 = vx * vx + vy * vy + vz * vz
                    if v2 >= C_LIGHT * C_LIGHT:
                        return False, f"HC Species 1 superluminal velocity v={math.sqrt(v2)}"
                    gamma_exact = 1.0 / math.sqrt(1.0 - v2 / (C_LIGHT * C_LIGHT))
                    e_out, g_out = rec[16], rec[17]
                    max_energy_rel_err = max(
                        max_energy_rel_err,
                        abs(e_out - (gamma_exact - 1.0) * mc2) / mc2,
                        abs(g_out - gamma_exact) / gamma_exact,
                    )

            if max_energy_rel_err >= 1e-5:
                return False, f"HC Species 1 energy/gamma mismatch: {max_energy_rel_err:.3e}"

            return True, (
                f"Passed (4 flush files, max_vx_err={max_vx_err:.2e}, "
                f"max_energy_rel_err={max_energy_rel_err:.2e})"
            )

        if test_name == "particle_tracker_amr":
            if n_real != 26 or cols != expected_cols_26:
                return False, f"Unexpected Header for AMR: nReal={n_real}, cols={cols}"

            files0, tr0 = _read_species_trajectories(tp_dir, 0, n_real)
            # n00000000 + n00000006 + n00000012 + n00000018 = 4 files
            if len(files0) != 4:
                return False, f"Expected 4 files (init + 3 flushes), got {len(files0)}"

            n_unique = len(tr0)
            # L0 has 16*16*4 = 1024 cells; L1 replaces central 128 L0 cells with 1024 L1 cells
            # minus interior body cells (r < 1.0).
            if not (1024 < n_unique < 1920):
                return (
                    False,
                    f"AMR leaf seeding particle count {n_unique} outside expected range (1024, 1920)",
                )

            max_energy_drift = 0.0
            max_vexb = 0.0
            n_multi = 0
            for recs in tr0.values():
                x0, y0, z0 = recs[0][1], recs[0][2], recs[0][3]
                r0 = math.sqrt(x0 * x0 + y0 * y0 + z0 * z0)
                if r0 < 0.99:
                    return False, f"Particle seeded inside body interior at r0={r0:.4f}"

                e0 = recs[0][24]
                if len(recs) >= 2:
                    n_multi += 1
                for rec in recs:
                    x = rec[1]
                    if x < -4.0 - 1e-4 or x > 4.0 + 1e-4:
                        return False, f"Particle escaped reflecting X boundary: x={x}"

                    vexb = math.sqrt(rec[19] ** 2 + rec[20] ** 2 + rec[21] ** 2)
                    max_vexb = max(max_vexb, vexb)

                    mu = rec[22]
                    pitch = rec[23]
                    energy = rec[24]
                    gamma = rec[25]
                    if not (math.isfinite(mu) and mu >= 0.0):
                        return False, f"Invalid mu={mu}"
                    if not (-1.0001 <= pitch <= 1.0001):
                        return False, f"Invalid pitch cos(alpha)={pitch}"
                    if not (math.isfinite(gamma) and gamma >= 1.0):
                        return False, f"Invalid gamma={gamma}"
                    if e0 > 0.0:
                        max_energy_drift = max(max_energy_drift, abs(energy - e0) / e0)

            if n_multi == 0:
                return False, "No multi-point trajectories recorded in AMR test"
            if max_energy_drift >= 1e-5:
                return False, f"AMR E=0 energy conservation error {max_energy_drift:.3e} >= 1e-5"
            if max_vexb >= 1e-6:
                return False, f"Nonzero v_exb={max_vexb:.3e} when E=0"

            return True, (
                f"Passed (nUnique={n_unique}, energy_drift={max_energy_drift:.2e}, "
                f"max_vexb={max_vexb:.2e})"
            )

        return False, f"Unknown particle_tracker variant '{test_name}'"
    except Exception as exc:
        return False, str(exc)
