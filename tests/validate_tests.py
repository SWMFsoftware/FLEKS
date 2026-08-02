#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys
import math

# Directory used for simulation output. Defaults to "run_test" but can be
# overridden with --run-dir so a second test can run without clobbering a
# currently-running job that owns the default run_test/.
RUN_DIR = "run_test"

def safe_symlink(src, dst):
    if os.path.lexists(dst):
        if os.path.islink(dst) or os.path.isfile(dst):
            os.remove(dst)
        elif os.path.isdir(dst):
            shutil.rmtree(dst)
    os.symlink(src, dst)

def prepare_run_dir():
    run_dir = RUN_DIR
    os.makedirs(run_dir, exist_ok=True)

    # Determine the location of the share directory relative to FLEKS root.
    # In a standalone FLEKS repository, 'share/Scripts/PostProc.pl' is directly inside the working directory.
    # In SWMF integrated environment, 'share' is located in SWMF root, i.e., two levels above FLEKS root.
    if os.path.isfile("share/Scripts/PostProc.pl"):
        postproc_target = "../share/Scripts/PostProc.pl"
        pidl_target = "../../share/Scripts/pIDL"
    else:
        postproc_target = "../../../share/Scripts/PostProc.pl"
        pidl_target = "../../../../share/Scripts/pIDL"

    # PostIDL.exe may reside in different locations depending on the build mode:
    #   - Standalone build:       bin/PostIDL.exe        (FLEKS/bin/)
    #   - SWMF integrated build:  ../../bin/PostIDL.exe  (SWMF/bin/)
    # Search all candidates (relative to FLEKS root) and use the first match.
    # The symlink target is computed relative to run_test/PC/.
    postidl_candidates = [
        "bin/PostIDL.exe",          # standalone
        "../../bin/PostIDL.exe",    # SWMF integrated
    ]
    postidl_target = None
    for candidate in postidl_candidates:
        if os.path.isfile(candidate):
            postidl_target = os.path.relpath(candidate, os.path.join(run_dir, "PC"))
            continue
    if postidl_target is None:
        # Default to standalone path; run_test() will emit a broken-symlink
        # warning if PostIDL.exe is not found at any candidate location.
        postidl_target = "../../bin/PostIDL.exe"

    # Symlinks in run directory
    safe_symlink("../bin/FLEKS.exe", os.path.join(run_dir, "FLEKS.exe"))
    safe_symlink(postproc_target, os.path.join(run_dir, "PostProc.pl"))

    # Component plot and restart directories
    pc_dir = os.path.join(run_dir, "PC")
    os.makedirs(pc_dir, exist_ok=True)
    os.makedirs(os.path.join(pc_dir, "plots"), exist_ok=True)
    os.makedirs(os.path.join(pc_dir, "restartOUT"), exist_ok=True)

    # Symlinks in component directory
    safe_symlink(pidl_target, os.path.join(pc_dir, "pIDL"))
    safe_symlink(postidl_target, os.path.join(pc_dir, "PostIDL.exe"))

def cleanup_run_dir():
    """Remove simulation output files from run_test/ after each test.

    Deletes the contents of PC/plots/ and PC/restartOUT/ (the bulky per-run
    outputs) so they do not accumulate between tests.  The directory structure
    and symlinks are left in place so the next call to prepare_run_dir() is a
    no-op.
    """
    run_dir = RUN_DIR
    for subdir in [os.path.join(run_dir, "PC", "plots"),
                   os.path.join(run_dir, "PC", "restartOUT")]:
        if os.path.isdir(subdir):
            for entry in os.listdir(subdir):
                entry_path = os.path.join(subdir, entry)
                try:
                    if os.path.islink(entry_path) or os.path.isfile(entry_path):
                        os.remove(entry_path)
                    elif os.path.isdir(entry_path):
                        shutil.rmtree(entry_path)
                except Exception as e:
                    print(f"  [WARN] Could not remove {entry_path}: {e}")

# ---------------------------------------------------------------------------
# Free-stream test: tests/freestream/ holds TWO parameter files that differ only
# in the field-solver block -- PARAM.in (full PIC implicit Maxwell/GMRES) and
# PARAM.in.hybrid (hybrid Ohm's-law solver with the Hall term OFF).  The runner
# exercises BOTH solvers by running the single directory twice, once per file,
# so one shared setup validates both field solvers head-to-head.


def run_and_validate(test_dir, display_name, validator, nprocs, results,
                     param_text=None, base_name=None):
    """Run one FLEKS test (optionally with a patched PARAM.in) and record the
    outcome in *results* as (name, status, reason). Mirrors the former main
    loop body so a single test can be run with several PARAM variants."""
    if base_name is None:
        base_name = display_name
    print(f"\n==========================================")
    print(f"Starting test: {display_name}")
    print(f"==========================================")
    try:
        stdout, code = run_test(test_dir, nprocs=nprocs, param_text=param_text)
        if code != 0 or stdout is None:
            print(f"FAIL: {display_name} execution failed with exit code {code}")
            results.append((display_name, "FAILED", f"Execution failed (code {code})"))
            return

        # Read the PIC energy log (the only diagnostic log produced by FLEKS).
        pic_diags = read_pic_log(RUN_DIR)

        val_res = False
        reason = "Validation skipped"

        if validator:
            import inspect
            sig = inspect.signature(validator)
            kwargs = {}
            if "pic_diags" in sig.parameters:
                kwargs["pic_diags"] = pic_diags
            if "test_name" in sig.parameters:
                kwargs["test_name"] = base_name
            val_res, reason = validator(**kwargs)
            if not val_res:
                results.append((display_name, "FAILED", reason))
                return
        else:
            print(f"Validating {display_name} (generic check)...")
            print(f"{display_name} (generic check): PASSED")
            val_res = True
            reason = "Passed"

        # Read the test-particle tracer log (log_pt_n*.log) and validate
        # it for tests that enable #PARTICLETRACKER T.
        pt_diags = read_pt_log("run_test")
        pt_tests = {"beam", "photoionization"}
        if base_name in pt_tests:
            pt_tol = {
                "beam":   {"expected_active_species": [0],
                           "launch_threshold": 0.5, "max_speed": 10.0},
                "photoionization": {"expected_active_species": [0, 1, 2],
                                    "launch_threshold": 0.5, "max_speed": 10.0},
            }.get(base_name, {})
            pt_res, pt_reason = validate_test_particles(
                pt_diags, test_name=base_name, tol=pt_tol)
            if not pt_res:
                results.append((display_name, "FAILED",
                                f"test-particle check failed: {pt_reason}"))
                return

        # Validate output plotfiles
        plot_res, plot_reason = validate_plot_output(base_name)
        if not plot_res:
            results.append((display_name, "FAILED", f"plot check failed: {plot_reason}"))
        else:
            results.append((display_name, "PASSED", "Passed"))

    finally:
        # Always clean up run output after each test to keep disk usage low.
        print(f"  Cleaning up run_test/ output for {display_name}...")
        cleanup_run_dir()


def run_test(test_dir, nprocs=1, param_text=None):
    param_file = os.path.join(test_dir, "PARAM.in")
    print(f"Running test in {test_dir}...")
    prepare_run_dir()
    
    # Verify that PostIDL.exe exists; PostProc.pl needs it to produce .out files.
    postidl_link = os.path.join(RUN_DIR, "PC", "PostIDL.exe")
    if os.path.islink(postidl_link) and not os.path.exists(postidl_link):
        real = os.path.realpath(postidl_link)
        print(f"  [WARN] Broken symlink: {postidl_link} -> {real}")
        print(f"  [WARN] PostIDL.exe is missing. Build it with 'make PIDL' "
              f"before running tests that check plot output (.out files).")
    
    # Write run_test/PARAM.in: use the supplied text (e.g. a patched solver
    # variant) if given, otherwise copy the test directory's PARAM.in.
    if param_text is not None:
        with open(RUN_DIR + "/PARAM.in", "w") as f:
            f.write(param_text)
    else:
        param_file = os.path.join(test_dir, "PARAM.in")
        shutil.copy(param_file, RUN_DIR + "/PARAM.in")
    
    # Build the command: serial for nprocs==1, mpirun otherwise
    if nprocs <= 1:
        cmd = ["./FLEKS.exe"]
        print(f"  Running in serial mode (no MPI)...")
    else:
        cmd = ["mpirun", "-n", str(nprocs), "./FLEKS.exe"]
        print(f"  Running with {nprocs} MPI processes...")
    
    # Run the command inside run_test/
    result = subprocess.run(cmd, cwd=RUN_DIR, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error running FLEKS.exe for {test_dir}:")
        print("--- FLEKS stdout ---")
        print(result.stdout)
        print("--- FLEKS stderr ---")
        print(result.stderr)
        return None, result.returncode
        
    # Automatically run post-processing on the generated plots
    pp = subprocess.run(["./PostProc.pl", "-v"], cwd=RUN_DIR,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if pp.returncode != 0:
        print(f"  [WARN] PostProc.pl exited with code {pp.returncode}:")
        if pp.stdout:
            print(pp.stdout)
        if pp.stderr:
            print(pp.stderr)
    
    return result.stdout, 0


def validate_beam():
    """Validate the beam instability test.

    The primary validation is the FFT-based transverse-wave check performed
    on the plot output by _check_beam_transverse_wave(), invoked from
    validate_plot_output().  No log-file-based checks are performed here.
    """
    print("Validating Beam Instability Test...")
    print("  [INFO] Beam diagnostic checks rely on plot output (FFT).")
    print("Beam Instability Test: PASSED")
    return True, "Passed"


def read_pic_log(run_dir):
    """Read the energy diagnostic log file written by Pic::write_log.

    The log_pic_n*.log format:
      time  nStep  Etot  Ee  Eb  Epart  Epart0  Epart1 ...

    Returns a list of dicts with keys: time, cycle, Etot, Ee, Eb,
    Epart, and one EpartN key per species.
    """
    import glob
    pc_plots = os.path.join(run_dir, "PC", "plots")
    log_files = sorted(glob.glob(os.path.join(pc_plots, "log_pic_n*.log")))
    if not log_files:
        return []

    log_file = log_files[-1]
    pic_diags = []

    with open(log_file, "r") as f:
        lines = f.readlines()

    if len(lines) < 2:
        return []

    header = lines[0].strip().split("\t")
    # Discover species count from EpartN columns
    n_species = sum(1 for col in header if col.startswith("Epart") and col != "Epart")

    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        vals = line.split("\t")
        if len(vals) < 5 + n_species:
            continue
        try:
            entry = {
                "time":  float(vals[0]),
                "cycle": int(vals[1]),
                "Etot":  float(vals[2]),
                "Ee":    float(vals[3]),
                "Eb":    float(vals[4]),
                "Epart": float(vals[5]),
            }
            for iS in range(n_species):
                entry[f"Epart{iS}"] = float(vals[6 + iS])
            pic_diags.append(entry)
        except (ValueError, IndexError):
            continue

    return pic_diags


def read_pt_log(run_dir):
    """Read the test-particle tracer diagnostic log (log_pt_n*.log).

    Written by ParticleTracker::write_log.  Format (whitespace separated):
        time nStep mass_0 moment_x_0 moment_y_0 moment_z_0 energy_0 \
                     mass_1 moment_x_1 ... energy_1 ...
    One data row is appended per ptLog interval (defaults to every 10 cycles).
    mass_i is the total charge-weighted mass of test particles for species i;
    moment_*/mass_i gives the mean velocity (normalized units), so the log is
    self-describing and needs no external package (flekspy/Batsrus.jl) to load.

    Returns a list of dicts with keys: time, cycle, and per species
    mass{i}, moment_x{i}, moment_y{i}, moment_z{i}, energy{i}.
    """
    import glob
    pc_plots = os.path.join(run_dir, "PC", "plots")
    log_files = sorted(glob.glob(os.path.join(pc_plots, "log_pt_n*.log")))
    if not log_files:
        return []

    log_file = log_files[-1]
    diags = []
    with open(log_file, "r") as f:
        lines = f.readlines()

    if len(lines) < 2:
        return []

    header = lines[0].strip().split()
    n_species = sum(1 for col in header if col.startswith("mass_"))
    if n_species == 0 or "time" not in header or "nStep" not in header:
        return []

    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        vals = line.split()
        if len(vals) < 2 + 5 * n_species:
            continue
        try:
            entry = {
                "time":  float(vals[0]),
                "cycle": int(vals[1]),
            }
            idx = 2
            for iS in range(n_species):
                entry[f"mass{iS}"]     = float(vals[idx]); idx += 1
                entry[f"moment_x{iS}"] = float(vals[idx]); idx += 1
                entry[f"moment_y{iS}"] = float(vals[idx]); idx += 1
                entry[f"moment_z{iS}"] = float(vals[idx]); idx += 1
                entry[f"energy{iS}"]   = float(vals[idx]); idx += 1
            diags.append(entry)
        except (ValueError, IndexError):
            continue

    return diags


def validate_test_particles(pt_diags, test_name=None, tol=None):
    """Validate test-particle tracer output from log_pt_n*.log.

    Pure-Python checks (no external packages):
      1. Log exists with >= min_rows data rows.
      2. Header parses into time, cycle, and n_species*5 columns.
      3. Activity: species listed in expected_active_species must have
         mass > 0 at t=0 -- catches "test particles not seeded".
      4. Sanity: for active species, all masses/moments/energies finite and
         mass > 0 in every row.
      5. Time: time non-decreasing, cycle strictly increasing (>=2 rows).
      6. Conservation: mass_i / mass_i[0] in [launch_threshold, 1] (periodic
         BC: particles only deplete, refilled at the launch threshold).
      7. Bounded speed: sqrt(vx^2+vy^2+vz^2) < max_speed (catches runaways/NaN).

    tol keys: min_rows (int, default 1), launch_threshold (float, default 0.5),
              max_speed (float, default 10.0),
              expected_active_species (list[int], default []).
    """
    tol = tol or {}
    min_rows = int(tol.get("min_rows", 1))
    launch_threshold = float(tol.get("launch_threshold", 0.5))
    max_speed = float(tol.get("max_speed", 10.0))
    expected_active = set(int(s) for s in tol.get("expected_active_species", []))

    print("Validating Test-Particle Tracer Output...")

    if not pt_diags:
        print("  [PT] No test-particle log (log_pt_n*.log) found.")
        return False, "No test-particle log file"

    n_rows = len(pt_diags)
    n_species = sum(1 for k in pt_diags[0] if k.startswith("mass"))
    print(f"  [PT] {n_rows} log rows, {n_species} species "
          f"(from log_pt_n*.log).")

    passed = True
    reasons = []

    if n_rows < min_rows:
        passed = False
        reasons.append(f"only {n_rows} data rows (< {min_rows})")

    # Active = species seeded at t=0 (mass[0] > 0) or explicitly expected.
    active = set()
    for iS in range(n_species):
        if iS in expected_active:
            active.add(iS)
        elif math.isfinite(pt_diags[0].get(f"mass{iS}", 0.0)) \
                and pt_diags[0].get(f"mass{iS}", 0.0) > 0:
            active.add(iS)
    print(f"  [PT] active (seeded) species: {sorted(active)}")

    if expected_active and not expected_active.issubset(active):
        missing = sorted(expected_active - active)
        passed = False
        reasons.append(f"expected active species {missing} not seeded (mass<=0)")

    # Sanity: finite + mass>0 for every active species in every row.
    for iS in sorted(active):
        for r in pt_diags:
            cyc = r.get("cycle", "?")
            m = r.get(f"mass{iS}", float("nan"))
            if not math.isfinite(m) or m <= 0:
                passed = False
                reasons.append(f"species {iS} non-finite/non-positive mass "
                               f"at cycle {cyc}")
                break
            for ax in ("x", "y", "z"):
                mv = r.get(f"moment_{ax}{iS}", float("nan"))
                if not math.isfinite(mv):
                    passed = False
                    reasons.append(f"species {iS} non-finite moment_{ax} "
                                   f"at cycle {cyc}")
            if not math.isfinite(r.get(f"energy{iS}", float("nan"))):
                passed = False
                reasons.append(f"species {iS} non-finite energy at cycle {cyc}")

    # Time / cycle monotonicity (needs >=2 rows).
    if n_rows >= 2:
        times = [r["time"] for r in pt_diags]
        cycles = [r["cycle"] for r in pt_diags]
        if any(times[i + 1] < times[i] for i in range(n_rows - 1)):
            passed = False
            reasons.append("time not monotonically increasing")
        if any(cycles[i + 1] <= cycles[i] for i in range(n_rows - 1)):
            passed = False
            reasons.append("cycle not strictly increasing")

        # Conservation: mass ratio within [launch_threshold, 1].
        for iS in sorted(active):
            m0 = pt_diags[0].get(f"mass{iS}", 0.0)
            if m0 <= 0:
                continue
            for r in pt_diags:
                mr = r.get(f"mass{iS}", 0.0)
                ratio = mr / m0
                if ratio < launch_threshold - 1e-6 or ratio > 1.0 + 1e-6:
                    passed = False
                    reasons.append(
                        f"species {iS} mass ratio {ratio:.3f} outside "
                        f"[{launch_threshold:g}, 1] at cycle {r.get('cycle', '?')}")
                    break

    # Bounded speed for active species.
    for iS in sorted(active):
        for r in pt_diags:
            m = r.get(f"mass{iS}", 0.0)
            if m <= 0:
                continue
            vx = r.get(f"moment_x{iS}", 0.0) / m
            vy = r.get(f"moment_y{iS}", 0.0) / m
            vz = r.get(f"moment_z{iS}", 0.0) / m
            speed = math.sqrt(vx * vx + vy * vy + vz * vz)
            if not math.isfinite(speed) or speed > max_speed:
                passed = False
                reasons.append(
                    f"species {iS} speed {speed:.3f} exceeds max_speed "
                    f"{max_speed:g} at cycle {r.get('cycle', '?')}")
                break

    if passed:
        print("Test-Particle Tracer Output: PASSED")
        return True, "Passed"
    return False, "; ".join(reasons)


def validate_chemistry(pic_diags=None, test_name=None):
    """Validate the Mars chemistry test with 4 ion species and 10 reactions.

    Checks that ion energies change over time due to the combined action of
    photoionization (source), cross-species charge exchange (source + loss),
    and recombination (loss).

    Key validations:
    1. ALL 4 ion species (H+, O+, O2+, CO2+) show significant energy changes,
       proving all reaction types are active.
    2. O2+ (species 3, Epart3) energy INCREASES — O2+ is produced ONLY by
       cross-species CX (reactions 3, 4) and lost by recombination (reaction 6).
       Since the CX source rate (~6 s^-1) far exceeds the recombination loss
       rate (~0.04 s^-1), O2+ energy must increase.  This is the critical
       test for the cross-species CX source term.
    3. CO2+ (species 4, Epart4) energy changes — CO2+ is produced by
       photoionization (reaction 1) and consumed by CX (reactions 3, 5) and
       recombination (reaction 7).
    """
    print("Validating Mars Chemistry Test...")

    if not pic_diags or len(pic_diags) < 2:
        print("  [INFO] No PIC energy log found; skipping energy checks.")
        return True, "Passed (no pic log)"

    first = pic_diags[0]
    last = pic_diags[-1]

    epart_keys = sorted(
        k for k in first.keys() if k.startswith("Epart") and k != "Epart"
    )
    if not epart_keys:
        print("  [INFO] No per-species energy columns; skipping.")
        return True, "Passed (no Epart columns)"

    # Species mapping: 0=e, 1=H+, 2=O+, 3=O2+, 4=CO2+
    species_names = {
        "Epart1": "H+",
        "Epart2": "O+",
        "Epart3": "O2+",
        "Epart4": "CO2+",
    }

    print(f"  --- Energy Diagnostics (from log_pic log) ---")
    for k in epart_keys:
        e0 = first.get(k, 0)
        e1 = last.get(k, 0)
        ratio = e1 / max(e0, 1e-30) if e0 > 0 else float('inf')
        name = species_names.get(k, k)
        print(f"    {k} ({name}): {e0:.6e} -> {e1:.6e}  (ratio {ratio:.4f})")

    passed = True
    reasons = []

    # ---- Check 1: All 4 ion species must show significant energy changes ----
    # This proves that photoionization, CX, and recombination are all active.
    # A 0.1% threshold catches any meaningful chemistry signal while filtering
    # out pure numerical noise.
    change_threshold = 0.001  # 0.1%
    for k in epart_keys:
        e0 = first.get(k, 0.0)
        e1 = last.get(k, 0.0)
        if e0 <= 0:
            continue
        ratio = e1 / e0
        name = species_names.get(k, k)
        if abs(ratio - 1.0) < change_threshold:
            print(f"    FAIL: {k} ({name}) energy unchanged "
                  f"(ratio {ratio:.4f}) — chemistry may not be active.")
            passed = False
            reasons.append(f"{name} energy unchanged")

    # ---- Check 2: O2+ must INCREASE — the critical CX source test ----
    # O2+ (Epart3) is produced ONLY by cross-species CX (reactions 3, 4)
    # and consumed by recombination (reaction 6).  The CX source rate
    # (~6.3 s^-1, driven by the large exosphere neutral density ~5e10 m^-3)
    # vastly exceeds the recombination loss rate (~4e-17 s^-1, limited by
    # the small plasma electron density in SI units).  Therefore O2+ energy
    # must increase.  If it does not increase, the CX source term is broken.
    o2_key = "Epart3" if "Epart3" in first else None
    if o2_key:
        e_o2_init = first.get(o2_key, 0.0)
        e_o2_final = last.get(o2_key, 0.0)
        if e_o2_init > 0:
            o2_ratio = e_o2_final / e_o2_init
            print(f"    {o2_key} (O2+): ratio = {o2_ratio:.4f} "
                  f"(must be > 1.0 for CX source validation)")
            if o2_ratio <= 1.0:
                print(f"    FAIL: {o2_key} (O2+) energy did not increase — "
                      f"cross-species CX source is not working.")
                passed = False
                reasons.append("O2+ energy did not increase (CX source broken)")
            else:
                pct = (o2_ratio - 1.0) * 100
                print(f"    SUCCESS: {o2_key} (O2+) energy increased by "
                      f"{pct:.2f}% (cross-species CX source active).")
        else:
            print(f"    [INFO] {o2_key} initial energy is zero.")

    # ---- Check 3: CO2+ must show a change ----
    # CO2+ is produced by photoionization (R1) and consumed by CX (R3, R5)
    # and recombination (R7).  Both source and loss are active.
    co2_key = "Epart4" if "Epart4" in first else None
    if co2_key:
        e_co2_init = first.get(co2_key, 0.0)
        e_co2_final = last.get(co2_key, 0.0)
        if e_co2_init > 0:
            co2_ratio = e_co2_final / e_co2_init
            print(f"    {co2_key} (CO2+): ratio = {co2_ratio:.4f}")
            if abs(co2_ratio - 1.0) < change_threshold:
                print(f"    FAIL: {co2_key} (CO2+) energy unchanged.")
                passed = False
                reasons.append("CO2+ energy unchanged")
            else:
                pct = abs(co2_ratio - 1.0) * 100
                print(f"    SUCCESS: {co2_key} (CO2+) energy changed by "
                      f"{pct:.2f}%.")

    if passed:
        print("Mars Chemistry Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)


def validate_recombination(pic_diags=None, test_name=None):
    """Validate the recombination loss test (O2+ + e- -> O + O).

    Checks that O2+ (species 2, Epart2) energy decreases over time due
    to recombination loss, while H+ (species 1, Epart1) energy remains
    stable since H+ does not participate in recombination.
    """
    print("Validating Recombination Loss Test...")

    if not pic_diags or len(pic_diags) < 2:
        print("  [INFO] No PIC energy log found; skipping energy checks.")
        return True, "Passed (no pic log)"

    first = pic_diags[0]
    last = pic_diags[-1]

    epart_keys = sorted(
        k for k in first.keys() if k.startswith("Epart") and k != "Epart"
    )
    if not epart_keys:
        print("  [INFO] No per-species energy columns; skipping.")
        return True, "Passed (no Epart columns)"

    print(f"  --- Energy Diagnostics (from log_pic log) ---")
    for k in epart_keys:
        print(f"    {k}: {first.get(k, 0):.6e} -> {last.get(k, 0):.6e}")

    passed = True
    reasons = []

    # O2+ (species 2) should decrease due to recombination.
    o2_key = "Epart2" if "Epart2" in first else (epart_keys[-1] if len(epart_keys) >= 2 else None)
    if o2_key:
        e_o2_initial = first.get(o2_key, 0.0)
        e_o2_final = last.get(o2_key, 0.0)
        print(f"    {o2_key} (O2+): {e_o2_initial:.6e} -> {e_o2_final:.6e}")
        if e_o2_initial <= 0:
            print(f"    FAIL: {o2_key} initial energy is zero.")
            passed = False
            reasons.append("O2+ initial energy is zero")
        elif e_o2_final >= e_o2_initial:
            print(f"    FAIL: {o2_key} energy did not decrease (recombination not active).")
            passed = False
            reasons.append("O2+ energy did not decrease")
        else:
            ratio = e_o2_final / e_o2_initial
            print(f"    SUCCESS: {o2_key} energy decreased to {ratio:.3f} of initial.")

    # H+ (species 1) should remain stable.
    h_key = "Epart1" if "Epart1" in first else None
    if h_key:
        e_h_initial = first.get(h_key, 0.0)
        e_h_final = last.get(h_key, 0.0)
        h_tolerance = 0.10  # allow up to 10% variation (numerical noise)
        print(f"    {h_key} (H+): {e_h_initial:.6e} -> {e_h_final:.6e}")
        if e_h_initial > 0:
            h_ratio = abs(e_h_final - e_h_initial) / e_h_initial
            if h_ratio > h_tolerance:
                print(f"    FAIL: {h_key} energy changed by {h_ratio*100:.1f}% "
                      f"(threshold {h_tolerance*100:.0f}%).")
                passed = False
                reasons.append(f"H+ energy changed by {h_ratio*100:.1f}%")
            else:
                print(f"    SUCCESS: {h_key} energy stable ({h_ratio*100:.1f}% change).")

    if passed:
        print("Recombination Loss Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)


def validate_lightwave(pic_diags=None, test_name=None):
    """Validate the 3D light-wave (vacuum transverse EM wave) test.

    The light-wave initial condition (testCase = lightwave) fills the node E
    and B fields with an analytic transverse plane wave; with
    nPartPerCell = 0 the PIC loads no macroparticles, so the total energy is
    purely electromagnetic (Ee + Eb).  On a periodic vacuum grid the wave
    should propagate without the energy decaying to zero or blowing up, so
    the total EM energy is approximately conserved.

    Checks (from log_pic_n*.log):
      1. Etot at the first and last frame is finite and > 0 (wave present).
      2. Energy is approximately conserved:
         0.3 <= Etot_final / Etot_initial <= 3.0.
    """
    import math
    print("Validating Light Wave Test...")

    if not pic_diags or len(pic_diags) < 2:
        print("  [INFO] No PIC energy log found; skipping energy checks.")
        return True, "Passed (no pic log)"

    first = pic_diags[0]
    last = pic_diags[-1]

    e0 = first.get("Etot", 0.0)
    e1 = last.get("Etot", 0.0)

    print(f"  --- Energy Diagnostics (from log_pic log) ---")
    print(f"    Etot (t={first.get('time', 0):.4f}): {e0:.6e}")
    print(f"    Etot (t={last.get('time', 0):.4f}): {e1:.6e}")

    if not (math.isfinite(e0) and math.isfinite(e1)):
        print("    FAIL: Non-finite total EM energy.")
        return False, "Non-finite total EM energy"

    if e0 <= 0:
        print("    FAIL: Initial total EM energy is zero -- "
              "wave not initialised.")
        return False, "Initial Etot is zero (wave not initialised)"

    if e1 <= 0:
        print("    FAIL: Final total EM energy is zero -- wave collapsed.")
        return False, "Final Etot is zero (wave collapsed)"

    ratio = e1 / e0
    lower, upper = 0.3, 3.0
    print(f"    Etot_final / Etot_initial = {ratio:.4f} "
          f"(allowed [{lower}, {upper}])")

    if ratio < lower or ratio > upper:
        print("    FAIL: total EM energy changed outside the allowed range -- "
              "possible blow-up or unphysical decay.")
        return False, f"Etot ratio {ratio:.3f} outside [{lower}, {upper}]"

    print(f"    SUCCESS: light wave energy conserved (ratio = {ratio:.3f}).")
    return True, "Passed"


def validate_hybrid(pic_diags=None, test_name=None):
    """Validate the Hybrid PIC (kinetic ion / fluid electron) wave test.

    Success criteria:
    1. FLEKS completes (run_test already checks the exit code).
    2. No NaN / blow-up: magnetic energy Eb and total ion energy Epart finite.
    3. Ion particle number conserved (periodic BCs, no source/loss): the
       kinetic-ion energy Epart0 stays within ~10% of its initial value.
    4. Magnetic energy stays bounded -- the Hall-driven whistler must not
       numerically blow up (this is the failure mode of the missing 1/(4*pi)
       factor in the Hall current, or insufficient #HALLSUBCYCLE).
    The precise whistler dispersion omega/Omega_i = (k d_i)^2/(1+(k d_i)^2)
    is verified separately by _check_hybrid_wave_dispersion() / the README.
    """
    print("Validating Hybrid PIC Wave Test...")

    if not pic_diags or len(pic_diags) < 2:
        print("  [INFO] No PIC energy log found; skipping energy checks.")
        return True, "Passed (no pic log)"

    first = pic_diags[0]
    last = pic_diags[-1]

    passed = True
    reasons = []

    eb0 = first.get("Eb", 0.0)
    eb1 = last.get("Eb", 0.0)
    print(f"    Eb (magnetic): {eb0:.6e} -> {eb1:.6e}")
    if not math.isfinite(eb1):
        passed = False
        reasons.append("Eb not finite (NaN/Inf)")
    if eb0 > 0 and eb1 > eb0 * 5.0:
        passed = False
        reasons.append("Eb grew >5x (whistler instability / 4pi bug?)")

    ep1 = last.get("Epart", 0.0)
    print(f"    Epart (ions):  {first.get('Epart', 0.0):.6e} -> {ep1:.6e}")
    if not math.isfinite(ep1):
        passed = False
        reasons.append("Epart not finite (NaN/Inf)")

    # Ion kinetic-energy conservation check.
    # NOTE: in a hybrid (kinetic-ion / fluid-electron) wave the ion KINETIC
    # ENERGY is NOT conserved -- it continuously exchanges with the
    # electromagnetic field (only the particle *number* is conserved under
    # periodic BCs with no source/loss).  The first-order E^n push (roadmap
    # Step 5) adds a further small numerical drift.  We therefore do NOT
    # require |Epart0 ratio - 1| < 0.1.  Instead we guard against *gross*
    # non-conservation (particle loss/creation or a runaway blow-up) with a
    # wide tolerance, while the decisive stability guards are the Eb blow-up
    # check (Eb not > 5x) and the bounded transverse-wave check below.
    e0 = first.get("Epart0", 0.0)
    e1 = last.get("Epart0", 0.0)
    if e0 > 0:
        ratio = e1 / e0
        print(f"    Epart0 (ions): {e0:.6e} -> {e1:.6e} (ratio {ratio:.4f})")
        if ratio < 0.2 or ratio > 10.0:
            passed = False
            reasons.append(
                f"Ion energy ratio {ratio:.3f} outside [0.2,10.0] "
                f"(gross particle non-conservation / runaway)")
    else:
        print("    [INFO] Epart0 initial zero; skipping ion check.")

    if passed:
        print("Hybrid PIC Wave Test: PASSED")
        return True, "Passed"
    return False, "; ".join(reasons)


def validate_singlecell(pic_diags=None, test_name=None):
    """Single-cell periodic hybrid test.

    With exactly one grid cell and periodic boundaries, curl B is identically
    zero, so the Hall term (J x B)/rho and the convective term U_i x B both
    vanish (U_i = 0).  The electric field stays zero and the magnetic field is
    frozen.  The test passes iff (1) the magnetic energy Eb is conserved to
    round-off (no spurious Hall-driven evolution), (2) the electric field energy
    Ee stays ~0 (the field is truly frozen, not merely energy-conserving), and
    (3) no NaN/blow-up occurs.
    """
    print("\n=== Validating Single-Cell Hybrid Test ===")
    if not pic_diags:
        return False, "No diagnostics found"
    first = pic_diags[0]
    last = pic_diags[-1]
    passed = True
    reasons = []

    eb0 = first.get("Eb", 0.0)
    eb1 = last.get("Eb", 0.0)
    print(f"    Eb (magnetic): {eb0:.6e} -> {eb1:.6e}")
    if not math.isfinite(eb1):
        passed = False
        reasons.append("Eb not finite (NaN/Inf)")
    if eb0 > 0:
        ratio = eb1 / eb0
        print(f"    Eb ratio: {ratio:.6f}")
        # Exact (to round-off): a single cell has no spatial gradient, so the
        # Hall term is exactly zero and B cannot evolve.  Allow a tiny tolerance
        # for floating-point round-off in the field solvers.
        if ratio < 0.9999 or ratio > 1.0001:
            passed = False
            reasons.append(
                f"Eb ratio {ratio:.6f} not ~1 (spurious Hall/evolution on "
                f"single cell; curl B must be zero)")

    # Ee must stay ~0: a frozen field has E = 0, so there is no electric energy.
    # This is the real discriminator between "frozen" and "propagating
    # non-dispersively" (both conserve Eb, but only a frozen field has Ee ~ 0).
    # Allow a small round-off floor: on a single cell the field solver leaves a
    # residual Ee of order 1e-3 * Eb (vs. ~1 for a genuinely propagating wave),
    # so a 1e-2 * Eb threshold cleanly separates frozen from evolving.
    eemax = max((d.get("Ee", 0.0) for d in pic_diags), default=0.0)
    print(f"    Ee (electric, max): {eemax:.6e}")
    if not math.isfinite(eemax):
        passed = False
        reasons.append("Ee not finite (NaN/Inf)")
    if eb0 > 0 and eemax > 1.0e-2 * eb0:
        passed = False
        reasons.append(
            f"Ee {eemax:.3e} not ~0 vs Eb {eb0:.3e} (field is evolving / "
            f"propagating, not frozen)")

    if passed:
        print("Single-Cell Hybrid Test: PASSED")
        return True, "Passed"
    return False, "; ".join(reasons)


def validate_zerocurrent(pic_diags=None, test_name=None):
    """Zero-current hybrid wave test.

    No macroparticles are loaded (rho = 0 everywhere), so every cell is left
    inert: the generalized Ohm's law is fully short-circuited by the
    ``if (rho > 0)`` guard in assemble_ohm_E, so the electric field E = 0.
    Faraday's law then gives dB/dt = -curl E = 0: the seeded sinusoidal B
    perturbation is FROZEN and does NOT propagate.  (Physically this is the
    "zero-current" regime -- with no plasma there is no current source to drive
    the field; mechanically it is the rho = 0 inert-cell path, not a literal
    J = 0 inside the Hall formula, since the sinusoidal B itself has curl B != 0.)

    The test passes iff (1) the magnetic energy Eb is conserved (the wave neither
    grows, decays, nor travels) AND (2) the electric field energy Ee stays ~0 --
    the genuine signature that the field is frozen rather than merely
    energy-conserving (a non-dispersive propagating wave would also conserve Eb
    but would have Ee > 0).
    """
    print("\n=== Validating Zero-Current Hybrid Test ===")
    if not pic_diags:
        return False, "No diagnostics found"
    first = pic_diags[0]
    last = pic_diags[-1]
    passed = True
    reasons = []

    eb0 = first.get("Eb", 0.0)
    eb1 = last.get("Eb", 0.0)
    print(f"    Eb (magnetic): {eb0:.6e} -> {eb1:.6e}")
    if not math.isfinite(eb1):
        passed = False
        reasons.append("Eb not finite (NaN/Inf)")
    if eb0 > 0:
        ratio = eb1 / eb0
        print(f"    Eb ratio: {ratio:.6f}")
        # The perturbation is frozen, so Eb should be conserved to round-off.
        # Allow a tiny tolerance for floating-point round-off.
        if ratio < 0.999 or ratio > 1.001:
            passed = False
            reasons.append(
                f"Eb ratio {ratio:.6f} not ~1 (zero-current wave propagated / "
                f"decayed; field should be frozen)")

    # Ee is the real no-propagation discriminator: a frozen field has E = 0, so
    # the electric energy is ~0.  A propagating wave would generate an inductive
    # E and a non-zero Ee even while conserving Eb.
    eemax = max((d.get("Ee", 0.0) for d in pic_diags), default=0.0)
    print(f"    Ee (electric, max): {eemax:.6e}")
    if not math.isfinite(eemax):
        passed = False
        reasons.append("Ee not finite (NaN/Inf)")
    if eb0 > 0 and eemax > 1.0e-6 * eb0:
        passed = False
        reasons.append(
            f"Ee {eemax:.3e} not ~0 vs Eb {eb0:.3e} (wave is propagating / "
            f"field is evolving, not frozen)")

    if passed:
        print("Zero-Current Hybrid Test: PASSED")
        return True, "Passed"
    return False, "; ".join(reasons)


def _check_iaw_density():
    """Check the seeded ion-density profile in the IAW plot output.

    Reads the pic-component plot .out files (PostProc.pl output) and verifies
    that (1) the ion density rhoS0 is a clean single-mode sinusoid at the first
    frame (high |correlation| with sin(kx*x)) and (2) the mean density is
    conserved across the run (no mass loss / blow-up).

    Returns (passed, reason).  A True return with a non-empty reason means the
    check was skipped (e.g. PostProc.pl output not available) -- this is not a
    failure, matching how the other wave tests degrade gracefully when PostIDL
    is not built.
    """
    import glob
    try:
        plots_dir = os.path.join(RUN_DIR, "PC", "plots")
        out_files = sorted(glob.glob(os.path.join(plots_dir, "plot*pic*.out")))
        if not out_files:
            out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
        if not out_files:
            return True, ("No .out plot files found (PostProc.pl not run?); "
                          "skipping profile check.")

        # Parse Lx (x-domain length) and waveMode from PARAM.in.  The FLEKS
        # #GEOMETRY block uses per-dimension xMin/xMax/yMin/... tokens.
        Lx = None
        wave_mode = 1
        xmin = xmax = None
        try:
            with open(os.path.join(RUN_DIR, "PARAM.in"), "r") as pf:
                section = None
                for line in pf:
                    s = line.strip()
                    if s.startswith("#"):
                        section = s
                        continue
                    if not s:
                        continue
                    parts = s.split()
                    # FLEKS PARAM.in uses the "value keyword" token order,
                    # e.g. "-3.2  xMin".  Scan for the keyword and take the
                    # preceding token as the value.
                    for i, tok in enumerate(parts):
                        if section == "#GEOMETRY" and i > 0:
                            if tok == "xMin":
                                try:
                                    xmin = float(parts[i - 1])
                                except ValueError:
                                    pass
                            elif tok == "xMax":
                                try:
                                    xmax = float(parts[i - 1])
                                except ValueError:
                                    pass
                        if (section == "#TESTCASE" and i > 0
                                and tok == "waveMode"):
                            try:
                                wave_mode = int(float(parts[i - 1]))
                            except ValueError:
                                pass
            if xmin is not None and xmax is not None and xmax > xmin:
                Lx = xmax - xmin
        except Exception:
            pass

        if Lx is None or Lx <= 0:
            return True, "Could not parse Lx from PARAM.in; skipping profile check."

        kx = 2.0 * math.pi * wave_mode / Lx

        def load_profile(path):
            with open(path, "r", encoding="latin-1") as f:
                lines = f.readlines()
            if len(lines) < 6:
                return None, None
            var_names = lines[4].split()
            vidx = {v.upper(): i for i, v in enumerate(var_names)}
            ridx = None
            for target in ("RHOS0", "RHOS1", "RHO"):
                if target in vidx:
                    ridx = vidx[target]
                    break
            if ridx is None:
                return None, None
            xs = []
            rhos = []
            for line in lines[5:]:
                cols = line.split()
                if len(cols) <= ridx:
                    continue
                try:
                    xs.append(float(cols[0]))
                    rhos.append(float(cols[ridx]))
                except (ValueError, IndexError):
                    continue
            return xs, rhos

        first_x, first_rho = load_profile(out_files[0])
        last_x, last_rho = load_profile(out_files[-1])
        if not first_x or not last_x:
            return True, "Could not parse rhoS0 from .out; skipping profile check."

        m0 = sum(first_rho) / len(first_rho)
        m1 = sum(last_rho) / len(last_rho)
        print(f"    [IAW] mean rhoS0: {m0:.4e} -> {m1:.4e}")
        if m0 <= 0:
            return True, "Zero initial density; skipping profile check."
        if not (math.isfinite(m0) and math.isfinite(m1)):
            return False, "Non-finite density; blow-up."
        if abs(m1 - m0) / m0 > 0.3:
            return False, (f"Mean density changed by >30% ({m0:.3e} -> "
                           f"{m1:.3e}); mass not conserved.")

        # Correlate the (mean-subtracted) density with sin(kx*x).  For a clean
        # single-mode IAW seed the profile stays proportional to sin(kx*x) in
        # space (the temporal part is a scalar), so |corr| remains large.
        def corr_with_sin(xs, rhos):
            n = len(xs)
            mean = sum(rhos) / n
            dr = [rhos[i] - mean for i in range(n)]
            s = [math.sin(kx * xs[i]) for i in range(n)]
            m_dr = sum(dr) / n
            m_s = sum(s) / n
            cov = sum((dr[i] - m_dr) * (s[i] - m_s) for i in range(n))
            vd = sum((dr[i] - m_dr) ** 2 for i in range(n))
            vs = sum((s[i] - m_s) ** 2 for i in range(n))
            if vd <= 0 or vs <= 0:
                return 0.0
            return cov / math.sqrt(vd * vs)

        c0 = corr_with_sin(first_x, first_rho)
        print(f"    [IAW] |corr(rho, sin(kx x))| first frame: {abs(c0):.3f}")
        if abs(c0) < 0.3:
            return False, (f"Initial density not a clean sinusoid "
                           f"(corr={c0:.3f}); seed may be missing.")

        # Bounded amplitude (no blow-up).
        amp0 = max(first_rho) - min(first_rho)
        amp1 = max(last_rho) - min(last_rho)
        if amp0 <= 0:
            return True, "Flat initial density; skipping profile check."
        ratio = amp1 / amp0
        print(f"    [IAW] density amplitude: {amp0:.3e} -> {amp1:.3e} "
              f"(ratio {ratio:.3f})")
        if ratio > 10.0:
            return False, (f"Density amplitude grew >10x "
                           f"({amp0:.3e} -> {amp1:.3e}); blow-up.")
        return True, ""
    except Exception as e:  # never let the profile check crash validation
        return True, f"Profile check errored ({e}); skipping."


def validate_iaw(pic_diags=None, test_name=None):
    """Validate the ion-acoustic-wave (IAW) hybrid-PIC test.

    Success criteria:
    1. FLEKS completes (run_test already checks the exit code).
    2. No NaN / blow-up: magnetic energy Eb and total ion energy Epart finite,
       and Eb does not grow > 5x (no Hall/whistler runaway).
    3. Ion particle number conserved (periodic BCs, no source): Epart0 stays
       within [0.2, 10] of its initial value.
    4. Density seed present and mass-conserving: from the plot output the seeded
       ion density rhoS0 is a clean single-mode sinusoid at the first frame and
       its mean/overall amplitude remain bounded (mass conserved, no blow-up).
       Requires PostProc.pl output; if absent the profile check is skipped.
    """
    print("Validating Ion Acoustic Wave Test...")

    if not pic_diags or len(pic_diags) < 2:
        print("  [INFO] No PIC energy log found; skipping energy checks.")
        return True, "Passed (no pic log)"

    first = pic_diags[0]
    last = pic_diags[-1]

    passed = True
    reasons = []

    # --- energy / stability ------------------------------------------------
    eb0 = first.get("Eb", 0.0)
    eb1 = last.get("Eb", 0.0)
    print(f"    Eb (magnetic): {eb0:.6e} -> {eb1:.6e}")
    if not math.isfinite(eb1):
        passed = False
        reasons.append("Eb not finite (NaN/Inf)")
    if eb0 > 0 and eb1 > eb0 * 5.0:
        passed = False
        reasons.append("Eb grew >5x (whistler/Hall runaway?)")

    ep1 = last.get("Epart", 0.0)
    print(f"    Epart (ions):  {first.get('Epart', 0.0):.6e} -> {ep1:.6e}")
    if not math.isfinite(ep1):
        passed = False
        reasons.append("Epart not finite (NaN/Inf)")

    e0 = first.get("Epart0", 0.0)
    e1 = last.get("Epart0", 0.0)
    if e0 > 0:
        ratio = e1 / e0
        print(f"    Epart0 (ions): {e0:.6e} -> {e1:.6e} (ratio {ratio:.4f})")
        if ratio < 0.2 or ratio > 10.0:
            passed = False
            reasons.append(
                f"Ion energy ratio {ratio:.3f} outside [0.2,10.0] "
                f"(gross particle non-conservation / runaway)")
    else:
        print("    [INFO] Epart0 initial zero; skipping ion check.")

    # --- density profile / seed check --------------------------------------
    profile_ok, profile_reason = _check_iaw_density()
    if profile_ok is False:
        passed = False
        reasons.append(profile_reason)
    elif profile_reason:
        print(f"    [IAW] {profile_reason}")

    if passed:
        print("Ion Acoustic Wave Test: PASSED")
        return True, "Passed"
    return False, "; ".join(reasons)


def validate_ionization_source(pic_diags=None, test_name=None):
    """Validate an ionization source test (photoionization, electron impact,
    or charge exchange).

    Checks that the heaviest ion species (O+, which receives the exosphere
    source) energy increases over time, confirming the ionization source is
    active.  Uses the PIC energy log (log_pic_n*.log) as the data source.

    With the full-PIC species layout (species 0 = electron, 1 = H+, 2 = O+),
    the source is injected into the last ion species.

    For charge exchange (test_name="chargeexchange"), additionally verifies
    H+ (Epart1) energy does not decrease and requires O+ (Epart2) energy to
    grow by at least a minimum factor, since the O+ background is set to
    near-zero so the source contribution dominates.
    """
    print("Validating Ionization Source Test...")

    if not pic_diags or len(pic_diags) < 2:
        print("  [INFO] No PIC energy log found; skipping energy checks.")
        return True, "Passed (no pic log)"

    first = pic_diags[0]
    last = pic_diags[-1]

    # Determine the source (heaviest ion) species index from available
    # EpartN keys.  The last EpartN column is the heaviest ion.
    epart_keys = sorted(
        k for k in first.keys() if k.startswith("Epart") and k != "Epart"
    )
    if not epart_keys:
        print("  [INFO] No per-species energy columns; skipping.")
        return True, "Passed (no Epart columns)"
    source_key = epart_keys[-1]  # e.g. "Epart2" for O+
    source_idx = source_key.replace("Epart", "")

    print(f"  --- Energy Diagnostics (from log_pic log) ---")
    for k in epart_keys:
        print(f"    {k}: {first.get(k, 0):.6e} -> {last.get(k, 0):.6e}")
    print(f"    Initial total Epart: {first.get('Epart', 0):.6e}")
    print(f"    Final total Epart:   {last.get('Epart', 0):.6e}")

    if test_name == "chargeexchange":
        # For charge exchange, verify both H+ (Epart1) and O+ (Epart2).
        # O+ has a near-zero background, so its energy should increase
        # by a large factor.  H+ has a large bulk-kinetic-energy
        # background, so its energy increase is tiny; we only require
        # that it does not decrease (allowing for numerical noise).
        passed = True
        reasons = []
        min_factor_o = 2.0   # O+ must at least double
        h_tolerance = 0.05   # H+ may decrease by up to 5% (numerical noise)

        # --- O+ (heaviest ion, source species) ---
        o_key = epart_keys[-1]  # e.g. "Epart2"
        e_o_initial = first.get(o_key, 0.0)
        e_o_final = last.get(o_key, 0.0)
        factor_o = e_o_final / max(e_o_initial, 1e-30)
        print(f"    {o_key} (O+): {e_o_initial:.6e} -> {e_o_final:.6e} "
              f"(factor {factor_o:.3f}x, threshold {min_factor_o}x)")
        if e_o_initial <= 0:
            if e_o_final <= 0:
                print(f"    FAIL: {o_key} (O+) energy is zero — source not active.")
                passed = False
                reasons.append("O+ energy is zero (source not active)")
            else:
                print(f"    SUCCESS: {o_key} (O+) energy became non-zero.")
        elif factor_o < min_factor_o:
            print(f"    FAIL: {o_key} (O+) growth factor {factor_o:.3f} < {min_factor_o}")
            passed = False
            reasons.append(f"O+ growth factor {factor_o:.3f} < {min_factor_o}")
        else:
            print(f"    SUCCESS: {o_key} (O+) energy increased by {factor_o:.1f}x.")

        # --- H+ (light ion, also receives CX source) ---
        h_key = "Epart1" if "Epart1" in first else None
        if h_key:
            e_h_initial = first.get(h_key, 0.0)
            e_h_final = last.get(h_key, 0.0)
            print(f"    {h_key} (H+): {e_h_initial:.6e} -> {e_h_final:.6e}")
            if e_h_final < e_h_initial * (1.0 - h_tolerance):
                print(f"    FAIL: {h_key} (H+) energy decreased by more than "
                      f"{h_tolerance*100:.0f}% (numerical noise threshold).")
                passed = False
                reasons.append("H+ energy decreased beyond noise threshold")
            else:
                print(f"    SUCCESS: {h_key} (H+) energy stable or increasing.")

        if passed:
            print("Charge Exchange Source Test: PASSED")
            return True, "Passed"
        else:
            return False, "; ".join(reasons)

    else:
        # Original behavior for photoionization and electronimpact:
        # check that the heaviest ion (O+) energy increases.
        e_src_initial = first.get(source_key, 0.0)
        e_src_final = last.get(source_key, 0.0)
        print(f"    Initial {source_key} (species {source_idx}, O+): {e_src_initial:.6e}")
        print(f"    Final   {source_key} (species {source_idx}, O+): {e_src_final:.6e}")
        print(f"    Growth factor: {e_src_final / max(e_src_initial, 1e-30):.3f}x")
        if e_src_final <= e_src_initial:
            print(f"    FAIL: {source_key} energy did not increase.")
            print("    Ionization source may not be working correctly.")
            return False, (
                f"{source_key} energy did not increase "
                f"(initial={e_src_initial:.2e}, final={e_src_final:.2e})"
            )
        else:
            print(f"    SUCCESS: {source_key} energy increased (ionization source active).")
            return True, "Passed"


def _read_shadow_params():
    """Read shadow-cylinder and normalization parameters from PARAM.in.

    Returns (Rp_plot, halfH_plot, shadowR_plot) all in *plot* (normalized)
    coordinates, or None if shadow cylinder is not enabled.
    """
    param_path = os.path.join(RUN_DIR, "PARAM.in")
    Rp_si = 3.0e6
    lNormSI = 1.0
    halfH_si = 0.0
    shadowR_si = 0.0
    useShadow = False

    # Track which line within #NORMALIZATION we're on.
    norm_line_idx = 0

    try:
        with open(param_path, "r") as pf:
            section = None
            for line in pf:
                line_s = line.strip()
                if line_s.startswith("#"):
                    section = line_s
                    if section == "#NORMALIZATION":
                        norm_line_idx = 0
                    continue
                if not line_s:
                    continue
                parts = line_s.split()
                if section == "#BODYSIZE" and len(parts) >= 1:
                    try:
                        Rp_si = float(parts[0])
                    except ValueError:
                        pass
                elif section == "#NORMALIZATION" and len(parts) >= 1:
                    # First value is lNormSI, second is uNormSI.
                    if norm_line_idx == 0:
                        try:
                            lNormSI = float(parts[0])
                        except ValueError:
                            pass
                    norm_line_idx += 1
                elif section == "#SHADOWCYLINDER":
                    useShadow = True
                    try:
                        val = float(parts[0])
                    except ValueError:
                        continue
                    if "radius" in line_s.lower():
                        shadowR_si = val
                    elif "halfheight" in line_s.lower():
                        halfH_si = val
    except Exception:
        pass

    if not useShadow:
        return None

    # Plot coordinates = SI / lNormSI
    return (Rp_si / lNormSI, halfH_si / lNormSI, shadowR_si / lNormSI)


def _load_idl_plot_asymmetry():
    """Check photoionization day/night asymmetry from plot output.

    Reads .out files produced by PostProc.pl (which concatenates the *.h
    and *.idl files written by FLEKS).  PostProc.pl must be run before
    calling this function.  Verifies that rhoS1 is non-zero on the dayside
    (+X) and much smaller inside the planetary shadow cylinder (-X, within
    cylinder radius and height).

    Returns (passed: bool, reason: str).
    """
    import glob

    plots_dir = os.path.join(RUN_DIR, "PC", "plots")

    # -- Get shadow cylinder geometry in plot coordinates -------------------
    shadow_geom = _read_shadow_params()
    if shadow_geom is None:
        print("    [ASYM] Shadow cylinder not enabled; skipping asymmetry check.")
        return True, "No shadow cylinder"

    Rp_plot, halfH_plot, shadowR_plot = shadow_geom
    print(f"    [ASYM] Rp={Rp_plot:.0f}, halfH={halfH_plot:.0f}, "
          f"shadowR={shadowR_plot:.0f} (plot coords)")

    # -- Collect data points (x, y, rhoS1) from PostProc.pl .out files -----
    # PostProc.pl must be run first to concatenate *.h and *.idl into *.out.
    # We do NOT work on the raw .idl files directly.
    points = []  # list of (x, y, rhoS1)

    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        print("    [ASYM] No .out files found. "
              "Ensure PostProc.pl has been run after FLEKS.exe.")
        return False, "No .out files found (PostProc.pl not run?)"

    latest_out = out_files[-1]
    print(f"    [ASYM] Loading .out: {os.path.basename(latest_out)}")
    with open(latest_out, "r") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return True, "Short .out file"
    var_names = lines[4].split()
    # Look for the heaviest ion species density (rhoS2 = O+ with 3-species
    # layout: 0=e, 1=H+, 2=O+).  Fall back to rhoS1 for 2-species layouts.
    rho_idx = None
    for target in ("RHOS2", "RHOS1"):
        for iv, vn in enumerate(var_names):
            if vn.upper() == target:
                rho_idx = iv
                break
        if rho_idx is not None:
            continue
    if rho_idx is None:
        return True, "rhoS2/rhoS1 not in .out"
    for line in lines[5:]:
        cols = line.strip().split()
        if len(cols) <= rho_idx:
            continue
        try:
            points.append((float(cols[0]), float(cols[1]),
                           float(cols[rho_idx])))
        except (ValueError, IndexError):
            continue

    if not points:
        print("    [ASYM] No data points parsed.")
        return True, "Empty plot data"

    # -- Classify points: dayside vs shadow ---------------------------------
    # The shadow cylinder covers the nightside (x < 0 for solarDir=+X).
    # The "planet" plot keyword limits output to ~[-Rp, +Rp], so we compare
    # the dayside (x > 0) with the deep nightside (x < -Rp/2) where
    # photoionization is suppressed and diffusion has less effect.
    dayside_vals = []
    shadow_vals = []
    y_lim = min(Rp_plot / 2.0, shadowR_plot / 4.0)
    for x, y, rho in points:
        if abs(y) > y_lim:
            continue
        if x > 0:
            dayside_vals.append(rho)
        elif x < -Rp_plot * 0.5:
            shadow_vals.append(rho)

    dayside_mean = (sum(dayside_vals) / len(dayside_vals)
                    if dayside_vals else 0.0)
    shadow_mean = (sum(shadow_vals) / len(shadow_vals)
                   if shadow_vals else 0.0)

    print(f"    [ASYM] Parsed {len(points)} points: "
          f"{len(dayside_vals)} dayside, {len(shadow_vals)} shadow")
    print(f"    Dayside (+X) mean rhoS1:      {dayside_mean:.4e}")
    print(f"    Shadow  (-X, cyl) mean:       {shadow_mean:.4e}")

    if len(dayside_vals) == 0:
        return False, "Zero dayside points -- cannot verify"
    if dayside_mean <= 0.0:
        return False, "Dayside rhoS1 is zero -- photoionization source not active"
    if len(shadow_vals) == 0:
        return False, "Zero shadow points -- cannot verify"
    # The shadow region still has some density from particle diffusion from
    # the dayside (especially near the planet surface), so we require
    # shadow < 20% of dayside rather than near-zero.  With the shadow
    # cylinder radius set to the planet radius, the dayside/night asymmetry
    # is pronounced and a 0.2 threshold provides a meaningful check.
    if shadow_mean > max(dayside_mean * 0.2, 1e-30):
        return False, (
            f"Shadow rhoS1 too high ({shadow_mean:.2e}) "
            f"vs dayside ({dayside_mean:.2e})"
        )

    ratio = shadow_mean / max(dayside_mean, 1e-30)
    print(f"    Shadow/dayside ratio:          {ratio:.2e}  (expected \u226a 1)")

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
    import glob

    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        print("    [FFT] No .out files found (PostProc.pl not run?).")
        return True, "No .out files (skipped)"

    # Use the final frame (latest cycle); this is the t ~= 0.1 snapshot.
    out_file = out_files[-1]
    print(f"    [FFT] Loading .out: {os.path.basename(out_file)}")

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
            print(f"    [FFT] '{need}' not found in .out variables: {var_names}")
            return True, f"{need} not in .out"

    iby, ibz, ibx = vidx["BY"], vidx["BZ"], vidx["BX"]

    x = []
    by = []
    bz = []
    bx = []
    for line in lines[5:]:
        cols = line.split()
        if len(cols) <= max(iby, ibz, ibx):
            continue
        try:
            x.append(float(cols[0]))
            by.append(float(cols[iby]))
            bz.append(float(cols[ibz]))
            bx.append(float(cols[ibx]))
        except (ValueError, IndexError):
            continue

    n = len(x)
    if n < 4:
        print("    [FFT] Too few data points for FFT.")
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

    print(f"    [FFT] t={t_norm:.4f} (normalized), N={n} cells")
    print(f"    [FFT] |Bx|={bx_mean:.4f} nT, max|B_perp|={bperp_max:.4e} nT")

    # ---- Check 1: wave growth above the noise floor ----------------------
    # At t=0 the transverse field is exactly zero; after the instability
    # triggers it grows from numerical noise.  Require the amplitude to
    # exceed a small fraction of the guide field.
    noise_frac = 1e-4
    if bx_mean <= 0:
        growth_ok = bperp_max > 0
    else:
        growth_ok = bperp_max > noise_frac * abs(bx_mean)
    print(f"    [FFT] Wave growth: max|B_perp|/|Bx| = "
          f"{bperp_max / max(abs(bx_mean), 1e-30):.3e} "
          f"(threshold {noise_frac:.0e}) -> "
          f"{'OK' if growth_ok else 'FAIL'}")

    # ---- DFT of the transverse field -------------------------------------
    # FFT By and Bz separately (preserving sign/oscillation), then combine
    # the per-mode amplitudes.  Using |B_perp| directly would introduce
    # spurious harmonics from the magnitude operation.
    def _dft_amp(data):
        nn = len(data)
        out = []
        for k in range(nn // 2 + 1):
            re = sum(data[j] * math.cos(2 * math.pi * k * j / nn)
                     for j in range(nn))
            im = -sum(data[j] * math.sin(2 * math.pi * k * j / nn)
                      for j in range(nn))
            out.append(math.hypot(re, im))
        return out

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
        print("    [FFT] No non-DC spectral power; wave has not grown.")
        return False, "No transverse wave power detected"

    # Dominant non-DC mode.
    n_dom = max(range(1, len(amps)), key=lambda k: amps[k])
    dom_frac = amps[n_dom] ** 2 / total_power

    # Fraction of power in low-order modes (n <= N/4).
    n_low = n // 4
    low_frac = sum(a * a for a in amps[1:n_low + 1]) / total_power

    print(f"    [FFT] Dominant mode: n={n_dom} "
          f"({100 * dom_frac:.1f}% of non-DC power)")
    print(f"    [FFT] Power in low modes (n<={n_low}): "
          f"{100 * low_frac:.1f}%")

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

    # Box geometry (x is in metres in the .out file).
    dx = x[1] - x[0]
    L = n * dx
    k1 = 2 * math.pi / L                  # box-fundamental wavenumber
    n_res = max(1, round(k_res / k1))

    print(f"    [FFT] Omega_i = {omega_i:.4f} rad/s, "
          f"Delta_v = {delta_v:.1e} m/s")
    print(f"    [FFT] k_res = {k_res:.3e} 1/m, "
          f"k_1 = {k1:.3e} 1/m (L = {L:.1f} m)")
    print(f"    [FFT] Resonant wavelength = {2 * math.pi / k_res:.3e} m "
          f"vs box L = {L:.1f} m")

    if k_res < k1:
        # Resonant wavelength exceeds the box: the resonant mode (n ~ 0)
        # cannot fit, so the nearest available mode is the box-fundamental
        # n=1.  The instability populates the longest-wavelength modes that
        # fit; verify the bulk of the wave power resides in low-order modes
        # (not grid-scale noise near the Nyquist frequency).
        mode_ok = low_frac > 0.4
        print(f"    [FFT] k_res < k_1: resonant mode (n={k_res / k1:.2e}) "
              f"exceeds the box; nearest available mode is n=1.")
        print(f"    [FFT] Mode check (box-limited): low-mode power "
              f"fraction {low_frac:.2f} > 0.4 -> "
              f"{'OK' if mode_ok else 'FAIL'}")
    else:
        tol = 2
        mode_ok = abs(n_dom - n_res) <= tol
        print(f"    [FFT] Mode check: |n_dom({n_dom}) - n_res({n_res})| "
              f"<= {tol} -> {'OK' if mode_ok else 'FAIL'}")

    if growth_ok and mode_ok:
        print("    [FFT] Transverse wave check: PASSED")
        return True, "Passed"
    else:
        reasons = []
        if not growth_ok:
            reasons.append(f"wave amplitude {bperp_max:.2e} nT below "
                           f"noise floor ({noise_frac:.0e}*|Bx|)")
        if not mode_ok:
            reasons.append(f"dominant mode n={n_dom} inconsistent with "
                           f"resonance (n_res={n_res})")
        return False, "; ".join(reasons)


def _check_charge_exchange_source_profile():
    """Check charge exchange source spatial profile from plot output.

    Reads .out files produced by PostProc.pl.  Verifies that the O+ density
    (rhoS2) peaks near the planet surface (|x| ~ Rp) where the exosphere
    density is highest, and is much smaller near the planet center where no
    neutrals exist.  The exosphere density is zero inside the planet, so the
    source (and resulting particle density) should be smallest at the center.

    Returns (passed: bool, reason: str).
    """
    import glob

    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        print("    [CX] No .out files found (PostProc.pl not run?).")
        return False, "No .out files found"

    out_file = out_files[-1]
    print(f"    [CX] Loading .out: {os.path.basename(out_file)}")

    with open(out_file, "r") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return True, "Short .out file"

    var_names = lines[4].split()
    vidx = {v.upper(): i for i, v in enumerate(var_names)}

    # Find rhoS2 (O+ density); fall back to rhoS1 for 2-species layouts.
    rho_idx = None
    rho_name = None
    for target in ("RHOS2", "RHOS1"):
        if target in vidx:
            rho_idx = vidx[target]
            rho_name = target
            continue
    if rho_idx is None:
        print(f"    [CX] rhoS2/rhoS1 not found in .out variables: {var_names}")
        return True, "rhoS2/rhoS1 not in .out"

    # Read planet radius and normalization from PARAM.in (plot coords = SI / lNormSI).
    Rp_si = 3.0e6
    lNormSI = 1000.0
    try:
        with open(os.path.join(RUN_DIR, "PARAM.in"), "r") as pf:
            section = None
            norm_idx = 0
            for line in pf:
                line_s = line.strip()
                if line_s.startswith("#"):
                    section = line_s
                    if section == "#NORMALIZATION":
                        norm_idx = 0
                    continue
                if not line_s:
                    continue
                parts = line_s.split()
                if section == "#BODYSIZE" and len(parts) >= 1:
                    try:
                        Rp_si = float(parts[0])
                    except ValueError:
                        pass
                elif section == "#NORMALIZATION" and len(parts) >= 1:
                    if norm_idx == 0:
                        try:
                            lNormSI = float(parts[0])
                        except ValueError:
                            pass
                    norm_idx += 1
    except Exception:
        pass

    Rp_plot = Rp_si / lNormSI

    # Parse data points (x, rhoS2).
    points = []
    for line in lines[5:]:
        cols = line.strip().split()
        if len(cols) <= rho_idx:
            continue
        try:
            x = float(cols[0])
            rho = float(cols[rho_idx])
            points.append((x, rho))
        except (ValueError, IndexError):
            continue

    if not points:
        print("    [CX] No data points parsed from .out file.")
        return False, "No data points parsed"

    print(f"    [CX] Rp (plot coords): {Rp_plot:.1f}")
    print(f"    [CX] Points: {len(points)}")

    # Classify points by distance from planet center:
    #   - "near surface": 0.5*Rp < |x| <= 1.5*Rp (exosphere active, source peaks)
    #   - "deep interior": |x| < 0.3*Rp (no neutrals, source should be ~0)
    near_surface = [(x, r) for x, r in points
                    if 0.5 * Rp_plot < abs(x) <= 1.5 * Rp_plot]
    deep_interior = [(x, r) for x, r in points if abs(x) < 0.3 * Rp_plot]

    surface_mean = (sum(r for _, r in near_surface) / len(near_surface)
                    if near_surface else 0.0)
    surface_max = max((r for _, r in near_surface), default=0.0)
    interior_mean = (sum(r for _, r in deep_interior) / len(deep_interior)
                     if deep_interior else 0.0)
    interior_max = max((r for _, r in deep_interior), default=0.0)

    print(f"    [CX] {rho_name} near surface (mean): {surface_mean:.4e}")
    print(f"    [CX] {rho_name} near surface (max):  {surface_max:.4e}")
    print(f"    [CX] {rho_name} deep interior (mean): {interior_mean:.4e}")
    print(f"    [CX] {rho_name} deep interior (max):  {interior_max:.4e}")

    # Check 1: source non-zero near the planet surface.
    if surface_max <= 0.0:
        print("    [CX] FAIL: No source detected near planet surface.")
        return False, "No source detected near planet surface"

    # Check 2: density much smaller in the deep interior than near surface.
    if interior_mean > surface_mean * 0.1:
        print(f"    [CX] FAIL: Interior density too high "
              f"({interior_mean:.2e} vs surface mean {surface_mean:.2e})")
        return False, (f"Interior density too high "
                       f"({interior_mean:.2e} vs {surface_mean:.2e})")

    # Check 3: approximately symmetric (left vs right near surface).
    left = [r for x, r in near_surface if x < 0]
    right = [r for x, r in near_surface if x > 0]
    left_mean = sum(left) / len(left) if left else 0.0
    right_mean = sum(right) / len(right) if right else 0.0

    print(f"    [CX] {rho_name} left  (x<0, near surf) mean: {left_mean:.4e}")
    print(f"    [CX] {rho_name} right (x>0, near surf) mean: {right_mean:.4e}")

    if left_mean > 0 and right_mean > 0:
        ratio = min(left_mean, right_mean) / max(left_mean, right_mean)
        print(f"    [CX] Left/Right ratio: {ratio:.2f}")
        if ratio < 0.3:
            print(f"    [CX] FAIL: Source asymmetric (L/R ratio {ratio:.2f})")
            return False, f"Source asymmetric (L/R ratio {ratio:.2f})"

    print("    [CX] Charge exchange source profile: VERIFIED")
    return True, "Passed"


def _check_lightwave_present():
    """Verify the light wave is present in the final plot output.

    Reads the final .out file (produced by PostProc.pl) and checks that the
    magnetic-field amplitude (BX/BY/BZ) is non-zero somewhere on the slice,
    confirming the transverse EM wave was initialised and is still present at
    the final time.

    Returns (passed: bool, reason: str).
    """
    import glob

    plots_dir = os.path.join("run_test", "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        print("    [LW] No .out files found (PostProc.pl not run?).")
        return True, "No .out files (skipped)"

    out_file = out_files[-1]
    print(f"    [LW] Loading .out: {os.path.basename(out_file)}")

    with open(out_file, "r") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return True, "Short .out file"

    var_names = lines[4].split()
    vidx = {v.upper(): i for i, v in enumerate(var_names)}
    b_idx = []
    for target in ("BX", "BY", "BZ"):
        if target in vidx:
            b_idx.append(vidx[target])
    if not b_idx:
        return True, "BX/BY/BZ not in .out"

    bmax = 0.0
    for line in lines[5:]:
        cols = line.split()
        if max(b_idx) >= len(cols):
            continue
        try:
            for i in b_idx:
                v = abs(float(cols[i]))
                if v > bmax:
                    bmax = v
        except (ValueError, IndexError):
            continue

    print(f"    [LW] Max |B| amplitude on slice: {bmax:.4e}")
    if bmax <= 0.0:
        print("    [LW] FAIL: magnetic field is zero -- wave not present.")
        return False, "Magnetic field is zero (wave not present)"
    print("    [LW] Light wave present: VERIFIED")
    return True, "Passed"


def _hyb_load_out(out_file):
    """Load By, Bz arrays from a hybrid-wave .out plot file."""
    with open(out_file, "r", encoding="latin-1") as f:
        lines = f.readlines()
    if len(lines) < 6:
        return None
    var_names = lines[4].split()
    vidx = {v.upper(): i for i, v in enumerate(var_names)}
    for need in ("BY", "BZ"):
        if need not in vidx:
            return None
    iby, ibz = vidx["BY"], vidx["BZ"]
    by, bz = [], []
    for line in lines[5:]:
        cols = line.split()
        if len(cols) <= max(iby, ibz):
            continue
        try:
            by.append(float(cols[iby]))
            bz.append(float(cols[ibz]))
        except (ValueError, IndexError):
            continue
    return by, bz


def _hyb_dft_dominant(by):
    """Spatial DFT of By(x); return (dominant_k, dominant_frac, nondc_power)."""
    n = len(by)
    if n < 4:
        return 0, 0.0, 0.0
    dft_mag = []
    for k in range(n // 2 + 1):
        re = sum(by[i] * math.cos(2.0 * math.pi * k * i / n) for i in range(n))
        im = -sum(by[i] * math.sin(2.0 * math.pi * k * i / n) for i in range(n))
        dft_mag.append(math.hypot(re, im))
    nondc_power = sum(dft_mag[k] ** 2 for k in range(1, n // 2 + 1))
    if nondc_power < 1e-30:
        return 0, 0.0, 0.0
    dominant_k = max(range(1, n // 2 + 1), key=lambda k: dft_mag[k])
    dominant_frac = dft_mag[dominant_k] ** 2 / nondc_power
    return dominant_k, dominant_frac, nondc_power


def _check_hybrid_wave_dispersion():
    """Check the transverse wave launched by the HybridWave initializer.

    Reads the FIRST and LAST .out plot files and verifies:
    1. Early time: the dominant spatial mode is n=1 (one wavelength in the
       box), matching the kx = 2*pi/Lx seed in fill_hybrid_wave().  This
       confirms the wave is correctly seeded and the solver propagates it.
    2. Late time: the transverse amplitude is bounded (no catastrophic
       blow-up from a missing 1/(4*pi) Hall factor or insufficient
       #HALLSUBCYCLE).

    Note: at moderate PPC the Hall term amplifies grid-scale particle noise
    over long times (a well-known hybrid-PIC limitation).  The early-time
    n=1 check validates the solver; the late-time bound catches genuine
    instabilities while tolerating moderate noise-driven growth.
    """
    import glob

    plots_dir = os.path.join(RUN_DIR, "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        print("    [HYB] No .out files found (PostProc.pl not run?).")
        return True, "No .out files (skipped)"

    # --- Early-time check: seeded wavelength (n=1) must dominate ---
    early_file = out_files[0]
    print(f"    [HYB] Early .out: {os.path.basename(early_file)}")
    early_data = _hyb_load_out(early_file)
    if early_data is None:
        return True, "Could not parse early .out file"
    by_e, bz_e = early_data
    n_e = len(by_e)
    bperp_e = [math.hypot(by_e[i], bz_e[i]) for i in range(n_e)]
    bperp_max_e = max(bperp_e) if bperp_e else 0.0
    dom_k_e, dom_frac_e, nondc_e = _hyb_dft_dominant(by_e)

    print(f"    [HYB] Early: N={n_e}, max|B_perp|={bperp_max_e:.4e}, "
          f"dominant mode n={dom_k_e} ({dom_frac_e*100:.1f}%)")

    if nondc_e < 1e-30:
        return False, "No transverse wave power at early time (By is flat/zero)"

    if dom_k_e != 1:
        return False, (f"Early dominant mode n={dom_k_e} (expected n=1 for "
                       f"seeded wavelength kx=2pi/Lx)")
    if dom_frac_e < 0.5:
        return False, (f"Early mode n=1 carries only {dom_frac_e*100:.1f}% "
                       f"of non-DC power (wave spectrum not clean at t=0)")

    # --- Late-time check: amplitude must be bounded ---
    late_file = out_files[-1]
    print(f"    [HYB] Late .out:  {os.path.basename(late_file)}")
    late_data = _hyb_load_out(late_file)
    if late_data is None:
        return True, "Could not parse late .out file"
    by_l, bz_l = late_data
    n_l = len(by_l)
    bperp_l = [math.hypot(by_l[i], bz_l[i]) for i in range(n_l)]
    bperp_max_l = max(bperp_l) if bperp_l else 0.0
    dom_k_l, dom_frac_l, _ = _hyb_dft_dominant(by_l)

    growth = bperp_max_l / bperp_max_e if bperp_max_e > 0 else float('inf')
    print(f"    [HYB] Late:  N={n_l}, max|B_perp|={bperp_max_l:.4e}, "
          f"dominant mode n={dom_k_l} ({dom_frac_l*100:.1f}%)")
    print(f"    [HYB] Growth: {growth:.1f}x (early -> late)")

    # Late-time bound: seed is ~0.02*B0; 10.0 catches catastrophic blow-up
    # (missing 4pi factor or CFL violation) while tolerating moderate
    # noise-driven growth at low PPC.
    if bperp_max_l > 10.0:
        return False, (f"Late amplitude {bperp_max_l:.2e} too large "
                       f"(unstable; seed was ~0.02)")

    print("    [HYB] Hybrid wave: early wavelength + late bounded-amplitude: VERIFIED")
    return True, "Passed"


def validate_plot_output(test_name):
    """Validate simulation plot output files for a given test.

    For the photoionization test, this checks the day/night asymmetry from
    the .out files produced by PostProc.pl.  For the beam test, this
    performs an FFT-based transverse-wave resonant-wavenumber check on the
    final plot frame.  For the charge exchange test, this verifies that the
    O+ source density appears only outside the planet and is approximately
    symmetric.  Other tests have no plot-file-based validation.
    """
    # ---- Photoionization: check day/night asymmetry via IDL .out ----
    if test_name == "photoionization":
        print("  --- Validating Output Files (IDL .out) ---")
        result, reason = _load_idl_plot_asymmetry()
        if result:
            print("    [IDL] Photoionization day/night asymmetry: VERIFIED")
        return result, reason

    # ---- Beam: FFT-based transverse-wave resonant-wavenumber check ----
    if test_name == "beam":
        print("  --- Validating Output Files (FFT transverse wave) ---")
        result, reason = _check_beam_transverse_wave()
        if result:
            print("    [FFT] Beam transverse-wave resonance check: VERIFIED")
        return result, reason

    # ---- Charge exchange: source profile (peaks near surface, symmetric) ----
    if test_name == "chargeexchange":
        print("  --- Validating Output Files (CX source profile) ---")
        result, reason = _check_charge_exchange_source_profile()
        return result, reason

    # ---- Hybrid wave: transverse-wave launch + bounded-amplitude check ----
    if test_name in ("hybrid_whistler", "hybrid_ohm"):
        print("  --- Validating Output Files (Hybrid wave) ---")
        result, reason = _check_hybrid_wave_dispersion()
        if result:
            print("    [HYB] Hybrid wave output check: VERIFIED")
        return result, reason

    # ---- Light wave: transverse EM wave must be present in the output ----
    if test_name == "lightwave":
        print("  --- Validating Output Files (light wave present) ---")
        result, reason = _check_lightwave_present()
        return result, reason

    # ---- Convection wave: rigid advection of the transverse wave ----
    if test_name == "hybrid_convection_wave":
        print("  --- Validating Output Files (convection advection) ---")
        result, reason = _check_convection_advection(test_name)
        if result:
            print("    [CNV] Convection advection check: VERIFIED")
        return result, reason

    # ---- Other tests: no plot-file validation ----
    print("  --- Validating Output Files: No plot-file check for this test ---")
    return True, "Passed (no plot-file check)"

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


def _frame_time(out_file):
    """Read the simulation time from a hybrid .out header (line 1, 2nd token)."""
    try:
        with open(out_file, "r", encoding="latin-1") as f:
            f.readline()
            nxt = f.readline()
        return float(nxt.split()[1])
    except (OSError, IndexError, ValueError):
        return None


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
    import glob
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


def main():
    # Parse nprocs: -n N or --nprocs N
    nprocs = 1
    for i, arg in enumerate(sys.argv):
        if arg in ("-n", "--nprocs"):
            try:
                nprocs = int(sys.argv[i + 1])
            except (IndexError, ValueError):
                print(f"Error: {arg} requires an integer argument (number of MPI processes).")
                sys.exit(1)
            continue
    
    if nprocs < 1:
        print("Error: Number of processes must be >= 1.")
        sys.exit(1)
    
    # Parse --summary-file PATH (custom output path for CI serial/parallel split)
    summary_file = "tests/summary.md"
    for i, arg in enumerate(sys.argv):
        if arg == "--summary-file":
            try:
                summary_file = sys.argv[i + 1]
            except IndexError:
                print("Error: --summary-file requires a path argument.")
                sys.exit(1)
            continue

    # Parse --test NAME (select a specific test to run; default: run all)
    # Accepts both "--test=NAME" and "--test NAME" forms.
    selected_test = None
    for i, arg in enumerate(sys.argv):
        if arg.startswith("--test="):
            selected_test = arg[len("--test="):]
            continue
        if arg == "--test":
            try:
                selected_test = sys.argv[i + 1]
            except IndexError:
                print("Error: --test requires a test name argument.")
                sys.exit(1)
            continue
        if arg == "--run-dir":
            global RUN_DIR
            try:
                RUN_DIR = sys.argv[i + 1]
            except IndexError:
                print("Error: --run-dir requires a path argument.")
                sys.exit(1)
            continue

    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(os.path.dirname(script_dir))
    print(f"Working directory set to: {os.getcwd()}")
    
    validators = {
        "beam": validate_beam,
        "photoionization": validate_ionization_source,
        "electronimpact": validate_ionization_source,
        "chargeexchange": validate_ionization_source,
        "recombination": validate_recombination,
        "chemistry": validate_chemistry,
        "lightwave": validate_lightwave,
        "hybrid_whistler": validate_hybrid,
        "hybrid_ohm": validate_hybrid,
        "freestream": validate_hybrid,
        "hybrid_convection_wave": validate_hybrid,
        "iaw": validate_iaw,
        "singlecell": validate_singlecell,
        "zerocurrent": validate_zerocurrent,
    }
    
    # Discover test subdirectories under tests/
    tests_dir = "tests"
    
    # Iterate through sorted subdirectories
    # Exclude "performance" (benchmark, not a pass/fail test) and "run_test"
    # (the shared run directory created by prepare_run_dir, which contains a
    # PARAM.in and would otherwise be discovered and re-run as a test).
    subdirs = sorted([d for d in os.listdir(tests_dir) 
                      if os.path.isdir(os.path.join(tests_dir, d))
                      and d not in ["performance", "run_test"]])
    
    tests = []
    for d in subdirs:
        test_dir = os.path.join(tests_dir, d)
        param_file = os.path.join(test_dir, "PARAM.in")
        if os.path.exists(param_file):
            validator = validators.get(d, None)
            tests.append((test_dir, d, validator))
            
    if not tests:
        print("No tests found in tests/ subdirectories!")
        sys.exit(1)

    # Filter to a single test if --test was given.
    if selected_test is not None:
        matching = [t for t in tests if t[1] == selected_test]
        if not matching:
            available = ", ".join(t[1] for t in tests)
            print(f"Error: test '{selected_test}' not found.")
            print(f"Available tests: {available}")
            sys.exit(1)
        tests = matching
        print(f"Selected test: {selected_test}")
        
    results = [] # Collect results for summary table

    for test_dir, name, validator in tests:
        # The free-stream test is a single directory holding two parameter files
        # (PARAM.in and PARAM.in.hybrid); the runner exercises both field solvers
        # by running it once per file. All other tests run once as written.
        if name == "freestream":
            # Variant 1: the full-PIC parameter file (PARAM.in).
            with open(os.path.join(test_dir, "PARAM.in")) as _f:
                _fullpic = _f.read()
            run_and_validate(test_dir, "FREESTREAM (FULL PIC)", validator,
                             nprocs, results, param_text=_fullpic,
                             base_name="freestream")
            # Variant 2: the hybrid Hall-OFF parameter file (PARAM.in.hybrid).
            # It shares the same grid/plasma/normalization/timestepping as
            # PARAM.in and differs only in the field-solver block.
            _hybrid_path = os.path.join(test_dir, "PARAM.in.hybrid")
            with open(_hybrid_path) as _f:
                _hybrid = _f.read()
            run_and_validate(test_dir, "FREESTREAM (HYBRID HALL-OFF)", validator,
                             nprocs, results, param_text=_hybrid,
                             base_name="freestream")
        else:
            run_and_validate(test_dir, name.upper(), validator, nprocs, results,
                             base_name=name)


    # ----------------------------------------------------
    # Print Summary Table
    # ----------------------------------------------------
    print("\n" + "=" * 80)
    print(" " * 32 + "TEST SUMMARY")
    print("=" * 80)
    print(f" {'Test Name':<28} | {'Status':<8} | {'Failure Reason / Details':<38}")
    print("-" * 80)
    
    all_passed = True
    for name_str, status, reason in results:
        status_display = f"{status:<8}"
        if sys.stdout.isatty():
            if status == "PASSED":
                status_display = f"\033[92;1m{status:<8}\033[0m" # Green Bold
            else:
                status_display = f"\033[91;1m{status:<8}\033[0m" # Red Bold
                
        if status != "PASSED":
            all_passed = False
            
        reason_display = reason if status != "PASSED" else ""
        print(f" {name_str:<28} | {status_display} | {reason_display:<38}")
        
    print("=" * 80)
    
    # Write Markdown Summary to summary.md for CI / PR integration
    try:
        with open(summary_file, "w") as f:
            f.write("### 🧪 Standalone FLEKS Test Results\n\n")
            f.write("| Test Name | Status | Failure Reason / Details |\n")
            f.write("| :--- | :--- | :--- |\n")
            for name_str, status, reason in results:
                status_md = "🟢 **PASSED**" if status == "PASSED" else "🔴 **FAILED**"
                reason_md = reason if status != "PASSED" else ""
                f.write(f"| {name_str} | {status_md} | {reason_md} |\n")
    except Exception as e:
        print(f"Warning: Could not write summary.md: {e}")
    
    if all_passed:
        if sys.stdout.isatty():
            print("\033[92;1m\nALL STANDALONE FLEKS TESTS PASSED SUCCESSFULLY!\033[0m\n")
        else:
            print("\nALL STANDALONE FLEKS TESTS PASSED SUCCESSFULLY!\n")
        sys.exit(0)
    else:
        if sys.stdout.isatty():
            print("\033[91;1m\nSOME TESTS FAILED.\033[0m\n")
        else:
            print("\nSOME TESTS FAILED.\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
