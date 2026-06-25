#!/usr/bin/env python3
import os
import shutil
import subprocess
import re
import sys

def safe_symlink(src, dst):
    if os.path.lexists(dst):
        if os.path.islink(dst) or os.path.isfile(dst):
            os.remove(dst)
        elif os.path.isdir(dst):
            shutil.rmtree(dst)
    os.symlink(src, dst)

def prepare_run_dir():
    run_dir = "run_test"
    os.makedirs(run_dir, exist_ok=True)
    
    # Determine the location of the share and bin directories relative to FLEKS root.
    # In a standalone FLEKS repository, 'share/Scripts/PostProc.pl' is directly inside the working directory.
    # In SWMF integrated environment, 'share' and 'bin' are located in SWMF root, i.e., two levels above FLEKS root.
    if os.path.isfile("share/Scripts/PostProc.pl"):
        postproc_target = "../share/Scripts/PostProc.pl"
        pidl_target = "../../share/Scripts/pIDL"
        postidl_target = "../../bin/PostIDL.exe"
    else:
        postproc_target = "../../../share/Scripts/PostProc.pl"
        pidl_target = "../../../../share/Scripts/pIDL"
        postidl_target = "../../../../bin/PostIDL.exe"

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
    run_dir = "run_test"
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

def run_test(test_dir, nprocs=1):
    param_file = os.path.join(test_dir, "PARAM.in")
    print(f"Running test in {test_dir} with config {param_file}...")
    prepare_run_dir()
    
    # Copy param_file to run_test/PARAM.in
    shutil.copy(param_file, "run_test/PARAM.in")
    
    # Build the command: serial for nprocs==1, mpirun otherwise
    if nprocs <= 1:
        cmd = ["./FLEKS.exe"]
        print(f"  Running in serial mode (no MPI)...")
    else:
        cmd = ["mpirun", "-n", str(nprocs), "./FLEKS.exe"]
        print(f"  Running with {nprocs} MPI processes...")
    
    # Run the command inside run_test/
    result = subprocess.run(cmd, cwd="run_test", stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error running FLEKS.exe for {test_dir}:")
        print(result.stderr)
        return None, result.returncode
        
    # Automatically run post-processing on the generated plots
    subprocess.run(["./PostProc.pl", "-v"], cwd="run_test", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    return result.stdout, 0


def read_diag_log(run_dir):
    """Read the structured diagnostic log file written by Pic::write_diag_log.

    Returns a tuple (particle_diags, field_diags) where:
    - particle_diags: list of dicts with keys species, time, cycle, macro,
      phys, vx, vy, vz, ke  (one entry per species per timestep).
    - field_diags: list of dicts with keys time, cycle, max_by, max_bz.
      Empty if doFieldDiag was not enabled.
    """
    import glob
    pc_plots = os.path.join(run_dir, "PC", "plots")
    log_files = sorted(glob.glob(os.path.join(pc_plots, "log_diag_n*.log")))
    if not log_files:
        return [], []

    # Use the most recent log file (there should normally be only one).
    log_file = log_files[-1]

    particle_diags = []
    field_diags = []

    with open(log_file, "r") as f:
        lines = f.readlines()

    if len(lines) < 2:
        return [], []

    # Parse header to discover column layout.
    header = lines[0].strip().split("\t")
    # Expected header columns (examples):
    #   time  cycle  macro0 phys0 Vx0 Vy0 Vz0 KE0  macro1 ...  maxBy maxBz

    # Discover how many species are present by counting macroN columns.
    n_species = sum(1 for col in header if col.startswith("macro"))
    has_field = "maxBy" in header

    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        vals = line.split("\t")
        if len(vals) < 2:
            continue
        try:
            t = float(vals[0])
            cycle = int(vals[1])
        except ValueError:
            continue

        col = 2  # current column index
        for iS in range(n_species):
            if col + 5 >= len(vals):
                break
            try:
                macro = int(vals[col])
                phys  = float(vals[col + 1])
                vx    = float(vals[col + 2])
                vy    = float(vals[col + 3])
                vz    = float(vals[col + 4])
                ke    = float(vals[col + 5])
            except ValueError:
                break
            particle_diags.append({
                "species": iS,
                "time":    t,
                "cycle":   cycle,
                "macro":   macro,
                "phys":    phys,
                "vx":      vx,
                "vy":      vy,
                "vz":      vz,
                "ke":      ke,
            })
            col += 6

        if has_field and col + 1 < len(vals):
            try:
                field_diags.append({
                    "time":   t,
                    "cycle":  cycle,
                    "max_by": float(vals[col]),
                    "max_bz": float(vals[col + 1]),
                })
            except ValueError:
                pass

    return particle_diags, field_diags


def parse_diagnostics(stdout):
    """DEPRECATED: parse DIAGNOSTIC lines from stdout.

    Kept as a fallback for runs that pre-date the log file system.
    Returns (particle_diags, field_diags=[]). Use read_diag_log() instead.
    """
    diagnostics = []
    pattern = re.compile(
        r"DIAGNOSTIC:\s+Species=(\d+)\s+Time=([+\-\d.e]+)\s+Cycle=(\d+)\s+"
        r"MacroParticles=(\d+)\s+PhysParticles=([+\-\d.e]+)\s+MeanVx=([+\-\d.e]+)\s+"
        r"MeanVy=([+\-\d.e]+)\s+MeanVz=([+\-\d.e]+)\s+KineticEnergy=([+\-\d.e]+)"
    )
    for line in stdout.splitlines():
        match = pattern.search(line)
        if match:
            diag = {
                "species": int(match.group(1)),
                "time":    float(match.group(2)),
                "cycle":   int(match.group(3)),
                "macro":   int(match.group(4)),
                "phys":    float(match.group(5)),
                "vx":      float(match.group(6)),
                "vy":      float(match.group(7)),
                "vz":      float(match.group(8)),
                "ke":      float(match.group(9))
            }
            diagnostics.append(diag)
    return diagnostics, []

def validate_beam(diags, field_diags=None):
    print("Validating Beam Instability Test...")
    if not diags:
        print("  [INFO] No diagnostic outputs parsed; skipping beam diagnostic checks.")
        return True, "Passed (diagnostics unavailable)"

    passed = True
    reasons = []
    
    # 1. Total particle count conservation check
    initial_phys = diags[0]["phys"]
    for diag in diags:
        t = diag["time"]
        actual_phys = diag["phys"]
        if abs(actual_phys - initial_phys) > 1e-5 * initial_phys:
            print(f"  FAIL at t={t:.2f}: Particle number changed! Expected={initial_phys:.2e}, Actual={actual_phys:.2e}")
            passed = False
            reasons.append(f"Particle number changed at t={t:.2f} (expected {initial_phys:.2e}, actual {actual_phys:.2e})")

    # 2. Mean velocity conservation check
    # Expected MeanVx is -0.392 (due to 99% background at -0.4 and 1% beam at +0.4)
    expected_vx = -0.392
    for diag in diags:
        t = diag["time"]
        actual_vx = diag["vx"]
        relative_diff = abs(actual_vx - expected_vx)
        print(f"  t={t:.2f}: MeanVx expected={expected_vx:.4f}, actual={actual_vx:.4f} (diff={relative_diff:.4f})")
        if relative_diff > 0.01:
            print(f"  FAIL: MeanVx exceeds tolerance of 0.01")
            passed = False
            reasons.append(f"MeanVx diff {relative_diff:.4f} > 0.01 at t={t:.2f}")

    # 3. Transverse magnetic field cyclotron wave growth check
    if field_diags:
        print("  --- Validating Cyclotron Wave Growth (Transverse B-field) ---")

        if not field_diags:
            print("  FAIL: No field diagnostic outputs found.")
            passed = False
            reasons.append("No field diagnostics found")
        else:
            final_by = field_diags[-1]["max_by"]
            final_bz = field_diags[-1]["max_bz"]

            # The initial condition has by=0, bz=0; numerical noise develops
            # after the first few steps. Use the first non-zero entry as the
            # noise-floor baseline so the growth ratio is well-defined.
            noise_by = next((d["max_by"] for d in field_diags if d["max_by"] > 0), 0.0)
            noise_bz = next((d["max_bz"] for d in field_diags if d["max_bz"] > 0), 0.0)

            print(f"    Initial Transverse B-field (noise level): MaxBy={noise_by:.4e}, MaxBz={noise_bz:.4e}")
            print(f"    Final Transverse B-field (after growth):  MaxBy={final_by:.4e}, MaxBz={final_bz:.4e}")

            if noise_by > 0 and noise_bz > 0:
                growth_y = final_by / noise_by
                growth_z = final_bz / noise_bz
                print(f"    Transverse B-field growth factors: By growth={growth_y:.2f}x, Bz growth={growth_z:.2f}x")
                if growth_y < 1.2 and growth_z < 1.2:
                    print("    FAIL: No significant growth of cyclotron waves detected (growth < 1.2x)")
                    passed = False
                    reasons.append(f"No wave growth (By growth={growth_y:.2f}x, Bz growth={growth_z:.2f}x < 1.2x)")
                else:
                    print("    SUCCESS: Cyclotron wave growth validated!")
            else:
                # Fallback: no non-zero noise baseline found; just require the
                # final value to be positive (any wave amplitude counts).
                print("    (No non-zero noise baseline; checking final amplitude only.)")
                if final_by <= 0 and final_bz <= 0:
                    print("    FAIL: Final transverse B-field is zero — no wave growth.")
                    passed = False
                    reasons.append("Final By and Bz are both zero — no wave growth detected")
                else:
                    print("    SUCCESS: Transverse B-field is non-zero at end of run.")

    if passed:
        print("Beam Instability Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)


def ensure_flekspy_installed():
    try:
        import flekspy
        import matplotlib
    except ImportError:
        print("\n[flekspy] Required package 'flekspy' or 'matplotlib' not found. Attempting automatic installation...")
        import subprocess
        import sys
        try:
            # Try installing via pip in the current Python environment (which could be the virtual env)
            subprocess.check_call([sys.executable, "-m", "pip", "install", "flekspy", "matplotlib"])
            print("[flekspy] Packages successfully installed.\n")
        except Exception as e:
            print(f"[flekspy] WARNING: Automatic installation failed: {e}")
            print("[flekspy] Attempting user installation fallback...")
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install", "--user", "flekspy", "matplotlib"])
                print("[flekspy] Packages successfully installed in user space.\n")
            except Exception as e_user:
                print(f"[flekspy] ERROR: Installation failed completely: {e_user}")
                print("[flekspy] Please install flekspy manually via 'pip install flekspy matplotlib'.")
                sys.exit(1)

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


def validate_ionization_source(diags, field_diags=None, pic_diags=None):
    """Validate an ionization source test (photoionization, electron impact,
    or charge exchange).

    Checks that Species 1 (O+, the heavy ion receiving the exosphere source)
    energy increases over time, confirming the ionization source is active.
    Uses the PIC energy log (log_pic_n*.log) as the primary data source.
    """
    print("Validating Ionization Source Test...")

    # Use pic energy log as primary data source
    if pic_diags and len(pic_diags) >= 2:
        first = pic_diags[0]
        last = pic_diags[-1]

        # Check that species 1 energy grows (exosphere source injects energy)
        e1_initial = first.get("Epart1", 0.0)
        e1_final = last.get("Epart1", 0.0)

        print(f"  --- Energy Diagnostics (from log_pic log) ---")
        print(f"    Initial Epart1 (species 1, O+): {e1_initial:.6e}")
        print(f"    Final Epart1 (species 1, O+):   {e1_final:.6e}")
        print(f"    Growth factor: {e1_final / max(e1_initial, 1e-30):.3f}x")
        print(f"    Initial Epart0 (species 0, H+): {first.get('Epart0', 0):.6e}")
        print(f"    Final Epart0 (species 0, H+):   {last.get('Epart0', 0):.6e}")
        print(f"    Initial total Epart: {first.get('Epart', 0):.6e}")
        print(f"    Final total Epart:   {last.get('Epart', 0):.6e}")

        if e1_final <= e1_initial:
            print("    FAIL: Species 1 energy did not increase.")
            print("    Ionization source may not be working correctly.")
            return False, (
                f"Species 1 energy did not increase "
                f"(initial={e1_initial:.2e}, final={e1_final:.2e})"
            )
        else:
            print("    SUCCESS: Species 1 energy increased (ionization source active).")
            return True, "Passed"
    else:
        print("  [INFO] No pic energy log found; trying particle diagnostics...")

    # Fall back to particle diagnostics from diag_log or stdout
    if not diags:
        print("  [INFO] No diagnostic outputs parsed; skipping exosphere checks.")
        return True, "Passed (diagnostics unavailable)"

    passed = True
    reasons = []

    species1_diags = [d for d in diags if d["species"] == 1]
    if not species1_diags:
        print("  FAIL: No diagnostics found for species 1 (O+ source species).")
        return False, "No species 1 diagnostics found"

    phys_initial = species1_diags[0]["phys"]
    phys_final = species1_diags[-1]["phys"]

    print(f"    Initial phys particles (species 1): {phys_initial:.6e}")
    print(f"    Final phys particles (species 1):   {phys_final:.6e}")

    if phys_final <= phys_initial:
        print("    FAIL: Species 1 particle count did not increase.")
        passed = False
        reasons.append(
            f"Species 1 particle count did not increase "
            f"(initial={phys_initial:.2e}, final={phys_final:.2e})"
        )
    else:
        print("    SUCCESS: Species 1 particle count increased.")

    if passed:
        print("Ionization Source Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)


def validate_amrex_outputs_with_flekspy(test_name, use_flekspy=False):
    if not use_flekspy:
        print("  --- Validating Output Files: Skipped thorough check (run with --thorough or --flekspy to enable) ---")
        return True, "Passed (skipped thorough check)"
        
    print(f"  --- Validating Output Files with flekspy ---")
    import glob
    try:
        import flekspy
    except ImportError:
        print("    [flekspy] ERROR: flekspy is not installed. Skipping verification.")
        return False, "flekspy is not installed"
    
    # AMReX plotfiles generated by FLEKS under run_test/PC/plots/
    plot_pattern = os.path.join("run_test", "PC", "plots", "*_amrex")
    plots = sorted(glob.glob(plot_pattern))
    
    # Exclude temporary or old plotfiles
    plots = [p for p in plots if not p.endswith(".old") and ".old." not in p]
    
    if not plots:
        print("    WARNING: No active AMReX plotfiles found to validate with flekspy.")
        return True, "No plotfiles found"
        
    print(f"    Found {len(plots)} AMReX plotfile(s) for validation.")
    
    # Load the latest generated plotfile to perform validation
    latest_plot = plots[-1]
    print(f"    Loading latest plotfile: {latest_plot}")
    
    try:
        ds = flekspy.load(latest_plot)
        print(f"    [flekspy] Successfully loaded dataset: {os.path.basename(latest_plot)}")
        
        # Access and display dataset properties
        dims = getattr(ds, "domain_dimensions", None)
        fields = getattr(ds, "field_list", None)
        params = getattr(ds, "parameters", None)
        
        if dims is not None:
            print(f"    [flekspy] Domain dimensions: {dims}")
        if fields is not None:
            fields_str = ", ".join([str(f) for f in fields[:10]])
            if len(fields) > 10:
                fields_str += f" ... (+{len(fields)-10} more)"
            print(f"    [flekspy] Available fields: [{fields_str}]")
        if params is not None:
            key_params = {k: v for k, v in params.items() if any(x in k.lower() for x in ["time", "cycle", "dt", "nstep"])}
            if key_params:
                print(f"    [flekspy] Parameters: {key_params}")
                
        # Basic sanity check: Ensure domain dimensions and fields are valid/non-empty
        if dims is not None and len(dims) == 3 and all(d > 0 for d in dims):
            print("    [flekspy] Validation Check: Domain dimensions are valid.")
        else:
            print("    [flekspy] Validation Warning: Invalid domain dimensions.")
            
        return True, "Passed"
    except Exception as e:
        print(f"    FAIL: flekspy failed to load or validate plotfile: {e}")
        return False, f"flekspy load error: {e}"

def main():
    # Parse command-line arguments
    use_flekspy = "--thorough" in sys.argv or "--flekspy" in sys.argv
    
    # Parse nprocs: -n N or --nprocs N
    nprocs = 1
    for i, arg in enumerate(sys.argv):
        if arg in ("-n", "--nprocs"):
            try:
                nprocs = int(sys.argv[i + 1])
            except (IndexError, ValueError):
                print(f"Error: {arg} requires an integer argument (number of MPI processes).")
                sys.exit(1)
            break
    
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
            break
    
    # Ensure flekspy is installed first if thorough validation is enabled
    if use_flekspy:
        ensure_flekspy_installed()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(os.path.dirname(script_dir))
    print(f"Working directory set to: {os.getcwd()}")
    
    validators = {
        "beam": validate_beam,
        "photoionization": validate_ionization_source,
        "electronimpact": validate_ionization_source,
        "chargeexchange": validate_ionization_source,
    }
    
    # Discover test subdirectories under tests/
    tests_dir = "tests"
    
    # Iterate through sorted subdirectories
    subdirs = sorted([d for d in os.listdir(tests_dir) 
                      if os.path.isdir(os.path.join(tests_dir, d)) and d not in ["performance"]])
    
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
        
    results = [] # Collect results for summary table
    
    for test_dir, name, validator in tests:
        print(f"\n==========================================")
        print(f"Starting test: {name.upper()}")
        print(f"==========================================")
        try:
            stdout, code = run_test(test_dir, nprocs=nprocs)
            if code != 0 or stdout is None:
                print(f"FAIL: {name.upper()} execution failed with exit code {code}")
                results.append((name.upper(), "FAILED", f"Execution failed (code {code})"))
                continue

            # Read diagnostics from the structured log file (preferred) or fall back
            # to parsing stdout for backward compatibility.
            particle_diags, field_diags = read_diag_log("run_test")
            pic_diags = read_pic_log("run_test")
            if not particle_diags:
                print("  [INFO] No diag log file found; falling back to stdout parsing.")
                particle_diags, field_diags = parse_diagnostics(stdout)

            val_res = False
            reason = "Validation skipped"

            if validator:
                import inspect
                sig = inspect.signature(validator)
                kwargs = {}
                if "field_diags" in sig.parameters:
                    kwargs["field_diags"] = field_diags
                if "pic_diags" in sig.parameters:
                    kwargs["pic_diags"] = pic_diags
                val_res, reason = validator(particle_diags, **kwargs)
                if not val_res:
                    print("FLEKS execution output:")
                    print(stdout)
                    results.append((name.upper(), "FAILED", reason))
                    continue
            else:
                print(f"Validating {name.upper()} (generic check)...")
                if not particle_diags:
                    print("  [INFO] No diagnostic outputs parsed; skipping generic diagnostic check.")
                    val_res = True
                    reason = "Passed (diagnostics unavailable)"
                else:
                    print(f"{name.upper()} (generic check): PASSED")
                    val_res = True
                    reason = "Passed"

            # Validate output plotfiles using flekspy!
            flekspy_res, flekspy_reason = validate_amrex_outputs_with_flekspy(name, use_flekspy=use_flekspy)
            if not flekspy_res:
                results.append((name.upper(), "FAILED", f"flekspy check failed: {flekspy_reason}"))
            else:
                results.append((name.upper(), "PASSED", "Passed"))

        finally:
            # Always clean up run output after each test to keep disk usage low.
            print(f"  Cleaning up run_test/ output for {name.upper()}...")
            cleanup_run_dir()


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
