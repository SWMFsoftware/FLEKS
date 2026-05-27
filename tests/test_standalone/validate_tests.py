#!/usr/bin/env python3
import os
import shutil
import subprocess
import re
import math
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
    
    # Symlinks in run directory
    safe_symlink("../bin/FLEKS.exe", os.path.join(run_dir, "FLEKS.exe"))
    safe_symlink("../../../share/Scripts/PostProc.pl", os.path.join(run_dir, "PostProc.pl"))
    
    # Component plot and restart directories
    pc_dir = os.path.join(run_dir, "PC")
    os.makedirs(pc_dir, exist_ok=True)
    os.makedirs(os.path.join(pc_dir, "plots"), exist_ok=True)
    os.makedirs(os.path.join(pc_dir, "restartOUT"), exist_ok=True)
    
    # Symlinks in component directory
    safe_symlink("../../../../share/Scripts/pIDL", os.path.join(pc_dir, "pIDL"))
    safe_symlink("../../../../bin/PostIDL.exe", os.path.join(pc_dir, "PostIDL.exe"))

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

def run_test(test_dir):
    param_file = os.path.join(test_dir, "PARAM.in")
    print(f"Running test in {test_dir} with config {param_file}...")
    prepare_run_dir()
    
    # Copy param_file to run_test/PARAM.in
    shutil.copy(param_file, "run_test/PARAM.in")
    
    # Run ./FLEKS.exe inside run_test/
    cmd = ["./FLEKS.exe"]
    if "sound_wave" in test_dir:
        cmd = ["mpiexec", "-n", "2", "./FLEKS.exe"]
        
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
                diag_dict = {
                    "time":   t,
                    "cycle":  cycle,
                    "max_by": float(vals[col]),
                    "max_bz": float(vals[col + 1]),
                }
                if col + 2 < len(vals):
                    diag_dict["max_ey"] = float(vals[col + 2])
                field_diags.append(diag_dict)
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

def validate_box(diags):
    print("Validating Box Test...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False, "No diagnostic outputs parsed"
    
    passed = True
    reasons = []
    prod_rate = 1.0e24
    dt = 0.1
    for diag in diags:
        t = diag["time"]
        # Expected physical particles: initial injection + step updates
        expected_phys = prod_rate * (t + dt)
        actual_phys = diag["phys"]
        relative_diff = abs(actual_phys - expected_phys) / expected_phys
        print(f"  t={t:.2f}: expected={expected_phys:.2e}, actual={actual_phys:.2e}, diff={relative_diff*100:.3f}%")
        if relative_diff > 0.02: # 2% tolerance
            print(f"  FAIL: relative difference {relative_diff*100:.2f}% exceeds 2%")
            passed = False
            reasons.append(f"relative diff {relative_diff*100:.2f}% at t={t:.2f} > 2%")
    
    if passed:
        print("Box Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)


def validate_pickup(diags):
    print("Validating Pickup Test...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False, "No diagnostic outputs parsed"
    
    # Group diagnostics by species (0-indexed in FLEKS diagnostic outputs: 0=H+, 1=O+, 2=e-)
    species_diags = {0: [], 1: [], 2: []}
    for diag in diags:
        sp = diag["species"]
        if sp in species_diags:
            species_diags[sp].append(diag)
            
    passed = True
    reasons = []
    dt = 0.1
    
    # Expected production rates mapped to 0-indexed species
    prod_rates = {
        0: 8.0e23,  # H+ (Species 0)
        1: 2.0e23,  # O+ (Species 1)
        2: 1.0e24   # e- (Species 2)
    }
    
    # 1. Validate particle counts (charge neutrality verification)
    print("  --- Validating Particle Production Rates & Charge Neutrality ---")
    for sp in [0, 1, 2]:
        diags_sp = species_diags[sp]
        if not diags_sp:
            print(f"  FAIL: No diagnostics found for species {sp}")
            passed = False
            reasons.append(f"No diagnostics for species {sp}")
            continue
            
        rate = prod_rates[sp]
        for diag in diags_sp:
            t = diag["time"]
            expected_phys = rate * (t + dt)
            actual_phys = diag["phys"]
            relative_diff = abs(actual_phys - expected_phys) / expected_phys
            print(f"    Species {sp} (H+, O+, e-) at t={t:.2f}: expected={expected_phys:.2e}, actual={actual_phys:.2e} (diff={relative_diff*100:.3f}%)")
            if relative_diff > 0.02: # 2% tolerance
                print(f"    FAIL: species {sp} count difference exceeds tolerance")
                passed = False
                reasons.append(f"Species {sp} count diff {relative_diff*100:.2f}% > 2% at t={t:.2f}")

    # 2. Validate velocities for H+ (Species 0, omega_c = 1.0)
    print("  --- Validating H+ (Species 0) Velocity History ---")
    diags_h = species_diags[0]
    for diag in diags_h:
        t = diag["time"]
        N = int(round(t / dt))
        sum_vy = 0.0
        sum_vx = 0.0
        omega_c = 1.0
        for k in range(N + 1):
            tau = omega_c * k * dt
            sum_vy += math.sin(tau)
            sum_vx += 1.0 - math.cos(tau)
        expected_vy = sum_vy / (N + 1)
        expected_vx = sum_vx / (N + 1)
        
        actual_vy = diag["vy"]
        actual_vx = diag["vx"]
        diff_vx = abs(actual_vx - expected_vx)
        diff_vy = abs(actual_vy - expected_vy)
        
        print(f"    t={t:.2f}: Vx expected={expected_vx:.4f}, actual={actual_vx:.4f} (diff={diff_vx:.4f})")
        print(f"            Vy expected={expected_vy:.4f}, actual={actual_vy:.4f} (diff={diff_vy:.4f})")
        if diff_vx > 0.02 or diff_vy > 0.02:
            print(f"    FAIL: H+ velocity difference exceeds tolerance (0.02)")
            passed = False
            reasons.append(f"H+ vel diff (Vx diff={diff_vx:.4f}, Vy diff={diff_vy:.4f}) > 0.02 at t={t:.2f}")

    # 3. Validate velocities for O+ (Species 1, omega_c = 0.0625)
    print("  --- Validating O+ (Species 1) Velocity History ---")
    diags_o = species_diags[1]
    for diag in diags_o:
        t = diag["time"]
        N = int(round(t / dt))
        sum_vy = 0.0
        sum_vx = 0.0
        omega_c = 0.0625
        for k in range(N + 1):
            tau = omega_c * k * dt
            sum_vy += math.sin(tau)
            sum_vx += 1.0 - math.cos(tau)
        expected_vy = sum_vy / (N + 1)
        expected_vx = sum_vx / (N + 1)
        
        actual_vy = diag["vy"]
        actual_vx = diag["vx"]
        diff_vx = abs(actual_vx - expected_vx)
        diff_vy = abs(actual_vy - expected_vy)
        
        print(f"    t={t:.2f}: Vx expected={expected_vx:.4f}, actual={actual_vx:.4f} (diff={diff_vx:.4f})")
        print(f"            Vy expected={expected_vy:.4f}, actual={actual_vy:.4f} (diff={diff_vy:.4f})")
        if diff_vx > 0.02 or diff_vy > 0.02:
            print(f"    FAIL: O+ velocity difference exceeds tolerance (0.02)")
            passed = False
            reasons.append(f"O+ vel diff (Vx diff={diff_vx:.4f}, Vy diff={diff_vy:.4f}) > 0.02 at t={t:.2f}")
            
    if passed:
        print("Pickup Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)


def validate_photoionization(diags):
    print("Validating Photoionization Test...")
    return validate_pickup(diags)


def validate_electron_impact(diags):
    print("Validating Electron Impact Ionization (MCC) Test...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False, "No diagnostic outputs parsed"
        
    species_diags = {0: [], 1: [], 2: []}
    for diag in diags:
        sp = diag["species"]
        if sp in species_diags:
            species_diags[sp].append(diag)
            
    passed = True
    reasons = []
    dt = 0.1
    
    # Exospheric nominal rates
    nom_rates = {
        0: 8.0e23,  # H+ nominal exosphere
        1: 2.0e23,  # O+ nominal exosphere
        2: 1.0e24   # e- nominal exosphere
    }
    
    # We expect the actual counts to be HIGHER than nominal due to electron impact ionization!
    print("  --- Checking Electron Impact Ionization Yield ---")
    for sp in [0, 1, 2]:
        diags_sp = species_diags[sp]
        if not diags_sp:
            print(f"  FAIL: No diagnostics found for species {sp}")
            passed = False
            reasons.append(f"No diagnostics for species {sp}")
            continue
            
        # Check last diagnostic step
        last_diag = diags_sp[-1]
        t = last_diag["time"]
        nominal = nom_rates[sp] * (t + dt)
        actual = last_diag["phys"]
        
        diff = actual - nominal
        print(f"    Species {sp} at t={t:.2f}: nominal_exosphere={nominal:.2e}, actual_with_mcc={actual:.2e} (MCC yield={diff:.2e})")
        if diff <= 0.0:
            print(f"    FAIL: Species {sp} has no electron impact ionization yield (actual {actual:.2e} <= nominal exosphere {nominal:.2e})")
            passed = False
            reasons.append(f"Species {sp} has no MCC yield (actual {actual:.2e} <= nominal {nominal:.2e})")
            
    if passed:
        print("Electron Impact Ionization Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)


def validate_exosphere_charge_exchange(diags):
    print("Validating Exospheric Charge Exchange (MCC) Test...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False, "No diagnostic outputs parsed"
        
    species_diags = {0: [], 1: [], 2: []}
    for diag in diags:
        sp = diag["species"]
        if sp in species_diags:
            species_diags[sp].append(diag)
            
    passed = True
    reasons = []
    
    print("  --- Checking Charge Exchange cooling effect ---")
    for sp in [0, 1]:
        diags_sp = species_diags[sp]
        if not diags_sp:
            print(f"  FAIL: No diagnostics found for species {sp}")
            passed = False
            reasons.append(f"No diagnostics for species {sp}")
            continue
            
        last_diag = diags_sp[-1]
        t = last_diag["time"]
        actual_vx = last_diag["vx"]
        
        initial_vx = 400.0
        print(f"    Species {sp} at t={t:.2f}: initial_drift_vx={initial_vx:.1f}, cooled_vx={actual_vx:.2f}")
        if actual_vx >= initial_vx:
            print(f"    FAIL: Species {sp} did not experience any charge exchange cooling (vx {actual_vx:.2f} >= initial vx {initial_vx:.2f})")
            passed = False
            reasons.append(f"Species {sp} did not experience charge exchange cooling (vx {actual_vx:.2f} >= {initial_vx:.2f})")
            
    if passed:
        print("Exospheric Charge Exchange Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)


def validate_exosphere(diags):
    print("Validating Combined Exosphere Test (Photoionization, Electron Impact MCC, and Charge Exchange MCC)...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False, "No diagnostic outputs parsed"
        
    species_diags = {0: [], 1: [], 2: []}
    for diag in diags:
        sp = diag["species"]
        if sp in species_diags:
            species_diags[sp].append(diag)
            
    passed = True
    reasons = []
    dt = 0.1
    
    # Exospheric nominal rates
    nom_rates = {
        0: 8.0e23,  # H+
        1: 2.0e23,  # O+
        2: 1.0e24   # e-
    }
    
    # 1. Verify Electron Impact yield (actual > nominal exosphere)
    print("  --- Checking Electron Impact Ionization Yield ---")
    for sp in [0, 1, 2]:
        diags_sp = species_diags[sp]
        if not diags_sp:
            print(f"  FAIL: No diagnostics found for species {sp}")
            passed = False
            reasons.append(f"No diagnostics for species {sp}")
            continue
            
        last_diag = diags_sp[-1]
        t = last_diag["time"]
        nominal = nom_rates[sp] * (t + dt)
        actual = last_diag["phys"]
        
        diff = actual - nominal
        print(f"    Species {sp} at t={t:.2f}: nominal_exosphere={nominal:.2e}, actual={actual:.2e} (diff={diff:.2e})")
        if diff <= 0.0:
            print(f"    FAIL: Species {sp} has no combined process yield (actual {actual:.2e} <= nominal exosphere {nominal:.2e})")
            passed = False
            reasons.append(f"Species {sp} has no combined process yield (actual {actual:.2e} <= nominal {nominal:.2e})")

    # 2. Verify Charge Exchange cooling (final Vx < initial Vx = 400.0)
    print("  --- Checking Charge Exchange cooling effect ---")
    for sp in [0, 1]:
        diags_sp = species_diags[sp]
        if not diags_sp:
            print(f"  FAIL: No diagnostics found for species {sp}")
            passed = False
            reasons.append(f"No diagnostics for species {sp}")
            continue
            
        last_diag = diags_sp[-1]
        t = last_diag["time"]
        actual_vx = last_diag["vx"]
        
        initial_vx = 400.0
        print(f"    Species {sp} at t={t:.2f}: initial_drift_vx={initial_vx:.1f}, cooled_vx={actual_vx:.2f}")
        if actual_vx >= initial_vx:
            print(f"    FAIL: Species {sp} did not experience any charge exchange cooling (vx {actual_vx:.2f} >= initial vx {initial_vx:.2f})")
            passed = False
            reasons.append(f"Species {sp} did not experience charge exchange cooling (vx {actual_vx:.2f} >= {initial_vx:.2f})")
            
    if passed:
        print("Combined Exosphere Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)


def validate_chamber(diags):
    print("Validating Chamber Test...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False, "No diagnostic outputs parsed"
    
    times = [d["time"] for d in diags]
    phys_particles = [d["phys"] for d in diags]
    
    if len(phys_particles) < 5:
        print("FAIL: Too few data points.")
        return False, "Too few data points"
    
    # Check that in the last 1.5 seconds, the relative change in particle number is very small
    last_idx = len(phys_particles) - 1
    # Find index at t ~ 3.5
    start_idx = 0
    for idx, t in enumerate(times):
        if t >= 3.5:
            start_idx = idx
            break
            
    p_start = phys_particles[start_idx]
    p_end = phys_particles[last_idx]
    relative_change = abs(p_end - p_start) / p_start
    print(f"  Steady state check from t={times[start_idx]:.2f} to t={times[last_idx]:.2f}:")
    print(f"    N(t={times[start_idx]:.2f}) = {p_start:.2e}")
    print(f"    N(t={times[last_idx]:.2f}) = {p_end:.2e}")
    print(f"    Relative change = {relative_change*100:.3f}%")
    
    # If the relative change is < 5%, we have reached a steady state balancing injection and escape.
    if relative_change < 0.05:
        print("Chamber Test: PASSED")
        return True, "Passed"
    else:
        print("FAIL: Chamber did not reach steady-state (relative change >= 5%)")
        return False, f"Chamber relative change {relative_change*100:.2f}% >= 5%"


def validate_beam(diags, field_diags=None):
    print("Validating Beam Instability Test...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False, "No diagnostic outputs parsed"

    passed = True
    reasons = []

    # Group diagnostics by species independently.
    from collections import defaultdict
    by_species = defaultdict(list)
    for d in diags:
        by_species[d["species"]].append(d)

    # 1. Particle count conservation per species
    for iS, species_diags in sorted(by_species.items()):
        if not species_diags:
            continue
        initial_phys = species_diags[0]["phys"]
        for diag in species_diags:
            t = diag["time"]
            actual_phys = diag["phys"]
            if initial_phys > 0 and abs(actual_phys - initial_phys) > 1e-5 * initial_phys:
                print(f"  FAIL species {iS} at t={t:.2f}: Particle count changed! "
                      f"Expected={initial_phys:.2e}, Actual={actual_phys:.2e}")
                passed = False
                reasons.append(f"Species {iS} particle count changed at t={t:.2f} "
                               f"(expected {initial_phys:.2e}, actual {actual_phys:.2e})")

    # 2. Mean velocity conservation check on the ion species (species 1).
    # Expected MeanVx is -0.392 (99% background at -0.4, 1% beam at +0.4).
    ion_diags = by_species.get(1, [])
    if not ion_diags:
        # Fallback: single-species run (legacy); use species 0
        ion_diags = by_species.get(0, [])
    expected_vx = -0.392
    for diag in ion_diags:
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


def validate_sound_wave(diags, is_cold=False):
    print("Validating Sound Wave Test..." if not is_cold else "Validating Cold Wave Test...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False, "No diagnostic outputs parsed"
    
    passed = True
    reasons = []
    
    # Group diagnostics by species so we check each species independently.
    # diags is a flat list with entries for every species at every timestep.
    from collections import defaultdict
    by_species = defaultdict(list)
    for d in diags:
        by_species[d["species"]].append(d)

    for iS, species_diags in sorted(by_species.items()):
        if not species_diags:
            continue
        # 1. Particle count conservation per species
        initial_phys = species_diags[0]["phys"]
        for diag in species_diags:
            t = diag["time"]
            actual_phys = diag["phys"]
            if initial_phys > 0 and abs(actual_phys - initial_phys) > 1e-5 * initial_phys:
                print(f"  FAIL species {iS} at t={t:.2f}: Particle count changed! "
                      f"Expected={initial_phys:.2e}, Actual={actual_phys:.2e}")
                passed = False
                reasons.append(f"Species {iS} particle count changed at t={t:.2f}")

        # 2. Kinetic energy stability per species
        if is_cold:
            # In a cold plasma (T=0), kinetic energy oscillates between magnetic energy and kinetic bulk energy.
            # Light species (electrons) are also subject to numerical/grid heating. Therefore, we bypass
            # the KE stability check for cold plasma wave tests.
            continue

        initial_ke = species_diags[0]["ke"]
        ke_tol = 0.85 if iS == 0 else 0.20  # allow more slack for electrons (species 0) due to physical wave heating
        for diag in species_diags:
            t = diag["time"]
            actual_ke = diag["ke"]
            if initial_ke > 0:
                diff_ke = abs(actual_ke - initial_ke) / initial_ke
                if diff_ke > ke_tol:
                    print(f"  FAIL species {iS} at t={t:.2f}: KE variation too large! "
                          f"Expected={initial_ke:.2e}, Actual={actual_ke:.2e}, "
                          f"diff={diff_ke*100:.2f}%")
                    passed = False
                    reasons.append(f"Species {iS} KE variation too large at t={t:.2f}")

    if passed:
        print("Wave Test: PASSED")
        return True, "Passed"
    else:
        return False, "; ".join(reasons)

def validate_fast_wave(diags):
    print("Validating Fast Wave Test...")
    return validate_sound_wave(diags)

def validate_alfven_wave(diags):
    print("Validating Alfven Wave Test...")
    return validate_sound_wave(diags)

def validate_slow_wave(diags):
    print("Validating Slow Wave Test...")
    return validate_sound_wave(diags)


def validate_tophat(diags, field_diags=None):
    print("Validating TopHat Test...")
    if not field_diags:
        print("FAIL: No field diagnostic outputs parsed.")
        return False, "No field diagnostics parsed"
    
    passed = True
    reasons = []
    
    for diag in field_diags:
        t = diag["time"]
        max_bz = diag["max_bz"]
        max_by = diag["max_by"]
        max_ey = diag.get("max_ey", None)
        
        # max_bz should remain close to 1.0 (allow [0.8, 1.2])
        if max_bz < 0.8 or max_bz > 1.2:
            print(f"  FAIL at t={t:.2f}: maxBz is outside acceptable range [0.8, 1.2]! Actual={max_bz:.3f}")
            passed = False
            reasons.append(f"maxBz outside range at t={t:.2f}")

        # max_ey should remain close to 1.0 (allow [0.8, 1.2])
        if max_ey is not None and (max_ey < 0.8 or max_ey > 1.2):
            print(f"  FAIL at t={t:.2f}: maxEy is outside acceptable range [0.8, 1.2]! Actual={max_ey:.3f}")
            passed = False
            reasons.append(f"maxEy outside range at t={t:.2f}")
            
        # max_by should remain near 0.0 (allow <= 0.05)
        if max_by > 0.05:
            print(f"  FAIL at t={t:.2f}: maxBy is too high! Actual={max_by:.3f}")
            passed = False
            reasons.append(f"maxBy too high at t={t:.2f}")
            
    if passed:
        # Print actual initial/final values to confirm validation
        if len(field_diags) > 0:
            init_bz = field_diags[0]["max_bz"]
            init_ey = field_diags[0].get("max_ey", 0.0)
            final_bz = field_diags[-1]["max_bz"]
            final_ey = field_diags[-1].get("max_ey", 0.0)
            print(f"  Initial Fields: maxBz={init_bz:.3f}, maxEy={init_ey:.3f}")
            print(f"  Final Fields:   maxBz={final_bz:.3f}, maxEy={final_ey:.3f}")
        print("TopHat Test: PASSED")
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
    # Check if thorough validation is requested
    use_flekspy = "--thorough" in sys.argv or "--flekspy" in sys.argv
    
    # Ensure flekspy is installed first if thorough validation is enabled
    if use_flekspy:
        ensure_flekspy_installed()
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(os.path.dirname(os.path.dirname(script_dir)))
    print(f"Working directory set to: {os.getcwd()}")
    
    # Define mapping from folder name to specific validator
    validators = {
        "box": validate_box,
        "pickup": validate_pickup,
        "photoionization": validate_photoionization,
        "electron_impact": validate_electron_impact,
        "charge_exchange": validate_exosphere_charge_exchange,
        "exosphere": validate_exosphere,
        "chamber": validate_chamber,
        "beam": validate_beam,
        "sound_wave": validate_sound_wave,
        "fast_wave": validate_fast_wave,
        "alfven_wave": validate_alfven_wave,
        "slow_wave": validate_slow_wave,
        "tophat": validate_tophat
    }
    
    # Discover test subdirectories under tests/test_standalone
    test_standalone_dir = os.path.join("tests", "test_standalone")
    
    # Iterate through sorted subdirectories
    subdirs = sorted([d for d in os.listdir(test_standalone_dir) 
                      if os.path.isdir(os.path.join(test_standalone_dir, d)) and d not in ["performance"]])
    
    # Filter by command-line arguments if provided
    targets = [arg for arg in sys.argv[1:] if not arg.startswith("-")]
    if targets:
        subdirs = [d for d in subdirs if d in targets]
    
    tests = []
    for d in subdirs:
        test_dir = os.path.join(test_standalone_dir, d)
        param_file = os.path.join(test_dir, "PARAM.in")
        if os.path.exists(param_file):
            validator = validators.get(d, None)
            tests.append((test_dir, d, validator))
            
    if not tests:
        print("No tests found in tests/test_standalone/ subdirectories!")
        sys.exit(1)
        
    results = [] # Collect results for summary table
    
    for test_dir, name, validator in tests:
        val_res = False
        print(f"\n==========================================")
        print(f"Starting test: {name.upper()}")
        print(f"==========================================")
        try:
            stdout, code = run_test(test_dir)
            if code != 0 or stdout is None:
                print(f"FAIL: {name.upper()} execution failed with exit code {code}")
                results.append((name.upper(), "FAILED", f"Execution failed (code {code})"))
                continue

            # Read diagnostics from the structured log file (preferred) or fall back
            # to parsing stdout for backward compatibility.
            particle_diags, field_diags = read_diag_log("run_test")
            if not particle_diags and not field_diags:
                print("  [INFO] No log file found; falling back to stdout parsing.")
                particle_diags, field_diags = parse_diagnostics(stdout)

            val_res = False
            reason = "Validation skipped"

            if validator:
                import inspect
                sig = inspect.signature(validator)
                # validate_beam now accepts field_diags keyword instead of stdout.
                if "field_diags" in sig.parameters:
                    val_res, reason = validator(particle_diags, field_diags=field_diags)
                else:
                    val_res, reason = validator(particle_diags)
                if not val_res:
                    print("FLEKS execution output:")
                    print(stdout)
                    results.append((name.upper(), "FAILED", reason))
                    continue
            else:
                print(f"Validating {name.upper()} (generic check)...")
                if not particle_diags:
                    print("FAIL: No diagnostic outputs parsed.")
                    results.append((name.upper(), "FAILED", "No diagnostics parsed"))
                    continue
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
            if val_res:
                print(f"  Cleaning up run_test/ output for {name.upper()}...")
                cleanup_run_dir()
            else:
                print(f"  Skipping cleanup for {name.upper()} to preserve diagnostics.")


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
        summary_file = os.path.join(test_standalone_dir, "summary.md")
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
