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

def run_test(test_dir):
    param_file = os.path.join(test_dir, "PARAM.in")
    print(f"Running test in {test_dir} with config {param_file}...")
    prepare_run_dir()
    
    # Copy param_file to run_test/PARAM.in
    shutil.copy(param_file, "run_test/PARAM.in")
    
    # Run ./FLEKS.exe inside run_test/
    result = subprocess.run(["./FLEKS.exe"], cwd="run_test", stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error running FLEKS.exe for {test_dir}:")
        print(result.stderr)
        return None, result.returncode
        
    # Automatically run post-processing on the generated plots
    subprocess.run(["./PostProc.pl", "-v"], cwd="run_test", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    return result.stdout, 0


def parse_diagnostics(stdout):
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
                "time": float(match.group(2)),
                "cycle": int(match.group(3)),
                "macro": int(match.group(4)),
                "phys": float(match.group(5)),
                "vx": float(match.group(6)),
                "vy": float(match.group(7)),
                "vz": float(match.group(8)),
                "ke": float(match.group(9))
            }
            diagnostics.append(diag)
    return diagnostics

def validate_box(diags):
    print("Validating Box Test...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False
    
    passed = True
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
    
    if passed:
        print("Box Test: PASSED")
    return passed

def validate_pickup(diags):
    print("Validating Pickup Test...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False
    
    # Group diagnostics by species (0-indexed in FLEKS diagnostic outputs: 0=H+, 1=O+, 2=e-)
    species_diags = {0: [], 1: [], 2: []}
    for diag in diags:
        sp = diag["species"]
        if sp in species_diags:
            species_diags[sp].append(diag)
            
    passed = True
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
            
    if passed:
        print("Pickup Test: PASSED")
    return passed


def validate_chamber(diags):
    print("Validating Chamber Test...")
    if not diags:
        print("FAIL: No diagnostic outputs parsed.")
        return False
    
    times = [d["time"] for d in diags]
    phys_particles = [d["phys"] for d in diags]
    
    if len(phys_particles) < 5:
        print("FAIL: Too few data points.")
        return False
    
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
        return True
    else:
        print("FAIL: Chamber did not reach steady-state (relative change >= 5%)")
        return False

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(os.path.dirname(os.path.dirname(script_dir)))
    print(f"Working directory set to: {os.getcwd()}")
    
    # Define mapping from folder name to specific validator
    validators = {
        "box": validate_box,
        "pickup": validate_pickup,
        "chamber": validate_chamber
    }
    
    # Discover test subdirectories under tests/test_standalone
    test_standalone_dir = os.path.join("tests", "test_standalone")
    
    # Iterate through sorted subdirectories (excluding exosphere which is merged into pickup)
    subdirs = sorted([d for d in os.listdir(test_standalone_dir) 
                      if os.path.isdir(os.path.join(test_standalone_dir, d)) and d != "exosphere"])
    
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
        
    all_passed = True
    for test_dir, name, validator in tests:
        print(f"\n==========================================")
        print(f"Starting test: {name.upper()}")
        print(f"==========================================")
        stdout, code = run_test(test_dir)
        if code != 0 or stdout is None:
            print(f"FAIL: {name.upper()} execution failed with exit code {code}")
            all_passed = False
            continue
            
        diags = parse_diagnostics(stdout)
        if validator:
            if not validator(diags):
                print("FLEKS execution output:")
                print(stdout)
                all_passed = False
        else:
            print(f"Validating {name.upper()} (generic check)...")
            if not diags:
                print("FAIL: No diagnostic outputs parsed.")
                all_passed = False
            else:
                print(f"{name.upper()} (generic check): PASSED")
                
    if all_passed:
        print("\nALL STANDALONE EXOSPHERE TESTS PASSED SUCCESSFULLY!")
        sys.exit(0)
    else:
        print("\nSOME TESTS FAILED.")
        sys.exit(1)

if __name__ == "__main__":
    main()
