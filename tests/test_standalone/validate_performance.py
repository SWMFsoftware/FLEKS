#!/usr/bin/env python3
import os
import shutil
import subprocess
import re
import sys
import time

def safe_symlink(src, dst):
    if os.path.lexists(dst):
        if os.path.islink(dst) or os.path.isfile(dst):
            os.remove(dst)
        elif os.path.isdir(dst):
            shutil.rmtree(dst)
    os.symlink(src, dst)

def prepare_run_dir(run_dir):
    os.makedirs(run_dir, exist_ok=True)
    safe_symlink(os.path.abspath("../../bin/FLEKS.exe"), os.path.join(run_dir, "FLEKS.exe"))
    
    # Component plots and restart directories
    pc_dir = os.path.join(run_dir, "PC")
    os.makedirs(pc_dir, exist_ok=True)
    os.makedirs(os.path.join(pc_dir, "plots"), exist_ok=True)
    os.makedirs(os.path.join(pc_dir, "restartOUT"), exist_ok=True)

def run_benchmark(n_proc, run_dir):
    print(f"Running benchmark in {run_dir} with {n_proc} MPI process(es)...")
    prepare_run_dir(run_dir)
    
    # Copy benchmark config
    shutil.copy("performance/PARAM.in", os.path.join(run_dir, "PARAM.in"))
    
    # Setup MPI command line
    cmd = ["mpirun", "-n", str(n_proc), "./FLEKS.exe"]
    
    start_time = time.perf_counter()
    result = subprocess.run(cmd, cwd=run_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    end_time = time.perf_counter()
    
    duration = end_time - start_time
    
    if result.returncode != 0:
        print(f"Error executing benchmark with {n_proc} processes:")
        print(result.stderr)
        return None, result.returncode
        
    # Parse total cycles and particles from stdout
    cycles = 0
    particles = 0
    
    # We parse the last DIAGNOSTIC output
    pattern = re.compile(r"DIAGNOSTIC:\s+Species=\d+\s+Time=[+\-\d.e]+\s+Cycle=(\d+)\s+MacroParticles=(\d+)")
    for line in result.stdout.splitlines():
        match = pattern.search(line)
        if match:
            cycles = int(match.group(1))
            particles = int(match.group(2))
            
    return {
        "duration": duration,
        "cycles": cycles,
        "particles": particles,
        "stdout": result.stdout
    }, 0

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    if not os.path.exists("../../bin/FLEKS.exe"):
        print("Error: Standalone FLEKS.exe not found. Please compile it first.")
        sys.exit(1)
        
    print("================================================================================")
    print("                      FLEKS PERFORMANCE REGRESSION TEST")
    print("================================================================================")
    
    # 1. Serial Run (1 MPI process)
    serial_res, code = run_benchmark(1, "run_perf_serial")
    if code != 0 or serial_res is None:
        print("FAIL: Serial benchmark run failed.")
        sys.exit(1)
        
    t_serial = serial_res["duration"]
    p_count = serial_res["particles"]
    cycles = serial_res["cycles"]
    
    if p_count == 0 or cycles == 0:
        print("Error: Could not parse diagnostic cycles/particles from serial output.")
        sys.exit(1)
        
    # Calculate Particle-Step Mover Time (PPS in microseconds)
    # PPS = (T * 1e6) / (Particles * Cycles)
    pps_serial = (t_serial * 1e6) / (p_count * cycles)
    print(f"Serial Timing: {t_serial:.3f} s, Particles: {p_count}, Cycles: {cycles}, Rate: {pps_serial:.3f} μs/part-step")
    
    # 2. Parallel Run (2 MPI processes)
    parallel_res, code = run_benchmark(2, "run_perf_parallel")
    if code != 0 or parallel_res is None:
        print("FAIL: Parallel benchmark run failed.")
        sys.exit(1)
        
    t_parallel = parallel_res["duration"]
    pps_parallel = (t_parallel * 1e6) / (p_count * cycles)
    speedup = t_serial / t_parallel
    print(f"Parallel Timing: {t_parallel:.3f} s, Rate: {pps_parallel:.3f} μs/part-step, Speedup: {speedup:.2f}x")
    
    # 3. Baselines verification
    serial_pps_baseline = 15.0  # 15 microseconds per particle-step is a very safe limit for virtualized runners
    speedup_baseline = 1.25     # Perfect scaling is 2.0x, but 1.25x is highly robust for virtual environments
    
    serial_passed = pps_serial <= serial_pps_baseline
    parallel_passed = speedup >= speedup_baseline
    
    print("\n" + "=" * 80)
    print(" " * 30 + "PERFORMANCE COMPARISON")
    print("=" * 80)
    print(f" {'Metric':<25} | {'Measured':<15} | {'Baseline Target':<18} | {'Status':<8}")
    print("-" * 80)
    
    status_serial = "PASSED" if serial_passed else "FAILED"
    status_parallel = "PASSED" if parallel_passed else "FAILED"
    
    status_serial_display = status_serial
    status_parallel_display = status_parallel
    
    if sys.stdout.isatty():
        status_serial_display = f"\033[92;1m{status_serial}\033[0m" if serial_passed else f"\033[91;1m{status_serial}\033[0m"
        status_parallel_display = f"\033[92;1m{status_parallel}\033[0m" if parallel_passed else f"\033[91;1m{status_parallel}\033[0m"
        
    print(f" {'Serial Rate (PPS)':<25} | {pps_serial:<15.3f} | <={serial_pps_baseline:<14.1f} μs | {status_serial_display}")
    print(f" {'Parallel Speedup':<25} | {speedup:<15.2f} | >={speedup_baseline:<14.2f} x  | {status_parallel_display}")
    print("=" * 80)
    
    # 4. Generate markdown report
    summary_path = "performance_summary.md"
    try:
        with open(summary_path, "w") as f:
            f.write("### ⚡ Standalone FLEKS Performance Report\n\n")
            f.write("| Performance Metric | Measured | Target Baseline | Status |\n")
            f.write("| :--- | :--- | :--- | :--- |\n")
            
            serial_status_md = "🟢 **PASSED**" if serial_passed else "🔴 **FAILED**"
            parallel_status_md = "🟢 **PASSED**" if parallel_passed else "🔴 **FAILED**"
            
            f.write(f"| Serial Rate (μs/part-step) | {pps_serial:.3f} μs | <= {serial_pps_baseline:.1f} μs | {serial_status_md} |\n")
            f.write(f"| Parallel Speedup (2 Cores) | {speedup:.2f}x | >= {speedup_baseline:.2f}x | {parallel_status_md} |\n\n")
            f.write(f"*Note: Test run on virtualized GitHub Actions runner. Base particles: {p_count}, cycles: {cycles}.*\n")
    except Exception as e:
        print(f"Warning: Could not write performance_summary.md: {e}")
        
    if serial_passed and parallel_passed:
        if sys.stdout.isatty():
            print("\033[92;1m\nALL PERFORMANCE TESTS PASSED SUCCESSFULLY!\033[0m\n")
        else:
            print("\nALL PERFORMANCE TESTS PASSED SUCCESSFULLY!\n")
        sys.exit(0)
    else:
        if sys.stdout.isatty():
            print("\033[91;1m\nPERFORMANCE REGRESSION DETECTED.\033[0m\n")
        else:
            print("\nPERFORMANCE REGRESSION DETECTED.\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
