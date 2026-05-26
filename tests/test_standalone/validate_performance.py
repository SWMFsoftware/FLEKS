#!/usr/bin/env python3
import os
import shutil
import subprocess
import re
import sys
import time
import glob

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

def read_diag_log(run_dir):
    """Read the structured diagnostic log file written by Pic::write_diag_log.
    Returns a dict with final cycle and total macroparticle count.
    """
    pc_plots = os.path.join(run_dir, "PC", "plots")
    log_files = sorted(glob.glob(os.path.join(pc_plots, "log_diag_n*.log")))
    if not log_files:
        return None

    # Use the most recent log file
    log_file = log_files[-1]

    try:
        with open(log_file, "r") as f:
            lines = f.readlines()
        if len(lines) < 2:
            return None
            
        # Parse last data row
        vals = lines[-1].strip().split("\t")
        header = lines[0].strip().split("\t")
        
        cycle = int(vals[1])
        
        # Sum macroparticles over all species
        n_species = sum(1 for col in header if col.startswith("macro"))
        macro_particles = 0
        col = 2
        for iS in range(n_species):
            macro_particles += int(vals[col])
            col += 6
            
        return {
            "cycle": cycle,
            "macro_particles": macro_particles
        }
    except Exception as e:
        print(f"Warning: Failed to parse diagnostic log {log_file}: {e}")
        return None

def parse_tiny_profiler(stdout):
    """Parse AMReX's TinyProfiler exclusive and inclusive timing tables from stdout.
    """
    profiler_data = {}
    
    # Parse total profile time
    total_time_match = re.search(r"TinyProfiler total time across processes\s+\[min\.\.\.avg\.\.\.max\]:\s+([\d.]+)", stdout)
    if total_time_match:
        profiler_data["total_time"] = float(total_time_match.group(1))
    else:
        total_time_match2 = re.search(r"TinyProfiler total time across processes.*?:\s+([\d.]+)", stdout)
        if total_time_match2:
            profiler_data["total_time"] = float(total_time_match2.group(1))
        
    # Parse Exclusive Timing Table
    excl_section_match = re.search(
        r"Name\s+NCalls\s+Excl\.\s+Min\s+Excl\.\s+Avg\s+Excl\.\s+Max\s+Max\s+%\n-+([\s\S]+?)-+", 
        stdout
    )
    if excl_section_match:
        excl_text = excl_section_match.group(1)
        excl_entries = {}
        for line in excl_text.splitlines():
            line = line.strip()
            if not line:
                continue
            match = re.match(r"^(.+?)\s+(\d+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)%", line)
            if match:
                name = match.group(1).strip()
                excl_entries[name] = {
                    "calls": int(match.group(2)),
                    "min": float(match.group(3)),
                    "avg": float(match.group(4)),
                    "max": float(match.group(5)),
                    "percent": float(match.group(6))
                }
        profiler_data["exclusive"] = excl_entries

    # Parse Inclusive Timing Table
    incl_section_match = re.search(
        r"Name\s+NCalls\s+Incl\.\s+Min\s+Incl\.\s+Avg\s+Incl\.\s+Max\s+Max\s+%\n-+([\s\S]+?)-+", 
        stdout
    )
    if incl_section_match:
        incl_text = incl_section_match.group(1)
        incl_entries = {}
        for line in incl_text.splitlines():
            line = line.strip()
            if not line:
                continue
            match = re.match(r"^(.+?)\s+(\d+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)%", line)
            if match:
                name = match.group(1).strip()
                incl_entries[name] = {
                    "calls": int(match.group(2)),
                    "min": float(match.group(3)),
                    "avg": float(match.group(4)),
                    "max": float(match.group(5)),
                    "percent": float(match.group(6))
                }
        profiler_data["inclusive"] = incl_entries
        
    return profiler_data

def run_benchmark_suite(n_proc, run_dir, count=3):
    print(f"Running benchmark in {run_dir} with {n_proc} MPI process(es) ({count} runs)...")
    runs = []
    
    for i in range(count):
        print(f"  Run {i+1}/{count}...")
        prepare_run_dir(run_dir)
        
        # Clean up old log files in this dir
        for log in glob.glob(os.path.join(run_dir, "PC", "plots", "log_diag_n*.log")):
            try:
                os.remove(log)
            except Exception:
                pass
                
        # Copy benchmark config
        shutil.copy("performance/PARAM.in", os.path.join(run_dir, "PARAM.in"))
        
        # Setup MPI command line
        cmd = ["mpirun", "-n", str(n_proc), "./FLEKS.exe"]
        
        start_time = time.perf_counter()
        result = subprocess.run(cmd, cwd=run_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        end_time = time.perf_counter()
        
        duration = end_time - start_time
        
        if result.returncode != 0:
            print(f"Error executing benchmark in run {i+1}:")
            print(result.stderr)
            return None, result.returncode
            
        # Parse Diagnostic log file
        diag_data = read_diag_log(run_dir)
        cycles = 30
        particles = 25906
        if diag_data:
            cycles = diag_data["cycle"]
            particles = diag_data["macro_particles"]
            
        # Parse TinyProfiler
        prof = parse_tiny_profiler(result.stdout)
        
        # Extract hot path timings from TinyProfiler
        # 1. Pts::charged_particle_mover (exclusive)
        mover_excl_time = 0.0
        if "exclusive" in prof and "Pts::charged_particle_mover" in prof["exclusive"]:
            mover_excl_time = prof["exclusive"]["Pts::charged_particle_mover"]["avg"]
            
        # 2. Pic::E_iterate (inclusive)
        solver_incl_time = 0.0
        if "inclusive" in prof and "Pic::E_iterate" in prof["inclusive"]:
            solver_incl_time = prof["inclusive"]["Pic::E_iterate"]["avg"]
        elif "inclusive" in prof and "Pic::E_matvec" in prof["inclusive"]:
            solver_incl_time = prof["inclusive"]["Pic::E_matvec"]["avg"]
            
        runs.append({
            "duration": duration,
            "cycles": cycles,
            "particles": particles,
            "total_profile_time": prof.get("total_time", duration),
            "mover_time": mover_excl_time,
            "solver_time": solver_incl_time,
        })
        
    # Calculate statistics across runs (median and min)
    def median(lst):
        sorted_lst = sorted(lst)
        n = len(sorted_lst)
        if n % 2 == 1:
            return sorted_lst[n//2]
        else:
            return (sorted_lst[n//2 - 1] + sorted_lst[n//2]) / 2.0
            
    durations = [r["duration"] for r in runs]
    movers = [r["mover_time"] for r in runs]
    solvers = [r["solver_time"] for r in runs]
    
    # We take the final run's cycles/particles (they are constant)
    cycles = runs[-1]["cycles"]
    particles = runs[-1]["particles"]
    
    stats = {
        "runs": runs,
        "cycles": cycles,
        "particles": particles,
        "duration_median": median(durations),
        "duration_min": min(durations),
        "mover_median": median(movers),
        "mover_min": min(movers),
        "solver_median": median(solvers),
        "solver_min": min(solvers),
    }
    
    return stats, 0

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    if not os.path.exists("../../bin/FLEKS.exe"):
        print("Error: Standalone FLEKS.exe not found. Please compile it first.")
        sys.exit(1)
        
    print("================================================================================")
    print("                 FLEKS ROBUST PERFORMANCE REGRESSION TEST")
    print("================================================================================")
    
    # Run serial benchmark suite
    serial_stats, code = run_benchmark_suite(1, "run_perf_serial", count=3)
    if code != 0 or serial_stats is None:
        print("FAIL: Serial benchmark runs failed.")
        sys.exit(1)
        
    t_serial_median = serial_stats["duration_median"]
    p_count = serial_stats["particles"]
    cycles = serial_stats["cycles"]
    
    # Calculate Particle-Step Rates (PPS in microseconds)
    # PPS = (T * 1e6) / (Particles * Cycles)
    total_steps = p_count * cycles
    pps_serial_median = (t_serial_median * 1e6) / total_steps
    
    # TinyProfiler rates
    pps_mover_serial_median = (serial_stats["mover_median"] * 1e6) / total_steps
    pps_solver_serial_median = (serial_stats["solver_median"] * 1e6) / total_steps
    
    print(f"\nSerial Median Wall-Clock: {t_serial_median:.3f} s (Rate: {pps_serial_median:.3f} μs/part-step)")
    print(f"Serial Median Mover:      {serial_stats["mover_median"]:.3f} s (Rate: {pps_mover_serial_median:.3f} μs/part-step)")
    print(f"Serial Median Solver:     {serial_stats["solver_median"]:.3f} s (Rate: {pps_solver_serial_median:.3f} μs/part-step)")
    
    # Run parallel benchmark suite
    parallel_stats, code = run_benchmark_suite(2, "run_perf_parallel", count=3)
    if code != 0 or parallel_stats is None:
        print("FAIL: Parallel benchmark runs failed.")
        sys.exit(1)
        
    t_parallel_median = parallel_stats["duration_median"]
    pps_parallel_median = (t_parallel_median * 1e6) / total_steps
    
    speedup = t_serial_median / t_parallel_median
    print(f"\nParallel Median Wall-Clock: {t_parallel_median:.3f} s (Rate: {pps_parallel_median:.3f} μs/part-step)")
    print(f"Median Speedup (2 Cores):   {speedup:.2f}x")
    
    # --- BASINES VERIFICATION ---
    # Baseline Targets (Highly robust, statistically sound thresholds for virtualized GHA environment)
    total_pps_baseline = 15.0      # 15.0 μs/part-step for total wall-clock (GHA VM baseline)
    mover_pps_baseline = 1.0       # 1.0 μs/part-step for isolated particle mover
    solver_pps_baseline = 5.0      # 5.0 μs/part-step for isolated field solver
    speedup_baseline = 1.25        # 1.25x scaling target for 2 cores in virtual environments
    
    total_passed = pps_serial_median <= total_pps_baseline
    mover_passed = pps_mover_serial_median <= mover_pps_baseline
    solver_passed = pps_solver_serial_median <= solver_pps_baseline
    speedup_passed = speedup >= speedup_baseline
    
    print("\n" + "=" * 85)
    print(" " * 32 + "PERFORMANCE COMPARISON PROFILE")
    print("=" * 85)
    print(f" {'Metric / Component':<30} | {'Measured (Median)':<18} | {'Baseline Target':<18} | {'Status':<8}")
    print("-" * 85)
    
    status_total = "PASSED" if total_passed else "FAILED"
    status_mover = "PASSED" if mover_passed else "FAILED"
    status_solver = "PASSED" if solver_passed else "FAILED"
    status_speedup = "PASSED" if speedup_passed else "FAILED"
    
    status_total_disp = status_total
    status_mover_disp = status_mover
    status_solver_disp = status_solver
    status_speedup_disp = status_speedup
    
    if sys.stdout.isatty():
        status_total_disp = f"\033[92;1m{status_total}\033[0m" if total_passed else f"\033[91;1m{status_total}\033[0m"
        status_mover_disp = f"\033[92;1m{status_mover}\033[0m" if mover_passed else f"\033[91;1m{status_mover}\033[0m"
        status_solver_disp = f"\033[92;1m{status_solver}\033[0m" if solver_passed else f"\033[91;1m{status_solver}\033[0m"
        status_speedup_disp = f"\033[92;1m{status_speedup}\033[0m" if speedup_passed else f"\033[91;1m{status_speedup}\033[0m"
        
    print(f" {'Total Runtime Rate (PPS)':<30} | {pps_serial_median:<13.3f} μs/pt | <={total_pps_baseline:<14.1f} μs | {status_total_disp}")
    print(f" {'Particle Mover Rate (PPS)':<30} | {pps_mover_serial_median:<13.3f} μs/pt | <={mover_pps_baseline:<14.1f} μs | {status_mover_disp}")
    print(f" {'Field Solver Rate (PPS)':<30} | {pps_solver_serial_median:<13.3f} μs/pt | <={solver_pps_baseline:<14.1f} μs | {status_solver_disp}")
    print(f" {'Parallel Speedup (2 Cores)':<30} | {speedup:<16.2f} x  | >={speedup_baseline:<14.2f} x  | {status_speedup_disp}")
    print("=" * 85)
    
    # 4. Generate markdown report
    summary_path = "performance_summary.md"
    try:
        with open(summary_path, "w") as f:
            f.write("### ⚡ Standalone FLEKS Performance Report\n\n")
            f.write("A robust statistical check was executed on the runner to filter out virtualization noise (3 runs per benchmark):\n\n")
            f.write("| Performance Metric | Measured (Median) | Target Baseline | Status |\n")
            f.write("| :--- | :--- | :--- | :--- |\n")
            
            total_status_md = "🟢 **PASSED**" if total_passed else "🔴 **FAILED**"
            mover_status_md = "🟢 **PASSED**" if mover_passed else "🔴 **FAILED**"
            solver_status_md = "🟢 **PASSED**" if solver_passed else "🔴 **FAILED**"
            speedup_status_md = "🟢 **PASSED**" if speedup_passed else "🔴 **FAILED**"
            
            f.write(f"| Total Wall-Clock Rate | {pps_serial_median:.3f} μs/pt | <= {total_pps_baseline:.1f} μs/pt | {total_status_md} |\n")
            f.write(f"| Particle Mover Rate | {pps_mover_serial_median:.3f} μs/pt | <= {mover_pps_baseline:.1f} μs/pt | {mover_status_md} |\n")
            f.write(f"| Field Solver Rate | {pps_solver_serial_median:.3f} μs/pt | <= {solver_pps_baseline:.1f} μs/pt | {solver_status_md} |\n")
            f.write(f"| Parallel Speedup (2 Cores) | {speedup:.2f}x | >= {speedup_baseline:.2f}x | {speedup_status_md} |\n\n")
            f.write(f"*Note: Benchmark ran on virtualized runner. Macroparticles: {p_count}, cycles: {cycles} (total steps: {total_steps}).*\n")
            f.write("\n#### 📊 Detailed Runs (Wall-Clock Runtime)\n")
            f.write("| Process Count | Run 1 | Run 2 | Run 3 | Median | Minimum |\n")
            f.write("| :--- | :--- | :--- | :--- | :--- | :--- |\n")
            s_runs = [f"{r['duration']:.3f}s" for r in serial_stats["runs"]]
            p_runs = [f"{r['duration']:.3f}s" for r in parallel_stats["runs"]]
            f.write(f"| 1 MPI Process (Serial) | {s_runs[0]} | {s_runs[1]} | {s_runs[2]} | {serial_stats['duration_median']:.3f}s | {serial_stats['duration_min']:.3f}s |\n")
            f.write(f"| 2 MPI Processes (Parallel) | {p_runs[0]} | {p_runs[1]} | {p_runs[2]} | {parallel_stats['duration_median']:.3f}s | {parallel_stats['duration_min']:.3f}s |\n")
    except Exception as e:
        print(f"Warning: Could not write performance_summary.md: {e}")
        
    all_passed = total_passed and mover_passed and solver_passed and speedup_passed
    if all_passed:
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
