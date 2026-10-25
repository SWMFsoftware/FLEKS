#!/usr/bin/env python3
"""Repeat alternating GPU-demo cases and retain timings, diagnostics and provenance."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shlex
import statistics
import subprocess
import time

from validate_demo import ROOT, compare_records, parse_diagnostics, require, validate_diagnostics


# Cases recorded before CGS2 became the default pin gmres_batched=0 so their
# meaning is unchanged; "combined" is the current recommended configuration.
CASES = {
    "baseline": {"gmres_batched": 0},
    "pairs": {"pair_symmetry": 1, "gmres_batched": 0},
    "pairs-fused": {"pair_symmetry": 1, "fused_push": 1, "gmres_batched": 0},
    "folded": {"pair_symmetry": 1, "fold_charge": 1, "gmres_batched": 0},
    "chunks": {"pair_symmetry": 1, "cell_chunk_size": 256, "gmres_batched": 0},
    "batched": {"gmres_batched": 1},
    "fused": {"fused_push": 1, "gmres_batched": 0},
    "combined": {"pair_symmetry": 1, "gmres_batched": 1, "fused_push": 1},
}


def capture(command):
    try:
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                                timeout=15, check=False)
        return {"returncode": result.returncode,
                "output": result.stdout.strip(), "stderr": result.stderr.strip()}
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"error": str(error)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", default=str(ROOT / "bin/main_gpu.exe"))
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=["baseline", "pairs"])
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--steps", type=int, default=16)
    parser.add_argument("--warmup-steps", type=int, default=3)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--launcher", default="mpirun")
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--timing", choices=("stages", "throughput"), default="throughput")
    parser.add_argument("--set", nargs="*", default=[], metavar="KEY=VALUE",
                        help="Common physics/grid overrides; case and timing controls stay fixed")
    args = parser.parse_args()
    require(args.repeats > 0 and args.ranks > 0 and args.timeout > 0,
            "Require positive repeats, ranks and timeout")
    require(0 <= args.warmup_steps < args.steps, "Require 0 <= warmup_steps < steps")
    executable = Path(args.executable).resolve()
    require(executable.is_file(), f"Executable not found: {executable}")
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    require(not any(output.iterdir()), "Use an empty output directory to preserve earlier results")
    parameters = dict(n_cell=64, max_grid_size=32, particles_per_cell=100,
                      initial_condition="langmuir", perturbation=0.01, b0=0.5,
                      dt=0.1, correct_div_e=1, solver_tolerance=1e-6,
                      response="cell", deposit="cell", seed=1)
    reserved = {"mode", "n_step", "warmup_steps", "diagnostic_interval",
                "benchmark_timing", "microbenchmarks"}
    reserved.update(key for case in CASES.values() for key in case)
    for setting in args.set:
        key, separator, value = setting.partition("=")
        require(separator and key and value and key not in reserved,
                f"Invalid or reserved common override: {setting}")
        parameters[key] = value
    parameters.update(mode="benchmark", n_step=args.steps, warmup_steps=args.warmup_steps,
                      diagnostic_interval=args.steps, benchmark_timing=args.timing,
                      microbenchmarks=0)
    cases = list(dict.fromkeys(["baseline", *args.cases]))
    results = dict(
        arguments=vars(args), common_parameters=parameters, trials=[], summary={},
        provenance={
            "host": platform.node(), "platform": platform.platform(),
            "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
            "revision": capture(["git", "rev-parse", "HEAD"]),
            "amrex_revision": capture(["git", "-C", "util/AMREX", "rev-parse", "HEAD"]),
            "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
            "working_tree": capture(["git", "status", "--short"]),
            "source_diff": capture(["git", "diff", "--", "src_gpu"]),
            "makefile_conf": (ROOT / "Makefile.conf").read_text(),
            "compiler": capture(["nvcc", "--version"]),
            "gpu": capture(["nvidia-smi", "--query-gpu=name,uuid,driver_version,memory.total,"
                            "clocks.sm,clocks.mem", "--format=csv"]),
            "affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
            "mpi_environment": {key: value for key, value in os.environ.items()
                                if key.startswith(("OMPI_MCA_", "UCX_")) or key == "OPAL_PREFIX"},
        })
    # Include untracked source content in provenance through an aggregate hash.
    digest = hashlib.sha256()
    for source in sorted((ROOT / "src_gpu").glob("*")):
        if source.suffix in (".cpp", ".h", ".py"):
            digest.update(source.name.encode())
            digest.update(source.read_bytes())
    results["provenance"]["source_sha256"] = digest.hexdigest()
    destination = output / "results.json"
    reference = None
    reference_warnings = None
    try:
        for repeat in range(args.repeats):
            # Reverse order every other trial to reduce monotonic drift bias.
            order = cases if repeat % 2 == 0 else list(reversed(cases))
            for case in order:
                settings = {**parameters, **CASES[case]}
                command = [str(executable), *(f"{key}={value}" for key, value in settings.items())]
                if args.launcher:
                    command = [*shlex.split(args.launcher), "-np", str(args.ranks), *command]
                else:
                    require(args.ranks == 1, "Multiple ranks require a launcher")
                label = f"{repeat + 1:02d}-{case}"
                trial = dict(case=case, repeat=repeat + 1, command=command,
                             load_before=os.getloadavg(), passed=False)
                results["trials"].append(trial)
                start = time.monotonic()
                with (output / f"{label}.log").open("w") as log:
                    process = subprocess.run(command, cwd=ROOT, stdout=log,
                                             stderr=subprocess.STDOUT, timeout=args.timeout)
                trial["process_seconds"] = time.monotonic() - start
                trial["returncode"] = process.returncode
                trial["load_after"] = os.getloadavg()
                text = (output / f"{label}.log").read_text()
                require(process.returncode == 0, f"{label}: failed; see log")
                require("FAIL" not in text, f"{label}: correctness failure; see log")
                records = parse_diagnostics(text)
                # The executable checks the configured relative true-residual
                # criterion and exits on failure. The correctness suite's
                # absolute 1e-7 bound is specific to its smaller test systems.
                validate_diagnostics(records, args.steps, max_residual=math.inf)
                if reference is None:
                    reference = records
                compare_records(reference, records, label)
                trial["diagnostics"] = records
                trial["arena_peak_mb"] = {
                    kind: {"min_rank": float(low), "max_rank": float(high)}
                    for kind, low, high in re.findall(
                        r"\[The\s+Arena\] max space \(MB\) (allocated|used)\s+"
                        r"spread across MPI: \[([\d.]+) \.\.\. ([\d.]+)\]", text)
                }
                timings = {}
                for line in text.splitlines():
                    if line.startswith("TIMING "):
                        fields = dict(item.split("=", 1) for item in line.split()[1:])
                        timings[fields["stage"]] = float(fields["seconds"])
                    if line.startswith("BENCHMARK "):
                        trial["benchmark"] = dict(item.split("=", 1) for item in line.split()[1:])
                require(trial.get("benchmark", {}).get("aborted") == "0", f"{label}: aborted window")
                warnings = int(trial["benchmark"]["correction_warnings"])
                if reference_warnings is None:
                    reference_warnings = warnings
                require(warnings <= reference_warnings,
                        f"{label}: extra correction-limit warnings ({warnings} vs {reference_warnings})")
                require(int(trial["benchmark"]["measured_steps"]) == args.steps - args.warmup_steps,
                        f"{label}: wrong number of measured steps")
                require(math.isfinite(timings.get("timestep", math.nan)) and timings["timestep"] > 0,
                        f"{label}: missing/invalid timestep timing")
                trial.update(timings=timings, passed=True)
                print(f"{label}: {timings['timestep']:.6f} s/step", flush=True)
                destination.write_text(json.dumps(results, indent=2) + "\n")
        for case in cases:
            values = [trial["timings"]["timestep"] for trial in results["trials"] if trial["case"] == case]
            results["summary"][case] = dict(median=statistics.median(values), min=min(values),
                                           max=max(values), stdev=statistics.stdev(values) if len(values) > 1 else 0)
        baseline = results["summary"]["baseline"]["median"]
        for summary in results["summary"].values():
            summary["speedup_vs_baseline"] = baseline / summary["median"]
        print(json.dumps(results["summary"], indent=2))
    finally:
        destination.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
