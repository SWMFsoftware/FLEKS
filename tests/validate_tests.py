#!/usr/bin/env python3
"""Common runner for the standalone FLEKS test suite.

This module holds only the shared infrastructure:

  * preparing / cleaning the run directory,
  * launching FLEKS.exe and PostProc.pl,
  * reading the PIC energy log and the test-particle tracer log,
  * the generic test-particle validator,
  * build-configuration pre-flight checks,
  * the test-discovery + execution loop and the summary table.

Each per-test validation (energy-log check, plot-output check, particle-log
tolerance) lives in a small ``validate.py`` module inside that test's own
directory (e.g. ``tests/beam/validate.py``).  This runner loads only the module
for the test currently being run, so it does not have to import the code for
every test up front.  Shared helpers (e.g. the hybrid family) live in
``tests/_shared/``.

Output is controlled with the standard :mod:`logging` levels:
  * INFO    -- essential per-test results and the summary (default),
  * DEBUG   -- intermediate diagnostics (enable with --verbose / -v),
  * WARNING / ERROR -- warnings and failures.
"""
import importlib
import logging
import os
import shutil
import subprocess
import sys

# Directory used for simulation output. Defaults to "run_test" but can be
# overridden with --run-dir so a second test can run without clobbering a
# currently-running job that owns the default run_test/.
RUN_DIR = "run_test"

logger = logging.getLogger("fleks.validate")


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
def setup_logging(verbose=False):
    """Configure the root logger used by the runner and per-test modules.

    Without *verbose* only INFO+ messages are shown (essential results, the
    summary, warnings and errors).  With *verbose* DEBUG messages (per-check
    diagnostics, energy numbers) are shown too.

    The root logger is configured (rather than a named logger) so that the
    per-test modules -- which each use ``logging.getLogger(__name__)`` -- print
    through the same handler and level without needing explicit configuration.
    """
    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger()
    root.setLevel(level)
    if not root.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter("%(message)s"))
        root.addHandler(handler)
    else:
        for h in root.handlers:
            h.setLevel(level)
    # Make sure child loggers propagate up to the configured root handler.
    for lname in list(logging.Logger.manager.loggerDict):
        logging.getLogger(lname).propagate = True


# ---------------------------------------------------------------------------
# Run-directory management
# ---------------------------------------------------------------------------
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
                    logger.warning("  Could not remove %s: %s", entry_path, e)


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------
def run_test(test_dir, nprocs=1, param_text=None):
    param_file = os.path.join(test_dir, "PARAM.in")
    logger.debug("Running test in %s...", test_dir)
    prepare_run_dir()

    # Verify that PostIDL.exe exists; PostProc.pl needs it to produce .out files.
    postidl_link = os.path.join(RUN_DIR, "PC", "PostIDL.exe")
    if os.path.islink(postidl_link) and not os.path.exists(postidl_link):
        real = os.path.realpath(postidl_link)
        logger.warning("  [WARN] Broken symlink: %s -> %s", postidl_link, real)
        logger.warning("  [WARN] PostIDL.exe is missing. Build it with 'make PIDL' "
                       "before running tests that check plot output (.out files).")

    # Write run_test/PARAM.in: use the supplied text (e.g. a patched solver
    # variant) if given, otherwise copy the test directory's PARAM.in.
    if param_text is not None:
        with open(RUN_DIR + "/PARAM.in", "w") as f:
            f.write(param_text)
    else:
        shutil.copy(param_file, RUN_DIR + "/PARAM.in")

    # Build the command: serial for nprocs==1, mpirun otherwise
    if nprocs <= 1:
        cmd = ["./FLEKS.exe"]
        logger.debug("  Running in serial mode (no MPI)...")
    else:
        cmd = ["mpirun", "-n", str(nprocs), "./FLEKS.exe"]
        logger.debug("  Running with %d MPI processes...", nprocs)

    # Run the command inside run_test/
    result = subprocess.run(cmd, cwd=RUN_DIR, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        logger.error("Error running FLEKS.exe for %s:", test_dir)
        logger.error("--- FLEKS stdout ---")
        logger.error(result.stdout)
        logger.error("--- FLEKS stderr ---")
        logger.error(result.stderr)
        return None, result.returncode

    # Automatically run post-processing on the generated plots
    pp = subprocess.run(["./PostProc.pl", "-v"], cwd=RUN_DIR,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if pp.returncode != 0:
        logger.warning("  [WARN] PostProc.pl exited with code %s:", pp.returncode)
        if pp.stdout:
            logger.warning(pp.stdout)
        if pp.stderr:
            logger.warning(pp.stderr)

    return result.stdout, 0


def run_and_validate(test_dir, display_name, validator, nprocs, results,
                     param_text=None, base_name=None):
    """Run one FLEKS test (optionally with a patched PARAM.in) and record the
    outcome in *results* as (name, status, reason).  *validator* is the object
    returned by load_validator() -- it exposes ``validate_log``, ``plot`` and
    ``particle_tol`` -- or None for a generic completion check.

    Mirrors the former main-loop body so a single test can be run with several
    PARAM variants (used by the free-stream test)."""
    if base_name is None:
        base_name = display_name
    logger.info("")
    logger.info("=" * 50)
    logger.info("Starting test: %s", display_name)
    logger.info("=" * 50)
    try:
        stdout, code = run_test(test_dir, nprocs=nprocs, param_text=param_text)
        if code != 0 or stdout is None:
            logger.error("FAIL: %s execution failed with exit code %s",
                         display_name, code)
            results.append((display_name, "FAILED", f"Execution failed (code {code})"))
            return

        # Read the PIC energy log (the only diagnostic log produced by FLEKS).
        pic_diags = read_pic_log(RUN_DIR)

        val_res = False
        reason = "Validation skipped"

        if validator is not None and validator.validate_log is not None:
            val_res, reason = validator.validate_log(pic_diags=pic_diags,
                                                     test_name=base_name)
            if not val_res:
                logger.error("%s: FAILED (%s)", display_name, reason)
                results.append((display_name, "FAILED", reason))
                return
        else:
            logger.debug("Validating %s (generic check)...", display_name)
            logger.info("%s (generic check): PASSED", display_name)
            val_res = True
            reason = "Passed"

        # Read the test-particle tracer log (log_pt_n*.log) and validate it for
        # tests that enable #PARTICLETRACKER T (declared via PARTICLE_TOL).
        if validator is not None and validator.particle_tol:
            pt_diags = read_pt_log(RUN_DIR)
            pt_res, pt_reason = validate_test_particles(
                pt_diags, test_name=base_name, tol=validator.particle_tol)
            if not pt_res:
                logger.error("%s: test-particle check failed (%s)",
                             display_name, pt_reason)
                results.append((display_name, "FAILED",
                                f"test-particle check failed: {pt_reason}"))
                return

        # Validate output plotfiles.
        plot_validator = validator.plot if validator is not None else None
        if plot_validator is None:
            plot_res, plot_reason = True, "Passed (no plot-file check)"
        else:
            plot_res, plot_reason = plot_validator(base_name)
        if not plot_res:
            logger.error("%s: plot check failed (%s)", display_name, plot_reason)
            results.append((display_name, "FAILED", f"plot check failed: {plot_reason}"))
        else:
            logger.info("%s: PASSED", display_name)
            results.append((display_name, "PASSED", "Passed"))

    finally:
        # Always clean up run output after each test to keep disk usage low.
        logger.debug("  Cleaning up run_test/ output for %s...", display_name)
        cleanup_run_dir()


# ---------------------------------------------------------------------------
# Build-configuration pre-flight checks
#
# The standalone suite expects a single binary built as:
#     ./Config.pl -lev=2 -u=Exo && make
# These helpers inspect the *configured* source (include/UserSource.h) and AMR
# level (include/Constants.h) so that a misconfigured build is caught before a
# test runs, instead of producing confusing physics failures (empty source) or
# an abort at start-up (nLevMax = 1 on lightwave).
# ---------------------------------------------------------------------------

# Tests that exercise the Exosphere ionization source.  They require the Exo
# user source (-u=Exo); with the default empty source they compile and run but
# no ionization source term is applied, so the physics checks fail.
EXO_SOURCE_TESTS = {
    "chargeexchange", "photoionization", "electronimpact",
    "recombination", "chemistry",
}

# Tests that require AMR (nLevMax >= 2).
AMR_TESTS = {"lightwave"}


def configured_user_source_is_exo():
    """Return True if include/UserSource.h is the Exo source.

    Reads the configured source file (the one copied by Config.pl -u=Exo from
    userfiles/ExoSource.h).  Returns False for the default empty source and on
    any I/O error (treat as "not Exo" so the check fails loudly).
    """
    try:
        with open("include/UserSource.h", "r") as f:
            content = f.read()
    except (OSError, IOError):
        return False
    # The Exo source guard is _EXOSOURCE_H_ and it sets useFluidSource = true.
    # The default empty source uses _USERSOURCE_H_ and leaves useFluidSource
    # false.  Checking either marker is enough to distinguish them.
    if "_EXOSOURCE_H_" in content or "useFluidSource = true" in content:
        return True
    return False


def configured_nlevmax():
    """Return the configured nLevMax from include/Constants.h (int).

    Returns None if it cannot be determined (file missing / malformed).
    """
    try:
        with open("include/Constants.h", "r") as f:
            content = f.read()
    except (OSError, IOError):
        return None
    import re
    m = re.search(r"nLevMax\s*=\s*(\d+)", content)
    return int(m.group(1)) if m else None


def preflight_check(test_name):
    """Fail fast if the configured binary cannot run `test_name`.

    Returns (ok: bool, reason: str|None).  When ok is False the caller should
    skip the test and report the reason (an actionable build instruction).
    """
    if test_name in EXO_SOURCE_TESTS and not configured_user_source_is_exo():
        return False, (
            f"Test '{test_name}' needs the Exosphere user source, but "
            "include/UserSource.h is the default empty source.\n"
            "  Rebuild with:  ./Config.pl -lev=2 -u=Exo && make"
        )

    if test_name in AMR_TESTS:
        nlev = configured_nlevmax()
        if nlev is None:
            return False, (
                f"Test '{test_name}' needs AMR (nLevMax >= 2), but "
                "include/Constants.h could not be read.\n"
                "  Rebuild with:  ./Config.pl -lev=2 -u=Exo && make"
            )
        if nlev < 2:
            return False, (
                f"Test '{test_name}' needs AMR (nLevMax >= 2), but the "
                f"binary is configured with nLevMax = {nlev}.\n"
                "  Rebuild with:  ./Config.pl -lev=2 -u=Exo && make"
            )

    return True, None


# ---------------------------------------------------------------------------
# Log readers
# ---------------------------------------------------------------------------
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
    import math
    tol = tol or {}
    min_rows = int(tol.get("min_rows", 1))
    launch_threshold = float(tol.get("launch_threshold", 0.5))
    max_speed = float(tol.get("max_speed", 10.0))
    expected_active = set(int(s) for s in tol.get("expected_active_species", []))

    logger.debug("Validating Test-Particle Tracer Output...")

    if not pt_diags:
        logger.debug("  [PT] No test-particle log (log_pt_n*.log) found.")
        return False, "No test-particle log file"

    n_rows = len(pt_diags)
    n_species = sum(1 for k in pt_diags[0] if k.startswith("mass"))
    logger.debug("  [PT] %d log rows, %d species (from log_pt_n*.log).",
                 n_rows, n_species)

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
    logger.debug("  [PT] active (seeded) species: %s", sorted(active))

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
        logger.info("Test-Particle Tracer Output: PASSED")
        return True, "Passed"
    return False, "; ".join(reasons)


# ---------------------------------------------------------------------------
# Per-test validator loading
# ---------------------------------------------------------------------------
def load_validator(test_name):
    """Dynamically import the ``validate.py`` module for *test_name*.

    Returns a small object with attributes ``validate_log``, ``plot`` and
    ``particle_tol`` (some may be None) so the runner can invoke it uniformly.
    If the test has no ``validate.py`` (e.g. ``hyper_resistivity``), returns
    None and the runner falls back to a generic no-op check.
    """
    module_name = f"tests.{test_name}.validate"
    try:
        module = importlib.import_module(module_name)
    except ImportError as exc:
        # Only swallow failures caused by a genuinely missing per-test module
        # (hyper_resistivity and friends).  Any other import error (a bug in
        # the module) should surface loudly.
        if exc.name != module_name and not str(exc).startswith("No module named"):
            raise
        logger.debug("No per-test module %s; using generic check.", module_name)
        return None

    class _Validator:
        pass

    v = _Validator()
    v.validate_log = getattr(module, "validate_log", None)
    # plot() expects a single positional arg (test_name).
    v.plot = getattr(module, "validate_plot", None)
    v.particle_tol = getattr(module, "PARTICLE_TOL", None)
    return v


def _sync_shared_run_dir():
    """Push the current RUN_DIR into the shared helpers used by plot checks."""
    import tests._shared.hybrid as _hyb
    _hyb.set_run_dir(RUN_DIR)


# ---------------------------------------------------------------------------
# Test discovery + main
# ---------------------------------------------------------------------------
def discover_tests(tests_dir="tests"):
    """Return a sorted list of (test_dir, name) for directories with PARAM.in.

    Excludes "performance" (a benchmark, not a pass/fail test) and "run_test"
    (the shared run directory created by prepare_run_dir, which contains a
    PARAM.in and would otherwise be discovered and re-run as a test).
    """
    subdirs = sorted([d for d in os.listdir(tests_dir)
                      if os.path.isdir(os.path.join(tests_dir, d))
                      and d not in ["performance", "run_test"]])
    tests = []
    for d in subdirs:
        test_dir = os.path.join(tests_dir, d)
        param_file = os.path.join(test_dir, "PARAM.in")
        if os.path.exists(param_file):
            tests.append((test_dir, d))
    return tests


def run_one_test(test_dir, name, nprocs, results):
    """Pre-flight check + run one test (or its PARAM variants) + record result."""
    # Fail fast if the configured binary cannot run this test (wrong user
    # source or AMR level), before spending time executing it.
    ok, preflight_reason = preflight_check(name)
    if not ok:
        logger.warning("SKIP: %s cannot run with the configured build.", name.upper())
        logger.warning(preflight_reason)
        results.append((name.upper(), "FAILED",
                        f"build misconfigured: {preflight_reason}"))
        return

    validator = load_validator(name)

    # The free-stream test is a single directory holding two parameter files
    # (PARAM.in and PARAM.in.hybrid); the runner exercises both field solvers
    # by running it once per file. All other tests run once as written.
    if name == "freestream":
        with open(os.path.join(test_dir, "PARAM.in")) as _f:
            _fullpic = _f.read()
        run_and_validate(test_dir, "FREESTREAM (FULL PIC)", validator, nprocs,
                         results, param_text=_fullpic, base_name="freestream")
        _hybrid_path = os.path.join(test_dir, "PARAM.in.hybrid")
        with open(_hybrid_path) as _f:
            _hybrid = _f.read()
        run_and_validate(test_dir, "FREESTREAM (HYBRID HALL-OFF)", validator,
                         nprocs, results, param_text=_hybrid,
                         base_name="freestream")
    else:
        run_and_validate(test_dir, name.upper(), validator, nprocs, results,
                         base_name=name)


def main():
    # Parse nprocs: -n N or --nprocs N
    nprocs = 1
    for i, arg in enumerate(sys.argv):
        if arg in ("-n", "--nprocs"):
            try:
                nprocs = int(sys.argv[i + 1])
            except (IndexError, ValueError):
                logger.error("Error: %s requires an integer argument (number of MPI processes).", arg)
                sys.exit(1)
            continue

    if nprocs < 1:
        logger.error("Error: Number of processes must be >= 1.")
        sys.exit(1)

    # Parse --summary-file PATH (custom output path for CI serial/parallel split)
    summary_file = "tests/summary.md"
    for i, arg in enumerate(sys.argv):
        if arg == "--summary-file":
            try:
                summary_file = sys.argv[i + 1]
            except IndexError:
                logger.error("Error: --summary-file requires a path argument.")
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
                logger.error("Error: --test requires a test name argument.")
                sys.exit(1)
            continue
        if arg == "--run-dir":
            global RUN_DIR
            try:
                RUN_DIR = sys.argv[i + 1]
            except IndexError:
                logger.error("Error: --run-dir requires a path argument.")
                sys.exit(1)
            continue

    # Parse --verbose / -v: show DEBUG diagnostics from each validator.
    verbose = "-v" in sys.argv or "--verbose" in sys.argv
    setup_logging(verbose=verbose)

    # Work from the FLEKS root so `include/...`, `tests/...` and `run_test/`
    # resolve relative to the repository, and so the `tests` package is
    # importable.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(script_dir)
    os.chdir(repo_root)
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    logger.debug("Working directory set to: %s", os.getcwd())

    _sync_shared_run_dir()

    tests = discover_tests()
    if not tests:
        logger.error("No tests found in tests/ subdirectories!")
        sys.exit(1)

    # Filter to a single test if --test was given.
    if selected_test is not None:
        matching = [t for t in tests if t[1] == selected_test]
        if not matching:
            available = ", ".join(t[1] for t in tests)
            logger.error("Error: test '%s' not found.", selected_test)
            logger.error("Available tests: %s", available)
            sys.exit(1)
        tests = matching
        logger.debug("Selected test: %s", selected_test)

    results = []  # Collect results for summary table

    for test_dir, name in tests:
        run_one_test(test_dir, name, nprocs, results)

    # ----------------------------------------------------
    # Print Summary Table
    # ----------------------------------------------------
    logger.info("")
    logger.info("=" * 80)
    logger.info(" " * 32 + "TEST SUMMARY")
    logger.info("=" * 80)
    logger.info(" %-28s | %-8s | %s", "Test Name", "Status", "Failure Reason / Details")
    logger.info("-" * 80)

    all_passed = True
    for name_str, status, reason in results:
        if status == "PASSED":
            status_display = "\033[92;1mPASSED\033[0m" if sys.stdout.isatty() else "PASSED"
        else:
            status_display = "\033[91;1mFAILED\033[0m" if sys.stdout.isatty() else "FAILED"
            all_passed = False
        reason_display = reason if status != "PASSED" else ""
        logger.info(" %-28s | %-8s | %s", name_str, status_display, reason_display)

    logger.info("=" * 80)

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
        logger.warning("Warning: Could not write summary.md: %s", e)

    if all_passed:
        if sys.stdout.isatty():
            logger.info("\033[92;1m\nALL STANDALONE FLEKS TESTS PASSED SUCCESSFULLY!\033[0m\n")
        else:
            logger.info("\nALL STANDALONE FLEKS TESTS PASSED SUCCESSFULLY!\n")
        sys.exit(0)
    else:
        if sys.stdout.isatty():
            logger.error("\033[91;1m\nSOME TESTS FAILED.\033[0m\n")
        else:
            logger.error("\nSOME TESTS FAILED.\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
