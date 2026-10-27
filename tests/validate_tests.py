#!/usr/bin/env python3
"""Common runner for the standalone FLEKS test suite.

This module holds only the shared infrastructure:

  * preparing / cleaning the run directory,
  * launching FLEKS.exe and PostProc.pl,
  * reading the PIC energy log and the test-particle tracer log,
  * the generic test-particle validator,
  * build-configuration pre-flight checks,
  * the test-discovery + execution loop and the summary table.

Each per-test validation lives in a small ``validate.py`` module
inside that test's own directory (e.g. ``tests/beam/validate.py``).
This runner loads only the module for the test currently being run,
so it does not have to import the code for every test up front.
Shared helpers (e.g. the hybrid family) live in ``tests/_shared/``.

Output is controlled with the standard :mod:`logging` levels:
  * INFO    -- essential per-test results and the summary (default),
  * DEBUG   -- intermediate diagnostics (enable with --verbose / -v),
  * WARNING / ERROR -- warnings and failures.
"""
import glob
import importlib
import logging
import os
import shutil
import subprocess
import sys

# Directory used for simulation output. Defaults can be overridden with
# --run-dir so a second test can run without clobbering a currently-running
# job that owns the default run_test/.
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


def ensure_postidl():
    """Ensure PostIDL.exe is compiled before running any test.

    If it is not present at any candidate location we compile it with
    ``make PIDL`` and re-search.  Returns True if PostIDL.exe is found.
    """
    postidl_candidates = [
        "bin/PostIDL.exe",          # standalone build
        "../../bin/PostIDL.exe",    # SWMF integrated build
    ]

    def _found():
        return any(os.path.isfile(c) for c in postidl_candidates)

    if _found():
        return True

    # Not built yet -- compile it. This keeps all tests (pcai, beam, ...) able
    # to produce plot output regardless of how the tree was set up.
    logger.warning("  [WARN] PostIDL.exe not found. Building it with 'make PIDL'...")
    try:
        build = subprocess.run(["make", "PIDL"], stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, timeout=900)
        build.stdout = (build.stdout or b"").decode("utf-8", errors="replace")
        build.stderr = (build.stderr or b"").decode("utf-8", errors="replace")
    except subprocess.TimeoutExpired:
        logger.warning("  [WARN] 'make PIDL' timed out; PostIDL.exe may be missing.")
        return False

    if build.returncode != 0 or not _found():
        logger.warning("  [WARN] 'make PIDL' failed; PostIDL.exe is missing. "
                       "Tests that check plot output (.out files) will not "
                       "produce frames.")
        if build.stderr:
            logger.warning(build.stderr)
        return False

    logger.info("  [INFO] PostIDL.exe built successfully.")
    return True


def prepare_run_dir(run_dir=None):
    if run_dir is None:
        run_dir = RUN_DIR
    os.makedirs(run_dir, exist_ok=True)

    ensure_postidl()

    # PostProc.pl / pIDL live in share/Scripts/; the standalone layout keeps
    # "share" next to FLEKS root, the SWMF layout two levels up.
    if os.path.isfile("share/Scripts/PostProc.pl"):
        postproc_src = os.path.abspath("share/Scripts/PostProc.pl")
        pidl_src = os.path.abspath("share/Scripts/pIDL")
    elif os.path.isfile("../../share/Scripts/PostProc.pl"):
        postproc_src = os.path.abspath("../../share/Scripts/PostProc.pl")
        pidl_src = os.path.abspath("../../share/Scripts/pIDL")
    else:
        postproc_src = os.path.abspath("share/Scripts/PostProc.pl")
        pidl_src = os.path.abspath("share/Scripts/pIDL")

    postidl_candidates = [
        os.path.abspath("bin/PostIDL.exe"),
        os.path.abspath("../../bin/PostIDL.exe"),
    ]
    postidl_src = None
    for candidate in postidl_candidates:
        if os.path.isfile(candidate):
            postidl_src = candidate
            break
    if postidl_src is None:
        postidl_src = os.path.abspath("../../bin/PostIDL.exe")

    fleks_src = os.path.abspath("bin/FLEKS.exe")
    safe_symlink(os.path.relpath(fleks_src, run_dir), os.path.join(run_dir, "FLEKS.exe"))
    safe_symlink(os.path.relpath(postproc_src, run_dir), os.path.join(run_dir, "PostProc.pl"))

    pc_dir = os.path.join(run_dir, "PC")
    os.makedirs(pc_dir, exist_ok=True)
    os.makedirs(os.path.join(pc_dir, "plots"), exist_ok=True)
    os.makedirs(os.path.join(pc_dir, "restartOUT"), exist_ok=True)

    safe_symlink(os.path.relpath(pidl_src, pc_dir), os.path.join(pc_dir, "pIDL"))
    safe_symlink(os.path.relpath(postidl_src, pc_dir), os.path.join(pc_dir, "PostIDL.exe"))


def cleanup_run_dir(run_dir=None):
    """Remove simulation output files from run directory after each test.

    Deletes the contents of PC/plots/ and PC/restartOUT/ while keeping the
    directory structure and symlinks.
    """
    if run_dir is None:
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
def run_test(test_dir, nprocs=1, param_text=None, run_dir=None):
    if run_dir is None:
        run_dir = RUN_DIR
    param_file = os.path.join(test_dir, "PARAM.in")
    logger.debug("Running test in %s...", test_dir)
    prepare_run_dir(run_dir)
    cleanup_run_dir(run_dir)

    postidl_link = os.path.join(run_dir, "PC", "PostIDL.exe")
    if os.path.islink(postidl_link) and not os.path.exists(postidl_link):
        real = os.path.realpath(postidl_link)
        logger.warning("  [WARN] Broken symlink: %s -> %s", postidl_link, real)
        logger.warning("  [WARN] PostIDL.exe is missing. Build it with 'make PIDL' "
                       "before running tests that check plot output (.out files).")

    if param_text is not None:
        with open(os.path.join(run_dir, "PARAM.in"), "w") as f:
            f.write(param_text)
    else:
        if os.path.isfile(param_file):
            shutil.copy(param_file, os.path.join(run_dir, "PARAM.in"))
        else:
            candidates = sorted(glob.glob(os.path.join(test_dir, "PARAM.in*")))
            if candidates:
                shutil.copy(candidates[0], os.path.join(run_dir, "PARAM.in"))

    # Copy the test's own data files (coefficient tables, lookup tables, ...)
    # into the run directory: the decks reference them by a relative name.
    for entry in sorted(os.listdir(test_dir)):
        if entry.startswith(".") or entry.startswith("PARAM.in") or \
                entry in ("validate.py", "README.md"):
            continue
        src = os.path.join(test_dir, entry)
        if os.path.isfile(src):
            shutil.copy(src, os.path.join(run_dir, entry))

    if nprocs <= 1:
        cmd = ["./FLEKS.exe"]
        logger.debug("  Running in serial mode (no MPI)...")
    else:
        cmd = ["mpirun", "-n", str(nprocs), "./FLEKS.exe"]
        logger.debug("  Running with %d MPI processes...", nprocs)

    result = subprocess.run(cmd, cwd=run_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result.stdout = (result.stdout or b"").decode("utf-8", errors="replace")
    result.stderr = (result.stderr or b"").decode("utf-8", errors="replace")
    # Keep FLEKS's stdout on disk. Several per-cycle diagnostics (notably the
    # AMR div(B) report) are only printed to stdout, so validators that assert
    # on them need a file to read.
    try:
        with open(os.path.join(run_dir, "fleks_stdout.log"), "w",
                  encoding="utf-8", errors="replace") as _f:
            _f.write(result.stdout)
    except OSError:
        pass
    if result.returncode != 0:
        logger.error("Error running FLEKS.exe for %s:", test_dir)
        logger.error("--- FLEKS stdout ---")
        logger.error(result.stdout)
        logger.error("--- FLEKS stderr ---")
        logger.error(result.stderr)
        return None, result.returncode

    pp = subprocess.run(["./PostProc.pl", "-v"], cwd=run_dir,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pp.stdout = (pp.stdout or b"").decode("utf-8", errors="replace")
    pp.stderr = (pp.stderr or b"").decode("utf-8", errors="replace")
    if pp.returncode != 0:
        logger.warning("  [WARN] PostProc.pl exited with code %s:", pp.returncode)
        if pp.stdout:
            logger.warning(pp.stdout)
        if pp.stderr:
            logger.warning(pp.stderr)

    return result.stdout, 0


def run_and_validate(test_dir, display_name, validator, nprocs, results,
                     param_text=None, base_name=None, run_dir=None):
    """Run one FLEKS test (optionally with a patched PARAM.in) and record the
    outcome in *results* as (name, status, reason, elapsed_seconds).  *validator*
    is the object returned by load_validator() -- it exposes ``validate_log``,
    ``plot`` and ``particle_tol`` -- or None for a generic completion check.

    Mirrors the former main-loop body so a single test can be run with several
    PARAM variants (used by the free-stream test)."""
    import time
    if run_dir is None:
        run_dir = RUN_DIR
    if base_name is None:
        base_name = display_name
    logger.info("")
    logger.info("=" * 50)
    logger.info("Starting test: %s", display_name)
    logger.info("=" * 50)
    t0 = time.time()
    try:
        # Plot validators locate their output relative to run_dir; keep the
        # per-test module in sync with it before anything is validated.
        _sync_shared_run_dir(run_dir)
        _sync_validator_run_dir(validator, run_dir)
        stdout, code = run_test(test_dir, nprocs=nprocs, param_text=param_text, run_dir=run_dir)
        dt = time.time() - t0
        if code != 0 or stdout is None:
            logger.error("FAIL: %s execution failed with exit code %s (%.2fs)",
                         display_name, code, dt)
            results.append((display_name, "FAILED", f"Execution failed (code {code})", dt))
            return

        # FLEKS's only diagnostic log is the PIC energy log.
        pic_diags = read_pic_log(run_dir)

        val_res = False
        reason = "Validation skipped"

        if validator is not None and validator.validate_log is not None:
            val_res, reason = validator.validate_log(pic_diags=pic_diags,
                                                     test_name=base_name)
            if not val_res:
                logger.error("%s: FAILED (%s) (%.2fs)", display_name, reason, dt)
                results.append((display_name, "FAILED", reason, dt))
                return
        else:
            logger.debug("Validating %s (generic check)...", display_name)
            logger.debug("%s (generic check): PASSED", display_name)
            val_res = True
            reason = "Passed"

        # Test-particle tracer log (log_pt_n*.log), validated when #PARTICLETRACKER T.
        if validator is not None and validator.particle_tol:
            has_pt = False
            active_param = os.path.join(run_dir, "PARAM.in")
            if os.path.isfile(active_param):
                with open(active_param, "r", errors="ignore") as pf:
                    ptext = pf.read()
                    if "#PARTICLETRACKER" in ptext:
                        pt_block = ptext.split("#PARTICLETRACKER", 1)[1].splitlines()
                        for ln in pt_block[1:]:
                            ln_s = ln.strip()
                            if not ln_s or ln_s.startswith("#"):
                                continue
                            if ln_s.startswith("T"):
                                has_pt = True
                            break
            if has_pt:
                pt_diags = read_pt_log(run_dir)
                pt_res, pt_reason = validate_test_particles(
                    pt_diags, test_name=base_name, tol=validator.particle_tol)
                if not pt_res:
                    logger.error("%s: test-particle check failed (%s) (%.2fs)",
                                 display_name, pt_reason, dt)
                    results.append((display_name, "FAILED",
                                    f"test-particle check failed: {pt_reason}", dt))
                    return

        # Validate output plotfiles.
        plot_validator = validator.plot if validator is not None else None
        if plot_validator is None:
            plot_res, plot_reason = True, "Passed (no plot-file check)"
        else:
            plot_res, plot_reason = plot_validator(base_name)
        dt = time.time() - t0
        if not plot_res:
            logger.error("%s: plot check failed (%s) (%.2fs)", display_name, plot_reason, dt)
            results.append((display_name, "FAILED", f"plot check failed: {plot_reason}", dt))
        else:
            logger.info("%s: PASSED (%.2fs)", display_name, dt)
            results.append((display_name, "PASSED", "Passed", dt))

    except ImportError as exc:
        dt = time.time() - t0
        # A check that imports a third-party package (e.g. numpy) inside a
        # function only finds out here that it is missing.
        package = missing_dependency(exc)
        if package is None:
            raise
        logger.warning("SKIP: %s needs the Python package '%s', which is not "
                       "installed. (%.2fs)", display_name, package, dt)
        results.append((display_name, "SKIPPED",
                        f"needs the Python package '{package}'", dt))

    finally:
        logger.debug("  Cleaning up output for %s...", display_name)
        cleanup_run_dir(run_dir)


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
AMR_TESTS = {
    "amr_equilibrium", "dynamic_amr", "lightwave", "reconnection_amr",
    "wave_amr_3d",
}

# Tests needing a true-3D AMReX library (real z-grid, nCellZ > 1).
AMREX3D_TESTS = {"lightwave", "wave_amr_3d"}

# Tests needing a true-2D AMReX library (built with -amrex2d).
AMREX2D_TESTS = {"photoionization", "body"}

# Tests excluded from the 2D suite.
AMREX2D_EXCLUDED_TESTS = {
    "chemistry",
}


def _read_amrex_spacedim(path):
    """Return AMREX_SPACEDIM from an AMReX_Config.H file, or None.

    Robust to a broken "InstallDir" symlink (which GitHub Actions cache does not
    reliably preserve): even if InstallDir/... is unreadable, we can still read
    the dimension-specific InstallDir2D/ or InstallDir3D/ directory directly.
    """
    import re
    try:
        with open(path, "r") as f:
            content = f.read()
    except (OSError, IOError):
        return None
    m = re.search(r"#define\s+AMREX_SPACEDIM\s+(\d+)", content)
    return int(m.group(1)) if m else None


def configured_amrex_dim():
    """AMReX dimension (2 or 3) the binary was built with, or None.

    Reads AMREX_SPACEDIM from the linked library's AMReX_Config.H in either the
    standalone (util/AMREX) or SWMF-component (../../util/AMREX) layout,
    checking the InstallDir symlink and its 2D/3D targets so a broken symlink
    does not prevent detection.
    """
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # (path-to-AMReX-root, candidates-to-check) per layout.
    layouts = [
        os.path.join(root, "util", "AMREX"),
        os.path.join(root, "..", "..", "util", "AMREX"),
    ]
    for base in layouts:
        for sub in ("InstallDir", "InstallDir3D", "InstallDir2D"):
            dim = _read_amrex_spacedim(
                os.path.join(base, sub, "include", "AMReX_Config.H"))
            if dim is not None:
                return dim
    return None


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

    Returns (ok, reason, skip).  ok=False means do not run: skip=True -> record
    as SKIPPED (needs a different AMReX dimension); skip=False -> FAILED.
    """
    if test_name in EXO_SOURCE_TESTS and not configured_user_source_is_exo():
        return False, (
            f"Test '{test_name}' needs the Exosphere user source, but "
            "include/UserSource.h is the default empty source.\n"
            "  Rebuild with:  ./Config.pl -lev=2 -u=Exo && make"
        ), False

    if test_name in AMR_TESTS:
        nlev = configured_nlevmax()
        if nlev is None:
            return False, (
                f"Test '{test_name}' needs AMR (nLevMax >= 2), but "
                "include/Constants.h could not be read.\n"
                "  Rebuild with:  ./Config.pl -lev=2 -u=Exo && make"
            ), False
        if nlev < 2:
            return False, (
                f"Test '{test_name}' needs AMR (nLevMax >= 2), but the "
                f"binary is configured with nLevMax = {nlev}.\n"
                "  Rebuild with:  ./Config.pl -lev=2 -u=Exo && make"
            ), False

    dim = configured_amrex_dim()
    required = None
    if test_name in AMREX2D_TESTS:
        required = 2
    elif test_name in AMREX3D_TESTS:
        required = 3
    if required is not None:
        if dim is None:
            return False, (
                f"Test '{test_name}' needs a {required}D AMReX library, but the "
                "configured AMReX dimension could not be determined "
                "(util/AMREX/InstallDir/include/AMReX_Config.H missing).\n"
                f"  Rebuild with:  ./Config.pl -amrex{required}d -lev=2 -u=Exo && make"
            ), False
        if dim != required:
            return False, f"needs a {required}D AMReX library (binary is {dim}D)", True

    # Skip tests excluded from the 2D suite under a 2D build.
    if dim == 2 and test_name in AMREX2D_EXCLUDED_TESTS:
        return False, "not part of the 2D AMReX suite", True

    return True, None, False


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

    header = [col.strip() for col in lines[0].strip().split("\t")]
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
            # Optional columns appended after the fixed ones (e.g. the #BODY
            # absorption tallies) are exposed under their header name.
            for iCol in range(6 + n_species, min(len(vals), len(header))):
                name = header[iCol]
                if name and name not in entry:
                    try:
                        entry[name] = float(vals[iCol])
                    except ValueError:
                        continue
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
              expected_active_species (list[int], default []),
              min_velocity_change (float, default None = disabled).
    """
    import math
    tol = tol or {}
    min_rows = int(tol.get("min_rows", 1))
    launch_threshold = float(tol.get("launch_threshold", 0.5))
    max_speed = float(tol.get("max_speed", 10.0))
    expected_active = set(int(s) for s in tol.get("expected_active_species", []))
    # Optional: require an active species to show a relative mean-velocity
    # change >= this (catches a silently disabled field gather).
    min_velocity_change = tol.get("min_velocity_change", None)
    if min_velocity_change is not None:
        min_velocity_change = float(min_velocity_change)

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

    # Field-interaction check: catch a zero-field gather (constant velocity).
    if min_velocity_change is not None and n_rows >= 2:
        got_change = False
        for iS in sorted(active):
            m0 = pt_diags[0].get(f"mass{iS}", 0.0)
            m1 = pt_diags[-1].get(f"mass{iS}", 0.0)
            if m0 <= 0 or m1 <= 0:
                continue
            v0 = [pt_diags[0].get(f"moment_{ax}{iS}", 0.0) / m0
                  for ax in ("x", "y", "z")]
            v1 = [pt_diags[-1].get(f"moment_{ax}{iS}", 0.0) / m1
                  for ax in ("x", "y", "z")]
            speed0 = math.sqrt(sum(c * c for c in v0))
            dv = math.sqrt(sum((a - b) ** 2 for a, b in zip(v1, v0)))
            if speed0 > 1e-30 and dv / speed0 >= min_velocity_change:
                got_change = True
                break
        if not got_change:
            passed = False
            reasons.append(
                f"no active species shows a mean-velocity change >= "
                f"{min_velocity_change:g} (relative); field gather may be "
                "silently disabled")

    if passed:
        logger.debug("Test-Particle Tracer Output: PASSED")
        return True, "Passed"
    return False, "; ".join(reasons)


# ---------------------------------------------------------------------------
# Per-test validator loading
# ---------------------------------------------------------------------------
class MissingDependency(Exception):
    """A validator needs a third-party Python package that is not installed."""

    def __init__(self, package):
        super().__init__(package)
        self.package = package


def missing_dependency(exc):
    """Return the top-level package name if *exc* is an ImportError for a
    third-party package (e.g. numpy), or None if the missing module belongs to
    the test tree itself (tests.*, or a module next to this script such as
    _shared), which is a bug that must surface."""
    if not isinstance(exc, ImportError) or not exc.name:
        return None
    top = exc.name.split(".")[0]
    here = os.path.dirname(os.path.abspath(__file__))
    if top == "tests" or os.path.isdir(os.path.join(here, top)) or \
            os.path.isfile(os.path.join(here, top + ".py")):
        return None
    return top


def load_validator(test_name):
    """Dynamically import the ``validate.py`` module for *test_name*.

    Returns a small object with attributes ``validate_log``, ``plot`` and
    ``particle_tol`` (some may be None) so the runner can invoke it uniformly.
    If the test has no ``validate.py`` (e.g. ``hyper_resistivity``), returns
    None and the runner falls back to a generic no-op check.  If the module
    needs a third-party package that is not installed (e.g. numpy), raises
    MissingDependency so the runner reports the test as SKIPPED instead of
    passing it with the generic check.
    """
    module_name = f"tests.{test_name}.validate"
    try:
        module = importlib.import_module(module_name)
    except ImportError as exc:
        # Only a genuinely missing per-test module (hyper_resistivity and
        # friends) falls back to the generic check.
        if exc.name in (module_name, f"tests.{test_name}"):
            logger.debug("No per-test module %s; using generic check.",
                         module_name)
            return None
        package = missing_dependency(exc)
        if package is not None:
            raise MissingDependency(package) from exc
        raise

    class _Validator:
        pass

    v = _Validator()
    v.validate_log = getattr(module, "validate_log", None)
    # plot() expects a single positional arg (test_name).
    v.plot = getattr(module, "validate_plot", None)
    v.particle_tol = getattr(module, "PARTICLE_TOL", None)
    # Kept so the run directory can be pushed into the module (see
    # _sync_validator_run_dir): several per-test modules own a module-level
    # RUN_DIR that they use to locate plot output.
    v.module = module
    return v


def _sync_shared_run_dir(run_dir=None):
    """Push the current RUN_DIR into the shared helpers used by plot checks."""
    if run_dir is None:
        run_dir = RUN_DIR
    try:
        import tests._shared.hybrid as _hyb
        _hyb.set_run_dir(run_dir)
    except ImportError:
        pass
    try:
        import tests._shared.run_dir as _rd
        _rd.set_run_dir(run_dir)
    except ImportError:
        pass


def _sync_validator_run_dir(validator, run_dir=None):
    """Push the current RUN_DIR into a per-test validator module.

    Modules that read plot output keep their own module-level ``RUN_DIR``
    (e.g. ``tests/shock/validate.py``).  Without this, ``--run-dir DIR`` would
    leave them pointing at the default ``run_test/`` and their plot check
    would silently find no frames and report "skipped" instead of validating.
    """
    if run_dir is None:
        run_dir = RUN_DIR
    module = getattr(validator, "module", None) if validator else None
    setter = getattr(module, "set_run_dir", None)
    if callable(setter):
        setter(run_dir)


# ---------------------------------------------------------------------------
# Test discovery + main
# ---------------------------------------------------------------------------
def discover_tests(tests_dir="tests"):
    """Return a sorted list of (test_dir, name) for directories with any
    PARAM.in variant (PARAM.in and/or PARAM.in.<suffix>).

    Excludes "performance" and "run_test" (the shared run directory, which
    would otherwise be re-run as a test); dev tests live outside tests/.
    """
    subdirs = sorted([d for d in os.listdir(tests_dir)
                      if os.path.isdir(os.path.join(tests_dir, d))
                      and d not in ["performance", "run_test"]])
    tests = []
    for d in subdirs:
        test_dir = os.path.join(tests_dir, d)
        if glob.glob(os.path.join(test_dir, "PARAM.in*")):
            tests.append((test_dir, d))
    return tests


# Variants or tests that are computationally expensive and excluded from the
# default test suite. They can be executed explicitly via `--test=NAME.VARIANT`
# or by passing `--all` / `--include-expensive`.
EXPENSIVE_VARIANTS = {
    # Only the force-free reconnection variants run in the default suite. The
    # Fadeev, GEM and asymmetric ones are research benchmarks: keep them
    # runnable explicitly (or with --all) so CI stays fast.
    "reconnection.fadeev_hybrid",
    "reconnection.fadeev_hybrid_long",
    "reconnection.fadeev_pic",
    "reconnection.gem_hybrid",
    "reconnection.gem_pic",
    "reconnection.asym_pic",
    "beam.instability",
    "iaw.landau",
}


def is_expensive(name, variant_token, test_dir=None):
    """Check whether a test variant is marked as expensive."""
    full_token = f"{name}.{variant_token}"
    if full_token in EXPENSIVE_VARIANTS or name in EXPENSIVE_VARIANTS:
        return True
    if test_dir is not None:
        expensive_marker = os.path.join(test_dir, ".expensive")
        if os.path.isfile(expensive_marker):
            try:
                with open(expensive_marker) as f:
                    marked = [
                        ln.strip()
                        for ln in f
                        if ln.strip() and not ln.startswith("#")
                    ]
                if not marked or variant_token in marked or name in marked:
                    return True
            except Exception:
                return True
    return False


def run_one_test(test_dir, name, nprocs, results, variant_filter=None,
                 include_expensive=False, run_dir=None):
    """Pre-flight check + run one test (or its PARAM variants) + record result.

    *variant_filter* (optional) is a PARAM suffix (without the leading dot, e.g.
    "hybrid") that restricts execution to that single variant of the test.  When
    None every non-expensive variant under *test_dir* runs (unless include_expensive=True).
    """
    if run_dir is None:
        run_dir = RUN_DIR

    ok, preflight_reason, skip = preflight_check(name)
    if not ok:
        if skip:
            # Dimension-mismatched test: skip silently, no detail needed.
            results.append((name.upper(), "SKIPPED", "", 0.0))
        else:
            logger.warning("SKIP: %s cannot run with the configured build.", name.upper())
            logger.warning(preflight_reason)
            results.append((name.upper(), "FAILED",
                            f"build misconfigured: {preflight_reason}", 0.0))
        return

    missing_package = None
    try:
        validator = load_validator(name)
    except MissingDependency as exc:
        validator = None
        missing_package = exc.package

    # A test directory may hold several parameter files (PARAM.in plus any
    # PARAM.in.<suffix>) that group related test configurations (e.g. the same
    # physics under different field solvers, or several wave kinds).  The
    # runner executes the test once per variant.  Each variant gets a distinct
    # base_name built from the suffix (PARAM.in -> name, PARAM.in.hybrid ->
    # name_hybrid, PARAM.in.foo -> name_foo) so the single per-directory
    # validate.py can branch on base_name and apply variant-appropriate checks.
    variants = []
    param_files = sorted(glob.glob(os.path.join(test_dir, "PARAM.in*")))
    for pf in param_files:
        suffix = os.path.basename(pf)[len("PARAM.in"):]  # "" | ".hybrid" | ...
        if not suffix:
            display_name = name.upper()
            base_name = name
        else:
            tag = suffix.lstrip(".").upper()
            display_name = f"{name.upper()} ({tag})"
            base_name = f"{name}{suffix.replace('.', '_')}"
        variants.append((pf, display_name, base_name))

    if name == "freestream":
        # Free-stream test variants: full (original GMRES), upwind, hybrid, 2D
        # hybrid, and the 2D hybrid inflow-layer regression test.
        variants = [
            (os.path.join(test_dir, "PARAM.in.full"),
             "FREESTREAM (FULL PIC GMRES)", "freestream_full"),
            (os.path.join(test_dir, "PARAM.in.upwind"),
             "FREESTREAM (FULL PIC UPWIND)", "freestream_upwind"),
            (os.path.join(test_dir, "PARAM.in.hybrid"),
             "FREESTREAM (HYBRID HALL-OFF)", "freestream_hybrid"),
            (os.path.join(test_dir, "PARAM.in.hybrid_2d"),
             "FREESTREAM (HYBRID 2D INFLOW)", "freestream_hybrid_2d"),
            (os.path.join(test_dir, "PARAM.in.hybrid_inflow"),
             "FREESTREAM (HYBRID INFLOW LAYER)", "freestream_hybrid_inflow"),
        ]

    # Restrict to a single PARAM variant when requested.  The variant token is
    # the filename suffix without the leading dot, with the base file (no
    # suffix) mapped to "full": `PARAM.in` -> "full", `PARAM.in.hybrid` ->
    # "hybrid".  This lets `--test=NAME.full` pick the base PARAM.in and
    # `--test=NAME.hybrid` pick the hybrid variant.
    def _variant_token(pf):
        suffix = os.path.basename(pf)[len("PARAM.in"):].lstrip(".")
        return suffix or "full"

    if name == "freestream" and variant_filter == "gmres":
        variant_filter = "full"

    if name == "whistler" and variant_filter == "hybrid":
        variant_filter = "rk4"

    if name == "bc_reflecting" and variant_filter is not None:
        # Backward compatibility for legacy --test=bc_reflecting.full and .hybrid
        if variant_filter == "full":
            variant_filter = "particles"
        elif variant_filter == "hybrid":
            variant_filter = "hybrid.particles"

    if variant_filter == "full":
        has_full = any(_variant_token(pf) == "full" for pf, _, _ in variants)
        has_implicit = any(_variant_token(pf) == "implicit" for pf, _, _ in variants)
        if not has_full and has_implicit:
            variant_filter = "implicit"
    elif variant_filter == "implicit":
        has_implicit = any(_variant_token(pf) == "implicit" for pf, _, _ in variants)
        has_full = any(_variant_token(pf) == "full" for pf, _, _ in variants)
        if not has_implicit and has_full:
            variant_filter = "full"

    if variant_filter is not None:
        kept = []
        for pf, display_name, base_name in variants:
            if _variant_token(pf) == variant_filter:
                kept.append((pf, display_name, base_name))
        if not kept:
            available = ", ".join(sorted({_variant_token(pf)
                                          for pf, _, _ in variants}))
            logger.error(
                "Error: no variant '%s' in test '%s'. Available variants: %s",
                variant_filter, name, available)
            results.append((name.upper(), "FAILED",
                            f"no variant '{variant_filter}'", 0.0))
            return
        variants = kept
        logger.debug("Restricted test '%s' to variant '%s'.", name,
                     variant_filter)
    elif not include_expensive:
        kept = []
        for pf, display_name, base_name in variants:
            v_tok = _variant_token(pf)
            if is_expensive(name, v_tok, test_dir):
                logger.info("  [INFO] Skipping expensive variant '%s.%s' (run with --test=%s.%s or --all)",
                            name, v_tok, name, v_tok)
            else:
                kept.append((pf, display_name, base_name))
        variants = kept

    for param_file, display_name, base_name in variants:
        if missing_package is not None:
            logger.warning("SKIP: %s needs the Python package '%s', which is "
                           "not installed.", display_name, missing_package)
            results.append((display_name, "SKIPPED",
                            f"needs the Python package '{missing_package}'", 0.0))
            continue
        with open(param_file) as _f:
            _param_text = _f.read()
        run_and_validate(test_dir, display_name, validator, nprocs, results,
                         param_text=_param_text, base_name=base_name, run_dir=run_dir)


def set_worker_run_dir(run_dir):
    """Mirror run_dir into runner globals and shared test helpers for a worker."""
    global RUN_DIR
    RUN_DIR = run_dir
    _sync_shared_run_dir(run_dir)


def _worker_execute_test(task):
    """Worker task executed in a subprocess by ProcessPoolExecutor."""
    test_idx, test_dir, name, nprocs, variant_filter, include_expensive, queue, verbose, base_run_dir = task
    slot = queue.get()
    worker_dir = f"{base_run_dir}_w{slot}"

    # Configure worker-local logger to capture output without interleaving
    import io
    log_stream = io.StringIO()
    handler = logging.StreamHandler(log_stream)
    handler.setFormatter(logging.Formatter("%(message)s"))
    worker_logger = logging.getLogger()
    worker_logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    worker_logger.handlers = [handler]

    set_worker_run_dir(worker_dir)

    res = []
    try:
        run_one_test(test_dir, name, nprocs, res, variant_filter=variant_filter,
                     include_expensive=include_expensive, run_dir=worker_dir)
    except Exception as e:
        res.append((name.upper(), "FAILED", f"Worker execution error: {e}", 0.0))
        import traceback
        traceback.print_exc(file=log_stream)
    finally:
        cleanup_run_dir(worker_dir)
        queue.put(slot)

    return test_idx, res, log_stream.getvalue()


def main():
    # Parse nprocs: -n N or --nprocs N
    nprocs = 1
    for i, arg in enumerate(sys.argv):
        if arg in ("-n", "--nprocs") and i + 1 < len(sys.argv):
            try:
                nprocs = int(sys.argv[i + 1])
            except ValueError:
                logger.error("Error: %s requires an integer argument (number of MPI processes).", arg)
                sys.exit(1)
            continue

    if nprocs < 1:
        logger.error("Error: Number of processes must be >= 1.")
        sys.exit(1)

    # Parse -j / --jobs: concurrent test execution
    jobs = 1
    for i, arg in enumerate(sys.argv):
        if arg.startswith("--jobs="):
            val = arg[len("--jobs="):]
            if val.lower() == "auto":
                jobs = max(1, (os.cpu_count() or 1) // nprocs)
            else:
                try:
                    jobs = int(val)
                except ValueError:
                    logger.error("Error: --jobs requires an integer or 'auto'.")
                    sys.exit(1)
        elif arg in ("-j", "--jobs"):
            if i + 1 < len(sys.argv) and not sys.argv[i + 1].startswith("-"):
                val = sys.argv[i + 1]
                if val.lower() == "auto":
                    jobs = max(1, (os.cpu_count() or 1) // nprocs)
                else:
                    try:
                        jobs = int(val)
                    except ValueError:
                        logger.error("Error: -j requires an integer or 'auto'.")
                        sys.exit(1)
            else:
                jobs = max(1, (os.cpu_count() or 1) // nprocs)

    if jobs < 1:
        logger.error("Error: Number of jobs must be >= 1.")
        sys.exit(1)

    # Parse --shard K/N: partition tests across CI jobs
    shard_str = None
    for i, arg in enumerate(sys.argv):
        if arg.startswith("--shard="):
            shard_str = arg[len("--shard="):]
        elif arg == "--shard" and i + 1 < len(sys.argv):
            shard_str = sys.argv[i + 1]

    # Parse --summary-file PATH (custom output path for CI serial/parallel split)
    summary_file = None
    for i, arg in enumerate(sys.argv):
        if arg == "--summary-file" and i + 1 < len(sys.argv):
            summary_file = os.path.abspath(sys.argv[i + 1])
            continue
        if arg.startswith("--summary-file="):
            summary_file = os.path.abspath(arg[len("--summary-file="):])
            continue

    # Parse --summary-title TITLE (custom heading for Markdown summary)
    summary_title = "Standalone FLEKS Test Results"
    for i, arg in enumerate(sys.argv):
        if arg == "--summary-title" and i + 1 < len(sys.argv):
            summary_title = sys.argv[i + 1]
            continue
        if arg.startswith("--summary-title="):
            summary_title = arg[len("--summary-title="):]
            continue

    # Parse --test NAME (select a specific test to run; default: run all)
    # Accepts both "--test=NAME" and "--test NAME" forms.
    selected_test = None
    for i, arg in enumerate(sys.argv):
        if arg.startswith("--test="):
            selected_test = arg[len("--test="):]
            continue
        if arg == "--test" and i + 1 < len(sys.argv):
            selected_test = sys.argv[i + 1]
            continue
        if arg == "--run-dir" and i + 1 < len(sys.argv):
            global RUN_DIR
            RUN_DIR = sys.argv[i + 1]
            continue
        if arg.startswith("--run-dir="):
            RUN_DIR = arg[len("--run-dir="):]
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

    if summary_file is None:
        summary_file = os.path.abspath("tests/summary.md")

    _sync_shared_run_dir()

    tests = discover_tests()
    if not tests:
        logger.error("No tests found in tests/ subdirectories!")
        sys.exit(1)

    # Filter to a single test if --test was given.
    variant_filter = None
    if selected_test is not None:
        if "." in selected_test:
            name_part, variant_filter = selected_test.split(".", 1)
        else:
            name_part = selected_test
        matching = [t for t in tests if t[1] == name_part]
        if not matching:
            available = ", ".join(t[1] for t in tests)
            logger.error("Error: test '%s' not found.", selected_test)
            logger.error("Available tests: %s", available)
            sys.exit(1)
        tests = matching
        logger.debug("Selected test: %s", selected_test)

    # Filter by --shard K/N if provided
    shard_info = ""
    if shard_str is not None:
        try:
            k_str, n_str = shard_str.split("/")
            shard_k, shard_n = int(k_str), int(n_str)
            if shard_k < 1 or shard_k > shard_n:
                raise ValueError()
        except ValueError:
            logger.error("Error: --shard must be in format K/N (e.g. 1/2) with 1 <= K <= N.")
            sys.exit(1)
        tests = [t for idx, t in enumerate(tests) if idx % shard_n == (shard_k - 1)]
        shard_info = f" (Shard {shard_k}/{shard_n})"
        logger.info("Running shard %d/%d (%d tests selected)", shard_k, shard_n, len(tests))

    # Parse --all / --include-expensive: run every variant including expensive benchmarks.
    include_expensive = "--all" in sys.argv or "--include-expensive" in sys.argv

    results = []  # Collect results for summary table

    if jobs > 1 and len(tests) > 1:
        logger.info("Running standalone tests concurrently with %d workers (nprocs=%d)...", jobs, nprocs)
        ensure_postidl()
        import multiprocessing
        from concurrent.futures import ProcessPoolExecutor, as_completed
        import shutil

        manager = multiprocessing.Manager()
        slot_queue = manager.Queue()
        for s in range(jobs):
            slot_queue.put(s)

        HEAVY_TESTS = {
            "shock", "whistler", "reconnection_amr", "reconnection",
            "langmuir", "freestream", "body", "opticaldepth"
        }

        indexed_tests = list(enumerate(tests))
        # Prioritize heavier tests first to maximize multi-core throughput
        scheduled_tests = sorted(
            indexed_tests,
            key=lambda item: (0 if item[1][1] in HEAVY_TESTS else 1, item[0])
        )

        tasks = [
            (idx, td, name, nprocs, variant_filter, include_expensive, slot_queue, verbose, RUN_DIR)
            for idx, (td, name) in scheduled_tests
        ]

        collected_results = {}
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            future_to_task = {executor.submit(_worker_execute_test, t): t[2] for t in tasks}
            for future in as_completed(future_to_task):
                test_idx, test_results, log_output = future.result()
                if log_output.strip():
                    print(log_output, end="", flush=True)
                collected_results[test_idx] = test_results

        # Flatten results in original deterministic order
        for idx, _ in indexed_tests:
            if idx in collected_results:
                results.extend(collected_results[idx])

        # Clean up worker directories
        for s in range(jobs):
            wdir = f"{RUN_DIR}_w{s}"
            if os.path.exists(wdir):
                shutil.rmtree(wdir, ignore_errors=True)
    else:
        for test_dir, name in tests:
            run_one_test(test_dir, name, nprocs, results, variant_filter=variant_filter,
                         include_expensive=include_expensive, run_dir=RUN_DIR)

    # ----------------------------------------------------
    # Print Summary Table
    # ----------------------------------------------------
    logger.info("")
    logger.info("=" * 80)
    logger.info(" " * 32 + f"TEST SUMMARY{shard_info}")
    logger.info("=" * 80)
    logger.info(" %-28s | %-8s | %-8s | %s", "Test Name", "Status", "Time", "Failure Reason / Details")
    logger.info("-" * 80)

    all_passed = True
    for item in results:
        name_str, status, reason = item[0], item[1], item[2]
        dt = item[3] if len(item) > 3 else 0.0
        time_display = f"{dt:.2f}s" if dt > 0 else "-"
        if status == "PASSED":
            status_display = "\033[92;1mPASSED\033[0m" if sys.stdout.isatty() else "PASSED"
        elif status == "SKIPPED":
            status_display = "\033[93;1mSKIPPED\033[0m" if sys.stdout.isatty() else "SKIPPED"
        else:
            status_display = "\033[91;1mFAILED\033[0m" if sys.stdout.isatty() else "FAILED"
            all_passed = False
        reason_display = reason if status != "PASSED" else ""
        logger.info(" %-28s | %-8s | %-8s | %s", name_str, status_display, time_display, reason_display)

    logger.info("=" * 80)

    # Write Markdown Summary to summary.md for CI / PR integration
    try:
        with open(summary_file, "w") as f:
            header_prefix = summary_title if summary_title.startswith("#") else f"### 🧪 {summary_title}"
            f.write(f"{header_prefix}{shard_info}\n\n")
            f.write("| Test Name | Status | Time | Failure Reason / Details |\n")
            f.write("| :--- | :--- | :--- | :--- |\n")
            for item in results:
                name_str, status, reason = item[0], item[1], item[2]
                dt = item[3] if len(item) > 3 else 0.0
                time_md = f"{dt:.2f}s" if dt > 0 else "-"
                if status == "PASSED":
                    status_md = "🟢 **PASSED**"
                elif status == "SKIPPED":
                    status_md = "🟡 **SKIPPED**"
                else:
                    status_md = "🔴 **FAILED**"
                reason_md = reason if status != "PASSED" else ""
                f.write(f"| {name_str} | {status_md} | {time_md} | {reason_md} |\n")
    except Exception as e:
        logger.warning("Warning: Could not write summary file %s: %s", summary_file, e)

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
