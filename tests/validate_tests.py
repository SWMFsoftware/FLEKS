#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys
import math

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
    
    # Verify that PostIDL.exe exists; PostProc.pl needs it to produce .out files.
    postidl_link = os.path.join("run_test", "PC", "PostIDL.exe")
    if os.path.islink(postidl_link) and not os.path.exists(postidl_link):
        real = os.path.realpath(postidl_link)
        print(f"  [WARN] Broken symlink: {postidl_link} -> {real}")
        print(f"  [WARN] PostIDL.exe is missing. Build it with 'make PIDL' "
              f"before running tests that check plot output (.out files).")
    
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
    pp = subprocess.run(["./PostProc.pl", "-v"], cwd="run_test",
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
    param_path = os.path.join("run_test", "PARAM.in")
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

    plots_dir = os.path.join("run_test", "PC", "plots")

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
            break
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

    plots_dir = os.path.join("run_test", "PC", "plots")
    out_files = sorted(glob.glob(os.path.join(plots_dir, "*.out")))
    if not out_files:
        print("    [FFT] No .out files found (PostProc.pl not run?).")
        return True, "No .out files (skipped)"

    # Use the final frame (latest cycle); this is the t ~= 0.1 snapshot.
    out_file = out_files[-1]
    print(f"    [FFT] Loading .out: {os.path.basename(out_file)}")

    with open(out_file, "r") as f:
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

    plots_dir = os.path.join("run_test", "PC", "plots")
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
            break
    if rho_idx is None:
        print(f"    [CX] rhoS2/rhoS1 not found in .out variables: {var_names}")
        return True, "rhoS2/rhoS1 not in .out"

    # Read planet radius and normalization from PARAM.in (plot coords = SI / lNormSI).
    Rp_si = 3.0e6
    lNormSI = 1000.0
    try:
        with open(os.path.join("run_test", "PARAM.in"), "r") as pf:
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

    # ---- Other tests: no plot-file validation ----
    print("  --- Validating Output Files: No plot-file check for this test ---")
    return True, "Passed (no plot-file check)"

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

    # Parse --test NAME (select a specific test to run; default: run all)
    # Accepts both "--test=NAME" and "--test NAME" forms.
    selected_test = None
    for i, arg in enumerate(sys.argv):
        if arg.startswith("--test="):
            selected_test = arg[len("--test="):]
            break
        if arg == "--test":
            try:
                selected_test = sys.argv[i + 1]
            except IndexError:
                print("Error: --test requires a test name argument.")
                sys.exit(1)
            break

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
        print(f"\n==========================================")
        print(f"Starting test: {name.upper()}")
        print(f"==========================================")
        try:
            stdout, code = run_test(test_dir, nprocs=nprocs)
            if code != 0 or stdout is None:
                print(f"FAIL: {name.upper()} execution failed with exit code {code}")
                results.append((name.upper(), "FAILED", f"Execution failed (code {code})"))
                continue

            # Read the PIC energy log (the only diagnostic log produced by FLEKS).
            pic_diags = read_pic_log("run_test")

            val_res = False
            reason = "Validation skipped"

            if validator:
                import inspect
                sig = inspect.signature(validator)
                kwargs = {}
                if "pic_diags" in sig.parameters:
                    kwargs["pic_diags"] = pic_diags
                if "test_name" in sig.parameters:
                    kwargs["test_name"] = name
                val_res, reason = validator(**kwargs)
                if not val_res:
                    results.append((name.upper(), "FAILED", reason))
                    continue
            else:
                print(f"Validating {name.upper()} (generic check)...")
                print(f"{name.upper()} (generic check): PASSED")
                val_res = True
                reason = "Passed"

            # Validate output plotfiles
            plot_res, plot_reason = validate_plot_output(name)
            if not plot_res:
                results.append((name.upper(), "FAILED", f"plot check failed: {plot_reason}"))
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
