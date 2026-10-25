# GPU PIC demo validation

The test executable returns a nonzero status when a check or field solve fails.
`validate_demo.py` runs the public modes in separate processes and checks the
resulting diagnostics, rather than relying only on printed PASS labels. It
uses the Python standard library.

## Build and run

```bash
make GPUTEST GPUTEST_REFERENCE -j8
python3 src_gpu/validate_demo.py --mpi --reference bin/gpu_reference.exe \
    --output-dir /tmp/fleks-gpu-validation
python3 src_gpu/validate_demo.py --physics \
    --output-dir /tmp/fleks-gpu-physics
```

Each case writes a log; `results.json` retains its exact command and diagnostics.
Without `--output-dir`, the runner uses and removes a temporary directory.
One-rank cases use `mpirun -np 1` by default, avoiding singleton-startup issues
in some MPI installations. Use `--direct` when launching directly is preferred,
`--mpiexec 'mpirun --oversubscribe'` to customize the launcher, and `--timeout`
to change the per-case time limit. All children inherit the calling environment.

On the development host, the relocated NVIDIA HPC SDK requires the installed
OpenMPI prefix; host-staged transport avoids device-buffer assumptions:

```bash
export OPAL_PREFIX=/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/comm_libs/12.9/hpcx/hpcx-2.22.1/ompi
export OMPI_MCA_btl=self,vader
export OMPI_MCA_pml=ob1
```

These are host-specific runtime settings, not build requirements. Running a
CUDA executable requires access to the NVIDIA driver and may need execution
outside a restrictive sandbox. Compilation alone is not GPU verification.

## Default suite

- Run `fields`, `particles`, and `deposit`, and optionally repeat on two ranks.
- Repeat particle and deposition checks with one and 27 particles per cell to
  cover different lattice sizes.
- Compare a one-step uniform plasma case against an actual production FLEKS
  CPU result in `production_one_step.json`.
- Compare global physical diagnostics over eight coupled steps between box
  sizes four and eight, and optionally between one and two MPI ranks.
- Repeat the box and MPI comparisons on a nonuniform longitudinal wave
  (`perturbation=0.01`) with div(E) correction enabled. The uniform case alone
  cannot detect stencil or seam errors.
- Rerun that nonuniform case with `response=cell` and require its diagnostics
  to match the per-particle atomic assembly. The `particles` mode also checks
  the cell-reduced `jHat` and mass matrix against the independent serial
  reference (oblique B, shared nodes, periodic corners), and `benchmark` mode
  compares all three assembly methods entry by entry.
- A div(E) solve that reaches `div_e_max_iter` prints `WARN` and the partial
  potential is applied, as in production `LinearSolver`; it is not a failure.
  An electric solve that does not converge still fails the run.
- Exercise uniform plasma with div(E) correction enabled and disabled.
- Run a nonuniform longitudinal perturbation with correction off/on and require
  the corrected Gauss residual to fall below `1e-6` times the uncorrected value
  plus `1e-12` (measured: about 1e-4 to 1e-14).
- Run 20 steps at `dt=0.5` (`c*dt/dx=8`, `n_cell=16`, `b0=0.5`) with
  correction enabled, the regime FLEKS targets. Require at least ten GMRES
  iterations in some step, so the case keeps exercising a nontrivial solve,
  and a corrected Gauss residual below `1e-7` after step 0. Optionally
  compare one and two ranks with relative tolerance `1e-6`.
- Require finite diagnostics, increasing steps/times, constant particle count,
  integrated neutral charge within `1e-9`, nonnegative energy, and the true
  electric-solver residual below `1e-7`.

Decomposition comparisons use relative tolerance `2e-7` and absolute tolerance
`1e-14`. These compare energy, charge, the minimum nodal `E_x`, particle
count, time, and step. Krylov
iteration counts and near-zero residual maxima may depend on accumulation order
and are intentionally not required to match bitwise.

## Exact production CPU operator comparison

`make GPUTEST_REFERENCE` generates `ProductionOperators.h` from the current
`src/GridUtility.cpp`. Extraction preserves the numerical function bodies and
records the source SHA256. Only dispatch changes: `ParallelFor` becomes
`LoopOnCpu`, and the device-only lambda annotation is removed. The reference
uses CPU-accessible pinned storage. Production source files remain untouched,
and the reference tool does not link `libFLEKS.a`.

`ProductionReference.cpp` compares the demo backend with those CPU bodies for:

1. Nodal-to-cell and cell-to-nodal gradients.
2. Nodal-to-cell and cell-to-nodal divergences.
3. Nodal-to-cell and cell-to-nodal curls.
4. The complete spatial electric-field operator, with zero particle response,
   at `coef_diff=0` and `coef_diff=0.1`.

The composed test calls the actual `DemoFields::apply_e_operator`; its CPU
reference uses production gradients/divergences, including the 3-by-3 transverse
averaging in production's extended cell-centered divergence. The manufactured
input contains mixed three-dimensional Fourier modes so cross terms matter.
The tolerance is `1e-12 * max(1, max(abs(reference)))`.

The tool writes paired files `<reference_prefix>.rankN.production.csv` and
`<reference_prefix>.rankN.demo.csv` with columns
`operator,box,i,j,k,component,value`. This preserves independently calculated
arrays for external inspection. Disable array output with `write_files=0`.
The six primitive operators passed on the development GPU with maximum absolute
error `7.11e-15`; the actual composed operator passed at both blend values with
maximum absolute error `1.67e-15`.

## Production CPU coupled one-step fixture

`production_one_step.json` contains an actual CPU FLEKS run's parameter deck,
binary SHA256, printed version header, and initial/final diagnostics. The deck
uses a four-cubed periodic mesh, a quiet two-by-two-by-two particle lattice,
ion/electron masses 25 and one, a uniform electron drift of `0.001`, zero ion
velocity and fields, `dt=0.01`, `theta=0.51`, `coefDiff=0.1`, and disabled div(E)
correction. Temperatures and resampling are zero/off. The production build is
identified by its own version header; this is not asserted to be a fresh build
of the current working tree.

The SI normalization inputs give number density `1.0000000259030656` in code
units rather than exactly one because the supplied decimal density is rounded.
The electric-field plot has eight significant digits. The comparison therefore
uses relative tolerance `2e-7` and absolute tolerance `1e-13` for the final
uniform electric field, field energy, and particle energy.

Observed production CPU one-step values:

| Quantity | Value |
| --- | --- |
| Initial kinetic energy | `5.000000129515328e-7` |
| Final kinetic energy | `4.993595549663868e-7` |
| Final electric energy | `6.278999854368504e-10` |
| Final uniform Ex | `1.2562184e-4` |

A uniform case tests normalization, response amplitude, charge sign, and time
centering. It does not exercise spatially varying particle-response couplings;
those require the separate tiny-case mass-matrix/reference checks. No binary
or large plot output is committed.

## Longer physics checks

`--physics` runs the uniform cold-plasma oscillator at `theta=0.5`, with
correction disabled, for 600 steps at `dt=0.01`. For normalized density one,
electron mass one, and ion mass 25, the expected frequency is
`sqrt(4*pi*(1+1/25))`. Successive maxima in electric energy occur at twice the
field frequency; the estimator accounts for that factor and interpolates peak
locations. It requires at least three peaks, frequency error below `1e-3`
(measured `1.09e-4`), and total energy drift below `1e-9` (measured
`2.2e-15`). A 1200-step run at half the timestep must reduce the frequency
error by at least a factor of three; second-order convergence predicts four,
and the measured ratio is 4.0.

The small-step `theta=0.51` smoke runs and `theta=0.5` energy run are distinct:
`theta=0.51` is dissipative. Neither run establishes stability for large plasma
timesteps. Nonuniform perturbations, clustered particles, magnetic fields,
solver-stress cases, race checking, and backend comparisons should be evaluated
separately before drawing production-port performance conclusions.

## Recorded development validation (2026-10-24)

The complete command above with `--mpi --physics --reference` completed with
exit status zero: 21 recorded case/results entries all passed. Artifacts were
written to `/tmp/fleks-gpu-validation/results.json` and adjacent per-case logs.
The observed largest true electric-solver residual was `3.47e-18`.

| Check | Observed result |
| --- | --- |
| Production CPU one-step Ex | Demo `1.2562184139005094e-4`; production ASCII `1.2562184e-4` |
| Production CPU one-step field energy | Demo `6.278999529185812e-10`; production `6.278999854368504e-10` |
| Nonuniform Gauss residual, correction off/on | `1.1326285671537255e-4` / `9.255066107916087e-15` |
| Cold plasma expected frequency | `3.615110709083961` |
| Measured frequency, dt=0.01 | `3.6147166079451716` (relative error `1.09e-4`) |
| Measured frequency, dt=0.005 | `3.6150122491711394` (relative error `2.72e-5`) |
| theta=0.5 relative total-energy drift, 600 steps | `2.22e-15` |

Halving the timestep reduced frequency error by approximately four, consistent
with second-order temporal convergence for this uniform case. These are
correctness results on the development CUDA backend, not a CPU/GPU speedup
comparison or a production-scale performance claim.
