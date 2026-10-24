# Particle implementation notes

## Production correspondence

- `DemoParticles` derives from the installed AMReX
  `AmrParticleContainer<4,1>` with the single-level vector-geometry constructor.
  The installed headers retain AoS containers without a deprecation attribute.
  Device compilation succeeds. The CUDA compiler reports pre-existing AMReX
  defaulted-function annotation and CUB `Sum` warnings, not AoS deprecation.
- Positions and `rdata(0..2)` velocities, `rdata(3)` signed macro-particle charge,
  and one integer species attribute match the production struct layout.
  Here `qp = species_charge * number_density * cell_volume / particles_per_cell`.
  Species mass and charge independently determine q/m; kinetic energy uses
  `qp * species_mass / species_charge`, which is positive for either sign.
- A quiet lattice is used for perfect-cube particle counts, with deterministic
  low-discrepancy samples otherwise. The seed permutes sample assignments to
  IDs. Global cell/sample/species indices determine IDs independently of boxes
  or ranks. The benchmark-only clustered option compresses positions to the
  middle fifth of the box and redistributes them.
- `initialize(A, drift)` assigns `vx = drift.x + A*cos(2*pi*x)`. It does not
  shift the position time level; the orchestrator calls `drift_positions(dt/2)`
  when adopting FLEKS's initial half-step positions.
- Gather and push are separate kernels with six per-particle scratch values
  stored outside the particle ABI. The fused path gathers and pushes directly.
  The default drift uses the updated velocity for dt, matching production
  `pushStage=0`; average-velocity drift is explicitly selectable for experiments.
  Periodic wrapping precedes AMReX redistribution.
- Moment components are charge density, three charge-weighted velocity moments,
  and six charge-weighted second moments (xx, yy, zz, xy, xz, yz), all divided
  by volume. This ten-component demo omits production's additional particle-count
  moment. These are current/charge moments, not mass moments or thermal pressure.
- `jHat` is deposited in raw charge units; `DemoFields::finish_response()` scales
  it by inverse volume. The matrix is already density-normalized. Directly
  scattering all 64 node pairs avoids production's host-side negative-z copying
  pass while retaining the same 27-neighbor x 3 x 3 storage and alpha coefficients.
- The cell-centered 27-component response uses the same CIC shape gradients,
  absolute macro charge, and inverse volume as `ParticlesMassMatrix.cpp`.
  Position correction uses sign(q)*4*pi*interpolated potential gradient and the
  production 0.1-cell displacement limiter. Periodic domains need no active-region
  rejection or wall guards. Redistribution follows every correction.

## Checks

`run_checks()` uses seven deliberately placed particles at nodes, cell centers,
periodic seams, arbitrary positions, and colliding deposition neighborhoods.
It checks count/weight normalization, ID-preserving multi-box/rank migration,
free streaming, uniform E acceleration, uniform B kinetic energy, rotation, and timestep convergence,
nonuniform CIC gather with an independent cross-product Boris reference,
fused/split equivalence, all ten nodal moments, jHat, all 243 matrix components,
cell charge, all 27 position-response components, and constant-potential position
invariance. Host references use independent nested-loop shape construction and
periodic indexing; no device CIC or Boris helper is reused. Snapshot copies are
explicit check-only host transfers, with synchronization before host buffer reuse.

## Performance experiments

`accumulate_response_buffered()` sorts particles by cell and constructs cell
ranges. Each cell thread sums its particle contributions into local 8-node
current and 64-pair 3x3 buffers, then flushes the buffers with global atomics.
The measured call includes sorting, range allocation/construction, and required
synchronization. This is an initial correctness-first alternative, not a tuned
shared-memory CUDA implementation: local arrays may spill and sparse cells may
lose to direct atomics. Timing must decide the result.

`benchmark_layouts()` copies particles into AoS and seven component SoA scratch
arrays. It measures identical fused gather/push kernels with positions frozen,
so interpolation remains box-local without migration. One warmup is excluded
from both paths; the maximum rank time is reported. Transpose/setup cost is
reported separately, as is the final velocity difference. This compares memory
access layouts only. It does not implement `ParticleContainerPureSoA` or establish
SoA container migration performance; a production conversion needs that separate
interface and lifecycle work.

## Verification results (2026-10-24)

- GPU particle checks passed on one and two MPI ranks with `n_cell=4`,
  `max_grid_size=2`, `particles_per_cell=8`. Both charge signs pass all checks,
  including the buffered matrix reference and gyromotion timestep refinement.
  The normalized buffered matrix error was `8.47e-22`; jHat error was `1.39e-17`.
- Compute Sanitizer memcheck passed with **0 errors** for the complete particle
  checks and a two-step nonuniform coupled implicit run with div(E) correction.
  In the latter run the Gauss maximum fell from `2.32e-4` initially to `1.51e-15`,
  with all 1024 particles retained and zero net charge.
- The installed HPCX MPI probes host-buffer pointer attributes using CUDA APIs
  that deliberately return invalid-value errors. A first sanitizer run reported
  these handled MPI API errors despite all application checks passing. The clean
  memory-check runs used `mpirun --mca opal_cuda_support 0` to disable MPI's CUDA
  pointer probes, consistent with the demo's host-staged MPI path. No application
  memory/API errors were suppressed. The invocation was:

  ```bash
  mpirun --mca opal_cuda_support 0 -np 1 compute-sanitizer \
    --tool memcheck --error-exitcode 99 bin/main_gpu.exe \
    mode=particles n_cell=4 max_grid_size=2 particles_per_cell=8 n_iter=1
  ```

- Racecheck was not used to claim global scatter safety: its shared-memory
  hazard checks do not verify global-memory atomic deposition. Independent
  collision/decomposition checks and memcheck cover the implemented paths.
