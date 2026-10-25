# Implementation progress

Branch: `gpu`. All commits remain local; do not push.

## 2026-10-24 — baseline and interfaces

- Read the revised implementation plan, including the div(E) correction and
  production-reference requirements.
- Preserved the uncommitted field-test baseline in commit `3619693` before
  moving the source with `git mv`.
- Split work into field/GMRES, particles/deposition, and validation tasks using
  the FLEKS skill's implementation-agent workflow.
- Established shared configuration/grid interfaces and began extracting the
  original field checks into `DemoChecks`.
- Validation is pending compilation and GPU runs. A successful CUDA build
  alone will not be recorded as runtime verification.

Further entries record completed work, validation evidence, and limitations as
implementation proceeds. The implementation plan remains the acceptance list.

## 2026-10-24 — stages 0–5 implemented; integration checks

- Added configuration/grid ownership, preserved baseline checks, dual-grid
  field kernels, production-layout AoS particles, CIC gather/ten-moment atomic
  scatter, Boris push, periodic migration, and deterministic initialization.
- Added custom 243-real mass response and jHat assembly, device-resident
  restarted GMRES, theta-stage E/B advance, and the scalar response solve and
  particle-position correction for div(E).
- Production audit corrected two assumptions in the plan: rhoTheta is 0.51,
  not 0.5, and the electric operator includes coefDiff=0.1 by default.
- Initial one-rank GPU tests passed all original field checks, manufactured
  field/solver checks, and both positive/negative species reference checks.
  Uniform and nonuniform coupled runs passed eight timesteps with div(E)
  correction enabled; the latter reached a Gauss residual around 2e-14.
- Independent production-operator comparison exposed transverse averaging in
  div_center_to_center that the first port had missed. Fixed the implementation;
  rerunning the complete composed operator comparison after rebuilding.
- Stage 6 now includes component-SoA scratch vs AoS gather/push, real cell-buffered
  vs atomic mass deposition, custom vs component matrix storage, fused-push
  selection, sorting intervals, clustered particles, and ghost-width timings.
- Full multi-rank, long plasma oscillation, production single-step comparison,
  profiling, and memory-check validation are still pending at this entry.

## 2026-10-24 — full functional regression passed

- Built both executables with nvcc/CUDA sm89 and formatted all new C++ with
  clang-format 20.1.8. Header dependency simulation confirmed that a changed
  DemoConfig.h rebuilds its field/particle/driver consumers. Compilation database
  entries point to src_gpu; production sources remain isolated.
- Full validation command: `python3 src_gpu/validate_demo.py --mpi --physics
  --reference bin/gpu_reference.exe --output-dir /tmp/fleks-gpu-validation`
  with the host-specific MPI environment documented in README. All 21 recorded
  results passed, including both GPUs, ppc 1/8/27, box/rank invariance, and the
  matched production CPU one-step fixture.
- Exact production primitive operators agreed within 7.11e-15 and the complete
  spatial E operator within 1.67e-15, including mixed three-dimensional fields.
- Nonuniform div(E) correction reduced Gauss error from 1.1326e-4 to 9.26e-15.
- Uniform plasma frequency errors were 1.09e-4 at dt=.01 and 2.72e-5 at dt=.005:
  the fourfold improvement demonstrates second-order time convergence.
  Theta=.5 total-energy drift was 2.22e-15 over 600 timesteps.
- Strengthened evolved-B divergence diagnostics with a checkerboard regression.
  Corrected the plan's scatter-cost estimate: full CIC response has 576 tensor
  component atomics/particle, while allocated matrix storage remains 243/node.
- Profiling and compute-sanitizer checks follow; no production speedup is
  inferred from these tiny correctness cases.

## 2026-10-24 — profiling, memory checking, and completion

- Compute-sanitizer memcheck returned zero errors for the comprehensive
  particle mode and for a corrected nonuniform implicit run. HPC-X's handled
  host-pointer CUDA probes were disabled with `--mca opal_cuda_support 0`,
  consistent with host-staged MPI; no application memory errors were suppressed.
  Reproduction details and log paths are in PARTICLE_NOTES.md.
- Ran eight profiling configurations: baseline, fused, sorted, buffered,
  clustered atomic/buffered, smaller boxes, and two GPUs. Every case passed.
  BENCHMARKS.md records the commands, measured costs, and interpretation.
- The basic per-cell buffering experiment is slower, particularly for clustered
  particles; it remains opt-in. Custom response exchange dominates the tiny
  baseline. The SoA comparison is explicitly a scratch kernel experiment, not
  a container/migration port. No production-scale speedup is claimed.
- Full C++ formatting check (clang-format 20.1.8), Python syntax checks,
  compilation database checks, production-library isolation, generated header
  dependency checks, incremental direct build, and `git diff --check` passed.
- Implementation and verification are complete for the single-level scope.
  All changes will be committed to local branch `gpu`; no push is authorized.

## 2026-10-24 — review follow-up: validation strength

- Tightened automated gates to roughly ten times the measured values:
  plasma frequency 6% -> 1e-3, theta=0.5 energy drift 2% -> 1e-9, div(E)
  reduction 5% -> six orders of magnitude, dt refinement must cut the error by
  at least 3x, decomposition absolute floor 2e-10 -> 1e-14, and minimum nodal
  E_x added to the compared quantities.
- Added box/MPI invariance on a nonuniform wave with div(E) correction, and a
  large-dt case (c*dt/dx=8, b0=0.5) that needs up to 28 GMRES iterations, plus
  its two-rank comparison.
- Added oblique-B Boris and mass-matrix references so every alpha entry,
  including the omega omega^T terms, is checked independently.
- Named production rhoTheta once in DemoFieldKernels.h and documented the
  ghost-fill ordering that the E operator relies on.
- Full suite (`--mpi --physics --reference`) passed on two RTX 4090s.
- Still open: the production CPU fixture compares only three global scalars
  for a spatially uniform state without div(E) correction. A nonuniform
  fixture with field dumps and correction enabled would require new
  production FLEKS runs.

## 2026-10-24 — GPU-optimized mass-matrix assembly

- Added `response=cell`: counting-sort bins plus one 64-thread block per cell
  that reduces the 64 node-pair tensors and jHat in registers before 576
  atomics per cell. Validated against the serial oblique-B reference, both
  species, one and two ranks, and the atomic path in benchmark mode and in a
  coupled nonuniform run (`validate_demo.py`).
- 64^3 cells, 100 particles per cell on one RTX 4090: assembly 0.536 s ->
  0.073 s (7.4x); timestep 1.064 s -> 0.603 s, now 3.6x faster than 32 CPU
  cores. See BENCHMARKS.md.
- div(E) solves that reach `div_e_max_iter` now warn instead of failing,
  matching production LinearSolver, which applies the partial potential.

## 2026-10-25 — cell-reduced particle deposits

- Implemented phases 1 and 2 of DEPOSIT_OPTIMIZATION_PLAN.md as
  `deposit=cell`: block-reduced div(E) charge/correction-matrix deposit and
  ten-moment deposit, with CPU fallback to the per-particle versions.
- Validated against the serial references on one and two ranks, per species
  against the atomic deposits in benchmark mode (uniform and clustered), and
  in a coupled nonuniform div(E) run (`cell-deposit-matches-atomic`). Full
  suite passed: 39 checks.
- Reference case timestep 0.599 s -> 0.387 s; moment deposit 7.5x, div(E)
  deposit 2.0x (now FP64-bound). The charge-only deposit stays per particle
  because the cell version measured slower. Phase 3 skipped by its criterion.

## 2026-10-25 — next performance plan implementation

- Added throughput-window timing, warmup/microbenchmark controls, complete-step
  accounting, solver work totals, profiling ranges, and a repeated alternating
  measurement runner with diagnostics, provenance, spread, and arena peaks.
- Implemented spatial pair reuse in both cell matrix deposits, including
  particle-sharing lane groups; added optional folded gradients and bounded
  cell chunks. Removed duplicate operator halo fills and fused normalization.
- Added optional CGS2 with device-local and MPI projection batches, fused
  vector updates, owner masking, and independently checked true residuals.
- Five alternating reference trials: recommended pair reuse plus fused push
  reduced median timestep from 0.41599 to 0.35493 s (14.7%), reaching
  147.7M particles/s. Folded gradients and batched GMRES were slower here;
  keep them optional. Chunking gave a small gain on the clustered comparison.
- Final GPU suite: 45 cases passed. CPU suite and refreshed solver/MPI checks
  passed, as did 12 timing edge cases. CUDA memcheck and racecheck were clean.
- Recorded an Nsight Systems trace and measured outcomes in BENCHMARKS.md.
  Migration deferral, preconditioning, adaptive selection, matrix-storage
  replacement, and cross-device scaling remain conditional follow-up work.

## 2026-10-25 — div(E) solver optimization

- Implemented DIVE_SOLVER_OPTIMIZATION_PLAN.md phases 0-2 and 4: solver
  counters and section timers (`SOLVER` lines, `solver_sections`), a
  synchronization-free batched CGS2 (one host wait per Arnoldi step on one
  rank, basis stored once, device-resident projections, fused launches,
  norm from the second-pass reduction), a zero-initial-guess shortcut and
  fused right-hand side, and opt-in `div_e_method=jacobi_gmres|pcg`.
  Phase 3 was not needed on one GPU.
- Scalar solves 33.8 -> 13.8 ms per step (8 boxes) and 35.6 -> 15.3 ms
  (64 boxes) with identical iterations and warnings; complete step -7%.
  CGS2 is now the default; `gmres_batched=0` keeps MGS2.
- Validated with unit checks, coupled runs on one and two ranks, the CPU
  build, compute-sanitizer memcheck and racecheck, and the full suite.

## 2026-10-25 — single-precision cell matrix sums

- Added opt-in `cell_precision=single`: the response and div(E) matrix
  kernels stage pair inputs and sum each cell in float, then add into the
  double-precision matrices. Charge, jHat, and moments stay double because
  they cancel between species.
- Reference case: response assembly 2.0x, correction deposit 1.36x, complete
  step 15.9% shorter (0.332 -> 0.279 s, throughput mode on GPU 1) with
  identical solver work.
- theta=0.5 energy drift rises from 2.2e-15 to 8.6e-9; the default stays
  double. Validation adds single-precision serial references, precision-aware
  benchmark comparisons, a single-precision energy/frequency gate (1e-7), and
  a single-precision nonuniform div(E) reduction check.
