# div(E) solver optimization plan

Status: phases 0, 1, 2, and 4 implemented; phase 3 not needed on one GPU.
Batched CGS2 is now the default (`gmres_batched=1`). See "Outcome" at the
end for measurements and differences from the expectations below.

## Goal

The scalar div(E) solve is the one timestep phase where one RTX 4090 loses
to the CPU (`CPU_GPU_COMPARISON.md`):

| Layout | CPU (ms/step) | GPU (ms/step) | CPU/GPU |
| --- | ---: | ---: | ---: |
| 64 boxes of 16^3, 64 CPU cores | 28.16 | 41.57 | 0.68x |
| 64 boxes of 16^3, 32 CPU cores | 34.63 | 39.03 | 0.89x |
| 8 boxes of 32^3, 8 CPU cores | 157.81 | 35.60 | 4.43x |

Make the GPU scalar solve at least 3x faster than now, under 12 ms per step on
the 8-box reference layout and under 15 ms on the 64-box layout, so that it
beats all 64 CPU cores. Do not change the production-matching solver result:
the same GMRES iterates up to roundoff, the same `div_e_tolerance=0.01`,
`div_e_max_iter=20`, three correction passes, and the same correction-limit
warnings. Algorithm changes that alter the partial solution at the iteration
limit are kept separate (phase 4) and opt-in.

Reference case: the `CPU_GPU_COMPARISON.md` workload (64^3 cells, 100
particles per cell per species, Langmuir wave, `b0=0.5`, `dt=0.1`,
`solver_tolerance=1e-6`) with the recommended flags `response=cell
deposit=cell pair_symmetry=1 fused_push=1`, on both the 8-box and the 64-box
layouts. Over five timed steps it performs 287 scalar iterations (about 19 per
solve, 15 solves) and 70 electric iterations.

## Diagnosis

Each scalar iteration costs about 0.68 ms (39 ms / 57.4 iterations per step).
The work per iteration is small: a vector is 2 MB (262,144 cells), the
27-point operator reads about 57 MB of coefficients, and orthogonalizing
against an average of 10.5 basis vectors twice reads about 90 MB. At roughly
1 TB/s that is about 0.15 ms. The solve is latency-bound, not bandwidth-bound.

The latency comes from host synchronization. In this AMReX version:

- `MultiFab::Dot` reduces on the device, copies the scalar to the host, and
  then performs an MPI all-reduce: one blocking round trip each.
- Fused whole-MultiFab vector operations (`Saxpy`, `LinComb`, `Copy`, ...)
  call `Gpu::streamSynchronize()` after their launch unless executed inside
  `Gpu::NoSyncRegion` (`AMReX_FabArray.H`, `FabArray<FAB>::Saxpy`).
- Every `MFIter` loop synchronizes the device when it finishes unless it is
  built with `MFItInfo().DisableDeviceSync()` (`AMReX_MFIter.cpp`). That
  includes `apply_div_e_operator` and the batched-solver helpers.
- Whole-MultiFab fusing applies only when a rank owns more than one box;
  otherwise operations launch one kernel per box.

The default twice-modified-Gram-Schmidt (MGS2) GMRES therefore performs, at
Arnoldi step `j`, 2(j+1) dot products and 2(j+1) Saxpys, each followed by a
host sync, plus the operator, the norm, and the normalization: about 46
blocking round trips per iteration on average. At about 15 us each, that is
the measured 0.68 ms.

On the CPU, the same synchronizations are a few-microsecond `MPI_Allreduce`
among ranks that each own 4,096 cells, which is why the CPU wins this phase.

The existing opt-in `gmres_batched=1` (twice classical Gram-Schmidt, CGS2)
already cuts the dot-product syncs to one per pass, but it measured slower for
the whole timestep (0.424 s versus 0.416 s). It has never been timed per
solver. Visible overheads in `DemoLinearSolver.cpp`:

- Each new basis vector is copied into `packedBasis` (`MultiFab::Copy`, which
  syncs), so the basis is stored twice.
- `add_basis` allocates a `Gpu::AsyncArray` and copies the coefficients from
  host to device for every box on every call, then synchronizes at the end of
  its `MFIter` loop. The coefficients are already on the device in
  `deviceProjection`.
- `project_basis` launches one kernel per box rather than one per vector set.
- The normalization `LinComb`, the norm `Dot`, and the operator still sync.

## Phase 0: per-solver baseline and instrumentation

1. Add stage-timing sub-phases inside `DemoLinearSolver::solve`: operator,
   orthogonalization, norm and normalization, solution update, and true
   residual. Report them separately for the scalar and electric solvers so a
   change to the shared solver can be judged on each.
2. Count blocking synchronizations per solve: wrap the solver's dot products
   and explicit syncs with a counter, and cross-check with Nsight Systems
   (`cudaStreamSynchronize` and device-to-host `cudaMemcpyAsync` calls inside
   an NVTX range around each solve).
3. Measure MGS2 and the existing `gmres_batched=1` per solver on both
   layouts, five alternating trials with `measure_performance.py`, on a quiet
   host with recorded load, and compare against the bandwidth bound above.

Deliverable: a table per solver of milliseconds per iteration, syncs per
iteration, and kernel launches per iteration for MGS2 and CGS2. This decides
whether phase 1 is enough or phase 3 is needed.

## Phase 1: remove synchronizations from batched CGS2

Keep the CGS2 algorithm unchanged; remove its implementation overheads.

1. Store the basis once. Allocate `packedBasis` with one ghost layer and make
   each `basis[j]` an AMReX alias of its components (`amrex::make_alias`), so
   the operator reads its input directly from the packed basis and no copy is
   made per iteration.
2. Keep projections on the device for the vector update. `add_basis` reads
   the negated coefficients from `deviceProjection` (negated by a one-thread
   kernel or by the update kernel itself). Only the Hessenberg column goes to
   the host, through the single pinned copy that already exists. Remove the
   per-box `AsyncArray`.
3. Launch once per operation for all boxes: rewrite `project_basis`,
   `add_basis`, the normalization, and `apply_div_e_operator` with
   whole-MultiFab `ParallelFor(mf, ...)` / `const_arrays()` launches, and
   wrap them in `Gpu::NoSyncRegion` or `MFItInfo().DisableDeviceSync()` so
   that ordering relies on the single stream, not on host waits.
4. Fuse the norm into the second CGS pass's reduction: compute `w . v_i` for
   all `i` and `w . w` in one reduction, then form the post-projection norm on
   the host as `sqrt(w.w - sum p_i^2)`. If that difference falls below a
   safety fraction of `w.w` (cancellation), recompute the norm with an
   explicit dot. Record how often the fallback triggers.

Expected syncs per iteration: two (one per CGS pass, the second also giving
the norm), down from about 46. The device work stays at about 0.15 ms, so the
target is about 0.25 ms per iteration: about 14 ms per step for the scalar
solve.

Acceptance:

- Scalar and electric solves take the same number of iterations and produce
  the same correction-limit warnings as MGS2 on the reference case; coupled
  diagnostics pass `compare_records` against the MGS2 run.
- Existing solver checks pass: manufactured potential, true residual, and
  the varying nonsymmetric batched-solver test in `DemoFields::run_checks`.
- The CPU build (host path) and two-rank runs pass; compute-sanitizer
  memcheck and racecheck are clean for the new kernels.
- If CGS2 is faster than MGS2 for both solvers on both layouts, make it the
  default and keep MGS2 as the reference option (`gmres_batched=0`).

## Phase 2: per-solve overhead

Each of the 15 scalar solves per measured run also pays fixed costs outside
the Arnoldi loop:

1. `solve_div_e` starts from `phi = 0`, so the initial true residual is the
   right-hand side. Skip the operator application and its dot when the
   caller states that the initial guess is zero. Production starts from zero
   as well.
2. Fuse the right-hand-side construction (divergence, `LinComb`, `Saxpy`,
   and mean subtraction) into one kernel plus one reduction, and fuse the
   final potential mean subtraction and `FillBoundary`.
3. Measure the `FillBoundary` in each operator application. On one rank it is
   a local periodic copy; keep it, but issue it without a host wait. With
   several ranks, overlap the halo exchange with the interior stencil
   (`FillBoundary_nowait` / `FillBoundary_finish`) only if phase 0 shows the
   exchange as a significant share.

Expected: about 2-4 ms per step. The true residual at the end of each
restart cycle stays, because validation reports it.

## Phase 3: one-synchronization orthogonalization (only if needed)

If phases 1 and 2 leave the scalar solve above target, replace CGS2 with a
low-synchronization variant that needs one global reduction per iteration
by lagging the normalization and fusing the reorthogonalization
(Swirydowicz, Langou, Ananthan, Yang, and Thomas, "Low synchronization
Gram-Schmidt and generalized minimal residual algorithms", Numer. Linear
Algebra Appl., 2021; and the DCGS2 variant of Bielich et al., Parallel
Computing, 2022). These are equivalent to CGS2 in exact arithmetic and have
comparable stability in the cited analyses, but they change rounding, so
they need the same acceptance tests as phase 1 plus an orthogonality check
(`||I - V^T V||`) on the varying nonsymmetric test operator.

Expected: one sync per iteration, about 0.2 ms per iteration.

## Phase 4: fewer iterations (changes results; production decision)

The solves average about 19 iterations, and six of the fifteen stop at the
20-iteration limit with the 1% tolerance unmet, so the cost is set mostly by
the limit rather than by convergence.
Options that reduce iterations change the partial potential that is applied,
and so change the particle correction:

- Jacobi (diagonal) preconditioning from the center entry of `centerMM`.
  It is cheap and may let more solves reach 1% within the limit.
- Conjugate gradients. The operator `4 pi^2 M` with `M` built from
  `|q| grad w_a . grad w_c` is symmetric positive semidefinite, with the
  constants as its null space, which the mean subtraction already removes.
  CG needs two reductions per iteration (one with the Chronopoulos-Gear
  formulation) independent of the iteration count, and no basis storage.

Evaluate them by Gauss-residual reduction per millisecond on the reference
case and on the nonuniform `validate_demo.py` cases, opt-in only. Adopting
either changes FLEKS physics output and must be decided for production
FLEKS, not just for the demo.

## Shared effects on the electric solve

The scalar and electric solves share `DemoLinearSolver`. The electric solve
(14 iterations per step, 47 ms) has a heavier operator (curl-curl, the
coefDiff term, and the 243-coefficient response, about 0.6 ms per
application), so it is less latency-bound and should gain less. Phases 1-3
must not slow it down; phase 0's per-solver table is the check.

## Validation summary

For every phase, in addition to the acceptance lists above:

- Full `validate_demo.py --mpi --physics --reference` suite.
- Reference case with identical electric and scalar iteration counts and
  warnings (phases 1-3).
- CPU build with `mode=all` through the host code path.
- Report the remaining GPU-to-CPU64 ratio for the scalar solve using the
  `CPU_GPU_COMPARISON.md` protocol (pinned cores, alternating trials,
  recorded host load).

## Expected outcome

| Phase | Syncs per iteration | Scalar solve, 8 boxes (ms/step) |
| --- | ---: | ---: |
| Now (MGS2) | about 46 | 35.6 |
| 1: synchronization-free CGS2 | 2 | about 14 |
| 2: per-solve overhead | 2 | about 11 |
| 3: one-sync orthogonalization (if needed) | 1 | about 9 |

Phase 4 is evaluated separately, because it changes the correction rather
than only its cost.

## Outcome (2026-10-25)

Measured on one RTX 4090 with the reference case and the recommended flags,
one MPI rank pinned to an idle core, three alternating trials per
configuration with `measure_performance.py` (`pairs-fused` = MGS2,
`combined` = the new CGS2). The host was shared; series disturbed by other
users' jobs were rerun, and the ranges below are from undisturbed series.

Stage mode (`benchmark_timing=stages`, no solver section timers):

| Stage (ms/step) | 8 boxes MGS2 | 8 boxes CGS2 | 64 boxes MGS2 | 64 boxes CGS2 |
| --- | ---: | ---: | ---: | ---: |
| Scalar div(E) solves | 33.80 | **13.79** | 35.58 | **15.30** |
| Electric solve | 46.28 | 42.78 | 46.78 | 43.53 |
| Stage-timed step | 343.6 | 319.8 | 396.2 | 372.5 |

Throughput mode (complete step, diagnostics included):

| Layout | MGS2 (s/step) | CGS2 (s/step) | Change |
| --- | ---: | ---: | ---: |
| 8 boxes of 32^3 | 0.3449 (0.3440-0.3452) | 0.3203 (0.3201-0.3204) | -7.1% |
| 64 boxes of 16^3 | 0.3999 (0.3951-0.4060) | 0.3728 (0.3722-0.3735) | -6.8% |

The 8-box CGS2 throughput is 163.7M particles per second. Every trial had
70 electric and 287 scalar iterations and 6 correction-limit warnings,
identical to MGS2, and no norm fallback triggered.

The scalar solve is now 2.0x faster than all 64 CPU cores (28.16 ms in
`CPU_GPU_COMPARISON.md`), on both layouts. The 3x target was reached for
neither layout (2.45x and 2.33x); the 64-box target of 15 ms is met within
measurement precision and the 8-box target of 12 ms is missed by 1.8 ms.

Phase 0 findings (solver counters, `solver_sections=1` for the section times):

- Host waits per scalar iteration fell from 43.8 to 1.10, and per electric
  iteration from 33.4 to 1.21. Nsight Systems confirms the count: over three
  steps the new solver issues 8,722 fewer `cudaStreamSynchronize` calls,
  against 8,696 predicted by the counters, and 20,316 -> 8,640 kernel
  launches.
- With CGS2, one scalar orthogonalization costs about 0.15 ms and one
  operator application about 0.06 ms, both within about 1.3x of the
  memory-bandwidth bound estimated in the diagnosis. The remaining gap to
  the 8-box target is bandwidth, not waiting.

Differences from the plan:

- Phase 1 already reaches one host wait per Arnoldi step on one rank: both
  CGS2 passes, their vector updates, and the norm run on the device, and the
  host reads the Hessenberg column once. Phase 3's purpose, one global
  reduction per step, is therefore met; low-synchronization variants would
  only help multi-rank runs, where each pass still needs an MPI reduction.
  Phase 3 was not implemented.
- The norm comes from the second pass's reduction as
  `sqrt(w.w - sum p_i^2)`, with an explicit dot product when more than 99%
  of `w.w` would cancel. The fallback never triggered in the reference case.
- Phase 2: the zero-initial-guess shortcut and the fused right-hand side are
  in place; their saving is included in the numbers above and was not
  measured separately.
- `MFIter` loops and AMReX vector operations inside `Gpu::NoSyncRegion` still
  spread boxes over up to four streams without a closing synchronization.
  The new code therefore uses only single-launch whole-MultiFab kernels in
  those regions; the electric operator keeps its own synchronizing loops.

Phase 4 (opt-in, `div_e_method`), same case, 8 boxes, CGS2:

| Method | Scalar solves (ms/step) | Iterations | Warnings | Post-correction Gauss residual, steps 1-6 |
| --- | ---: | ---: | ---: | --- |
| `gmres` (production) | 13.8 | 287 | 6 | 1.5e-8 to 4.1e-8 |
| `jacobi_gmres` | 14.1 | 283 | 7 | 1.7e-8 to 5.2e-8 |
| `pcg` | 4.9 | 271 | 7 | 6.4e-9 to 8.7e-8 |

Jacobi preconditioning does not help: in a near-uniform plasma the
diagonal of the correction operator is nearly constant. Jacobi PCG costs
2.8x less per step, because its iteration cost does not grow with a Krylov
basis, at a similar iteration count; its Gauss residual is better on some
steps and worse on others. Adopting PCG changes the applied correction and
needs evaluation in production FLEKS on longer, more strongly nonuniform
runs; it remains opt-in.

Validation: solver unit checks and coupled runs for both GMRES variants and
all three methods on one and two ranks; compute-sanitizer memcheck (0
errors) and racecheck (0 hazards) for `gmres` and `pcg`; CPU build on four
ranks with every method; `validate_demo.py` now pins the reference runs to
MGS2 so its CGS2 comparisons stay meaningful, and checks that `jacobi_gmres`
and `pcg` achieve the same nonuniform Gauss-residual reduction as GMRES.
