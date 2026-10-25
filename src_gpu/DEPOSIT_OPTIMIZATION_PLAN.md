# Particle deposit optimization plan

Status: phases 1 and 2 implemented (`deposit=cell`); phase 3 skipped by its
own criterion. See "Outcome" at the end for measured results and where they
differ from the expectations below.
Development concluded on 2026-10-25.

## Goal

After the cell-reduced mass-matrix assembly (`response=cell`, commit
`767a754`), the remaining particle-to-grid deposits dominate the GPU timestep.
Apply the same per-cell reduction to them and cut the reference timestep from
0.60 s to about 0.3 s without changing any result beyond floating-point
summation order.

Reference case for every measurement in this plan: one RTX 4090, `n_cell=64`,
`particles_per_cell=100` per species (52.4M particles), `max_grid_size=32`,
`initial_condition=langmuir perturbation=0.01 b0=0.5 dt=0.1 correct_div_e=1
solver_tolerance=1e-6 response=cell`, `mode=benchmark n_step=6 n_iter=3`.

Current GPU cost per timestep (benchmark mode, stage-synchronized):

| Stage | Seconds | Share | Atomics per particle |
| --- | ---: | ---: | ---: |
| div(E) charge + matrix deposit, 3 passes | 0.220 | 36% | 8 + 64 per pass |
| Moment deposit, 10 moments | 0.121 | 20% | 80 |
| Final charge deposit | 0.013 | 2% | 8 |
| Mass-matrix assembly (`response=cell`) | 0.073 | 12% | 576 per cell |
| Everything else | 0.176 | 29% | |
| **Timestep** | **0.603** | | |

## Diagnosis

All three deposits use one thread per particle and add each contribution
directly to global memory. Particles are stored in cell order, so the 32
lanes of a warp usually belong to the same cell and add to the same
addresses; the hardware serializes those atomics. The charge-only final
deposit (8 atomics per particle) costs 13 ms per call, while each matrix pass
(72 atomics per particle) costs 73 ms. The cost tracks the number of
contended atomics, not arithmetic.

The mass-matrix work showed that reducing within a cell first, and then
issuing one atomic per output value per cell, removes this contention:
assembly went from 0.536 s to 0.073 s.

## Shared design rules

These apply to every new kernel in this plan.

- Keep the existing per-particle deposits as the default and as the CPU
  fallback (`#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)`), as
  `accumulate_response_cell` does.
- Select the new path with one option, `deposit=atomic|cell`, covering the
  div(E) charge/matrix deposit, the final charge deposit, and the moment
  deposit. Keep `response` separate for the mass matrix.
- Bin with `amrex::DenseBins` (index permutation; particles do not move).
  Bin by the clamped floor cell of the deposit's own shape function and
  evaluate the CIC fractions relative to that cell inside the block, exactly
  as `accumulate_response_cell` does, so particles on box faces get exact
  weights.
- Launch one block per bin. Empty bins return before any `__syncthreads()`.
  Stage per-particle quantities in shared memory in chunks equal to the block
  size, accumulate in registers, and issue global atomics only after the cell
  is finished.
- Skip particles with `id() <= 0` by zeroing their staged weights.
- The caller-side contract does not change: kernels add to existing storage;
  `DemoPic` zeroes, calls `SumBoundary` once after all species, and fills
  ghosts.

## Phase 1: div(E) charge and correction matrix (largest gain)

Add `DemoParticles::deposit_center_charge_cell(MultiFab& rho, MultiFab*
centerMatrix)`, matching `deposit_center_charge`.

Binning. The center shape uses `shift = 0.5`, so its base index is
`floor(u - 0.5)`, which ranges over `[small - 1, big]` in each direction for
particles in the box. Bin over `Box(small - 1, big)`, i.e. `(n + 1)^3` bins
for an `n^3` box; the low layer writes into ghost cells and `SumBoundary`
handles it, as today.

Kernel. One 64-thread block per center cell. Thread `t` owns pair
`(a, c) = (t / 8, t % 8)` of the eight surrounding cell centers and accumulates
one scalar, `|q| / V * dot(grad w_a, grad w_c)`, into a register. Threads
`t < 8` also accumulate `q / V * w_t` for `rho`. Stage per particle in shared
memory: 8 weights and 8 x 3 weight gradients (32 doubles; 16 KB for a
64-particle chunk), plus `q`. Write with the existing layout
`offset = (ii - i + 1) * 9 + (jj - j + 1) * 3 + (kk - k + 1)`. That is 64 + 8
atomics per cell instead of per particle.

When `centerMatrix` is null (the final charge deposit) launch the same kernel
with only the `rho` threads doing work, or a smaller 8-thread-per-cell variant;
measure which is faster.

Optional refinement after the measurement: the matrix contribution is
symmetric in `a` and `c`, so 36 unique pairs suffice. It halves arithmetic,
not atomics; implement it only if the kernel is arithmetic-bound.

Expected: the three passes fall from 0.220 s to roughly 0.03-0.04 s, of which
about 8 ms per pass is the counting sort.

## Phase 2: moment deposit

Add `DemoParticles::deposit_cell(MultiFab& moments)`, matching `deposit`.

Production calls `sum_moments(true)` every PIC step and uses the moments to
update `dt`, so the deposit stays every step; skipping it would make the demo
cheaper than production for the wrong reason. The demo uses the moments only
for diagnostics, so record that difference in `README.md` without exploiting
it.

Kernel. Nodal shape, the same binning as `accumulate_response_cell`. One
block of 96 threads per cell; thread `t < 80` owns node `t / 10` and moment
`t % 10` and accumulates `w_node * value_moment` over the cell's particles.
Stage 8 weights and the 10 moment values per particle (18 doubles). That is 80
atomics per cell instead of per particle.

Document in the code, as the plan required earlier, that all ten moments are
charge-weighted (`q`, `q u`, `q u u`, ...).

Expected: 0.121 s down to roughly 0.015 s.

## Phase 3: share bins between deposits (optional, measure first)

Every bin build costs about 0.5 ms per box and species, about 8 ms per build
for the reference case. With phases 1 and 2 a step builds bins six times:
mass matrix, three div(E) passes, final charge, moments. Particle positions
are identical for the final charge deposit, the moment deposit, and the next
step's mass-matrix assembly.

Option: bin once by half-cell index, `bin = 8 * cell + subcell`, where
`subcell` records which half of the cell the particle occupies in each
direction. A nodal cell is then the contiguous range
`[offsets[8 * cell], offsets[8 * cell + 8])`, and a center cell is the union of
eight sub-bins from neighboring nodal cells. Cache the bins with a validity
flag cleared on every position change (`drift_positions`, `push`,
`correct_positions`, `sort_by_cell`, `set_particles`, `initialize`), as
`gatheredValid` does for gathered fields.

This saves about three builds (about 24 ms) per step. Implement it only if
phases 1 and 2 leave binning as a top-three cost.

## Alternative considered: warp-aggregated atomics

Reduce within a warp with `__match_any_sync` and shuffles, then let one lane
per distinct address issue the atomic. It needs no binning and fewer lines of
code, but its benefit depends on the particle order (which degrades after
migration unless particles are sorted), and it is CUDA-specific. Keep it as a
fallback experiment if binning turns out to dominate after phase 3.

## Validation

Each phase must pass the following before it is measured.

- Serial references in `DemoParticles::run_checks`: extend "cell charge CIC
  serial reference", "divE position-response matrix serial reference", and
  "atomic ten-moment CIC scatter serial reference" to the new kernels, on one
  and two ranks, including particles at periodic corners and box faces.
- Benchmark mode: compare each new deposit entry by entry with the atomic
  deposit after `SumBoundary` (tolerance 1e-10 absolute, as for the mass
  matrix), for uniform and `clustered=1` particle distributions.
- `validate_demo.py`: add a coupled nonuniform div(E) run with `deposit=cell`
  and require `compare_records` agreement with the atomic run, as
  `cell-response-matches-atomic` does now.
- CPU build (scratchpad CPU AMReX) passes `mode=all` with `deposit=cell`
  through the fallback.
- Full suite `--mpi --physics --reference` passes.

## Measurement and reporting

- Rerun the reference case with `deposit=atomic` and `deposit=cell`; report
  every stage, the timestep, and particles per second, plus the 32-CPU-core
  baseline (2.163 s) for context.
- Use Nsight Systems to split counting-sort time from kernel time. Nsight
  Compute is unavailable on this host (`RmProfilingAdminOnly=1`), so state any
  arithmetic-versus-atomic conclusion as an estimate.
- Record results in `BENCHMARKS.md` and a dated entry in `PROGRESS.md`.

## Expected outcome

| Stage | Now (s) | Target (s) |
| --- | ---: | ---: |
| div(E) charge + matrix, 3 passes | 0.220 | 0.03-0.04 |
| Moment deposit | 0.121 | about 0.015 |
| Final charge deposit | 0.013 | about 0.01 |
| Mass matrix | 0.073 | 0.073 |
| Everything else | 0.176 | 0.176 |
| **Timestep** | **0.603** | **about 0.31** |

At that point the remaining costs are spread over the mass matrix (pair
symmetry), the two GMRES solvers (latency-bound global reductions; batch the
Gram-Schmidt dot products), and particle migration after each correction
pass. Those belong in a later plan.

## Outcome (2026-10-25)

Measured on the reference case with the host quiet (load average about 2),
the process pinned to one idle core, `response=cell`:

| Stage | `deposit=atomic` (s) | `deposit=cell` (s) | Speedup | Plan target (s) |
| --- | ---: | ---: | ---: | ---: |
| div(E) charge + matrix, 3 passes | 0.219 | 0.111 | 2.0x | 0.03-0.04 |
| Moment deposit | 0.122 | 0.016 | 7.5x | about 0.015 |
| Final charge deposit | 0.013 | 0.013 | 1.0x | about 0.01 |
| **Timestep** | **0.599** | **0.387** | **1.55x** | **about 0.31** |

Throughput rose from 87.5M to 135.5M particles per second (32 CPU cores:
24.2M). Coupled diagnostics after six steps match the atomic run to about
1e-15 relative, and the per-species deposit fields agree within 4e-15
(1e-14 for clustered particles).

Differences from the plan:

- Phase 1 reached 2.0x, not about 6x. Nsight Systems puts the counting sort
  at 0.22 ms and the kernel at 2.1 ms per box and species. Once the atomics
  are reduced per cell, the per-pair dot products are limited by FP64
  throughput on this consumer GPU (1/64 of FP32); a throughput model gives
  about 1.6 ms. The plan assumed the cost would stay atomic-bound.
- The charge-only cell kernel was slower than the per-particle deposit
  (16.4 ms versus 13.0 ms), because 8 atomics per particle cost less than
  binning plus staging. `deposit_center_charge_cell` therefore routes the
  charge-only call to the per-particle deposit.
- Clustered particles (`clustered=1`, about 1,250 per occupied cell) gain
  nothing from the cell deposits (matrix 0.033 s versus 0.030 s): only about
  2,100 cells are occupied, too few blocks to fill the GPU, each with a long
  serial loop. The mass-matrix kernel still wins there because its
  per-particle atomic cost is much larger.
- The comparison against the atomic deposit is done per species. Summed over
  a neutral plasma, the charge-weighted fields nearly cancel, so roundoff
  becomes a large relative error (7e-12 clustered) that says nothing about
  correctness.
- Phase 3 was skipped: binning totals about 18 ms per step (under 5%), so it
  is not a top-three cost.

Remaining options for the div(E) matrix pass, in order of expected payoff:
use the pair symmetry with particle-split thread groups (36 unique pairs,
up to about 40% less pair arithmetic), fold `|q|` into the staged gradients
(one FP64 instruction less per pair and particle), and split oversized cells
across several blocks for clustered distributions. On a data-center GPU with
full-rate FP64 the arithmetic limit largely disappears.
