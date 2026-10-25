# div(E) deposition optimization plan

Status: phase 1 implemented (`center_kernel=moments`, now the default);
phases 2 and 3 not started. See "Outcome" at the end.

## Goal

The div(E) charge and correction-matrix deposit (`center_deposit`) is the
largest phase of the GPU timestep: 95.8 ms of 318 ms per step in double
precision (30%) and 70.5 ms of 268 ms with `cell_precision=single` (26%),
over three correction passes. Cut it by about 2.5x, to about 35-40 ms per
step, **in double precision** and without changing results beyond
floating-point summation order.

Reference case: the `GPU_DEVELOPMENT_SUMMARY.md` workload (64^3 cells, 100
particles per cell per species, 8 boxes of 32^3, recommended flags
`response=cell deposit=cell pair_symmetry=1 fused_push=1`), one RTX 4090.

## Where the time goes (measured)

Nsight Systems trace of one correction pass, double precision:

| Part of one pass | Time | Share |
| --- | ---: | ---: |
| Cell kernel, 2 species x 8 boxes (1.8 ms per launch) | 28.6 ms | 89% |
| Counting sort (`DenseBins`), 16 builds | 4.0 ms | 12% |
| Clearing and boundary exchange of rho and the 27-entry matrix | 0.15 ms | 0.5% |
| **Pass total (x3 per step)** | **32 ms** | |

The kernel is close to the RTX 4090's FP64 limit (1/64 of its FP32 rate).
Per cell and species it executes about 780 FP64 warp instructions:

- staging, about 70 FP64 operations per particle: 8 CIC weights and 24
  weight gradients, stored in shared memory;
- the pair loop, about 5 FP64 instructions per particle for each of the 36
  unique node pairs (gradient dot product, charge factor, accumulate).

At the 4090's FP64 rate that predicts about 1.3 ms per launch, against 1.8 ms
measured. The kernel uses 120 registers and 17 KB of shared memory per
64-thread block. Only a reduction in arithmetic will make it materially
faster.

## Key idea: separable moments

CIC weights factor by dimension, w_a = X_a(f_x) Y_a(f_y) Z_a(f_z) with
X ∈ {1-f, f}, and each gradient component replaces one factor by
±1/dx. Every matrix entry therefore becomes

    M(a,c) = (1/V) sum_d (1/dx_d^2) s_a^d s_c^d T_d[i(a,c,e1)][j(a,c,e2)]

where s = ±1 is the node's side in direction d, e1 and e2 are the other two
directions, and the per-cell sums are

    T_d[i][j] = sum over particles of |q| P_i(f_e1) P_j(f_e2),
    P_0 = (1-f)^2,  P_1 = f(1-f),  P_2 = f^2.

That is **27 sums per cell** instead of 36 pair accumulators, and each costs
**one product per particle** instead of a three-term gradient dot product.
The charge density needs 8 more sums, q X Y Z. No gradients or weights are
staged.

Checked numerically (random particles, unequal dx): the 27-sum form matches
the current per-particle dot products to 3e-16 relative, and keeps the zero
row sums (constants in the null space) that the div(E) operator relies on.
All 27 sums add non-negative terms, so the accuracy is the same as today's.

Per particle the new work is about 60 FP64 operations in total: the
fractions, 9 basis values, 9 charge-scaled values, 27 accumulations for the
matrix and 8 for the charge. That is about 2 FP64 warp instructions per
particle with one thread per particle, against about 7.8 now.

## Phase 1: moment kernel (main gain)

Add a moment-based `deposit_center_charge_cell` kernel, keeping the current
pair kernel selectable for comparison (`center_kernel=moments|pairs`,
default `moments` once validated).

- **Binning:** unchanged: the same counting sort over half-cell-shifted base
  cells, including the low ghost layer.
- **Mapping:** a group of G lanes per base cell (G = 8 or 16, tuned; 128-thread
  blocks hold 128/G cells). Lanes stride through the cell's particles; each
  keeps 35 double accumulators (27 matrix sums + 8 charge sums) in registers.
  No shared-memory staging and no barriers in the particle loop.
- **Reduction:** shuffle within the G lanes (log2 G steps), with the 35 values
  spread across lanes so each lane finishes a few of them.
- **Output:** the lanes convert the 27 sums to the 64 matrix entries with a
  constant (direction, i, j, sign) table and issue 64 + 8 atomics per cell,
  as today. Same caller contract: add to existing storage; `SumBoundary`
  once after both species.
- **Dense cells:** keep `cell_chunk_size` working through the existing
  `cell_task_range`; each task adds its partial sums.
- **Precision:** double only. With about 4x less FP64 work,
  `cell_precision=single` should no longer pay off for this kernel; measure,
  and if it does not, restrict that option to the response kernel and say so
  in README.

Expected: the kernel drops from 1.8 ms toward its memory floor (reading
about 210 MB of particle data per launch, roughly 0.2-0.3 ms) at an
estimated 0.4-0.6 ms per launch. Per pass: about 8 ms of kernel plus 4 ms
of sorting, about 12 ms instead of 32 ms.

## Phase 2: deterministic gather output (optional, measure first)

After phase 1, the 72 atomics per cell (about 40M per pass) may become
visible. If removing them (a test build with output atomics disabled) saves
more than about 10% of the pass:

- write the 35 sums per base cell into a cell-based array with no atomics,
  accumulating both species into the same array (each base cell has one
  owner block per species);
- run `SumBoundary` on that array (35 components, replacing today's 28);
- one gather kernel builds each cell's 27 matrix entries and its charge from
  the up to 8 base cells that share it.

That makes the deposit bitwise reproducible from run to run, which atomics
cannot provide.

## Phase 3: sorting (optional)

The counting-sort count kernel already runs at about memory bandwidth
(0.21 ms to read 3.3M particles of 64 bytes each), so its contention is not
the issue. After phase 1 it may be a third of the remaining pass. Options,
only if that is the case:

- have the sort's permutation pass also write a compact, bin-ordered copy of
  the four values the moment kernel needs (three fractions and q), so the
  kernel reads coalesced data instead of gathering 64-byte particles;
- longer term, compute the next pass's bin index inside
  `correct_positions`, which already rewrites every position.

## Not in scope

- The mass-matrix kernel. Its weights are separable too, but it needs 243
  sums per cell (27 weight patterns x 9 tensor components) against 324
  pair-symmetric accumulators today: a much smaller gain, and too many
  registers for the one-thread-per-particle mapping.
- The moment and charge-only deposits: they would need their own sort,
  which costs more than they could save.
- Position correction and particle migration (50 ms per step): a separate
  phase.

## Validation

- Serial references in `DemoParticles::run_checks` for the new kernel
  (charge and divE matrix, at the existing double-precision tolerance),
  on one and two ranks, with periodic corners and box faces.
- Benchmark mode: per-species entry-by-entry comparison with the atomic
  deposit at the existing 1e-12 relative tolerance, uniform and
  `clustered=1`, with and without `cell_chunk_size`.
- Coupled reference run: identical electric and div(E) iteration counts and
  warnings; `validate_demo.py` comparison of a nonuniform div(E) run against
  the atomic path; the full `--mpi --physics --reference` suite.
- CPU build unchanged (per-particle fallback); compute-sanitizer memcheck
  and racecheck on the new kernel.

## Measurement

- Nsight Systems per-pass breakdown (kernel, sorting, exchange) before and
  after.
- Stage and throughput modes, three alternating trials of
  `center_kernel=pairs` versus `moments`, on 8-box and 64-box layouts, both
  precisions, on an idle GPU with a pinned core and recorded host load.
- Record results in `BENCHMARKS.md`, `PROGRESS.md`, and this plan's outcome
  section; update `GPU_DEVELOPMENT_SUMMARY.md`.

## Expected outcome

| | Now, double | Now, single | Target (double) |
| --- | ---: | ---: | ---: |
| div(E) deposit, 3 passes | 95.8 ms | 70.5 ms | about 35-40 ms |
| Whole step (stage mode) | 318 ms | 268 ms | about 260 ms |

The target makes the default double-precision step faster than today's
single-precision step, without its energy-conservation cost. These are
estimates from operation counts and memory traffic; phase 1's measurement
decides whether phases 2 and 3 are worth doing.

## Outcome (2026-10-25)

Phase 1 is implemented as `DemoParticles::deposit_center_charge_moments`
and selected by `center_kernel=moments` (default; `pairs` keeps the
node-pair kernel). Groups of G lanes share a base cell, keep 35 double
accumulators in registers, combine them with an xor shuffle, and issue the
same 72 atomics per cell. Tuning on the reference case:

| Lanes per cell G | 4 | 8 | 16 | 32 |
| --- | ---: | ---: | ---: | ---: |
| div(E) deposit, 3 passes (ms/step, GPU 1) | 44.3 | **39.6** | 44.2 | 56.5 |

G = 8 with 128-thread blocks is used; 64- and 256-thread blocks were within
1%.

Reference case, stage mode, one RTX 4090 (GPU 0), two undisturbed trials
per configuration (another user's GPU job disturbed the others; an earlier
three-trial series on GPU 1 agrees: 99.5 -> 44.1 ms at G = 16):

| | Pairs, double | Moments, double | Pairs, single | Moments, single |
| --- | ---: | ---: | ---: | ---: |
| div(E) deposit, 3 passes | 96.1 ms | **38.9 ms** | 70.8 ms | 38.9 ms |
| Whole step | 319.2 ms | **262.0 ms** | 268.4 ms | 236.6 ms |

- The deposit is **2.47x** faster in double precision and the step 18%
  shorter (1.22x). The default double-precision step (262 ms) is now faster
  than the previous single-precision step (268 ms).
- `cell_precision=single` no longer affects this deposit: the moment kernel
  is double precision only, and its single-precision gain now comes from
  the response kernel alone (262.0 -> 236.6 ms).
- Throughput mode (one undisturbed pair): 0.3233 -> 0.2660 s per step.
  64-box layout: deposit 106.3 -> 44.4 ms, step 371.5 -> 309.7 ms.
- Nsight Systems, per pass: kernel 1.8 -> 0.69 ms per launch at G = 16
  (before tuning to G = 8), sorting 4.0-4.5 ms unchanged. The kernel is
  2.6x faster rather than the estimated 3-4x: the 35-value shuffle
  reduction and the particle gathers through the sort permutation remain.
- Clustered particles benefit too: with `cell_chunk_size=256` the deposit
  takes 1.12 ms against 3.24 ms for the atomic path (32^3 cells, 20
  particles per cell, `clustered=1`); the pair kernel gained nothing there.
- Solver work is identical (70 electric and 287 div(E) iterations,
  6 warnings).

Validation: serial references at the double-precision tolerance on one and
two ranks with and without chunking (errors about 1e-17); per-species
comparison with the atomic deposit within 1.5e-15 (uniform) and 1.2e-14
(clustered, with and without chunking); coupled comparisons of both kernels
against the atomic path in `validate_demo.py`; CPU build on four ranks;
compute-sanitizer memcheck and racecheck clean.

Next candidates: sorting is now about 4 ms of the remaining ~13 ms per
pass, so phase 3 (a compact bin-ordered copy, or binning inside
`correct_positions`) is the larger remaining item; phase 2 (atomic-free
gather) also makes the deposit bitwise reproducible.
