# Initial GPU measurements — 2026-10-24

These are short exploratory measurements on NVIDIA GeForce RTX 4090 GPUs
(24 GiB each), AMReX `45a1fc5e1e5a`, nvcc sm89, double precision, and host-staged
MPI. They are correctness-preserving comparisons, not production scaling
results or statistically established speedups. No other demo or sanitizer
process was running during this series.

All eight cases passed their conservation and storage/layout comparisons.
Common arguments were:

```text
mode=benchmark n_cell=16 max_grid_size=8 particles_per_cell=8
n_step=8 n_iter=10 correct_div_e=0 diagnostic_interval=8
```

Each case has 65,536 particles, two species, and a uniform cold electron-drift
initial condition, except that the clustered cases compress positions into
the central region. Coupled timings exclude the first step and diagnostics,
include charge-history bookkeeping, and synchronize at measured phase
boundaries. Two-rank times use the maximum rank elapsed time.

| Case | Changed arguments | Mean timestep (ms) |
| --- | --- | ---: |
| Uniform baseline | none | 6.072 |
| Uniform fused gather/push | `fused_push=1` | 5.706 |
| Uniform sorting each step | `sort_interval=1` | 6.520 |
| Uniform buffered response | `buffered_response=1` | 14.242 |
| Clustered atomic | `clustered=1` | 7.653 |
| Clustered buffered | `clustered=1 buffered_response=1` | 388.876 |
| Smaller boxes | `max_grid_size=4` | 98.612 |
| Two GPUs | `mpirun -np 2` | 68.085 |

Baseline phase means were 0.768 ms for response assembly, 3.520 ms for matrix
exchange, 0.958 ms for the electric solve, and 0.364 ms for separate particle
gather plus push/migration. The tiny problem is dominated by communication and
launch/control overhead, so the two-GPU case provides a communication test,
not evidence of useful strong scaling.

For the baseline state's microbenchmarks:

| Operation | Mean time (ms) |
| --- | ---: |
| Atomic response assembly including clearing | 0.764 |
| Cell-buffered response including sorting/setup/clearing | 10.017 |
| Custom 243-real matrix application | 0.0570 |
| Component MultiFab matrix application | 0.0434 |
| Ghost exchange, width 1 | 0.00749 |
| Ghost exchange, width 3 | 0.00770 |
| AoS frozen-position gather/push, electron species | 0.0257 |
| Component-SoA scratch gather/push, electron species | 0.0279 |

The matrix layouts produced matching results; the particle layout test's
velocity difference was zero. SoA setup/transpose was measured separately
(0.0736 ms for that species) and excluded from kernel time. This is a scratch
layout test, not an AMReX PureSoA container/migration benchmark.

The cell-buffered implementation serially accumulates all particles in a cell
before flushing 64 tensor blocks. It reduces global atomics but introduces
large thread-local storage and poor parallelism in densely occupied cells.
The clustered result is a concrete reason to keep the atomic path as default;
this buffering strategy is an experimental comparison, not an optimization to
transfer directly into FLEKS. Each thread holds 600 doubles (4.8 KB) of
private buffers, which spill from registers to local memory, and each call
allocates range arrays and synchronizes the stream once per box. The result
therefore measures this implementation, not buffering in general; it does not
rule out block- or warp-level reduction in shared memory. Similarly, the small fused-kernel difference
needs longer repeated trials before attributing a speedup to fusion.

Logs and structured results from this run are in
`/tmp/fleks-gpu-benchmarks/`; rerun the documented commands to reproduce them.
The most useful next profiling target is custom matrix communication at larger
per-box workloads, followed by cooperative buffering with bounded work per
thread rather than a serial per-cell accumulation.

## Cell-reduced mass-matrix assembly — 2026-10-24

`response=cell` (`DemoParticles::accumulate_response_cell`) bins particles by
cell with an AMReX `DenseBins` counting sort (an index permutation; particles
are not moved) and launches one 64-thread block per cell. Thread `t` owns the
node pair `(t/8, t%8)` and accumulates its 3x3 tensor in registers over every
particle in the cell; particle weights, `alpha` and `jHat` terms are staged in
shared memory in chunks of 64. Each cell then issues 576 atomics, with
distinct addresses inside a warp, instead of 576 atomics per particle with
all lanes of a warp hitting the same addresses. The CPU build falls back to
the per-particle path.

One RTX 4090, 64^3 cells, 100 particles per cell per species (52.4M
particles), 8 boxes of 32^3, Langmuir wave, `b0=0.5`, `dt=0.1`, div(E)
correction on, `solver_tolerance=1e-6`; steps 2-6 timed:

| Assembly | Seconds per assembly | Timestep (s) | Particles/s |
| --- | ---: | ---: | ---: |
| `atomic` | 0.536 | 1.064 | 49.3M |
| `cell` | 0.073 | 0.603 | 87.0M |
| `buffered` | 0.672 | (not run) | |
| 32 CPU cores, `atomic` | 0.443 | 2.163 | 24.2M |

All three GPU assemblies agree entry by entry within 7e-17, and the coupled
diagnostics after six steps match the atomic run to 2e-15. Clustered particles
(`clustered=1`, 10 particles per cell, about 1,250 per occupied cell) remain
favorable: 0.050 s for `cell` versus 0.112 s for `atomic`.

Nsight Systems attributes about 0.5 ms per box and species to the counting
sort and about 5.9 ms to the block kernel. Hardware counters are restricted to
administrators on this host (`RmProfilingAdminOnly=1`), so Nsight Compute
could not run. A simple throughput model (FP64 on this consumer GPU is 1/64 of
FP32) estimates about 5 ms: roughly 3.6 ms for the ten FP64 instructions per
particle and node pair, plus 1.5 ms of per-particle staging. The kernel is
therefore likely FP64-bound. Exploiting the symmetry of the pair weights
(36 unique pairs instead of 64) could remove about 40% of the pair work; at
this point, however, the divE charge/matrix deposit (0.22 s) and the moment
deposit (0.12 s) cost more per step and would benefit from the same per-cell
reduction.

## Cell-reduced particle deposits — 2026-10-25

`deposit=cell` applies the per-cell block reduction to the div(E)
charge/correction-matrix deposit (`deposit_center_charge_cell`, bins over the
half-cell-shifted base cells including one low ghost layer) and to the
ten-moment deposit (`deposit_cell`). Same reference case as above, with
`response=cell`, the process pinned to one idle core:

| Stage | `deposit=atomic` (s) | `deposit=cell` (s) | Speedup |
| --- | ---: | ---: | ---: |
| div(E) charge + matrix, 3 passes | 0.219 | 0.111 | 2.0x |
| Moment deposit | 0.122 | 0.016 | 7.5x |
| Final charge deposit (per particle in both) | 0.013 | 0.013 | 1.0x |
| **Timestep** | **0.599** | **0.387** | **1.55x** |

Throughput is 135.5M particles per second, 5.6x the 32-core CPU run. The
div(E) kernel is FP64-bound after the atomic reduction, and clustered
particles do not benefit; `DEPOSIT_OPTIMIZATION_PLAN.md` ("Outcome") gives
the profile and the remaining options.

A competing 16-rank CPU job started during the first measurements and
doubled the timestep (1.28 s, with the electric solve at 0.195 s instead of
0.047 s): the host drives every Krylov iteration, so GPU timings on this
shared machine are sensitive to CPU load. The numbers above are from a quiet
period and were reproduced twice.
