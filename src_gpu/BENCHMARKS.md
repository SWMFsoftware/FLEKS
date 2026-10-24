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
