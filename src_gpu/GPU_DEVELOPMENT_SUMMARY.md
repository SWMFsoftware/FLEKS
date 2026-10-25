# GPU development summary, 2026-10-23 to 2026-10-25

This summarizes the GPU work on branch `gpu` over about 48 hours: what each
change did and how much it sped up the code. Details live in the documents
linked at the end.

## Result in one line

One RTX 4090 now runs the reference semi-implicit PIC step in **0.318 s**
(default double precision) or **0.268 s** (opt-in single-precision matrix
sums), down from **1.066 s** for the first GPU version: **3.4x** and **4.0x**
faster. It is **4.8x** (double) and **5.7x** (single) faster than all 64 CPU
cores of the same workstation running the same code (1.525 s).

## Reference case and how speed was measured

- Hardware: NVIDIA GeForce RTX 4090 (24 GB, FP64 at 1/64 of FP32) and an
  AMD Ryzen Threadripper 7980X (64 cores), one shared workstation.
- Workload: 3D periodic grid of 64^3 cells, two kinetic species with 100
  particles per cell each (52.4M particles), Langmuir wave, B0 = 0.5,
  dt = 0.1 (c dt/dx = 6.4), three div(E) correction passes per step,
  electric-solver tolerance 1e-6 (production default).
- GPU: one MPI rank, 8 boxes of 32^3 unless noted. CPU: one MPI rank per
  physical core, 64 boxes of 16^3, the same source built against a CPU-only
  AMReX.
- "Step" is the mean time per timestep over steps 2-6 in benchmark stage
  mode, which synchronizes the GPU around each phase. Throughput mode, which
  avoids those barriers, gives slightly larger absolute times on this case;
  each comparison below uses one protocol only.
- The machine was shared. Runs disturbed by other users' jobs were discarded
  and repeated; reported trials agree within about 1-2%.

## Timeline of changes and gains

### 1. Baseline AMReX GPU test (`ca0766a`, `3619693`)

A standalone AMReX test of GPU field operations (fills, reductions, ghost
exchange, stencils, copies) with host reference checks, preserved as the
starting point. No PIC yet; no speed figure.

### 2. Standalone semi-implicit PIC demo (`3bde4eb`)

A self-contained single-level port of the FLEKS semi-implicit algorithm in
`src_gpu/`, built with `make GPUTEST` and independent of `libFLEKS.a`:
FLEKS dual-grid field operators, AoS particles with the production layout,
CIC gather and Boris push, the 243-component mass matrix and jHat, a
device-resident GMRES electric solve, the B update, and the div(E) position
correction. Every kernel was checked line by line against production FLEKS
during review.

- First GPU run of the reference case: **1.066 s per step**.
- Same code on 32 CPU cores: 2.163 s (GPU 2.0x faster).
- Profile: 82% of the GPU step was particle-to-grid deposition. Mass-matrix
  assembly alone took 0.536 s, slower than on the CPU, because every
  particle issued 576 double-precision atomic adds and the 32 lanes of a
  warp hit the same addresses.

### 3. Stronger validation (`4173abb`)

No speed change. Automated gates tightened to about ten times the measured
values (plasma frequency, theta = 0.5 energy drift, dt refinement, div(E)
reduction, decomposition invariance); added nonuniform and large-dt
(c dt/dx = 8) cases and oblique-B references that cover every term of the
Boris rotation.

### 4. Cell-reduced mass-matrix assembly (`767a754`)

`response=cell`: a counting sort bins particles by cell, and one 64-thread
block per cell accumulates the 64 node-pair 3x3 tensors in registers before
issuing 576 atomics per cell instead of per particle.

| | Before | After | Gain |
| --- | ---: | ---: | ---: |
| Mass-matrix assembly | 0.536 s | 0.073 s | **7.4x** |
| Step | 1.064 s | 0.603 s | **1.76x** |

Also: a div(E) solve that reaches its iteration limit now warns instead of
failing, as in production FLEKS.

### 5. Cell-reduced particle deposits (`d4ef247`)

`deposit=cell`: the same per-cell reduction for the div(E) charge and
correction-matrix deposit and for the ten-moment deposit. The charge-only
deposit stays per particle, because the cell version measured slower.

| | Before | After | Gain |
| --- | ---: | ---: | ---: |
| div(E) charge + matrix deposit, 3 passes | 0.219 s | 0.111 s | **2.0x** |
| Moment deposit | 0.122 s | 0.016 s | **7.5x** |
| Step | 0.599 s | 0.387 s | **1.55x** |

The div(E) deposit became limited by FP64 arithmetic rather than by
atomics, which is why it gained less than the moments.

### 6. Pair reuse, fused push, and benchmark tooling (`422512c`)

Spatial pair symmetry in the cell kernels (`pair_symmetry=1`, 36 unique node
pairs with cooperative lane groups), fused gather and push (`fused_push=1`),
optional cell chunking and charge folding, throughput-mode timing, and the
repeatable `measure_performance.py` runner.

| | Before | After | Gain |
| --- | ---: | ---: | ---: |
| Step, throughput mode | 0.416 s | 0.355 s | **1.17x** |

Recommended flags since then: `response=cell deposit=cell pair_symmetry=1
fused_push=1`. In stage mode the recommended configuration ran at
0.344 s per step.

### 7. GMRES without host waits (`d9b7af2`)

The div(E) solve was the one phase where the GPU lost to the CPU (35.6 ms
versus 28.2 ms on 64 cores). AMReX dot products, vector operations, and
`MFIter` loops each made the host wait for the GPU: about 46 waits per
GMRES iteration. The batched twice-classical-Gram-Schmidt solver was
rewritten so both projection passes, their updates, and the norm stay on
the GPU, leaving one host wait per iteration. It is now the default
(`gmres_batched=1`). Also: zero-initial-guess shortcut, fused right-hand
side, per-solver counters, and opt-in `div_e_method=jacobi_gmres|pcg`.

| | Before | After | Gain |
| --- | ---: | ---: | ---: |
| Host waits per div(E) iteration | 43.8 | 1.1 | 40x fewer |
| div(E) solves, 3 passes (8 boxes) | 33.8 ms | 13.8 ms | **2.45x** |
| div(E) solves, 3 passes (64 boxes) | 35.6 ms | 15.3 ms | **2.33x** |
| Electric solve | 46.3 ms | 42.8 ms | 1.08x |
| Step, stage mode | 0.344 s | 0.320 s | **1.07x** |
| Step, throughput mode | 0.345 s | 0.320 s | **1.08x** |

Iteration counts and warnings are identical to the previous solver. The
opt-in `pcg` (Jacobi-preconditioned conjugate gradients) makes the div(E)
solves 2.8x cheaper again (4.9 ms) but changes the applied correction, so it
is a production decision and stays off.

### 8. Single-precision cell matrix sums (`d4e6584`, opt-in)

`cell_precision=single`: the response and div(E) correction-matrix kernels
stage their pair inputs and sum each cell in float, then add into the
double-precision matrices. Charge, jHat, and moments stay in double
precision because they nearly cancel between electrons and ions.

| | Double | Single | Gain |
| --- | ---: | ---: | ---: |
| Mass-matrix assembly | 50.6 ms | 25.2 ms | **2.0x** |
| div(E) charge + matrix deposit | 95.8 ms | 70.5 ms | **1.36x** |
| Step, stage mode | 0.318 s | 0.268 s | **1.19x** |
| Step, throughput mode (GPU 1) | 0.332 s | 0.279 s | **1.19x** |

Trade-off: the theta = 0.5 energy drift over 600 steps rises from 2.2e-15
to 8.6e-9, because float sums break the exact consistency between the mass
matrix and the particle push. Default stays double. Tensor cores were
assessed and do not help on this card (no FP64 support; TF32 is no faster
than FP32 here).

## Cumulative speedup

Stage-mode step time, reference case, one RTX 4090, 8 boxes:

| After change | Step (s) | Gain over previous | Gain over first GPU version |
| --- | ---: | ---: | ---: |
| 2. First PIC demo | 1.066 | | 1.0x |
| 4. Cell-reduced mass matrix | 0.603 | 1.77x | 1.8x |
| 5. Cell-reduced deposits | 0.387 | 1.56x | 2.8x |
| 6. Pair reuse + fused push | 0.344 | 1.13x | 3.1x |
| 7. GMRES without host waits | 0.318 | 1.08x | **3.4x** (default) |
| 8. Single-precision matrix sums | 0.268 | 1.19x | **4.0x** (opt-in) |

## Current speed of each part, GPU versus CPU

Fresh runs of the current code, medians of three trials. CPU: one rank per
physical core, 64 boxes of 16^3. GPU: 8 boxes of 32^3, recommended flags.

| Part (ms per step) | CPU 64 cores | CPU 32 cores | GPU double | GPU single | GPU single vs 64 cores |
| --- | ---: | ---: | ---: | ---: | ---: |
| div(E) charge + matrix deposit, 3 passes | 425.1 | 667.7 | 95.8 | 70.6 | 6.0x |
| Mass-matrix assembly | 318.2 | 358.5 | 50.6 | 25.2 | 12.6x |
| div(E) position correction + migration, 3 passes | 268.7 | 286.3 | 50.2 | 50.2 | 5.4x |
| Electric-field GMRES solve | 240.6 | 297.8 | 42.7 | 42.7 | 5.6x |
| Particle push + migration | 90.0 | 91.4 | 24.3 | 24.3 | 3.7x |
| Mass-matrix exchange | 63.7 | 61.8 | 11.7 | 11.7 | 5.4x |
| Moment deposit | 56.6 | 90.6 | 16.2 | 16.2 | 3.5x |
| Final charge deposit | 32.7 | 39.2 | 13.0 | 13.0 | 2.5x |
| div(E) GMRES solves, 3 passes | 23.5 | 33.3 | 13.7 | 13.7 | 1.7x |
| Magnetic-field update | 1.3 | 1.6 | 0.2 | 0.2 | 6.9x |
| **Step** | **1525** | **1930** | **318.5** | **267.9** | **5.7x** |

On the CPU's own layout (64 boxes of 16^3) the GPU takes 370 ms (double) and
325 ms (single), 4.1x and 4.7x faster than 64 cores; small boxes mainly
slow the GPU's mass-matrix exchange. Doubling the CPU from 32 to 64 cores
gains only 1.27x.

## Correctness and validation

Every change kept the physics: solver iteration counts and warnings are
identical across the double-precision optimizations, and coupled
diagnostics match the reference paths to roundoff. The validation suite
(`validate_demo.py --mpi --physics --reference`) passes and now covers:

- serial references for every particle kernel, including oblique B and both
  precisions;
- entry-by-entry comparisons of every optimized kernel with the atomic path;
- box and MPI decomposition invariance on nonuniform and large-dt cases;
- plasma frequency, second-order dt convergence, and theta = 0.5 energy
  conservation;
- div(E) correction effectiveness for every solver method;
- a comparison against production FLEKS CPU operators and a one-step CPU
  fixture.

compute-sanitizer memcheck and racecheck are clean for the new kernels, and
a CPU build of the same source passes all checks on four ranks.

## Remaining bottlenecks and options

- The div(E) correction loop is half of the step (deposit, solve, and
  position correction with particle migration, three passes each).
- The div(E) deposit is limited by its per-particle setup arithmetic in
  double precision; moving that setup to float is the next candidate.
- Particle migration after each correction pass (about 50 ms) is the
  largest phase not yet optimized.
- Opt-in options that change results and need a production decision:
  `div_e_method=pcg` and `cell_precision=single`.
- On A100/H100-class GPUs, FP64 runs at half the FP32 rate, so the
  single-precision gain would largely disappear; FP64 tensor cores there
  could run the per-cell matrix sums in full precision.

## Defaults and recommended flags

- Defaults: `response=atomic deposit=atomic gmres_batched=1
  div_e_method=gmres cell_precision=double` (production-matching results).
- Recommended for speed with production-matching results:
  `response=cell deposit=cell pair_symmetry=1 fused_push=1`.
- Fastest, with the accuracy trade-off above: add `cell_precision=single`.

## Further reading

- [README.md](README.md): build, run, and all options.
- [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md): scope of the demo.
- [DEPOSIT_OPTIMIZATION_PLAN.md](DEPOSIT_OPTIMIZATION_PLAN.md) and
  [DIVE_SOLVER_OPTIMIZATION_PLAN.md](DIVE_SOLVER_OPTIMIZATION_PLAN.md):
  plans with measured outcomes.
- [BENCHMARKS.md](BENCHMARKS.md): all measurements and their protocols.
- [VALIDATION.md](VALIDATION.md): what each check verifies.
- [PROGRESS.md](PROGRESS.md): dated log of the work.
