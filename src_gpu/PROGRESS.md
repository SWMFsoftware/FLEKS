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
