# Multi-Solver Domain Coupling (`#REGIONSOLVER`)

Tests spatial partitioning of a shared FLEKS domain across field solvers (`Implicit`, `Explicit`, `Hybrid`) and wave transmission across the solver seam using the `HybridWave` parallel right-hand whistler-Alfven eigenmode (`lNormSI = 1e5 m`, `uNormSI = 5e5 m/s`, `c / v_A ~ 5`, `Lx = 12.8 d_i`, `k d_i ~ 0.49`).

## Variants

1. `PARAM.in` (`solver_coupling`): Implicit Full-PIC ($x < 0$) coupled to RK4 Hybrid-PIC ($x \ge 0$) with evolved electron pressure (`#ELECTRONPRESSURE`), refreshed kinetic-electron collar (`>= 2` cells), masked `div(E)` cleaning, and unified Faraday `B` advance.
2. `PARAM.in.explicit_hybrid` (`solver_coupling_explicit_hybrid`): Explicit Predictor-Corrector Full-PIC ($x < 0$) coupled to RK4 Hybrid-PIC ($x \ge 0$).
3. `PARAM.in.implicit_explicit` (`solver_coupling_implicit_explicit`): Implicit Full-PIC ($x < 0$) coupled to Explicit Leapfrog Full-PIC ($x \ge 0$).
4. `PARAM.in.blended_amr` (`solver_coupling_blended_amr`): 2D periodic run with cosine blending across the seam (`nBlendCell = 2`) and an AMR refinement patch straddling $x = 0$.
5. `PARAM.in.three_way` (`solver_coupling_three_way`): 1D periodic three-way coupling across `x in [-9.6, 9.6]` with Implicit Full-PIC (`x < -3.2`), Explicit Leapfrog Full-PIC (`-3.2 <= x < 3.2`), and RK4 Hybrid-PIC with `#ELECTRONPRESSURE` (`x >= 3.2`).

## Validation Criteria

- **Energy stability**: finite non-negative `Etot`, `Ee`, `Eb`, `Epart`, `|Etot[-1] - Etot[0]| / Etot[0] < 0.15`, and `0.80 <= Eb[-1] / Eb[0] <= 1.20`.
- **Solver mask & divergence**: `solver` diagnostic matches each subdomain's `PicSolverMode` in the interior, `max |divB| < 1e-9` everywhere across all variants, and `max |divBc| < 1e-9` on the 1D runs.
- **Species & electron pressure continuity**: `rhoS0`, `rhoS1`, `pS0`, `pS1` positive across all subdomains (verifying synthesized Hybrid electron moments), and `Pe` positive and continuous within 30% when `#ELECTRONPRESSURE` is active.
- **Wave transmission**: balanced transverse wave amplitude across subdomains (`0.65 <= rms(|B_perp|_L) / rms(|B_perp|_R) <= 1.55` for 2-way, `max/min <= 1.65` for 3-way), mode-1 power fraction `>= 0.75`, and positive unwrapped phase advance across frames.
