# Restart Regression Tests (`tests/restart/`)

Standalone checkpoint/restart regression tests verifying bitwise/round-off
trajectory continuity across all three FLEKS field solvers:

| Variant | Input Deck | Field Solver | Key Exercised Restart State |
| :--- | :--- | :--- | :--- |
| Implicit | `PARAM.in.implicit` | `#FIELDSOLVER implicit` | Nodal `E`/`B`, cell-centered `B`, multi-species particles |
| Explicit | `PARAM.in.explicit` | `#FIELDSOLVER explicit` (`leapfrog`) | Nodal `E`/`B`, cell-centered `B`, multi-species particles |
| Hybrid | `PARAM.in.hybrid` | `#FIELDSOLVER hybrid` (`rk4` + `#HALLTERM` + `#ELECTRONPRESSURE`) | Nodal `E`/`B`, cell-centered `B`, `nodePlasmaPrev`, smoothed moments, `centerPeState`, ion particles |

## Validation Protocol

For each variant, `tests/validate_tests.py` runs an uninterrupted 10-step
simulation split across two sessions (`#STOP 5` with `#SAVERESTART T`, then `#RUN`
to `#STOP 10` with `#SAVERESTART F`). `validate.py` then:

1. Verifies the cycle-5 checkpoint header (`PC/restartOUT/FLEKS1_restart.H`,
   `FLEKS1_amrex_restart.H`) and solver MultiFab/particle checkpoint files.
2. Copies `PC/restartOUT` to `PC/restartIN` and launches a fresh `FLEKS.exe`
   restart run driven by `#INCLUDE PC/restartIN/FLEKS1_restart.H` from cycle 5 to
   cycle 10 (matching the initial run's MPI rank count).
3. Asserts relative difference `< 1e-10` between the restarted and uninterrupted
   runs across overlapping diagnostic log energies (`Etot`, `Eb`, `Ee`, `Epart`
   for cycles `5..10`) and across every field/moment variable in the cycle-10
   `real8` plot output.
