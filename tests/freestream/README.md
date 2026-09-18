# Free-stream — uniform-equilibrium field-solver test

**Test directory:** `tests/freestream/`

A 1D periodic **free-stream** uniform equilibrium test:
- High-speed bulk flow $u_x = 100\text{ km/s}$ ($u_{\text{code}} = 1.0 = v_A$, Mach 1 Alfvenic flow) across a 3D oblique magnetic field $\vec{B} = (8, 6, 3) \times 10^{-9}\text{ T}$ ($|\vec{B}| \approx 10.44\text{ nT}$).
- Convective electric field $\vec{E} = -\vec{u} \times \vec{B} = (0, 3.0\times 10^{-4}, -6.0\times 10^{-4})\text{ V/m}$ balancing the Lorentz force so the plasma drifts in steady state without spurious charge separation or acceleration.
- Plasma beta $\beta \approx 1.0$ ($T = 314{,}000\text{ K}$): thermal pressure balances magnetic pressure, enabling rigorous validation of temperature and pressure stability.
- In **Full PIC (Original GMRES)** (`PARAM.in.full`), two kinetic species are modeled: ions ($H^+$, $q=1, m=1$) and electrons ($e^-$, $q=-1, m=0.04$). The standard semi-implicit GMRES solver is used without comoving-frame current or artificial viscosity.
- In **Full PIC (Upwind)** (`PARAM.in.upwind`), two kinetic species are modeled with the comoving frame current solver (`#COMOVING`), Lax-Friedrichs upwind viscosity (`#UPWINDE`, `#UPWINDB`), and current smoothing (`#SMOOTHJ`) to suppress high-frequency oscillations.
- In **Hybrid PIC** (`PARAM.in.hybrid`), kinetic ions and massless fluid electrons are simulated with the Ohm's law + Faraday advance.

## Configuration

The directory holds three parameter files that share the identical grid, normalization, timestepping, and initial equilibrium:

| File | Field solver | Particle Species | Description |
| --- | --- | --- | --- |
| `PARAM.in.full` | Full PIC implicit Maxwell/GMRES (`#SOLVEEM T`, `#HYBRIDPIC F`) | Species 0: ions ($m=1, q=1$), Species 1: electrons ($m=0.04, q=-1$) | Canonical semi-implicit GMRES solve without comoving/upwind terms |
| `PARAM.in.upwind` | Full PIC implicit Maxwell/GMRES (`#SOLVEEM T`, `#HYBRIDPIC F`) | Species 0: ions ($m=1, q=1$), Species 1: electrons ($m=0.04, q=-1$) | Full kinetic E-B solve with comoving current & upwind damping |
| `PARAM.in.hybrid` | Hybrid (Ohm + Faraday) (`#SOLVEEM F`, `#HYBRIDPIC T`, `#HALLTERM F`) | Species 0: ions ($m=1, q=1$), massless fluid electrons | Convective Ohm's-law advance |

## Validation

All variants are checked by `tests/freestream/validate.py`:
- **Energy log:** Total particle kinetic energy `Epart`, per-species energies `Epart0`/`Epart1`, magnetic energy `Eb`, and electric energy `Ee` are conserved within tight tolerances ($\pm 2\%$ for particle energy, $\pm 5\%$ for field energies).
- **Plot output:**
  - Bulk velocity: $u_{xS0}$ (and $u_{xS1}$ for full PIC) is preserved to within 5%; transverse velocities remain negligible noise.
  - Temperature / pressure: Mean scalar pressure $\langle p_{S0} \rangle$ (and $\langle p_{S1} \rangle$ for full PIC) is conserved to within 5%, ensuring no numerical heating or cooling.
  - Magnetic & electric fields: Normal field $B_x$ is strictly uniform and preserved; mean field components match the theoretical convective equilibrium, with spatial spread $\le 10\%$ in Full PIC.

## Physics scales

- Normalization: $l_{\text{NormSI}} = 10^5\text{ m}$ ($= d_i$), $u_{\text{NormSI}} = 10^5\text{ m/s}$ ($= v_A$), so $t_{\text{Norm}} = 1\text{ s}$ ($= \Omega_i^{-1}$).
- Ion density $n_i = 5\text{ cm}^{-3}$, electron density $n_e = 5\text{ cm}^{-3}$ (charge neutral).
- 3D oblique magnetic field: $B_x = 8.0\text{ nT}, B_y = 6.0\text{ nT}, B_z = 3.0\text{ nT}$ ($|\vec{B}| \approx 10.44\text{ nT}$; $\Omega_i \approx 1.0\text{ rad/s}$, gyroperiod $T_{ci} \approx 6.28\text{ s}$).
- Convective electric field: $E_x = 0, E_y = 0.30\text{ mV/m}, E_z = -0.60\text{ mV/m}$.
- Temperature $T = 314{,}000\text{ K}$ sets total plasma $\beta = (p_i + p_e) / (B^2 / 2\mu_0) \approx 1.00$.
- Grid: $L_x = 6.4 d_i$ with 32 cells (periodic in X, Y, Z, 1 cell in y/z).
- Fixed $dt = 0.02\text{ s}$ for $t_{\text{max}} = 10\text{ s}$ ($10\text{ code units} \approx 1.6 T_{ci}$). Fast bulk flow $u_x = 1.0$ traverses $1.56\times$ the periodic box.

## Relationship to other tests

- `whistler/`, `ohm/` — seeded wave tests built on related equilibria.

## Running

All three variants run together under the test name:

```bash
python3 tests/validate_tests.py --test=freestream          # all 3 variants
```

To run a **single variant**, append its token to the test name:

```bash
python3 tests/validate_tests.py --test=freestream.full      # full PIC original GMRES (PARAM.in.full)
python3 tests/validate_tests.py --test=freestream.upwind    # full PIC upwind (PARAM.in.upwind)
python3 tests/validate_tests.py --test=freestream.hybrid    # hybrid only (PARAM.in.hybrid)
```