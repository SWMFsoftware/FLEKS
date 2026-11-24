# Evolved & Anisotropic Electron Pressure (`tests/electron_pressure`)

Verification suite for the hybrid-PIC electron pressure models in FLEKS:

1. **Isotropic scalar PDE (`#ELECTRONPRESSURE`)**:
   ```
   dPe/dt + div(u_e Pe) + (gamma_e - 1) Pe div(u_e)
       = (gamma_e - 1) [ eta |J|^2 - div(kappa_hat . grad(Te)) + H_ei ]
   ```
   with `u_e = (J_i - J) / rho_q`, SSP-RK2 advection + exponential compression,
   conservative face-centered heat conduction (`rkl2` default, `point-implicit`,
   `subcycle`), Ohmic heating (`#ELECTRONHEATING`, default `T`), and
   point-implicit electron–ion equilibration (`#ELECTRONCOLLISION`).

2. **Gyrotropic / anisotropic electron pressure (`#ELECTRONANISO`)**:
   - **`cgl`**: evolves `(Pe, PePar)` with CGL parallel/perpendicular compression,
     anisotropic parallel Spitzer conduction + conservative mirror-force exchange,
     marginal-stability firehose/whistler clamping, and `-div(P_e)/rho_q` in Ohm's law.
   - **`le2009`**: algebraic Le et al. (2009) trapped–passing EOS bridging
     Boltzmann and CGL limits as a function of `n/n0` and `B/B0`.

## Decks

| Deck | Model / Physics | Validation Criterion |
|------|-----------------|----------------------|
| `adiabatic` | Scalar PDE (`kappa0 = 0`, `B = 0`) | `Pe / rho^gamma` flat vs `Pe / rho` |
| `conduction` | Scalar PDE (`point-implicit`, `B = 0`) | `Pe / rho` flat vs `Pe / rho^gamma` |
| `crossfield` | Scalar PDE (`point-implicit`, `By != 0`) | Cross-field conduction suppressed (`Pe / rho^gamma` flat) |
| `conduction_rkl2` | Scalar PDE (`rkl2` default, `B = 0`) | Flat `Te`, `sum(Pe)` conserved to `< 1e-5`, zero checkerboard |
| `ohmic` | `#ELECTRONHEATING T` (`rk4`) | Resistive decay heats `Pe` and conserves `B^2/(8*pi) + Pe/(gamma-1)` |
| `pc` | `#ELECTRONHEATING T` (`rk4-pc`) | Predictor-corrector gives identical non-doubled Ohmic heating |
| `collision` | `#ELECTRONCOLLISION` | Matches analytical Braginskii relaxation ODE to `< 5%` |
| `biermann_2d` | 2D crossed `grad(ne) x grad(Te)` | Generates out-of-plane `Bz` via `-div(Pe)/rho_q` |
| `amr` | 2-level refined mesh (`refRatio = 2`) | Coarse/fine adiabatic `Pe / rho^gamma` consistency |
| `cgl_compression` | `#ELECTRONANISO cgl`, `Bx != 0`, `kappa0 = 0` | `PePar/PePerp` tracks `(rho/rho0)^2` and satisfies trace identity |
| `aniso_firehose` | `#ELECTRONANISO cgl`, `Tperp/Tpar0 = 0.05` | Conserves `Pe` and clamps `(PePar - PePerp)/P_B` to `2.0` |
| `aniso_whistler` | `#ELECTRONANISO cgl`, `Tperp/Tpar0 = 10.0` | Conserves `Pe` and clamps `PePerp/PePar - 1` to `0.36/beta_par^0.55` |
| `aniso_wave` | `#ELECTRONANISO cgl`, `HybridWave`, `Tperp/Tpar0 = 3` | Matches gyrotropic whistler-Alfvén dispersion shift `alpha = 1 - (beta_par-beta_perp)/2` |
| `le2009` | `#ELECTRONANISO le2009`, `IonAcousticWave` | Matches closed-form Le et al. (2009) `PePar/PePerp` and `TePar/TePerp` pointwise |

## Running

```bash
python3 tests/validate_tests.py --test=electron_pressure -n 2 --verbose
```
