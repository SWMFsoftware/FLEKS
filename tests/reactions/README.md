# Reaction Tests (analytic source and loss terms)

## Description

Uniform and multi-species plasmas in which reactions of the exospheric
source (`-u=Exo`) run against analytic solutions and physical balances. They
exercise the plasma-state coupling of the stand-alone source terms: every rate
that depends on the plasma (recombination, ion-neutral chemistry, charge
exchange, electron impact) is evaluated from the particles and plasma state.

| Variant | Solver | Reaction / Physics | Validation check |
|---|---|---|---|
| `ionneutral` | hybrid | A+ + N -> B+, nu dt = 0.01 | N_A = N_A0 exp(-nu t); N_A + N_B conserved; B+ keeps the 100 km/s drift; no electron particles |
| `stiff` | hybrid | A+ + N -> B+, nu dt = 3 | N_A + N_B conserved exactly (exponential loss + stiff product factor) |
| `recombination` | full PIC | O2+ + e -> O + O | n = n0/(1 + k n0 t), k = k0 (Tref/Te)^0.5 with Te from particles; electrons follow ions |
| `recombination_hybrid` | hybrid | O2+ + e -> O + O | same, n_e and Te from hybrid electron closure |
| `cx` | hybrid | H+ + H -> H + H+ (`#CHARGEEXCHANGE`) | mass conserved; u = u0/(1 + n sigma u0 t); no electrons created |
| `weightfloor` | hybrid | photoionization below `#SOURCEPARTICLES` minWeightRatio | injected mass S V t on average |
| `chemistry` | full PIC | 10-reaction Mars chemistry network (H+, O+, O2+, CO2+) | All 4 ion species show energy change; O2+ cross-species CX source growth (>1); CO2+ energy change |
| `electronimpact` | full PIC | Electron-impact ionization source | O+ energy growth (source active) |
| `cx_multi` | full PIC | 2-component Chamberlain neutral atmosphere CX | O+ energy growth, H+ stability, and surface-peaked source profile in plot output |

## Running

```bash
python3 tests/validate_tests.py --test=reactions           # all variants
python3 tests/validate_tests.py --test=reactions.chemistry # single variant
```
