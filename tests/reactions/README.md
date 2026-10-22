# Reaction Tests (analytic source and loss terms)

## Description

Uniform, periodic 1D plasmas in which a single reaction of the exosphere
source (`-u=Exo`) runs against its analytic solution. They exercise the
plasma-state coupling of the stand-alone source terms: every rate that depends
on the plasma (recombination, ion-neutral chemistry, charge exchange) has to be
evaluated from the particles, not from the frozen `#UNIFORMSTATE`.

| Variant | Solver | Reaction | Analytic check |
|---|---|---|---|
| `ionneutral` | hybrid | A+ + N -> B+, nu dt = 0.01 | N_A = N_A0 exp(-nu t); N_A + N_B conserved; B+ keeps the 100 km/s drift; no electron particles |
| `stiff` | hybrid | A+ + N -> B+, nu dt = 3 | N_A + N_B conserved exactly (exponential loss + stiff product factor) |
| `recombination` | full PIC | O2+ + e -> O + O | n = n0/(1 + k n0 t), k = k0 (Tref/Te)^0.5 with Te from the particles; electrons follow the ions |
| `recombination_hybrid` | hybrid | O2+ + e -> O + O | same, n_e and Te from the hybrid electron closure |
| `cx` | hybrid | H+ + H -> H + H+ (`#CHARGEEXCHANGE`) | mass conserved; u = u0/(1 + n sigma u0 t); no electrons created |
| `weightfloor` | hybrid | photoionization below `#SOURCEWEIGHTFLOOR` | injected mass S V t on average |

The checks read the per-species particle mass (`MpartN`) and energy
(`EpartN`) columns of `log_pic_*.log`.

## Running

```bash
python3 tests/validate_tests.py --test=reactions            # all variants
python3 tests/validate_tests.py --test=reactions.ionneutral # one variant
```
