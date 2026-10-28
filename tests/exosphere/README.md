# Exosphere Tests (spatial density profiles, radiation & attenuation)

## Description

Tests exercising spatial exospheric neutral profiles (exponential, Chamberlain,
power-law), UV/photoionization radiation attenuation (optical depth and Chapman
functions), geometric shadow cylinders, and dual-cell averaging for thin layers:

| Variant | Atmosphere / Radiation | Validation check |
|---|---|---|
| `aligned` | Exponential atmosphere $H = \Delta x / 33$ | Dual-cell averaging prevents grid aliasing; injected mass matches volume integral within 10% |
| `shifted` | Shifted by $\Delta x / 3$ | Invariance under grid translation within 10% |
| `corona` | Chamberlain corona with optical depth | Geometric shadow lit check; photoproduction outside body matches within 5% |
| `sourceppc` | Constant PPV (#CONSTANTPPV) + source | Injects exactly target particles per cell and step (#SOURCEPARTICLES) |
| `opticaldepth` | Plain $\exp(-\tau/\mu)$ attenuation | Radial $O^+$ and $H^+$ profile along subsolar line matches model; nightside production suppressed |
| `chapman` | Chapman function attenuation | Radial profile matches Chapman function model along subsolar line |
| `shadow` | 2D/3D photoionization with `#SHADOWCYLINDER` | Nightside shadow suppression ($< 20\%$ of dayside), photoion temperature $\sim 300\text{ K}$, and block boundary smoothness |

## Running

```bash
python3 tests/validate_tests.py --test=exosphere              # all variants
python3 tests/validate_tests.py --test=exosphere.opticaldepth # single variant
```
