# FLEKS Tools and Post-Processing

## Scripts

| Script | Purpose |
|---|---|
| `tools/amrex2tec.py` | AMReX plot → Tecplot |
| `tools/amrex2vtk.sh` | AMReX plot → VTK (ParaView) |
| `tools/tec2vtk.sh` | Tecplot → VTK |
| `tools/converter.py` | general data conversion |
| `tools/clean_dat.py` | clean up `.dat` output |
| `tools/format_all.py` | bulk reformat (C++ + Fortran), required before a PR |
| `tools/generate_compile_commands.py` | regenerate `compile_commands.json` |
| `tools/check_docs.py` | verify the agent documentation tree — runs in CI |

## Reading output

FLEKS writes AMReX block-structured output (plus IDL/`.h` for SWMF
post-processing). Do not use generic NetCDF loaders:

- Python: `flekspy` (`pip install flekspy`), integrates with Matplotlib and YT.
- Julia: `Batsrus.jl` (`using Batsrus; load("filename.out")`).
- ParaView: reads HDF5 natively, but handles AMReX block boundaries better with
  the BATSRUSReader plugin shipped by `flekspy`; converting to `.vtm` or
  Tecplot with the scripts above is often simpler.
