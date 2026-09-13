# FLEKS ↔ SWMF Coupling Reference

> Canonical: `doc/DEVELOPING.md` §4 (Modify the SWMF interface) and
> `.agent/workflows/add-coupling-var.md`. This page owns the full entry-point
> tables and the Fortran/C++ interoperability pitfalls.

## Three-layer pattern

```
SWMF coupler (Fortran)
  → PC_wrapper.f90 / PT_wrapper.f90   (Fortran wrappers)
    → FleksInterface.cpp              (extern "C" entry points)
      → Domain / Pic / Particles      (core FLEKS C++)
```

GM → PC carries the MHD state (B, E, rho, u, p); PC → GM carries PIC moments
(rho, u, p, J).

## Files (`srcInterface/`)

| File | Language | Purpose |
|---|---|---|
| `FleksInterface.cpp` | C++ | C-linkage entry points; owns the global `fleksDomains` object |
| `PC_wrapper.f90` | Fortran | SWMF wrapper for the **PC** (PIC) component |
| `PT_wrapper.f90` | Fortran | SWMF wrapper for the **PT** (particle tracker) component |

## Entry points (`FleksInterface.cpp`)

| Function | Called by | Purpose |
|---|---|---|
| `fleks_init_mpi_` | `PC_set_param` | Initialize the MPI communicator |
| `fleks_init_` | `PC_init_session` | Initialize FLEKS with the simulation time |
| `fleks_read_param_` | `PC_set_param` | Pass `PARAM.in` content to FLEKS |
| `fleks_from_gm_init_` | `PC_put_from_gm_init` | Receive grid / normalization info from GM |
| `fleks_finalize_init_` | `PC_put_from_gm` | Finalize initialization after the first IC set |
| `fleks_run_` | `PC_run` | Advance to the target time |
| `fleks_set_state_var_` | `PC_put_from_gm` | Receive state variables from GM |
| `fleks_get_state_var_` | `PC_get_for_gm` | Provide moments to GM |
| `fleks_find_points_` | `PC_find_points` | Find the MPI ranks owning given coordinates |
| `fleks_set_dt_` | `PC_put_from_gm_dt` | Receive the coupler time step |
| `fleks_save_restart_` | `PC_save_restart` | Write restart files |
| `fleks_get_grid_` | `PC_get_grid` | Return grid node positions |
| `fleks_set_grid_info_` | `PC_put_from_gm_grid_info` | Receive AMR grid structure from GM |
| `fleks_end_` | `PC_finalize` | Clean up and finalize |

## Wrapper subroutines (`PC_wrapper.f90`)

| Subroutine | SWMF action | Purpose |
|---|---|---|
| `PC_set_param` | READ/CHECK | Read parameters, initialize MPI |
| `PC_init_session` | INIT | Initialize the session |
| `PC_run` | RUN | Advance to `TimeSimulationLimit` |
| `PC_save_restart` | SAVE | Save restart data |
| `PC_finalize` | FINALIZE | Release resources |
| `PC_put_from_gm` | COUPLE GM→PC | Receive the MHD state |
| `PC_get_for_gm` | COUPLE PC→GM | Send moments back to GM |
| `PC_put_from_gm_init` | INIT COUPLE | Receive initial grid parameters |
| `PC_put_from_gm_grid_info` | GRID COUPLE | Receive AMR grid updates |
| `PC_put_from_gm_dt` | DT COUPLE | Receive the time step |
| `PC_find_points` | POINT QUERY | Find the owning processors |
| `PC_get_grid_info` | GRID QUERY | Return grid dimension info |

## Fortran ↔ C++ pitfalls

1. **Indexing & layout:** C++ is 0-indexed and row-major, Fortran is 1-indexed
   and column-major. Check every `Data_VI(iVar, iPoint)` offset.
2. **Type mismatches:** `double` in C++ ↔ `real(8)` in Fortran; watch
   "possible change of value in conversion" warnings.
3. **Array rank remapping:** passing multidimensional arrays where 1-D
   contiguous memory is expected requires explicit contiguity or rank-1
   slicing, otherwise gfortran/ifort raise rank-mismatch errors.
4. **Name mangling:** C functions meant for Fortran use `extern "C"` with a
   trailing underscore (`fleks_run_`) unless `bind(C, name="...")` is explicit.
5. **ISO_C_BINDING:** prefer `use iso_c_binding` in interface blocks when
   adding new bindings.
6. **Units / nVar counts:** FLEKS is CGS internally while GM may be SI or
   normalized — convert in `FluidInterface`. Update `nVar` on *both* the
   Fortran wrapper and the C++ side when adding a variable.
7. **Query points:** FLEKS sends **all** of its grid nodes (valid + ghost +
   shared) to GM, including nodes inside `r < rBody`. GM must handle points
   outside its radial domain (clamp to the inner boundary) — a `-1` block
   result must never be used as a buffer index.

## Adding an interface function

1. Add the C++ function in `FleksInterface.cpp` with `extern "C"` linkage:

   ```cpp
   extern "C" int fleks_new_function_(int *param) { /* ... */ return 0; }
   ```

2. Declare it in `include/FleksInterface.h`.
3. Add the interface block and wrapper subroutine in `PC_wrapper.f90`:

   ```fortran
   use iso_c_binding, only: c_int
   interface
     integer(c_int) function fleks_new_function(param) bind(C)
       integer(c_int), intent(in) :: param
     end function
   end interface
   ```

4. Rebuild (`make LIB -j8`) and validate with `make test16_3d` from the SWMF
   root.

See also `.agent/workflows/add-coupling-var.md` for adding an exchanged
variable end-to-end.
