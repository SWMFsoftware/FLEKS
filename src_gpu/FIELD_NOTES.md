# Field and solver implementation notes

## Production algorithm audit

The standalone demo uses the same dual grid as FLEKS: nodal E and interpolated
B, cell-centered evolved B, and quarter-averaged transverse finite differences.
`DemoFieldKernels.h` implements the direction-generic form of the expressions in
`src/GridUtility.cpp`. The difference is second order for smooth periodic
fields; it is not the component-wise staggering of a Yee mesh.

The demo diagnoses div(B) directly from evolved cell-centered B onto nodes.
This deliberately avoids the interpolation/smoothing used by production's
existing diagnostics: interpolating a longitudinal checkerboard B to nodes
first can erase its divergence. A checkerboard regression checks this case.

The electric operator is

```
A(E) = E + (theta dt)^2 [curl_cn(curl_nc(E))
       + coefDiff grad_cn(div_cc(avg_nc(E)) - div_nc(E))]
       + 4 pi theta dt M E.
```

The first curl composition equals `grad_cn(div_nc(E)) - lap_nn(E)` for these
commuting uniform-grid derivatives, so this reproduces `Pic::update_E_matvec`
without allocating a separate scalar gradient per electric component. The
cell-to-cell divergence includes production's 3x3 transverse averaging of
central differences; using a plain central divergence would change the
operator for mixed-dimensional modes. The additional divergence mixture is necessary: production `FieldSolver` defaults
to `coefDiff = 0.1`, not zero. The demo exposes `coef_diff`, defaulting to 0.1.
The RHS is `Eold + theta dt [curl_cn(Bold) - 4 pi jHat]`. The solution is the
theta-stage field, followed by `Enew = (Eth - (1-theta) Eold)/theta` and
`Bnew = Bold - dt curl_nc(Eth)`.

Normalized Gaussian units retain FLEKS's explicit `4 pi` factors and use c=1.
Particle response `M` already includes inverse cell volume, macro charge, and
`q dt/(2m)` in its coefficients. `jHat` is deposited as raw macro charge times
rotated velocity; `finish_response()` sums contributions and divides by cell
volume exactly once. Field energy is integrated as `(E^2+B^2)/(8 pi)`.

## div(E) correction

An audit found `Pic::rhoTheta = 0.51`, independently of electric-field theta.
The correction RHS is `(div(Enew)-4 pi rhoN)/rhoTheta`, where
`rhoN=(1-rhoTheta)rhoOld+rhoTheta rhoNew`. The 27-component scalar response
operator is `16 pi^2 centerMM phi` with the production component ordering
`(di+1)*9+(dj+1)*3+dk+1`. The periodic constant mode is removed from the RHS and
potential. Particle position correction and response deposition live in
`DemoParticles` and follow `Particles::divE_correct_position` and
`accumulate_mass_matrix_contribution`.

The separate cell-centered GMRES shares the implementation used for nodal E.
Its tolerance and iteration count are separate settings, because production
intentionally uses a looser inner correction solve. Nonconvergence remains a
reported failure rather than silently accepting a potential.

## GPU storage and communication findings

Production `Arr1D` in `include/Array1D.h` lacks GPU method annotations. The demo
keeps a local layout-compatible `DemoRealMM` with device-callable indexing,
addition and initialization; production headers are untouched. It stores 243
reals per nodal point in `FabArray<BaseFab<DemoRealMM>>`.

Installed AMReX's `HasAtomicAdd<DemoRealMM>` is false. Its GPU FabArray
communication automatically selects the deterministic, masked component-object
addition path for custom types. No specialization of AMReX's atomic traits or
unsafe object-size atomics is required. Particle scatter atomically adds
individual Real components. The communication smoke check scatters each cell's
1/8 contribution to its eight nodes for all 243 components and verifies the
uniform periodic result in valid and ghost nodes after SumBoundary and
FillBoundary, including box/rank seams.

Nodal fields computed from stencils use OverrideSync, then FillBoundary.
Deposits use SumBoundary once, then OverrideSync and FillBoundary. Scalar
centered response uses the corresponding cell-centered sum/fill operations.
Do not repeat SumBoundary on a field that already holds its globally summed
nodal value.

## Solver interface findings

Production `LinearSolver` stores flattened host `double*` vectors. The demo
keeps restarted GMRES basis vectors, residuals and matvec scratch in persistent
MultiFabs on the device. Only the Hessenberg matrix, Givens rotations and
convergence scalars live on the host. Twice-modified Gram-Schmidt uses globally
reduced owner-mask dot products; unique periodic nodal ownership avoids box and
periodic-face duplication. A true device matvec residual is recomputed after
each restart/convergence candidate.

The production port needs an operator-on-device-vector interface and device
vector BLAS/reductions; changing only `MATVEC` while retaining host vector
storage would continue to incur whole-vector transfers each iteration.

## Built-in checks

`DemoFields::run_checks()` verifies both manufactured curls, cell-to-node B
interpolation, unique periodic nodal counts, custom243 communication, a nonzero
mass response in the complete electric operator, nodal GMRES solution and true
residual, and a manufactured scalar potential with a seven-point negative
Laplacian encoded through the general 27-component response storage. Analytic
Fourier responses use exact discrete sine factors rather than continuum
approximations.
