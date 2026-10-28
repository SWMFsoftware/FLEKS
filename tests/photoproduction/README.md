# Unresolved Photoproduction Test

## Description

An exponential neutral atmosphere with a scale height of 1 km around a body
of radius 150 km, on 33 km cells (the Mars ionosphere: 6.5 km vs 849 km),
photoionized at a uniform rate (hybrid, no fields, 4 steps). Sampled at the
nodes, the source aliases with the grid: the production depends on where the
nodes fall relative to the thin layer. The source terms use averages over the
dual cell of every node and deliver the production of the nodes at the body
surface to the plasma cells, so the injected mass per step must match the
analytic column `nu n0 4 pi (R^2 H + 2 R H^2 + 2 H^3)` within 10%, for a grid
aligned with the planet (`aligned`) and one shifted by 0.3 cells (`shifted`).

The `corona` variant uses a `Mixed` #EXOSPHERE: a transparent Chamberlain
corona `n0 exp(-H0 (1/R - 1/r))` (local scale height 15 km at the body, H0 =
10 R) next to an exponential component that is not ionized, with #OPTICALDEPTH
on. Only the geometric shadow of the planet behind it is dark, so the
production must equal `nu` times the integral of the corona over the box
outside the body and the shadow cylinder, within 5% (measured 0.998; a dark
nightside, the earlier behavior, gives 0.66).

Needs the Exo user source and a 3D AMReX build.
