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

Needs the Exo user source and a 3D AMReX build.
