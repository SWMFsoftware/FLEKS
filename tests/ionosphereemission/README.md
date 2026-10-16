# Ionospheric Emission Shell Test (`#IONOSPHEREEMISSION`)

## Description

An absorbing spherical body (`#BODY`, radius 1.5) sits inside an emission
sphere of radius 2 (`#IONOSPHEREEMISSION`, after the rhybrid
IonosphereInjector): H+ at 1e16 /s and O+ at 2e15 /s, T = 1e5 K, outward
half-Maxwellian, noon/night profile 1/0.1. Hybrid PIC with B = 0 and Te = 0,
so E = 0 and the ions fly ballistically into vacuum (3D, 24^3 cells, 2 s).

## Validation

From the PIC log (`mBodyInject`, `mBodyAbsorb`, `Mpart0`, `Mpart1`), with the
code mass of a proton `e^2 / (1e7 lNormSI m_p)`:

* the injected mass rate equals the requested rates (relative 1e-6), i.e. the
  shell weights are normalized over the plasma cells and the SI -> code
  conversion of the rate is right;
* the body absorbs nothing (the emission is outward);
* each species' mass in the box grows at its own rate (relative 1e-3).

Needs a 3D AMReX build.
