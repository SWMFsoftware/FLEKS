#DESCRIPTION
Hybrid-PIC upwind (Rusanov) dissipation of the Ohm's-law E in a frozen plasma
(no particles) with the MC-limited reconstruction and a fixed signal speed.
The limited scheme is nonlinear; its damping rate must lie between the
unlimited second-order and the first-order rates.
#INITFROMSWMF
F                       initFromSWMF

#GEOMETRY
-3.2                    xMin  [code units; Lx = 6.4 d_i]
 3.2                    xMax
-0.05                   yMin
 0.05                   yMax
-0.05                   zMin
 0.05                   zMax

#NCELL
32                      nCellX
1                       nCellY
1                       nCellZ

#MAXBLOCKSIZE
8                       maxBlockSizeX
1                       maxBlockSizeY
1                       maxBlockSizeZ

#PERIODICITY
T                       isPeriodicX
T                       isPeriodicY
T                       isPeriodicZ

#SOLVEEM
F                       solveEM

#HYBRIDPIC
T                       useHybridPIC

#HYBRIDUPWIND
T                       useHybridUpwind
mc                      upwindLimiter
1.0                     upwindCoefWhistler
0.5                     upwindFixedSpeed  [code units]
1.0                     upwindSensorFloor

#RESISTIVITY
0.0                     etaResistivity  [m^2/s]

#ELECTRONTEMPERATURE
0.0                     electronTemperature  [eV]
1.0                     electronGamma
1.0                     electronDensity0

#BSUBCYCLE
1                       nBSubcycle

#FIELDINTEGRATOR
rk4                     fieldIntegrator

#PLASMA
1                       nSpecies
1.0                     mass
1.0                     charge

#UNIFORMSTATE
5.0                     rho [/cc]
0.0                     ux [km/s]
0.0                     uy [km/s]
0.0                     uz [km/s]
10000.0                 T [K]
1.0439e-8               bx [T]          Guide field Bx0
0.0                     by [T]
0.0                     bz [T]
0.0                     ex [V/m]
0.0                     ey [V/m]
0.0                     ez [V/m]

#TESTCASE
waveic

#WAVEIC
1                       seedB
1                       guideField
4                       waveMode
0.5                     frac

#TIMESTEPPING
T                       useFixedDt
0.002                   dt

#NORMALIZATION
1.0e5                   lNormSI   [m]
1.0e5                   uNormSI   [m/s]

#MONITOR
100                     dnReport

#SAVELOG
1                       dnSavePic
1                       dnSavePT

#PARTICLES
0                       nParticle per cell in X
0                       nParticle per cell in Y
0                       nParticle per cell in Z

#SAVEPLOT
1                       nPlot
z=0 fluid ascii pic     plotString
-1                      dn
2.0                     dt
1                       dx

#STOP
-1                      MaxIter
2.0                     TimeMax
