#DESCRIPTION
BC_DISCONTINUITY (RD): Rotational discontinuity propagating along +x.

Normal field Bx = 3 nT and a 180-degree circular rotation of the tangential
field By: +4 nT -> -4 nT at constant |B| = 5 nT and exact Walen velocity jump
du = +(B2 - B1)/sqrt(mu0 * rho_m).

#INITFROMSWMF
F                       initFromSWMF

#NORMALIZATION
1.0e5                   lNorm [m] (100 km)
1.0e6                   uNorm [m/s]

#GEOMETRY
0.0                     xMin
10.0                    xMax
0.0                     yMin
1.0                     yMax
0.0                     zMin
0.5                     zMax

#NCELL
40                      nCellX
2                       nCellY
1                       nCellZ

#MAXBLOCKSIZE
40                      maxBlockSizeX
2                       maxBlockSizeY
1                       maxBlockSizeZ

#PERIODICITY
F                       isPeriodicX
T                       isPeriodicY
T                       isPeriodicZ

#FIELDBOXBOUNDARY
inflow                  fieldBoxBoundaryLoX
outflow                 fieldBoxBoundaryHiX
periodic                fieldBoxBoundaryLoY
periodic                fieldBoxBoundaryHiY
periodic                fieldBoxBoundaryLoZ
periodic                fieldBoxBoundaryHiZ

#PARTICLEBOXBOUNDARY
inflow                  particleBoxBoundaryLoX
outflow                 particleBoxBoundaryHiX
periodic                particleBoxBoundaryLoY
periodic                particleBoxBoundaryHiY
periodic                particleBoxBoundaryLoZ
periodic                particleBoxBoundaryHiZ

#PLASMA
1                       nS
1.0                     mass
1.0                     charge

#PARTICLES
8                       npcelx
2                       npcely
1                       npcelz

#RANDOMPARTICLESLOCATION
F                       isParticleLocationRandom

#UNIFORMSTATE
5.0                     rho [/cc]
150.0                   ux [km/s]
39.01855                uy [km/s]
0.0                     uz [km/s]
50000.0                 T [K]
3.0e-9                  bx [T]
4.0e-9                  by [T]
0.0                     bz [T]
0.0                     ex [V/m]
0.0                     ey [V/m]
-4.8294435e-4           ez [V/m]

#DISCONTINUITY
rd                      typeDiscontinuity
1.0                     normalX
0.0                     normalY
0.0                     normalZ
0.8                     width [code units]
0.0                     tStart [s]
5.0                     rho [/cc]
150.0                   ux [km/s]
-39.01855               uy [km/s]
0.0                     uz [km/s]
50000.0                 T [K]
3.0e-9                  bx [T]
-4.0e-9                 by [T]
0.0                     bz [T]
-1.0                    vn [km/s]

#FIELDSOLVER
hybrid                  fieldSolverMode

#HALLTERM
F                       useHallTerm

#ELECTRONTEMPERATURE
0.0                     electronTemperature [eV]
1.0                     electronGamma
1.0                     electronDensity0

#FIELDSUBCYCLE
2                       nSubcycle
F                       isAutoSubcycle

#TIMESTEPPING
T                       useFixedDt
0.05                    dt [s]

#SAVELOG
10                      dnSavePic
-1                      dnSavePT

#STOP
120                     nStep
6.0                     tMax

#SAVEPLOT
1                       nPlot
z=0 fluid ascii si      plotString
120                     dnSavePlot
-1                      dtSavePlot
-1                      dx
