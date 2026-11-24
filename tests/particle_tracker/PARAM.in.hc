#DESCRIPTION
Test particle tracker: relativistic Higuera-Cary (2017) crossed E x B drift
and multi-flush extended TPOUTPUT benchmark.

Identical setup to PARAM.in, selecting TPPUSHER hc and runtime TPOUTPUT
with PBE vexb energy gamma and nRecordFlush = 8 over 24 steps (3 flushes).

#INITFROMSWMF
F                       initFromSWMF

#PIC
F                       usePIC

#GEOMETRY
-10.0                   xMin
 10.0                   xMax
-10.0                   yMin
 10.0                   yMax
-1.0                    zMin
 1.0                    zMax

#NCELL
32                      nCellX
8                       nCellY
1                       nCellZ

#MAXBLOCKSIZE
16                      maxBlockSizeX
8                       maxBlockSizeY
1                       maxBlockSizeZ

#PERIODICITY
T                       isPeriodicX
T                       isPeriodicY
T                       isPeriodicZ

#SOLVEEM
F                       solveEM

#HYBRIDPIC
F                       useHybridPIC

#NORMALIZATION
2.99792458e8            lNormSI
2.99792458e8            uNormSI

#PLASMA
2                       nSpecies
1.0                     mass
1.0                     charge
1.0                     mass
1.0                     charge

#UNIFORMSTATE
1.0                     rho [/cc]
239833.9664             ux [km/s]
0.0                     uy [km/s]
0.0                     uz [km/s]
0.0                     T [K]
1.0                     rho [/cc]
0.0                     ux [km/s]
0.0                     uy [km/s]
0.0                     uz [km/s]
0.0                     T [K]
0.0                     bx [T]
0.0                     by [T]
1.0e-7                  bz [T]
0.0                     ex [V/m]
23.98339664             ey [V/m]
0.0                     ez [V/m]

#TIMESTEPPING
T                       useFixedDt
0.25                    dt

#MONITOR
10                      dnReport

#SAVELOG
10                      dnSavePic
1                       dnSavePT

#PARTICLETRACKER
T                       usePT

#TPRELATIVISTIC
T                       isRelativistic

#TPPUSHER
hc                      pusherType

#TPPARTICLES
1                       npcelx
1                       npcely
1                       npcelz

#TPREGION
uniform                 region

#TPCELLINTERVAL
4                       nIntervalX
4                       nIntervalY
1                       nIntervalZ

#TPOUTPUT
PBE vexb energy gamma   outputVars
8                       nRecordFlush

#TPSAVE
0                       iSpecies
si                      IOUnit
1                       dnSave
0.0                     launchThreshold

#TPSAVE
1                       iSpecies
si                      IOUnit
1                       dnSave
0.0                     launchThreshold

#SAVEPLOT
0                       nPlot

#STOP
24                      MaxIter
-1.0                    TimeMax
