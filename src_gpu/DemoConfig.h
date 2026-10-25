#ifndef DEMO_CONFIG_H_
#define DEMO_CONFIG_H_

#include <string>

#include <AMReX_REAL.H>

struct DemoConfig {
  int nCell = 16;
  int maxGridSize = 8;
  int nIter = 10;
  int nStep = 20;
  int particlesPerCell = 8;
  int seed = 1;
  int diagnosticInterval = 1;
  int gmresRestart = 20;
  int gmresMaxIter = 200;
  int divEMaxIter = 20;
  int nDivECorrection = 3;
  int productionGhost = 3;
  int sortInterval = 0;
  int warmupSteps = 1;
  int cellChunkSize = 0;
  amrex::Real dt = 0.01;
  amrex::Real theta = 0.51;
  amrex::Real massRatio = 25;
  amrex::Real perturbation = 0.001;
  amrex::Real b0 = 0;
  amrex::Real solverTolerance = 1.e-9;
  amrex::Real divETolerance = 0.01;
  amrex::Real coefDiff = 0.1;
  bool correctDivE = true;
  bool fusedPush = false;
  bool bufferedResponse = false;
  bool clustered = false;
  bool microbenchmarks = true;
  bool pairSymmetry = false;
  bool foldCharge = false;
  // Twice classical Gram-Schmidt with one host wait per Arnoldi step (see
  // DemoLinearSolver.cpp); 0 selects the reference twice modified
  // Gram-Schmidt.
  bool batchedGmres = true;
  // Time solver sections (operator, orthogonalization, ...). Each section
  // boundary synchronizes, so this inflates the stage-mode solve times.
  bool solverSections = false;
  // div(E) correction solve: gmres (production), jacobi_gmres (right Jacobi
  // preconditioning), or pcg (Jacobi-preconditioned conjugate gradients).
  // The last two change the partial potential at the iteration limit.
  std::string divEMethod = "gmres";
  std::string benchmarkTiming = "stages";
  std::string mode = "fields";
  // Mass-matrix assembly: atomic (per particle), cell (block-reduced per
  // cell, GPU), or buffered (serial per-cell thread buffers).
  std::string response = "atomic";
  // div(E) charge/matrix, final charge, and moment deposits: atomic (per
  // particle) or cell (block-reduced per cell, GPU).
  std::string deposit = "atomic";
  std::string initialCondition = "uniform";
  std::string referencePrefix;

  DemoConfig();
};

#endif
