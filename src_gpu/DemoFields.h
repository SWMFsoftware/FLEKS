#ifndef DEMO_FIELDS_H
#define DEMO_FIELDS_H

#include "DemoConfig.h"
#include "DemoFieldKernels.h"
#include "DemoGrid.h"
#include "DemoLinearSolver.h"

class DemoFields {
public:
  DemoFields(const DemoGrid&, const DemoConfig&);
  amrex::MultiFab nodeE, nodeEth, nodeB, centerB, jHat, moments, rho, phi,
      centerMM, rhoOld;
  DemoNodeMMFab nodeMM;
  void clear_response();
  void finish_response();
  void refresh_node_b();
  void advance_b();
  bool solve_e();
  bool solve_div_e();
  void apply_e_operator(amrex::MultiFab&, amrex::MultiFab&);
  void apply_div_e_operator(amrex::MultiFab&, amrex::MultiFab&);
  // Prepared inputs must already have synchronized valid nodes and ghosts.
  // Output ghosts are optional when the next consumer uses valid data only.
  void curl_node_to_center(amrex::MultiFab&, amrex::MultiFab&,
                           bool inputPrepared = false, bool fillOutput = true);
  void curl_center_to_node(amrex::MultiFab&, amrex::MultiFab&,
                           bool inputPrepared = false, bool fillOutput = true);
  void divergence(amrex::MultiFab&, amrex::MultiFab&,
                  bool inputPrepared = false, bool fillOutput = true);
  int run_checks();
  amrex::Real field_energy();
  amrex::Real gauss_error();
  amrex::Real div_b_error();
  double lastResidual = 0, lastDivResidual = 0;
  int lastIterations = 0, lastDivIterations = 0;

private:
  const DemoGrid& grid;
  const DemoConfig& config;
  amrex::MultiFab rhs, curlCenter, curlNode, divE, scalarRhs, centerAverage,
      alternativeDiv;
  std::unique_ptr<amrex::iMultiFab> owner;
  DemoLinearSolver eSolver, divSolver;
};
#endif
