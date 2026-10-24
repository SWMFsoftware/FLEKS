#ifndef DEMO_LINEAR_SOLVER_H
#define DEMO_LINEAR_SOLVER_H

#include <functional>
#include <memory>
#include <vector>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

struct DemoSolveResult {
  bool converged = false;
  int iterations = 0;
  amrex::Real residual = 0, initialResidual = 0;
};

// The host manages Arnoldi/Hessenberg data, while all full-size vectors stay
// in persistent GPU MultiFabs. Supports both nodal E and cell-centered phi.
class DemoLinearSolver {
public:
  using Operator = std::function<void(amrex::MultiFab&, amrex::MultiFab&)>;
  DemoLinearSolver(const amrex::BoxArray&, const amrex::DistributionMapping&,
                   int ncomp, int restart, const amrex::Geometry&);
  DemoSolveResult solve(amrex::MultiFab& x, const amrex::MultiFab& rhs,
                        const Operator&, amrex::Real tolerance, int maxIter);
  amrex::Real dot(const amrex::MultiFab&, const amrex::MultiFab&) const;

private:
  int nComp, restartLength;
  std::unique_ptr<amrex::iMultiFab> owner;
  std::vector<std::unique_ptr<amrex::MultiFab> > basis;
  amrex::MultiFab residual, work;
};
#endif
