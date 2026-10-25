#ifndef DEMO_LINEAR_SOLVER_H
#define DEMO_LINEAR_SOLVER_H

#include <functional>
#include <memory>
#include <vector>

#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

struct DemoSolveResult {
  bool converged = false;
  int iterations = 0;
  amrex::Real residual = 0, initialResidual = 0;
};

// Counters always accumulate; section times only while timing is enabled
// (benchmark stage mode), because they synchronize the device at each
// boundary and add overhead. hostWaits
// counts the blocking host waits the solver itself issues (reductions
// returned to the host and synchronizing vector operations); waits inside
// the caller's operator are not included.
struct DemoSolverStats {
  long solves = 0, iterations = 0, hostWaits = 0, normFallbacks = 0;
  double operatorSeconds = 0, orthogonalizeSeconds = 0, updateSeconds = 0,
         residualSeconds = 0;
};

// The host manages Arnoldi/Hessenberg data, while all full-size vectors stay
// in persistent GPU MultiFabs. Supports both nodal E and cell-centered phi.
class DemoLinearSolver {
public:
  using Operator = std::function<void(amrex::MultiFab&, amrex::MultiFab&)>;
  DemoLinearSolver(const amrex::BoxArray&, const amrex::DistributionMapping&,
                   int ncomp, int restart, const amrex::Geometry&,
                   bool batched = false);
  // zeroInitialGuess: the caller guarantees x == 0, so the initial residual
  // is rhs and the first operator application is skipped.
  DemoSolveResult solve(amrex::MultiFab& x, const amrex::MultiFab& rhs,
                        const Operator&, amrex::Real tolerance, int maxIter,
                        bool zeroInitialGuess = false);
  amrex::Real dot(const amrex::MultiFab&, const amrex::MultiFab&) const;
  void set_timing(bool on) { timing = on; }
  const DemoSolverStats& stats() const { return statistics; }
  void reset_stats() { statistics = DemoSolverStats(); }

  // The helpers below launch device lambdas and are public only because nvcc
  // forbids extended lambdas in private or protected member functions.
  // A local box in the flattened cell index used by the fused reductions.
  struct BoxRange {
    int lo[3], length[3];
    amrex::Long offset;
  };

  template <class F> void section(double DemoSolverStats::*field, F&& work);
  amrex::Real waited_dot(const amrex::MultiFab&, const amrex::MultiFab&);
  // Batched (CGS2) helpers. Projections live in deviceProjection; slot
  // `pass` holds `count` values, plus work.work when withNorm is set.
  void project(int count, int pass, bool withNorm);
  void subtract_projection(int count, int pass);
  void scale_into(amrex::MultiFab& destination, const amrex::MultiFab& source,
                  amrex::Real scale);
  void residual_from(const amrex::MultiFab& rhs);
  void add_solution(amrex::MultiFab& x, const std::vector<amrex::Real>& y,
                    int used);
  void solve_batched_column(int j, std::vector<amrex::Real>& column,
                            amrex::Real& next);


private:
  int nComp, restartLength;
  bool batched, timing = false;
  DemoSolverStats statistics;
  // CGS2 stores the basis once, component-packed with one ghost layer;
  // basis[j] aliases components [j*nComp, (j+1)*nComp).
  std::unique_ptr<amrex::MultiFab> packedBasis;
  amrex::Gpu::DeviceVector<amrex::Real> deviceProjection, deviceCoefficients;
  amrex::Gpu::PinnedVector<amrex::Real> hostProjection, hostCoefficients;
  amrex::Gpu::DeviceVector<BoxRange> boxRanges;
  amrex::Long localCells = 0;
  std::unique_ptr<amrex::iMultiFab> owner;
  std::vector<std::unique_ptr<amrex::MultiFab> > basis;
  amrex::MultiFab residual, work;
};
#endif
