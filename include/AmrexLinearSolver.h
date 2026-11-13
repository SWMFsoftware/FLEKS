#ifndef AMREX_LINEAR_SOLVER_H_
#define AMREX_LINEAR_SOLVER_H_

#include <cstring>

#include <AMReX_GMRES.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include "GpuBackend.h"
#include "LinearSolver.h"

class Pic;

enum class ESolverBackend { LegacyGMRES, AmrexGMRES };

enum class EPrecondType {
  Auto, // Adaptive selection (None for <5 iters or JFNK, BlockJacobi otherwise)
  None, // P = I (unpreconditioned)
  Jacobi,      // Scalar diagonal scaling (3 components per node)
  BlockJacobi, // Point-block 3x3 exact inverse per node
  Smoother     // Multi-step weighted block-Jacobi smoother using full matvec
};

namespace amrex {

/**
 * @brief Shared 0-ghost MultiFab vector operations for amrex::GMRES wrappers.
 *
 * All Krylov vectors are 0-ghost MultiFabs whose inactive DOFs are maintained
 * at zero, enabling contiguous SIMD inner products and norms without per-element
 * mask branching.
 */
class AmrexMultiFabOps {
public:
  using RT = Real;
  using V = MultiFab;

  void assign(V& lhs, V const& rhs) {
    if constexpr (fleksUseGpu) {
      MultiFab::Copy(lhs, rhs, 0, 0, lhs.nComp(), 0);
    } else {
      const int nFabs = lhs.local_size();
      for (int li = 0; li < nFabs; ++li) {
        Real* AMREX_RESTRICT pL = lhs.atLocalIdx(li).dataPtr();
        const Real* AMREX_RESTRICT pR = rhs.atLocalIdx(li).dataPtr();
        const Long n = lhs.atLocalIdx(li).size();
        std::memcpy(pL, pR, static_cast<std::size_t>(n) * sizeof(Real));
      }
    }
  }

  RT dotProduct(V const& v1, V const& v2) {
    if constexpr (fleksUseGpu) {
      return MultiFab::Dot(v1, 0, v2, 0, v1.nComp(), 0);
    } else {
      RT sm = RT(0.0);
      const int nFabs = v1.local_size();
      for (int li = 0; li < nFabs; ++li) {
        const Real* AMREX_RESTRICT p1 = v1.atLocalIdx(li).dataPtr();
        const Real* AMREX_RESTRICT p2 = v2.atLocalIdx(li).dataPtr();
        const Long n = v1.atLocalIdx(li).size();
        RT localSm = RT(0.0);
#pragma omp simd reduction(+ : localSm)
        for (Long i = 0; i < n; ++i) {
          localSm += p1[i] * p2[i];
        }
        sm += localSm;
      }
      ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
      return sm;
    }
  }

  void increment(V& lhs, V const& rhs, RT a) {
    if constexpr (fleksUseGpu) {
      MultiFab::Saxpy(lhs, a, rhs, 0, 0, lhs.nComp(), 0);
    } else {
      const int nFabs = lhs.local_size();
      for (int li = 0; li < nFabs; ++li) {
        Real* AMREX_RESTRICT pL = lhs.atLocalIdx(li).dataPtr();
        const Real* AMREX_RESTRICT pR = rhs.atLocalIdx(li).dataPtr();
        const Long n = lhs.atLocalIdx(li).size();
#pragma omp simd
        for (Long i = 0; i < n; ++i) {
          pL[i] += a * pR[i];
        }
      }
    }
  }

  void linComb(V& lhs, RT a, V const& rhs_a, RT b, V const& rhs_b) {
    if constexpr (fleksUseGpu) {
      MultiFab::LinComb(lhs, a, rhs_a, 0, b, rhs_b, 0, 0, lhs.nComp(), 0);
    } else {
      const int nFabs = lhs.local_size();
      for (int li = 0; li < nFabs; ++li) {
        Real* AMREX_RESTRICT pL = lhs.atLocalIdx(li).dataPtr();
        const Real* AMREX_RESTRICT pA = rhs_a.atLocalIdx(li).dataPtr();
        const Real* AMREX_RESTRICT pB = rhs_b.atLocalIdx(li).dataPtr();
        const Long n = lhs.atLocalIdx(li).size();
#pragma omp simd
        for (Long i = 0; i < n; ++i) {
          pL[i] = a * pA[i] + b * pB[i];
        }
      }
    }
  }

  RT norm2(V const& v) {
    if constexpr (fleksUseGpu) {
      return std::sqrt(MultiFab::Dot(v, 0, v, 0, v.nComp(), 0));
    } else {
      RT sm = RT(0.0);
      const int nFabs = v.local_size();
      for (int li = 0; li < nFabs; ++li) {
        const Real* AMREX_RESTRICT p = v.atLocalIdx(li).dataPtr();
        const Long n = v.atLocalIdx(li).size();
        RT localSm = RT(0.0);
#pragma omp simd reduction(+ : localSm)
        for (Long i = 0; i < n; ++i) {
          localSm += p[i] * p[i];
        }
        sm += localSm;
      }
      ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
      return std::sqrt(sm);
    }
  }

  void scale(V& v, RT fac) {
    if constexpr (fleksUseGpu) {
      v.mult(fac, 0, v.nComp(), 0);
    } else {
      const int nFabs = v.local_size();
      for (int li = 0; li < nFabs; ++li) {
        Real* AMREX_RESTRICT p = v.atLocalIdx(li).dataPtr();
        const Long n = v.atLocalIdx(li).size();
#pragma omp simd
        for (Long i = 0; i < n; ++i) {
          p[i] *= fac;
        }
      }
    }
  }

  void setToZero(V& v) {
    if constexpr (fleksUseGpu) {
      v.setVal(RT(0.0));
    } else {
      const int nFabs = v.local_size();
      for (int li = 0; li < nFabs; ++li) {
        Real* AMREX_RESTRICT p = v.atLocalIdx(li).dataPtr();
        const Long n = v.atLocalIdx(li).size();
        std::memset(p, 0, static_cast<std::size_t>(n) * sizeof(Real));
      }
    }
  }
};

/**
 * @brief Linear operator and preconditioner wrapper for amrex::GMRES
 * operating directly on nodal 3-component MultiFabs for the FLEKS implicit E
 * solver.
 */
class AmrexLinearSolver : public AmrexMultiFabOps {
public:
  AmrexLinearSolver() = default;
  AmrexLinearSolver(Pic* pic, int iLev) : m_pic(pic), m_iLev(iLev) {}

  void define(Pic* pic, int iLev, EPrecondType precondType = EPrecondType::Auto,
              int precondIters = 2, Real smootherWeight = Real(2.0 / 3.0),
              FieldSolverMode mode = FieldSolverMode::GMRES);

  void set_newton_base(const MultiFab& baseE, const MultiFab& baseMatvec);

  void mask_inactive_nodes(MultiFab& mf) const;

  void solve(MultiFab& sol, const MultiFab& rhs, Real tolRel, Real tolAbs,
             int maxIters, int restartLen, int verbose, int& outIters,
             Real& outResNorm);

  // Required by amrex::GMRES
  void apply(V& Ax, V const& x);

  V makeVecRHS();
  V makeVecLHS();

  void precond(V& lhs, V const& rhs);

  const iMultiFab& dot_mask() const { return m_dotMask; }

private:
  void ensure_mask();
  void apply_linear_matvec(V& Ax, V const& x);
  void apply_point_precond(V& dst, V const& src, Real scaleFactor = 1.0) const;

  Pic* m_pic = nullptr;
  int m_iLev = -1;
  EPrecondType m_precondType = EPrecondType::None;
  int m_lastIters = 0;
  int m_precondIters = 2;
  Real m_smootherWeight = Real(2.0 / 3.0);
  FieldSolverMode m_mode = FieldSolverMode::GMRES;

  // Persistent GMRES instance to avoid re-initializing Hessenberg tables
  GMRES<MultiFab, AmrexLinearSolver> m_gmres;
  bool m_gmresDefined = false;

  // Active-DOF mask: 1 on owner nodes that belong to the linear system, 0
  // on shared duplicate boundary/periodic nodes and excluded body nodes.
  iMultiFab m_dotMask;

  // Precomputed inverse diagonal:
  // - 3 components for EPrecondType::Jacobi
  // - 9 components (row-major 3x3 inverse) for BlockJacobi and Smoother
  MultiFab m_invDiagMF;

  // Scratch MultiFabs for multi-step smoother (0-ghost)
  MultiFab m_smoothAx;
  MultiFab m_smoothRes;

  // 0-ghost solution and RHS buffers for GMRES::solve
  MultiFab m_sol0;
  MultiFab m_rhs0;

  // Scratch MultiFabs for Jacobian-Free Newton-Krylov (0-ghost)
  MultiFab m_nkBase;
  MultiFab m_nkBaseMatvec;
  Real m_nkNormBase = 0.0;
};

/**
 * @brief Linear operator and preconditioner wrapper for amrex::GMRES
 * operating directly on cell-centered 1-component MultiFabs for FLEKS div(E)
 * cleaning.
 */
class AmrexDivESolver : public AmrexMultiFabOps {
public:
  AmrexDivESolver() = default;
  AmrexDivESolver(Pic* pic, int iLev) : m_pic(pic), m_iLev(iLev) {}

  void define(Pic* pic, int iLev, EPrecondType precondType = EPrecondType::Auto);

  void mask_inactive_cells(MultiFab& mf) const;

  void solve(MultiFab& sol, const MultiFab& rhs, Real tolRel, Real tolAbs,
             int maxIters, int restartLen, int verbose, int& outIters,
             Real& outResNorm);

  // Required by amrex::GMRES
  void apply(V& Ax, V const& x);

  V makeVecRHS();
  V makeVecLHS();

  void precond(V& lhs, V const& rhs);

  const iMultiFab& dot_mask() const { return m_dotMask; }

private:
  void ensure_mask();

  Pic* m_pic = nullptr;
  int m_iLev = -1;
  EPrecondType m_precondType = EPrecondType::None;
  int m_lastIters = 0;

  GMRES<MultiFab, AmrexDivESolver> m_gmres;
  bool m_gmresDefined = false;

  // Active-DOF mask: 1 on active cells, 0 on skipped AMR/body cells.
  iMultiFab m_dotMask;

  // Precomputed scalar inverse diagonal (1 component) for Jacobi preconditioning.
  MultiFab m_invDiagMF;

  // 0-ghost solution and RHS buffers for GMRES::solve
  MultiFab m_sol0;
  MultiFab m_rhs0;
};

} // namespace amrex

#endif
