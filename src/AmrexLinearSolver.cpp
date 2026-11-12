#include "AmrexLinearSolver.h"
#include "Bit.h"
#include "Pic.h"

namespace amrex {

void AmrexLinearSolver::ensure_mask() {
  const BoxArray& ba = m_pic->node_box_array(m_iLev);
  const DistributionMapping& dm = m_pic->DistributionMap(m_iLev);
  if (!m_dotMask.empty() && m_dotMask.boxArray() == ba &&
      m_dotMask.DistributionMap() == dm) {
    return;
  }

  m_dotMask.define(ba, dm, 1, 0);
  const bool skipAllBody = m_pic->is_body_linetied();
  const bool skipInterior = m_pic->is_body_conducting();
  const auto& statusMF = m_pic->node_status(m_iLev);

  for (MFIter mfi(m_dotMask, doTiling); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const Array4<int>& mskArr = m_dotMask[mfi].array();
    const Array4<int const>& nodeArr = statusMF[mfi].const_array();
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int st = nodeArr(i, j, k);
      const bool isActive = bit::is_owner(st) &&
                            !(skipAllBody && bit::is_body(st)) &&
                            !(skipInterior && bit::is_body_interior(st));
      mskArr(i, j, k) = isActive ? 1 : 0;
    });
  }
}

void AmrexLinearSolver::mask_inactive_nodes(MultiFab& mf) const {
#ifdef AMREX_USE_GPU
  for (MFIter mfi(mf, doTiling); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const Array4<Real>& arr = mf[mfi].array();
    const Array4<int const>& msk = m_dotMask[mfi].const_array();
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const Real m = static_cast<Real>(msk(i, j, k));
      arr(i, j, k, 0) *= m;
      arr(i, j, k, 1) *= m;
      arr(i, j, k, 2) *= m;
    });
  }
#else
  const int nFabs = mf.local_size();
  for (int li = 0; li < nFabs; ++li) {
    Real* AMREX_RESTRICT p0 = mf.atLocalIdx(li).dataPtr(0);
    Real* AMREX_RESTRICT p1 = mf.atLocalIdx(li).dataPtr(1);
    Real* AMREX_RESTRICT p2 = mf.atLocalIdx(li).dataPtr(2);
    const int* AMREX_RESTRICT pM = m_dotMask.atLocalIdx(li).dataPtr();
    const Long nPts = mf.atLocalIdx(li).box().numPts();
    for (Long i = 0; i < nPts; ++i) {
      const Real m = static_cast<Real>(pM[i]);
      p0[i] *= m;
      p1[i] *= m;
      p2[i] *= m;
    }
  }
#endif
}

void AmrexLinearSolver::define(Pic* pic, int iLev, EPrecondType precondType,
                               int precondIters, Real smootherWeight,
                               FieldSolverMode mode) {
  m_pic = pic;
  m_iLev = iLev;
  m_precondIters = std::max(1, precondIters);
  m_smootherWeight = smootherWeight;
  m_mode = mode;

  if (precondType == EPrecondType::Auto) {
    if (m_mode == FieldSolverMode::NewtonKrylov || m_lastIters < 5) {
      m_precondType = EPrecondType::None;
    } else {
      m_precondType = EPrecondType::BlockJacobi;
    }
  } else {
    m_precondType = precondType;
  }

  ensure_mask();

  if (m_precondType != EPrecondType::None) {
    m_pic->build_E_preconditioner(m_invDiagMF, m_iLev, m_precondType);
  }

  if (m_precondType == EPrecondType::Smoother && m_precondIters > 1) {
    const BoxArray& ba = m_pic->node_box_array(m_iLev);
    const DistributionMapping& dm = m_pic->DistributionMap(m_iLev);
    if (m_smoothAx.empty() || m_smoothAx.boxArray() != ba ||
        m_smoothAx.DistributionMap() != dm) {
      m_smoothAx.define(ba, dm, 3, 0);
      m_smoothAx.setVal(0.0);
      m_smoothRes.define(ba, dm, 3, 0);
      m_smoothRes.setVal(0.0);
    }
  }

  if (!m_gmresDefined) {
    m_gmres.define(*this);
    m_gmresDefined = true;
  }
}

void AmrexLinearSolver::set_newton_base(const MultiFab& baseE,
                                        const MultiFab& baseMatvec) {
  const BoxArray& ba = m_pic->node_box_array(m_iLev);
  const DistributionMapping& dm = m_pic->DistributionMap(m_iLev);
  if (m_nkBase.empty() || m_nkBase.boxArray() != ba ||
      m_nkBase.DistributionMap() != dm) {
    m_nkBase.define(ba, dm, 3, 0);
    m_nkBaseMatvec.define(ba, dm, 3, 0);
  }
  MultiFab::Copy(m_nkBase, baseE, 0, 0, 3, 0);
  mask_inactive_nodes(m_nkBase);
  m_nkNormBase = norm2(m_nkBase);

  MultiFab::Copy(m_nkBaseMatvec, baseMatvec, 0, 0, 3, 0);
  mask_inactive_nodes(m_nkBaseMatvec);
}

void AmrexLinearSolver::solve(MultiFab& sol, const MultiFab& rhs, Real tolRel,
                              Real tolAbs, int maxIters, int restartLen,
                              int verbose, int& outIters, Real& outResNorm) {
  const BoxArray& ba = m_pic->node_box_array(m_iLev);
  const DistributionMapping& dm = m_pic->DistributionMap(m_iLev);
  if (m_sol0.empty() || m_sol0.boxArray() != ba ||
      m_sol0.DistributionMap() != dm) {
    m_sol0.define(ba, dm, 3, 0);
    m_rhs0.define(ba, dm, 3, 0);
  }
  setToZero(m_sol0);
  MultiFab::Copy(m_rhs0, rhs, 0, 0, 3, 0);
  mask_inactive_nodes(m_rhs0);

  m_gmres.setRestartLength(restartLen);
  m_gmres.setMaxIters(maxIters);
  m_gmres.setVerbose(verbose);
  m_gmres.solve(m_sol0, m_rhs0, tolRel, tolAbs);
  outIters = m_gmres.getNumIters();
  outResNorm = m_gmres.getResidualNorm();
  m_lastIters = outIters;

  sol.setVal(0.0);
  MultiFab::Copy(sol, m_sol0, 0, 0, 3, 0);
}

void AmrexLinearSolver::apply_linear_matvec(V& Ax, V const& x) {
  m_pic->update_E_matvec(x, Ax, m_iLev, true);
  mask_inactive_nodes(Ax);
}

void AmrexLinearSolver::apply(V& Ax, V const& x) {
  if (m_mode == FieldSolverMode::NewtonKrylov) {
    const Real normDir = norm2(x);
    const Real eps =
        fleks_jfnk::finite_difference_epsilon(m_nkNormBase, normDir);
    if (eps == 0.0) {
      setToZero(Ax);
      return;
    }
    // Use Ax temporarily to hold perturbed state base + eps * x (0-ghost,
    // masked on inactive nodes), then evaluate nonlinear F(base + eps * x).
    linComb(Ax, 1.0, m_nkBase, eps, x);
    m_pic->update_E_matvec(Ax, Ax, m_iLev, false);
    mask_inactive_nodes(Ax);
    const Real invEps = 1.0 / eps;
    linComb(Ax, invEps, Ax, -invEps, m_nkBaseMatvec);
  } else {
    apply_linear_matvec(Ax, x);
  }
}

MultiFab AmrexLinearSolver::makeVecRHS() {
  MultiFab mf(m_pic->node_box_array(m_iLev), m_pic->DistributionMap(m_iLev), 3,
              0);
  setToZero(mf);
  return mf;
}

MultiFab AmrexLinearSolver::makeVecLHS() {
  MultiFab mf(m_pic->node_box_array(m_iLev), m_pic->DistributionMap(m_iLev), 3,
              0);
  setToZero(mf);
  return mf;
}

void AmrexLinearSolver::apply_point_precond(V& dst, V const& src,
                                            Real scaleFactor) const {
#ifdef AMREX_USE_GPU
  if (m_precondType == EPrecondType::Jacobi) {
    for (MFIter mfi(dst, doTiling); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      const Array4<Real>& d = dst[mfi].array();
      const Array4<Real const>& s = src[mfi].const_array();
      const Array4<Real const>& inv = m_invDiagMF[mfi].const_array();
      ParallelFor(
          bx, 3, [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
            d(i, j, k, c) = scaleFactor * inv(i, j, k, c) * s(i, j, k, c);
          });
    }
  } else {
    for (MFIter mfi(dst, doTiling); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      const Array4<Real>& d = dst[mfi].array();
      const Array4<Real const>& s = src[mfi].const_array();
      const Array4<Real const>& inv = m_invDiagMF[mfi].const_array();
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real sx = s(i, j, k, 0);
        const Real sy = s(i, j, k, 1);
        const Real sz = s(i, j, k, 2);
        d(i, j, k, 0) =
            scaleFactor * (inv(i, j, k, 0) * sx + inv(i, j, k, 1) * sy +
                           inv(i, j, k, 2) * sz);
        d(i, j, k, 1) =
            scaleFactor * (inv(i, j, k, 3) * sx + inv(i, j, k, 4) * sy +
                           inv(i, j, k, 5) * sz);
        d(i, j, k, 2) =
            scaleFactor * (inv(i, j, k, 6) * sx + inv(i, j, k, 7) * sy +
                           inv(i, j, k, 8) * sz);
      });
    }
  }
#else
  const int nFabs = dst.local_size();
  if (m_precondType == EPrecondType::Jacobi) {
    for (int li = 0; li < nFabs; ++li) {
      Real* AMREX_RESTRICT pD = dst.atLocalIdx(li).dataPtr();
      const Real* AMREX_RESTRICT pS = src.atLocalIdx(li).dataPtr();
      const Real* AMREX_RESTRICT pInv = m_invDiagMF.atLocalIdx(li).dataPtr();
      const Long n = dst.atLocalIdx(li).size();
      for (Long i = 0; i < n; ++i) {
        pD[i] = scaleFactor * pInv[i] * pS[i];
      }
    }
  } else {
    for (int li = 0; li < nFabs; ++li) {
      Real* AMREX_RESTRICT d0 = dst.atLocalIdx(li).dataPtr(0);
      Real* AMREX_RESTRICT d1 = dst.atLocalIdx(li).dataPtr(1);
      Real* AMREX_RESTRICT d2 = dst.atLocalIdx(li).dataPtr(2);
      const Real* AMREX_RESTRICT s0 = src.atLocalIdx(li).dataPtr(0);
      const Real* AMREX_RESTRICT s1 = src.atLocalIdx(li).dataPtr(1);
      const Real* AMREX_RESTRICT s2 = src.atLocalIdx(li).dataPtr(2);
      const Real* AMREX_RESTRICT inv00 = m_invDiagMF.atLocalIdx(li).dataPtr(0);
      const Real* AMREX_RESTRICT inv01 = m_invDiagMF.atLocalIdx(li).dataPtr(1);
      const Real* AMREX_RESTRICT inv02 = m_invDiagMF.atLocalIdx(li).dataPtr(2);
      const Real* AMREX_RESTRICT inv10 = m_invDiagMF.atLocalIdx(li).dataPtr(3);
      const Real* AMREX_RESTRICT inv11 = m_invDiagMF.atLocalIdx(li).dataPtr(4);
      const Real* AMREX_RESTRICT inv12 = m_invDiagMF.atLocalIdx(li).dataPtr(5);
      const Real* AMREX_RESTRICT inv20 = m_invDiagMF.atLocalIdx(li).dataPtr(6);
      const Real* AMREX_RESTRICT inv21 = m_invDiagMF.atLocalIdx(li).dataPtr(7);
      const Real* AMREX_RESTRICT inv22 = m_invDiagMF.atLocalIdx(li).dataPtr(8);
      const Long nPts = dst.atLocalIdx(li).box().numPts();
      for (Long i = 0; i < nPts; ++i) {
        const Real sx = s0[i];
        const Real sy = s1[i];
        const Real sz = s2[i];
        d0[i] = scaleFactor * (inv00[i] * sx + inv01[i] * sy + inv02[i] * sz);
        d1[i] = scaleFactor * (inv10[i] * sx + inv11[i] * sy + inv12[i] * sz);
        d2[i] = scaleFactor * (inv20[i] * sx + inv21[i] * sy + inv22[i] * sz);
      }
    }
  }
#endif

  if (m_pic->is_body_conducting()) {
    m_pic->project_body_E(dst, m_iLev);
  }
}

void AmrexLinearSolver::precond(V& lhs, V const& rhs) {
  if (m_precondType == EPrecondType::None) {
    assign(lhs, rhs);
    return;
  }

  if (m_precondType == EPrecondType::Jacobi ||
      m_precondType == EPrecondType::BlockJacobi || m_precondIters <= 1) {
    apply_point_precond(lhs, rhs, 1.0);
    return;
  }

  // Multi-step weighted block-Jacobi smoother:
  //   z^(1) = omega * D^-1 * rhs
  //   z^(m+1) = z^(m) + omega * D^-1 * (rhs - A * z^(m))
  apply_point_precond(lhs, rhs, m_smootherWeight);

  const Real omega = m_smootherWeight;
#ifndef AMREX_USE_GPU
  const int nFabs = lhs.local_size();
#endif
  for (int iter = 1; iter < m_precondIters; ++iter) {
    apply_linear_matvec(m_smoothAx, lhs);
    linComb(m_smoothRes, 1.0, rhs, -1.0, m_smoothAx);

#ifdef AMREX_USE_GPU
    for (MFIter mfi(lhs, doTiling); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      const Array4<Real>& d = lhs[mfi].array();
      const Array4<Real const>& r = m_smoothRes[mfi].const_array();
      const Array4<Real const>& inv = m_invDiagMF[mfi].const_array();
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const Real rx = r(i, j, k, 0);
        const Real ry = r(i, j, k, 1);
        const Real rz = r(i, j, k, 2);
        d(i, j, k, 0) += omega * (inv(i, j, k, 0) * rx + inv(i, j, k, 1) * ry +
                                  inv(i, j, k, 2) * rz);
        d(i, j, k, 1) += omega * (inv(i, j, k, 3) * rx + inv(i, j, k, 4) * ry +
                                  inv(i, j, k, 5) * rz);
        d(i, j, k, 2) += omega * (inv(i, j, k, 6) * rx + inv(i, j, k, 7) * ry +
                                  inv(i, j, k, 8) * rz);
      });
    }
#else
    for (int li = 0; li < nFabs; ++li) {
      Real* AMREX_RESTRICT d0 = lhs.atLocalIdx(li).dataPtr(0);
      Real* AMREX_RESTRICT d1 = lhs.atLocalIdx(li).dataPtr(1);
      Real* AMREX_RESTRICT d2 = lhs.atLocalIdx(li).dataPtr(2);
      const Real* AMREX_RESTRICT r0 = m_smoothRes.atLocalIdx(li).dataPtr(0);
      const Real* AMREX_RESTRICT r1 = m_smoothRes.atLocalIdx(li).dataPtr(1);
      const Real* AMREX_RESTRICT r2 = m_smoothRes.atLocalIdx(li).dataPtr(2);
      const Real* AMREX_RESTRICT inv00 = m_invDiagMF.atLocalIdx(li).dataPtr(0);
      const Real* AMREX_RESTRICT inv01 = m_invDiagMF.atLocalIdx(li).dataPtr(1);
      const Real* AMREX_RESTRICT inv02 = m_invDiagMF.atLocalIdx(li).dataPtr(2);
      const Real* AMREX_RESTRICT inv10 = m_invDiagMF.atLocalIdx(li).dataPtr(3);
      const Real* AMREX_RESTRICT inv11 = m_invDiagMF.atLocalIdx(li).dataPtr(4);
      const Real* AMREX_RESTRICT inv12 = m_invDiagMF.atLocalIdx(li).dataPtr(5);
      const Real* AMREX_RESTRICT inv20 = m_invDiagMF.atLocalIdx(li).dataPtr(6);
      const Real* AMREX_RESTRICT inv21 = m_invDiagMF.atLocalIdx(li).dataPtr(7);
      const Real* AMREX_RESTRICT inv22 = m_invDiagMF.atLocalIdx(li).dataPtr(8);
      const Long nPts = lhs.atLocalIdx(li).box().numPts();
      for (Long i = 0; i < nPts; ++i) {
        const Real rx = r0[i];
        const Real ry = r1[i];
        const Real rz = r2[i];
        d0[i] += omega * (inv00[i] * rx + inv01[i] * ry + inv02[i] * rz);
        d1[i] += omega * (inv10[i] * rx + inv11[i] * ry + inv12[i] * rz);
        d2[i] += omega * (inv20[i] * rx + inv21[i] * ry + inv22[i] * rz);
      }
    }
#endif

    if (m_pic->is_body_conducting()) {
      m_pic->project_body_E(lhs, m_iLev);
    }
  }
}

} // namespace amrex
