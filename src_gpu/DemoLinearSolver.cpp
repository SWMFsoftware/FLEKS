#include <algorithm>
#include <cmath>
#include <limits>

#include "DemoLinearSolver.h"
using namespace amrex;

DemoLinearSolver::DemoLinearSolver(const BoxArray& ba,
                                   const DistributionMapping& dm, int nc,
                                   int restart, const Geometry& geom)
    : nComp(nc),
      restartLength(restart),
      residual(ba, dm, nc, 1),
      work(ba, dm, nc, 1) {
  owner = residual.OwnerMask(geom.periodicity());
  for (int n = 0; n <= restart; ++n)
    basis.push_back(std::make_unique<MultiFab>(ba, dm, nc, 1));
}
Real DemoLinearSolver::dot(const MultiFab& a, const MultiFab& b) const {
  return MultiFab::Dot(*owner, a, 0, b, 0, nComp, 0);
}
DemoSolveResult DemoLinearSolver::solve(MultiFab& x, const MultiFab& rhs,
                                        const Operator& apply, Real tolerance,
                                        int maxIter) {
  DemoSolveResult result;
  const Real rhsNorm = std::sqrt(dot(rhs, rhs));
  const Real target = tolerance * std::max(rhsNorm, Real(1e-14));
  auto true_residual = [&]() {
    apply(x, work);
    MultiFab::LinComb(residual, 1, rhs, 0, -1, work, 0, 0, nComp, 0);
    return std::sqrt(dot(residual, residual));
  };
  result.initialResidual = true_residual();
  result.residual = result.initialResidual;
  while (result.iterations < maxIter && result.residual > target) {
    const int m = std::min(restartLength, maxIter - result.iterations);
    std::vector<Real> h((m + 1) * m, 0), cs(m), sn(m), g(m + 1, 0), y(m, 0);
    auto H = [&](int i, int j) -> Real& { return h[j * (m + 1) + i]; };
    MultiFab::Copy(*basis[0], residual, 0, 0, nComp, 0);
    basis[0]->mult(1 / result.residual, 0, nComp, 0);
    g[0] = result.residual;
    int used = 0;
    for (int j = 0; j < m; ++j) {
      apply(*basis[j], work);
      // Twice modified Gram-Schmidt prevents false convergence for ill-scaled
      // response matrices and costs only scalar global reductions.
      for (int pass = 0; pass < 2; ++pass)
        for (int i = 0; i <= j; ++i) {
          const Real projection = dot(work, *basis[i]);
          H(i, j) += projection;
          MultiFab::Saxpy(work, -projection, *basis[i], 0, 0, nComp, 0);
        }
      H(j + 1, j) = std::sqrt(dot(work, work));
      const Real next = H(j + 1, j);
      if (next > 0) {
        MultiFab::Copy(*basis[j + 1], work, 0, 0, nComp, 0);
        basis[j + 1]->mult(1 / next, 0, nComp, 0);
      }
      for (int i = 0; i < j; ++i) {
        const Real t = cs[i] * H(i, j) + sn[i] * H(i + 1, j);
        H(i + 1, j) = -sn[i] * H(i, j) + cs[i] * H(i + 1, j);
        H(i, j) = t;
      }
      const Real norm = std::hypot(H(j, j), H(j + 1, j));
      cs[j] = norm > 0 ? H(j, j) / norm : 1;
      sn[j] = norm > 0 ? H(j + 1, j) / norm : 0;
      H(j, j) = norm;
      H(j + 1, j) = 0;
      g[j + 1] = -sn[j] * g[j];
      g[j] *= cs[j];
      used = j + 1;
      ++result.iterations;
      if (std::abs(g[j + 1]) <= target ||
          next <= std::numeric_limits<Real>::epsilon() * norm)
        break;
    }
    bool singular = false;
    for (int i = used - 1; i >= 0; --i) {
      Real r = g[i];
      for (int j = i + 1; j < used; ++j)
        r -= H(i, j) * y[j];
      if (std::abs(H(i, i)) <= std::numeric_limits<Real>::min()) {
        singular = true;
        break;
      }
      y[i] = r / H(i, i);
    }
    if (singular)
      break;
    for (int j = 0; j < used; ++j)
      MultiFab::Saxpy(x, y[j], *basis[j], 0, 0, nComp, 0);
    result.residual = true_residual();
    if (!std::isfinite(result.residual))
      break;
  }
  result.converged =
      std::isfinite(result.residual) && result.residual <= target;
  return result;
}
