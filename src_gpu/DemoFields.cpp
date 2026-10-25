#include <algorithm>
#include <cmath>
#include <limits>

#include <AMReX.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_ParReduce.H>
#include <AMReX_Print.H>

#include "DemoFields.h"
using namespace amrex;
using namespace demo_field;

DemoFields::DemoFields(const DemoGrid& g, const DemoConfig& c)
    : nodeE(g.nodeBa, g.dm, 3, 1),
      nodeEth(g.nodeBa, g.dm, 3, 1),
      nodeB(g.nodeBa, g.dm, 3, 1),
      centerB(g.ba, g.dm, 3, 1),
      jHat(g.nodeBa, g.dm, 3, 1),
      moments(g.nodeBa, g.dm, 10, 1),
      rho(g.ba, g.dm, 1, 1),
      phi(g.ba, g.dm, 1, 1),
      centerMM(g.ba, g.dm, 27, 1),
      rhoOld(g.ba, g.dm, 1, 1),
      nodeMM(g.nodeBa, g.dm, 1, 1),
      grid(g),
      config(c),
      rhs(g.nodeBa, g.dm, 3, 1),
      curlCenter(g.ba, g.dm, 3, 1),
      curlNode(g.nodeBa, g.dm, 3, 1),
      divE(g.ba, g.dm, 1, 1),
      scalarRhs(g.ba, g.dm, 1, 1),
      centerAverage(g.ba, g.dm, 3, 1),
      alternativeDiv(g.ba, g.dm, 1, 1),
      jacobiInverse(g.ba, g.dm, 1, 0),
      scalarP(g.ba, g.dm, 1, 1),
      scalarQ(g.ba, g.dm, 1, 0),
      scalarR(g.ba, g.dm, 1, 0),
      scalarZ(g.ba, g.dm, 1, 0),
      scalarU(g.ba, g.dm, 1, 1),
      eSolver(g.nodeBa, g.dm, 3, c.gmresRestart, g.geom, c.batchedGmres),
      divSolver(g.ba, g.dm, 1, c.gmresRestart, g.geom, c.batchedGmres) {
  owner = nodeE.OwnerMask(grid.geom.periodicity());
  nodeE.setVal(0);
  nodeEth.setVal(0);
  nodeB.setVal(0);
  centerB.setVal(0);
  centerB.setVal(c.b0, 2, 1, 1);
  nodeB.setVal(c.b0, 2, 1, 1);
  moments.setVal(0);
  rho.setVal(0);
  rhoOld.setVal(0);
  phi.setVal(0);
  centerMM.setVal(0);
  clear_response();
}
void DemoFields::clear_response() {
  jHat.setVal(0);
  nodeMM.setVal(DemoRealMM(0));
}
void DemoFields::finish_response() {
  const auto period = grid.geom.periodicity();
  {
    BL_PROFILE("response.jHat.SumBoundary");
    jHat.SumBoundary(period);
  }
  {
    BL_PROFILE("response.jHat.OverrideSync");
    jHat.OverrideSync(period);
  }
  const auto dxi = grid.geom.InvCellSizeArray();
  jHat.mult(dxi[0] * dxi[1] * dxi[2], 0, 3, 0);
  {
    BL_PROFILE("response.jHat.FillBoundary");
    jHat.FillBoundary(period);
  }
  {
    BL_PROFILE("response.nodeMM.SumBoundary");
    nodeMM.SumBoundary(period);
  }
  {
    BL_PROFILE("response.nodeMM.OverrideSync");
    nodeMM.OverrideSync(period);
  }
  {
    BL_PROFILE("response.nodeMM.FillBoundary");
    nodeMM.FillBoundary(period);
  }
}
void DemoFields::curl_node_to_center(MultiFab& input, MultiFab& output,
                                     bool inputPrepared, bool fillOutput) {
  if (!inputPrepared) {
    input.OverrideSync(grid.geom.periodicity());
    input.FillBoundary(grid.geom.periodicity());
  }
  const auto invDx = grid.geom.InvCellSizeArray();
  for (MFIter mfi(output); mfi.isValid(); ++mfi) {
    const auto a = input.const_array(mfi);
    const auto b = output.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  b(i, j, k, c) = curl(a, i, j, k, c, invDx, true);
                });
  }
  if (fillOutput)
    output.FillBoundary(grid.geom.periodicity());
}
void DemoFields::curl_center_to_node(MultiFab& input, MultiFab& output,
                                     bool inputPrepared, bool fillOutput) {
  if (!inputPrepared)
    input.FillBoundary(grid.geom.periodicity());
  const auto invDx = grid.geom.InvCellSizeArray();
  for (MFIter mfi(output); mfi.isValid(); ++mfi) {
    const auto a = input.const_array(mfi);
    const auto b = output.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  b(i, j, k, c) = curl(a, i, j, k, c, invDx, false);
                });
  }
  output.OverrideSync(grid.geom.periodicity());
  if (fillOutput)
    output.FillBoundary(grid.geom.periodicity());
}
void DemoFields::divergence(MultiFab& input, MultiFab& output,
                            bool inputPrepared, bool fillOutput) {
  if (!inputPrepared) {
    input.OverrideSync(grid.geom.periodicity());
    input.FillBoundary(grid.geom.periodicity());
  }
  const auto invDx = grid.geom.InvCellSizeArray();
  for (MFIter mfi(output); mfi.isValid(); ++mfi) {
    const auto a = input.const_array(mfi);
    const auto b = output.array(mfi);
    ParallelFor(mfi.validbox(),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  b(i, j, k) = derivative(a, i, j, k, 0, 0, invDx, true) +
                               derivative(a, i, j, k, 1, 1, invDx, true) +
                               derivative(a, i, j, k, 2, 2, invDx, true);
                });
  }
  if (fillOutput)
    output.FillBoundary(grid.geom.periodicity());
}
void DemoFields::refresh_node_b() {
  centerB.FillBoundary(grid.geom.periodicity());
  for (MFIter mfi(nodeB); mfi.isValid(); ++mfi) {
    auto a = centerB.const_array(mfi);
    auto b = nodeB.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  Real sum = 0;
                  for (int p = 0; p < 2; ++p)
                    for (int q = 0; q < 2; ++q)
                      for (int r = 0; r < 2; ++r)
                        sum += a(i - p, j - q, k - r, c);
                  b(i, j, k, c) = .125 * sum;
                });
  }
  nodeB.OverrideSync(grid.geom.periodicity());
  nodeB.FillBoundary(grid.geom.periodicity());
}
void DemoFields::apply_e_operator(MultiFab& input, MultiFab& output) {
  // Production I - delt^2 lap + delt^2 grad(div) + 4pi delt M.
  // curlcurl reproduces the coefDiff=0 part exactly.
  // curl_node_to_center synchronizes input and fills its ghosts; the M*x
  // stencil below reads those ghosts, so keep this call first.
  curl_node_to_center(input, curlCenter);
  curl_center_to_node(curlCenter, curlNode, true, false);
  divergence(input, divE, true, false);
  for (MFIter mfi(centerAverage); mfi.isValid(); ++mfi) {
    auto x = input.const_array(mfi);
    auto avg = centerAverage.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  Real sum = 0;
                  for (int p = 0; p < 2; ++p)
                    for (int q = 0; q < 2; ++q)
                      for (int r = 0; r < 2; ++r)
                        sum += x(i + p, j + q, k + r, c);
                  avg(i, j, k, c) = .125 * sum;
                });
  }
  centerAverage.FillBoundary(grid.geom.periodicity());
  const auto invDx = grid.geom.InvCellSizeArray();
  for (MFIter mfi(alternativeDiv); mfi.isValid(); ++mfi) {
    auto a = centerAverage.const_array(mfi);
    auto standard = divE.const_array(mfi);
    auto alt = alternativeDiv.array(mfi);
    ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j,
                                                     int k) noexcept {
      alt(i, j, k) = center_divergence(a, i, j, k, invDx) - standard(i, j, k);
    });
  }
  alternativeDiv.FillBoundary(grid.geom.periodicity());
  const Real delt = config.theta * config.dt, delt2 = delt * delt,
             coef = fourPi * delt, diff = config.coefDiff;
  for (MFIter mfi(output); mfi.isValid(); ++mfi) {
    auto x = input.const_array(mfi);
    auto cc = curlNode.const_array(mfi);
    auto m = nodeMM.const_array(mfi);
    auto y = output.array(mfi);
    auto alt = alternativeDiv.const_array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  Real response = 0;
                  for (int dk = -1; dk <= 1; ++dk)
                    for (int dj = -1; dj <= 1; ++dj)
                      for (int di = -1; di <= 1; ++di) {
                        const int index = (dk + 1) * 81 + (dj + 1) * 27 +
                                          (di + 1) * 9 + c * 3;
                        for (int e = 0; e < 3; ++e)
                          response += m(i, j, k)[index + e] *
                                      x(i + di, j + dj, k + dk, e);
                      }
                  y(i, j, k, c) = x(i, j, k, c) +
                                  delt2 * (cc(i, j, k, c) +
                                           diff * derivative(alt, i, j, k, 0, c,
                                                             invDx, false)) +
                                  coef * response;
                });
  }
  output.OverrideSync(grid.geom.periodicity());
}
bool DemoFields::solve_e() {
  curl_center_to_node(centerB, curlNode);
  MultiFab::Copy(rhs, nodeE, 0, 0, 3, 0);
  MultiFab::Saxpy(rhs, config.theta * config.dt, curlNode, 0, 0, 3, 0);
  MultiFab::Saxpy(rhs, -fourPi * config.theta * config.dt, jHat, 0, 0, 3, 0);
  MultiFab::Copy(nodeEth, nodeE, 0, 0, 3, 0);
  const auto result = eSolver.solve(
      nodeEth, rhs,
      [this](MultiFab& x, MultiFab& y) { apply_e_operator(x, y); },
      config.solverTolerance, config.gmresMaxIter);
  lastResidual = result.residual;
  lastIterations = result.iterations;
  if (!result.converged)
    return false;
  MultiFab::LinComb(nodeE, -(1 - config.theta) / config.theta, nodeE, 0,
                    1 / config.theta, nodeEth, 0, 0, 3, 0);
  nodeE.OverrideSync(grid.geom.periodicity());
  nodeE.FillBoundary(grid.geom.periodicity());
  nodeEth.OverrideSync(grid.geom.periodicity());
  nodeEth.FillBoundary(grid.geom.periodicity());
  return true;
}
void DemoFields::advance_b() {
  curl_node_to_center(nodeEth, curlCenter, false, false);
  MultiFab::Saxpy(centerB, -config.dt, curlCenter, 0, 0, 3, 0);
  refresh_node_b();
}
void DemoFields::apply_div_e_operator(MultiFab& input, MultiFab& output) {
  // One launch over all boxes on the current stream, with a non-waiting
  // local halo copy: the solver consumes the output on the same stream, and
  // any later MFIter consumer synchronizes the stream before it starts.
  Gpu::NoSyncRegion noSync;
  input.FillBoundary(grid.geom.periodicity());
  const auto x = input.const_arrays();
  const auto m = centerMM.const_arrays();
  const auto y = output.arrays();
  ParallelFor(output, IntVect(0),
              [=] AMREX_GPU_DEVICE(int box, int i, int j, int k) noexcept {
                Real sum = 0;
                for (int di = -1; di <= 1; ++di)
                  for (int dj = -1; dj <= 1; ++dj)
                    for (int dk = -1; dk <= 1; ++dk)
                      sum += m[box](i, j, k,
                                    (di + 1) * 9 + (dj + 1) * 3 + dk + 1) *
                             x[box](i + di, j + dj, k + dk);
                y[box](i, j, k) = fourPi * fourPi * sum;
              });
}
bool DemoFields::solve_div_e() {
  divergence(nodeE, divE);
  // One kernel forms the right-hand side, one reduction returns its mean,
  // and the mean subtraction runs without a host wait.
  const Real cells = static_cast<Real>(grid.geom.Domain().numPts());
  {
    Gpu::NoSyncRegion noSync;
    const Real a = 1 / rhoTheta, b = -fourPi * (1 - rhoTheta) / rhoTheta,
               c = -fourPi;
    const auto out = scalarRhs.arrays();
    const auto div = divE.const_arrays();
    const auto old = rhoOld.const_arrays();
    const auto now = rho.const_arrays();
    ParallelFor(scalarRhs, IntVect(0),
                [=] AMREX_GPU_DEVICE(int box, int i, int j, int k) noexcept {
                  out[box](i, j, k) = a * div[box](i, j, k) +
                                      b * old[box](i, j, k) +
                                      c * now[box](i, j, k);
                });
  }
  const Real mean = scalarRhs.sum(0) / cells;
  {
    Gpu::NoSyncRegion noSync;
    const auto out = scalarRhs.arrays();
    const auto potential = phi.arrays();
    ParallelFor(scalarRhs, IntVect(0),
                [=] AMREX_GPU_DEVICE(int box, int i, int j, int k) noexcept {
                  out[box](i, j, k) -= mean;
                  potential[box](i, j, k) = 0;
                });
  }
  DemoSolveResult result;
  if (config.divEMethod == "gmres") {
    result = divSolver.solve(
        phi, scalarRhs,
        [this](MultiFab& x, MultiFab& y) { apply_div_e_operator(x, y); },
        config.divETolerance, config.divEMaxIter, true);
  } else if (config.divEMethod == "jacobi_gmres") {
    // Right preconditioning: solve (A D^-1) u = b, then phi = D^-1 u. The
    // residual of the preconditioned system equals that of A phi = b.
    update_div_e_jacobi();
    {
      Gpu::NoSyncRegion noSync;
      scalarU.setVal(0);
    }
    result = divSolver.solve(
        scalarU, scalarRhs,
        [this](MultiFab& x, MultiFab& y) {
          scale_by_jacobi(x, scalarP);
          apply_div_e_operator(scalarP, y);
        },
        config.divETolerance, config.divEMaxIter, true);
    scale_by_jacobi(scalarU, phi);
  } else {
    update_div_e_jacobi();
    result = solve_div_e_pcg();
  }
  phi.plus(-phi.sum(0) / cells, 0, 1, 0);
  phi.FillBoundary(grid.geom.periodicity());
  lastDivResidual = result.residual;
  lastDivIterations = result.iterations;
  return result.converged;
}
void DemoFields::update_div_e_jacobi() {
  // Inverse of the operator diagonal 4 pi^2 M_center. Cells without
  // particles have a zero row; they keep unit scaling.
  Gpu::NoSyncRegion noSync;
  const auto inverse = jacobiInverse.arrays();
  const auto m = centerMM.const_arrays();
  ParallelFor(jacobiInverse, IntVect(0),
              [=] AMREX_GPU_DEVICE(int box, int i, int j, int k) noexcept {
                const Real diagonal = fourPi * fourPi * m[box](i, j, k, 13);
                inverse[box](i, j, k) = diagonal > 0 ? 1 / diagonal : 1;
              });
}
void DemoFields::scale_by_jacobi(const MultiFab& source, MultiFab& destination) {
  Gpu::NoSyncRegion noSync;
  const auto inverse = jacobiInverse.const_arrays();
  const auto src = source.const_arrays();
  const auto dst = destination.arrays();
  ParallelFor(destination, IntVect(0),
              [=] AMREX_GPU_DEVICE(int box, int i, int j, int k) noexcept {
                dst[box](i, j, k) = inverse[box](i, j, k) * src[box](i, j, k);
              });
}
// Jacobi-preconditioned conjugate gradients for the symmetric positive
// semidefinite correction operator (null space: constants, removed from the
// right-hand side). Two reductions per iteration; vector updates are fused
// and do not wait. phi must be zero on entry. The stopping test uses the
// recurrence residual; the reported residual is recomputed from A phi.
DemoSolveResult DemoFields::solve_div_e_pcg() {
  DemoSolveResult result;
  auto sum2 = [](const MultiFab& a, const MultiFab& b, const MultiFab& c,
                 const MultiFab& d) {
    const auto aa = a.const_arrays(), bb = b.const_arrays();
    const auto cc = c.const_arrays(), dd = d.const_arrays();
    auto value = ParReduce(
        TypeList<ReduceOpSum, ReduceOpSum>{}, TypeList<Real, Real>{}, a,
        IntVect(0),
        [=] AMREX_GPU_DEVICE(int box, int i, int j, int k)
            -> GpuTuple<Real, Real> {
          return { aa[box](i, j, k) * bb[box](i, j, k),
                   cc[box](i, j, k) * dd[box](i, j, k) };
        });
    Real out[2] = { amrex::get<0>(value), amrex::get<1>(value) };
    ParallelDescriptor::ReduceRealSum(out, 2);
    return std::make_pair(out[0], out[1]);
  };
  {
    Gpu::NoSyncRegion noSync;
    const auto b = scalarRhs.const_arrays();
    const auto inverse = jacobiInverse.const_arrays();
    const auto r = scalarR.arrays(), z = scalarZ.arrays(), p = scalarP.arrays();
    ParallelFor(scalarR, IntVect(0),
                [=] AMREX_GPU_DEVICE(int box, int i, int j, int k) noexcept {
                  r[box](i, j, k) = b[box](i, j, k);
                  z[box](i, j, k) = inverse[box](i, j, k) * b[box](i, j, k);
                  p[box](i, j, k) = z[box](i, j, k);
                });
  }
  auto [rz, rr] = sum2(scalarR, scalarZ, scalarR, scalarR);
  const Real target = config.divETolerance * std::max(std::sqrt(rr), 1e-14);
  result.initialResidual = std::sqrt(rr);
  Real recurrence = result.initialResidual;
  while (result.iterations < config.divEMaxIter && recurrence > target) {
    apply_div_e_operator(scalarP, scalarQ);
    const Real pq = sum2(scalarP, scalarQ, scalarP, scalarQ).first;
    if (!(pq > 0))
      break;
    const Real alpha = rz / pq;
    {
      Gpu::NoSyncRegion noSync;
      const auto inverse = jacobiInverse.const_arrays();
      const auto x = phi.arrays(), r = scalarR.arrays(), z = scalarZ.arrays();
      const auto p = scalarP.const_arrays(), q = scalarQ.const_arrays();
      ParallelFor(scalarR, IntVect(0),
                  [=] AMREX_GPU_DEVICE(int box, int i, int j, int k) noexcept {
                    x[box](i, j, k) += alpha * p[box](i, j, k);
                    r[box](i, j, k) -= alpha * q[box](i, j, k);
                    z[box](i, j, k) = inverse[box](i, j, k) * r[box](i, j, k);
                  });
    }
    const auto [rzNew, rrNew] = sum2(scalarR, scalarZ, scalarR, scalarR);
    recurrence = std::sqrt(rrNew);
    ++result.iterations;
    const Real beta = rzNew / rz;
    rz = rzNew;
    {
      Gpu::NoSyncRegion noSync;
      const auto p = scalarP.arrays();
      const auto z = scalarZ.const_arrays();
      ParallelFor(scalarP, IntVect(0),
                  [=] AMREX_GPU_DEVICE(int box, int i, int j, int k) noexcept {
                    p[box](i, j, k) = z[box](i, j, k) + beta * p[box](i, j, k);
                  });
    }
  }
  // True residual b - A phi, reported like the GMRES solves.
  apply_div_e_operator(phi, scalarQ);
  {
    Gpu::NoSyncRegion noSync;
    const auto b = scalarRhs.const_arrays();
    const auto q = scalarQ.arrays();
    ParallelFor(scalarQ, IntVect(0),
                [=] AMREX_GPU_DEVICE(int box, int i, int j, int k) noexcept {
                  q[box](i, j, k) = b[box](i, j, k) - q[box](i, j, k);
                });
  }
  result.residual = std::sqrt(sum2(scalarQ, scalarQ, scalarQ, scalarQ).first);
  result.converged =
      std::isfinite(result.residual) && result.residual <= target;
  return result;
}
Real DemoFields::field_energy() {
  const auto dx = grid.geom.CellSizeArray();
  const Real e = MultiFab::Dot(*owner, nodeE, 0, nodeE, 0, 3, 0);
  const Real b = MultiFab::Dot(centerB, 0, centerB, 0, 3, 0);
  return (e + b) * dx[0] * dx[1] * dx[2] / (2 * fourPi);
}
Real DemoFields::gauss_error() {
  divergence(nodeE, divE);
  MultiFab::Saxpy(divE, -fourPi * rhoTheta, rho, 0, 0, 1, 0);
  MultiFab::Saxpy(divE, -fourPi * (1 - rhoTheta), rhoOld, 0, 0, 1, 0);
  return divE.norm0(0);
}
Real DemoFields::div_b_error() {
  // Diagnose the evolved B directly: interpolating to nodes first can erase
  // a longitudinal checkerboard error before the divergence is measured.
  centerB.FillBoundary(grid.geom.periodicity());
  const auto invDx = grid.geom.InvCellSizeArray();
  for (MFIter mfi(curlNode); mfi.isValid(); ++mfi) {
    const auto input = centerB.const_array(mfi);
    const auto output = curlNode.array(mfi);
    ParallelFor(mfi.validbox(),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  Real value = 0;
                  for (int d = 0; d < 3; ++d)
                    value += derivative(input, i, j, k, d, d, invDx, false);
                  output(i, j, k, 0) = value;
                });
  }
  return curlNode.norm0(0);
}

int DemoFields::run_checks() {
  int failures = 0;
  auto check = [&](const char* name, Real error, Real tolerance) {
    bool pass = std::isfinite(error) && error <= tolerance;
    if (!pass)
      ++failures;
    Print() << "  " << (pass ? "PASS  " : "FAIL  ") << name
            << ": error = " << error << " (tol = " << tolerance << ")\n";
  };
  const Real tol =
      std::max(Real(5e-11), Real(500) * std::numeric_limits<Real>::epsilon());
  const int nc = config.nCell;
  const Real fac = 2 * nc * std::sin(pi / nc);
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    auto a = nodeE.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  int idx[3] = { i, j, k };
                  a(i, j, k, c) = std::sin(2 * pi * idx[(c + 1) % 3] / nc);
                });
  }
  curl_node_to_center(nodeE, curlCenter);
  for (MFIter mfi(curlCenter); mfi.isValid(); ++mfi) {
    auto a = curlCenter.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  int idx[3] = { i, j, k };
                  a(i, j, k, c) +=
                      fac * std::cos(2 * pi * (idx[(c + 2) % 3] + .5) / nc);
                });
  }
  check("FLEKS node-to-center curl", curlCenter.norm0(0, 3, IntVect(0)), tol);
  for (MFIter mfi(centerB); mfi.isValid(); ++mfi) {
    auto a = centerB.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  int idx[3] = { i, j, k };
                  a(i, j, k, c) =
                      std::sin(2 * pi * (idx[(c + 1) % 3] + .5) / nc);
                });
  }
  curl_center_to_node(centerB, curlNode);
  for (MFIter mfi(curlNode); mfi.isValid(); ++mfi) {
    auto a = curlNode.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  int idx[3] = { i, j, k };
                  a(i, j, k, c) +=
                      fac * std::cos(2 * pi * idx[(c + 2) % 3] / nc);
                });
  }
  check("FLEKS center-to-node curl", curlNode.norm0(0, 3, IntVect(0)), tol);
  refresh_node_b();
  for (MFIter mfi(nodeB); mfi.isValid(); ++mfi) {
    auto a = nodeB.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  int idx[3] = { i, j, k };
                  a(i, j, k, c) -= std::cos(pi / nc) *
                                   std::sin(2 * pi * idx[(c + 1) % 3] / nc);
                });
  }
  check("Cell-to-node magnetic interpolation", nodeB.norm0(0, 3, IntVect(0)),
        tol);
  const Real nodeCount = nc * static_cast<Real>(nc) * nc;
  rhs.setVal(1);
  check("Unique periodic nodal ownership",
        std::abs(eSolver.dot(rhs, rhs) - 3 * nodeCount), tol * nodeCount);

  // Each cell scatters 1/8 to eight corners. Every periodic node must receive
  // exactly one after local+MPI SumBoundary, for all 243 custom components.
  nodeMM.setVal(DemoRealMM(0));
  for (MFIter mfi(centerB); mfi.isValid(); ++mfi) {
    auto mm = nodeMM.array(mfi);
    ParallelFor(mfi.validbox(),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  for (int p = 0; p < 2; ++p)
                    for (int q = 0; q < 2; ++q)
                      for (int r = 0; r < 2; ++r)
                        for (int c = 0; c < 243; ++c)
                          Gpu::Atomic::Add(&mm(i + p, j + q, k + r)[c],
                                           Real(.125) * (c + 1));
                });
  }
  nodeMM.SumBoundary(grid.geom.periodicity());
  nodeMM.FillBoundary(grid.geom.periodicity());
  for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
    auto mm = nodeMM.const_array(mfi);
    auto err = rhs.array(mfi);
    ParallelFor(mfi.growntilebox(1), [=] AMREX_GPU_DEVICE(int i, int j,
                                                          int k) noexcept {
      Real error = 0;
      for (int c = 0; c < 243; ++c)
        error = amrex::max(error, std::abs(mm(i, j, k)[c] - (c + 1)));
      err(i, j, k) = error;
    });
  }
  check("Custom RealMM atomic scatter/SumBoundary/FillBoundary",
        rhs.norm0(0, 1), tol);

  // A nonzero diagonal response exercises the 243-value storage application.
  nodeMM.setVal(DemoRealMM(0));
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    auto mm = nodeMM.array(mfi);
    ParallelFor(mfi.validbox(),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  for (int c = 0; c < 3; ++c)
                    mm(i, j, k)[13 * 9 + 3 * c + c] = .2;
                });
  }
  nodeMM.FillBoundary(grid.geom.periodicity());
  apply_e_operator(nodeE, rhs);
  // For this transverse manufactured field, curlcurl is fac^2 E.
  MultiFab::Copy(nodeEth, nodeE, 0, 0, 3, 0);
  const Real scale = 1 + std::pow(config.theta * config.dt * fac, 2) +
                     fourPi * config.theta * config.dt * .2;
  nodeEth.mult(scale, 0, 3, 0);
  MultiFab::Subtract(nodeEth, rhs, 0, 0, 3, 0);
  check("Implicit E operator manufactured Fourier response",
        nodeEth.norm0(0, 3, IntVect(0)), tol);
  nodeEth.setVal(0);
  const auto solved = eSolver.solve(
      nodeEth, rhs,
      [this](MultiFab& x, MultiFab& y) { apply_e_operator(x, y); },
      config.solverTolerance, config.gmresMaxIter);
  MultiFab::Subtract(nodeEth, nodeE, 0, 0, 3, 0);
  check("Device GMRES manufactured solution", nodeEth.norm0(0, 3, IntVect(0)),
        std::max(tol, 10 * config.solverTolerance));
  check(
      "Device GMRES independently recomputed residual",
      solved.converged ? solved.residual : Real(1),
      std::max(tol, config.solverTolerance * std::sqrt(eSolver.dot(rhs, rhs))));
  // Mixed-dimensional plane wave catches the transverse smoothing in the
  // production coefDiff term; axis-aligned transverse waves cannot detect it.
  for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
    auto e = nodeE.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  e(i, j, k, c) = (c + 1) * std::sin(2 * pi * (i + j + k) / nc);
                });
  }
  apply_e_operator(nodeE, rhs);
  const Real cosine = std::cos(pi / nc), q = fac * cosine * cosine;
  const Real transverse = (1 + 2 * std::cos(2 * pi / nc)) / 3;
  const Real delta = config.theta * config.dt,
             mixture = cosine * cosine * transverse * transverse - 1;
  const Real diff = config.coefDiff;
  for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
    auto value = rhs.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  const Real a = c + 1;
                  const Real coefficient =
                      a +
                      delta * delta * q * q * (3 * a - 6 - diff * mixture * 6) +
                      fourPi * delta * .2 * a;
                  value(i, j, k, c) -=
                      coefficient * std::sin(2 * pi * (i + j + k) / nc);
                });
  }
  check("Full 3D coefDiff operator Fourier response",
        rhs.norm0(0, 3, IntVect(0)), tol);

  // Independently test the cell-centered solver, null constant mode, and
  // full27 scalar stencil indexing using a manufactured negative Laplacian.
  centerMM.setVal(0);
  const Real scalarScale = nc * static_cast<Real>(nc) / (fourPi * fourPi);
  for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
    auto m = centerMM.array(mfi);
    auto x = phi.array(mfi);
    ParallelFor(mfi.validbox(),
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  m(i, j, k, 13) = 6 * scalarScale;
                  m(i, j, k, 4) = -scalarScale;
                  m(i, j, k, 22) = -scalarScale;
                  m(i, j, k, 10) = -scalarScale;
                  m(i, j, k, 16) = -scalarScale;
                  m(i, j, k, 12) = -scalarScale;
                  m(i, j, k, 14) = -scalarScale;
                  x(i, j, k) = std::sin(2 * pi * (i + .5) / nc);
                });
  }
  apply_div_e_operator(phi, scalarRhs);
  MultiFab::Copy(divE, phi, 0, 0, 1, 0);
  phi.setVal(0);
  const auto scalarSolved = divSolver.solve(
      phi, scalarRhs,
      [this](MultiFab& x, MultiFab& y) { apply_div_e_operator(x, y); },
      config.solverTolerance, config.gmresMaxIter);
  MultiFab::Subtract(phi, divE, 0, 0, 1, 0);
  check("Cell-centered GMRES manufactured potential", phi.norm0(0),
        std::max(tol, 10 * config.solverTolerance));
  check("Cell-centered GMRES true residual",
        scalarSolved.converged ? scalarSolved.residual : Real(1),
        std::max(tol, config.solverTolerance *
                          std::sqrt(divSolver.dot(scalarRhs, scalarRhs))));
  // The optional preconditioned methods must reach the correction tolerance
  // on the same manufactured potential, with a recomputed true residual.
  const Real scalarRhsNorm = std::sqrt(divSolver.dot(scalarRhs, scalarRhs));
  update_div_e_jacobi();
  phi.setVal(0);
  const auto pcgSolved = solve_div_e_pcg();
  check("Jacobi PCG manufactured potential true relative residual",
        pcgSolved.converged ? pcgSolved.residual / scalarRhsNorm : Real(1),
        config.divETolerance);
  scalarU.setVal(0);
  const auto jacobiSolved = divSolver.solve(
      scalarU, scalarRhs,
      [this](MultiFab& x, MultiFab& y) {
        scale_by_jacobi(x, scalarP);
        apply_div_e_operator(scalarP, y);
      },
      config.divETolerance, config.divEMaxIter, true);
  check("Jacobi GMRES manufactured potential true relative residual",
        jacobiSolved.converged ? jacobiSolved.residual / scalarRhsNorm
                               : Real(1),
        config.divETolerance);
  phi.setVal(0);
  centerMM.setVal(0);
  // A varying, nonsymmetric block operator exercises multiple Arnoldi
  // columns and owner-masked batched projections, unlike a single Fourier
  // eigenmode. Run the optional solver's check even with the default MGS.
  {
    DemoLinearSolver batchedSolver(grid.nodeBa, grid.dm, 3, 10, grid.geom,
                                   true);
    MultiFab exact(grid.nodeBa, grid.dm, 3, 0);
    MultiFab answer(grid.nodeBa, grid.dm, 3, 1);
    MultiFab forcing(grid.nodeBa, grid.dm, 3, 0);
    for (MFIter mfi(exact); mfi.isValid(); ++mfi) {
      const auto value = exact.array(mfi);
      ParallelFor(mfi.validbox(), 3,
                  [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                    value(i, j, k, c) =
                        (c + 1) * std::cos(2 * pi * (i + j + k) / nc);
                  });
    }
    auto nonnormal = [&](MultiFab& input, MultiFab& output) {
      for (MFIter mfi(output); mfi.isValid(); ++mfi) {
        const auto x = input.const_array(mfi);
        const auto y = output.array(mfi);
        ParallelFor(
            mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              const Real d = 1 + .5 * std::sin(2 * pi * i / nc);
              y(i, j, k, 0) = d * x(i, j, k, 0) + .5 * x(i, j, k, 1);
              y(i, j, k, 1) = (d + .3) * x(i, j, k, 1) + .25 * x(i, j, k, 2);
              y(i, j, k, 2) = (d + .6) * x(i, j, k, 2);
            });
      }
      output.OverrideSync(grid.geom.periodicity());
    };
    nonnormal(exact, forcing);
    answer.setVal(0);
    const auto batchResult =
        batchedSolver.solve(answer, forcing, nonnormal, 1.e-10, 200);
    nonnormal(answer, exact);
    MultiFab::Subtract(exact, forcing, 0, 0, 3, 0);
    check("Batched CGS2 nonnormal true residual",
          batchResult.converged ? std::sqrt(batchedSolver.dot(exact, exact))
                                : Real(1),
          1.e-10 * std::sqrt(batchedSolver.dot(forcing, forcing)));
  }
  if (nc % 2 == 0) {
    centerB.setVal(0);
    for (MFIter mfi(centerB); mfi.isValid(); ++mfi) {
      const auto b = centerB.array(mfi);
      ParallelFor(mfi.validbox(),
                  [=] AMREX_GPU_DEVICE(int i, int j, int k)
                      noexcept { b(i, j, k, 0) = i % 2 == 0 ? 1 : -1; });
    }
    check("Evolved-B divergence detects checkerboard",
          std::abs(div_b_error() - 2 * nc), tol);
  }
  clear_response();
  nodeE.setVal(0);
  nodeEth.setVal(0);
  centerB.setVal(0);
  centerB.setVal(config.b0, 2, 1, 1);
  refresh_node_b();
  return failures;
}
