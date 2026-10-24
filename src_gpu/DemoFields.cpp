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
      eSolver(g.nodeBa, g.dm, 3, c.gmresRestart, g.geom),
      divSolver(g.ba, g.dm, 1, c.gmresRestart, g.geom) {
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
  jHat.SumBoundary(period);
  jHat.OverrideSync(period);
  const auto dxi = grid.geom.InvCellSizeArray();
  jHat.mult(dxi[0] * dxi[1] * dxi[2], 0, 3, 0);
  jHat.FillBoundary(period);
  nodeMM.SumBoundary(period);
  nodeMM.OverrideSync(period);
  nodeMM.FillBoundary(period);
}
void DemoFields::curl_node_to_center(MultiFab& input, MultiFab& output) {
  input.OverrideSync(grid.geom.periodicity());
  input.FillBoundary(grid.geom.periodicity());
  const auto invDx = grid.geom.InvCellSizeArray();
  for (MFIter mfi(output); mfi.isValid(); ++mfi) {
    const auto a = input.const_array(mfi);
    const auto b = output.array(mfi);
    ParallelFor(mfi.validbox(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) noexcept {
                  b(i, j, k, c) = curl(a, i, j, k, c, invDx, true);
                });
  }
  output.FillBoundary(grid.geom.periodicity());
}
void DemoFields::curl_center_to_node(MultiFab& input, MultiFab& output) {
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
  output.FillBoundary(grid.geom.periodicity());
}
void DemoFields::divergence(MultiFab& input, MultiFab& output) {
  input.OverrideSync(grid.geom.periodicity());
  input.FillBoundary(grid.geom.periodicity());
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
  curl_node_to_center(input, curlCenter);
  curl_center_to_node(curlCenter, curlNode);
  divergence(input, divE);
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
  curl_node_to_center(nodeEth, curlCenter);
  MultiFab::Saxpy(centerB, -config.dt, curlCenter, 0, 0, 3, 0);
  refresh_node_b();
}
void DemoFields::apply_div_e_operator(MultiFab& input, MultiFab& output) {
  input.FillBoundary(grid.geom.periodicity());
  for (MFIter mfi(output); mfi.isValid(); ++mfi) {
    auto x = input.const_array(mfi);
    auto m = centerMM.const_array(mfi);
    auto y = output.array(mfi);
    ParallelFor(
        mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          Real sum = 0;
          for (int di = -1; di <= 1; ++di)
            for (int dj = -1; dj <= 1; ++dj)
              for (int dk = -1; dk <= 1; ++dk)
                sum += m(i, j, k, (di + 1) * 9 + (dj + 1) * 3 + dk + 1) *
                       x(i + di, j + dj, k + dk);
          y(i, j, k) = fourPi * fourPi * sum;
        });
  }
}
bool DemoFields::solve_div_e() {
  divergence(nodeE, divE);
  // Production rhoTheta is fixed at .51 (independent of field theta).
  constexpr Real rhoTheta = .51;
  MultiFab::LinComb(scalarRhs, 1 / rhoTheta, divE, 0,
                    -fourPi * (1 - rhoTheta) / rhoTheta, rhoOld, 0, 0, 1, 0);
  MultiFab::Saxpy(scalarRhs, -fourPi, rho, 0, 0, 1, 0);
  const Real cells = static_cast<Real>(grid.geom.Domain().numPts());
  const Real mean = scalarRhs.sum(0) / cells;
  scalarRhs.plus(-mean, 0, 1, 0);
  phi.setVal(0);
  const auto result = divSolver.solve(
      phi, scalarRhs,
      [this](MultiFab& x, MultiFab& y) { apply_div_e_operator(x, y); },
      config.divETolerance, config.divEMaxIter);
  phi.plus(-phi.sum(0) / cells, 0, 1, 0);
  phi.FillBoundary(grid.geom.periodicity());
  lastDivResidual = result.residual;
  lastDivIterations = result.iterations;
  return result.converged;
}
Real DemoFields::field_energy() {
  const auto dx = grid.geom.CellSizeArray();
  const Real e = MultiFab::Dot(*owner, nodeE, 0, nodeE, 0, 3, 0);
  const Real b = MultiFab::Dot(centerB, 0, centerB, 0, 3, 0);
  return (e + b) * dx[0] * dx[1] * dx[2] / (2 * fourPi);
}
Real DemoFields::gauss_error() {
  divergence(nodeE, divE);
  MultiFab::Saxpy(divE, -fourPi * .51, rho, 0, 0, 1, 0);
  MultiFab::Saxpy(divE, -fourPi * .49, rhoOld, 0, 0, 1, 0);
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
  phi.setVal(0);
  centerMM.setVal(0);
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
