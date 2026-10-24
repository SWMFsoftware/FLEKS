// Exact production FLEKS field-operator bodies run on CPU pinned memory,
// compared with the demo kernels on the configured AMReX backend. The bodies
// are generated from src/GridUtility.cpp by validate_demo.py at build time.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_Loop.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "DemoFieldKernels.h"
#include "DemoFields.h"
#include "ProductionOperators.h"

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
reference_value(int i, int j, int k, int c, Real dx, bool nodal) noexcept {
  const Real shift = nodal ? 0 : .5;
  const Real x = (i + shift) * dx;
  const Real y = (j + shift) * dx;
  const Real z = (k + shift) * dx;
  return std::sin(2 * demo_field::pi * ((c + 1) * x + y)) +
         .3 * std::cos(2 * demo_field::pi * (y + (c + 1) * z)) +
         .2 * std::sin(2 * demo_field::pi * (z - (c + 1) * x));
}

int check_operator(int op, const BoxArray& cells, const DistributionMapping& dm,
                   const Geometry& geom, std::ofstream& productionFile,
                   std::ofstream& demoFile) {
  const bool toCenter = op % 2 == 0;
  const bool curl = op >= 4;
  const bool gradient = op < 2;
  const int nComp = gradient || curl ? 3 : 1;
  BoxArray nodes = convert(cells, IntVect(1));
  const BoxArray& inputBA = toCenter ? nodes : cells;
  const BoxArray& outputBA = toCenter ? cells : nodes;
  MFInfo hostInfo;
  hostInfo.SetArena(The_Pinned_Arena());
  MultiFab hostInput(inputBA, dm, 3, 2, hostInfo);
  MultiFab hostOutput(outputBA, dm, nComp, 0, hostInfo);
  MultiFab deviceInput(inputBA, dm, 3, 2);
  MultiFab deviceOutput(outputBA, dm, nComp, 0);
  const Real dx = geom.CellSize(0);
  const auto invDx = geom.InvCellSizeArray();
  for (MFIter mfi(hostInput); mfi.isValid(); ++mfi) {
    const auto host = hostInput.array(mfi);
    const auto device = deviceInput.array(mfi);
    const Box box = mfi.fabbox();
    LoopOnCpu(box, 3, [=](int i, int j, int k, int c) {
      host(i, j, k, c) = reference_value(i, j, k, c, dx, toCenter);
    });
    ParallelFor(box, 3, [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) {
      device(i, j, k, c) = reference_value(i, j, k, c, dx, toCenter);
    });
  }
  // These calls execute the current production numerical expressions on CPU.
  const Real* hInvDx = geom.InvCellSize();
  if (op == 0)
    production_reference::grad_node_to_center(hostInput, hostOutput, hInvDx);
  else if (op == 1)
    production_reference::grad_center_to_node(hostInput, hostOutput, hInvDx);
  else if (op == 2)
    production_reference::div_node_to_center(hostInput, hostOutput, hInvDx);
  else if (op == 3)
    production_reference::div_center_to_node(hostInput, hostOutput, hInvDx);
  else if (op == 4)
    production_reference::curl_node_to_center(hostInput, hostOutput, hInvDx);
  else
    production_reference::curl_center_to_node(hostInput, hostOutput, hInvDx);
  for (MFIter mfi(deviceOutput); mfi.isValid(); ++mfi) {
    const auto input = deviceInput.const_array(mfi);
    const auto output = deviceOutput.array(mfi);
    ParallelFor(mfi.validbox(), nComp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) {
                  if (curl) {
                    output(i, j, k, c) =
                        demo_field::curl(input, i, j, k, c, invDx, toCenter);
                  } else if (gradient) {
                    output(i, j, k, c) = demo_field::derivative(
                        input, i, j, k, 0, c, invDx, toCenter);
                  } else {
                    Real result = 0;
                    for (int d = 0; d < 3; ++d)
                      result += demo_field::derivative(input, i, j, k, d, d,
                                                       invDx, toCenter);
                    output(i, j, k, c) = result;
                  }
                });
  }
  Real error = 0, scale = 0;
  for (MFIter mfi(deviceOutput); mfi.isValid(); ++mfi) {
    const auto& device = deviceOutput[mfi];
    FArrayBox hostCopy(device.box(), nComp, The_Pinned_Arena());
    Gpu::dtoh_memcpy_async(hostCopy.dataPtr(), device.dataPtr(),
                           device.nBytes());
    Gpu::streamSynchronize();
    const auto actual = hostCopy.const_array();
    const auto expected = hostOutput.const_array(mfi);
    const int boxId = mfi.index();
    LoopOnCpu(mfi.validbox(), nComp, [&](int i, int j, int k, int c) {
      error =
          std::max(error, std::abs(actual(i, j, k, c) - expected(i, j, k, c)));
      scale = std::max(scale, std::abs(expected(i, j, k, c)));
      if (productionFile.is_open()) {
        productionFile << op << ',' << boxId << ',' << i << ',' << j << ',' << k
                       << ',' << c << ',' << expected(i, j, k, c) << '\n';
        demoFile << op << ',' << boxId << ',' << i << ',' << j << ',' << k
                 << ',' << c << ',' << actual(i, j, k, c) << '\n';
      }
    });
  }
  ParallelDescriptor::ReduceRealMax(error);
  ParallelDescriptor::ReduceRealMax(scale);
  const Real tolerance = 1e-12 * std::max(Real(1), scale);
  const char* names[] = { "grad-node-to-center", "grad-center-to-node",
                          "div-node-to-center",  "div-center-to-node",
                          "curl-node-to-center", "curl-center-to-node" };
  Print() << (error <= tolerance ? "PASS " : "FAIL ") << "production-CPU/"
          << names[op] << " error=" << error << " tolerance=" << tolerance
          << '\n';
  return error > tolerance;
}

// Compare the actual demo matrix-vector operator with production's composed
// I - delt^2 lap + delt^2 grad(blended div), with zero particle response.
int check_composed_operator(DemoFields& fields, const DemoGrid& grid,
                            const DemoConfig& config, int op,
                            std::ofstream& productionFile,
                            std::ofstream& demoFile) {
  MFInfo hostInfo;
  hostInfo.SetArena(The_Pinned_Arena());
  MultiFab input(grid.nodeBa, grid.dm, 3, 3, hostInfo);
  MultiFab lap(grid.nodeBa, grid.dm, 3, 0, hostInfo);
  MultiFab gradient(grid.nodeBa, grid.dm, 3, 0, hostInfo);
  MultiFab vectorScratch(grid.ba, grid.dm, 3, 1, hostInfo);
  MultiFab divergence(grid.ba, grid.dm, 1, 1, hostInfo);
  MultiFab extendedDiv(grid.ba, grid.dm, 1, 1, hostInfo);
  BoxArray grownCells = grid.ba;
  grownCells.grow(2);
  MultiFab averaged(grownCells, grid.dm, 3, 0, hostInfo);
  MultiFab actual(grid.nodeBa, grid.dm, 3, 0);
  const Real dx = grid.geom.CellSize(0);
  const Real* invDx = grid.geom.InvCellSize();
  for (MFIter mfi(input); mfi.isValid(); ++mfi) {
    const auto a = input.array(mfi);
    LoopOnCpu(mfi.fabbox(), 3, [=](int i, int j, int k, int c) {
      a(i, j, k, c) = reference_value(i, j, k, c, dx, true);
    });
    const auto device = fields.nodeE.array(mfi);
    ParallelFor(fields.nodeE[mfi].box(), 3,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int c) {
                  device(i, j, k, c) = reference_value(i, j, k, c, dx, true);
                });
  }
  for (int c = 0; c < 3; ++c) {
    MultiFab inputComponent(input, make_alias, c, 1);
    MultiFab outputComponent(lap, make_alias, c, 1);
    production_reference::grad_node_to_center(inputComponent, vectorScratch,
                                              invDx);
    production_reference::div_center_to_node(vectorScratch, outputComponent,
                                             invDx);
  }
  production_reference::div_node_to_center(input, divergence, invDx);
  // Grow the valid region of the CPU average so the production extended
  // divergence can read its full transverse stencil without communication.
  production_reference::average_node_to_center(input, averaged);
  production_reference::div_center_to_center(averaged, extendedDiv, invDx);
  const Real blend = config.coefDiff;
  for (MFIter mfi(divergence); mfi.isValid(); ++mfi) {
    auto standard = divergence.array(mfi);
    auto extended = extendedDiv.const_array(mfi);
    LoopOnCpu(mfi.fabbox(), [=](int i, int j, int k) {
      standard(i, j, k) =
          (1 - blend) * standard(i, j, k) + blend * extended(i, j, k);
    });
  }
  production_reference::grad_center_to_node(divergence, gradient, invDx);
  fields.apply_e_operator(fields.nodeE, actual);
  Real error = 0, scale = 0;
  const Real delt2 = config.theta * config.dt * config.theta * config.dt;
  for (MFIter mfi(actual); mfi.isValid(); ++mfi) {
    const auto& device = actual[mfi];
    FArrayBox copy(device.box(), 3, The_Pinned_Arena());
    Gpu::dtoh_memcpy_async(copy.dataPtr(), device.dataPtr(), device.nBytes());
    Gpu::streamSynchronize();
    const auto values = copy.const_array();
    const auto original = input.const_array(mfi);
    const auto laplacian = lap.const_array(mfi);
    const auto gradDiv = gradient.const_array(mfi);
    const int box = mfi.index();
    LoopOnCpu(mfi.validbox(), 3, [&](int i, int j, int k, int c) {
      const Real expected =
          original(i, j, k, c) +
          delt2 * (-laplacian(i, j, k, c) + gradDiv(i, j, k, c));
      error = std::max(error, std::abs(values(i, j, k, c) - expected));
      scale = std::max(scale, std::abs(expected));
      if (productionFile.is_open()) {
        productionFile << op << ',' << box << ',' << i << ',' << j << ',' << k
                       << ',' << c << ',' << expected << '\n';
        demoFile << op << ',' << box << ',' << i << ',' << j << ',' << k << ','
                 << c << ',' << values(i, j, k, c) << '\n';
      }
    });
  }
  ParallelDescriptor::ReduceRealMax(error);
  ParallelDescriptor::ReduceRealMax(scale);
  const Real tolerance = 1e-12 * std::max(Real(1), scale);
  Print() << (error <= tolerance ? "PASS " : "FAIL ")
          << "production-CPU/composed-E-operator coef_diff=" << config.coefDiff
          << " error=" << error << " tolerance=" << tolerance << '\n';
  return error > tolerance;
}

int main(int argc, char* argv[]) {
  amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, []() {
    ParmParse pp("amrex");
    if (!pp.contains("use_gpu_aware_mpi"))
      pp.add("use_gpu_aware_mpi", 0);
    if (!pp.contains("the_arena_init_size"))
      pp.add("the_arena_init_size", 64 * 1024 * 1024);
  });
  int failures = 0;
  {
    ParmParse pp;
    int nCell = 8, maxGridSize = 4, writeFiles = 1;
    std::string prefix = "gpu_reference";
    pp.query("n_cell", nCell);
    pp.query("max_grid_size", maxGridSize);
    pp.query("reference_prefix", prefix);
    pp.query("write_files", writeFiles);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        nCell >= 4 && maxGridSize > 0,
        "Reference mesh requires n_cell>=4 and max_grid_size>0");
    Box domain(IntVect(0), IntVect(nCell - 1));
    RealBox realBox({ 0., 0., 0. }, { 1., 1., 1. });
    Geometry geom(domain, realBox, CoordSys::cartesian,
                  Array<int, 3>{ 1, 1, 1 });
    BoxArray cells(domain);
    cells.maxSize(maxGridSize);
    DistributionMapping dm(cells);
    std::ofstream productionFile, demoFile;
    if (writeFiles) {
      const std::string rank =
          ".rank" + std::to_string(ParallelDescriptor::MyProc());
      productionFile.open(prefix + rank + ".production.csv");
      demoFile.open(prefix + rank + ".demo.csv");
      AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
          productionFile && demoFile,
          "Cannot write production reference arrays");
      productionFile << std::setprecision(17)
                     << "operator,box,i,j,k,component,value\n";
      demoFile << std::setprecision(17)
               << "operator,box,i,j,k,component,value\n";
    }
    for (int op = 0; op < 6; ++op)
      failures += check_operator(op, cells, dm, geom, productionFile, demoFile);
    DemoConfig config;
    config.nCell = nCell;
    config.maxGridSize = maxGridSize;
    DemoGrid grid(config);
    DemoFields fields(grid, config);
    config.coefDiff = 0;
    failures += check_composed_operator(fields, grid, config, 6, productionFile,
                                        demoFile);
    config.coefDiff = .1;
    failures += check_composed_operator(fields, grid, config, 7, productionFile,
                                        demoFile);
    Print() << (failures ? "FAILED PRODUCTION OPERATOR CHECKS\n"
                         : "ALL PRODUCTION OPERATOR CHECKS PASSED\n");
  }
  amrex::Finalize();
  return failures ? 1 : 0;
}
