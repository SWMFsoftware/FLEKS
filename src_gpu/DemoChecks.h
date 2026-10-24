#ifndef DEMO_CHECKS_H_
#define DEMO_CHECKS_H_

#include <cmath>
#include <string>

#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParReduce.H>

#include "DemoConfig.h"
#include "DemoGrid.h"

class DemoChecks {
public:
  // Read the parameters and create the grid and the MultiFabs
  DemoChecks(const DemoGrid& grid, const DemoConfig& config);

  // Run all tests. Return the number of failed checks.
  int run();

  void test_fill();
  void test_reductions();
  void test_fill_boundary();
  void test_arithmetic();
  void test_laplacian();
  void test_parallel_copy();
  void test_device_to_host();
  void time_kernels();

  void apply_laplacian();

  // Global max |mfIn - expected(i,j,k)| over the valid cells and nGhostIn
  // ghost cells
  template <typename F>
  amrex::Real max_error(const amrex::MultiFab& mfIn, int nGhostIn,
                        F const& expected) const {
    auto const& arrs = mfIn.const_arrays();
    amrex::Real err = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpMax>{}, amrex::TypeList<amrex::Real>{},
        mfIn, amrex::IntVect(nGhostIn),
        [=] AMREX_GPU_DEVICE(int iBox, int i, int j,
                             int k) -> amrex::GpuTuple<amrex::Real> {
          return { std::abs(arrs[iBox](i, j, k) - expected(i, j, k)) };
        });
    amrex::ParallelDescriptor::ReduceRealMax(err);
    return err;
  }

  // Periodic test function on the unit cube, evaluated at the cell center of
  // (i,j,k). It is also valid for ghost cells outside the domain.
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real exact_value(
      int i, int j, int k, amrex::GpuArray<amrex::Real, 3> const& dx) noexcept {
    const amrex::Real x = (i + 0.5) * dx[0];
    const amrex::Real y = (j + 0.5) * dx[1];
    const amrex::Real z = (k + 0.5) * dx[2];
    return 1 + std::sin(twoPi * x) * std::cos(twoPi * y) +
           0.5 * std::sin(twoPi * z);
  }

  // Exact 7-point discrete Laplacian of exact_value() for a uniform grid
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real exact_laplacian(
      int i, int j, int k, amrex::GpuArray<amrex::Real, 3> const& dx) noexcept {
    const amrex::Real x = (i + 0.5) * dx[0];
    const amrex::Real y = (j + 0.5) * dx[1];
    const amrex::Real z = (k + 0.5) * dx[2];
    const amrex::Real dx2 = dx[0] * dx[0];
    const amrex::Real factor = (2 * std::cos(twoPi * dx[0]) - 2) / dx2;
    return 2 * factor * std::sin(twoPi * x) * std::cos(twoPi * y) +
           0.5 * factor * std::sin(twoPi * z);
  }

private:
  void report(const std::string& name, amrex::Real err, amrex::Real tol);

  void print_device_info() const;

  static constexpr amrex::Real twoPi = 2 * 3.14159265358979323846;

  static constexpr int nGhost = 1;

  int nCell = 128;
  int maxGridSize = 32;
  int nIter = 100;

  int nFail = 0;

  // Periodic unit cube
  const amrex::Box& domain;
  const amrex::Geometry& geom;
  const amrex::BoxArray& ba;
  const amrex::DistributionMapping& dm;

  // Test function with ghost cells, and its Laplacian
  amrex::MultiFab mf;
  amrex::MultiFab lap;
};

#endif
