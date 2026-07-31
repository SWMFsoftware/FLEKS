#include "LightWaveIC.h"

#include <AMReX.H>
#include <cmath>

using namespace amrex;

// Replicates the kernel Pic::fill_lightwaves(48.0) exactly: an oblique
// (kx:ky = 0.6:0.8) transverse plane wave filling both E and B on the node
// grids and B on the cell-centred grid, wavelength 48.0, time 0.0.
void LightWaveIC::set_fields(PicICFields& f) const {
  const Real wavelength = 48.0;
  const int EorB = -1;  // -1: both E and B
  const Real time = 0.0;
  const int lev = -1;   // -1: all levels

  const int nLev = f.n_lev();

  for (int iLev = 0; iLev < nLev; iLev++) {
    auto& nodeE = f.node_E(iLev);
    auto& nodeB = f.node_B(iLev);
    nodeE.setVal(0.0);
    nodeB.setVal(0.0);
    f.center_B(iLev).setVal(0.0);

    if (lev != -1 && iLev != lev) continue;

    for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
      FArrayBox& fab = nodeE[mfi];
      FArrayBox& fab2 = nodeB[mfi];

      const Box& box = mfi.fabbox();
      const Array4<Real>& arrE = fab.array();
      const Array4<Real>& arrB = fab2.array();
      const auto& prob_lo = f.geom(iLev).ProbLo();
      const auto& dx = f.geom(iLev).CellSize();

      ParallelFor(box, [&](int i, int j, int k) {
        IntVect ijk = {AMREX_D_DECL(i, j, k)};
        const Real phase =
            (2.0 * (dPI) *
             ((prob_lo[0] + dx[0] * i) * 0.6 +
              (prob_lo[1] + dx[1] * j) * 0.8 - time)) /
            wavelength;
        if (EorB == -1 || EorB == 0) {
          arrE(ijk, ix_) = -0.8 * std::sin(phase);
          arrE(ijk, iy_) = 0.6 * std::sin(phase);
          arrE(ijk, iz_) = -std::cos(phase);
        }
        if (EorB == -1 || EorB == 1) {
          arrB(ijk, ix_) = -0.8 * std::cos(phase);
          arrB(ijk, iy_) = 0.6 * std::cos(phase);
          arrB(ijk, iz_) = std::sin(phase);
        }
      });
    }
    nodeE.FillBoundary(f.geom(iLev).periodicity());
    nodeB.FillBoundary(f.geom(iLev).periodicity());
  }

  for (int iLev = 0; iLev < nLev; iLev++) {
    auto& centerB = f.center_B(iLev);
    for (MFIter mfi(centerB); mfi.isValid(); ++mfi) {
      FArrayBox& fab = centerB[mfi];
      const Box& box = mfi.fabbox();
      const Array4<Real>& arrcB = fab.array();
      const auto& prob_lo = f.geom(iLev).ProbLo();
      const auto& dx = f.geom(iLev).CellSize();

      ParallelFor(box, [&](int i, int j, int k) {
        IntVect ijk = {AMREX_D_DECL(i, j, k)};
        const Real phase =
            (2.0 * (dPI) *
             (((prob_lo[0] + dx[0] * (i + 0.5)) * 0.6 +
               (prob_lo[1] + dx[1] * (j + 0.5)) * 0.8 - time))) /
            wavelength;
        if (EorB == -1 || EorB == 1) {
          arrcB(ijk, ix_) = -0.8 * std::cos(phase);
          arrcB(ijk, iy_) = 0.6 * std::cos(phase);
          arrcB(ijk, iz_) = std::sin(phase);
        }
      });
    }
    centerB.FillBoundary(f.geom(iLev).periodicity());
  }
}
