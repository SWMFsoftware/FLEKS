
#include <AMReX_BC_TYPES.H>
#include <AMReX_PhysBCFunct.H>
#include <Constants.h>

#include "BC.h"

using namespace amrex;

// Use FabArray.setDomainBoundary()!!!!!!!!!!!!!!!!!!!
void zero_boundary_cpu(amrex::Box const& bx,
                       amrex::Array4<amrex::Real> const& arr, const int iStart,
                       const int nComp, amrex::GeometryData const& geom,
                       const amrex::Real time, const amrex::BCRec* bcr,
                       const int bcomp, const int orig_comp) {
  //"g" means "global".
  const Box& gbx = geom.Domain();
  const auto glo = gbx.loVect(); // Do not include ghost cells.
  const auto ghi = gbx.hiVect();

  int igMin = glo[ix_], igMax = ghi[ix_];
  int jgMin = glo[iy_], jgMax = ghi[iy_];
  int kgMin = glo[iz_], kgMax = ghi[iz_];

  // Include ghost cells.
  const auto lo = bx.loVect();
  const auto hi = bx.hiVect();

  int iMin = lo[ix_], iMax = hi[ix_];
  int jMin = lo[iy_], jMax = hi[iy_];
  int kMin = lo[iz_], kMax = hi[iz_];

  const Real val = 0;

  // x left
  if (bcr->lo(ix_) == BCType::ext_dir) {
    for (int iVar = iStart; iVar < nComp; iVar++)
      for (int k = kMin; k <= kMax; k++)
        for (int j = jMin; j <= jMax; j++)
          for (int i = iMin; i <= igMin - 1 + nVirGst; i++) {
            arr(i, j, k, iVar) = val;
          }
  }

  // x right
  if (bcr->hi(ix_) == BCType::ext_dir) {
    for (int iVar = iStart; iVar < nComp; iVar++)
      for (int k = kMin; k <= kMax; k++)
        for (int j = jMin; j <= jMax; j++)
          for (int i = igMax + 1 - nVirGst; i <= iMax; i++) {
            arr(i, j, k, iVar) = val;
          }
  }

  // y left
  if (bcr->lo(iy_) == BCType::ext_dir) {
    for (int iVar = iStart; iVar < nComp; iVar++)
      for (int k = kMin; k <= kMax; k++)
        for (int j = jMin; j <= jgMin - 1 + nVirGst; j++)
          for (int i = iMin; i <= iMax; i++) {
            arr(i, j, k, iVar) = val;
          }
  }

  // y right
  if (bcr->lo(iy_) == BCType::ext_dir) {
    for (int iVar = iStart; iVar < nComp; iVar++)
      for (int k = kMin; k <= kMax; k++)
        for (int j = jgMax + 1 - nVirGst; j <= jMax; j++)
          for (int i = iMin; i <= iMax; i++) {
            arr(i, j, k, iVar) = val;
          }
  }

  // z left
  if (bcr->lo(iz_) == BCType::ext_dir) {
    for (int iVar = iStart; iVar < nComp; iVar++)
      for (int k = kMin; k <= kgMin - 1 + nVirGst; k++)
        for (int j = jMin; j <= jMax; j++)
          for (int i = iMin; i <= iMax; i++) {
            arr(i, j, k, iVar) = val;
          }
  }

  // z right
  if (bcr->hi(iz_) == BCType::ext_dir) {
    for (int iVar = iStart; iVar < nComp; iVar++)
      for (int k = kgMax + 1 - nVirGst; k <= kMax; k++)
        for (int j = jMin; j <= jMax; j++)
          for (int i = iMin; i <= iMax; i++) {
            arr(i, j, k, iVar) = val;
          }
  }
}

void apply_zero_boundary(MultiFab& mf, const Geometry& geom) {

  if (geom.isAllPeriodic())
    return;
  if (mf.nGrow() == 0)
    return;

  Vector<BCRec> bc(mf.nComp());
  for (int n = 0; n < mf.nComp(); ++n) {
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      bc[n].setLo(idim, BCType::ext_dir);
      bc[n].setHi(idim, BCType::ext_dir);
    }
  }

  CpuBndryFuncFab bc_func(zero_boundary_cpu);

  PhysBCFunct<CpuBndryFuncFab> physbcf(geom, bc, bc_func);
  physbcf.FillBoundary(mf, 0, mf.nComp(), 0.0, 0);
}