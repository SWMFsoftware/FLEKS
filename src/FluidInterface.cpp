#include <AMReX_MultiFabUtil.H>

#include "FluidInterface.h"
#include "Utility.h"

using namespace amrex;

void FluidInterface::init() {
  set_myrank(ParallelDescriptor::MyProc());
  set_nProcs(ParallelDescriptor::NProcs());
  int iCycle = 0;
  setCycle(iCycle);
}

void FluidInterface::receive_info_from_gm(const int* const paramInt,
                                          const double* const gridDim,
                                          const double* const paramDouble,
                                          const std::string& paramString) {
  std::stringstream* ss = nullptr;
  ReadFromGMinit(paramInt, gridDim, paramDouble, ss);
  readParam = paramString;
}

void FluidInterface::make_grid(const amrex::DistributionMapping& dmIn,
                               const amrex::Geometry& geomIn,
                               const amrex::BoxArray& centerBAIn,
                               const amrex::BoxArray& nodeBAIn,
                               const int nGst) {
  dm = dmIn;
  geom = geomIn;
  centerBA = centerBAIn;
  nodeBA = nodeBAIn;

  nodeFluid.define(nodeBA, dm, nVarCoupling, nGst);
  nodeFluid.setVal(0);

  centerB.define(centerBA, dm, nDimMax, nGst);
  centerB.setVal(0);

  Array<int, 3> period;

  // As a interface between PC and MHD. It can not be periodic unless MHD is 2D.
  for (int i = 0; i < nDimMax; i++) {
    if (i < getnDim()) {
      period[i] = 0;
    } else {
      period[i] = 1;
    }
  }

  geom.setPeriodicity(period);
}

int FluidInterface::loop_through_node(std::string action, double* const pos_DI,
                                      const double* const data,
                                      const int* const index) {
  bool doCount = false;
  bool doGetLoc = false;
  bool doFill = false;

  if (action == "count") {
    doCount = true;
  } else if (action == "loc") {
    doGetLoc = true;
  } else if (action == "fill") {
    doFill = true;
  } else {
    amrex::Abort("Error: unknown action!\n");
  }

  const Real* dx = geom.CellSize();
  const auto plo = geom.ProbLo();
  const auto phi = geom.ProbHi();

  const double no2siL = getNo2SiL();

  int nIdxCount = 0;
  int nCount = 0;
  for (MFIter mfi(nodeFluid); mfi.isValid(); ++mfi) {
    const Box& box = mfi.fabbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const Array4<Real>& arr = nodeFluid[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {
          if (doCount) {
            nCount++;
          } else if (doGetLoc) {
            pos_DI[nCount++] = (i * dx[ix_] + plo[ix_]) * no2siL;
            pos_DI[nCount++] = (j * dx[iy_] + plo[iy_]) * no2siL;
            if (getnDim() > 2)
              pos_DI[nCount++] = (k * dx[iz_] + plo[iz_]) * no2siL;
          } else if (doFill) {
            for (int iVar = 0; iVar < nVarFluid; iVar++) {
              int idx;
              idx = iVar + nVarFluid * (index[nIdxCount] - 1);
              arr(i, j, k, iVar) = data[idx];
            }
            nIdxCount++;
          }

        } // for k
  }

  return nCount;
}

int FluidInterface::count_couple_node_number() {
  return loop_through_node("count");
}

void FluidInterface::get_couple_node_loc(double* const pos_DI) {
  int tmp = loop_through_node("loc", pos_DI);
}

void FluidInterface::set_couple_node_value(const double* const data,
                                           const int* const index) {
  int tmp = loop_through_node("fill", nullptr, data, index);

  calc_current();
  normalize_fluid_variables();
  convert_moment_to_velocity();

  MultiFab currentMF(nodeFluid, make_alias, iJx, nDimMax);
}

void FluidInterface::calc_current() {
  average_node_to_cellcenter(centerB, 0, nodeFluid, iBx, centerB.nComp(),
                             centerB.nGrow());

  // currentMF is just an alias of current components of nodeFluid.
  MultiFab currentMF(nodeFluid, make_alias, iJx, nDimMax);

  curl_center_to_node(centerB, currentMF, geom.InvCellSize());
  currentMF.mult(1.0 / (getNo2SiL() * fourPI * 1e-7), currentMF.nGrow());
  currentMF.FillBoundary();

  // The current in the ghost cells can not be calculated from the nodeB. So
  // fill in the ghost cell current with float boundary condition.
  apply_float_boundary(currentMF, geom, 0, currentMF.nComp());
}

void FluidInterface::normalize_fluid_variables() {
  for (int i = 0; i < nodeFluid.nComp(); ++i) {
    MultiFab tmpMF(nodeFluid, make_alias, i, 1);
    tmpMF.mult(Si2No_V[i], tmpMF.nGrow());
  }

  centerB.mult(Si2NoB, centerB.nGrow());
}

void FluidInterface::convert_moment_to_velocity() {

  for (MFIter mfi(nodeFluid); mfi.isValid(); ++mfi) {
    const Box& box = mfi.fabbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const Array4<Real>& arr = nodeFluid[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {

          if (useMultiSpecies) {
            double Rhot = 0;
            for (int iIon = 0; iIon < nIon; ++iIon) {
              // Rho = sum(Rhoi) + Rhoe;
              Rhot += arr(i, j, k, iRho_I[iIon]) *
                      (1 + MoMi0_S[0] / MoMi0_S[iIon + 1]);
            } // iIon

            arr(i, j, k, iUx_I[0]) /= Rhot;
            arr(i, j, k, iUy_I[0]) /= Rhot;
            arr(i, j, k, iUz_I[0]) /= Rhot;
          } else {
            for (int iFluid = 0; iFluid < nFluid; ++iFluid) {
              const double& rho = arr(i, j, k, iRho_I[iFluid]);
              arr(i, j, k, iUx_I[iFluid]) /= rho;
              arr(i, j, k, iUy_I[iFluid]) /= rho;
              arr(i, j, k, iUz_I[iFluid]) /= rho;
            } // iFluid
          }   // else
        }
  }
}
