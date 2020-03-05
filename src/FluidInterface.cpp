#include <AMReX_MultiFabUtil.H>

#include "FluidInterface.h"
#include "Utility.h"
#include "GridUtility.h"

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

void FluidInterface::regrid(const amrex::BoxArray& centerBAIn,
                            const amrex::DistributionMapping& dmIn) {
  std::string nameFunc = "FluidInterface::regrid";
  Print() << nameFunc << " is runing..." << std::endl;

  if (centerBAIn == centerBA)
    return;

  centerBA = centerBAIn;
  nodeBA = convert(centerBA, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });
  dm = dmIn;

  const bool doCopy = true;
  distribute_FabArray(nodeFluid, nodeBA, dm, nVarCoupling, nGst, doCopy);
  distribute_FabArray(centerB, centerBA, dm, nDimMax, nGst, doCopy);
}

void FluidInterface::set_geom(const int nGstIn, const amrex::Geometry& geomIn) {
  nGst = nGstIn;

  geom = geomIn;

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
  int ifab = 0;
  for (MFIter mfi(nodeFluid); mfi.isValid(); ++mfi) {
    ifab++;
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

  Print() << "action = " << action << " nCount = " << nCount << std::endl;
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
  // All centerB, including all ghost cells are accurate.
  average_node_to_cellcenter(centerB, 0, nodeFluid, iBx, centerB.nComp(),
                             centerB.nGrow());

  // currentMF is just an alias of current components of nodeFluid.
  MultiFab currentMF(nodeFluid, make_alias, iJx, nDimMax);

  // The outmost layer of currentMF can not be calculated from centerB
  curl_center_to_node(centerB, currentMF, geom.InvCellSize());
  currentMF.mult(1.0 / (getNo2SiL() * fourPI * 1e-7), currentMF.nGrow());

  currentMF.FillBoundary(geom.periodicity());

  /*
  Q: The outmost layer of currentMF is not accurate. Why not use
  apply_float_boundary to fill in first-order estimation?

  A: If the whole domain is just ONE block, it will work. Otherwise, it will
  not. For example, For a 2D simulation domain of 6x3 with 2 blocks. In the
  x-direction, block-1 convers cell 0 (c+0) to cell 2 (c+2), and block-2 covers
  c+3 to c+5. The node (n+5, n-1) is the corner ghost node for the block-1, and
  it is the face ghost node for the block-2. On block-2, this node can be
  calculated from curl_center_to_node(centerB, currentMF....). However, block-1
  does not know how to calculate it, and FillBoundary will also NOT copy this
  node from block-2 to block-1, because this node is a boundary node and it is
  not covered by any physical node.

  So, we should keep in mind, the variables that are directly received from the
  MHD side, such as the magnetic fields, are accurate on all ghost nodes, but
  the current releated variables (current, plasma velocities, and electric
  field) are unknown at the outmost boundary node layer.
  */
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

void FluidInterface::set_plasma_charge_and_mass(amrex::Real qomEl) {

  if (!useElectronFluid) {
    MoMi_S[0] = QoQi_S[0] / qomEl;
  }

  SumMass = 0.0;
  for (int is = 0; is < nSIn; is++)
    SumMass += MoMi_S[is];

  // Fix the values of MoMi_S and QoQi_S;
  MoMi0_S = new double[nSIn];
  QoQi0_S = new double[nSIn];
  for (int i = 0; i < nSIn; i++) {
    MoMi0_S[i] = MoMi_S[i];
    QoQi0_S[i] = QoQi_S[i];
  }

  if (doSplitSpecies) {
    int idx;
    for (int i = 0; i < nS; i++) {
      idx = iSPic2Mhd_I[i];
      QoQi_S[i] = QoQi0_S[idx];
      MoMi_S[i] = MoMi0_S[idx];
    }
  }

  Print() << "=========== Plasma species ============" << std::endl;
  for (int is = 0; is < nS; is++) {
    Print() << "Q/Qi[" << is << "] = " << QoQi_S[is] << std::endl;

    Print() << "M/Mi[" << is << "] = " << MoMi_S[is] << std::endl;
  }
  if (!useMhdPe && !useElectronFluid)
    Print() << "Pe/Ptotal = " << PeRatio << std::endl;
  Print() << "===================================" << std::endl;
}

void FluidInterface::load_balance(const DistributionMapping& dmIn) {
  dm = dmIn;

  redistribute_FabArray(nodeFluid, dm); // false?
  redistribute_FabArray(centerB, dm);   // false?
}