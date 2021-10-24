#include <AMReX_MultiFabUtil.H>

#include "FluidInterface.h"
#include "GridUtility.h"
#include "Utility.h"

using namespace amrex;
using namespace std;

void FluidInterface::init(int domainIDIn) {
  myrank = ParallelDescriptor::MyProc();
  domainID = domainIDIn;
  {
    std::stringstream ss;
    ss << "FLEKS" << domainID;
    domainName = ss.str();
    printPrefix = domainName + ": ";
  }
}

void FluidInterface::receive_info_from_gm(const int* const paramInt,
                                          const double* const gridDim,
                                          const double* const paramDouble,
                                          const std::string& paramString) {
  std::stringstream* ss = nullptr;
  read_from_GM(paramInt, gridDim, paramDouble, ss);
  readParam = paramString;
}

void FluidInterface::regrid(const amrex::BoxArray& centerBAIn,
                            const amrex::DistributionMapping& dmIn) {
  std::string nameFunc = "FluidInterface::regrid";

  // Why need 'isGridInitialized'? See the explaination in Domain::regrid().
  if (centerBAIn == centerBA && isGridInitialized) {
    // The interface grid does not change.
    return;
  }

  isGridEmpty = centerBAIn.empty();

  centerBA = centerBAIn;
  nodeBA = convert(centerBA, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });
  dm = dmIn;

  const bool doCopy = true;
  distribute_FabArray(nodeFluid, nodeBA, dm, nVarCoupling, nGst, doCopy);
  distribute_FabArray(centerB, centerBA, dm, nDimMax, nGst, doCopy);

  isGridInitialized = true;
}

void FluidInterface::set_geom(const int nGstIn, const amrex::Geometry& geomIn) {
  nGst = nGstIn;

  geom = geomIn;

  Array<int, 3> period;

  // As a interface between PC and MHD. It can not be periodic unless MHD is 2D.
  for (int i = 0; i < nDimMax; i++) {
    if (i < get_fluid_dimension()) {
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

  const double no2siL = get_No2SiL();

  int nIdxCount = 0;
  int nCount = 0;
  int ifab = 0;
  if (!nodeFluid.empty())
    for (MFIter mfi(nodeFluid); mfi.isValid(); ++mfi) {
      ifab++;
      // For each block, looping through all nodes, including ghost nodes.
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
              if (get_fluid_dimension() > 2)
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

  // Print() << "action = " << action << " nCount = " << nCount << std::endl;
  return nCount;
}

int FluidInterface::count_couple_node_number() {
  return loop_through_node("count");
}

void FluidInterface::get_couple_node_loc(double* const pos_DI) {
  loop_through_node("loc", pos_DI);
}

void FluidInterface::set_couple_node_value(const double* const data,
                                           const int* const index) {
  loop_through_node("fill", nullptr, data, index);

  calc_current();
  normalize_fluid_variables();
  convert_moment_to_velocity();

  MultiFab currentMF(nodeFluid, make_alias, iJx, nDimMax);
}

void FluidInterface::calc_current() {
  if (nodeFluid.empty())
    return;

  // All centerB, including all ghost cells are accurate.
  average_node_to_cellcenter(centerB, 0, nodeFluid, iBx, centerB.nComp(),
                             centerB.nGrow());

  // currentMF is just an alias of current components of nodeFluid.
  MultiFab currentMF(nodeFluid, make_alias, iJx, nDimMax);

  // The outmost layer of currentMF can not be calculated from centerB
  curl_center_to_node(centerB, currentMF, geom.InvCellSize());
  currentMF.mult(1.0 / (get_No2SiL() * fourPI * 1e-7), currentMF.nGrow());

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
                      (1 + MoMi_S[0] / MoMi_S[iIon + 1]);
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
  for (int is = 0; is < nS; is++)
    SumMass += MoMi_S[is];

  invSumMass = 1. / SumMass;
}

void FluidInterface::load_balance(const DistributionMapping& dmIn) {
  dm = dmIn;

  redistribute_FabArray(nodeFluid, dm); // false?
  redistribute_FabArray(centerB, dm);   // false?
}

//-------------------------------------------------------------------------

void FluidInterface::calc_normalized_units() {

  // normalization variables
  double RHOnorm, Bnorm, Jnorm, Pnorm;

  RHOnorm = Mnorm / (Lnorm * Lnorm * Lnorm);

  /* How to calculate Bnorm?
     1. Method 1
      In CGS unit, we have:
      1) [B] = [E]
      2) div(E) ~ rhoq, where rhoq is the charge density -> [E] = [rhoq]*[L]
      3) moment equation: d(rho*u)/dt ~ rhoq*u x B/c     ->
                         [RHO]*[U]/[T] = [rhoq]*[B]
        From the three equations above, we obtain: [B]=[E]=sqrt[RHO]*[U]

    2. Method 2
      In CGS, [P] = [RHO]*[U]*[U] = [B]*[B] -> [B]=sqrt[RHO]*[U]
  */
  Bnorm = sqrt(RHOnorm) * Unorm;

  // CGS Gauss's law: div(E) ~ rhoq -> [B] = [E] = [rhoq]*[L] = [Q]/[L]^2
  // -> [Q] = sqrt[RHO]*[U]*[L]^2 = sqrt[M*L]*[U].
  Qnorm = sqrt(Mnorm * Lnorm) * Unorm;
  Pnorm = RHOnorm * Unorm * Unorm;
  Jnorm = Qnorm * Unorm / (Lnorm * Lnorm * Lnorm);

  // SI -> CGS conversion
  Si2NoRho = 0.001; // [kg/m^3] -> [g/cm^3] and get numnber desity
  Si2NoV = 100.0;   // [m/s] -> [cm/s]
  Si2NoB = 1.0e4;   // [Tesla] -> [gauss]
  Si2NoP = 10.0;    // [Pa] -> [Ba]
  Si2NoL = 100;     // [m] ->[cm]

  /*
     1) SI: mu0*J_SI = curl(B_SI) with unit T/m
     2) CGS: 4pi/c_CGS*J_CGS = curl(B_CGS) with unit G/cm
     3) 1 T/m = 100 G/cm
     From 1) to 3) -> 100*mu0*J_SI = 4*pi/c_CGS*J_CGS ->
     J_CGS = 100*[ 4pi*10^(-7) ]/(4pi)*c_CGS = 10^(-5)*c_CGS,
     and c_CGS = Unorm.
  */
  Si2NoJ = 1.0e-5 * Unorm;

  /*
    1) E_SI = -U_SI x B_SI with unit m/s*T = 1e6 cm/s*G
    2) E_CGS*C_CGS = - U_CGS x B_CGS with unit cm/s*G
    -> E_CGS*C_CGS = 1e6*E_SI
    -> E_CGS = 1e6*E_SI/C_CGS = 1e6*E_SI/Unorm.
   */
  Si2NoE = 1e6 / Unorm;

  // Normalization: CGS -> non dimensional cgs
  Si2NoRho /= RHOnorm;
  Si2NoV /= Unorm;
  Si2NoB /= Bnorm;
  Si2NoP /= Pnorm;
  Si2NoJ /= Jnorm;
  Si2NoL /= Lnorm;
  Si2NoE /= Bnorm;

  Si2NoM = Si2NoRho * Si2NoV;

  No2SiV = 1. / Si2NoV;
  No2SiL = 1. / Si2NoL;

  Si2No_V[iBx] = Si2NoB;
  Si2No_V[iBy] = Si2NoB;
  Si2No_V[iBz] = Si2NoB;

  Si2No_V[iJx] = Si2NoJ;
  Si2No_V[iJy] = Si2NoJ;
  Si2No_V[iJz] = Si2NoJ;

  if (useElectronFluid) {
    Si2No_V[iEx] = Si2NoE;
    Si2No_V[iEy] = Si2NoE;
    Si2No_V[iEz] = Si2NoE;
  }

  if (useMhdPe)
    Si2No_V[iPe] = Si2NoP;
  if (useMultiSpecies)
    Si2No_V[iRhoTotal] = Si2NoRho;

  int iMax;
  iMax = nFluid;
  if (useMultiSpecies)
    iMax = nIon;

  for (int i = 0; i < iMax; ++i) {
    Si2No_V[iRho_I[i]] = Si2NoRho;
    Si2No_V[iRhoUx_I[i]] = Si2NoM;
    Si2No_V[iRhoUy_I[i]] = Si2NoM;
    Si2No_V[iRhoUz_I[i]] = Si2NoM;
    Si2No_V[iP_I[i]] = Si2NoP;
    if (useAnisoP)
      Si2No_V[iPpar_I[i]] = Si2NoP;
  }

  // Get back to SI units
  for (int iVar = 0; iVar < nVarCoupling; iVar++)
    No2Si_V[iVar] = 1.0 / Si2No_V[iVar];
}
//-------------------------------------------------------------------------

void FluidInterface::normalize_length() {
  // Normalization
  for (int i = 0; i < 3; i++) {
    dx_D[i] *= Si2NoL;
    lenPhy_D[i] *= Si2NoL;
    phyMin_D[i] *= Si2NoL;
    phyMax_D[i] *= Si2NoL;
    if (i >= nDimFluid) {
      // If MHD is 2D.
      phyMin_D[i] = 0;
      phyMax_D[i] = dx_D[i];
    }
  }
}
//-------------------------------------------------------------------------

/** Get nomal and pendicular vector to magnetic field */
void FluidInterface::calc_mag_base_vector(const double Bx, const double By,
                                          const double Bz,
                                          MDArray<double>& norm_DD) const {
  double inv;
  int Norm_, Perp1_, Perp2_, X_, Y_, Z_;
  Norm_ = 0;
  Perp1_ = 1;
  Perp2_ = 2;
  X_ = 0;
  Y_ = 1;
  Z_ = 2;

  inv = 1.0 / sqrt(Bx * Bx + By * By + Bz * Bz);
  norm_DD(Norm_, X_) = Bx * inv;
  norm_DD(Norm_, Y_) = By * inv;
  norm_DD(Norm_, Z_) = Bz * inv;

  if (norm_DD(Norm_, Z_) < 0.5) {
    norm_DD(Perp1_, X_) = norm_DD(Norm_, Y_);
    norm_DD(Perp1_, Y_) = -norm_DD(Norm_, X_);
    norm_DD(Perp1_, Z_) = 0.0;
    norm_DD(Perp2_, X_) = norm_DD(Norm_, Z_) * norm_DD(Norm_, X_);
    norm_DD(Perp2_, Y_) = norm_DD(Norm_, Z_) * norm_DD(Norm_, Y_);
    norm_DD(Perp2_, Z_) =
        -pow(norm_DD(Norm_, X_), 2) - pow(norm_DD(Norm_, Y_), 2);
  } else {
    norm_DD(Perp1_, X_) = 0.0;
    norm_DD(Perp1_, Y_) = norm_DD(Norm_, Z_);
    norm_DD(Perp1_, Z_) = -norm_DD(Norm_, Y_);
    norm_DD(Perp2_, X_) =
        -pow(norm_DD(Norm_, Y_), 2) - pow(norm_DD(Norm_, Z_), 2);
    norm_DD(Perp2_, Y_) = norm_DD(Norm_, Y_) * norm_DD(Norm_, X_);
    norm_DD(Perp2_, Z_) = norm_DD(Norm_, Z_) * norm_DD(Norm_, X_);
  }

  inv = 1.0 / sqrt(norm_DD(Perp1_, X_) * norm_DD(Perp1_, X_) +
                   norm_DD(Perp1_, Y_) * norm_DD(Perp1_, Y_) +
                   norm_DD(Perp1_, Z_) * norm_DD(Perp1_, Z_));
  norm_DD(Perp1_, X_) *= inv;
  norm_DD(Perp1_, Y_) *= inv;
  norm_DD(Perp1_, Z_) *= inv;

  inv = 1.0 / sqrt(norm_DD(Perp2_, X_) * norm_DD(Perp2_, X_) +
                   norm_DD(Perp2_, Y_) * norm_DD(Perp2_, Y_) +
                   norm_DD(Perp2_, Z_) * norm_DD(Perp2_, Z_));
  norm_DD(Perp2_, X_) *= inv;
  norm_DD(Perp2_, Y_) *= inv;
  norm_DD(Perp2_, Z_) *= inv;
}

// Data recived from SWMF coupler
void FluidInterface::read_from_GM(const int* const paramint,
                                  const double* const ParamRealRegion,
                                  const double* const ParamRealComm,
                                  const stringstream* const ss) {

  nDimFluid = paramint[0];
  nVarFluid = paramint[2];
  nFluid = paramint[3];
  nSpeciesFluid = paramint[4];

  // c++ index starts from 0. So, minus 1.
  iPe = paramint[5] - 1;
  iBx = paramint[6] - 1;
  iBy = iBx + 1;
  iBz = iBy + 1;

  iEx = paramint[7] - 1;
  iEy = iEx + 1;
  iEz = iEy + 1;

  nCellPerPatch = paramint[8];

  useElectronFluid = iEx > 1;

  if (useElectronFluid) {
    nIonFluid = -1; // Do not distinguish between electrons and ions.
    nIon = -1;
    nS = nFluid;
  } else {
    nIonFluid = nFluid;
    nIon = nFluid + nSpeciesFluid - 1; // Assuming one electron species.
    nS = nIon + 1;                     // + electron
  }

  useMultiFluid = nIonFluid > 1;
  useMultiSpecies = nSpeciesFluid > 1;

  nVarCoupling = nVarFluid + 3; // nVarFluid + (Jx, Jy, Jz)

  iRho_I.resize(nS);
  iRhoUx_I.resize(nS);
  iRhoUy_I.resize(nS);
  iRhoUz_I.resize(nS);
  iUx_I.resize(nS);
  iUy_I.resize(nS);
  iUz_I.resize(nS);
  iPpar_I.resize(nS);
  iP_I.resize(nS);

  int n = 9;
  if (useMultiSpecies) {
    // MultiSpecies. Densities of each species are known. Total velocity
    // and total pressure are known.
    iRhoTotal = paramint[n++] - 1;
    iRho_I[0] = iRhoTotal + 1;
    iRhoUx_I[0] = paramint[n++] - 1;
    iUx_I[0] = iRhoUx_I[0];
    iRhoUy_I[0] = iRhoUx_I[0] + 1;
    iUy_I[0] = iRhoUy_I[0];
    iRhoUz_I[0] = iRhoUx_I[0] + 2;
    iUz_I[0] = iRhoUz_I[0];
    iPpar_I[0] = paramint[n++] - 1;
    iP_I[0] = paramint[n++] - 1;

    for (int iIon = 1; iIon < nIon; ++iIon) {
      iRho_I[iIon] = iRho_I[0] + iIon;
      iRhoUx_I[iIon] = iRhoUx_I[0];
      iUx_I[iIon] = iUx_I[0];
      iRhoUy_I[iIon] = iRhoUy_I[0];
      iUy_I[iIon] = iUy_I[0];
      iRhoUz_I[iIon] = iRhoUz_I[0];
      iUz_I[iIon] = iUz_I[0];
      iPpar_I[iIon] = iPpar_I[0];
      iP_I[iIon] = iP_I[0];
    }
  } else {
    // Not multi-species
    for (int iFluid = 0; iFluid < nFluid; ++iFluid)
      iRho_I[iFluid] = paramint[n++] - 1;
    for (int iFluid = 0; iFluid < nFluid; ++iFluid) {
      iRhoUx_I[iFluid] = paramint[n++] - 1;
      iUx_I[iFluid] = iRhoUx_I[iFluid];
      iRhoUy_I[iFluid] = iRhoUx_I[iFluid] + 1;
      iUy_I[iFluid] = iRhoUy_I[iFluid];
      iRhoUz_I[iFluid] = iRhoUx_I[iFluid] + 2;
      iUz_I[iFluid] = iRhoUz_I[iFluid];
    }

    for (int iFluid = 0; iFluid < nFluid; ++iFluid)
      iPpar_I[iFluid] = paramint[n++] - 1;
    for (int iFluid = 0; iFluid < nFluid; ++iFluid)
      iP_I[iFluid] = paramint[n++] - 1;
  }

  int nVec = nFluid + 1;
  if (useElectronFluid)
    nVec++; // + E field.

  vecIdx_I.resize(nVec);

  for (int iVec = 0; iVec < nFluid; iVec++)
    vecIdx_I[iVec] = iRhoUx_I[iVec];
  vecIdx_I[nFluid] = iBx;
  if (useElectronFluid)
    vecIdx_I[nFluid + 1] = iEx;

  // See GM/BATSRUS/src/ModExtraVariables.f90.
  useAnisoP = iPpar_I[0] != 0;
  useMhdPe = iPe != 0;

  iJx = nVarFluid;
  iJy = iJx + 1;
  iJz = iJx + 2;

  Si2No_V.resize(nVarCoupling);
  No2Si_V.resize(nVarCoupling);

  for (int i = 0; i < nVarCoupling; i++)
    Si2No_V[i] = 1;

  n = 0;
  for (int i = 0; i < 3; i++) {
    phyMin_D[i] = ParamRealRegion[n++]; // Lmin
    phyMax_D[i] = phyMin_D[i] + ParamRealRegion[n++];
    dx_D[i] = ParamRealRegion[n++]; // dx
    lenPhy_D[i] = phyMax_D[i] - phyMin_D[i];
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      R_DD[i][j] = ParamRealRegion[n++];
    }
  }

  // Normalization parameters.
  Lnorm = ParamRealRegion[n++];
  Unorm = ParamRealRegion[n++];
  Mnorm = ParamRealRegion[n++];

  doRotate = false;
  double csmall = 1e-7;
  for (int i = 0; i < nDimFluid; i++)
    if (fabs(R_DD[i][i] - 1) > csmall) {
      doRotate = true;
    }

  for (int i = 0; i < 3; i++)
    nPhyCell_D[i] = (int)(lenPhy_D[i] / dx_D[i] + 0.5);

  QoQi_S.resize(nS);
  MoMi_S.resize(nS);

  /** Do not change the order of the following lines. */
  n = 0;
  if (useElectronFluid) {
    for (int i = 0; i < nS; ++i) {
      QoQi_S[i] = ParamRealComm[n++];
      MoMi_S[i] = ParamRealComm[n++];
    }
  } else {
    QoQi_S[0] = -1.0;
    for (int i = 1; i < nS; ++i) {
      QoQi_S[i] = ParamRealComm[n++];
      MoMi_S[i] = ParamRealComm[n++];
    }
  }

  // Electron pressure ratio: Pe/Ptotal
  PeRatio = ParamRealComm[n++];

  rPlanetSi = ParamRealComm[n++];
  MhdNo2SiL = ParamRealComm[n++];
  /** Do not change the order of above lines. */

  // Normalization units converted [SI] -> [cgs]
  if (Lnorm > 0) {
    Lnorm *= 100.0;
  } else {
    Lnorm = 1.0;
  }
  if (Unorm > 0) {
    Unorm *= 100.0;
  } else {
    Unorm = 1.0;
  }
  if (Mnorm > 0) {
    Mnorm *= 1000.0;
  } else {
    Mnorm = 1.0;
  }

  calc_normalized_units();
  normalize_length();

  if (useMultiFluid && !useMhdPe) {
    cout << printPrefix
         << " Use multi-fluid but do not use electron pressure. This case is "
            "not supported so far!!!"
         << endl;
    abort();
  }
}

/** print info for coupling */
void FluidInterface::print_info() const {

  if (myrank == 0) {
    cout << endl;
    cout.precision(15);
    cout << printPrefix << "Number of PIC species     = " << nS << endl;
    cout << printPrefix << "Total mass of all species = " << SumMass << endl;
    cout << printPrefix << "useMultiFluid  = " << (useMultiFluid ? "T" : "F")
         << endl;

    if (useMultiFluid)
      cout << printPrefix << "nFluid = " << nFluid << endl;
    cout << printPrefix << "useMultiSpecies = " << (useMultiSpecies ? "T" : "F")
         << endl;
    if (useMultiSpecies)
      cout << printPrefix << "nSpeciesFluid =" << nSpeciesFluid << endl;
    cout << printPrefix
         << "useElectronFluid = " << (useElectronFluid ? "T" : "F") << endl;
    for (int is = 0; is < nS; is++) {
      cout << printPrefix << "Q/Qi[" << is << "] = " << QoQi_S[is] << endl;
      cout << printPrefix << "M/Mi[" << is << "] = " << MoMi_S[is] << endl;
    }
    if (!useMhdPe)
      cout << printPrefix << "Pe/Ptotal = " << PeRatio << endl;
    cout << printPrefix
         << "============= Normalization factors =============" << endl;
    cout << printPrefix << "Mnorm    = " << Mnorm << endl;
    cout << printPrefix << "Unorm    = " << Unorm << endl;
    cout << printPrefix << "Qnorm    = " << Qnorm << endl;
    cout << printPrefix << "Lnorm    = " << Lnorm << endl;
    cout << printPrefix
         << "========== Unit conversion factors ==============" << endl;
    cout << printPrefix << "Si2NoRho = " << Si2NoRho << endl;
    cout << printPrefix << "Si2NoV   = " << Si2NoV << endl;
    cout << printPrefix << "Si2NoB   = " << Si2NoB << endl;
    cout << printPrefix << "Si2NoE   = " << Si2NoE << endl;
    cout << printPrefix << "Si2NoP   = " << Si2NoP << endl;
    cout << printPrefix << "Si2NoJ   = " << Si2NoJ << endl;
    cout << printPrefix << "Si2NoL   = " << Si2NoL << endl;
    cout << printPrefix
         << "===================================================" << endl;
  }
}

void FluidInterface::calc_fluid_state(const double* dataPIC_I,
                                      double* data_I) const {
  /* Input: dataPIC_I
     Output: data_I
     Function:
     dataPIC_I contains the information collected from EACH PIC species,
     and the E and B field. dataPIC_I have the normalized PIC units, and
     the vectors are in PIC coordinates. This function collects the
     MHD values, data_I, from dataPIC_I. data_I is in SI units and the
     vectors are in MHD coordinates.
   */

  double Rhoi, Pi;
  double Rhoe;
  double Rho;
  double PeXX, PeYY, PeZZ, PeXY, PeXZ, PeYZ, PiXX, PiYY, PiZZ, PiXY, PiXZ, PiYZ;
  double PtXX, PtYY, PtZZ, PtXY, PtXZ, PtYZ, BX, BY, BZ;
  double PitXX, PitYY, PitZZ, PitXY, PitXZ, PitYZ;
  double Mx, My, Mz;    // Total momentum.
  double Mix, Miy, Miz; // Momentum of i species.

  // (rho + 3*Moment + 6*p) = 10 variables per PIC species.
  const int nVarSpecies = 10;
  // For example, the rho index for PIC species iSpecies is
  // iRhoPIC+iSpecies*nVarSpecies
  const int iRhoPIC = 0, iMxPIC = 1, iMyPIC = 2, iMzPIC = 3, iPxxPIC = 4,
            iPyyPIC = 5, iPzzPIC = 6, iPxyPIC = 7, iPxzPIC = 8, iPyzPIC = 9;
  int iBxPIC = nS * nVarSpecies, iByPIC = iBxPIC + 1, iBzPIC = iByPIC + 1,
      iExPIC = iBzPIC + 1, iEyPIC = iExPIC + 1, iEzPIC = iEyPIC + 1;

  for (int i = 0; i < nVarFluid; i++) {
    // The variable for hyperbolic clean, is not known by iPIC3D,
    // but it is needed to be passed back. So data_I need to be
    // initilized.
    data_I[i] = 0;
  }

  BX = dataPIC_I[iBxPIC];
  BY = dataPIC_I[iByPIC];
  BZ = dataPIC_I[iBzPIC];

  data_I[iBx] = BX;
  data_I[iBy] = BY;
  data_I[iBz] = BZ;

  if (useElectronFluid) {

    data_I[iEx] = dataPIC_I[iExPIC];
    data_I[iEy] = dataPIC_I[iEyPIC];
    data_I[iEz] = dataPIC_I[iEzPIC];

    for (int iFluid = 0; iFluid < nFluid; ++iFluid) {
      Rhoi = dataPIC_I[iRhoPIC + iFluid * nVarSpecies];

      Mix = dataPIC_I[iMxPIC + iFluid * nVarSpecies];
      Miy = dataPIC_I[iMyPIC + iFluid * nVarSpecies];
      Miz = dataPIC_I[iMzPIC + iFluid * nVarSpecies];

      PiXX = dataPIC_I[iPxxPIC + iFluid * nVarSpecies];
      PiYY = dataPIC_I[iPyyPIC + iFluid * nVarSpecies];
      PiZZ = dataPIC_I[iPzzPIC + iFluid * nVarSpecies];
      Pi = (PiXX + PiYY + PiZZ) / 3.0;
      PiXY = dataPIC_I[iPxyPIC + iFluid * nVarSpecies];
      PiXZ = dataPIC_I[iPxzPIC + iFluid * nVarSpecies];
      PiYZ = dataPIC_I[iPyzPIC + iFluid * nVarSpecies];

      data_I[iRho_I[iFluid]] = Rhoi;

      data_I[iRhoUx_I[iFluid]] = Mix;
      data_I[iRhoUy_I[iFluid]] = Miy;
      data_I[iRhoUz_I[iFluid]] = Miz;

      data_I[iP_I[iFluid]] = Pi;
      if (useAnisoP) {
        data_I[iPpar_I[iFluid]] =
            (BX * PiXX * BX + BY * PiYY * BY + BZ * PiZZ * BZ +
             2.0 * BX * PiXY * BY + 2.0 * BX * PiXZ * BZ +
             2.0 * BY * PiYZ * BZ) /
            (BX * BX + BY * BY + BZ * BZ + 1e-40);
      }
    }
  } else {

    Rhoe = 0;
    PeXX = 0;
    PeYY = 0;
    PeZZ = 0;
    PeXY = 0;
    PeXZ = 0;
    PeYZ = 0;
    Mx = 0;
    My = 0;
    Mz = 0;
    for (int iSpecies = 0; iSpecies < nS; iSpecies++) {
      int iMHD = iSpecies;
      if (iMHD == 0) {
        // Electron

        Rhoe += dataPIC_I[iRhoPIC + iSpecies * nVarSpecies];
        Mx += dataPIC_I[iMxPIC + iSpecies * nVarSpecies];
        My += dataPIC_I[iMyPIC + iSpecies * nVarSpecies];
        Mz += dataPIC_I[iMzPIC + iSpecies * nVarSpecies];

        PeXX += dataPIC_I[iPxxPIC + iSpecies * nVarSpecies];
        PeYY += dataPIC_I[iPyyPIC + iSpecies * nVarSpecies];
        PeZZ += dataPIC_I[iPzzPIC + iSpecies * nVarSpecies];

        PeXY += dataPIC_I[iPxyPIC + iSpecies * nVarSpecies];
        PeXZ += dataPIC_I[iPxzPIC + iSpecies * nVarSpecies];
        PeYZ += dataPIC_I[iPyzPIC + iSpecies * nVarSpecies];
      }
    } // iSpecies

    if (useMhdPe)
      data_I[iPe] = (PeXX + PeYY + PeZZ) / 3.0;

    Rho = Rhoe;
    PtXX = PeXX;
    PtYY = PeYY;
    PtZZ = PeZZ;
    PtXY = PeXY;
    PtXZ = PeXZ;
    PtYZ = PeYZ;

    for (int iSpecies = 0; iSpecies < nS; ++iSpecies) {
      int iIon;
      iIon = iSpecies - 1; // The first species is electron;

      if (iIon >= 0) {
        Rhoi = dataPIC_I[iRhoPIC + iSpecies * nVarSpecies];
        Mix = dataPIC_I[iMxPIC + iSpecies * nVarSpecies];
        Miy = dataPIC_I[iMyPIC + iSpecies * nVarSpecies];
        Miz = dataPIC_I[iMzPIC + iSpecies * nVarSpecies];

        PiXX = dataPIC_I[iPxxPIC + iSpecies * nVarSpecies];
        PiYY = dataPIC_I[iPyyPIC + iSpecies * nVarSpecies];
        PiZZ = dataPIC_I[iPzzPIC + iSpecies * nVarSpecies];
        PiXY = dataPIC_I[iPxyPIC + iSpecies * nVarSpecies];
        PiXZ = dataPIC_I[iPxzPIC + iSpecies * nVarSpecies];
        PiYZ = dataPIC_I[iPyzPIC + iSpecies * nVarSpecies];

        // Sum to total density/pressure.
        Rho += Rhoi;
        Mx += Mix;
        My += Miy;
        Mz += Miz;

        PtXX += PiXX;
        PtYY += PiYY;
        PtZZ += PiZZ;
        PtXY += PiXY;
        PtXZ += PiXZ;
        PtYZ += PiYZ;

        // Density
        if (useMultiFluid || useMultiSpecies) {
          data_I[iRho_I[iIon]] += Rhoi;
        }

        // Pressure.
        if (useMultiFluid) {
          // ONLY works for iso pressure so far!!!!!
          data_I[iP_I[iIon]] += (PiXX + PiYY + PiZZ) / 3;
          if (useAnisoP) {
            cout << printPrefix
                 << "Multi-fluid model can not work with aniso pressure "
                    "now!!"
                 << endl;
            abort();
          }
        }

        // Momentum.
        if (useMultiFluid) {
          data_I[iRhoUx_I[iIon]] += Mix;
          data_I[iRhoUy_I[iIon]] += Miy;
          data_I[iRhoUz_I[iIon]] += Miz;
        }
      } // if(iIon > 0)
    }   // iSpecies

    if (!(useMultiFluid || useMultiSpecies)) {
      int iIon = 0;
      // Only one ion species.Rho = Rhoi + Rhoe
      data_I[iRho_I[iIon]] = Rho;
    }

    // Do not includes electron density. The total density passed
    // in MHD side is useless, so it doesnot matter whether Rho include
    // electron or not. -- Yuxi
    if (useMultiSpecies)
      data_I[iRhoTotal] = Rho - Rhoe;

    // Momentum
    if (!useMultiFluid && !useElectronFluid) {
      // Include electron momentum.
      data_I[iRhoUx_I[0]] = Mx;
      data_I[iRhoUy_I[0]] = My;
      data_I[iRhoUz_I[0]] = Mz;
    }

    // Sum of ion pressure.
    PitXX = PtXX - PeXX;
    PitYY = PtYY - PeYY;
    PitZZ = PtZZ - PeZZ;
    PitXY = PtXY - PeXY;
    PitXZ = PtXZ - PeXZ;
    PitYZ = PtYZ - PeYZ;

    // Pressure
    if (!useMultiFluid && !useElectronFluid) {
      // AnisoP
      if (useAnisoP) {
        if (useMhdPe)
          data_I[iPpar_I[0]] = (BX * PitXX * BX + BY * PitYY * BY +
                                BZ * PitZZ * BZ + 2.0 * BX * PitXY * BY +
                                2.0 * BX * PitXZ * BZ + 2.0 * BY * PitYZ * BZ) /
                               (BX * BX + BY * BY + BZ * BZ + 1e-40);
        else
          data_I[iPpar_I[0]] = (BX * PtXX * BX + BY * PtYY * BY +
                                BZ * PtZZ * BZ + 2.0 * BX * PtXY * BY +
                                2.0 * BX * PtXZ * BZ + 2.0 * BY * PtYZ * BZ) /
                               (BX * BX + BY * BY + BZ * BZ + 1e-40);
      } // useAnisoP

      // Isotropic Pressure.
      if (useMhdPe)
        data_I[iP_I[0]] = (PitXX + PitYY + PitZZ) / 3;
      else
        data_I[iP_I[0]] = (PtXX + PtYY + PtZZ) / 3;
    }
  } // useElectronFluid is true or false

  // Convert to SI units
  for (int iVar = 0; iVar < nVarFluid; ++iVar) {
    data_I[iVar] *= get_No2Si_V(iVar);
  }

  // Convert the vectors from PIC coordinates to MHD coordinates.
  double mhd_D[3], pic_D[3];
  for (int iVec = 0; iVec < vecIdx_I.size(); iVec++) {
    int idx0;
    idx0 = vecIdx_I[iVec];
    for (int iVar = idx0; iVar < idx0 + nDimFluid; iVar++)
      pic_D[iVar - idx0] = data_I[iVar];
    pic_to_Mhd_Vec(pic_D, mhd_D, true);
    for (int iVar = idx0; iVar < idx0 + nDimFluid; iVar++)
      data_I[iVar] = mhd_D[iVar - idx0];
  } // iVec
}
void FluidInterface::pic_to_Mhd_Vec(double const* vecIn_D, double* vecOut_D,
                                    bool isZeroOrigin) const {
  /** 1) Change a vector in coupling PIC coordinates to MHD coordinates.
      If not isZeroOrigin, then shifting the origin of the coupling
      PIC coordinates to INxRange_I.
   **/

  for (int iDim = 0; iDim < nDimFluid; iDim++) {
    vecOut_D[iDim] = 0;
    for (int jDim = 0; jDim < nDimFluid; jDim++) {
      vecOut_D[iDim] += R_DD[iDim][jDim] * vecIn_D[jDim];
    }
    if (!isZeroOrigin)
      vecOut_D[iDim] += phyMin_D[iDim];
  }
}

void FluidInterface::mhd_to_Pic_Vec(double const* vecIn_D, double* vecOut_D,
                                    bool isZeroOrigin) const {
  /** 1) Change a vector in MHD coordinates to coupling PIC coordinates.
      If not isZeroOrigin, then shifting the origin of the coupling
      PIC coordinates to phyMin_D.
  **/

  double vec_D[3];
  if (!isZeroOrigin) {
    for (int iDim = 0; iDim < nDimFluid; iDim++)
      vec_D[iDim] = vecIn_D[iDim] - phyMin_D[iDim];
  } else {
    for (int iDim = 0; iDim < nDimFluid; iDim++)
      vec_D[iDim] = vecIn_D[iDim];
  }

  for (int iDim = 0; iDim < nDimFluid; iDim++) {
    vecOut_D[iDim] = 0;
    for (int jDim = 0; jDim < nDimFluid; jDim++) {
      vecOut_D[iDim] += R_DD[jDim][iDim] * vec_D[jDim];
    }
  } // iDim
}

void FluidInterface::update_nodeFluid(const MultiFabFLEKS& nodeIn,
                                      const double dt) {

  double No2MhdNoL = No2SiL * (1.0 / MhdNo2SiL);
  double dtSI = dt * get_No2SiT();

  nodeFluid.setVal(0);

  for (MFIter mfi(nodeFluid); mfi.isValid(); ++mfi) {
    const Box& box = mfi.fabbox();
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    const Array4<Real>& arr = nodeFluid[mfi].array();
    const Array4<const Real>& arrIn = nodeIn[mfi].array();

    for (int k = lo.z; k <= hi.z; ++k)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int i = lo.x; i <= hi.x; ++i) {

          Real z = geom.CellCenter(k, iz_) * No2MhdNoL;
          Real y = geom.CellCenter(j, iy_) * No2MhdNoL;
          Real x = geom.CellCenter(i, ix_) * No2MhdNoL;

          if (useMultiSpecies) {
            // double Rhot = 0;
            // for (int iIon = 0; iIon < nIon; ++iIon) {
            //   // Rho = sum(Rhoi) + Rhoe;
            //   Rhot += arr(i, j, k, iRho_I[iIon]) *
            //           (1 + MoMi_S[0] / MoMi_S[iIon + 1]);
            // } // iIon

            // arr(i, j, k, iUx_I[0]) /= Rhot;
            // arr(i, j, k, iUy_I[0]) /= Rhot;
            // arr(i, j, k, iUz_I[0]) /= Rhot;
          } else {

            Real x0 = 5, y0 = 5, z0 = 0, r0 = 2;
            x -= x0;
            y -= y0;
            z -= z0;
            Real r = sqrt(x * x + y * y + z * z);
            Real ratio = exp(-r / r0) * 0.1;
            for (int iFluid = 0; iFluid < nFluid; ++iFluid) {

              arr(i, j, k, iRho_I[iFluid]) =
                  ratio * arrIn(i, j, k, iRho_I[iFluid]);

              arr(i, j, k, iUx_I[iFluid]) = arrIn(i, j, k, iUx_I[iFluid]);
              arr(i, j, k, iUy_I[iFluid]) = arrIn(i, j, k, iUy_I[iFluid]);
              arr(i, j, k, iUz_I[iFluid]) = arrIn(i, j, k, iUz_I[iFluid]);

              arr(i, j, k, iP_I[iFluid]) = ratio * arrIn(i, j, k, iP_I[iFluid]);

            } // iFluid
          }
        }
  }
}
