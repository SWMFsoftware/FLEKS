#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include "FluidInterface.h"
#include "GridUtility.h"
#include "Utility.h"

using namespace amrex;
using namespace std;

void FluidInterface::post_process_param() {
  if (initFromSWMF)
    return;

  nDimFluid = (Geom(0).Domain().length(iz_) == 1) ? 2 : 3;

  const Real protonMassPerChargeSI = cProtonMassSI / cUnitChargeSI;

  mNormSI = 1e7 * lNormSI * pow(protonMassPerChargeSI * ScalingFactor, 2);

  {
    iRho_I.resize(nS);
    iRhoUx_I.resize(nS);
    iRhoUy_I.resize(nS);
    iRhoUz_I.resize(nS);
    iUx_I.resize(nS);
    iUy_I.resize(nS);
    iUz_I.resize(nS);
    iPpar_I.resize(nS);
    iP_I.resize(nS);

    int idx = 0;
    for (int i = 0; i < nS; i++) {
      iRho_I[i] = idx++;
      varNames.push_back("rho" + to_string(i));
      iUx_I[i] = idx++;
      varNames.push_back("ux" + to_string(i));
      iUy_I[i] = idx++;
      varNames.push_back("uy" + to_string(i));
      iUz_I[i] = idx++;
      varNames.push_back("uz" + to_string(i));
      iPpar_I[i] = idx++;
      varNames.push_back("ppar" + to_string(i));
      iP_I[i] = idx++;
      varNames.push_back("p" + to_string(i));
    }
    iBx = idx++;
    varNames.push_back("bx");
    iBy = idx++;
    varNames.push_back("by");
    iBz = idx++;
    varNames.push_back("bz");
    iEx = idx++;
    varNames.push_back("ex");
    iEy = idx++;
    varNames.push_back("ey");
    iEz = idx++;
    varNames.push_back("ez");
    if (nVarCoupling > nVarFluid) {
      iJx = nVarFluid;
      varNames.push_back("jx");
      iJy = iJx + 1;
      varNames.push_back("jy");
      iJz = iJx + 2;
      varNames.push_back("jz");
    }
    iRhoUx_I = iUx_I;
    iRhoUy_I = iUy_I;
    iRhoUz_I = iUz_I;
  }
  calc_normalized_units();
}

FluidInterface::FluidInterface(Geometry const& gm, AmrInfo const& amrInfo,
                               int nGst, int id, std::string tag,
                               const amrex::Vector<int>& iParam,
                               const amrex::Vector<double>& norm,
                               const amrex::Vector<double>& paramComm)
    : Grid(gm, amrInfo, nGst, id, tag) {

  initFromSWMF = true;

  if (iParam.empty() || norm.empty() || paramComm.empty())
    amrex::Abort("Error: one of the input vector is empty!\n");

  nDimFluid = (Geom(0).Domain().length(iz_) == 1) ? 2 : 3;
  nVarFluid = iParam[2];
  nFluid = iParam[3];
  nSpeciesFluid = iParam[4];

  // c++ index starts from 0. So, minus 1.
  iPe = iParam[5] - 1;
  iBx = iParam[6] - 1;
  iBy = iBx + 1;
  iBz = iBy + 1;

  iEx = iParam[7] - 1;
  iEy = iEx + 1;
  iEz = iEy + 1;

  nCellPerPatch = iParam[8];

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

  nVarCoupling = nVarFluid;
  if (myType == PICFluid)
    nVarCoupling += 3; // nVarFluid + (Jx, Jy, Jz)

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
    iRhoTotal = iParam[n++] - 1;
    iRho_I[0] = iRhoTotal + 1;
    iRhoUx_I[0] = iParam[n++] - 1;
    iUx_I[0] = iRhoUx_I[0];
    iRhoUy_I[0] = iRhoUx_I[0] + 1;
    iUy_I[0] = iRhoUy_I[0];
    iRhoUz_I[0] = iRhoUx_I[0] + 2;
    iUz_I[0] = iRhoUz_I[0];
    iPpar_I[0] = iParam[n++] - 1;
    iP_I[0] = iParam[n++] - 1;

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
      iRho_I[iFluid] = iParam[n++] - 1;
    for (int iFluid = 0; iFluid < nFluid; ++iFluid) {
      iRhoUx_I[iFluid] = iParam[n++] - 1;
      iUx_I[iFluid] = iRhoUx_I[iFluid];
      iRhoUy_I[iFluid] = iRhoUx_I[iFluid] + 1;
      iUy_I[iFluid] = iRhoUy_I[iFluid];
      iRhoUz_I[iFluid] = iRhoUx_I[iFluid] + 2;
      iUz_I[iFluid] = iRhoUz_I[iFluid];
    }

    for (int iFluid = 0; iFluid < nFluid; ++iFluid)
      iPpar_I[iFluid] = iParam[n++] - 1;
    for (int iFluid = 0; iFluid < nFluid; ++iFluid)
      iP_I[iFluid] = iParam[n++] - 1;
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

  n = 0;
  // Normalization parameters.
  lNormSI = norm[n++];
  uNormSI = norm[n++];
  mNormSI = norm[n++];
  ScalingFactor = norm[n++];

  QoQi_S.resize(nS);
  MoMi_S.resize(nS);

  /** Do not change the order of the following lines. */
  n = 0;
  if (useElectronFluid) {
    for (int i = 0; i < nS; ++i) {
      QoQi_S[i] = paramComm[n++];
      MoMi_S[i] = paramComm[n++];
    }
  } else {
    QoQi_S[0] = -1.0;
    for (int i = 1; i < nS; ++i) {
      QoQi_S[i] = paramComm[n++];
      MoMi_S[i] = paramComm[n++];
    }
  }

  // Electron pressure ratio: Pe/Ptotal
  PeRatio = paramComm[n++];

  rPlanetSi = paramComm[n++];
  MhdNo2SiL = paramComm[n++];
  /** Do not change the order of above lines. */

  calc_normalized_units();

  if (useMultiFluid && !useMhdPe) {
    cout << printPrefix
         << " Use multi-fluid but do not use electron pressure. This "
            "case is "
            "not supported so far!!!"
         << endl;
    abort();
  }
}

//==========================================================
void FluidInterface::read_param(const std::string& command, ReadParam& param) {
  if (command == "#NORMALIZATION") {
    param.read_var("lNorm", lNormSI);
    param.read_var("uNorm", uNormSI);
  } else if (command == "#SCALINGFACTOR") {
    param.read_var("scaling", ScalingFactor);
  } else if (command == "#BODYSIZE") {
    param.read_var("scaling", rPlanetSi);
  } else if (command == "#PLASMA") {
    param.read_var("nS", nS);
    QoQi_S.resize(nS);
    MoMi_S.resize(nS);
    for (int i = 0; i < nS; i++) {
      param.read_var("mass", MoMi_S[i]);
      param.read_var("charge", QoQi_S[i]);
    }
    nFluid = nS;
    // (rho, vx, vy, vz, ppar, p)*nFluid + B + E
    nVarFluid = 6 * nFluid + 3 + 3;
    nVarCoupling = nVarFluid + 3; // nVarFluid + (Jx, Jy, Jz)
  } else if (command == "#UNIFORMSTATE") {
    if (nS <= 0) {
      amrex::Abort("Error: number of species <=0! Use #PLASMA command "
                   "to set plasma species information.");
    }
    double tmp;
    for (int i = 0; i < nS; i++) {
      double rho;
      param.read_var("rho", rho);
      uniformState.push_back(rho);
      param.read_var("ux", tmp);
      uniformState.push_back(tmp * rho);
      param.read_var("uy", tmp);
      uniformState.push_back(tmp * rho);
      param.read_var("uz", tmp);
      uniformState.push_back(tmp * rho);
      param.read_var("ppar", tmp);
      uniformState.push_back(tmp);
      param.read_var("p", tmp);
      uniformState.push_back(tmp);
    }
    param.read_var("bx", tmp);
    uniformState.push_back(tmp);
    param.read_var("by", tmp);
    uniformState.push_back(tmp);
    param.read_var("bz", tmp);
    uniformState.push_back(tmp);
    param.read_var("ex", tmp);
    uniformState.push_back(tmp);
    param.read_var("ey", tmp);
    uniformState.push_back(tmp);
    param.read_var("ez", tmp);
    uniformState.push_back(tmp);
  }
}

void FluidInterface::regrid(const amrex::BoxArray& centerBAIn,
                            const amrex::DistributionMapping& dmIn) {
  std::string nameFunc = "FluidInterface::regrid";

  // Why need 'isGridInitialized'? See the explaination in
  // Domain::regrid().
  if (centerBAIn == cGrid && isGridInitialized) {
    // The interface grid does not change.
    return;
  }

  isGridEmpty = centerBAIn.empty();

  cGrid = centerBAIn;
  nGrid = convert(cGrid, amrex::IntVect{ AMREX_D_DECL(1, 1, 1) });

  if (!cGrid.empty()) {
    // This method will call MakeNewLevelFromScratch() and
    // PostProcessBaseGrids()
    InitFromScratch(0.0);
    SetDistributionMap(0, dmIn);
  }

  const bool doCopy = true;
  distribute_FabArray(nodeFluid, nGrid, DistributionMap(0), nVarCoupling, nGst,
                      doCopy);
  distribute_FabArray(centerB, cGrid, DistributionMap(0), nDimMax, nGst,
                      doCopy);

  isGridInitialized = true;
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

  const Real* dx = Geom(0).CellSize();
  const auto plo = Geom(0).ProbLo();

  // Global NODE box.
  const Box gbx = convert(Geom(0).Domain(), { 1, 1, 1 });

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
              int idx[nDimMax] = { i, j, k };
              for (int iDim = 0; iDim < nDimMax; iDim++) {
                if (Geom(0).isPeriodic(iDim)) {
                  idx[iDim] = shift_periodic_index(
                      idx[iDim], gbx.smallEnd(iDim), gbx.bigEnd(iDim));
                }
              }

              for (int iDim = 0; iDim < get_fluid_dimension(); iDim++) {
                pos_DI[nCount++] = (idx[iDim] * dx[iDim] + plo[iDim]) * no2siL;
              }
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
  loop_through_node("loc", pos_DI);
}

void FluidInterface::set_node_fluid(const double* const data,
                                    const int* const index) {
  if (isGridEmpty)
    return;

  loop_through_node("fill", nullptr, data, index);

  calc_current();
  normalize_fluid_variables();
  convert_moment_to_velocity();

  if (nodeFluid.size() > iJx) {
    MultiFab currentMF(nodeFluid, make_alias, iJx, nDimMax);
  }

  // save_amrex_file();
}

void FluidInterface::set_node_fluid() {
  if (isGridEmpty)
    return;

  if (uniformState.empty()) {
    Abort("Error: use #UNIFORMSTATE command to set the initail state.");
  }

  for (int i = 0; i < nVarFluid; i++)
    nodeFluid.setVal(uniformState[i], i, 1, nodeFluid.nGrow());

  calc_current();
  normalize_fluid_variables();
  convert_moment_to_velocity();

  if (nodeFluid.size() > iJx) {
    MultiFab currentMF(nodeFluid, make_alias, iJx, nDimMax);
  }

  // save_amrex_file();
}

void FluidInterface::set_node_fluid(const FluidInterface& other) {
  if (isGridEmpty)
    return;

  MultiFab::Copy(nodeFluid, other.nodeFluid, 0, 0, nodeFluid.nComp(),
                 nodeFluid.nGrow());

  calc_current();
  normalize_fluid_variables();
  convert_moment_to_velocity();

  if (nodeFluid.size() > iJx) {
    MultiFab currentMF(nodeFluid, make_alias, iJx, nDimMax);
  }

  // save_amrex_file();
}

void FluidInterface::calc_current() {
  if (isGridEmpty)
    return;

  if (nVarCoupling == nVarFluid)
    return;

  // All centerB, including all ghost cells are accurate.
  average_node_to_cellcenter(centerB, 0, nodeFluid, iBx, centerB.nComp(),
                             centerB.nGrow());

  // currentMF is just an alias of current components of nodeFluid.
  MultiFab currentMF(nodeFluid, make_alias, iJx, nDimMax);

  // The outmost layer of currentMF can not be calculated from centerB
  curl_center_to_node(centerB, currentMF, Geom(0).InvCellSize());
  currentMF.mult(1.0 / (get_No2SiL() * fourPI * 1e-7), currentMF.nGrow());

  currentMF.FillBoundary(Geom(0).periodicity());

  /*
  Q: The outmost layer of currentMF is not accurate. Why not use
  apply_float_boundary to fill in first-order estimation?

  A: If the whole domain is just ONE block, it will work. Otherwise, it
  will not. For example, For a 2D simulation domain of 6x3 with 2 blocks.
  In the x-direction, block-1 convers cell 0 (c+0) to cell 2 (c+2), and
  block-2 covers c+3 to c+5. The node (n+5, n-1) is the corner ghost node
  for the block-1, and it is the face ghost node for the block-2. On
  block-2, this node can be calculated from curl_center_to_node(centerB,
  currentMF....). However, block-1 does not know how to calculate it, and
  FillBoundary will also NOT copy this node from block-2 to block-1,
  because this node is a boundary node and it is not covered by any
  physical node.

  So, we should keep in mind, the variables that are directly received
  from the MHD side, such as the magnetic fields, are accurate on all
  ghost nodes, but the current releated variables (current, plasma
  velocities, and electric field) are unknown at the outmost boundary node
  layer.
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
  SetDistributionMap(0, dmIn);

  redistribute_FabArray(nodeFluid, DistributionMap(0)); // false?
  redistribute_FabArray(centerB, DistributionMap(0));   // false?
}

//-------------------------------------------------------------------------

void FluidInterface::calc_normalized_units() {

  // Normalization units converted [SI] -> [cgs]
  if (lNormSI > 0) {
    Lnorm = 100.0 * lNormSI;
  } else {
    Lnorm = 1.0;
  }
  if (uNormSI > 0) {
    Unorm = 100.0 * uNormSI;
  } else {
    Unorm = 1.0;
  }
  if (mNormSI > 0) {
    Mnorm = 1000.0 * mNormSI;
  } else {
    Mnorm = 1.0;
  }

  // normalization variables
  double RHOnorm, Bnorm, Jnorm, Pnorm;

  RHOnorm = Mnorm / (Lnorm * Lnorm * Lnorm);

  /* How to calculate Bnorm?
     1. Method 1
      In CGS unit, we have:
      1) [B] = [E]
      2) div(E) ~ rhoq, where rhoq is the charge density -> [E] =
    [rhoq]*[L] 3) moment equation: d(rho*u)/dt ~ rhoq*u x B/c     ->
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

  Si2No_V.resize(nVarCoupling);
  No2Si_V.resize(nVarCoupling);

  for (int i = 0; i < nVarCoupling; i++)
    Si2No_V[i] = 1;

  Si2No_V[iBx] = Si2NoB;
  Si2No_V[iBy] = Si2NoB;
  Si2No_V[iBz] = Si2NoB;

  if (Si2No_V.size() > iJz) {
    Si2No_V[iJx] = Si2NoJ;
    Si2No_V[iJy] = Si2NoJ;
    Si2No_V[iJz] = Si2NoJ;
  }

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

/** print info for coupling */
void FluidInterface::print_info() const {

  if (ParallelDescriptor::MyProc() == 0) {
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

    cout << printPrefix << "useResistivity = " << (useResist ? "T" : "F")
         << "\n etaSI (magnetic diffusivity with unit m^2/s) = " << etaSI
         << "\n etaNO (magnetic diffusivity)                 = " << etaNO
         << endl;
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
  double pic_D[3];
  for (int iVec = 0; iVec < vecIdx_I.size(); iVec++) {
    int idx0;
    idx0 = vecIdx_I[iVec];
    for (int iVar = idx0; iVar < idx0 + nDimFluid; iVar++)
      pic_D[iVar - idx0] = data_I[iVar];

    for (int iVar = idx0; iVar < idx0 + nDimFluid; iVar++)
      data_I[iVar] = pic_D[iVar - idx0];
  } // iVec
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

          Real z = Geom(0).CellCenter(k, iz_) * No2MhdNoL;
          Real y = Geom(0).CellCenter(j, iy_) * No2MhdNoL;
          Real x = Geom(0).CellCenter(i, ix_) * No2MhdNoL;

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

void FluidInterface::save_amrex_file() {
  string filename = component + "/plots/" + tag;
  Print() << "Writing FluidInterface file " << filename;

  for (int i = 0; i < nodeFluid.nComp(); i++) {
    Real no2out = No2Si_V[i];
    nodeFluid.mult(no2out, i, 1, nodeFluid.nGrow());
  }

  WriteSingleLevelPlotfile(filename, nodeFluid, varNames, Geom(0), 0, 0);

  for (int i = 0; i < nodeFluid.nComp(); i++) {
    Real out2no = Si2No_V[i];
    nodeFluid.mult(out2no, i, 1, nodeFluid.nGrow());
  }
}

void FluidInterface::get_for_points(const int nDim, const int nPoint,
                                    const double* const xyz_I,
                                    double* const data_I, const int nVar,
                                    const double coef) {
  std::string nameFunc = "FI::get_for_points";

  const RealBox& range = Geom(0).ProbDomain();
  for (int iPoint = 0; iPoint < nPoint; iPoint++) {
    double xyz[3];
    for (int iDim = 0; iDim < nDim; iDim++) {
      xyz[iDim] = xyz_I[iPoint * nDim + iDim] * get_Si2NoL();
    }

    const Real xp = xyz[0];
    const Real yp = (nDim > 1) ? xyz[1] : 0.0;
    const Real zp = (nDim > 2) ? xyz[2] : 0.0;

    // Check if this point is inside this FLEKS domain.
    if (!range.contains(RealVect(xp, yp, zp), 1e-10))
      continue;

    const int iStart = iPoint * nVar;
    for (int iVar = 0; iVar <= nVar; iVar++) {
      data_I[iStart + iVar] =
          get_value_at_loc(nodeFluid, Geom(0), xp, yp, zp, iVar) * coef;
    }
  }
}
