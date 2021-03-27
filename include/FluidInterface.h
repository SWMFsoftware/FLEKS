#ifndef _FLUIDINTERFACE_H_
#define _FLUIDINTERFACE_H_

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>

#include "BC.h"
#include "Constants.h"
#include "FluidPicInterface.h"
#include "MDArray.h"
#include "ReadParam.h"
#include "Utility.h"
#include "Writer.h"

class FluidInterface {

protected:
  static const int nDimMax = 3;
  int myrank;

  int nCellPerPatch;

  int nDim; // number of dimentions

  // Min and Max of the physical domain in normalized PIC units.
  double phyMin_D[3], phyMax_D[3];

  // The length of the computational domain in normalized PIC units, including
  // the ghost cell layers.
  double lenPhy_D[3];

  // Cell Size
  double dx_D[3];

  // Rotation matrix.
  double R_DD[3][3];

  bool doRotate;

  // Number of nodes in each direction, including the ghost cell layers.
  int nPhyCell_D[nDimMax];

  // Number of variables passing between MHD and PIC.
  int nVarFluid;

  // Number of fluid at the MHD side. One 'fluid' has its own density,
  // velocity and pressure. Electron can be one fluid.
  int nFluid;

  // Number of ion fluid at the MHD side.
  int nIonFluid;

  // Number of species at the MHD side. One 'species' only has its own density.
  int nSpecies;

  // Total number of ion/electron species exit in the fluid code.
  int nIon;

  int nVarCoupling;

  bool useMultiSpecies, useMultiFluid, useElectronFluid;

  static const int NG = 1; // number of ghost cell

  bool useAnisoP; // Use anisotripic pressure

  bool useMhdPe;

  double rPlanetSi;

  double* Si2No_V; // array storing unit conversion factors
  double* No2Si_V; // array storing inverse unit conversion factors
  double Si2NoM, Si2NoV, Si2NoRho, Si2NoB, Si2NoP, Si2NoJ, Si2NoL, Si2NoE;
  double No2SiV, No2SiL;
  double MhdNo2SiL; // Length in BATSRUS normalized unit -> Si
  double Lnorm, Unorm, Mnorm,
      Qnorm; // normalization units for length, velocity, mass and charge
             // Normalized q/m ==1 for proton in CGS units

  //-------------------------------------------------------------------
  // nSIn is the number of species exists at the MHD side. PIC may use
  // two or more species to represent one MHD species. nS is the nuber of the
  // PIC species.
  long nS;         // number of particle species
  int nSIn;        // number of particle species before splitting.
  double* MoMi_S;  // masses for the particles species
  double* QoQi_S;  // charge for each particle species
  double* MoMi0_S; // masses for the particles species before splitting
  double* QoQi0_S; // charge for each particle species before splitting
  //-------------------------------------------------------------------

  double PeRatio; // temperature ratio for electrons: PeRatio = Pe/Ptotal
  double SumMass; // Sum of masses of each particle species
public:
  ReadParam readParam;

protected:
  // Variables for IDL format output.

public:
  // These variables are also used in PSKOutput.h
  int *iRho_I, *iRhoUx_I, *iRhoUy_I, *iRhoUz_I, iBx, iBy, iBz, iEx, iEy, iEz,
      iPe, *iPpar_I, *iP_I, iJx, iJy, iJz, *iUx_I, *iUy_I, *iUz_I, iRhoTotal;

  // At most 10 vectors are supported during the coupling.
  static const int nVecMax = 10;
  int vecIdxStart_I[nVecMax], nVec;

private:
  // ------Grid info----------
  amrex::DistributionMapping dm;
  amrex::Geometry geom;
  amrex::BoxArray centerBA;
  amrex::BoxArray nodeBA;
  //------------------------

  amrex::MultiFab nodeFluid;
  amrex::MultiFab centerB;

  double invSumMass;

  int nGst;

  int domainID;
  std::string domainName;
  std::string printPrefix;

  bool isGridInitialized = false;
  bool isGridEmpty = false;

public:
  FluidInterface() {}
  ~FluidInterface() = default;
  void init(int domainIDIn);
  void receive_info_from_gm(const int* const paramInt,
                            const double* const gridDim,
                            const double* const paramDouble,
                            const std::string& paramString);

  void set_geom(const int nGstIn, const amrex::Geometry& geomIn);

  void regrid(const amrex::BoxArray& centerBAIn,
              const amrex::DistributionMapping& dmIn);

  int count_couple_node_number();

  int loop_through_node(std::string action, double* const pos_DI = nullptr,
                        const double* const data = nullptr,
                        const int* const index = nullptr);

  void get_couple_node_loc(double* const pos_DI);

  void set_couple_node_value(const double* const data, const int* const index);

  void calc_current();

  void normalize_fluid_variables();

  void convert_moment_to_velocity();

  void set_plasma_charge_and_mass(amrex::Real qomEl);

  void load_balance(const amrex::DistributionMapping& dmIn);

  void InitData();

  void ReNormLength();

  void ReadFromGMinit(const int* const paramint,
                      const double* const ParamRealRegion,
                      const double* const ParamRealComm,
                      const std::stringstream* const ss);

  void checkParam();

  /** Get nomal and pendicular vector to magnetic field */
  void MagneticBaseVectors(const double Bx, const double By, const double Bz,
                           MDArray<double>& norm_DD) const;

  void CalcFluidState(const double* dataPIC_I, double* dataFluid_I) const;

  void mhd_to_Pic_Vec(double const* vecIn_D, double* vecOut_D,
                      bool isZeroOrigin = false) const;
  void pic_to_Mhd_Vec(double const* vecIn_D, double* vecOut_D,
                      bool isZeroOrigin = false) const;
  void PrintFluidPicInterface();

  int get_nCellPerPatch() const { return nCellPerPatch; }
  /** Use anisotropisc pressure when seting upt the particle distribution */
  bool getUseAnisoP() const { return (useAnisoP); }
  // void setAnisoP(bool useAnisoPIn){useAnisoP = useAnisoPIn;}

  bool get_useElectronFluid() const { return useElectronFluid; }

  /** Get convertion factor to from IPIC3D internal units */
  inline double getNo2Si_V(int idx) const { return (No2Si_V[idx]); }

  /** Get convertion factor to from IPIC3D internal units */
  inline double getSi2No_V(int idx) const { return (Si2No_V[idx]); }

  // The begining 'physical' point of this IPIC region. Assume there is one
  // layer PIC ghost cell.
  double getphyMin(int i) const { return phyMin_D[i]; }
  double getphyMax(int i) const { return phyMax_D[i]; }

  int getnDim() const { return (nDim); }

  int get_nS() const { return nS; }
  int get_nFluid() const { return (nFluid); }

  double getSi2NoL() const { return (Si2NoL); }
  double getSi2NoT() const { return Si2NoL / Si2NoV; }
  double getNo2SiL() const { return (No2SiL); }
  double getNo2SiRho() const { return (1. / Si2NoRho); }
  double getNo2SiV() const { return (1. / Si2NoV); }
  double getNo2SiB() const { return (1. / Si2NoB); }
  double getNo2SiP() const { return (1. / Si2NoP); }
  double getNo2SiJ() const { return (1. / Si2NoJ); }
  double getNo2SiT() const { return Si2NoV / Si2NoL; }

  double getMiSpecies(int i) const { return MoMi_S[i]; };
  double getQiSpecies(int i) const { return QoQi_S[i]; };
  double get_qom(int is) const { return QoQi_S[is] / MoMi_S[is]; }

  double get_cLight_SI() const { return Unorm / 100; /*Unorm is in cgs unit*/ };
  // return planet radius in SI unit.
  inline double get_rPlanet_SI() const { return (rPlanetSi); }
  // return MhdNo2SiL
  inline double getMhdNo2SiL() const { return (MhdNo2SiL); }
  // BATSRUS normalized unit -> PIC normalized unit;
  inline double getMhdNo2NoL() const { return (MhdNo2SiL * Si2NoL); }

  inline double getFluidNxc() const { return nPhyCell_D[ix_]; }
  inline double getFluidNyc() const { return nPhyCell_D[iy_]; }  
  inline double getFluidNzc() const { return nPhyCell_D[iz_]; }

  void save_restart_data() {
    if (isGridEmpty)
      return;

    std::string restartDir = "PC/restartOUT/";
    amrex::VisMF::Write(nodeFluid,
                        restartDir + domainName + "_Interface_nodeFluid");
    amrex::VisMF::Write(centerB,
                        restartDir + domainName + "_Interface_centerB");
  };

  void read_restart() {
    std::string restartDir = "PC/restartIN/";
    amrex::VisMF::Read(nodeFluid,
                       restartDir + domainName + "_Interface_nodeFluid");
    amrex::VisMF::Read(centerB, restartDir + domainName + "_Interface_centerB");
  }

  // ---------Functions to read/interpolate value from nodeFluid.
  // Begin------------
  const amrex::MultiFab& get_nodeFluid() const { return nodeFluid; }

  amrex::Real get_center_b(const amrex::MFIter& mfi, const int i, const int j,
                           const int k, const int iDir) const {
    const auto& arr = centerB[mfi].array();
    return arr(i, j, k, iDir);
  }

  amrex::Real get_value(const amrex::MFIter& mfi, const int i, const int j,
                        const int k, const int iVar) const {
    const auto& arr = nodeFluid[mfi].array();
    return arr(i, j, k, iVar);
  }

  amrex::Real get_value(const amrex::MFIter& mfi, const amrex::Real x,
                        const amrex::Real y, const amrex::Real z,
                        const int iVar) const {
    return get_value_at_loc(nodeFluid, mfi, geom, x, y, z, iVar);
  }

  template <typename T>
  amrex::Real get_number_density(const amrex::MFIter& mfi, const T x, const T y,
                                 const T z, const int is) const {
    amrex::Real Rho, NumDens;

    if (useElectronFluid) {
      Rho = get_value(mfi, x, y, z, iRho_I[is]);
      // TODO: change division to multiplication.
      NumDens = Rho / MoMi0_S[is];
    } else if (useMultiFluid || useMultiSpecies) {
      if (is == 0) {
        // Electron
        NumDens = 0;
        for (int iIon = 0; iIon < nIon; ++iIon) {
          Rho = get_value(mfi, x, y, z, iRho_I[iIon]);
          // TODO: change division to multiplication.
          NumDens += Rho / MoMi0_S[iIon + 1];
        }
      } else {
        // Ion
        Rho = get_value(mfi, x, y, z, iRho_I[is - 1]);
        // TODO: change division to multiplication.
        NumDens = Rho / MoMi0_S[is];
      }
    } else {
      // Electrons and iones have same density, ignoring is
      Rho = get_value(mfi, x, y, z, iRho_I[0]);
      NumDens = Rho * invSumMass;
    }
    return (NumDens);
  }

  template <typename Type>
  amrex::Real get_u(const amrex::MFIter& mfi, const Type x, const Type y,
                    const Type z, const int is, const int* iU_I,
                    const int iJ) const {

    amrex::Real U, J, Rhoit, Qit, Rhot;

    if (useElectronFluid) {
      U = get_value(mfi, x, y, z, iU_I[is]);
    } else if (useMultiFluid) {
      if (is == 0) {
        // Electron
        /** Ue = (J - sum(ni*qi*Ui))/(ne*qe)
               = J/(ne*qe) + sum(ni*Ui)/ne */
        amrex::Real Ui, ni, ne;
        J = get_value(mfi, x, y, z, iJ);
        ne = get_number_density(mfi, x, y, z, 0);
        U = J / (QoQi0_S[0] * ne);

        for (int iIon = 0; iIon < nIon; ++iIon) {
          Ui = get_u(mfi, x, y, z, iIon + 1, iU_I, iJ);
          ni = get_number_density(mfi, x, y, z, iIon + 1);
          U += ni * Ui / ne;
        }
      } else {
        // Ion
        U = get_value(mfi, x, y, z, iU_I[is - 1]);
      }
    } else {
      // Single fluid or multi-species.

      // Ui = U_{MHD} - me/qe*J/Rho;
      // where Rho is total density include electrons.
      // Ue = U_{MHD} - Rhoit/Qit*J/Rho;
      // where Rhoit = sum(ni*Mi), Qit = sum(ni*Qi).

      double moq, Numi;
      if (is == 0) {
        Rhoit = 0;
        Qit = 0;
        for (int iIon = 0; iIon < nIon; ++iIon) {
          Numi = get_number_density(mfi, x, y, z, iIon + 1);
          Rhoit += Numi * MoMi0_S[iIon + 1];
          Qit += Numi * QoQi0_S[iIon + 1];
        }
        moq = 0;
        if (Qit != 0)
          moq = Rhoit / Qit;
      } else
        moq = MoMi0_S[0] / QoQi0_S[0];

      Rhot = 0;
      for (int is0 = 0; is0 < nSIn; ++is0) {
        Rhot += MoMi0_S[is0] * get_number_density(mfi, x, y, z, is0);
      }

      U = get_value(mfi, x, y, z, iU_I[0]);
      J = get_value(mfi, x, y, z, iJ);

      if (Rhot != 0)
        U -= moq * J / Rhot;
    }
    return U;
  }

  template <typename Type>
  amrex::Real get_ux(const amrex::MFIter& mfi, const Type x, const Type y,
                     const Type z, const int is) const {
    return get_u(mfi, x, y, z, is, iUx_I, iJx);
  }

  template <typename Type>
  amrex::Real get_uy(const amrex::MFIter& mfi, const Type x, const Type y,
                     const Type z, const int is) const {
    return get_u(mfi, x, y, z, is, iUy_I, iJy);
  }

  template <typename Type>
  amrex::Real get_uz(const amrex::MFIter& mfi, const Type x, const Type y,
                     const Type z, const int is) const {
    return get_u(mfi, x, y, z, is, iUz_I, iJz);
  }

  template <typename Type>
  amrex::Real get_ppar(const amrex::MFIter& mfi, const Type x, const Type y,
                       const Type z, const int is) const {
    amrex::Real P;
    if (useMultiSpecies || useMultiFluid) {
      std::cout << " getFluidPpar has not implemented for "
                   "multifluid/multispecies!!"
                << std::endl;
      abort();
    }

    if (useElectronFluid) {
      P = get_value(mfi, x, y, z, iPpar_I[is]);
    } else if (useMhdPe) {
      if (is == 0)
        P = get_value(mfi, x, y, z, iPe); // Electron
      if (is == 1)
        P = get_value(mfi, x, y, z, iPpar_I[0]); // Ion
    } else {
      P = get_value(mfi, x, y, z, iPpar_I[0]);
      if (is == 0)
        P *= PeRatio;
      else if (is == 1)
        P *= (1 - PeRatio);
    }

    return P;
  }

  template <typename Type>
  amrex::Real get_p(const amrex::MFIter& mfi, const Type x, const Type y,
                    const Type z, const int is) const {
    amrex::Real P;

    if (useElectronFluid) {
      P = get_value(mfi, x, y, z, iP_I[is]);
    } else if (useMultiFluid) {
      // Multi-fluid.
      if (is == 0)
        P = get_value(mfi, x, y, z, iPe); // Electron
      else
        P = get_value(mfi, x, y, z, iP_I[is - 1]); // Ion
    } else {
      // Single-fluid and multi-species.
      if (!useMhdPe) {
        P = get_value(mfi, x, y, z, iP_I[0]);
        if (is == 0)
          P *= PeRatio;
        else if (is > 0)
          P *= (1 - PeRatio);
      } else {
        if (is == 0)
          P = get_value(mfi, x, y, z, iPe); // Electron
        else if (is > 0)
          P = get_value(mfi, x, y, z, iP_I[0]); // Ion
      }

      // Split pressure among ions.
      if (useMultiSpecies && is > 0) {
        amrex::Real Numit;
        Numit = 0; // Number of all ions.
        for (int iIon = 0; iIon < nIon; ++iIon)
          Numit += get_number_density(mfi, x, y, z, iIon + 1);

        P *= get_number_density(mfi, x, y, z, is) / Numit;
      }
    }
    return P;
  }

  template <typename Type>
  amrex::Real get_pxx(const amrex::MFIter& mfi, const Type x, const Type y,
                      const Type z, const int is) const {
    amrex::Real Pxx;
    if (useAnisoP) {
      amrex::Real Bx = get_value(mfi, x, y, z, iBx);
      amrex::Real By = get_value(mfi, x, y, z, iBy);
      amrex::Real Bz = get_value(mfi, x, y, z, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(mfi, x, y, z, is);
      amrex::Real P = get_p(mfi, x, y, z, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pxx = Pperp + (Ppar - Pperp) * Bx * Bx / Bt2;

    } else {
      Pxx = get_p(mfi, x, y, z, is);
    }
    return (QoQi0_S[is] *
            (Pxx / MoMi0_S[is] + get_number_density(mfi, x, y, z, is) *
                                     pow(get_ux(mfi, x, y, z, is), 2)));
  }

  template <typename Type>
  amrex::Real get_pyy(const amrex::MFIter& mfi, const Type x, const Type y,
                      const Type z, const int is) const {
    amrex::Real Pyy;
    if (useAnisoP) {
      amrex::Real Bx = get_value(mfi, x, y, z, iBx);
      amrex::Real By = get_value(mfi, x, y, z, iBy);
      amrex::Real Bz = get_value(mfi, x, y, z, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(mfi, x, y, z, is);
      amrex::Real P = get_p(mfi, x, y, z, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pyy = Pperp + (Ppar - Pperp) * By * By / Bt2;

    } else {
      Pyy = get_p(mfi, x, y, z, is);
    }
    return (QoQi0_S[is] *
            (Pyy / MoMi0_S[is] + get_number_density(mfi, x, y, z, is) *
                                     pow(get_uy(mfi, x, y, z, is), 2)));
  }

  template <typename Type>
  amrex::Real get_pzz(const amrex::MFIter& mfi, const Type x, const Type y,
                      const Type z, const int is) const {
    amrex::Real Pzz;
    if (useAnisoP) {
      amrex::Real Bx = get_value(mfi, x, y, z, iBx);
      amrex::Real By = get_value(mfi, x, y, z, iBy);
      amrex::Real Bz = get_value(mfi, x, y, z, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(mfi, x, y, z, is);
      amrex::Real P = get_p(mfi, x, y, z, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pzz = Pperp + (Ppar - Pperp) * Bz * Bz / Bt2;

    } else {
      Pzz = get_p(mfi, x, y, z, is);
    }
    return (QoQi0_S[is] *
            (Pzz / MoMi0_S[is] + get_number_density(mfi, x, y, z, is) *
                                     pow(get_uz(mfi, x, y, z, is), 2)));
  }

  template <typename Type>
  amrex::Real get_pxy(const amrex::MFIter& mfi, const Type x, const Type y,
                      const Type z, const int is) const {
    amrex::Real Pxy = 0;
    if (useAnisoP) {
      amrex::Real Bx = get_value(mfi, x, y, z, iBx);
      amrex::Real By = get_value(mfi, x, y, z, iBy);
      amrex::Real Bz = get_value(mfi, x, y, z, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(mfi, x, y, z, is);
      amrex::Real P = get_p(mfi, x, y, z, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pxy = (Ppar - Pperp) * Bx * By / Bt2;
    }

    return QoQi0_S[is] *
           (Pxy / MoMi0_S[is] + get_number_density(mfi, x, y, z, is) *
                                    get_ux(mfi, x, y, z, is) *
                                    get_uy(mfi, x, y, z, is));
  }

  template <typename Type>
  amrex::Real get_pxz(const amrex::MFIter& mfi, const Type x, const Type y,
                      const Type z, const int is) const {
    amrex::Real Pxz = 0;
    if (useAnisoP) {
      amrex::Real Bx = get_value(mfi, x, y, z, iBx);
      amrex::Real By = get_value(mfi, x, y, z, iBy);
      amrex::Real Bz = get_value(mfi, x, y, z, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(mfi, x, y, z, is);
      amrex::Real P = get_p(mfi, x, y, z, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pxz = (Ppar - Pperp) * Bx * Bz / Bt2;
    }

    return QoQi0_S[is] *
           (Pxz / MoMi0_S[is] + get_number_density(mfi, x, y, z, is) *
                                    get_ux(mfi, x, y, z, is) *
                                    get_uz(mfi, x, y, z, is));
  }

  template <typename Type>
  amrex::Real get_pyz(const amrex::MFIter& mfi, const Type x, const Type y,
                      const Type z, const int is) const {
    amrex::Real Pyz = 0;
    if (useAnisoP) {
      amrex::Real Bx = get_value(mfi, x, y, z, iBx);
      amrex::Real By = get_value(mfi, x, y, z, iBy);
      amrex::Real Bz = get_value(mfi, x, y, z, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(mfi, x, y, z, is);
      amrex::Real P = get_p(mfi, x, y, z, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pyz = (Ppar - Pperp) * By * Bz / Bt2;
    }

    return QoQi0_S[is] *
           (Pyz / MoMi0_S[is] + get_number_density(mfi, x, y, z, is) *
                                    get_uy(mfi, x, y, z, is) *
                                    get_uz(mfi, x, y, z, is));
  }

  template <typename Type>
  amrex::Real get_uth_iso(const amrex::MFIter& mfi, const Type x, const Type y,
                          const Type z, const int is) const {
    amrex::Real Uth = 0, p, ni;
    p = get_p(mfi, x, y, z, is);
    ni = get_number_density(mfi, x, y, z, is);
    if (ni > 0)
      Uth = sqrt(p / (ni * MoMi0_S[is]));
    return Uth;
  }

  template <typename Type>
  void set_particle_uth_iso(const amrex::MFIter& mfi, const Type x,
                            const Type y, const Type z, double* u, double* v,
                            double* w, const double rand1, const double rand2,
                            const double rand3, const double rand4,
                            const int is) const {
    double harvest, prob, theta, Uth;

    // u = X velocity
    harvest = rand1;
    prob = sqrt(-2.0 * log(1.0 - .999999999 * harvest));
    harvest = rand2;
    theta = 2.0 * M_PI * harvest;
    Uth = get_uth_iso(mfi, x, y, z, is);

    (*u) = Uth * prob * cos(theta);
    // v = Y velocity
    (*v) = Uth * prob * sin(theta);
    // w = Z velocity
    harvest = rand3;
    prob = sqrt(-2.0 * log(1.0 - .999999999 * harvest));
    harvest = rand4;
    theta = 2.0 * M_PI * harvest;
    (*w) = Uth * prob * cos(theta);
  }

  template <typename Type>
  void set_particle_uth_aniso(const amrex::MFIter& mfi, const Type x,
                              const Type y, const Type z, double* u, double* v,
                              double* w, const double rand1, const double rand2,
                              const double rand3, const double rand4,
                              const int is) const {
    amrex::Real Bx, By, Bz, P, Ppar, Pperp, Uthperp, Uthpar, Uthperp1, Uthperp2,
        prob, theta;
    MDArray<double> norm_DD;
    // indexes for the norm_DD matix
    int Norm_, Perp1_, Perp2_, X_, Y_, Z_;

    if (useMultiFluid || useMultiSpecies) {
      std::cout << " setFluidanisoUth has not implemented for "
                   "multifluid/multispecies!!!"
                << std::endl;
      abort();
    }

    Norm_ = 0;
    Perp1_ = 1;
    Perp2_ = 2;
    X_ = 0;
    Y_ = 1;
    Z_ = 2;

    // Get number density and B at the particle position
    double ni = get_number_density(mfi, x, y, z, is);
    Bx = get_value(mfi, x, y, z, iBx);
    By = get_value(mfi, x, y, z, iBy);
    Bz = get_value(mfi, x, y, z, iBz);

    // Get Parallel and perpendicular presure
    Ppar = get_ppar(mfi, x, y, z, is);
    P = get_p(mfi, x, y, z, is);
    Pperp = 0.5 * (3.0 * P - Ppar);

    // Get 3 vertors spaning the vector space
    norm_DD.init(3, 3);
    MagneticBaseVectors(Bx, By, Bz, norm_DD);

    // Get the thermal verlocities
    prob = sqrt(-2.0 * log(1.0 - .999999999 * rand1));
    theta = 2.0 * M_PI * rand2;
    Uthpar = sqrt(Ppar / (MoMi0_S[is] * ni)) * prob * cos(theta);

    prob = sqrt(-2.0 * log(1.0 - .999999999 * rand3));
    theta = 2.0 * M_PI * rand4;
    Uthperp = sqrt(Pperp / (MoMi0_S[is] * ni)) * prob;
    Uthperp1 = Uthperp * cos(theta);
    Uthperp2 = Uthperp * sin(theta);

    // Set particle thermal velocity
    (*u) = Uthpar * norm_DD(Norm_, X_) + Uthperp1 * norm_DD(Perp1_, X_) +
           Uthperp2 * norm_DD(Perp2_, X_);
    (*v) = Uthpar * norm_DD(Norm_, Y_) + Uthperp1 * norm_DD(Perp1_, Y_) +
           Uthperp2 * norm_DD(Perp2_, Y_);
    (*w) = Uthpar * norm_DD(Norm_, Z_) + Uthperp1 * norm_DD(Perp1_, Z_) +
           Uthperp2 * norm_DD(Perp2_, Z_);
  }

  template <typename Type>
  amrex::Real get_bx(const amrex::MFIter& mfi, const Type x, const Type y,
                     const Type z) const {
    return get_value(mfi, x, y, z, iBx);
  }

  template <typename Type>
  amrex::Real get_by(const amrex::MFIter& mfi, const Type x, const Type y,
                     const Type z) const {
    return get_value(mfi, x, y, z, iBy);
  }

  template <typename Type>
  amrex::Real get_bz(const amrex::MFIter& mfi, const Type x, const Type y,
                     const Type z) const {
    return get_value(mfi, x, y, z, iBz);
  }

  template <typename Type>
  amrex::Real get_ex(const amrex::MFIter& mfi, const Type x, const Type y,
                     const Type z) const {
    amrex::Real Ex;
    if (useElectronFluid) {
      Ex = get_value(mfi, x, y, z, iEx);
    } else {
      Ex = get_uz(mfi, x, y, z, 0) * get_by(mfi, x, y, z) -
           get_uy(mfi, x, y, z, 0) * get_bz(mfi, x, y, z);
    }
    return Ex;
  }

  template <typename Type>
  amrex::Real get_ey(const amrex::MFIter& mfi, const Type x, const Type y,
                     const Type z) const {
    amrex::Real Ey;
    if (useElectronFluid) {
      Ey = get_value(mfi, x, y, z, iEy);
    } else {
      Ey = get_ux(mfi, x, y, z, 0) * get_bz(mfi, x, y, z) -
           get_uz(mfi, x, y, z, 0) * get_bx(mfi, x, y, z);
    }
    return Ey;
  }

  template <typename Type>
  amrex::Real get_ez(const amrex::MFIter& mfi, const Type x, const Type y,
                     const Type z) const {
    amrex::Real Ez;
    if (useElectronFluid) {
      Ez = get_value(mfi, x, y, z, iEz);
    } else {
      Ez = get_uy(mfi, x, y, z, 0) * get_bx(mfi, x, y, z) -
           get_ux(mfi, x, y, z, 0) * get_by(mfi, x, y, z);
    }
    return Ez;
  }
  // ---------Functions to read/interpolate value from nodeFluid.
  // End------------
};
#endif
