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

#include "Constants.h"
#include "Grid.h"
#include "MDArray.h"
#include "ReadParam.h"
#include "Utility.h"
#include "Writer.h"

// amrex::MultiFab does not provide an assignment operator. MultiFabFLEKS is
// derived from MultiFab and provides an empty assignment operator.
class MultiFabFLEKS : public amrex::MultiFab {
public:
  MultiFabFLEKS() {}
  ~MultiFabFLEKS() = default;
  MultiFabFLEKS& operator=(const MultiFabFLEKS& other) { return *this; };
};

class FluidInterface : public Grid {
private:
  static const int OhmUe_ = 1, OhmUi_ = 2, OhmUMHD_ = 3;

  FluidType myType = PICFluid;

  MultiFabFLEKS nodeFluid;
  MultiFabFLEKS centerB;

  bool isGridInitialized = false;
  bool isGridEmpty = false;

  int nCellPerPatch = 1;

  int nDimFluid;

  // Number of variables passing between MHD and PIC.
  int nVarFluid;

  // nVarFluid + (Jx, Jy, Jz)
  int nVarCoupling;

  // Number of fluid at the MHD side. One 'fluid' has its own density,
  // velocity and pressure. Electron can be one fluid.
  int nFluid;

  // Number of ion fluid at the MHD side.
  int nIonFluid = -1;

  // Number of species at the MHD side. One 'species' only has its own density.
  int nSpeciesFluid = 0;

  // Total number of ion/electron species exit in the fluid code.
  int nIon = -1;

  // These default flags are set for stand-alone PIC initialization
  bool useMultiSpecies = false;
  bool useMultiFluid = false;
  bool useElectronFluid = true;
  bool useAnisoP = true;
  bool useMhdPe = false;

  //-------------------------------------------------------------------
  int nS;                       // number of particle species
  amrex::Vector<double> MoMi_S; // masses for the particles species
  amrex::Vector<double> QoQi_S; // charge for each particle species
  //-------------------------------------------------------------------

  // temperature ratio for electrons: PeRatio = Pe/Ptotal
  double PeRatio;

  // Sum of masses of each particle species
  double SumMass, invSumMass;

  amrex::Vector<int> iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I, iPpar_I, iP_I, iUx_I,
      iUy_I, iUz_I;

  int iBx, iBy, iBz, iEx, iEy, iEz, iPe, iJx, iJy, iJz, iRhoTotal;

  double rPlanetSi = 1;

  int ScalingFactor = 1;

  // normalization units for length, velocity, mass and charge
  // Normalized q/m ==1 for proton in CGS units
  double Lnorm, Unorm, lNormSI = 1, uNormSI = 1;
  double Mnorm, Qnorm, mNormSI;

  amrex::Vector<double> Si2No_V, No2Si_V;
  double Si2NoM, Si2NoV, Si2NoRho, Si2NoB, Si2NoP, Si2NoJ, Si2NoL, Si2NoE;
  double No2SiV, No2SiL;

  // Variable names of nodeFluid.
  amrex::Vector<std::string> varNames;

  amrex::Vector<double> uniformState;

  // Length in BATSRUS normalized unit -> Si
  double MhdNo2SiL;

  amrex::Vector<int> vecIdx_I;

  bool useResist = false;
  double etaSI = 0, etaNO = 0;
  int OhmU = OhmUe_;

  bool initFromSWMF;

public:
  FluidInterface(amrex::Geometry const& gm, amrex::AmrInfo const& amrInfo,
                 int nGst, int id, std::string tag,
                 const amrex::Vector<int>& iParam,
                 const amrex::Vector<double>& norm,
                 const amrex::Vector<double>& paramComm);

  FluidInterface(amrex::Geometry const& gm, amrex::AmrInfo const& amrInfo,
                 int nGst, int id, std::string tag, FluidType typeIn = PICFluid)
      : Grid(gm, amrInfo, nGst, id, tag), myType(typeIn) {

    initFromSWMF = false;

    if (myType != PICFluid)
      nS = 0;
  }

  // Initialization from other FluidInterface
  FluidInterface(amrex::Geometry const& gm, amrex::AmrInfo const& amrInfo,
                 int nGst, int id, std::string tag,
                 const FluidInterface* const other,
                 FluidType typeIn = PICFluid);

  ~FluidInterface() = default;

  FluidType my_type() { return myType; };

  void save_amrex_file();

  void read_param(const std::string& command, ReadParam& param);

  void post_process_param();

  void set_var_idx();

  void regrid(const amrex::BoxArray& centerBAIn,
              const amrex::DistributionMapping& dmIn);

  void update_nodeFluid(const MultiFabFLEKS& nodeIn, const double dt);

  int count_couple_node_number();

  int loop_through_node(std::string action, double* const pos_DI = nullptr,
                        const double* const data = nullptr,
                        const int* const index = nullptr);

  void get_couple_node_loc(double* const pos_DI);

  void set_node_fluid(const double* const data, const int* const index);

  void set_node_fluid();

  void calc_current();

  void normalize_fluid_variables();

  void convert_moment_to_velocity();

  void set_plasma_charge_and_mass(amrex::Real qomEl);

  void load_balance(const amrex::DistributionMapping& dmIn);

  void calc_normalized_units();

  /** Get nomal and pendicular vector to magnetic field */
  void calc_mag_base_vector(const double Bx, const double By, const double Bz,
                            MDArray<double>& norm_DD) const;

  void calc_fluid_state(const double* dataPIC_I, double* dataFluid_I) const;

  void print_info() const;

  int get_nCellPerPatch() const { return nCellPerPatch; }

  bool get_UseAnisoP() const { return (useAnisoP); }

  bool get_useElectronFluid() const { return useElectronFluid; }

  int get_fluid_dimension() const { return (nDimFluid); }

  int get_nS() const { return nS; }

  double get_No2Si_V(int idx) const { return (No2Si_V[idx]); }
  double get_Si2No_V(int idx) const { return (Si2No_V[idx]); }
  double get_Si2NoL() const { return (Si2NoL); }
  double get_Si2NoT() const { return Si2NoL / Si2NoV; }
  double get_No2SiL() const { return (No2SiL); }
  double get_No2SiRho() const { return (1. / Si2NoRho); }
  double get_No2SiV() const { return (1. / Si2NoV); }
  double get_Si2NoV() const { return Si2NoV; }
  double get_No2SiB() const { return (1. / Si2NoB); }
  double get_No2SiP() const { return (1. / Si2NoP); }
  double get_No2SiJ() const { return (1. / Si2NoJ); }
  double get_No2SiT() const { return Si2NoV / Si2NoL; }

  double get_species_mass(int i) const { return MoMi_S[i]; };
  double get_species_charge(int i) const { return QoQi_S[i]; };

  double get_lnorm_si() const { return lNormSI; }
  double get_unorm_si() const { return uNormSI; }
  double get_mnorm_si() const { return mNormSI; };

  double get_cLight_SI() const { return uNormSI; }

  double get_rPlanet_SI() const { return rPlanetSi; }

  int get_scaling_factor() const { return ScalingFactor; }

  // return MhdNo2SiL
  double get_MhdNo2SiL() const { return (MhdNo2SiL); }
  // BATSRUS normalized unit -> PIC normalized unit;
  double get_MhdNo2NoL() const { return (MhdNo2SiL * Si2NoL); }

  void set_resistivity(double etaSIIn) {
    // In SI unit R = u_si*L_si/eta_si, where eta_si is magnetic diffusivity
    // with unit m^2/s. In normalized CGS unit R = u_pic*L_pic/(eta_pic/4pi),
    // where eta_pic is also magnetic diffusivity. Magnetic Reynolds number R
    // should not change, so eta_pic = 4*pi*eta_si*Si2NoV*Si2NoL.
    etaSI = etaSIIn;
    useResist = etaSI > 0;
    if (useResist)
      etaNO = fourPI * etaSI * Si2NoV * Si2NoL;
  }

  void set_ohm_u(std::string ss) {
    if (ss.find("ue") != std::string::npos) {
      OhmU = OhmUe_;
    } else if (ss.find("ui") != std::string::npos) {
      OhmU = OhmUi_;
    } else if (ss.find("umhd") != std::string::npos) {
      OhmU = OhmUMHD_;
    } else {
      amrex::Print() << "Error: unknown velocity type: " << ss
                     << " It should be 'ue', 'ui' or 'umhd' " << std::endl;
      amrex::Abort(" ");
    }
  }

  void save_restart_data() {
    if (isGridEmpty)
      return;

    std::string restartDir = component + "/restartOUT/";
    amrex::VisMF::Write(nodeFluid,
                        restartDir + gridName + "_Interface_nodeFluid");
    amrex::VisMF::Write(centerB, restartDir + gridName + "_Interface_centerB");
  };

  void read_restart() {
    std::string restartDir = component + "/restartIN/";
    amrex::VisMF::Read(nodeFluid,
                       restartDir + gridName + "_Interface_nodeFluid");
    amrex::VisMF::Read(centerB, restartDir + gridName + "_Interface_centerB");
  }

  // ---------Functions to read/interpolate value from nodeFluid.
  // Begin------------
  const MultiFabFLEKS& get_nodeFluid() const { return nodeFluid; }

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
    return get_value_at_loc(nodeFluid, mfi, Geom(0), x, y, z, iVar);
  }

  template <typename T>
  amrex::Real get_number_density(const amrex::MFIter& mfi, const T x, const T y,
                                 const T z, const int is) const {
    amrex::Real Rho, NumDens;

    if (useElectronFluid) {
      Rho = get_value(mfi, x, y, z, iRho_I[is]);
      // TODO: change division to multiplication.
      NumDens = Rho / MoMi_S[is];
    } else if (useMultiFluid || useMultiSpecies) {
      if (is == 0) {
        // Electron
        NumDens = 0;
        for (int iIon = 0; iIon < nIon; ++iIon) {
          Rho = get_value(mfi, x, y, z, iRho_I[iIon]);
          // TODO: change division to multiplication.
          NumDens += Rho / MoMi_S[iIon + 1];
        }
      } else {
        // Ion
        Rho = get_value(mfi, x, y, z, iRho_I[is - 1]);
        // TODO: change division to multiplication.
        NumDens = Rho / MoMi_S[is];
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
                    const Type z, const int is, const amrex::Vector<int>& iU_I,
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
        U = J / (QoQi_S[0] * ne);

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
          Rhoit += Numi * MoMi_S[iIon + 1];
          Qit += Numi * QoQi_S[iIon + 1];
        }
        moq = 0;
        if (Qit != 0)
          moq = Rhoit / Qit;
      } else
        moq = MoMi_S[0] / QoQi_S[0];

      Rhot = 0;
      for (int is0 = 0; is0 < nS; ++is0) {
        Rhot += MoMi_S[is0] * get_number_density(mfi, x, y, z, is0);
      }

      U = get_value(mfi, x, y, z, iU_I[0]);
      J = get_value(mfi, x, y, z, iJ);

      if (Rhot != 0)
        U -= moq * J / Rhot;
    }
    return U;
  }

  template <typename Type>
  amrex::Real get_fluid_ux(const amrex::MFIter& mfi, const Type x, const Type y,
                           const Type z, const int is) const {
    return get_value(mfi, x, y, z, iUx_I[is]);
  }

  template <typename Type>
  amrex::Real get_fluid_uy(const amrex::MFIter& mfi, const Type x, const Type y,
                           const Type z, const int is) const {
    return get_value(mfi, x, y, z, iUy_I[is]);
  }

  template <typename Type>
  amrex::Real get_fluid_uz(const amrex::MFIter& mfi, const Type x, const Type y,
                           const Type z, const int is) const {
    return get_value(mfi, x, y, z, iUz_I[is]);
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
    return (QoQi_S[is] *
            (Pxx / MoMi_S[is] + get_number_density(mfi, x, y, z, is) *
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
    return (QoQi_S[is] *
            (Pyy / MoMi_S[is] + get_number_density(mfi, x, y, z, is) *
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
    return (QoQi_S[is] *
            (Pzz / MoMi_S[is] + get_number_density(mfi, x, y, z, is) *
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

    return QoQi_S[is] *
           (Pxy / MoMi_S[is] + get_number_density(mfi, x, y, z, is) *
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

    return QoQi_S[is] *
           (Pxz / MoMi_S[is] + get_number_density(mfi, x, y, z, is) *
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

    return QoQi_S[is] *
           (Pyz / MoMi_S[is] + get_number_density(mfi, x, y, z, is) *
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
      Uth = sqrt(p / (ni * MoMi_S[is]));
    return Uth;
  }

  template <typename Type>
  void set_particle_uth_iso(const amrex::MFIter& mfi, const Type x,
                            const Type y, const Type z, double* u, double* v,
                            double* w, const double rand1, const double rand2,
                            const double rand3, const double rand4,
                            const int is, const double uthIn = -1) const {
    double harvest, prob, theta, Uth;

    // u = X velocity
    harvest = rand1;
    prob = sqrt(-2.0 * log(1.0 - .999999999 * harvest));
    harvest = rand2;
    theta = 2.0 * M_PI * harvest;
    Uth = (uthIn >= 0 ? uthIn : get_uth_iso(mfi, x, y, z, is));

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
                              const int is, const double uthParIn = -1,
                              const double uthPerpIn = -1) const {
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
    calc_mag_base_vector(Bx, By, Bz, norm_DD);

    // Get the thermal verlocities
    prob = sqrt(-2.0 * log(1.0 - .999999999 * rand1));
    theta = 2.0 * M_PI * rand2;
    Uthpar = uthParIn >= 0 ? uthParIn
                           : sqrt(Ppar / (MoMi_S[is] * ni)) * prob * cos(theta);

    prob = sqrt(-2.0 * log(1.0 - .999999999 * rand3));
    theta = 2.0 * M_PI * rand4;
    Uthperp =
        uthPerpIn >= 0 ? uthPerpIn : sqrt(Pperp / (MoMi_S[is] * ni)) * prob;
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
      const bool UseGradPe = false;

      amrex::Real uz, uy;

      if (OhmU == OhmUe_) {
        uz = get_uz(mfi, x, y, z, 0);
        uy = get_uy(mfi, x, y, z, 0);
      } else if (OhmU == OhmUi_) {
        uz = get_uz(mfi, x, y, z, 1);
        uy = get_uy(mfi, x, y, z, 1);
      } else if (OhmU == OhmUMHD_) {
        const amrex::Real r0 = MoMi_S[0] / (MoMi_S[0] + MoMi_S[1]);
        const amrex::Real r1 = 1 - r0;
        uz = r0 * get_uz(mfi, x, y, z, 0) + r1 * get_uz(mfi, x, y, z, 1);
        uy = r0 * get_uy(mfi, x, y, z, 0) + r1 * get_uy(mfi, x, y, z, 1);
      }

      Ex = uz * get_by(mfi, x, y, z) - uy * get_bz(mfi, x, y, z);

      if (UseGradPe) {
        amrex::Real ne = get_number_density(mfi, x, y, z, 0);
        amrex::Real gradpe = get_grad_pe_x(mfi, x, y, z);
        Ex -= gradpe / fabs(QoQi_S[0] * ne);
      }

      if (useResist) {
        Ex += etaNO * get_value(mfi, x, y, z, iJx);
      }
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
      const bool UseGradPe = false;

      amrex::Real ux, uz;

      if (OhmU == OhmUe_) {
        uz = get_uz(mfi, x, y, z, 0);
        ux = get_ux(mfi, x, y, z, 0);
      } else if (OhmU == OhmUi_) {
        uz = get_uz(mfi, x, y, z, 1);
        ux = get_ux(mfi, x, y, z, 1);
      } else if (OhmU == OhmUMHD_) {
        const amrex::Real r0 = MoMi_S[0] / (MoMi_S[0] + MoMi_S[1]);
        const amrex::Real r1 = 1 - r0;
        uz = r0 * get_uz(mfi, x, y, z, 0) + r1 * get_uz(mfi, x, y, z, 1);
        ux = r0 * get_ux(mfi, x, y, z, 0) + r1 * get_ux(mfi, x, y, z, 1);
      }

      Ey = ux * get_bz(mfi, x, y, z) - uz * get_bx(mfi, x, y, z);

      if (UseGradPe) {
        amrex::Real ne = get_number_density(mfi, x, y, z, 0);
        amrex::Real gradpe = get_grad_pe_y(mfi, x, y, z);
        Ey -= gradpe / fabs(QoQi_S[0] * ne);
      }

      if (useResist) {
        Ey += etaNO * get_value(mfi, x, y, z, iJy);
      }
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
      const bool UseGradPe = false;

      amrex::Real ux, uy;

      if (OhmU == OhmUe_) {
        ux = get_ux(mfi, x, y, z, 0);
        uy = get_uy(mfi, x, y, z, 0);
      } else if (OhmU == OhmUi_) {
        ux = get_ux(mfi, x, y, z, 1);
        uy = get_uy(mfi, x, y, z, 1);
      } else if (OhmU == OhmUMHD_) {
        const amrex::Real r0 = MoMi_S[0] / (MoMi_S[0] + MoMi_S[1]);
        const amrex::Real r1 = 1 - r0;
        ux = r0 * get_ux(mfi, x, y, z, 0) + r1 * get_ux(mfi, x, y, z, 1);
        uy = r0 * get_uy(mfi, x, y, z, 0) + r1 * get_uy(mfi, x, y, z, 1);
      }

      Ez = uy * get_bx(mfi, x, y, z) - ux * get_by(mfi, x, y, z);

      if (UseGradPe) {
        amrex::Real ne = get_number_density(mfi, x, y, z, 0);
        amrex::Real gradpe = get_grad_pe_z(mfi, x, y, z);
        Ez -= gradpe / fabs(QoQi_S[0] * ne);
      }

      if (useResist) {
        Ez += etaNO * get_value(mfi, x, y, z, iJz);
      }
    }
    return Ez;
  }

  // Calculate grad(pe) at node (x,y,z). If this node is at the boundary of
  // the fab, it will return grad(pe) at a node that is one cell away from the
  // boundary. Only works when useElectronFluid is False.
  amrex::Real get_grad_pe_x(const amrex::MFIter& mfi, const int x, const int y,
                            const int z) const {

    const amrex::Box& box = mfi.fabbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int xCenter = x;

    if (x == lo.x) {
      xCenter = x + 1;
    } else if (x == hi.x) {
      xCenter = x - 1;
    }

    amrex::Real gradpe =
        0.5 * Geom(0).InvCellSize(ix_) *
        (get_p(mfi, xCenter + 1, y, z, 0) - get_p(mfi, xCenter - 1, y, z, 0));

    return gradpe;
  }

  amrex::Real get_grad_pe_y(const amrex::MFIter& mfi, const int x, const int y,
                            const int z) const {

    const amrex::Box& box = mfi.fabbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int yCenter = y;

    if (y == lo.y) {
      yCenter = y + 1;
    } else if (y == hi.y) {
      yCenter = y - 1;
    }

    amrex::Real gradpe =
        0.5 * Geom(0).InvCellSize(iy_) *
        (get_p(mfi, x, yCenter + 1, z, 0) - get_p(mfi, x, yCenter - 1, z, 0));

    return gradpe;
  }

  amrex::Real get_grad_pe_z(const amrex::MFIter& mfi, const int x, const int y,
                            const int z) const {

    const amrex::Box& box = mfi.fabbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    int zCenter = z;

    if (z == lo.z) {
      zCenter = z + 1;
    } else if (z == hi.z) {
      zCenter = z - 1;
    }

    amrex::Real gradpe =
        0.5 * Geom(0).InvCellSize(iz_) *
        (get_p(mfi, x, y, zCenter + 1, 0) - get_p(mfi, x, y, zCenter - 1, 0));

    return gradpe;
  }

  // ---------Functions to read/interpolate value from nodeFluid.
  // End------------
};
#endif
