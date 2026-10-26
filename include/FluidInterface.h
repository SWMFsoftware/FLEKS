#ifndef _FLUIDINTERFACE_H_
#define _FLUIDINTERFACE_H_

#include <array>
#include <memory>

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
#include "GpuBackend.h"
#include "Grid.h"
#include "GridAccess.h"
#include "GridUtility.h"
#include "ReadParam.h"

// Immutable SI<->normalized conversion factors, computed once from a
// FluidInterface and shared read-only (shared_ptr<const>) by all interfaces.
class FluidInterface; // forward declaration; defined later in this header

using Basis3x3 = std::array<std::array<amrex::Real, 3>, 3>;

class NormalizationParams {
public:
  NormalizationParams() = default;
  // Build the derived normalization from fi.
  NormalizationParams(const FluidInterface& fi, bool scalarOnly = false);

  // Derived temperature conversions.
  double si2noT() const { return Si2NoL / Si2NoV; }
  double no2siT() const { return Si2NoV / Si2NoL; }

  // normalization units for length, velocity, mass and charge
  // Normalized q/m == 1 for proton in CGS units
  double Lnorm = 1, Unorm = 1;
  double Mnorm = 0, Qnorm = 0;

  amrex::Vector<double> Si2No_V, No2Si_V;
  double Si2NoM = 0, Si2NoV = 0, Si2NoRho = 0, Si2NoB = 0, Si2NoP = 0,
         Si2NoJ = 0, Si2NoL = 0, Si2NoE = 0;
  double No2SiV = 0, No2SiL = 0;

private:
  void calc_normalization_units(double lNormSI, double uNormSI, double mNormSI);
  void compute_var_conversions(const FluidInterface& fi);
};

// The node-fluid data of one FAB and the geometry that locates a point in it.
// Device-copyable, so a kernel can evaluate the fluid state (see
// FluidStateParams). In a GPU build `fluid` points to device memory.
struct FluidData {
  amrex::Array4<amrex::Real const> fluid;
  amrex::Real plo[3] = { 0, 0, 0 };
  amrex::Real invDx[3] = { 0, 0, 0 };

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real value(
      const amrex::IntVect& ijk, const int iVar) const {
    return fluid(ijk, iVar);
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real value(
      const amrex::RealVect& xyz, const int iVar) const {
    return value_at_loc(fluid, plo, invDx, xyz, iVar);
  }
};

// The parameters of the fluid state and the formulas that derive the plasma
// state (densities, velocities, pressures, E, ...) from the node fluid. It is
// a plain struct of POD members, so a kernel can capture a copy of it; the
// formulas are written once, here, and are used both by the host accessors of
// FluidInterface and by the kernels (see FluidInterface::fluid_view()).
class FluidStateParams {
protected:
  static const int OhmUe_ = 1, OhmUi_ = 2, OhmUMHD_ = 3;

  // Total number of ion/electron species exist in the fluid code.
  int nIon = -1;

  // These default flags are set for stand-alone PIC initialization
  bool useMultiSpecies = false;
  bool useMultiFluid = false;
  bool useElectronFluid = true;
  bool useAnisoP = false;
  bool useMhdPe = false;
  bool useHybridPIC = false;

  KernelArray<int> iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I, iPpar_I, iP_I, iUx_I,
      iUy_I, iUz_I;

  int iBx = -1, iBy = -1, iBz = -1;
  int iEx = -1, iEy = -1, iEz = -1;
  int iPe = -1;

  //-------------------------------------------------------------------
  int nS = 0;                 // number of particle species
  KernelArray<double> MoMi_S; // masses for the particles species
  KernelArray<double> QoQi_S; // charge for each particle species
  //-------------------------------------------------------------------

  // temperature ratio for electrons: PeRatio = Pe/Ptotal
  double PeRatio = 0;

  // Sum of masses of each particle species
  double SumMass = 0, invSumMass = 0;

  int iJx, iJy, iJz;

  bool useResist = false;
  double etaSI = 0, etaNO = 0;
  int OhmU = OhmUe_;

  bool initFromSWMF;

public:
  /** Get normal and perpendicular unit vectors relative to magnetic field */
  AMREX_GPU_HOST_DEVICE static void calc_mag_base_vector(const amrex::Real Bx,
                                                         const amrex::Real By,
                                                         const amrex::Real Bz,
                                                         Basis3x3& norm_DD) {
    enum { Norm_ = 0, Perp1_ = 1, Perp2_ = 2, X_ = 0, Y_ = 1, Z_ = 2 };

    amrex::Real inv = 1.0 / std::sqrt(Bx * Bx + By * By + Bz * Bz);
    norm_DD[Norm_][X_] = Bx * inv;
    norm_DD[Norm_][Y_] = By * inv;
    norm_DD[Norm_][Z_] = Bz * inv;

    if (norm_DD[Norm_][Z_] < 0.5) {
      norm_DD[Perp1_][X_] = norm_DD[Norm_][Y_];
      norm_DD[Perp1_][Y_] = -norm_DD[Norm_][X_];
      norm_DD[Perp1_][Z_] = 0.0;
    } else {
      norm_DD[Perp1_][X_] = 0.0;
      norm_DD[Perp1_][Y_] = norm_DD[Norm_][Z_];
      norm_DD[Perp1_][Z_] = -norm_DD[Norm_][Y_];
    }

    inv = 1.0 / std::sqrt(norm_DD[Perp1_][X_] * norm_DD[Perp1_][X_] +
                          norm_DD[Perp1_][Y_] * norm_DD[Perp1_][Y_] +
                          norm_DD[Perp1_][Z_] * norm_DD[Perp1_][Z_]);
    norm_DD[Perp1_][X_] *= inv;
    norm_DD[Perp1_][Y_] *= inv;
    norm_DD[Perp1_][Z_] *= inv;

    // Perp2 = Norm x Perp1 (cross product of orthonormal unit vectors is
    // already unit length)
    norm_DD[Perp2_][X_] = norm_DD[Norm_][Y_] * norm_DD[Perp1_][Z_] -
                          norm_DD[Norm_][Z_] * norm_DD[Perp1_][Y_];
    norm_DD[Perp2_][Y_] = norm_DD[Norm_][Z_] * norm_DD[Perp1_][X_] -
                          norm_DD[Norm_][X_] * norm_DD[Perp1_][Z_];
    norm_DD[Perp2_][Z_] = norm_DD[Norm_][X_] * norm_DD[Perp1_][Y_] -
                          norm_DD[Norm_][Y_] * norm_DD[Perp1_][X_];
  }

  // The fluid state at the node ijk (IntVect) or the location xyz (RealVect)
  // of the FAB `d`.
  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_pe(const FluidData& d,
                                           const Type xyz) const {
    if (iPe >= 0) {
      return d.value(xyz, iPe);
    } else if (!iP_I.empty() && iP_I[0] >= 0 && PeRatio > 0.0) {
      return PeRatio * d.value(xyz, iP_I[0]);
    }
    return 0.0;
  }

  template <typename T>
  AMREX_GPU_HOST_DEVICE amrex::Real get_number_density(const FluidData& d,
                                                       const T xyz,
                                                       const int is) const {
    amrex::Real Rho, NumDens;

    if (!initFromSWMF && is < static_cast<int>(iRho_I.size())) {
      Rho = d.value(xyz, iRho_I[is]);
      return Rho / MoMi_S[is];
    }

    if (useHybridPIC) {
      if (useMultiFluid || useMultiSpecies) {
        if (is < static_cast<int>(iRho_I.size()) && iRho_I[is] >= 0) {
          Rho = d.value(xyz, iRho_I[is]);
          return (is < static_cast<int>(MoMi_S.size()) && MoMi_S[is] > 0.0)
                     ? (Rho / MoMi_S[is])
                     : Rho;
        }
      } else {
        if (!iRho_I.empty() && iRho_I[0] >= 0) {
          Rho = d.value(xyz, iRho_I[0]);
          return (is < static_cast<int>(MoMi_S.size()) && MoMi_S[is] > 0.0)
                     ? (Rho / MoMi_S[is])
                     : Rho;
        }
      }
      return 0.0;
    }

    if (useElectronFluid) {
      Rho = d.value(xyz, iRho_I[is]);
      NumDens = Rho / MoMi_S[is];
    } else if (useMultiFluid || useMultiSpecies) {
      if (is == 0) {
        // Electron
        NumDens = 0;
        for (int iIon = 0; iIon < nIon; ++iIon) {
          Rho = d.value(xyz, iRho_I[iIon]);
          NumDens += Rho / MoMi_S[iIon + 1];
        }
      } else {
        // Ion
        Rho = d.value(xyz, iRho_I[is - 1]);
        NumDens = Rho / MoMi_S[is];
      }
    } else {
      // Electrons and iones have same density, ignoring is
      Rho = d.value(xyz, iRho_I[0]);
      NumDens = Rho * invSumMass;
    }
    return (NumDens);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_u(const FluidData& d, const Type xyz,
                                          const int is,
                                          const KernelArray<int>& iU_I,
                                          const int iJ) const {

    amrex::Real U, J, Rhoit, Qit, Rhot;

    if (useHybridPIC) {
      if (useMultiFluid) {
        if (is < static_cast<int>(iU_I.size()) && iU_I[is] >= 0) {
          return d.value(xyz, iU_I[is]);
        }
      } else {
        if (!iU_I.empty() && iU_I[0] >= 0) {
          return d.value(xyz, iU_I[0]);
        }
      }
      return 0.0;
    }

    if (useElectronFluid) {
      U = d.value(xyz, iU_I[is]);
    } else if (useMultiFluid) {
      if (is == 0) {
        // Electron
        /** Ue = (J - sum(ni*qi*Ui))/(ne*qe)
               = J/(ne*qe) + sum(ni*Ui)/ne */
        amrex::Real Ui, ni, ne;
        J = d.value(xyz, iJ);
        ne = get_number_density(d, xyz, 0);
        U = J / (QoQi_S[0] * ne);

        for (int iIon = 0; iIon < nIon; ++iIon) {
          Ui = get_u(d, xyz, iIon + 1, iU_I, iJ);
          ni = get_number_density(d, xyz, iIon + 1);
          U += ni * Ui / ne;
        }
      } else {
        // Ion
        U = d.value(xyz, iU_I[is - 1]);
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
          Numi = get_number_density(d, xyz, iIon + 1);
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
        Rhot += MoMi_S[is0] * get_number_density(d, xyz, is0);
      }

      U = d.value(xyz, iU_I[0]);
      J = d.value(xyz, iJ);

      if (Rhot != 0)
        U -= moq * J / Rhot;
    }
    return U;
  }

  template <typename T>
  AMREX_GPU_HOST_DEVICE amrex::Real get_fluid_mass_density(const FluidData& d,
                                                           const T xyz,
                                                           const int is) const {
    return d.value(xyz, iRho_I[is]);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_fluid_p(const FluidData& d,
                                                const Type xyz,
                                                const int is) const {
    return d.value(xyz, iP_I[is]);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_fluid_uth(const FluidData& d,
                                                  const Type xyz,
                                                  const int is) const {
    // 'uth' returned by this method is defined as: uth = sqrt(kT/m)
    // PV = nkT
    // kT/m = PV/(nm) = P/rho
    // So, uth = sqrt(P/rho)
    amrex::Real uth = 0, p, rho;
    p = get_fluid_p(d, xyz, is);
    rho = get_fluid_mass_density(d, xyz, is);
    if (rho > 0)
      uth = sqrt(p / rho);
    return uth;
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_fluid_ux(const FluidData& d,
                                                 const Type xyz,
                                                 const int is) const {
    return d.value(xyz, iUx_I[is]);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_fluid_uy(const FluidData& d,
                                                 const Type xyz,
                                                 const int is) const {
    return d.value(xyz, iUy_I[is]);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_fluid_uz(const FluidData& d,
                                                 const Type xyz,
                                                 const int is) const {
    return d.value(xyz, iUz_I[is]);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_ux(const FluidData& d, const Type xyz,
                                           const int is) const {
    return get_u(d, xyz, is, iUx_I, iJx);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_uy(const FluidData& d, const Type xyz,
                                           const int is) const {
    return get_u(d, xyz, is, iUy_I, iJy);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_uz(const FluidData& d, const Type xyz,
                                           const int is) const {
    return get_u(d, xyz, is, iUz_I, iJz);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_ppar(const FluidData& d, const Type xyz,
                                             const int is) const {
    amrex::Real P;
    if (useMultiSpecies || useMultiFluid) {
      amrex::Abort("get_ppar has not been implemented for "
                   "multifluid/multispecies!");
    }

    if (useHybridPIC) {
      if (useMhdPe) {
        return d.value(xyz, iPpar_I[0]);
      } else {
        return d.value(xyz, iPpar_I[0]) * (1.0 - PeRatio);
      }
    }

    if (useElectronFluid) {
      P = d.value(xyz, iPpar_I[is]);
    } else if (useMhdPe) {
      if (is == 0)
        P = d.value(xyz, iPe); // Electron
      if (is == 1)
        P = d.value(xyz, iPpar_I[0]); // Ion
    } else {
      P = d.value(xyz, iPpar_I[0]);
      if (is == 0)
        P *= PeRatio;
      else if (is == 1)
        P *= (1 - PeRatio);
    }

    return P;
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_p(const FluidData& d, const Type xyz,
                                          const int is) const {
    if (!initFromSWMF && is < static_cast<int>(iP_I.size())) {
      return d.value(xyz, iP_I[is]);
    }

    if (useHybridPIC) {
      if (useMultiFluid) {
        if (is < static_cast<int>(iP_I.size()) && iP_I[is] >= 0) {
          return d.value(xyz, iP_I[is]);
        }
        return 0.0;
      } else {
        amrex::Real pTotIon = 0.0;
        if (!iP_I.empty() && iP_I[0] >= 0) {
          pTotIon = d.value(xyz, iP_I[0]);
          if (!useMhdPe) {
            pTotIon *= (1.0 - PeRatio);
          }
        }
        if (useMultiSpecies && nIon > 1) {
          amrex::Real numTot = 0.0;
          for (int iIon = 0; iIon < nIon; ++iIon) {
            numTot += get_number_density(d, xyz, iIon);
          }
          if (numTot > 0.0) {
            return pTotIon * (get_number_density(d, xyz, is) / numTot);
          }
        }
        return pTotIon;
      }
    }

    amrex::Real P = 0.0;

    if (useElectronFluid) {
      P = d.value(xyz, iP_I[is]);
    } else if (useMultiFluid) {
      // Multi-fluid.
      if (is == 0)
        P = d.value(xyz, iPe); // Electron
      else
        P = d.value(xyz, iP_I[is - 1]); // Ion
    } else {
      // Single-fluid and multi-species.
      if (!useMhdPe) {
        P = d.value(xyz, iP_I[0]);
        if (is == 0)
          P *= PeRatio;
        else if (is > 0)
          P *= (1 - PeRatio);
      } else {
        if (is == 0)
          P = d.value(xyz, iPe); // Electron
        else if (is > 0)
          P = d.value(xyz, iP_I[0]); // Ion
      }

      // Split pressure among ions.
      if (useMultiSpecies && is > 0) {
        amrex::Real Numit;
        Numit = 0; // Number of all ions.
        for (int iIon = 0; iIon < nIon; ++iIon)
          Numit += get_number_density(d, xyz, iIon + 1);

        P *= get_number_density(d, xyz, is) / Numit;
      }
    }
    return P;
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_pxx(const FluidData& d, const Type xyz,
                                            const int is) const {
    amrex::Real Pxx;
    if (useAnisoP) {
      amrex::Real Bx = d.value(xyz, iBx);
      amrex::Real By = d.value(xyz, iBy);
      amrex::Real Bz = d.value(xyz, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(d, xyz, is);
      amrex::Real P = get_p(d, xyz, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pxx = Pperp + (Ppar - Pperp) * Bx * Bx / Bt2;

    } else {
      Pxx = get_p(d, xyz, is);
    }
    const amrex::Real ux = get_ux(d, xyz, is);
    return (QoQi_S[is] *
            (Pxx / MoMi_S[is] + get_number_density(d, xyz, is) * ux * ux));
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_pyy(const FluidData& d, const Type xyz,
                                            const int is) const {
    amrex::Real Pyy;
    if (useAnisoP) {
      amrex::Real Bx = d.value(xyz, iBx);
      amrex::Real By = d.value(xyz, iBy);
      amrex::Real Bz = d.value(xyz, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(d, xyz, is);
      amrex::Real P = get_p(d, xyz, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pyy = Pperp + (Ppar - Pperp) * By * By / Bt2;

    } else {
      Pyy = get_p(d, xyz, is);
    }
    const amrex::Real uy = get_uy(d, xyz, is);
    return (QoQi_S[is] *
            (Pyy / MoMi_S[is] + get_number_density(d, xyz, is) * uy * uy));
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_pzz(const FluidData& d, const Type xyz,
                                            const int is) const {
    amrex::Real Pzz;
    if (useAnisoP) {
      amrex::Real Bx = d.value(xyz, iBx);
      amrex::Real By = d.value(xyz, iBy);
      amrex::Real Bz = d.value(xyz, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(d, xyz, is);
      amrex::Real P = get_p(d, xyz, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pzz = Pperp + (Ppar - Pperp) * Bz * Bz / Bt2;

    } else {
      Pzz = get_p(d, xyz, is);
    }
    const amrex::Real uz = get_uz(d, xyz, is);
    return (QoQi_S[is] *
            (Pzz / MoMi_S[is] + get_number_density(d, xyz, is) * uz * uz));
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_pxy(const FluidData& d, const Type xyz,
                                            const int is) const {
    amrex::Real Pxy = 0;
    if (useAnisoP) {
      amrex::Real Bx = d.value(xyz, iBx);
      amrex::Real By = d.value(xyz, iBy);
      amrex::Real Bz = d.value(xyz, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(d, xyz, is);
      amrex::Real P = get_p(d, xyz, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pxy = (Ppar - Pperp) * Bx * By / Bt2;
    }

    return QoQi_S[is] *
           (Pxy / MoMi_S[is] + get_number_density(d, xyz, is) *
                                   get_ux(d, xyz, is) * get_uy(d, xyz, is));
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_pxz(const FluidData& d, const Type xyz,
                                            const int is) const {
    amrex::Real Pxz = 0;
    if (useAnisoP) {
      amrex::Real Bx = d.value(xyz, iBx);
      amrex::Real By = d.value(xyz, iBy);
      amrex::Real Bz = d.value(xyz, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(d, xyz, is);
      amrex::Real P = get_p(d, xyz, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pxz = (Ppar - Pperp) * Bx * Bz / Bt2;
    }

    return QoQi_S[is] *
           (Pxz / MoMi_S[is] + get_number_density(d, xyz, is) *
                                   get_ux(d, xyz, is) * get_uz(d, xyz, is));
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_pyz(const FluidData& d, const Type xyz,
                                            const int is) const {
    amrex::Real Pyz = 0;
    if (useAnisoP) {
      amrex::Real Bx = d.value(xyz, iBx);
      amrex::Real By = d.value(xyz, iBy);
      amrex::Real Bz = d.value(xyz, iBz);
      amrex::Real Bt2 = Bx * Bx + By * By + Bz * Bz;

      amrex::Real Ppar = get_ppar(d, xyz, is);
      amrex::Real P = get_p(d, xyz, is);
      amrex::Real Pperp = 0.5 * (3.0 * P - Ppar);

      Pyz = (Ppar - Pperp) * By * Bz / Bt2;
    }

    return QoQi_S[is] *
           (Pyz / MoMi_S[is] + get_number_density(d, xyz, is) *
                                   get_uy(d, xyz, is) * get_uz(d, xyz, is));
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_uth_iso(const FluidData& d,
                                                const Type xyz,
                                                const int is) const {
    amrex::Real uth = 0, p, ni;
    p = get_p(d, xyz, is);
    ni = get_number_density(d, xyz, is);
    if (ni > 0)
      uth = sqrt(p / (ni * MoMi_S[is]));
    return uth;
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE void set_particle_uth_iso(
      const FluidData& d, const Type xyz, double* u, double* v, double* w,
      const double rand1, const double rand2, const double rand3,
      const double rand4, const int is, const double uthIn = -1) const {
    double harvest, prob, theta, uth;

    // u = X velocity
    harvest = rand1;
    prob = sqrt(-2.0 * log(1.0 - .999999999 * harvest));
    harvest = rand2;
    theta = 2.0 * M_PI * harvest;
    uth = (uthIn >= 0 ? uthIn : get_uth_iso(d, xyz, is));

    (*u) = uth * prob * cos(theta);
    // v = Y velocity
    (*v) = uth * prob * sin(theta);
    // w = Z velocity
    harvest = rand3;
    prob = sqrt(-2.0 * log(1.0 - .999999999 * harvest));
    harvest = rand4;
    theta = 2.0 * M_PI * harvest;
    (*w) = uth * prob * cos(theta);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE void set_particle_uth_aniso(
      const FluidData& d, const Type xyz, double* u, double* v, double* w,
      const double rand1, const double rand2, const double rand3,
      const double rand4, const int is, const double uthParIn = -1,
      const double uthPerpIn = -1) const {
    amrex::Real Bx, By, Bz, P, Ppar, Pperp, Uthperp, Uthpar, Uthperp1, Uthperp2,
        prob, theta;
    Basis3x3 norm_DD;
    // indices for the norm_DD matrix
    enum { Norm_ = 0, Perp1_ = 1, Perp2_ = 2, X_ = 0, Y_ = 1, Z_ = 2 };

    if (useMultiFluid || useMultiSpecies) {
      amrex::Abort("setFluidanisoUth has not been implemented for "
                   "multifluid/multispecies!");
    }

    // Get number density and B at the particle position
    double ni = get_number_density(d, xyz, is);
    Bx = d.value(xyz, iBx);
    By = d.value(xyz, iBy);
    Bz = d.value(xyz, iBz);

    // Get Parallel and perpendicular pressure
    Ppar = get_ppar(d, xyz, is);
    P = get_p(d, xyz, is);
    Pperp = 0.5 * (3.0 * P - Ppar);

    // Get 3 unit vectors spanning the vector space (stack-allocated, zero heap)
    calc_mag_base_vector(Bx, By, Bz, norm_DD);

    // Get the thermal velocities
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
    (*u) = Uthpar * norm_DD[Norm_][X_] + Uthperp1 * norm_DD[Perp1_][X_] +
           Uthperp2 * norm_DD[Perp2_][X_];
    (*v) = Uthpar * norm_DD[Norm_][Y_] + Uthperp1 * norm_DD[Perp1_][Y_] +
           Uthperp2 * norm_DD[Perp2_][Y_];
    (*w) = Uthpar * norm_DD[Norm_][Z_] + Uthperp1 * norm_DD[Perp1_][Z_] +
           Uthperp2 * norm_DD[Perp2_][Z_];
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_bx(const FluidData& d,
                                           const Type xyz) const {
    return d.value(xyz, iBx);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_by(const FluidData& d,
                                           const Type xyz) const {
    return d.value(xyz, iBy);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_bz(const FluidData& d,
                                           const Type xyz) const {
    return d.value(xyz, iBz);
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_ex(const FluidData& d,
                                           const Type xyz) const {
    amrex::Real Ex;
    if (useElectronFluid) {
      Ex = d.value(xyz, iEx);
    } else {

      amrex::Real uz = 0, uy = 0;

      if (useHybridPIC) {
        uz = get_uz(d, xyz, 0);
        uy = get_uy(d, xyz, 0);
      } else if (OhmU == OhmUe_) {
        uz = get_uz(d, xyz, 0);
        uy = get_uy(d, xyz, 0);
      } else if (OhmU == OhmUi_) {
        uz = get_uz(d, xyz, 1);
        uy = get_uy(d, xyz, 1);
      } else if (OhmU == OhmUMHD_) {
        const amrex::Real r0 = MoMi_S[0] / (MoMi_S[0] + MoMi_S[1]);
        const amrex::Real r1 = 1 - r0;
        uz = r0 * get_uz(d, xyz, 0) + r1 * get_uz(d, xyz, 1);
        uy = r0 * get_uy(d, xyz, 0) + r1 * get_uy(d, xyz, 1);
      }

      Ex = uz * get_by(d, xyz) - uy * get_bz(d, xyz);

      if (useResist) {
        Ex += etaNO * d.value(xyz, iJx);
      }
    }
    return Ex;
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_ey(const FluidData& d,
                                           const Type xyz) const {
    amrex::Real Ey;
    if (useElectronFluid) {
      Ey = d.value(xyz, iEy);
    } else {

      amrex::Real ux = 0, uz = 0;

      if (useHybridPIC) {
        uz = get_uz(d, xyz, 0);
        ux = get_ux(d, xyz, 0);
      } else if (OhmU == OhmUe_) {
        uz = get_uz(d, xyz, 0);
        ux = get_ux(d, xyz, 0);
      } else if (OhmU == OhmUi_) {
        uz = get_uz(d, xyz, 1);
        ux = get_ux(d, xyz, 1);
      } else if (OhmU == OhmUMHD_) {
        const amrex::Real r0 = MoMi_S[0] / (MoMi_S[0] + MoMi_S[1]);
        const amrex::Real r1 = 1 - r0;
        uz = r0 * get_uz(d, xyz, 0) + r1 * get_uz(d, xyz, 1);
        ux = r0 * get_ux(d, xyz, 0) + r1 * get_ux(d, xyz, 1);
      }

      Ey = ux * get_bz(d, xyz) - uz * get_bx(d, xyz);

      if (useResist) {
        Ey += etaNO * d.value(xyz, iJy);
      }
    }

    return Ey;
  }

  template <typename Type>
  AMREX_GPU_HOST_DEVICE amrex::Real get_ez(const FluidData& d,
                                           const Type xyz) const {
    amrex::Real Ez;
    if (useElectronFluid) {
      Ez = d.value(xyz, iEz);
    } else {

      amrex::Real ux = 0, uy = 0;

      if (useHybridPIC) {
        ux = get_ux(d, xyz, 0);
        uy = get_uy(d, xyz, 0);
      } else if (OhmU == OhmUe_) {
        ux = get_ux(d, xyz, 0);
        uy = get_uy(d, xyz, 0);
      } else if (OhmU == OhmUi_) {
        ux = get_ux(d, xyz, 1);
        uy = get_uy(d, xyz, 1);
      } else if (OhmU == OhmUMHD_) {
        const amrex::Real r0 = MoMi_S[0] / (MoMi_S[0] + MoMi_S[1]);
        const amrex::Real r1 = 1 - r0;
        ux = r0 * get_ux(d, xyz, 0) + r1 * get_ux(d, xyz, 1);
        uy = r0 * get_uy(d, xyz, 0) + r1 * get_uy(d, xyz, 1);
      }

      Ez = uy * get_bx(d, xyz) - ux * get_by(d, xyz);

      if (useResist) {
        Ez += etaNO * d.value(xyz, iJz);
      }
    }
    return Ez;
  }
};

// A copy of the fluid parameters together with the node fluid of one FAB:
// everything a kernel needs to evaluate the fluid state.
struct FluidView : public FluidStateParams {
  FluidData d;
};

class FluidInterfaceParameters : public FluidStateParams {
protected:
  FluidInterfaceParameters() = default;
  FluidInterfaceParameters(const FluidInterfaceParameters& fip) = default;

public:
  FluidInterfaceParameters& operator=(const FluidInterfaceParameters& fip) =
      default;

protected:
  double tStartSI;

  int nCellPerPatch = 1;

  // If true, nodeFluid contains (Jx, Jy, Jz)
  bool useCurrent = true;

  int nDimFluid;

  // Number of variables passing between MHD and PIC.
  int nVarFluid = 0;

  // Number of fluid at the MHD side. One 'fluid' has its own density,
  // velocity and pressure. Electron can be one fluid.
  int nFluid = 1;

  // Number of species at the MHD side. One 'species' only has its own density.
  int nSpeciesFluid = 0;

  int iRhoTotal = -1, iLevSet = -1;

  // Variable names of nodeFluid.
  amrex::Vector<std::string> varNames;

  // Normalization base scalars, read from SWMF iParam/norm or from
  // #NORMALIZATION/#SCALINGFACTOR/#PLANETRADIUS.
  double lNormSI = 1.0, uNormSI = 1.0, mNormSI = 1.0;
  double rPlanetSi = 1.0, ScalingFactor = 1.0, MhdNo2SiL = 1.0;

  // Read-only once published; changable via finalize_normalization().
  std::shared_ptr<const NormalizationParams> normParams =
      std::make_shared<const NormalizationParams>();

  amrex::Vector<double> uniformState;

public:
  // Inflow boundary upstream state per species in code units.
  struct InflowVel {
    double nDens = -1.0;
    double vth = 0.0; // sqrt(kT/m)
    double ux = 0.0, uy = 0.0, uz = 0.0;
  };
  amrex::Vector<InflowVel> inflowState;
  bool inflowDefined = false;

protected:
  std::string restartOut = component + "/restartOUT/";

  // Exosphere neutral profile parameters
  std::string exosphereType = "None";
  int nExoComponent = 0;
  amrex::Vector<double> exoN0;
  amrex::Vector<double> exoH0;
  amrex::Vector<double> exoT0;
  amrex::Vector<double> exoK0;

public:
  amrex::Vector<std::string> get_var_names() const { return varNames; }
  void set_var_names(const amrex::Vector<std::string>& names) {
    varNames = names;
  }

  void set_restart_out_dir(std::string dir) {
    restartOut = dir.empty() ? component + "/restartOUT/" : dir;
  }
  std::string get_restart_out_dir() const { return restartOut; }

  bool get_inflow_defined() const { return inflowDefined; }
  const InflowVel* get_inflow_vel(const int iS) const {
    if (!inflowDefined || iS < 0 || iS >= inflowState.size())
      return nullptr;
    return &inflowState[iS];
  }
  void set_inflow_state(const amrex::Vector<InflowVel>& s) { inflowState = s; }
  void set_inflow_defined(bool v) { inflowDefined = v; }
  const amrex::Vector<double>& get_uniform_state() const {
    return uniformState;
  }
};

class FluidInterface : public GridAccess, public FluidInterfaceParameters {
  // NormalizationParams is a derived, immutable snapshot of fi's parameters.
  friend class NormalizationParams;

protected:
  // The shared Grid reference and the mesh queries live in GridAccess.
  FluidType myType = PICFluid;

  std::string tag;
  std::string gridName;
  std::string printPrefix;
  int gridID;

  amrex::Vector<amrex::MultiFab> nodeFluid;
  amrex::Vector<amrex::MultiFab> centerB;

  bool isnodeFluidReady = false;

public:
  // Shared mesh queries (n_lev, Geom(iLev), DistributionMap, cell_status,
  // get_base_grid, lev_string, get_finest_lev, ...) are inherited from
  // GridAccess.

  FluidInterface(Grid& gridIn, std::string tag,
                 const amrex::Vector<int>& iParam,
                 const amrex::Vector<double>& norm,
                 const amrex::Vector<double>& paramComm,
                 bool useHybridIn = false);

  FluidInterface(Grid& gridIn, std::string tag, bool useHybridIn = false,
                 FluidType typeIn = PICFluid)
      : GridAccess(gridIn), myType(typeIn), tag(tag), gridID(gridIn.get_id()) {
    gridName = std::string("FLEKS") + std::to_string(gridID);
    printPrefix = tag.empty() ? gridName + ": " : gridName + " " + tag + ": ";
    initFromSWMF = false;
    useHybridPIC = useHybridIn;

    if (myType != PICFluid)
      nS = 0;
  }

  // Initialization from other FluidInterface
  FluidInterface(Grid& gridIn, const FluidInterface& other, int id,
                 std::string tag, FluidType typeIn = PICFluid)
      : GridAccess(gridIn),
        FluidInterfaceParameters(other),
        myType(typeIn),
        tag(tag),
        gridID(id) {
    gridName = std::string("FLEKS") + std::to_string(gridID);
    printPrefix = tag.empty() ? gridName + ": " : gridName + " " + tag + ": ";
  }

  virtual ~FluidInterface() = default;

  FluidType my_type() { return myType; }

  // Copy the FluidInterfaceParameters slice from another FluidInterface.
  // Used by Domain to keep the source / test-particle grids' parameters in
  // sync with the primary fluid interface: they are copy-constructed from
  // *fi before read_param populates *fi, so their snapshot is stale until
  // this sync runs.  Only the FluidInterfaceParameters members are copied;
  // the Grid data and SourceInterface / UserSource members are preserved.
  void sync_fluid_interface_params(const FluidInterface& other) {
    static_cast<FluidInterfaceParameters&>(*this) =
        static_cast<const FluidInterfaceParameters&>(other);
  }

  void set_period_start_si(double t) { tStartSI = t; }

  double get_period_start_si() const { return tStartSI; }

  void save_amrex_file();

  void read_param(const std::string& command, ReadParam& param);

  void post_process_param(bool receiveICOnly);

  virtual void post_regrid() {
    // Array redistribution creates fine nodes before their fluid state can be
    // interpolated from the retained coarse level.
    distribute_arrays();
    fill_new_cells();
  }

  void distribute_arrays();
  void fill_new_cells();

  int count_couple_node_number();

  int loop_through_node(std::string action, double* const pos_DI = nullptr,
                        const double* const data = nullptr,
                        const int* const index = nullptr);

  void find_mpi_rank_for_points(const int nPoint, const double* const xyz_I,
                                int* const rank_I);

  void get_couple_node_loc(double* const pos_DI);

  void set_node_fluid(const double* const data, const int* const index,
                      const std::vector<std::string>& names);

  void set_node_fluid();

  void set_node_fluid(const FluidInterface& other);

  void set_node_fluid_to_zero() {
    for (int iLev = 0; iLev < n_lev(); ++iLev)
      nodeFluid[iLev].setVal(0.0);
  }

  void calc_current();

  void normalize_fluid_variables();

  // Print the nodes where a fluid density is zero (see
  // convert_moment_to_velocity).
  void report_zero_density(int iLev, bool phyNodeOnly) const;

  void convert_moment_to_velocity(bool phyNodeOnly = false, bool doWarn = true);

  void set_plasma_charge_and_mass(amrex::Real qomEl);

  // Publish the immutable normParams; scalarOnly skips per-variable factors.
  void finalize_normalization(bool scalarOnly = false);

  void analyze_var_names(bool useNeutral = false);

  void calc_fluid_state(const double* dataPIC_I, double* dataFluid_I) const;
  void calc_fluid_state_hybrid(const double* dataPIC_I,
                               double* dataFluid_I) const;

  void print_info() const;

  void get_for_points(const int nDim, const int nPoint,
                      const double* const xyz_I, double* const data_I,
                      const int nVar, const double coef = 1,
                      const amrex::Vector<int>& idxMap = {});

  void get_moments_for_points(const int nDim, const int nPoint,
                              const double* const xyz_I, double* const data_I,
                              const int nVar, const double coef) {
    amrex::Vector<int> idxMap;
    idxMap.reserve(5 * nFluid);
    for (int i = 0; i < nFluid; i++) {
      idxMap.push_back(iRho_I[i]);
      idxMap.push_back(iRhoUx_I[i]);
      idxMap.push_back(iRhoUy_I[i]);
      idxMap.push_back(iRhoUz_I[i]);
      idxMap.push_back(iP_I[i]);
    }

    if (nVar != idxMap.size()) {
      amrex::Print() << "nVar = " << nVar
                     << " idxMap.size() = " << idxMap.size() << std::endl;
      amrex::Abort("Error: nVar != idxMap.size(). It is likely BATSRUS was "
                   "compiled with incorrect equation module.");
    }
    get_for_points(nDim, nPoint, xyz_I, data_I, nVar, coef, idxMap);
  }

  int get_nCellPerPatch() const { return nCellPerPatch; }

  bool get_UseAnisoP() const { return (useAnisoP); }

  bool get_useElectronFluid() const { return useElectronFluid; }

  bool get_useHybridPIC() const { return useHybridPIC; }

  bool get_useMhdPe() const { return useMhdPe; }
  int get_iPe() const { return iPe; }

  double get_MoMi(int is) const {
    return (is >= 0 && is < static_cast<int>(MoMi_S.size())) ? MoMi_S[is] : 1.0;
  }
  double get_QoQi(int is) const {
    return (is >= 0 && is < static_cast<int>(QoQi_S.size())) ? QoQi_S[is] : 1.0;
  }

  template <typename Type>
  amrex::Real get_pe(const amrex::MFIter& mfi, const Type xyz,
                     const int iLev = 0) const {
    return FluidStateParams::get_pe(fluid_data(mfi, iLev), xyz);
  }

  int get_fluid_dimension() const { return (nDimFluid); }

  int get_nS() const { return nS; }

  int get_nFluid() const { return nFluid; }

  const amrex::Vector<std::string>& get_var_names() const { return varNames; }

  // Read-only accessors for the shared, frozen normalization state.
  double get_Si2NoL() const { return (normParams->Si2NoL); }
  double get_Si2NoM() const { return (1. / mNormSI); }
  double get_Si2NoRho() const { return normParams->Si2NoRho; }
  double get_Si2NoV() const { return normParams->Si2NoV; }
  double get_Si2NoB() const { return normParams->Si2NoB; }
  double get_Si2NoP() const { return normParams->Si2NoP; }

  double get_lnorm_si() const { return lNormSI; }
  double get_unorm_si() const { return uNormSI; }
  double get_mnorm_si() const { return mNormSI; }

  double get_cLight_SI() const { return uNormSI; }

  double get_rPlanet_SI() const { return rPlanetSi; }

  int get_scaling_factor() const { return ScalingFactor; }

  double get_MhdNo2SiL() const { return (MhdNo2SiL); }

  double get_Si2No_V(int idx) const { return normParams->Si2No_V[idx]; }

  double get_No2Si_V(int idx) const { return normParams->No2Si_V[idx]; }
  double get_No2SiL() const { return (normParams->No2SiL); }
  double get_No2SiRho() const { return (1. / normParams->Si2NoRho); }
  double get_No2SiV() const { return (1. / normParams->Si2NoV); }
  double get_No2SiB() const { return (1. / normParams->Si2NoB); }
  double get_No2SiE() const {
    return (normParams->Si2NoE > 0.0 ? 1. / normParams->Si2NoE : 1.0);
  }
  double get_No2SiP() const { return (1. / normParams->Si2NoP); }
  double get_No2SiJ() const { return (1. / normParams->Si2NoJ); }
  double get_No2SiM() const { return mNormSI; }

  double get_species_mass(int i) const { return MoMi_S[i]; }
  double get_species_charge(int i) const { return QoQi_S[i]; }

  // Derived temperature + MHD-length conversions.
  double get_Si2NoT() const { return normParams->si2noT(); }
  double get_No2SiT() const { return normParams->no2siT(); }
  double get_MhdNo2NoL() const { return (MhdNo2SiL * normParams->Si2NoL); }

  void sum_boundary() {
    timing_func("FI::sum_boundary");
    for (int iLev = 0; iLev < n_lev(); ++iLev)
      ::sum_boundary(nodeFluid[iLev], Geom(iLev).periodicity());

    for (int iLev = n_lev() - 2; iLev >= 0; iLev--) {
      sum_two_lev_interface_node(
          nodeFluid[iLev], nodeFluid[iLev + 1], 0, nodeFluid[iLev].nComp(),
          refRatio(iLev), Geom(iLev), Geom(iLev + 1), node_status(iLev + 1));
    }

    // Correct domain edge nodes
    for (int iLev = 0; iLev < n_lev() - 1; iLev++) {
      interp_from_coarse_to_fine_for_domain_edge(
          nodeFluid[iLev], nodeFluid[iLev + 1], 0, nodeFluid[iLev].nComp(),
          refRatio(iLev), Geom(iLev), Geom(iLev + 1), node_status(iLev + 1));
    }
  }

  void fill_boundary() {
    timing_func("FI::fill_boundary");
    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      if (!nodeFluid[iLev].empty())
        nodeFluid[iLev].FillBoundary(Geom(iLev).periodicity());
    }
  }

  virtual int get_neu_source_region(const amrex::MFIter& mfi,
                                    const amrex::IntVect ijk,
                                    const int iLev) const {
    return -1;
  }

  void set_resistivity(double etaSIIn) {
    // In SI unit R = u_si*L_si/eta_si, where eta_si is magnetic
    // diffusivity with unit m^2/s. In normalized CGS unit R =
    // u_pic*L_pic/(eta_pic/4pi), where eta_pic is also magnetic
    // diffusivity. Magnetic Reynolds number R should not change, so
    // eta_pic = 4*pi*eta_si*Si2NoV*Si2NoL.
    etaSI = etaSIIn;
    useResist = etaSI > 0;
    if (useResist)
      etaNO = fourPI * etaSI * normParams->Si2NoV * normParams->Si2NoL;
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
    if (is_grid_empty())
      return;

    std::string restartDir = get_restart_out_dir();

    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      amrex::VisMF::Write(nodeFluid[iLev], restartDir + gridName +
                                               "_Interface_nodeFluid" +
                                               lev_string(iLev));
      amrex::VisMF::Write(centerB[iLev], restartDir + gridName +
                                             "_Interface_centerB" +
                                             lev_string(iLev));
    }
  }

  void read_restart() {
    if (is_grid_empty())
      return;

    std::string restartDir = component + "/restartIN/";

    for (int iLev = 0; iLev < n_lev(); ++iLev) {
      amrex::VisMF::Read(nodeFluid[iLev], restartDir + gridName +
                                              "_Interface_nodeFluid" +
                                              lev_string(iLev));
      amrex::VisMF::Read(centerB[iLev], restartDir + gridName +
                                            "_Interface_centerB" +
                                            lev_string(iLev));
    }
    isnodeFluidReady = true;
  }

  void add_rho_to_loc(const amrex::Real& val, const amrex::MFIter& mfi,
                      const amrex::RealVect xyz, const int iFluid,
                      const int iLev) {
    add_to_mf(val, nodeFluid[iLev], mfi, Geom(iLev), xyz, iRho_I[iFluid]);
  }

  void add_mx_to_loc(const amrex::Real& val, const amrex::MFIter& mfi,
                     const amrex::RealVect xyz, const int iFluid,
                     const int iLev) {
    add_to_mf(val, nodeFluid[iLev], mfi, Geom(iLev), xyz, iRhoUx_I[iFluid]);
  }

  void add_my_to_loc(const amrex::Real& val, const amrex::MFIter& mfi,
                     const amrex::RealVect xyz, const int iFluid,
                     const int iLev) {
    add_to_mf(val, nodeFluid[iLev], mfi, Geom(iLev), xyz, iRhoUy_I[iFluid]);
  }

  void add_mz_to_loc(const amrex::Real& val, const amrex::MFIter& mfi,
                     const amrex::RealVect xyz, const int iFluid,
                     const int iLev) {
    add_to_mf(val, nodeFluid[iLev], mfi, Geom(iLev), xyz, iRhoUz_I[iFluid]);
  }

  void add_p_to_loc(const amrex::Real& val, const amrex::MFIter& mfi,
                    const amrex::RealVect xyz, const int iFluid,
                    const int iLev) {
    add_to_mf(val, nodeFluid[iLev], mfi, Geom(iLev), xyz, iP_I[iFluid]);
  }

  void add_to_loc(const amrex::Real& val, const amrex::MFIter& mfi,
                  const amrex::RealVect xyz, const int iVar, const int iLev) {
    add_to_mf(val, nodeFluid[iLev], mfi, Geom(iLev), xyz, iVar);
  }

  const amrex::MultiFab& get_center_b(const int iLev) const {
    return centerB[iLev];
  }

  amrex::Real get_center_b(const amrex::MFIter& mfi, const amrex::IntVect ijk,
                           const int iDir, const int iLev = 0) const {
    const auto& arr = centerB[iLev][mfi].array();
    return arr(ijk, iDir);
  }

  // The node fluid of the FAB mfi of level iLev and its geometry.
  // The node fluid of level iLev (e.g. to stage it on the host).
  amrex::MultiFab& get_node_fluid(int iLev) { return nodeFluid[iLev]; }

  FluidData fluid_data(const amrex::MFIter& mfi, const int iLev) const {
    FluidData d;
    d.fluid = nodeFluid[iLev].const_array(mfi);
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      d.plo[i] = Geom(iLev).ProbLo(i);
      d.invDx[i] = Geom(iLev).InvCellSize(i);
    }
    return d;
  }

  // Everything a kernel needs to evaluate the fluid state of the FAB mfi of
  // level iLev. Capture it by value.
  FluidView fluid_view(const amrex::MFIter& mfi, const int iLev) const {
    FluidView v;
    static_cast<FluidStateParams&>(v) =
        static_cast<const FluidStateParams&>(*this);
    v.d = fluid_data(mfi, iLev);
    return v;
  }

  amrex::Real get_value(const amrex::MFIter& mfi, const amrex::IntVect ijk,
                        const int iVar, const int iLev = 0) const {
    return fluid_data(mfi, iLev).value(ijk, iVar);
  }

  amrex::Real get_value(const amrex::MFIter& mfi, const amrex::RealVect xyz,
                        const int iVar, const int iLev = 0) const {
    return fluid_data(mfi, iLev).value(xyz, iVar);
  }

  // Host accessors of the fluid state. They read nodeFluid on the host, so
  // a GPU build may only call them where nodeFluid is host accessible; kernels
  // use fluid_view() and the FluidStateParams methods instead.
  template <typename T>
  amrex::Real get_number_density(const amrex::MFIter& mfi, const T xyz,
                                 const int is, const int iLev = 0) const {
    return FluidStateParams::get_number_density(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_u(const amrex::MFIter& mfi, const Type xyz, const int is,
                    const KernelArray<int>& iU_I, const int iJ,
                    const int iLev = 0) const {
    return FluidStateParams::get_u(fluid_data(mfi, iLev), xyz, is, iU_I, iJ);
  }

  template <typename T>
  amrex::Real get_fluid_mass_density(const amrex::MFIter& mfi, const T xyz,
                                     const int is, const int iLev = 0) const {
    return FluidStateParams::get_fluid_mass_density(fluid_data(mfi, iLev), xyz,
                                                    is);
  }

  template <typename Type>
  amrex::Real get_fluid_p(const amrex::MFIter& mfi, const Type xyz,
                          const int is, const int iLev = 0) const {
    return FluidStateParams::get_fluid_p(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_fluid_uth(const amrex::MFIter& mfi, const Type xyz,
                            const int is, const int iLev = 0) const {
    return FluidStateParams::get_fluid_uth(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_fluid_ux(const amrex::MFIter& mfi, const Type xyz,
                           const int is, const int iLev = 0) const {
    return FluidStateParams::get_fluid_ux(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_fluid_uy(const amrex::MFIter& mfi, const Type xyz,
                           const int is, const int iLev = 0) const {
    return FluidStateParams::get_fluid_uy(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_fluid_uz(const amrex::MFIter& mfi, const Type xyz,
                           const int is, const int iLev = 0) const {
    return FluidStateParams::get_fluid_uz(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_ux(const amrex::MFIter& mfi, const Type xyz, const int is,
                     const int iLev = 0) const {
    return FluidStateParams::get_ux(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_uy(const amrex::MFIter& mfi, const Type xyz, const int is,
                     const int iLev = 0) const {
    return FluidStateParams::get_uy(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_uz(const amrex::MFIter& mfi, const Type xyz, const int is,
                     const int iLev = 0) const {
    return FluidStateParams::get_uz(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_ppar(const amrex::MFIter& mfi, const Type xyz, const int is,
                       const int iLev = 0) const {
    return FluidStateParams::get_ppar(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_p(const amrex::MFIter& mfi, const Type xyz, const int is,
                    const int iLev = 0) const {
    return FluidStateParams::get_p(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_pxx(const amrex::MFIter& mfi, const Type xyz, const int is,
                      const int iLev = 0) const {
    return FluidStateParams::get_pxx(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_pyy(const amrex::MFIter& mfi, const Type xyz, const int is,
                      const int iLev = 0) const {
    return FluidStateParams::get_pyy(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_pzz(const amrex::MFIter& mfi, const Type xyz, const int is,
                      const int iLev = 0) const {
    return FluidStateParams::get_pzz(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_pxy(const amrex::MFIter& mfi, const Type xyz, const int is,
                      const int iLev = 0) const {
    return FluidStateParams::get_pxy(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_pxz(const amrex::MFIter& mfi, const Type xyz, const int is,
                      const int iLev = 0) const {
    return FluidStateParams::get_pxz(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_pyz(const amrex::MFIter& mfi, const Type xyz, const int is,
                      const int iLev = 0) const {
    return FluidStateParams::get_pyz(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  amrex::Real get_uth_iso(const amrex::MFIter& mfi, const Type xyz,
                          const int is, const int iLev = 0) const {
    return FluidStateParams::get_uth_iso(fluid_data(mfi, iLev), xyz, is);
  }

  template <typename Type>
  void set_particle_uth_iso(const int iLev, const amrex::MFIter& mfi,
                            const Type xyz, double* u, double* v, double* w,
                            const double rand1, const double rand2,
                            const double rand3, const double rand4,
                            const int is, const double uthIn = -1) const {
    FluidStateParams::set_particle_uth_iso(fluid_data(mfi, iLev), xyz, u, v, w,
                                           rand1, rand2, rand3, rand4, is, -1);
  }

  template <typename Type>
  void set_particle_uth_aniso(const int iLev, const amrex::MFIter& mfi,
                              const Type xyz, double* u, double* v, double* w,
                              const double rand1, const double rand2,
                              const double rand3, const double rand4,
                              const int is, const double uthParIn = -1,
                              const double uthPerpIn = -1) const {
    FluidStateParams::set_particle_uth_aniso(fluid_data(mfi, iLev), xyz, u, v,
                                             w, rand1, rand2, rand3, rand4, is,
                                             -1, -1);
  }

  template <typename Type>
  amrex::Real get_bx(const amrex::MFIter& mfi, const Type xyz,
                     const int iLev = 0) const {
    return FluidStateParams::get_bx(fluid_data(mfi, iLev), xyz);
  }

  template <typename Type>
  amrex::Real get_by(const amrex::MFIter& mfi, const Type xyz,
                     const int iLev = 0) const {
    return FluidStateParams::get_by(fluid_data(mfi, iLev), xyz);
  }

  template <typename Type>
  amrex::Real get_bz(const amrex::MFIter& mfi, const Type xyz,
                     const int iLev = 0) const {
    return FluidStateParams::get_bz(fluid_data(mfi, iLev), xyz);
  }

  template <typename Type>
  amrex::Real get_ex(const amrex::MFIter& mfi, const Type xyz,
                     const int iLev = 0) const {
    return FluidStateParams::get_ex(fluid_data(mfi, iLev), xyz);
  }

  template <typename Type>
  amrex::Real get_ey(const amrex::MFIter& mfi, const Type xyz,
                     const int iLev = 0) const {
    return FluidStateParams::get_ey(fluid_data(mfi, iLev), xyz);
  }

  template <typename Type>
  amrex::Real get_ez(const amrex::MFIter& mfi, const Type xyz,
                     const int iLev = 0) const {
    return FluidStateParams::get_ez(fluid_data(mfi, iLev), xyz);
  }

  // Calculate grad(pe) at node (x,y,z). If this node is at the boundary of
  // the fab, it will return grad(pe) at a node that is one cell away from the
  // boundary. Only works when useElectronFluid is False.
  amrex::Real get_grad_pe_x(const amrex::MFIter& mfi, amrex::IntVect ijk,
                            const int iLev) const {
    const amrex::Box& box = mfi.fabbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    if (ijk[ix_] == lo.x) {
      ijk[ix_] += 1;
    } else if (ijk[ix_] == hi.x) {
      ijk[ix_] -= 1;
    }

    amrex::IntVect ijkL = ijk;
    amrex::IntVect ijkR = ijk;

    ijkL[ix_] += 1;
    ijkR[ix_] -= 1;

    amrex::Real gradpe =
        0.5 * Geom(iLev).InvCellSize(ix_) *
        (get_p(mfi, ijkL, 0, iLev) - get_p(mfi, ijkR, 0, iLev));

    return gradpe;
  }

  amrex::Real get_grad_pe_y(const amrex::MFIter& mfi, amrex::IntVect ijk,
                            const int iLev) const {

    const amrex::Box& box = mfi.fabbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    if (ijk[iy_] == lo.y) {
      ijk[iy_] += 1;
    } else if (ijk[iy_] == hi.y) {
      ijk[iy_] -= 1;
    }

    amrex::IntVect ijkL = ijk;
    amrex::IntVect ijkR = ijk;

    ijkL[iy_] += 1;
    ijkR[iy_] -= 1;

    amrex::Real gradpe =
        0.5 * Geom(iLev).InvCellSize(iy_) *
        (get_p(mfi, ijkL, 0, iLev) - get_p(mfi, ijkR, 0, iLev));

    return gradpe;
  }

  amrex::Real get_grad_pe_z(const amrex::MFIter& mfi, amrex::IntVect ijk,
                            const int iLev) const {
    if (nDim <= 2)
      return 0.0;

    const amrex::Box& box = mfi.fabbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    if (ijk[iz_] == lo.z) {
      ijk[iz_] += 1;
    } else if (ijk[iz_] == hi.z) {
      ijk[iz_] -= 1;
    }

    amrex::IntVect ijkL = ijk;
    amrex::IntVect ijkR = ijk;

    ijkL[iz_] += 1;
    ijkR[iz_] -= 1;

    amrex::Real gradpe =
        0.5 * Geom(iLev).InvCellSize(iz_) *
        (get_p(mfi, ijkL, 0, iLev) - get_p(mfi, ijkR, 0, iLev));

    return gradpe;
  }
};
#endif
