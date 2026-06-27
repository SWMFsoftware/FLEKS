#ifndef _EXOSOURCE_H_
#define _EXOSOURCE_H_

#include "SourceInterface.h"

// #define _EXOSPHERE_

#ifdef _EXOSPHERE_
// Get source from MHD side.
extern "C" {
void get_source_wrapper(double xyzSI[3], double sourceSI[6]);
}
#else
inline void get_source_wrapper(double xyzSI[3], double sourceSI[6]) {
  for (int i = 0; i < 6; ++i) sourceSI[i] = 0.0;
}
#endif

class UserSource : public SourceInterface {
public:
  UserSource(const FluidInterface& other, int id, std::string tag,
             FluidType typeIn = SourceFluid)
      : SourceInterface(other, id, tag, typeIn) {
    info = "Exosphere Source";
    useFluidSource = true;
  }

  // void sum_to_single_source() override {
  //   for (int iLev = 0; iLev < n_lev(); iLev++) {
  //     if (!nodeFluid[iLev].empty()) {
  //       for (amrex::MFIter mfi(nodeFluid[iLev]); mfi.isValid(); ++mfi) {
  //         const amrex::Box& box = mfi.fabbox();
  //         const auto lo = lbound(box);
  //         const auto hi = ubound(box);

  //         const amrex::Array4<amrex::Real>& arr = nodeFluid[iLev][mfi].array();

  //         for (int k = lo.z; k <= hi.z; ++k)
  //           for (int j = lo.y; j <= hi.y; ++j)
  //             for (int i = lo.x; i <= hi.x; ++i) {

  //               auto sum_moment = [&, this](amrex::Vector<int>& idx) {
  //                 amrex::Real sum = 0;
  //                 for (int ii = 0; ii < idx.size(); ++ii) {
  //                   sum += arr(i, j, k, idx[ii]);
  //                 }
  //                 for (int ii = 0; ii < idx.size(); ++ii) {
  //                   arr(i, j, k, idx[ii]) = sum;
  //                 }
  //               };

  //               sum_moment(iRho_I);
  //               sum_moment(iRhoUx_I);
  //               sum_moment(iRhoUy_I);
  //               sum_moment(iRhoUz_I);
  //               sum_moment(iP_I);
  //             }
  //       }
  //     }
  //   }
  // }

  // ---- Exosphere density profiles ----

  double get_exosphere_density(double r) const override {
    if (exosphereType == "None") return 0.0;
    if (r < rPlanetSi) return 0.0;

    double sum = 0.0;
    if (exosphereType == "Exponential") {
      for (int i = 0; i < nExoComponent; ++i) {
        if (exoH0[i] > 0.0) {
          sum += exoN0[i] * exp(-(r - rPlanetSi) / exoH0[i]);
        }
      }
    } else if (exosphereType == "Power-Law") {
      for (int i = 0; i < nExoComponent; ++i) {
        if (r > 0.0) {
          sum += exoN0[i] * pow(rPlanetSi / r, exoK0[i]);
        }
      }
    } else if (exosphereType == "Chamberlain") {
      for (int i = 0; i < nExoComponent; ++i) {
        if (rPlanetSi > 0.0 && r > 0.0) {
          sum += exoN0[i] * exp(-exoH0[i] * (1.0 / rPlanetSi - 1.0 / r));
        }
      }
    }
    return sum;
  }

  double get_exosphere_component_density(double r,
                                         int iC) const override {
    if (exosphereType == "None") return 0.0;
    if (iC < 0 || iC >= nExoComponent) return 0.0;
    if (r < rPlanetSi) return 0.0;

    if (exosphereType == "Exponential") {
      if (exoH0[iC] > 0.0) {
        return exoN0[iC] * exp(-(r - rPlanetSi) / exoH0[iC]);
      }
    } else if (exosphereType == "Power-Law") {
      if (r > 0.0) {
        return exoN0[iC] * pow(rPlanetSi / r, exoK0[iC]);
      }
    } else if (exosphereType == "Chamberlain") {
      if (rPlanetSi > 0.0 && r > 0.0) {
        return exoN0[iC] * exp(-exoH0[iC] * (1.0 / rPlanetSi - 1.0 / r));
      }
    }
    return 0.0;
  }

  //-------------------------------------------------------------------
  // Voronov 1997 electron impact ionization rate coefficient.
  // Input: Te in eV.  Returns <sigma*v> in [m^3/s].
  amrex::Real impact_rate(amrex::Real Te_eV, int iC) const {
    if (Te_eV <= 0.0) return 0.0;
    double Ei = impactEIon[iC];  // ionization energy [eV]
    double A = impactA[iC];      // Voronov A [cm^3/s]
    double K = impactK[iC];      // Voronov K
    double X = impactX[iC];      // Voronov X
    double t = Te_eV / Ei;
    // <sigma v> in cm^3/s, convert to m^3/s
    return A * pow(t, K) / (X + t) * exp(-Ei / Te_eV) * 1e-6;
  }

  //-------------------------------------------------------------------
  // Charge exchange rate coefficient (constant cross-section model).
  // u_mag_SI: relative ion-neutral speed in [m/s].
  // Returns <sigma*v> in [m^3/s].
  amrex::Real cx_rate(amrex::Real u_mag_SI, int iC) const {
    // sigmaCX in [cm^2], convert to [m^2]
    return cxSigma[iC] * 1e-4 * u_mag_SI;
  }

  //-------------------------------------------------------------------
  // Check whether a point at (x, y, z) [m] relative to planet center
  // lies inside the planetary shadow cylinder.
  bool is_in_shadow(double x, double y, double z) const {
    if (!useShadowCylinder) return false;
    // Project position vector onto solar direction.
    // solarDir points from planet center toward the Sun.
    double proj = x * solarDir[0] + y * solarDir[1] + z * solarDir[2];
    if (proj >= 0.0) return false;  // dayside — not in shadow
    // Finite-height check: only shadow points within the cylinder height.
    if (-proj > shadowCylinderHalfHeight)
      return false;
    // Perpendicular distance from the sun–planet axis.
    double r2 = x * x + y * y + z * z;
    double perp2 = r2 - proj * proj;
    return perp2 <= shadowCylinderRadius * shadowCylinderRadius;
  }

  //-------------------------------------------------------------------
  // Photoionation rate at position (x, y, z) [m] from planet center.
  // Returns zero inside the shadow cylinder (nightside).
  amrex::Real photo_rate(double x, double y, double z, int iC) const {
    if (is_in_shadow(x, y, z)) return 0.0;
    double r2 = x * x + y * y + z * z;
    double r_m = sqrt(r2);
    double ratio = rPlanetSi / r_m;
    return photoNu0[iC] * ratio * ratio;
  }

  //-------------------------------------------------------------------
  // Compute electron temperature in eV from PIC-normalized pressure
  // and number density.  pe: electron pressure (PIC units).
  // ne: electron number density (PIC units).
  amrex::Real electron_temperature(amrex::Real pe, amrex::Real ne) const {
    if (ne <= 0.0) return 0.0;
    // T_eV = (m_p/e) * uNorm^2 * (pe/ne)
    const amrex::Real protonMassPerCharge =
        cProtonMassSI / cUnitChargeSI;   // [kg/C]
    const amrex::Real ur2 =
        get_unorm_si() * get_unorm_si(); // [m^2/s^2]
    return protonMassPerCharge * ur2 * pe / ne;
  }

  //-------------------------------------------------------------------
  // Read ionization-related parameter commands from PARAM.in.
  void read_param(const std::string& command, ReadParam& param) override {
    if (command == "#PHOTOIONIZATION") {
      usePhotoIonization = true;
      param.read_var("nComponent", nPhotoComponent);
      photoNu0.resize(nPhotoComponent);
      for (int i = 0; i < nPhotoComponent; ++i) {
        param.read_var("nuPhoto0", photoNu0[i]);
      }
    } else if (command == "#ELECTRONIMPACT") {
      useElectronImpact = true;
      param.read_var("nComponent", nImpactComponent);
      impactEIon.resize(nImpactComponent);
      impactA.resize(nImpactComponent);
      impactK.resize(nImpactComponent);
      impactX.resize(nImpactComponent);
      for (int i = 0; i < nImpactComponent; ++i) {
        param.read_var("eIon", impactEIon[i]);
        param.read_var("Acoeff", impactA[i]);
        param.read_var("Kcoeff", impactK[i]);
        param.read_var("Xcoeff", impactX[i]);
      }
    } else if (command == "#CHARGEEXCHANGE") {
      useChargeExchange = true;
      param.read_var("nComponent", nCXComponent);
      cxSigma.resize(nCXComponent);
      for (int i = 0; i < nCXComponent; ++i) {
        param.read_var("sigmaCX", cxSigma[i]);
      }
    } else if (command == "#SHADOWCYLINDER") {
      useShadowCylinder = true;
      param.read_var("solarDirX", solarDir[0]);
      param.read_var("solarDirY", solarDir[1]);
      param.read_var("solarDirZ", solarDir[2]);
      param.read_var("radius", shadowCylinderRadius);
      param.read_var("halfHeight", shadowCylinderHalfHeight);
      // Normalize the solar direction vector.
      double norm = sqrt(solarDir[0] * solarDir[0] +
                         solarDir[1] * solarDir[1] +
                         solarDir[2] * solarDir[2]);
      if (norm > 0.0) {
        solarDir[0] /= norm;
        solarDir[1] /= norm;
        solarDir[2] /= norm;
      }
    }
  }

  //-------------------------------------------------------------------
  // Validate consistency between exosphere and ionization commands.
  void post_process_param() override {
    if (exosphereType == "None") return;

    // The source code assumes species 0 is the electron (used for ne and Te
    // in electron impact ionization, and skipped in the exosphere-to-ion
    // mapping where component iC -> species iC+1).  Enforce this here.
    if (nS < 1) {
      amrex::Print() << printPrefix << "Error: no plasma species defined. "
                     << "Use #PLASMA to set species.\n";
      std::abort();
    }
    if (QoQi_S[0] >= 0.0) {
      amrex::Print() << printPrefix << "Error: species 0 must be the electron "
                     << "(negative charge). Got Q/Qi[0] = " << QoQi_S[0]
                     << ". Reorder #PLASMA so the electron is first.\n";
      std::abort();
    }
    if (nExoComponent > nS - 1) {
      amrex::Print() << printPrefix << "Error: #EXOSPHERE nComponent ("
                     << nExoComponent << ") exceeds the number of ion species ("
                     << nS - 1 << "). Add more ion species in #PLASMA.\n";
      std::abort();
    }

    if (usePhotoIonization && nPhotoComponent != nExoComponent) {
      amrex::Print() << printPrefix << "Error: #PHOTOIONIZATION nComponent ("
                     << nPhotoComponent << ") != #EXOSPHERE nComponent ("
                     << nExoComponent << ")\n";
      std::abort();
    }
    if (useElectronImpact && nImpactComponent != nExoComponent) {
      amrex::Print() << printPrefix << "Error: #ELECTRONIMPACT nComponent ("
                     << nImpactComponent << ") != #EXOSPHERE nComponent ("
                     << nExoComponent << ")\n";
      std::abort();
    }
    if (useChargeExchange && nCXComponent != nExoComponent) {
      amrex::Print() << printPrefix << "Error: #CHARGEEXCHANGE nComponent ("
                     << nCXComponent << ") != #EXOSPHERE nComponent ("
                     << nExoComponent << ")\n";
      std::abort();
    }
  }

  //-------------------------------------------------------------------
  // Set nodeFluid from plasma-state-dependent ionization processes.
  void set_source(const FluidInterface& other) override {
    std::string nameFunc = "FS:get_source_from_fluid";
    amrex::Print() << nameFunc << " is called.";

    set_node_fluid(other);

    // Global NODE box.
    const amrex::Box gbx = convert(Geom(0).Domain(), { AMREX_D_DECL(1, 1, 1) });

    const bool doPhoto = usePhotoIonization;
    const bool doImpact = useElectronImpact;
    const bool doCX = useChargeExchange;

    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (!nodeFluid[iLev].empty()) {
        for (amrex::MFIter mfi(nodeFluid[iLev]); mfi.isValid(); ++mfi) {
          const amrex::Real* dx = Geom(iLev).CellSize();
          const auto plo = Geom(iLev).ProbLo();

          // For each block, looping through all nodes, including ghost nodes.
          const amrex::Box& box = mfi.fabbox();
          const auto lo = lbound(box);
          const auto hi = ubound(box);

          const amrex::Array4<amrex::Real>& arr = nodeFluid[iLev][mfi].array();

          for (int k = lo.z; k <= hi.z; ++k)
            for (int j = lo.y; j <= hi.y; ++j)
              for (int i = lo.x; i <= hi.x; ++i) {
                amrex::IntVect idx = { AMREX_D_DECL(i, j, k) };
                for (int iDim = 0; iDim < nDim; iDim++) {
                  if (Geom(iLev).isPeriodic(iDim)) {
                    idx[iDim] = shift_periodic_index(
                        idx[iDim], gbx.smallEnd(iDim), gbx.bigEnd(iDim));
                  }
                }
                double xyz[3] = { 0, 0, 0 };
                for (int iDim = 0; iDim < get_fluid_dimension(); iDim++) {
                  xyz[iDim] = idx[iDim] * dx[iDim] + plo[iDim];
                }

                double source[6] = {0.0};
#ifdef _EXOSPHERE_
                get_source_wrapper(xyz, source);
#else
                double r_val = 0.0;
                for (int d = 0; d < 3; ++d) {
                  r_val += xyz[d] * xyz[d];
                }
                r_val = sqrt(r_val);

                // Pre-fetch plasma state once per cell (lazy evaluation)
                double ne = 0, pe = 0, Te_eV = 0;
                bool plasmaFetched = false;
                auto fetch_electron_plasma = [&, &other=other]() {
                  if (plasmaFetched) return;
                  ne = other.get_number_density(mfi, idx, 0, iLev);
                  pe = other.get_p(mfi, idx, 0, iLev);
                  Te_eV = electron_temperature(pe, ne);
                  plasmaFetched = true;
                };

                // Per-ion-species source accumulators.
                // Species 0 = electron (skipped), 1 = H+, 2 = O+, ...
                // Exosphere component iC maps to ion species (iC + 1).
                // srcRho[iS], srcP[iS] in SI [kg m^-3 s^-1, Pa s^-1].
                const int nIonS = nS - 1;  // number of ion species
                std::vector<double> srcRho(nIonS + 1, 0.0);
                std::vector<double> srcP(nIonS + 1, 0.0);

                for (int iC = 0; iC < nExoComponent; ++iC) {
                  // Map exosphere component iC to ion species (iC + 1).
                  // Species 0 = electron, so ion species index = iC + 1.
                  const int iSp = iC + 1;
                  if (iSp > nIonS) break;  // safety: not enough ion species

                  double dens_i =
                      get_exosphere_component_density(r_val, iC);
                  double nu_tot = 0.0;

                  // ---- Photoionization ----
                  if (doPhoto) {
                    nu_tot += photo_rate(xyz[0], xyz[1], xyz[2], iC);
                  }

                  // ---- Electron impact ionization ----
                  if (doImpact) {
                    fetch_electron_plasma();
                    if (ne > 0 && Te_eV > 0) {
                      nu_tot += ne * impact_rate(Te_eV, iC);
                    }
                  }

                  // ---- Charge exchange ----
                  // Use the *same* ion species that this component maps to.
                  if (doCX) {
                    double ni_i = other.get_number_density(mfi, idx, iSp, iLev);
                    double ux_i = other.get_ux(mfi, idx, iSp, iLev);
                    double uy_i = other.get_uy(mfi, idx, iSp, iLev);
                    double uz_i = other.get_uz(mfi, idx, iSp, iLev);
                    double u_mag_SI =
                        sqrt(ux_i * ux_i + uy_i * uy_i + uz_i * uz_i) *
                        get_unorm_si();
                    if (ni_i > 0 && u_mag_SI > 0) {
                      nu_tot += ni_i * cx_rate(u_mag_SI, iC);
                    }
                  }

                  // Number density production rate [m^-3 s^-1]
                  double S_n = dens_i * nu_tot;
                  // Convert to mass density production rate [kg m^-3 s^-1]
                  // using the mass of the corresponding ion species.
                  double mass_amu = get_species_mass(iSp);
                  srcRho[iSp] += S_n * mass_amu * cProtonMassSI;
                  // Pressure source at T_source = 100 K
                  srcP[iSp] += S_n * mass_amu * cProtonMassSI * cBoltzmannSI * 100.0;
                }

#endif

                for (int iFluid = 0; iFluid < nFluid; iFluid++) {
                  arr(i, j, k, iRho_I[iFluid]) = 0;
                  arr(i, j, k, iUx_I[iFluid]) = 0;
                  arr(i, j, k, iUy_I[iFluid]) = 0;
                  arr(i, j, k, iUz_I[iFluid]) = 0;
                  arr(i, j, k, iP_I[iFluid]) = 0;
                }

                // Write per-species source terms.  Zero momentum (source
                // particles are born at rest in the planet frame).
                bool anySource = false;
                for (int iSp = 1; iSp <= nIonS; ++iSp) {
                  if (srcRho[iSp] > 0) {
                    anySource = true;
                    arr(i, j, k, iRho_I[iSp]) =
                        srcRho[iSp] * Si2NoRho / get_Si2NoT();
                    arr(i, j, k, iUx_I[iSp]) = 0.0;
                    arr(i, j, k, iUy_I[iSp]) = 0.0;
                    arr(i, j, k, iUz_I[iSp]) = 0.0;
                    arr(i, j, k, iP_I[iSp]) =
                        srcP[iSp] * Si2NoP / get_Si2NoT();
                  }
                }
                // Electron pressure source: sum of all ion sources.
                if (anySource && iPe >= 0) {
                  double srcPe = 0.0;
                  for (int iSp = 1; iSp <= nIonS; ++iSp)
                    srcPe += srcP[iSp];
                  arr(i, j, k, iPe) = srcPe * Si2NoP / get_Si2NoT();
                }
              } // for k
        }
      }
    }

    if (!isGridEmpty && useCurrent) {
      // The current stored in nodeFluid are used to initialize particle
      // velocities. Since the source particles only contribute little to the
      // total current, set current to zero.
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        amrex::MultiFab currentMF(nodeFluid[iLev], amrex::make_alias, iJx, 3);
        currentMF.setVal(0, currentMF.nGrow());
      }
    }
  }
};

#endif
