#ifndef _EXOSOURCE_H_
#define _EXOSOURCE_H_

#include "SourceInterface.h"

class UserSource : public SourceInterface {
public:
  UserSource(const FluidInterface& other, int id, std::string tag,
             FluidType typeIn = SourceFluid)
      : SourceInterface(other, id, tag, typeIn) {
    info = "Exosphere Source";
    useFluidSource = true;
  }

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

  //=================================================================
  // Per-cell source term functions.
  //
  // Each of the three ionization processes is implemented as a function
  // operating on a single spatial cell.  They return the ionization
  // frequency nu [s^-1] contributed by that process; multiplying by the
  // neutral number density gives the production rate S_n [m^-3 s^-1].
  //=================================================================

  //-------------------------------------------------------------------
  // Photoionization frequency [s^-1] at position xyz [m] (measured from
  // the planet center) for neutral component iC.  Returns zero inside
  // the planetary shadow cylinder (nightside).
  amrex::Real photoionization_rate(const double xyz[3], int iC) const {
    // Zero inside the planetary shadow cylinder (nightside).
    if (is_in_shadow(xyz[0], xyz[1], xyz[2])) return 0.0;
    double r2 = xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
    double r_m = sqrt(r2);
    double ratio = rPlanetSi / r_m;
    // Geometric dilution: nu = nu0 * (r_planet / r)^2
    return photoNu0[iC] * ratio * ratio;
  }

  //-------------------------------------------------------------------
  // Electron-impact ionization frequency [s^-1] for neutral component
  // iC.  ne: electron number density [m^-3]; Te_eV: electron temperature
  // [eV].  Returns ne * <sigma*v>_impact, or zero when the plasma state
  // is not positive.
  amrex::Real impact_ionization_rate(amrex::Real ne, amrex::Real Te_eV,
                                     int iC) const {
    if (ne <= 0.0 || Te_eV <= 0.0) return 0.0;
    // Voronov 1997 fit: <sigma*v> in cm^3/s, converted to m^3/s.
    double Ei = impactEIon[iC];  // ionization energy [eV]
    double A = impactA[iC];      // Voronov A [cm^3/s]
    double K = impactK[iC];      // Voronov K
    double X = impactX[iC];      // Voronov X
    double t = Te_eV / Ei;
    double sigma_v = A * pow(t, K) / (X + t) * exp(-Ei / Te_eV) * 1e-6;
    return ne * sigma_v;
  }

  //-------------------------------------------------------------------
  // Charge-exchange frequency [s^-1] for neutral component iC reacting
  // with ion species iSp.  The ion density and bulk velocity are fetched
  // from the plasma state (other) at cell (mfi, idx, iLev).  The
  // cross-section is looked up from the matrix cxSigma[iC, iIon] where
  // iIon = iSp - 1 is the 0-based ion index.
  //
  // TODO: Use the ion-neutral *relative* velocity (accounting for the
  // neutral bulk flow) instead of the ion speed alone.
  amrex::Real charge_exchange_rate(const FluidInterface& other,
                                   const amrex::MFIter& mfi,
                                   const amrex::IntVect& idx, int iLev,
                                   int iSp, int iC) const {
    int iIon = iSp - 1;  // 0-based ion index for matrix lookup
    double ni = other.get_number_density(mfi, idx, iSp, iLev);
    double ux_i = other.get_ux(mfi, idx, iSp, iLev);
    double uy_i = other.get_uy(mfi, idx, iSp, iLev);
    double uz_i = other.get_uz(mfi, idx, iSp, iLev);
    double u_mag_SI =
        sqrt(ux_i * ux_i + uy_i * uy_i + uz_i * uz_i) * get_unorm_si();
    if (ni <= 0.0 || u_mag_SI <= 0.0) return 0.0;
    // Constant cross-section model: sigmaCX in [cm^2], convert to [m^2].
    double sigma = cxSigma[iC * nCXIonSpecies + iIon];
    return ni * sigma * 1e-4 * u_mag_SI;
  }

  //-------------------------------------------------------------------
  // Require #EXOSPHERE to have been read (with nComponent > 0) before an
  // ionization command is parsed.  Domain syncs nExoComponent from fi into
  // this object just before dispatching ionization commands, so a missing or
  // out-of-order #EXOSPHERE leaves nExoComponent == 0 and is caught here.
  void require_exosphere(const std::string& command) const {
    if (nExoComponent > 0) return;
    amrex::Abort(printPrefix + "Error: " + command + " requires "
                 + "#EXOSPHERE with nComponent > 0 to be specified before "
                 + command + ".");
  }

  //-------------------------------------------------------------------
  // Read ionization-related parameter commands from PARAM.in.
  void read_param(const std::string& command, ReadParam& param) override {
    if (command == "#PHOTOIONIZATION") {
      usePhotoIonization = true;
      // The neutral component count is inherited from #EXOSPHERE.
      require_exosphere(command);
      photoNu0.resize(nExoComponent);
      for (int i = 0; i < nExoComponent; ++i) {
        param.read_var("nuPhoto0", photoNu0[i]);
      }
    } else if (command == "#ELECTRONIMPACT") {
      useElectronImpact = true;
      // The neutral component count is inherited from #EXOSPHERE.
      require_exosphere(command);
      impactEIon.resize(nExoComponent);
      impactA.resize(nExoComponent);
      impactK.resize(nExoComponent);
      impactX.resize(nExoComponent);
      for (int i = 0; i < nExoComponent; ++i) {
        param.read_var("eIon", impactEIon[i]);
        param.read_var("Acoeff", impactA[i]);
        param.read_var("Kcoeff", impactK[i]);
        param.read_var("Xcoeff", impactX[i]);
      }
    } else if (command == "#CHARGEEXCHANGE") {
      useChargeExchange = true;
      // The neutral component count is inherited from #EXOSPHERE.
      require_exosphere(command);
      param.read_var("nIonSpecies", nCXIonSpecies);
      cxSigma.resize(nExoComponent * nCXIonSpecies);
      for (int iC = 0; iC < nExoComponent; ++iC) {
        for (int iIon = 0; iIon < nCXIonSpecies; ++iIon) {
          param.read_var("sigmaCX", cxSigma[iC * nCXIonSpecies + iIon]);
        }
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
      amrex::Abort(printPrefix + "Error: no plasma species defined. "
                   + "Use #PLASMA to set species.");
    }
    if (QoQi_S[0] >= 0.0) {
      amrex::Abort(printPrefix + "Error: species 0 must be the electron "
                   + "(negative charge). Got Q/Qi[0] = "
                   + std::to_string(QoQi_S[0])
                   + ". Reorder #PLASMA so the electron is first.");
    }
    if (nExoComponent > nS - 1) {
      amrex::Abort(printPrefix + "Error: #EXOSPHERE nComponent ("
                   + std::to_string(nExoComponent)
                   + ") exceeds the number of ion species ("
                   + std::to_string(nS - 1)
                   + "). Add more ion species in #PLASMA.");
    }

    if (useChargeExchange) {
      // nExoComponent is used directly as the neutral component count for
      // charge exchange (no separate nComponent in #CHARGEEXCHANGE).  The
      // cross-section matrix was sized to nExoComponent * nCXIonSpecies in
      // read_param; only the ion-species count needs validating here.
      if (nCXIonSpecies != nS - 1) {
        amrex::Abort(printPrefix + "Error: #CHARGEEXCHANGE nIonSpecies ("
                     + std::to_string(nCXIonSpecies)
                     + ") != number of ion species ("
                     + std::to_string(nS - 1) + ")");
      }
    }

    // The neutral temperature exoT0[iC] is used as the source temperature
    // for the pressure source term in set_source.
    for (int iC = 0; iC < nExoComponent; ++iC) {
      if (exoT0[iC] <= 0.0) {
        amrex::Abort(printPrefix + "Error: #EXOSPHERE T0 for component "
                     + std::to_string(iC) + " must be positive (got "
                     + std::to_string(exoT0[iC])
                     + " K). It is used as the ionization source temperature.");
      }
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

                  // Sum the per-cell ionization frequency nu [s^-1] from
                  // each enabled process.  Each *_rate() below operates on
                  // this single spatial cell.
                  double nu_tot = 0.0;
                  if (doPhoto)
                    nu_tot += photoionization_rate(xyz, iC);
                  if (doImpact) {
                    fetch_electron_plasma();
                    nu_tot += impact_ionization_rate(ne, Te_eV, iC);
                  }
                  if (doCX) {
                    // Sum the charge-exchange frequency over all ion
                    // species; each neutral can exchange with every ion.
                    for (int iSpCX = 1; iSpCX <= nIonS; ++iSpCX) {
                      nu_tot += charge_exchange_rate(other, mfi, idx, iLev,
                                                     iSpCX, iC);
                    }
                  }

                  // Number density production rate [m^-3 s^-1]
                  double S_n = dens_i * nu_tot;
                  // Convert to mass density production rate [kg m^-3 s^-1]
                  // using the mass of the corresponding ion species.
                  double mass_amu = get_species_mass(iSp);
                  srcRho[iSp] += S_n * mass_amu * cProtonMassSI;
                  // Pressure source using the neutral temperature
                  // exoT0[iC] [K] from the #EXOSPHERE profile.
                  srcP[iSp] +=
                      S_n * mass_amu * cProtonMassSI * cBoltzmannSI *
                      exoT0[iC];
                }

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
