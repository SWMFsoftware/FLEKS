#ifndef _EXOSOURCE_H_
#define _EXOSOURCE_H_

#include "SourceInterface.h"

class UserSource : public SourceInterface {
public:
  UserSource(const FluidInterface& other, int id, std::string tag,
             FluidType typeIn, const DomainParameters& dp)
      : SourceInterface(other, id, tag, typeIn, dp) {
    info = "Exosphere Source";
    useFluidSource = true;
  }

  // ---- Exosphere density profiles ----

  amrex::Real get_exosphere_density(amrex::Real r) const override {
    if (exosphereType == "None") return 0.0;
    if (r < get_rPlanet_SI()) return 0.0;

    amrex::Real sum = 0.0;
    if (exosphereType == "Exponential") {
      for (int i = 0; i < nExoComponent; ++i) {
        if (exoH0[i] > 0.0) {
          sum += exoN0[i] * exp(-(r - get_rPlanet_SI()) / exoH0[i]);
        }
      }
    } else if (exosphereType == "Power-Law") {
      for (int i = 0; i < nExoComponent; ++i) {
        if (r > 0.0) {
          sum += exoN0[i] * pow(get_rPlanet_SI() / r, exoK0[i]);
        }
      }
    } else if (exosphereType == "Chamberlain") {
      for (int i = 0; i < nExoComponent; ++i) {
        if (get_rPlanet_SI() > 0.0 && r > 0.0) {
          sum += exoN0[i] * exp(-exoH0[i] * (1.0 / get_rPlanet_SI() - 1.0 / r));
        }
      }
    }
    return sum;
  }

  amrex::Real get_exosphere_component_density(amrex::Real r,
                                              int iC) const override {
    if (exosphereType == "None") return 0.0;
    if (iC < 0 || iC >= nExoComponent) return 0.0;
    if (r < get_rPlanet_SI()) return 0.0;

    if (exosphereType == "Exponential") {
      if (exoH0[iC] > 0.0) {
        return exoN0[iC] * exp(-(r - get_rPlanet_SI()) / exoH0[iC]);
      }
    } else if (exosphereType == "Power-Law") {
      if (r > 0.0) {
        return exoN0[iC] * pow(get_rPlanet_SI() / r, exoK0[iC]);
      }
    } else if (exosphereType == "Chamberlain") {
      if (get_rPlanet_SI() > 0.0 && r > 0.0) {
        return exoN0[iC] * exp(-exoH0[iC] * (1.0 / get_rPlanet_SI() - 1.0 / r));
      }
    }
    return 0.0;
  }

  // Check whether (x, y, z) [m] relative to planet center is in shadow.
  bool is_in_shadow(amrex::Real x, amrex::Real y, amrex::Real z) const {
    if (!useShadowCylinder) return false;
    amrex::Real proj = x * solarDir[0] + y * solarDir[1] + z * solarDir[2];
    if (proj >= 0.0) return false;
    if (-proj > shadowCylinderHalfHeight)
      return false;
    amrex::Real r2 = x * x + y * y + z * z;
    amrex::Real perp2 = r2 - proj * proj;
    return perp2 <= shadowCylinderRadius * shadowCylinderRadius;
  }

  // Electron temperature in eV from PIC-normalized pressure and density.
  amrex::Real electron_temperature(amrex::Real pe, amrex::Real ne) const {
    if (ne <= 0.0) return 0.0;
    const amrex::Real protonMassPerCharge =
        cProtonMassSI / cUnitChargeSI;   // [kg/C]
    const amrex::Real ur2 =
        get_unorm_si() * get_unorm_si(); // [m^2/s^2]
    return protonMassPerCharge * ur2 * pe / ne;
  }

  // Photoionization frequency [s^-1] for neutral component iC.
  amrex::Real photoionization_rate(const amrex::Real xyz[3], int iC,
                                   amrex::Real photoDilution = -1.0) const {
    if (photoDilution >= 0.0) {
      return photoNu0[iC] * photoDilution;
    }
    if (is_in_shadow(xyz[0], xyz[1], xyz[2])) return 0.0;
    amrex::Real r2 = xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2];
    amrex::Real r_m = sqrt(r2);
    if (r_m <= 0.0) return 0.0;
    amrex::Real ratio = get_rPlanet_SI() / r_m;
    return photoNu0[iC] * ratio * ratio;
  }

  // Electron-impact ionization frequency [s^-1] for neutral component iC.
  amrex::Real impact_ionization_rate(amrex::Real ne, amrex::Real Te_eV,
                                     int iC) const {
    if (ne <= 0.0 || Te_eV <= 0.0) return 0.0;
    amrex::Real ne_SI = ne / (get_Si2NoRho() * cProtonMassSI);
    amrex::Real Ei = impactEIon[iC];
    amrex::Real A = impactA[iC];
    amrex::Real K = impactK[iC];
    amrex::Real X = impactX[iC];
    amrex::Real t = Te_eV / Ei;
    amrex::Real sigma_v = A * pow(t, K) / (X + t) * exp(-Ei / Te_eV) * 1e-6;
    return ne_SI * sigma_v;
  }

  // Charge-exchange frequency [s^-1] for neutral component iC with ion iSp.
  amrex::Real charge_exchange_rate(const FluidInterface& other,
                                   const amrex::MFIter& mfi,
                                   const amrex::IntVect& idx, int iLev,
                                   int iSp, int iC) const {
    int iIon = iSp - 1;
    amrex::Real ni = other.get_number_density(mfi, idx, iSp, iLev);
    amrex::Real ux_i = other.get_ux(mfi, idx, iSp, iLev);
    amrex::Real uy_i = other.get_uy(mfi, idx, iSp, iLev);
    amrex::Real uz_i = other.get_uz(mfi, idx, iSp, iLev);
    amrex::Real u_mag_SI =
        sqrt(ux_i * ux_i + uy_i * uy_i + uz_i * uz_i) * get_unorm_si();
    if (ni <= 0.0 || u_mag_SI <= 0.0) return 0.0;
    amrex::Real ni_SI = ni / (get_Si2NoRho() * cProtonMassSI);
    amrex::Real sigma = cxSigma[iC * nCXIonSpecies + iIon];
    return ni_SI * sigma * 1e-4 * u_mag_SI;
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
      amrex::Real norm = sqrt(solarDir[0] * solarDir[0] +
                              solarDir[1] * solarDir[1] +
                              solarDir[2] * solarDir[2]);
      if (norm > 0.0) {
        solarDir[0] /= norm;
        solarDir[1] /= norm;
        solarDir[2] /= norm;
      }
    } else if (command == "#RECOMBINATION") {
      useRecombination = true;
      int nRecomb;
      param.read_var("nReactions", nRecomb);
      recombIonIndex.resize(nRecomb);
      recombRate0.resize(nRecomb);
      recombTempExp.resize(nRecomb);
      recombRefTemp.resize(nRecomb);
      for (int i = 0; i < nRecomb; ++i) {
        param.read_var("ionSpecies", recombIonIndex[i]);
        param.read_var("rateCoef", recombRate0[i]);
        param.read_var("tempExponent", recombTempExp[i]);
        param.read_var("refTemp", recombRefTemp[i]);
      }
    } else if (command == "#CHEMISTRY") {
      useChemistry = true;
      int nRxns;
      param.read_var("nReactions", nRxns);
      chemReactions.resize(nRxns);
      for (int i = 0; i < nRxns; ++i) {
        param.read_var("reactantIon", chemReactions[i].reactantIon);
        param.read_var("productIon", chemReactions[i].productIon);
        param.read_var("neutralComp", chemReactions[i].neutralComp);
        param.read_var("rateType", chemReactions[i].rateType);
        param.read_var("rateCoef", chemReactions[i].rateCoef);
        param.read_var("tempExp", chemReactions[i].tempExp);
        param.read_var("refTemp", chemReactions[i].refTemp);
      }
    }
  }

  // Validate consistency between exosphere and ionization commands.
  void post_process_param() override {
    // Recombination validation.
    if (useRecombination) {
      if (nS < 1) {
        amrex::Abort(printPrefix + "Error: #RECOMBINATION requires "
                     + "plasma species. Use #PLASMA to set species.");
      }
      if (!useElectronFluid) {
        amrex::Abort(printPrefix + "Error: #RECOMBINATION requires "
                     + "useElectronFluid = true (species 0 must be the "
                     + "electron). Set #PLASMA with an electron species.");
      }
      const int nIonS = nS - 1;
      for (int i = 0; i < static_cast<int>(recombIonIndex.size()); ++i) {
        int iSp = recombIonIndex[i];
        if (iSp < 1 || iSp > nIonS) {
          amrex::Abort(printPrefix + "Error: #RECOMBINATION ionSpecies "
                       + std::to_string(iSp) + " is out of range [1, "
                       + std::to_string(nIonS) + "].");
        }
        if (recombRate0[i] <= 0.0) {
          amrex::Abort(printPrefix + "Error: #RECOMBINATION rateCoef must "
                       + "be positive (got "
                       + std::to_string(recombRate0[i]) + ").");
        }
      }
    }

    // Chemistry validation.
    if (useChemistry) {
      if (usePhotoIonization || useChargeExchange || useRecombination) {
        amrex::Abort(printPrefix + "Error: #CHEMISTRY cannot be combined with "
                     + "#PHOTOIONIZATION, #CHARGEEXCHANGE, or #RECOMBINATION. "
                     + "Use either the unified #CHEMISTRY table or individual commands.");
      }
      if (nS < 2) {
        amrex::Abort(printPrefix + "Error: #CHEMISTRY requires at least "
                     + "2 plasma species (electron + 1 ion).");
      }
      const int nIonS = nS - 1;
      bool needsNeutral = false;
      bool needsElectron = false;
      for (int i = 0; i < static_cast<int>(chemReactions.size()); ++i) {
        const auto& rxn = chemReactions[i];
        if (rxn.reactantIon < 0 || rxn.reactantIon > nIonS) {
          amrex::Abort(printPrefix + "Error: #CHEMISTRY reactantIon "
                       + std::to_string(rxn.reactantIon) + " out of range "
                       + "[0, " + std::to_string(nIonS) + "].");
        }
        if (rxn.productIon < 0 || rxn.productIon > nIonS) {
          amrex::Abort(printPrefix + "Error: #CHEMISTRY productIon "
                       + std::to_string(rxn.productIon) + " out of range "
                       + "[0, " + std::to_string(nIonS) + "].");
        }
        if (rxn.reactantIon == 0 && rxn.productIon == 0) {
          amrex::Abort(printPrefix + "Error: #CHEMISTRY reaction "
                       + std::to_string(i) + " has both reactantIon and "
                       + "productIon = 0 (no-op).");
        }
        if (rxn.rateType == 1 && rxn.neutralComp < 0) {
          amrex::Abort(printPrefix + "Error: #CHEMISTRY photoionization "
                       + "reaction " + std::to_string(i)
                       + " requires a neutral component.");
        }
        if (rxn.rateType == 0 && rxn.neutralComp < 0 &&
            rxn.reactantIon == 0) {
          amrex::Abort(printPrefix + "Error: #CHEMISTRY thermal reaction "
                       + std::to_string(i) + " with no neutral and no "
                       + "reactant ion is invalid.");
        }
        if (rxn.neutralComp >= 0) {
          needsNeutral = true;
          if (nExoComponent > 0 && rxn.neutralComp >= nExoComponent) {
            amrex::Abort(printPrefix + "Error: #CHEMISTRY reaction "
                         + std::to_string(i) + " neutralComp "
                         + std::to_string(rxn.neutralComp)
                         + " >= nExoComponent "
                         + std::to_string(nExoComponent) + ".");
          }
        }
        if (rxn.productIon == 0 || rxn.tempExp != 0.0) {
          needsElectron = true;
        }
      }
      if (needsNeutral && nExoComponent <= 0) {
        amrex::Abort(printPrefix + "Error: #CHEMISTRY reactions with neutralComp >= 0 "
                     + "require #EXOSPHERE with nComponent > 0 to be specified.");
      }
      if (needsElectron && !useElectronFluid) {
        amrex::Abort(printPrefix + "Error: #CHEMISTRY reactions with recombination "
                     + "or temperature dependence require useElectronFluid = true. "
                     + "Set #PLASMA with an electron species.");
      }
    }

    if (exosphereType == "None") return;

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
      if (nCXIonSpecies != nS - 1) {
        amrex::Abort(printPrefix + "Error: #CHARGEEXCHANGE nIonSpecies ("
                     + std::to_string(nCXIonSpecies)
                     + ") != number of ion species ("
                     + std::to_string(nS - 1) + ")");
      }
    }

    for (int iC = 0; iC < nExoComponent; ++iC) {
      if (exoT0[iC] <= 0.0) {
        amrex::Abort(printPrefix + "Error: #EXOSPHERE T0 for component "
                     + std::to_string(iC) + " must be positive (got "
                     + std::to_string(exoT0[iC])
                     + " K). It is used as the ionization source temperature.");
      }
    }
  }

  // Reaction frequency [s^-1] for a chemistry reaction.
  amrex::Real chem_reaction_rate(const ChemistryReaction& rxn,
                                 const amrex::Real xyz[3], amrex::Real r_val,
                                 amrex::Real ne, amrex::Real Te_eV,
                                 amrex::Real photoDilution = -1.0) const {
    if (rxn.rateType == 1) {
      if (photoDilution >= 0.0) {
        return rxn.rateCoef * photoDilution;
      }
      if (is_in_shadow(xyz[0], xyz[1], xyz[2])) return 0.0;
      if (r_val <= 0.0) return 0.0;
      amrex::Real ratio = get_rPlanet_SI() / r_val;
      return rxn.rateCoef * ratio * ratio;
    }

    // Thermal rate coefficient k(Te) [cm^3/s] -> [m^3/s]
    amrex::Real k_si = rxn.rateCoef * 1e-6;
    if (rxn.tempExp != 0.0 && Te_eV > 0.0) {
      amrex::Real Te_K = Te_eV * cUnitChargeSI / cBoltzmannSI;
      k_si *= pow(rxn.refTemp / Te_K, rxn.tempExp);
    }

    if (rxn.neutralComp >= 0) {
      amrex::Real n_neutral =
          get_exosphere_component_density(r_val, rxn.neutralComp);
      return k_si * n_neutral;
    }

    return k_si * ne / (get_Si2NoRho() * cProtonMassSI);
  }

  //-------------------------------------------------------------------
  // Apply chemistry reaction source term to accumulators.
  void chem_apply_source(const ChemistryReaction& rxn, amrex::Real rate,
                         const FluidInterface& other,
                         const amrex::MFIter& mfi,
                         const amrex::IntVect& idx, int iLev,
                         int nIonS, amrex::Real r_val,
                         std::vector<amrex::Real>& srcRho,
                         std::vector<amrex::Real>& srcP,
                         std::vector<amrex::Real>& srcRhoUx,
                         std::vector<amrex::Real>& srcRhoUy,
                         std::vector<amrex::Real>& srcRhoUz) const {
    if (rxn.productIon <= 0 || rxn.productIon > nIonS) return;

    int iSpProd = rxn.productIon;
    amrex::Real mass_prod = get_species_mass(iSpProd);

    if (rxn.reactantIon > 0 && rxn.reactantIon <= nIonS) {
      // Cross-species CX: product inherits reactant velocity and temperature.
      int iSpReac = rxn.reactantIon;
      amrex::Real rho_reac_norm =
          other.get_value(mfi, idx, iRho_I[iSpReac], iLev);
      if (rho_reac_norm <= 0.0) return;

      amrex::Real mass_reac = get_species_mass(iSpReac);
      amrex::Real n_reac_norm = rho_reac_norm / mass_reac;
      amrex::Real n_reac_si =
          n_reac_norm / (get_Si2NoRho() * cProtonMassSI);

      amrex::Real srcRho_si =
          rate * n_reac_si * mass_prod * cProtonMassSI;
      srcRho[iSpProd] += srcRho_si;

      amrex::Real ux_reac = other.get_ux(mfi, idx, iSpReac, iLev);
      amrex::Real uy_reac = other.get_uy(mfi, idx, iSpReac, iLev);
      amrex::Real uz_reac = other.get_uz(mfi, idx, iSpReac, iLev);
      srcRhoUx[iSpProd] += srcRho_si * ux_reac;
      srcRhoUy[iSpProd] += srcRho_si * uy_reac;
      srcRhoUz[iSpProd] += srcRho_si * uz_reac;

      amrex::Real p_reac_norm =
          other.get_value(mfi, idx, iP_I[iSpReac], iLev);
      amrex::Real p_reac_si = p_reac_norm / get_Si2NoP();
      srcP[iSpProd] += rate * p_reac_si;
    } else {
      // Photoionization: source from neutral at rest (zero velocity).
      amrex::Real n_neutral =
          get_exosphere_component_density(r_val, rxn.neutralComp);
      amrex::Real S_n = n_neutral * rate;
      srcRho[iSpProd] += S_n * mass_prod * cProtonMassSI;
      if (rxn.neutralComp < nExoComponent) {
        srcP[iSpProd] += S_n * cBoltzmannSI * exoT0[rxn.neutralComp];
      }

      // Co-create neutralizing electrons in full-PIC mode.
      if (nS > 0 && get_species_charge(0) < 0) {
        amrex::Real mass_e = get_species_mass(0);
        srcRho[0] += S_n * mass_e * cProtonMassSI;
        if (rxn.neutralComp < nExoComponent) {
          srcP[0] += S_n * cBoltzmannSI * exoT0[rxn.neutralComp];
        }
      }
    }
  }

  // Apply chemistry reaction loss term to loss array.
  void chem_apply_loss(const ChemistryReaction& rxn, amrex::Real rate,
                       const FluidInterface& other,
                       const amrex::MFIter& mfi,
                       const amrex::IntVect& idx, int iLev,
                       int nIonS, int i, int j, int k,
                       amrex::Array4<amrex::Real>& lossArr) const {
    if (rxn.reactantIon <= 0 || rxn.reactantIon > nIonS) return;

    int iSpReac = rxn.reactantIon;
    amrex::Real rho_reac_norm =
        other.get_value(mfi, idx, iRho_I[iSpReac], iLev);
    if (rho_reac_norm <= 0.0) return;

    amrex::Real lossRho_norm = rate * rho_reac_norm / get_Si2NoT();
    if (lossRho_norm > 0.0) {
      lossArr(i, j, k, iSpReac) += lossRho_norm;

      if (rxn.productIon == 0 && nS > 0 && get_species_charge(0) < 0) {
        amrex::Real mass_reac = get_species_mass(iSpReac);
        amrex::Real mass_e = get_species_mass(0);
        lossArr(i, j, k, 0) += lossRho_norm * (mass_e / mass_reac);
      }
    }
  }

  // Apply #RECOMBINATION loss terms.
  void apply_recombination_loss(const FluidInterface& other,
                                const amrex::MFIter& mfi,
                                const amrex::IntVect& idx, int iLev,
                                int nIonS, amrex::Real ne, amrex::Real Te_eV,
                                int i, int j, int k,
                                amrex::Array4<amrex::Real>& lossArr) const {
    for (int iR = 0; iR < static_cast<int>(recombIonIndex.size()); ++iR) {
      int iSp = recombIonIndex[iR];
      if (iSp < 1 || iSp > nIonS) continue;

      amrex::Real k_si = recombRate0[iR] * 1e-6;
      if (Te_eV > 0.0 && recombTempExp[iR] != 0.0) {
        amrex::Real Te_K = Te_eV * cUnitChargeSI / cBoltzmannSI;
        k_si *= pow(recombRefTemp[iR] / Te_K, recombTempExp[iR]);
      }

      amrex::Real rho_ion_norm =
          other.get_value(mfi, idx, iRho_I[iSp], iLev);
      if (rho_ion_norm <= 0.0) continue;

      amrex::Real lossRho_norm = k_si * ne * rho_ion_norm /
          (get_Si2NoRho() * cProtonMassSI * get_Si2NoT());
      if (lossRho_norm > 0.0) {
        lossArr(i, j, k, iSp) += lossRho_norm;

        if (nS > 0 && get_species_charge(0) < 0) {
          amrex::Real mass_ion = get_species_mass(iSp);
          amrex::Real mass_e = get_species_mass(0);
          lossArr(i, j, k, 0) += lossRho_norm * (mass_e / mass_ion);
        }
      }
    }
  }

  // Set nodeFluid from plasma-state-dependent ionization processes.
  void set_source(const FluidInterface& other) override {
#ifdef AMREX_USE_GPU
    const_cast<FluidInterface&>(other).sync_host_fluid();
#endif
    set_node_fluid(other);
    set_node_loss_fluid_to_zero();

    const amrex::Box gbx = convert(Geom(0).Domain(), { AMREX_D_DECL(1, 1, 1) });

    const bool doPhoto = usePhotoIonization;
    const bool doImpact = useElectronImpact;
    const bool doCX = useChargeExchange;
    const bool doRecomb = useRecombination;
    const bool doChem = useChemistry;

    const amrex::Real rhoNormPerT = get_Si2NoRho() / get_Si2NoT();
    const amrex::Real pNormPerT = get_Si2NoP() / get_Si2NoT();
    const amrex::Real rPlanet = get_rPlanet_SI();
    const int nIonS = nS - 1;

    for (int iLev = 0; iLev < n_lev(); iLev++) {
#ifdef AMREX_USE_GPU
      auto& targetFluid = h_nodeFluid;
      auto& targetLossFluid = h_nodeLossFluid;
#else
      auto& targetFluid = nodeFluid;
      auto& targetLossFluid = nodeLossFluid;
#endif
      if (!targetFluid[iLev].empty()) {
        for (amrex::MFIter mfi(targetFluid[iLev]); mfi.isValid(); ++mfi) {
          const amrex::Real* dx = Geom(iLev).CellSize();
          const auto plo = Geom(iLev).ProbLo();

          const amrex::Box& box = mfi.validbox();
          const auto lo = lbound(box);
          const auto hi = ubound(box);

          const amrex::Array4<amrex::Real>& arr = targetFluid[iLev][mfi].array();

          amrex::Array4<amrex::Real> lossArr;
          if (doRecomb || doChem) {
            lossArr = targetLossFluid[iLev][mfi].array();
          }

          // Source accumulators in SI units.
          std::vector<amrex::Real> srcRho(nIonS + 1, 0.0);
          std::vector<amrex::Real> srcP(nIonS + 1, 0.0);
          std::vector<amrex::Real> srcRhoUx(nIonS + 1, 0.0);
          std::vector<amrex::Real> srcRhoUy(nIonS + 1, 0.0);
          std::vector<amrex::Real> srcRhoUz(nIonS + 1, 0.0);

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
                amrex::Real xyz[3] = { 0, 0, 0 };
                for (int iDim = 0; iDim < get_fluid_dimension(); iDim++) {
                  xyz[iDim] = idx[iDim] * dx[iDim] + plo[iDim];
                }

                amrex::Real r_val = 0.0;
                for (int d = 0; d < 3; ++d) {
                  r_val += xyz[d] * xyz[d];
                }
                r_val = sqrt(r_val);

                const bool inShadow = is_in_shadow(xyz[0], xyz[1], xyz[2]);
                const amrex::Real photoDilution =
                    (inShadow || r_val <= 0.0) ? 0.0 : (rPlanet / r_val) * (rPlanet / r_val);

                amrex::Real ne = 0, pe = 0, Te_eV = 0;
                bool plasmaFetched = false;
                auto fetch_electron_plasma = [&, &other=other]() {
                  if (plasmaFetched) return;
                  ne = other.get_number_density(mfi, idx, 0, iLev);
                  pe = other.get_p(mfi, idx, 0, iLev);
                  Te_eV = electron_temperature(pe, ne);
                  plasmaFetched = true;
                };

                std::fill(srcRho.begin(), srcRho.end(), 0.0);
                std::fill(srcP.begin(), srcP.end(), 0.0);
                std::fill(srcRhoUx.begin(), srcRhoUx.end(), 0.0);
                std::fill(srcRhoUy.begin(), srcRhoUy.end(), 0.0);
                std::fill(srcRhoUz.begin(), srcRhoUz.end(), 0.0);

                // Exosphere-based ionization sources.
                for (int iC = 0; iC < nExoComponent; ++iC) {
                  const int iSp = iC + 1;
                  if (iSp > nIonS) break;

                  amrex::Real dens_i =
                      get_exosphere_component_density(r_val, iC);

                  amrex::Real nu_tot = 0.0;
                  if (doPhoto)
                    nu_tot += photoNu0[iC] * photoDilution;
                  if (doImpact) {
                    fetch_electron_plasma();
                    nu_tot += impact_ionization_rate(ne, Te_eV, iC);
                  }
                  if (doCX) {
                    for (int iSpCX = 1; iSpCX <= nIonS; ++iSpCX) {
                      nu_tot += charge_exchange_rate(other, mfi, idx, iLev,
                                                     iSpCX, iC);
                    }
                  }

                  amrex::Real S_n = dens_i * nu_tot;
                  amrex::Real mass_amu = get_species_mass(iSp);
                  srcRho[iSp] += S_n * mass_amu * cProtonMassSI;
                  srcP[iSp] += S_n * cBoltzmannSI * exoT0[iC];

                  // Co-create neutralizing electrons in full-PIC mode.
                  if (nS > 0 && get_species_charge(0) < 0) {
                    amrex::Real mass_e = get_species_mass(0);
                    srcRho[0] += S_n * mass_e * cProtonMassSI;
                    srcP[0] += S_n * cBoltzmannSI * exoT0[iC];
                  }
                }

                // Chemistry sources.
                if (doChem) {
                  fetch_electron_plasma();
                  for (int iR = 0;
                       iR < static_cast<int>(chemReactions.size()); ++iR) {
                    const auto& rxn = chemReactions[iR];
                    amrex::Real rate = chem_reaction_rate(rxn, xyz, r_val,
                                                          ne, Te_eV, photoDilution);
                    if (rate <= 0.0) continue;
                    chem_apply_source(rxn, rate, other, mfi, idx, iLev,
                                      nIonS, r_val, srcRho, srcP,
                                      srcRhoUx, srcRhoUy, srcRhoUz);
                  }
                }

                // Write accumulated sources to nodeFluid.
                for (int iFluid = 0; iFluid < nFluid; iFluid++) {
                  arr(i, j, k, iRho_I[iFluid]) = 0;
                  arr(i, j, k, iUx_I[iFluid]) = 0;
                  arr(i, j, k, iUy_I[iFluid]) = 0;
                  arr(i, j, k, iUz_I[iFluid]) = 0;
                  arr(i, j, k, iP_I[iFluid]) = 0;
                }

                bool anySource = false;
                if (nS > 0 && iRho_I[0] >= 0 && srcRho[0] > 0) {
                  arr(i, j, k, iRho_I[0]) = srcRho[0] * rhoNormPerT;
                  arr(i, j, k, iP_I[0]) = srcP[0] * pNormPerT;
                  arr(i, j, k, iUx_I[0]) = 0.0;
                  arr(i, j, k, iUy_I[0]) = 0.0;
                  arr(i, j, k, iUz_I[0]) = 0.0;
                }
                for (int iSp = 1; iSp <= nIonS; ++iSp) {
                  if (srcRho[iSp] > 0) {
                    anySource = true;
                    amrex::Real rho_norm = srcRho[iSp] * rhoNormPerT;
                    arr(i, j, k, iRho_I[iSp]) = rho_norm;
                    arr(i, j, k, iUx_I[iSp]) = srcRhoUx[iSp] * rhoNormPerT;
                    arr(i, j, k, iUy_I[iSp]) = srcRhoUy[iSp] * rhoNormPerT;
                    arr(i, j, k, iUz_I[iSp]) = srcRhoUz[iSp] * rhoNormPerT;
                    arr(i, j, k, iP_I[iSp]) = srcP[iSp] * pNormPerT;
                  }
                }
                if (anySource && iPe >= 0) {
                  amrex::Real srcPe = 0.0;
                  for (int iSp = 1; iSp <= nIonS; ++iSp)
                    srcPe += srcP[iSp];
                  arr(i, j, k, iPe) = srcPe * pNormPerT;
                }

                // Loss terms.
                if (doRecomb) {
                  fetch_electron_plasma();
                  apply_recombination_loss(other, mfi, idx, iLev, nIonS,
                                           ne, Te_eV, i, j, k, lossArr);
                }

                if (doChem) {
                  fetch_electron_plasma();
                  for (int iR = 0;
                       iR < static_cast<int>(chemReactions.size()); ++iR) {
                    const auto& rxn = chemReactions[iR];
                    amrex::Real rate = chem_reaction_rate(rxn, xyz, r_val,
                                                          ne, Te_eV, photoDilution);
                    if (rate <= 0.0) continue;
                    chem_apply_loss(rxn, rate, other, mfi, idx, iLev,
                                    nIonS, i, j, k, lossArr);
                  }
                }
              } // for k
        }
      }
    }

#ifdef AMREX_USE_GPU
    for (int iLev = 0; iLev < n_lev(); iLev++) {
      if (!nodeFluid[iLev].empty()) {
        nodeFluid[iLev].ParallelCopy(h_nodeFluid[iLev]);
      }
      if ((doRecomb || doChem) && !nodeLossFluid[iLev].empty()) {
        nodeLossFluid[iLev].ParallelCopy(h_nodeLossFluid[iLev]);
      }
    }
    amrex::Gpu::streamSynchronize();
#endif

    if (!isGridEmpty && useCurrent) {
      for (int iLev = 0; iLev < n_lev(); iLev++) {
        amrex::MultiFab currentMF(nodeFluid[iLev], amrex::make_alias, iJx, 3);
        currentMF.setVal(0, currentMF.nGrow());
      }
    }
  }
};

#endif
