#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <AMReX_MultiFabUtil.H>

#include "GridUtility.h"
#include "Pic.h"
#include "Timer.h"

using namespace amrex;

//==========================================================
void Pic::read_param(const std::string& command, ReadParam& param) {

  if (command == "#PIC") {
    param.read_var("usePIC", usePIC);
  } else if (command == "#SOLVEEM") {
    param.read_var("solveEM", solveEM);
  } else if (command == "#PARTMODE") {
    std::string s;
    param.read_var("partMode", s);
    if (s == "SEP")
      pMode = PartMode::SEP;
    else if (s == "PIC")
      pMode = PartMode::PIC;
    else
      Abort("Error: wrong input for partMode.");
  } else if (command == "#PARTICLEBOXBOUNDARY") {
    int iSpecies;
    std::string lo, hi;
    param.read_var("iSpecies", iSpecies);
    if (iSpecies < 0)
      Abort("Error: negative species index in #PARTICLEBOXBOUNDARY.");
    // nSpecies is only known in post_process_param(), so pBCs grows on
    // demand here and is padded to nSpecies there.
    if (iSpecies >= static_cast<int>(pInfo.pBCs.size())) {
      pInfo.pBCs.resize(iSpecies + 1);
      pInfo.pBCsSet.resize(iSpecies + 1, 0);
    }
    pInfo.pBCsSet[iSpecies] = 1;

    for (int i = 0; i < nDim; ++i) {
      param.read_var("particleBoxBoundaryLo", lo);
      param.read_var("particleBoxBoundaryHi", hi);
      pInfo.pBCs[iSpecies].set(i, 0, ParticleBC::parse(lo));
      pInfo.pBCs[iSpecies].set(i, 1, ParticleBC::parse(hi));
    }
  } else if (command == "#FIELDBOXBOUNDARY" ||
             command == "#BFIELDBOXBOUNDARY") {
    if (command == "#BFIELDBOXBOUNDARY")
      add_bc_warning("#BFIELDBOXBOUNDARY is deprecated; use "
                     "#FIELDBOXBOUNDARY instead.");
    fieldBCSet_ = true;

    std::string lo, hi;
    for (int i = 0; i < nDim; ++i) {
      param.read_var("fieldBoxBoundaryLo", lo);
      param.read_var("fieldBoxBoundaryHi", hi);
      bcField.set(i, 0, FieldBC::parse(lo));
      bcField.set(i, 1, FieldBC::parse(hi));
    }
  } else if (command == "#ABSORB") {
    param.read_var("charSpeed", absorbCharSpeed);
  } else if (command == "#INFLOW") {
    double tmp;
    param.read_var("rho", tmp);
    inflowRho_ = tmp; // [amu/cc]
    param.read_var("ux", tmp);
    inflowUx_ = tmp; // [km/s]
    param.read_var("uy", tmp);
    inflowUy_ = tmp; // [km/s]
    param.read_var("uz", tmp);
    inflowUz_ = tmp; // [km/s]
    param.read_var("T", tmp);
    inflowT_ = tmp; // [K]
    inflowDefined_ = true;
  } else if (command == "#WAVEBC") {
    waveBC.read_param(param, fi);
  } else if (command == "#MEMORY") {
    param.read_var("dnMemory", dnMemory);
  } else if (command == "#RANDOMPARTICLESLOCATION") {
    param.read_var("isParticleLocationRandom", pInfo.isParticleLocationRandom);
  } else if (command == "#CONSTANTPPV") {
    param.read_var("isPPVconstant", pInfo.isPPVconstant);
  } else if (command == "#PRESPLITTING") {
    param.read_var("doPreSplitting", pInfo.doPreSplitting);
  } else if (command == "#DIVE") {
    param.read_var("doCorrectDivE", doCorrectDivE);
    if (doCorrectDivE) {
      param.read_var("nDivECorrection", nDivECorrection);
    }
  } else if (command == "#EXPLICITPIC") {
    param.read_var("useExplicitPIC", useExplicitPIC);
  } else if (command == "#EFIELDSOLVER") {
    Real tol;
    int nIter;
    param.read_var("tol", tol);
    param.read_var("nIter", nIter);
    eSolver.set_tol(tol);
    eSolver.set_nIter(nIter);
  } else if (command == "#PARTICLES") {
    param.read_var("npcelx", pInfo.nPartPerCell[ix_]);
    param.read_var("npcely", pInfo.nPartPerCell[iy_]);
    if (nDim == 3)
      param.read_var("npcelz", pInfo.nPartPerCell[iz_]);
  } else if (command == "#SOURCEPARTICLES") {
    param.read_var("npcelx", nSourcePPC[ix_]);
    param.read_var("npcely", nSourcePPC[iy_]);
    if (nDim == 3)
      param.read_var("npcelz", nSourcePPC[iz_]);
  } else if (command == "#KINETICSOURCE") {
    param.read_var("kineticSource", kineticSource);
  } else if (command == "#ELECTRON") {
    param.read_var("qom", qomEl);
  } else if (command == "#DISCRETIZE" || command == "#DISCRETIZATION") {
    param.read_var("theta", fsolver.theta);
    param.read_var("coefDiff", fsolver.coefDiff);
  } else if (command == "#COMOVING") {
    param.read_var("solveFieldInCoMov", solveFieldInCoMov);
    param.read_var("nSmoothBackGroundU", nSmoothBackGroundU);
  } else if (command == "#UPWINDE") {
    param.read_var("useUpwindE", useUpwindE);
    param.read_var("limiterThetaE", limiterThetaE);
  } else if (command == "#LAGGEDLIMITER") {
    param.read_var("useLaggedLimiter", fsolver.useLaggedLimiter);
  } else if (command == "#CMAXE") {
    param.read_var("cMaxE", cMaxE);
  } else if (command == "#SMOOTHE") {
    param.read_var("doSmoothE", doSmoothE);
    if (doSmoothE) {
      param.read_var("nSmoothE", nSmoothE);
    }
  } else if (command == "#SMOOTHJ") {
    param.read_var("doSmoothJ", doSmoothJ);
    if (doSmoothJ) {
      param.read_var("nSmoothJ", nSmoothJ);
      param.read_var("coefSmoothJ", coefSmoothJ);
    }
  } else if (command == "#SMOOTHMOMENTS") {
    param.read_var("doSmoothMoments", doSmoothMoments);
    if (doSmoothMoments) {
      param.read_var("nSmoothMoments", nSmoothMoments);
      param.read_var("coefSmoothMoments", coefSmoothMoments);
    }
  } else if (command == "#UPWINDB") {
    param.read_var("useUpwindB", useUpwindB);
    param.read_var("theta", limiterThetaB);
    if (useUpwindB) {
      useHyperbolicCleaning = true;
    }
    param.read_optional("fixedUpwindVel", fixedUpwindVel);
  } else if (command == "#FIXEDUMAX") {
    param.read_optional("fixedUMax", fixedUMax);
  } else if (command == "#DIVB") {
    param.read_var("useHyperbolicCleaning", useHyperbolicCleaning);
    if (useHyperbolicCleaning) {
      param.read_var("hypDecay", hypDecay);
    }
  } else if (command == "#RESAMPLING") {
    param.read_var("doReSampling", doReSampling);
    if (doReSampling) {
      param.read_var("reSamplingLowLimit", reSamplingLowLimit);
      param.read_var("reSamplingHighLimit", reSamplingHighLimit);
      param.read_var("maxWeightRatio", maxWeightRatio);
    }
  } else if (command == "#FASTMERGE") {
    param.read_var("fastMerge", pInfo.fastMerge);
    if (pInfo.fastMerge) {
      param.read_var("nMergeOld", pInfo.nPartCombine);
      param.read_var("nMergeNew", pInfo.nPartNew);
      param.read_var("nMergeTry", pInfo.nMergeTry);
      param.read_var("mergeRatioMax", pInfo.mergeRatioMax);
    }
  } else if (command == "#ADAPTIVESOURCEPPC") {
    param.read_var("adaptiveSourcePPC", adaptiveSourcePPC);
  } else if (command == "#MERGELIGHT") {
    param.read_var("mergeLight", pInfo.mergeLight);
    if (pInfo.mergeLight) {
      param.read_var("mergePartRatioMax", pInfo.mergePartRatioMax);
    }
  } else if (command == "#VACUUM") {
    param.read_var("vacuum", pInfo.vacuumIO);
  } else if (command == "#PARTICLELEVRATIO") {
    param.read_var("particleLevRatio", pInfo.pLevRatio);
  } else if (command == "#OHION") {
    param.read_var("rAnalytic", pInfo.ionOH.rAnalytic);
    param.read_var("doGetFromOH", pInfo.ionOH.doGetFromOH);

    if (!pInfo.ionOH.doGetFromOH) {
      param.read_var("rCutoff", pInfo.ionOH.rCutoff);
      param.read_var("swRho", pInfo.ionOH.swRho);
      param.read_var("swT", pInfo.ionOH.swT);
      param.read_var("swU", pInfo.ionOH.swU);
    }
  } else if (command == "#SUPID") {
    int n = 0;
    param.read_var("nSpecies", n);
    for (int i = 0; i < n; ++i) {
      int supid;
      param.read_var("supid", supid);
      pInfo.supIDs.push_back(supid);
    }
  } else if (command == "#MAXCHARGEEXCHANGERATE") {
    param.read_var("maxChargeExchangeRate", maxExchangeRatioLimit);
  } else if (command == "#TESTCASE") {
    std::string testcase;
    param.read_var("testCase", testcase);

    ic_ = ICRegistry::instance().create(testcase);
    if (!ic_) {
      std::string known;
      for (const auto& n : ICRegistry::instance().names()) {
        if (!known.empty())
          known += ", ";
        known += n;
      }
      amrex::Abort("Unknown #TESTCASE name '" + testcase +
                   "'. Registered names: " + known + ".");
    }
    ic_->read_param(param);
  } else if (command == "#WAVEIC") {
    if (!ic_) {
      amrex::Abort("The #WAVEIC block must follow a #TESTCASE that selects a "
                   "wave initial condition (waveic / lightwave / hybridwave / "
                   "convectionwave / ionacousticwave).");
    }
    ic_->read_param(param);
  } else if (command == "#FADEEVIC") {
    if (!ic_ || std::string(ic_->name()) != "fadeev") {
      amrex::Abort("The #FADEEVIC block must follow a #TESTCASE that selects "
                   "the fadeev (magnetic reconnection) initial condition.");
    }
    ic_->read_param(param);
  } else if (command == "#HYBRIDPIC") {
    param.read_var("useHybridPIC", useHybridPIC);
  } else if (command == "#RESISTIVITY") {
    param.read_var("etaResistivity", etaResistivitySI);
  } else if (command == "#ELECTRONTEMPERATURE") {
    param.read_var("electronTemperature", electronTemperatureEV);
    param.read_var("electronGamma", electronGamma);
    param.read_var("electronDensity0", electronDensity0In);
  } else if (command == "#BSUBCYCLE") {
    param.read_var("nBSubcycle", nBSubcycle);
  } else if (command == "#HALLTERM") {
    param.read_var("useHallTerm", useHallTerm);
  } else if (command == "#HYPERRESISTIVITY") {
    param.read_var("etaHyperSI", etaHyperSI);
    param.read_var("etaHyperMode", etaHyperMode);
    param.read_var("etaHyperCh", etaHyperCh);
  } else if (command == "#MINIMUMDENSITY") {
    param.read_var("rhoMinOhm", rhoMinOhm);
  } else if (command == "#FIELDINTEGRATOR") {
    param.read_var("fieldIntegrator", fieldIntegrator);
  } else if (command == "#AVGFIELDB") {
    param.read_var("useAvgFieldB", useAvgFieldB);
    param.read_var("nAvgFieldB", nAvgFieldB);
  } else if (command == "#SELECTPARTICLE") {
    param.read_var("doSelectParticle", doSelectParticle);
    if (doSelectParticle) {
      param.read_var("selectParticleInputFile", selectParticleInputFile);
    }
  }
}

//==========================================================
// Report de-duplicated boundary-condition warnings and abort.
void Pic::report_bc_warnings(const std::string& context) {
  if (bcWarnings_.empty())
    return;

  std::string msg;
  for (const std::string& w : bcWarnings_)
    msg += "\n  - " + w;
  amrex::Abort("Error: " + context + " boundary conditions:" + msg);
}

//==========================================================
// Autofill periodic boundaries from Geometry for field and particle BCs,
// warning on conflicting user configurations.
void Pic::apply_periodicity_autofill(const Geometry& gm) {
  static const char* const dimNames[3] = { "x", "y", "z" };
  const int nSpeciesBC = static_cast<int>(pInfo.pBCs.size());
  const int nBCsSet = static_cast<int>(pInfo.pBCsSet.size());

  for (int d = 0; d < nDim; ++d) {
    if (!gm.isPeriodic(d))
      continue;

    for (int side = 0; side < 2; ++side) {
      if (fieldBCSet_) {
        const int type = bcField.face(d, side);
        if (type != FieldBC::periodic)
          add_bc_warning(
              std::string("dimension ") + dimNames[d] +
              " is periodic (#PERIODICITY), but the field boundary was set "
              "to '" +
              FieldBC::to_string(static_cast<FieldBC::Type>(type)) +
              "'; using 'periodic'.");
      }
      bcField.set(d, side, FieldBC::periodic);

      for (int i = 0; i < nSpeciesBC; ++i) {
        if (i < nBCsSet && pInfo.pBCsSet[i] != 0) {
          const int type = pInfo.pBCs[i].face(d, side);
          if (type != ParticleBC::periodic)
            add_bc_warning(
                std::string("dimension ") + dimNames[d] +
                " is periodic (#PERIODICITY), but the particle boundary of "
                "species " +
                std::to_string(i) + " was set to '" +
                ParticleBC::to_string(static_cast<ParticleBC::Type>(type)) +
                "'; using 'periodic'.");
        }
        pInfo.pBCs[i].set(d, side, ParticleBC::periodic);
      }
    }
  }
}

//==========================================================
// Validate boundary condition consistency across field, particles, and
// geometry. Field vs particle checks are skipped if field BC was not explicitly
// set.
void Pic::validate_bc_pairing(const Geometry& gm) {
  static const char* const dimNames[3] = { "x", "y", "z" };
  static const char* const faceNames[3][2] = { { "x-lo", "x-hi" },
                                               { "y-lo", "y-hi" },
                                               { "z-lo", "z-hi" } };

  const bool isStandalone = domainParameters.isStandalone;
  const bool hasNonPeriodic = !gm.isAllPeriodic();

  if (isStandalone && hasNonPeriodic) {
    if (!fieldBCSet_)
      amrex::Abort(
          "Error: #FIELDBOXBOUNDARY command is required when there are "
          "non-periodic boundaries in standalone mode.");
    if (usePIC && (pInfo.pBCsSet.empty() || pInfo.pBCsSet[0] == 0))
      amrex::Abort(
          "Error: #PARTICLEBOXBOUNDARY command is required when there are "
          "non-periodic boundaries in standalone mode.");
  }

  const int nSpeciesBC = static_cast<int>(pInfo.pBCs.size());
  const int nBCsSet = static_cast<int>(pInfo.pBCsSet.size());
  const int nCheck = std::min(nSpeciesBC, nBCsSet);

  for (int d = 0; d < nDim; ++d) {
    const bool dimPeriodic = gm.isPeriodic(d);
    const int fLo = bcField.face(d, 0);
    const int fHi = bcField.face(d, 1);

    if (!dimPeriodic) {
      if (fLo == FieldBC::periodic)
        add_bc_warning(std::string("field boundary ") + faceNames[d][0] +
                       " is 'periodic' but #PERIODICITY is F for that "
                       "dimension.");
      if (fHi == FieldBC::periodic)
        add_bc_warning(std::string("field boundary ") + faceNames[d][1] +
                       " is 'periodic' but #PERIODICITY is F for that "
                       "dimension.");
    }
    if ((fLo == FieldBC::periodic) != (fHi == FieldBC::periodic))
      add_bc_warning(std::string("field boundary ") + dimNames[d] +
                     " is periodic on one side only; #PERIODICITY applies to "
                     "a whole dimension.");

    for (int side = 0; side < 2; ++side) {
      const auto fType = static_cast<FieldBC::Type>(bcField.face(d, side));
      const char* face = faceNames[d][side];

      if (fType == FieldBC::inflow && !inflowDefined_)
        add_bc_warning(std::string("field boundary ") + face +
                       " is 'inflow' but no #INFLOW block was given; the face "
                       "falls back to a zero-gradient copy.");

      if (fType == FieldBC::wave && !waveBC.active)
        add_bc_warning(std::string("field boundary ") + face +
                       " is 'wave' but no #WAVEBC block was given; the face "
                       "carries no wave source.");

      for (int i = 0; i < nCheck; ++i) {
        if (pInfo.pBCsSet[i] == 0)
          continue;

        const auto pType =
            static_cast<ParticleBC::Type>(pInfo.pBCs[i].face(d, side));

        auto what = [&]() {
          return "species " + std::to_string(i) + " particle boundary " + face;
        };

        if (pType == ParticleBC::periodic && !dimPeriodic)
          add_bc_warning(what() + " is 'periodic' but #PERIODICITY is F for "
                                  "that dimension.");

        if (!fieldBCSet_)
          continue; // Field is default coupled: skip field-pairing checks.

        if (pType == ParticleBC::inflow && fType != FieldBC::inflow &&
            fType != FieldBC::fixed)
          add_bc_warning(what() + " is 'inflow' but the field boundary is '" +
                         FieldBC::to_string(fType) +
                         "'; the injected flux has no upstream field to "
                         "match.");

        if (fType == FieldBC::conducting &&
            (pType == ParticleBC::outflow || pType == ParticleBC::vacuum ||
             pType == ParticleBC::absorb))
          add_bc_warning(what() + " is '" + ParticleBC::to_string(pType) +
                         "' on a 'conducting' wall: particles leave through "
                         "the wall.");

        if (pType == ParticleBC::absorb && fType != FieldBC::absorb)
          add_bc_warning(what() + " is 'absorb' but the field boundary is '" +
                         FieldBC::to_string(fType) +
                         "'; only the particles are absorbed.");
      }
    }
  }

  // In hybrid PIC, only centerB is evolved; wall BCs constrain B and close
  // the ghost ring for E.
  if (useHybridPIC) {
    bool hasWall = false;
    for (int d = 0; d < nDim && !hasWall; ++d) {
      for (int side = 0; side < 2; ++side) {
        const auto t = static_cast<FieldBC::Type>(bcField.face(d, side));
        if (t == FieldBC::conducting || t == FieldBC::wave) {
          hasWall = true;
          break;
        }
      }
    }
    if (hasWall)
      amrex::Print() << "  Note: hybrid solver: only centerB is evolved, so a "
                     << "conducting / wave field boundary constrains "
                     << "B; for the Ohm's-law E it only closes the ghost ring "
                     << "(it is not an independent constraint).\n";
  }
}

//==========================================================
void Pic::post_process_param() {
  fi->set_plasma_charge_and_mass(qomEl);
  nSpecies = fi->get_nS();
  // Species without a #PARTICLEBOXBOUNDARY block keep the default (coupled),
  // or inherit from species 0 if species 0 was specified.
  if (static_cast<int>(pInfo.pBCs.size()) < nSpecies) {
    pInfo.pBCs.resize(nSpecies);
    pInfo.pBCsSet.resize(nSpecies, 0);
  }

  const bool hasSpecies0 = (!pInfo.pBCsSet.empty() && pInfo.pBCsSet[0] != 0);
  if (hasSpecies0) {
    for (int i = 1; i < nSpecies; ++i) {
      if (pInfo.pBCsSet[i] == 0) {
        pInfo.pBCs[i] = pInfo.pBCs[0];
        pInfo.pBCsSet[i] = 1;
      }
    }
  }

  fsolver.mode = (!fsolver.useLaggedLimiter && limiterThetaE != 0)
                     ? FieldSolverMode::NewtonKrylov
                     : FieldSolverMode::GMRES;

  // Classify species: negative charge -> electron, otherwise kinetic ion.
  kineticSpecies_.clear();
  iElectron_ = -1;
  for (int i = 0; i < nSpecies; ++i) {
    // Guard: parts may not be fully populated yet.
    if (i < (int)parts.size() && parts[i] && parts[i]->get_charge() < 0) {
      if (iElectron_ < 0)
        iElectron_ = i;
    } else {
      kineticSpecies_.push_back(i);
    }
  }

  // Hybrid and implicit solver are mutually exclusive.
  if (useHybridPIC)
    solveEM = false;

  // Convert input units to normalized code units.
  if (useHybridPIC) {
    // Resistivity SI->code conversions are deferred to
    // finalize_units_conversion().

    useRK4 = (fieldIntegrator == "rk4");
    if (fieldIntegrator != "rk4" && fieldIntegrator != "ssprk3") {
      amrex::Print() << "  WARNING: unknown #FIELDINTEGRATOR '"
                     << fieldIntegrator << "'; defaulting to 'rk4'\n";
      fieldIntegrator = "rk4";
      useRK4 = true;
    }
    amrex::Print() << "  fieldIntegrator: " << fieldIntegrator << "\n";
    amrex::Print() << "  useAvgFieldB: " << useAvgFieldB
                   << "   nAvgFieldB: " << nAvgFieldB << "\n";
    if (nAvgFieldB < 1)
      nAvgFieldB = 1;
    if (electronTemperatureEV > 0) {
      // Te_code = Te_eV * e / (mp * uNorm_SI^2)
      double unormSI = fi->get_unorm_si();
      electronTemperature = electronTemperatureEV * cUnitChargeSI /
                            (cProtonMassSI * unormSI * unormSI);
      amrex::Print() << "  electronTemperature: " << electronTemperatureEV
                     << " [eV] -> " << electronTemperature << " [code units]\n";
    }

    // Conversion to code units deferred until convert_electron_density0()
    // (Si2NoRho is not yet available here).
    if (rhoMinOhm <= 0)
      rhoMinOhm = 0.0; // resolved to 1e-6*electronDensity0 on first advance
  }
  // report_bc_warnings() is called by Domain after BC autofill and validation.
}

//==========================================================
void Pic::finalize_units_conversion() {
  // Convert input units to code units using normalization factors from fi.
  convert_resistivity();
  convert_electron_density0();
  convert_inflow_state();
}

//==========================================================
void Pic::convert_resistivity() {
  if (!useHybridPIC)
    return;

  const Real Si2NoV = fi->get_Si2NoV();
  const Real Si2NoL = fi->get_Si2NoL();

  // Resistive term eta*J: [eta] = [U]*[L], so
  // eta_code = 4*pi * eta_SI * Si2NoV * Si2NoL.
  if (etaResistivitySI > 0) {
    etaResistivity = fourPI * etaResistivitySI * Si2NoV * Si2NoL;
    amrex::Print() << "  etaResistivity: " << etaResistivitySI << " [m^2/s] -> "
                   << etaResistivity << " [code units]"
                   << "  (Si2NoV = " << Si2NoV << ", Si2NoL = " << Si2NoL
                   << ")\n";
  }

  // Hyper-resistive term eta_h*nabla^2 J: [eta_h] = [U]*[L]^3, so
  // eta_h_code = 4*pi * eta_h_SI * Si2NoV * Si2NoL^3. A single physical value
  // is used on every level (the same choice as grid mode in update_B_hybrid).
  if (etaHyperSI > 0 && etaHyperMode == "si") {
    const Real etaHyper = fourPI * etaHyperSI * Si2NoV * std::pow(Si2NoL, 3);
    for (int iLev = 0; iLev < n_lev_max(); ++iLev)
      etaHyperLev[iLev] = etaHyper;
    amrex::Print() << "  etaHyper: " << etaHyperSI << " [m^4/s, si] -> "
                   << etaHyper << " [code units]\n";
  }

  // Guard against uninitialized normalization producing non-positive
  // coefficients.
  if ((etaResistivitySI > 0 && !(etaResistivity > 0)) ||
      (etaHyperSI > 0 && etaHyperMode == "si" &&
       (etaHyperLev.empty() || !(etaHyperLev[0] > 0)))) {
    amrex::Abort("Pic::convert_resistivity: the SI->code conversion produced a "
                 "non-positive resistivity. Check the normalization "
                 "(#NORMALIZATION lNormSI / uNormSI).");
  }
}

//==========================================================
void Pic::convert_electron_density0() {
  // Convert electron density from amu/cc to code units.
  electronDensity0 =
      electronDensity0In * 1.0e6 * cProtonMassSI * fi->get_Si2NoRho();

  // Auto density floor in code units.
  if (rhoMinOhm <= 0)
    rhoMinOhm = 1.0e-6 * electronDensity0;

  amrex::Print() << "  electronDensity0: " << electronDensity0In
                 << " [amu/cc] -> " << electronDensity0
                 << " [code units]  (Si2NoRho = " << fi->get_Si2NoRho()
                 << ")\n";
}

//==========================================================
void Pic::convert_inflow_state() {
  if (!inflowDefined_)
    return;

  // Convert density from amu/cc to code units.
  const double Si2NoRho = fi->get_Si2NoRho();
  inflowRho_ *= 1.0e6 * cProtonMassSI * Si2NoRho;

  // Convert velocity from km/s to code units.
  const double vFactor = 1.0e3 * fi->get_Si2NoV();
  inflowUx_ *= vFactor;
  inflowUy_ *= vFactor;
  inflowUz_ *= vFactor;

  // Convert temperature T [K] to code units (kT / (m_p * uNorm^2)).
  const double unormSI = fi->get_unorm_si();
  inflowT_ = cBoltzmannSI * inflowT_ / (cProtonMassSI * unormSI * unormSI);

  // Publish converted state to FluidInterface for boundary particle injection.
  FluidInterfaceParameters::InflowVel baseVel;
  baseVel.nDens = inflowRho_;
  baseVel.ux = inflowUx_;
  baseVel.uy = inflowUy_;
  baseVel.uz = inflowUz_;
  baseVel.vth = 0.0;

  amrex::Vector<FluidInterfaceParameters::InflowVel> stateVec(nSpecies,
                                                              baseVel);
  const int nParts = static_cast<int>(parts.size());
  for (int iS = 0; iS < nSpecies; ++iS) {
    const double mass_i =
        (iS < nParts && parts[iS]) ? parts[iS]->get_mass() : 1.0;
    if (inflowT_ > 0 && mass_i > 0)
      stateVec[iS].vth = std::sqrt(inflowT_ / mass_i);
  }
  fi->set_inflow_state(stateVec);
  fi->set_inflow_defined(true);

  amrex::Print() << "  #INFLOW state (code units):"
                 << " n=" << inflowRho_ << " u=(" << inflowUx_ << ","
                 << inflowUy_ << "," << inflowUz_ << ")"
                 << " vth=" << (inflowT_ > 0 ? std::sqrt(inflowT_) : 0.0)
                 << "  (Si2NoRho=" << Si2NoRho
                 << ", Si2NoV=" << fi->get_Si2NoV() << ")\n";
}
