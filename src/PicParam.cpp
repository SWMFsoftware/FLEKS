#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <memory>
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
    if (pInfo.pBCs.empty()) {
      pInfo.pBCs.resize(1);
      pInfo.pBCsSet.resize(1, 0);
    }
    pInfo.pBCsSet[0] = 1;

    std::string lo, hi;
    for (int i = 0; i < nDim; ++i) {
      param.read_var("particleBoxBoundaryLo", lo);
      param.read_var("particleBoxBoundaryHi", hi);
      pInfo.pBCs[0].set(i, 0, ParticleBC::parse(lo));
      pInfo.pBCs[0].set(i, 1, ParticleBC::parse(hi));
    }
    for (int s = 1; s < static_cast<int>(pInfo.pBCs.size()); ++s) {
      pInfo.pBCs[s] = pInfo.pBCs[0];
      pInfo.pBCsSet[s] = 1;
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
    update_bc_flags();
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
  } else if (command == "#DIPOLE") {
    // Static dipole of the inner body. Everything is SI; convert_intrinsic_B()
    // turns it into code units once the normalization is known.
    double bEq = 0.0, theta = 0.0, phi = 0.0, rRef = -1.0;
    param.read_var("strength", bEq);
    param.read_var("theta", theta);
    param.read_var("phi", phi);
    param.read_var("rRef", rRef);

    if (intrinsicB_ == nullptr)
      intrinsicB_ = std::make_unique<IntrinsicBField>();
    intrinsicB_->read_dipole(bEq, theta, phi, rRef);

    Print() << "  intrinsic field: dipole, strength = " << bEq
            << " [nT] at rRef, tilt = (" << theta << ", " << phi << ") [deg], "
            << "rRef = " << rRef << " [m, <0 = from #BODY/#PLANETRADIUS]\n";
  } else if (command == "#CRUSTALFIELD") {
    // Static crustal field from a BATSRUS spherical harmonic coefficient file.
    std::string fileName;
    int nMax = 0;
    param.read_var("fileName", fileName);
    param.read_var("nMax", nMax);

    if (nMax < 2)
      Abort("Error: #CRUSTALFIELD nMax must be at least 2 (the number of "
            "harmonic degrees, i.e. the BATSRUS NNm), got " +
            std::to_string(nMax) + ".");

    if (intrinsicB_ == nullptr)
      intrinsicB_ = std::make_unique<IntrinsicBField>();
    intrinsicB_->read_crustal(fileName, nMax);

    Print() << "  intrinsic field: crustal, file = " << fileName
            << ", nMax = " << nMax << "\n";
  } else if (command == "#BODY") {
    std::string type;
    param.read_var("type", type);
    if (type != "sphere")
      Abort("Error: #BODY type '" + type +
            "' is not supported. Only 'sphere' is implemented.");

    Real radius;
    param.read_var("radius", radius);

    Real center[nDim];
    for (int i = 0; i < nDim; i++)
      param.read_var("center", center[i]);

    // The geometry is in code units (one length unit is lNormSI metres),
    // like #REGION, and unlike #PLANETRADIUS which is in SI.
    set_body(center, radius);

    Print() << "  inner body: sphere, radius = " << bodyRadius
            << ", center = (";
    for (int i = 0; i < nDim; i++)
      Print() << (i > 0 ? ", " : "") << bodyCenter[i];
    Print() << ") [code units]\n";
  } else if (command == "#BODYBOUNDARY") {
    std::string particle, field;
    param.read_var("particleBoundary", particle);
    param.read_var("fieldBoundary", field);

    bodyParticleBC = ParticleBC::parse(particle);
    if (bodyParticleBC != ParticleBC::absorb &&
        bodyParticleBC != ParticleBC::reflect)
      Abort("Error: #BODYBOUNDARY particleBoundary '" + particle +
            "' is not supported for the inner body. Accepted values: "
            "absorb, reflect.");

    bodyFieldBC = BodyFieldBC::parse(field);
    bodyBoundarySet_ = true;

    Print() << "  inner body BC: particles = "
            << ParticleBC::to_string(bodyParticleBC)
            << ", fields = " << BodyFieldBC::to_string(bodyFieldBC) << "\n";
  } else if (command == "#REGIONRESISTIVITY") {
    RegionalResistivityConfig cfg;
    param.read_var("region", cfg.regionStr);
    param.read_var("etaResistivitySI", cfg.etaSI);
    regionalResistivityConfigs.push_back(cfg);
  } else if (command == "#REGIONHYPERRESISTIVITY") {
    RegionalHyperResistivityConfig cfg;
    param.read_var("region", cfg.regionStr);
    param.read_var("etaHyperSI", cfg.etaSI);
    param.read_var("etaHyperMode", cfg.mode);
    param.read_var("etaHyperCh", cfg.ch);
    regionalHyperResistivityConfigs.push_back(cfg);
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
      Abort("Unknown #TESTCASE name '" + testcase +
            ". Registered names: " + known + ".");
    }
    ic_->read_param(param);
  } else if (command == "#WAVEIC") {
    if (!ic_) {
      Abort("The #WAVEIC block must follow a #TESTCASE that selects a "
            "wave initial condition (waveic / lightwave / hybridwave / "
            "alfvenpulse / convectionwave / ionacousticwave).");
    }
    ic_->read_param(param);
  } else if (command == "#FADEEVIC") {
    if (!ic_ || std::string(ic_->name()) != "fadeev") {
      Abort("The #FADEEVIC block must follow a #TESTCASE that selects "
            "the fadeev (magnetic reconnection) initial condition.");
    }
    ic_->read_param(param);
  } else if (command == "#GEMIC") {
    if (!ic_ || std::string(ic_->name()) != "gem") {
      Abort("The #GEMIC block must follow a #TESTCASE that selects "
            "the gem initial condition.");
    }
    ic_->read_param(param);
  } else if (command == "#FORCEFREEIC") {
    if (!ic_ || std::string(ic_->name()) != "forcefree") {
      Abort("The #FORCEFREEIC block must follow a #TESTCASE that selects "
            "the forcefree (magnetic reconnection) initial condition.");
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
  } else if (command == "#ELECTRONPRESSURE") {
    param.read_var("useElectronPressureEq", useElectronPressureEq);
    if (useElectronPressureEq) {
      param.read_var("heatCondKappa0SI", heatCondKappa0SI);
      param.read_var("coulombLog", coulombLog);
      param.read_var("fieldAlignedConduction", fieldAlignedConduction);
      param.read_var("fieldAlignedFraction", fieldAlignedFraction);
      param.read_var("fieldAlignedBMinSI", fieldAlignedBMinSI);
      param.read_var("heatFluxLimiter", heatFluxLimiter);
      param.read_var("peMin", peMin);
      param.read_var("ambipolarInStages", ambipolarInStages);
    }
  } else if (command == "#ELECTRONADVECTION") {
    param.read_var("peAdvectionLimiter", peAdvectionLimiter);
    param.read_var("peCompressionScheme", peCompressionScheme);
  } else if (command == "#ELECTRONCONDUCTION") {
    param.read_var("heatCondMethod", heatCondMethod);
    param.read_var("nCondIter", nCondIter);
    param.read_var("nCondSubcycleMax", nCondSubcycleMax);
  } else if (command == "#ELECTRONCOLLISION") {
    param.read_var("useHeatExchange", useHeatExchange);
    if (useHeatExchange) {
      param.read_var("collisionFactor", collisionFactor);
    }
  } else if (command == "#FIELDINTEGRATOR") {
    param.read_var("fieldIntegrator", fieldIntegrator);
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
  Abort("Error: " + context + " boundary conditions:" + msg);
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
      Abort("Error: #FIELDBOXBOUNDARY command is required when there are "
            "non-periodic boundaries in standalone mode.");
    if (usePIC && (pInfo.pBCsSet.empty() || pInfo.pBCsSet[0] == 0))
      Abort("Error: #PARTICLEBOXBOUNDARY command is required when there are "
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
      Print() << "  Note: hybrid solver: only centerB is evolved, so a "
              << "conducting / wave field boundary constrains "
              << "B; for the Ohm's-law E it only closes the ghost ring "
              << "(it is not an independent constraint).\n";
  }

  update_bc_flags();
}

//==========================================================
void Pic::post_process_param() {
  // Validate raw user input before deriving defaults or converting units.
  // These values are used in denominators and dispatch solver branches, so
  // silently correcting them can produce a run with different physics than
  // the input deck describes.
  if (nBSubcycle < 1)
    Abort("Invalid #BSUBCYCLE: nBSubcycle must be at least 1.");
  if (electronGamma <= 0)
    Abort("Invalid #ELECTRONTEMPERATURE: electronGamma must be > 0.");
  if (electronDensity0In <= 0)
    Abort("Invalid #ELECTRONTEMPERATURE: electronDensity0 must be > 0.");
  if (etaHyperSI < 0)
    Abort("Invalid #HYPERRESISTIVITY: etaHyperSI must be non-negative.");
  if (etaHyperCh < 0)
    Abort("Invalid #HYPERRESISTIVITY: etaHyperCh must be non-negative.");
  if (rhoMinOhm < 0)
    Abort("Invalid #MINIMUMDENSITY: rhoMinOhm must be non-negative.");
  if (useElectronPressureEq) {
    if (!useHybridPIC)
      Abort("Invalid #ELECTRONPRESSURE: the evolved electron pressure "
            "equation is only supported by the hybrid PIC solver "
            "(#HYBRIDPIC T).");
    // gamma_e - 1 multiplies both the pdV term and every source term, so an
    // isothermal index would silently reduce the equation to pure advection.
    if (electronGamma <= 1.0)
      Abort("Invalid #ELECTRONPRESSURE: #ELECTRONTEMPERATURE electronGamma "
            "must be > 1 for the evolved electron pressure equation "
            "(gamma_e = 1 removes the pdV term and the heat-flux source).");
    if (electronTemperatureEV <= 0)
      Abort("Invalid #ELECTRONPRESSURE: #ELECTRONTEMPERATURE "
            "electronTemperature must be > 0 (it sets the initial Pe).");
    if (heatCondKappa0SI < 0)
      Abort("Invalid #ELECTRONPRESSURE: heatCondKappa0SI must be "
            "non-negative.");
    if (coulombLog <= 0)
      Abort("Invalid #ELECTRONPRESSURE: coulombLog must be positive.");
    if (fieldAlignedBMinSI < 0)
      Abort("Invalid #ELECTRONPRESSURE: fieldAlignedBMinSI must be "
            "non-negative.");
    if (fieldAlignedFraction < 0 || fieldAlignedFraction > 1)
      Abort("Invalid #ELECTRONPRESSURE: fieldAlignedFraction must be "
            "between 0 (isotropic) and 1 (field aligned).");
    if (heatFluxLimiter < 0)
      Abort("Invalid #ELECTRONPRESSURE: heatFluxLimiter must be "
            "non-negative.");
    if (peMin < 0)
      Abort("Invalid #ELECTRONPRESSURE: peMin must be non-negative.");

    if (peAdvectionLimiter == "upwind1") {
      peLimiterType = 0;
    } else if (peAdvectionLimiter == "minmod") {
      peLimiterType = 1;
    } else if (peAdvectionLimiter == "vanleer") {
      peLimiterType = 2;
    } else if (peAdvectionLimiter == "mc") {
      peLimiterType = 3;
    } else {
      Abort("Invalid #ELECTRONADVECTION: peAdvectionLimiter must be "
            "'upwind1', 'minmod', 'vanleer', or 'mc'.");
    }

    if (peCompressionScheme == "exponential") {
      peCompressionExp = true;
    } else if (peCompressionScheme == "explicit") {
      peCompressionExp = false;
    } else {
      Abort("Invalid #ELECTRONADVECTION: peCompressionScheme must be "
            "'exponential' or 'explicit'.");
    }

    if (heatCondMethod != "point-implicit" && heatCondMethod != "subcycle") {
      Abort("Invalid #ELECTRONCONDUCTION: heatCondMethod must be "
            "'point-implicit' or 'subcycle'.");
    }
    if (nCondIter < 1)
      Abort("Invalid #ELECTRONCONDUCTION: nCondIter must be at least 1.");
    if (nCondSubcycleMax < 1)
      Abort(
          "Invalid #ELECTRONCONDUCTION: nCondSubcycleMax must be at least 1.");

    if (collisionFactor < 0)
      Abort(
          "Invalid #ELECTRONCOLLISION: collisionFactor must be non-negative.");

    // Sized here rather than in the constructor: the command has not been
    // read when Pic is constructed.
    centerPeState.resize(n_lev_max());
    centerPeRho.resize(n_lev_max());
    centerPeTe.resize(n_lev_max());
    nodePeVec.resize(n_lev_max());
    nodePeRho.resize(n_lev_max());
    nodePeAux.resize(n_lev_max());
  }
  if (fieldIntegrator != "rk4" && fieldIntegrator != "ssprk3")
    Abort("Invalid #FIELDINTEGRATOR '" + fieldIntegrator +
          "'. Expected 'rk4' or 'ssprk3'.");
  if (etaHyperMode != "si" && etaHyperMode != "grid")
    Abort("Invalid #HYPERRESISTIVITY etaHyperMode '" + etaHyperMode +
          "'. Expected 'si' or 'grid'.");

  if (useBody) {
    if (bodyRadius <= 0)
      Abort("Invalid #BODY: radius must be positive.");

    if (n_lev_max() > 1 && refineRegions && !refineRegions->empty())
      Print() << "  Warning: #BODY has not been verified with AMR "
              << "(refinement regions are defined).\n";

    // The body has to be strictly inside the domain: a body crossing a
    // domain face (or a periodic face) would need a mask that is consistent
    // across the periodic images, which is not implemented. The invariant
    // direction of a fake-2D run (one cell) is not checked, because the body
    // necessarily extends beyond it.
    const auto plo = Geom(0).ProbLo();
    const auto phi = Geom(0).ProbHi();
    const auto& dom = Geom(0).Domain();
    for (int i = 0; i < nDim; i++) {
      if (dom.length(i) <= 1)
        continue;
      if (bodyCenter[i] - bodyRadius <= plo[i] ||
          bodyCenter[i] + bodyRadius >= phi[i])
        Abort("Invalid #BODY: the body must be strictly inside the "
              "simulation domain.");
    }
  } else if (bodyBoundarySet_) {
    // #BODYBOUNDARY only has a meaning together with #BODY.
    Print() << "  Warning: #BODYBOUNDARY is ignored because no #BODY is "
            << "defined.\n";
    bodyParticleBC = ParticleBC::absorb;
    bodyFieldBC = BodyFieldBC::linetied;
    bodyBoundarySet_ = false;
  }

  // A conducting inner body pins the radial component of the *evolved* field
  // to zero on the body surface (project_body_B). A planetary intrinsic field
  // has a radial component that necessarily crosses the surface, so the two
  // conditions contradict each other. Reject the combination instead of
  // silently dropping one of them; use 'linetied' or 'insulating' instead.
  if (useBody && bodyFieldBC == BodyFieldBC::conducting &&
      intrinsicB_ != nullptr && intrinsicB_->is_active()) {
    Abort("Invalid combination: #BODYBOUNDARY fieldBoundary 'conducting' "
          "cannot be used together with an intrinsic magnetic field "
          "(#DIPOLE / #CRUSTALFIELD). Use 'linetied' or 'insulating'.");
  }

  hasRegionalResistivity_ = !regionalResistivityConfigs.empty();
  hasRegionalHyper_ = !regionalHyperResistivityConfigs.empty();

  if (hasRegionalResistivity_ || hasRegionalHyper_) {
    if (!useHybridPIC) {
      Abort("Invalid configuration: #REGIONRESISTIVITY and "
            "#REGIONHYPERRESISTIVITY "
            "are only supported for the hybrid PIC solver (#HYBRIDPIC T).");
    }
  }

  if (useHybridPIC) {
    if (hasRegionalResistivity_)
      nodeEtaRegional.resize(n_lev_max());
    if (hasRegionalHyper_)
      nodeEtaHyperRegional.resize(n_lev_max());
  }

  for (const auto& cfg : regionalResistivityConfigs) {
    if (cfg.regionStr.empty())
      Abort("Invalid #REGIONRESISTIVITY: region expression cannot be empty.");
    if (cfg.etaSI < 0.0)
      Abort("Invalid #REGIONRESISTIVITY: etaResistivitySI must be >= 0.");
  }

  for (const auto& cfg : regionalHyperResistivityConfigs) {
    if (cfg.regionStr.empty())
      Abort("Invalid #REGIONHYPERRESISTIVITY: region expression cannot be "
            "empty.");
    if (cfg.mode != "si" && cfg.mode != "grid")
      Abort("Invalid #REGIONHYPERRESISTIVITY etaHyperMode '" + cfg.mode +
            "'. Expected 'si' or 'grid'.");
    if (cfg.etaSI < 0.0)
      Abort("Invalid #REGIONHYPERRESISTIVITY: etaHyperSI must be >= 0.");
    if (cfg.ch < 0.0)
      Abort("Invalid #REGIONHYPERRESISTIVITY: etaHyperCh must be >= 0.");
  }

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
    Print() << "  fieldIntegrator: " << fieldIntegrator << "\n";
    if (electronTemperatureEV > 0) {
      // Te_code = Te_eV * e / (mp * uNorm_SI^2)
      double unormSI = fi->get_unorm_si();
      electronTemperature = electronTemperatureEV * cUnitChargeSI /
                            (cProtonMassSI * unormSI * unormSI);
      Print() << "  electronTemperature: " << electronTemperatureEV
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
  convert_electron_heat_conduction();
  convert_electron_collision();
  convert_inflow_state();
  convert_intrinsic_B();
}

//==========================================================
void Pic::convert_intrinsic_B() {
  if (intrinsicB_ == nullptr)
    return;

  if (!intrinsicB_->is_active()) {
    // Neither model was actually requested; drop the object so that every
    // consumer keeps its fast path.
    intrinsicB_.reset();
    return;
  }

  const double bodyCenterTmp[3] = { bodyCenter[ix_], bodyCenter[iy_],
                                    bodyCenter[iz_] };

  intrinsicB_->convert_units(fi->get_Si2NoB(), fi->get_Si2NoL(),
                             fi->get_rPlanet_SI(), bodyCenterTmp, bodyRadius,
                             useBody, nDim);

  Print() << intrinsicB_->describe() << "\n";
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
    Print() << "  etaResistivity: " << etaResistivitySI << " [m^2/s] -> "
            << etaResistivity << " [code units]"
            << "  (Si2NoV = " << Si2NoV << ", Si2NoL = " << Si2NoL << ")\n";
  }

  // Hyper-resistive term eta_h*nabla^2 J: [eta_h] = [U]*[L]^3, so
  // eta_h_code = 4*pi * eta_h_SI * Si2NoV * Si2NoL^3. A single physical value
  // is used on every level (the same choice as grid mode in update_B_hybrid).
  if (etaHyperSI > 0 && etaHyperMode == "si") {
    const Real etaHyper = fourPI * etaHyperSI * Si2NoV * std::pow(Si2NoL, 3);
    for (int iLev = 0; iLev < n_lev_max(); ++iLev)
      etaHyperLev[iLev] = etaHyper;
    Print() << "  etaHyper: " << etaHyperSI << " [m^4/s, si] -> " << etaHyper
            << " [code units]\n";
  }

  // Regional resistivity
  for (auto& cfg : regionalResistivityConfigs) {
    if (cfg.etaSI > 0) {
      cfg.etaCode = fourPI * cfg.etaSI * Si2NoV * Si2NoL;
      Print() << "  regionalResistivity [" << cfg.regionStr
              << "]: " << cfg.etaSI << " [m^2/s] -> " << cfg.etaCode
              << " [code units]\n";
    }
  }

  // Regional hyper-resistivity (si mode)
  for (auto& cfg : regionalHyperResistivityConfigs) {
    cfg.etaLev.resize(n_lev_max(), 0.0);
    if (cfg.etaSI > 0 && cfg.mode == "si") {
      const Real etaHyper = fourPI * cfg.etaSI * Si2NoV * std::pow(Si2NoL, 3);
      for (int iLev = 0; iLev < n_lev_max(); ++iLev)
        cfg.etaLev[iLev] = etaHyper;
      Print() << "  regionalHyper [" << cfg.regionStr << "]: " << cfg.etaSI
              << " [m^4/s, si] -> " << etaHyper << " [code units]\n";
    }
  }

  // Guard against uninitialized normalization producing non-positive
  // coefficients.
  if ((etaResistivitySI > 0 && !(etaResistivity > 0)) ||
      (etaHyperSI > 0 && etaHyperMode == "si" &&
       (etaHyperLev.empty() || !(etaHyperLev[0] > 0)))) {
    Abort("Pic::convert_resistivity: the SI->code conversion produced a "
          "non-positive resistivity. Check the normalization "
          "(#NORMALIZATION lNormSI / uNormSI).");
  }

  for (const auto& cfg : regionalResistivityConfigs) {
    if (cfg.etaSI > 0 && !(cfg.etaCode > 0)) {
      Abort("Pic::convert_resistivity: #REGIONRESISTIVITY produced a "
            "non-positive resistivity. Check the normalization "
            "(#NORMALIZATION).");
    }
  }
  for (const auto& cfg : regionalHyperResistivityConfigs) {
    if (cfg.etaSI > 0 && cfg.mode == "si" &&
        (cfg.etaLev.empty() || !(cfg.etaLev[0] > 0))) {
      Abort("Pic::convert_resistivity: #REGIONHYPERRESISTIVITY produced a "
            "non-positive hyper-resistivity. Check the normalization "
            "(#NORMALIZATION).");
    }
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

  Print() << "  electronDensity0: " << electronDensity0In << " [amu/cc] -> "
          << electronDensity0
          << " [code units]  (Si2NoRho = " << fi->get_Si2NoRho() << ")\n";
}

//==========================================================
// SI -> code conversion of the Spitzer electron heat conduction coefficient.
//
// The heat flux is q = -kappa * grad(Te) with kappa = kappa0 * Te^2.5, so
//   [kappa] = [energy flux] / ([T]^3.5 / [L]) = [EnergyDens]*[U]*[L]/[T]^3.5
// (the same expression BATSRUS uses in ModHeatConduction.f90). In FLEKS code
// units the energy density is rho*U^2 and the temperature is carried as
// k_B*T/m_p in units of uNorm^2 (see the electronTemperature conversion
// above), hence
//   Si2No(EnergyDens) = Si2NoRho * Si2NoV^2
//   Si2No(T)          = cBoltzmannSI / (cProtonMassSI * uNormSI^2)
void Pic::convert_electron_heat_conduction() {
  if (!useElectronPressureEq)
    return;

  const Real Si2NoV = fi->get_Si2NoV();
  const Real Si2NoL = fi->get_Si2NoL();
  const Real Si2NoRho = fi->get_Si2NoRho();
  const Real uNormSI = fi->get_unorm_si();

  const Real Si2NoEnergyDens = Si2NoRho * Si2NoV * Si2NoV;
  const Real Si2NoTemperature =
      cBoltzmannSI / (cProtonMassSI * uNormSI * uNormSI);

  const Real kappa0SI = heatCondKappa0SI * (coulombLog / 20.0);
  heatCondKappa0 = kappa0SI * Si2NoEnergyDens * Si2NoV * Si2NoL /
                   std::pow(Si2NoTemperature, 3.5);

  // heatCondKappa0SI = 0 is legitimate: it switches the electron heat flux
  // off and leaves advection + pdV. Only a positive input that converts to a
  // non-positive code value is an error.
  // Minimum field strength for the field-aligned direction to be meaningful.
  // Without it a globally unmagnetized run still carries round-off level B,
  // whose normalised direction is arbitrary and would silently steer the heat
  // flux along a noise direction.
  fieldAlignedBMin = fieldAlignedBMinSI * fi->get_Si2NoB();

  if (heatCondKappa0SI > 0 && !(heatCondKappa0 > 0))
    Abort("Pic::convert_electron_heat_conduction: the SI->code conversion "
          "produced a non-positive heat conduction coefficient. Check the "
          "normalization (#NORMALIZATION lNormSI / uNormSI).");

  Print() << "  electron heat conduction: kappa0 = " << heatCondKappa0SI
          << " [W/(m K^(7/2))] * (coulombLog/20 = " << coulombLog / 20.0
          << ") -> " << heatCondKappa0 << " [code units]\n";
  Print() << "    gamma_e = " << electronGamma
          << (fieldAlignedConduction
                  ? std::string(", field aligned (fraction ") +
                        std::to_string(fieldAlignedFraction) + ")"
                  : std::string(", isotropic"))
          << (heatFluxLimiter > 0 ? ", free-streaming limiter " +
                                        std::to_string(heatFluxLimiter)
                                  : ", no heat-flux limiter")
          << "\n";
}

//==========================================================
// SI -> code conversion of the electron-ion thermal equilibration coefficient
// following the BATSRUS Braginskii formulation.
void Pic::convert_electron_collision() {
  if (!useElectronPressureEq || !useHeatExchange)
    return;

  const Real Si2NoV = fi->get_Si2NoV();
  const Real Si2NoL = fi->get_Si2NoL();
  const Real Si2NoRho = fi->get_Si2NoRho();
  const Real uNormSI = fi->get_unorm_si();

  const Real reducedMassSI =
      (cElectronMassSI * cProtonMassSI) / (cElectronMassSI + cProtonMassSI);
  const Real twoPiKB = 2.0 * dPI * cBoltzmannSI;
  const Real e2OverEps = (cUnitChargeSI * cUnitChargeSI) / cEps0SI;

  const Real coefSI = coulombLog * std::sqrt(reducedMassSI / cProtonMassSI) *
                      (e2OverEps * e2OverEps) / (3.0 * std::pow(twoPiKB, 1.5));

  const Real Si2NoT = Si2NoL / Si2NoV;
  const Real No2SiN = 1.0 / (Si2NoRho * cProtonMassSI);
  const Real Si2NoTemperature =
      cBoltzmannSI / (cProtonMassSI * uNormSI * uNormSI);

  collisionCoefEi = 2.0 * collisionFactor * coefSI * No2SiN *
                    std::pow(Si2NoTemperature, 1.5) / Si2NoT;

  Print() << "  electron-ion collision (heat exchange): factor = "
          << collisionFactor << ", coulombLog = " << coulombLog
          << " -> collisionCoefEi = " << collisionCoefEi << " [code units]\n";
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

  Vector<FluidInterfaceParameters::InflowVel> stateVec(nSpecies, baseVel);
  const int nParts = static_cast<int>(parts.size());
  for (int iS = 0; iS < nSpecies; ++iS) {
    const double mass_i =
        (iS < nParts && parts[iS]) ? parts[iS]->get_mass() : 1.0;
    if (inflowT_ > 0 && mass_i > 0)
      stateVec[iS].vth = std::sqrt(inflowT_ / mass_i);
  }
  fi->set_inflow_state(stateVec);
  fi->set_inflow_defined(true);

  Print() << "  #INFLOW state (code units):"
          << " n=" << inflowRho_ << " u=(" << inflowUx_ << "," << inflowUy_
          << "," << inflowUz_ << ")"
          << " vth=" << (inflowT_ > 0 ? std::sqrt(inflowT_) : 0.0)
          << "  (Si2NoRho=" << Si2NoRho << ", Si2NoV=" << fi->get_Si2NoV()
          << ")\n";

  const auto& unif = fi->get_uniform_state();
  if (nSpecies > 1 && !unif.empty()) {
    const double rawInflowN = inflowRho_ / (1.0e6 * cProtonMassSI * Si2NoRho);
    for (int iS = 0; iS < nSpecies; ++iS) {
      if (iS * 5 < static_cast<int>(unif.size()) && iS < fi->get_nS() &&
          fi->get_species_mass(iS) > 0.0) {
        const double speciesN =
            unif[iS * 5] / (fi->get_species_mass(iS) * cProtonMassSI * 1.0e6);
        if (std::abs(speciesN - rawInflowN) >
            1e-4 * std::max(speciesN, rawInflowN)) {
          Print()
              << "  Warning: #INFLOW supplies a single uniform number density "
              << "(n=" << rawInflowN << " /cc) for all species, but species "
              << iS << " has #UNIFORMSTATE density n=" << speciesN << " /cc.\n";
          break;
        }
      }
    }
  }
}
