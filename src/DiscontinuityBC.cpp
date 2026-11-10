#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>

#include <AMReX_Print.H>

#include "Constants.h"
#include "DiscontinuityBC.h"

namespace {

inline amrex::Real dot3(const amrex::Real a[3], const amrex::Real b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline void cross3(const amrex::Real a[3], const amrex::Real b[3],
                   amrex::Real out[3]) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}

inline amrex::Real norm3(const amrex::Real a[3]) {
  return std::sqrt(dot3(a, a));
}

// Compute mass-weighted bulk velocity from per-species states.
inline void compute_bulk_velocity(
    const std::vector<DiscontinuitySpeciesState>& states, bool isHybridPIC,
    amrex::Real uBulk[3]) {
  amrex::Real rhoSum = 0.0;
  uBulk[0] = 0.0;
  uBulk[1] = 0.0;
  uBulk[2] = 0.0;
  for (const auto& s : states) {
    if (s.nDens <= 0.0 || s.mass <= 0.0 || (isHybridPIC && s.charge < 0.0))
      continue;
    const amrex::Real rhoS = s.nDens * s.mass;
    rhoSum += rhoS;
    uBulk[0] += rhoS * s.ux;
    uBulk[1] += rhoS * s.uy;
    uBulk[2] += rhoS * s.uz;
  }
  if (rhoSum > 0.0) {
    uBulk[0] /= rhoSum;
    uBulk[1] /= rhoSum;
    uBulk[2] /= rhoSum;
  }
}

// Compute total mass density rho_m = sum_s n_s * m_s in code units.
inline amrex::Real compute_mass_density(
    const std::vector<DiscontinuitySpeciesState>& states, bool isHybridPIC) {
  amrex::Real rhoSum = 0.0;
  for (const auto& s : states) {
    if (s.nDens > 0.0 && s.mass > 0.0 && !(isHybridPIC && s.charge < 0.0))
      rhoSum += s.nDens * s.mass;
  }
  return rhoSum;
}

// Compute hybrid fluid electron pressure P_e in code units from ion states.
inline amrex::Real compute_hybrid_pe(
    const std::vector<DiscontinuitySpeciesState>& states, bool useHybridPe,
    amrex::Real te0Code, amrex::Real ne0Code, amrex::Real gammaE) {
  if (!useHybridPe || te0Code <= 0.0 || ne0Code <= 0.0)
    return 0.0;
  amrex::Real ne = 0.0;
  for (const auto& s : states) {
    if (s.nDens > 0.0 && s.charge > 0.0)
      ne += s.charge * s.nDens;
  }
  if (ne <= 0.0)
    return 0.0;
  const amrex::Real pe0 = ne0Code * te0Code;
  return pe0 * std::pow(ne / ne0Code, gammaE);
}

// Compute total thermal + magnetic pressure in FLEKS CGS-like code units:
//   P_tot = sum_s n_s * m_s * vth_s^2 + P_e(hybrid) + |B|^2 / (8*pi).
inline amrex::Real compute_total_pressure(
    const std::vector<DiscontinuitySpeciesState>& states,
    const amrex::Real b[3], bool isHybridPIC, bool useHybridPe,
    amrex::Real te0Code, amrex::Real ne0Code, amrex::Real gammaE) {
  amrex::Real pTh = 0.0;
  for (const auto& s : states) {
    if (s.nDens > 0.0 && s.mass > 0.0 && !(isHybridPIC && s.charge < 0.0))
      pTh += s.nDens * s.mass * s.vth * s.vth;
  }
  pTh += compute_hybrid_pe(states, useHybridPe, te0Code, ne0Code, gammaE);
  const amrex::Real bSq = dot3(b, b);
  const amrex::Real pMag = bSq / (2.0 * fourPI);
  return pTh + pMag;
}

} // namespace

void DiscontinuityBoundaryManager::read_param(ReadParam& param) {
  active = true;
  isCodeUnits = false;

  rho2In.clear();
  ux2In.clear();
  uy2In.clear();
  uz2In.clear();
  t2In.clear();

  param.read_var("typeDiscontinuity", typeStr);
  for (char& c : typeStr)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

  if (typeStr == "td") {
    type = DiscontinuityType::TD;
    rdRotationSign = 1;
  } else if (typeStr == "rd" || typeStr == "rd+") {
    type = DiscontinuityType::RD;
    rdRotationSign = 1;
  } else if (typeStr == "rd-") {
    type = DiscontinuityType::RD;
    rdRotationSign = -1;
  } else if (typeStr == "general") {
    type = DiscontinuityType::General;
    rdRotationSign = 1;
  } else {
    amrex::Abort("Invalid #DISCONTINUITY: typeDiscontinuity must be 'td', "
                 "'rd', 'rd+', 'rd-', or 'general', got '" +
                 typeStr + "'.");
  }

  param.read_var("normalX", normal[0]);
  param.read_var("normalY", normal[1]);
  param.read_var("normalZ", normal[2]);

  param.read_var("width", width);
  if (width < 0.0) {
    amrex::Abort("Invalid #DISCONTINUITY: width must be >= 0.");
  }

  param.read_var("tStart", tStartSI);

  // Read at least one species State 2 block (rho [1/cc], ux/uy/uz [km/s], T
  // [K]), then any optional per-species blocks for additional species.
  amrex::Real rho = 0.0, ux = 0.0, uy = 0.0, uz = 0.0, temp = 0.0;
  param.read_var("rho", rho);
  param.read_var("ux", ux);
  param.read_var("uy", uy);
  param.read_var("uz", uz);
  param.read_var("T", temp);
  rho2In.push_back(rho);
  ux2In.push_back(ux);
  uy2In.push_back(uy);
  uz2In.push_back(uz);
  t2In.push_back(temp);

  while (param.read_optional("rho", rho)) {
    param.read_var("ux", ux);
    param.read_var("uy", uy);
    param.read_var("uz", uz);
    param.read_var("T", temp);
    rho2In.push_back(rho);
    ux2In.push_back(ux);
    uy2In.push_back(uy);
    uz2In.push_back(uz);
    t2In.push_back(temp);
  }

  param.read_var("bx", b2In[0]);
  param.read_var("by", b2In[1]);
  param.read_var("bz", b2In[2]);

  vnInput = -1.0;
  param.read_optional("vn", vnInput);
}

void DiscontinuityBoundaryManager::convert_and_validate(
    const FluidInterface* fi, const amrex::Geometry& geom, int nDim,
    bool isFake2D, bool useHybridPIC, amrex::Real electronTemperatureCode,
    amrex::Real electronDensity0Code, amrex::Real electronGammaIn) {
  if (!active || isCodeUnits)
    return;

  if (!fi) {
    amrex::Abort("DISCONTINUITY: FluidInterface is null.");
  }

  const amrex::Real normLen = norm3(normal);
  if (!(normLen > 0.0)) {
    amrex::Abort("Invalid #DISCONTINUITY: normal vector (normalX, normalY, "
                 "normalZ) must be non-zero.");
  }
  for (int d = 0; d < 3; ++d)
    normal[d] /= normLen;

  if ((nDim < 3 || isFake2D) && std::abs(normal[2]) > 1.0e-12) {
    amrex::Abort("Invalid #DISCONTINUITY: normalZ must be 0 in 2D or fake-2D "
                 "simulations.");
  }

  const auto& uniformState = fi->get_uniform_state();
  const int nS = fi->get_nS();
  if (nS <= 0 || static_cast<int>(uniformState.size()) < 5 * nS + 3) {
    amrex::Abort("Invalid #DISCONTINUITY: #UNIFORMSTATE must be defined to "
                 "provide the pre-discontinuity background (State 1).");
  }
  if (static_cast<int>(fi->inflowState.size()) < nS) {
    amrex::Abort("Invalid #DISCONTINUITY: State 1 inflow state is not "
                 "initialized.");
  }

  isHybridPIC = useHybridPIC;
  useHybridElectronP = (useHybridPIC && electronTemperatureCode > 0.0);
  te0Code = electronTemperatureCode;
  ne0Code = electronDensity0Code;
  gammaE = electronGammaIn;

  // Conversion factors: rho [1/cc], u [km/s], T [K], B [T], t [s] -> code units
  const amrex::Real si2NoRho = fi->get_Si2NoRho();
  const amrex::Real rhoFactor = 1.0e6 * cProtonMassSI * si2NoRho;
  const amrex::Real vFactor = 1.0e3 * fi->get_Si2NoV();
  const amrex::Real uNormSI = fi->get_unorm_si();
  const amrex::Real tFactor =
      cBoltzmannSI / (cProtonMassSI * uNormSI * uNormSI);
  const amrex::Real si2NoB = fi->get_Si2NoB();
  const amrex::Real si2NoT = fi->get_Si2NoT();

  tStart = tStartSI * si2NoT;

  // Populate State 1 (pre-discontinuity) from fi->inflowState and uniformState
  // B1
  state1.resize(nS);
  for (int iS = 0; iS < nS; ++iS) {
    const amrex::Real mass_i =
        (fi->get_species_mass(iS) > 0.0) ? fi->get_species_mass(iS) : 1.0;
    const amrex::Real charge_i = fi->get_species_charge(iS);
    state1[iS].nDens = fi->inflowState[iS].nDens;
    state1[iS].ux = fi->inflowState[iS].ux;
    state1[iS].uy = fi->inflowState[iS].uy;
    state1[iS].uz = fi->inflowState[iS].uz;
    state1[iS].vth = fi->inflowState[iS].vth;
    state1[iS].mass = mass_i;
    state1[iS].charge = charge_i;
  }
  for (int d = 0; d < 3; ++d) {
    b1[d] = uniformState[5 * nS + d] * si2NoB;
    b2[d] = b2In[d] * si2NoB;
  }

  // Populate State 2 (post-discontinuity) from #DISCONTINUITY blocks
  const int nBlock = static_cast<int>(rho2In.size());
  if (nBlock == 0) {
    amrex::Abort("Invalid #DISCONTINUITY: at least one species State 2 block "
                 "must be specified.");
  }

  state2.resize(nS);
  for (int iS = 0; iS < nS; ++iS) {
    const int iB = std::min(iS, nBlock - 1);
    if (rho2In[iB] <= 0.0) {
      amrex::Abort("Invalid #DISCONTINUITY: State 2 density rho must be > 0.");
    }
    if (t2In[iB] < 0.0) {
      amrex::Abort(
          "Invalid #DISCONTINUITY: State 2 temperature T must be >= 0.");
    }
    const amrex::Real mass_i = state1[iS].mass;
    const amrex::Real charge_i = state1[iS].charge;
    const amrex::Real nCode = rho2In[iB] * rhoFactor;
    const amrex::Real tempCode = t2In[iB] * tFactor;

    state2[iS].nDens = nCode;
    state2[iS].ux = ux2In[iB] * vFactor;
    state2[iS].uy = uy2In[iB] * vFactor;
    state2[iS].uz = uz2In[iB] * vFactor;
    state2[iS].vth = (tempCode > 0.0) ? std::sqrt(tempCode / mass_i) : 0.0;
    state2[iS].mass = mass_i;
    state2[iS].charge = charge_i;
  }

  compute_bulk_velocity(state1, isHybridPIC, uBulk1);
  compute_bulk_velocity(state2, isHybridPIC, uBulk2);

  // First-contact corner projection dMin = min_{x in Domain}(n . x)
  const amrex::Real* probLo = geom.ProbLo();
  const amrex::Real* probHi = geom.ProbHi();
  dMin = 0.0;
  for (int d = 0; d < nDim; ++d) {
    if (isFake2D && d == 2)
      continue;
    dMin +=
        (normal[d] >= 0.0) ? (normal[d] * probLo[d]) : (normal[d] * probHi[d]);
  }

  // Validate Rankine-Hugoniot / MHD discontinuity jump conditions
  const amrex::Real b1Mag = norm3(b1);
  const amrex::Real b2Mag = norm3(b2);
  const amrex::Real bScale = std::max({ b1Mag, b2Mag, 1.0e-30 });
  const amrex::Real bn1 = dot3(b1, normal);
  const amrex::Real bn2 = dot3(b2, normal);

  // 1. Divergence-free normal magnetic field continuity: B1.n == B2.n
  if (std::abs(bn2 - bn1) > 1.0e-2 * bScale) {
    std::ostringstream oss;
    oss << "Invalid #DISCONTINUITY: div(B) = 0 requires B1.n == B2.n within 1% "
        << "of |B|.\n"
        << "  normal = (" << normal[0] << ", " << normal[1] << ", " << normal[2]
        << ")\n"
        << "  B1.n = " << bn1 / si2NoB << " [T], B2.n = " << bn2 / si2NoB
        << " [T] (|B2.n - B1.n| / max(|B1|,|B2|) = "
        << std::abs(bn2 - bn1) / bScale << ").";
    amrex::Abort(oss.str());
  }
  // Remove any residual sub-tolerance normal difference so div(B) = 0 is exact
  for (int d = 0; d < 3; ++d) {
    b2[d] -= (bn2 - bn1) * normal[d];
  }
  bn = bn1;

  const amrex::Real un1 = dot3(uBulk1, normal);
  const amrex::Real un2 = dot3(uBulk2, normal);
  const amrex::Real uScale =
      std::max({ norm3(uBulk1), norm3(uBulk2), state1[0].vth, 1.0e-30 });

  const amrex::Real rhoM1 = compute_mass_density(state1, isHybridPIC);
  invSqrtFourPiRho1 = (rhoM1 > 0.0) ? (1.0 / std::sqrt(fourPI * rhoM1)) : 0.0;

  if (type == DiscontinuityType::TD) {
    // Tangential discontinuity checks:
    // (a) B1.n == 0 and B2.n == 0
    if (std::abs(bn1) > 1.0e-2 * bScale) {
      std::ostringstream oss;
      oss << "Invalid #DISCONTINUITY (type 'td'): tangential discontinuity "
          << "requires B.n = 0.\n"
          << "  Got |B1.n| / |B| = " << std::abs(bn1) / bScale
          << " (B1.n = " << bn1 / si2NoB << " [T]).";
      amrex::Abort(oss.str());
    }
    // Project out tiny sub-tolerance Bn so B.n == 0 is exact for TD
    for (int d = 0; d < 3; ++d) {
      b1[d] -= bn1 * normal[d];
      b2[d] -= bn1 * normal[d];
    }
    bn = 0.0;

    // (b) Normal velocity continuity: (u2 - u1).n == 0
    if (std::abs(un2 - un1) > 1.0e-2 * uScale) {
      std::ostringstream oss;
      oss << "Invalid #DISCONTINUITY (type 'td'): tangential discontinuity "
          << "requires normal velocity continuity u1.n == u2.n.\n"
          << "  u1.n = " << un1 / vFactor << " [km/s], u2.n = " << un2 / vFactor
          << " [km/s].";
      amrex::Abort(oss.str());
    }
    // Remove any sub-tolerance normal velocity difference so u_n == V_n is
    // exact
    const amrex::Real dUn = un2 - un1;
    for (int d = 0; d < 3; ++d) {
      uBulk2[d] -= dUn * normal[d];
    }
    for (int iS = 0; iS < nS; ++iS) {
      const amrex::Real u1sN = state1[iS].ux * normal[0] +
                               state1[iS].uy * normal[1] +
                               state1[iS].uz * normal[2];
      const amrex::Real u2sN = state2[iS].ux * normal[0] +
                               state2[iS].uy * normal[1] +
                               state2[iS].uz * normal[2];
      const amrex::Real dUsN = u2sN - u1sN;
      state2[iS].ux -= dUsN * normal[0];
      state2[iS].uy -= dUsN * normal[1];
      state2[iS].uz -= dUsN * normal[2];
    }

    // (c) Total pressure balance: P_tot,1 == P_tot,2
    const amrex::Real pTot1 = compute_total_pressure(
        state1, b1, isHybridPIC, useHybridElectronP, te0Code, ne0Code, gammaE);
    const amrex::Real pTot2 = compute_total_pressure(
        state2, b2, isHybridPIC, useHybridElectronP, te0Code, ne0Code, gammaE);
    const amrex::Real pScale = std::max({ pTot1, pTot2, 1.0e-30 });
    if (std::abs(pTot2 - pTot1) > 1.5e-2 * pScale) {
      const amrex::Real pMag1 = dot3(b1, b1) / (2.0 * fourPI);
      const amrex::Real pMag2 = dot3(b2, b2) / (2.0 * fourPI);
      const amrex::Real pPe2 = compute_hybrid_pe(state2, useHybridElectronP,
                                                 te0Code, ne0Code, gammaE);
      const amrex::Real pIonTarget2 = pTot1 - pMag2 - pPe2;
      const amrex::Real tTargetK =
          (pIonTarget2 > 0.0 && state2[0].nDens > 0.0)
              ? (pIonTarget2 / (state2[0].nDens * tFactor))
              : -1.0;
      std::ostringstream oss;
      oss << "Invalid #DISCONTINUITY (type 'td'): total pressure is not "
          << "balanced across the tangential discontinuity "
          << "(|P_tot,2 - P_tot,1| / P_tot = "
          << std::abs(pTot2 - pTot1) / pScale << " > 1.5%).\n"
          << "  State 1: P_tot,1 = " << pTot1 << " (P_mag,1 = " << pMag1
          << ") [code units]\n"
          << "  State 2: P_tot,2 = " << pTot2 << " (P_mag,2 = " << pMag2
          << ") [code units]";
      if (tTargetK > 0.0) {
        oss << "\n  Hint: for rho2 = " << rho2In[0]
            << " [/cc] and |B2| = " << norm3(b2In)
            << " [T], single-ion pressure balance requires T2 ~= " << tTargetK
            << " [K].";
      }
      amrex::Abort(oss.str());
    }

    if (vnInput > 0.0) {
      const amrex::Real vnOverride = vnInput * vFactor;
      if (std::abs(vnOverride - un1) > 1.0e-2 * uScale) {
        std::ostringstream oss;
        oss << "Invalid #DISCONTINUITY (type 'td'): override vn = " << vnInput
            << " [km/s] is inconsistent with the tangential discontinuity "
            << "normal flow speed u1.n = " << un1 / vFactor
            << " [km/s] (required by Faraday's law dB/dt + curl(E) = 0).";
        amrex::Abort(oss.str());
      }
    }
    vn = un1;

  } else if (type == DiscontinuityType::RD) {
    // Rotational discontinuity checks:
    // (a) Non-zero normal magnetic field B_n != 0
    if (std::abs(bn) <= 1.0e-6 * bScale) {
      amrex::Abort("Invalid #DISCONTINUITY (type 'rd'): rotational "
                   "discontinuity requires a non-zero normal magnetic field "
                   "(|B.n| > 1e-6 |B|).");
    }

    // (b) Constant magnetic field magnitude |B1| == |B2|
    if (std::abs(b2Mag - b1Mag) > 1.0e-2 * bScale) {
      std::ostringstream oss;
      oss << "Invalid #DISCONTINUITY (type 'rd'): rotational discontinuity "
          << "requires constant magnetic field magnitude |B1| == |B2|.\n"
          << "  |B1| = " << b1Mag / si2NoB << " [T], |B2| = " << b2Mag / si2NoB
          << " [T] (||B2| - |B1|| / |B| = " << std::abs(b2Mag - b1Mag) / bScale
          << ").";
      amrex::Abort(oss.str());
    }

    // (c) Constant density and temperature across RD for each species
    for (int iS = 0; iS < nS; ++iS) {
      if (isHybridPIC && state1[iS].charge < 0.0)
        continue;
      const amrex::Real nScale = std::max(state1[iS].nDens, 1.0e-30);
      if (std::abs(state2[iS].nDens - state1[iS].nDens) > 1.0e-2 * nScale) {
        std::ostringstream oss;
        oss << "Invalid #DISCONTINUITY (type 'rd'): rotational discontinuity "
            << "requires constant density n2 == n1 for species " << iS
            << " (got n1 = " << state1[iS].nDens / rhoFactor
            << " [/cc], n2 = " << state2[iS].nDens / rhoFactor << " [/cc]).";
        amrex::Abort(oss.str());
      }
      const amrex::Real vthScale = std::max(state1[iS].vth, 1.0e-30);
      if (std::abs(state2[iS].vth - state1[iS].vth) > 1.0e-2 * vthScale) {
        amrex::Abort("Invalid #DISCONTINUITY (type 'rd'): rotational "
                     "discontinuity requires constant temperature T2 == T1 for "
                     "each species.");
      }
    }

    // (d) Walen relation: du = u2 - u1 == sigma * (B2 - B1) / sqrt(4*pi*rho_m1)
    amrex::Real dB[3] = { b2[0] - b1[0], b2[1] - b1[1], b2[2] - b1[2] };
    amrex::Real du[3] = { uBulk2[0] - uBulk1[0], uBulk2[1] - uBulk1[1],
                          uBulk2[2] - uBulk1[2] };
    amrex::Real dVa[3] = { dB[0] * invSqrtFourPiRho1, dB[1] * invSqrtFourPiRho1,
                           dB[2] * invSqrtFourPiRho1 };
    const amrex::Real dVaMag = norm3(dVa);

    amrex::Real errPlus = 0.0, errMinus = 0.0;
    for (int d = 0; d < 3; ++d) {
      const amrex::Real dp = du[d] - dVa[d];
      const amrex::Real dm = du[d] + dVa[d];
      errPlus += dp * dp;
      errMinus += dm * dm;
    }
    errPlus = std::sqrt(errPlus);
    errMinus = std::sqrt(errMinus);

    if (errPlus <= errMinus) {
      walenSign = 1.0;
    } else {
      walenSign = -1.0;
    }
    const amrex::Real walenErr = std::min(errPlus, errMinus);
    const amrex::Real walenScale = std::max({ dVaMag, uScale, 1.0e-30 });
    if (walenErr > 1.5e-2 * walenScale) {
      std::ostringstream oss;
      oss << "Invalid #DISCONTINUITY (type 'rd'): State 2 velocity violates "
          << "the Walen relation u2 - u1 = +/- (B2 - B1) / sqrt(4*pi*rho_m).\n"
          << "  Got du = (" << du[0] / vFactor << ", " << du[1] / vFactor
          << ", " << du[2] / vFactor << ") [km/s]\n"
          << "  Expected du = +/- (" << dVa[0] / vFactor << ", "
          << dVa[1] / vFactor << ", " << dVa[2] / vFactor << ") [km/s].";
      amrex::Abort(oss.str());
    }

    // Validate per-species Walen velocity jump so multi-species opposing
    // offsets cannot cancel in uBulk while violating species continuity.
    for (int iS = 0; iS < nS; ++iS) {
      if (isHybridPIC && state1[iS].charge < 0.0)
        continue;
      const amrex::Real duS[3] = { state2[iS].ux - state1[iS].ux,
                                   state2[iS].uy - state1[iS].uy,
                                   state2[iS].uz - state1[iS].uz };
      amrex::Real errS = 0.0;
      for (int d = 0; d < 3; ++d) {
        const amrex::Real diff = duS[d] - walenSign * dVa[d];
        errS += diff * diff;
      }
      errS = std::sqrt(errS);
      if (errS > 1.5e-2 * walenScale) {
        std::ostringstream oss;
        oss << "Invalid #DISCONTINUITY (type 'rd'): species " << iS
            << " State 2 velocity violates the Walen relation u2 - u1 = "
            << (walenSign > 0.0 ? "+" : "-")
            << " (B2 - B1) / sqrt(4*pi*rho_m).\n"
            << "  Got du_s = (" << duS[0] / vFactor << ", " << duS[1] / vFactor
            << ", " << duS[2] / vFactor << ") [km/s].";
        amrex::Abort(oss.str());
      }
    }

    // Build orthonormal basis (e1, e2) in the plane perpendicular to normal
    amrex::Real b1tVec[3] = { b1[0] - bn * normal[0], b1[1] - bn * normal[1],
                              b1[2] - bn * normal[2] };
    amrex::Real b2tVec[3] = { b2[0] - bn * normal[0], b2[1] - bn * normal[1],
                              b2[2] - bn * normal[2] };
    bt1 = norm3(b1tVec);
    bt2 = norm3(b2tVec);

    if (bt1 > 1.0e-12 * bScale) {
      for (int d = 0; d < 3; ++d)
        e1[d] = b1tVec[d] / bt1;
    } else {
      // Pick arbitrary unit vector perpendicular to normal
      amrex::Real ref[3] = { 0.0, 0.0, 0.0 };
      const int minAxis =
          (std::abs(normal[0]) <= std::abs(normal[1]) &&
           std::abs(normal[0]) <= std::abs(normal[2]))
              ? 0
              : ((std::abs(normal[1]) <= std::abs(normal[2])) ? 1 : 2);
      ref[minAxis] = 1.0;
      cross3(normal, ref, e1);
      const amrex::Real e1Len = norm3(e1);
      for (int d = 0; d < 3; ++d)
        e1[d] /= e1Len;
    }
    cross3(normal, e1, e2);

    theta1 = 0.0;
    const amrex::Real c2 = dot3(b2tVec, e1);
    const amrex::Real s2 = dot3(b2tVec, e2);
    const amrex::Real theta2 =
        (bt2 > 1.0e-12 * bScale) ? std::atan2(s2, c2) : 0.0;

    // Shortest angular displacement in (-pi, pi]
    amrex::Real dTh = theta2 - theta1;
    while (dTh <= -dPI)
      dTh += 2.0 * dPI;
    while (dTh > dPI)
      dTh -= 2.0 * dPI;

    // For a 180-degree rotation (|dTh| ~= pi), use rdRotationSign (+pi vs -pi).
    // For general angles, if the user explicitly requested 'rd-', traverse the
    // complementary arc when dTh > 0 (or force negative rotation).
    if (std::abs(std::abs(dTh) - dPI) < 1.0e-3) {
      dTheta = static_cast<amrex::Real>(rdRotationSign) * dPI;
    } else if (rdRotationSign < 0 && dTh > 0.0) {
      dTheta = dTh - 2.0 * dPI;
    } else {
      dTheta = dTh;
    }

    const amrex::Real vAn = bn * invSqrtFourPiRho1;
    const amrex::Real vnExpected = un1 - walenSign * vAn;
    if (vnInput > 0.0) {
      const amrex::Real vnOverride = vnInput * vFactor;
      if (std::abs(vnOverride - vnExpected) > 1.5e-2 * walenScale) {
        std::ostringstream oss;
        oss << "Invalid #DISCONTINUITY (type 'rd'): override vn = " << vnInput
            << " [km/s] is inconsistent with the Walen-frame propagation "
            << "speed u1.n -/+ v_An = " << vnExpected / vFactor
            << " [km/s] (required by Faraday's law dB/dt + curl(E) = 0).";
        amrex::Abort(oss.str());
      }
    }
    vn = vnExpected;

  } else {
    // DiscontinuityType::General
    vn = (vnInput > 0.0) ? (vnInput * vFactor) : un1;
  }

  if (!(vn > 0.0)) {
    std::ostringstream oss;
    oss << "Invalid #DISCONTINUITY: normal propagation speed V_n = "
        << vn / vFactor
        << " [km/s] must be > 0 so the discontinuity enters the domain.\n"
        << "  Check the sign of (normalX, normalY, normalZ) or specify a "
        << "positive override 'vn' at the end of #DISCONTINUITY.";
    amrex::Abort(oss.str());
  }

  isCodeUnits = true;

  amrex::Print() << "  #DISCONTINUITY (" << typeStr << "):\n"
                 << "    normal = (" << normal[0] << ", " << normal[1] << ", "
                 << normal[2] << "), dMin = " << dMin << ", width = " << width
                 << " [code]\n"
                 << "    tStart = " << tStartSI << " [s] (" << tStart
                 << " [code]), V_n = " << vn / vFactor << " [km/s] (" << vn
                 << " [code])\n"
                 << "    B1 = (" << b1[0] << ", " << b1[1] << ", " << b1[2]
                 << "), B2 = (" << b2[0] << ", " << b2[1] << ", " << b2[2]
                 << ") [code]\n";
}

amrex::Real DiscontinuityBoundaryManager::phase_distance(const amrex::Real* pos,
                                                         amrex::Real t) const {
  const amrex::Real ndotx =
      normal[0] * pos[0] + normal[1] * pos[1] + normal[2] * pos[2];
  const amrex::Real pad = (width > 0.0) ? (3.0 * width) : 0.0;
  const amrex::Real dtActive = (t > tStart) ? (t - tStart) : 0.0;
  return (ndotx - dMin) + pad - vn * dtActive;
}

amrex::Real DiscontinuityBoundaryManager::alpha_at(const amrex::Real* pos,
                                                   amrex::Real t) const {
  if (!active || t <= tStart)
    return 0.0;
  const amrex::Real s = phase_distance(pos, t);
  if (width <= 0.0) {
    return (s < 0.0) ? 1.0 : 0.0;
  }
  const amrex::Real threeW = 3.0 * width;
  if (s >= threeW)
    return 0.0;
  if (s <= -threeW)
    return 1.0;
  static const amrex::Real invTanh3 = 1.0 / std::tanh(3.0);
  const amrex::Real a = 0.5 * (1.0 - std::tanh(s / width) * invTanh3);
  return std::clamp(a, 0.0, 1.0);
}

void DiscontinuityBoundaryManager::eval_B(const amrex::Real* pos, amrex::Real t,
                                          amrex::Real bOut[3]) const {
  const amrex::Real alpha = alpha_at(pos, t);
  if (alpha <= 0.0) {
    for (int d = 0; d < 3; ++d)
      bOut[d] = b1[d];
    return;
  }
  if (alpha >= 1.0) {
    for (int d = 0; d < 3; ++d)
      bOut[d] = b2[d];
    return;
  }

  if (type == DiscontinuityType::RD) {
    const amrex::Real bt = (1.0 - alpha) * bt1 + alpha * bt2;
    const amrex::Real theta = theta1 + alpha * dTheta;
    const amrex::Real c = std::cos(theta);
    const amrex::Real s = std::sin(theta);
    for (int d = 0; d < 3; ++d) {
      bOut[d] = bn * normal[d] + bt * (c * e1[d] + s * e2[d]);
    }
  } else {
    for (int d = 0; d < 3; ++d) {
      bOut[d] = (1.0 - alpha) * b1[d] + alpha * b2[d];
    }
  }
}

void DiscontinuityBoundaryManager::eval_u_bulk(const amrex::Real* pos,
                                               amrex::Real t,
                                               amrex::Real uOut[3]) const {
  const amrex::Real alpha = alpha_at(pos, t);
  if (alpha <= 0.0) {
    for (int d = 0; d < 3; ++d)
      uOut[d] = uBulk1[d];
    return;
  }
  if (alpha >= 1.0) {
    for (int d = 0; d < 3; ++d)
      uOut[d] = uBulk2[d];
    return;
  }

  if (type == DiscontinuityType::RD) {
    amrex::Real bCur[3];
    eval_B(pos, t, bCur);
    for (int d = 0; d < 3; ++d) {
      const amrex::Real dVaCur =
          walenSign * (bCur[d] - b1[d]) * invSqrtFourPiRho1;
      const amrex::Real dVaTot =
          walenSign * (b2[d] - b1[d]) * invSqrtFourPiRho1;
      const amrex::Real resid = (uBulk2[d] - uBulk1[d]) - dVaTot;
      uOut[d] = uBulk1[d] + dVaCur + alpha * resid;
    }
  } else {
    for (int d = 0; d < 3; ++d) {
      uOut[d] = (1.0 - alpha) * uBulk1[d] + alpha * uBulk2[d];
    }
  }
}

void DiscontinuityBoundaryManager::eval_E(const amrex::Real* pos, amrex::Real t,
                                          amrex::Real eOut[3]) const {
  amrex::Real bCur[3];
  amrex::Real uCur[3];
  eval_B(pos, t, bCur);
  eval_u_bulk(pos, t, uCur);
  // E = -u_bulk x B = B x u_bulk
  cross3(bCur, uCur, eOut);
}

void DiscontinuityBoundaryManager::eval_species(int iS, const amrex::Real* pos,
                                                amrex::Real t,
                                                amrex::Real& nOut,
                                                amrex::Real uOut[3],
                                                amrex::Real& vthOut) const {
  if (iS < 0 || iS >= static_cast<int>(state1.size())) {
    nOut = 0.0;
    uOut[0] = 0.0;
    uOut[1] = 0.0;
    uOut[2] = 0.0;
    vthOut = 0.0;
    return;
  }

  const auto& s1 = state1[iS];
  const auto& s2 = state2[iS];
  const amrex::Real alpha = alpha_at(pos, t);

  if (alpha <= 0.0) {
    nOut = s1.nDens;
    uOut[0] = s1.ux;
    uOut[1] = s1.uy;
    uOut[2] = s1.uz;
    vthOut = s1.vth;
    return;
  }
  if (alpha >= 1.0) {
    nOut = s2.nDens;
    uOut[0] = s2.ux;
    uOut[1] = s2.uy;
    uOut[2] = s2.uz;
    vthOut = s2.vth;
    return;
  }

  nOut = (1.0 - alpha) * s1.nDens + alpha * s2.nDens;

  if (type == DiscontinuityType::RD) {
    amrex::Real bCur[3];
    eval_B(pos, t, bCur);
    const amrex::Real s1u[3] = { s1.ux, s1.uy, s1.uz };
    const amrex::Real s2u[3] = { s2.ux, s2.uy, s2.uz };
    for (int d = 0; d < 3; ++d) {
      const amrex::Real dVaCur =
          walenSign * (bCur[d] - b1[d]) * invSqrtFourPiRho1;
      const amrex::Real dVaTot =
          walenSign * (b2[d] - b1[d]) * invSqrtFourPiRho1;
      const amrex::Real resid = (s2u[d] - s1u[d]) - dVaTot;
      uOut[d] = s1u[d] + dVaCur + alpha * resid;
    }
  } else {
    uOut[0] = (1.0 - alpha) * s1.ux + alpha * s2.ux;
    uOut[1] = (1.0 - alpha) * s1.uy + alpha * s2.uy;
    uOut[2] = (1.0 - alpha) * s1.uz + alpha * s2.uz;
  }

  // Base thermal pressure interpolation P_s(s) = (1 - alpha)*P_1,s +
  // alpha*P_2,s
  const amrex::Real p1s = s1.nDens * s1.mass * s1.vth * s1.vth;
  const amrex::Real p2s = s2.nDens * s2.mass * s2.vth * s2.vth;
  amrex::Real pCurS = (1.0 - alpha) * p1s + alpha * p2s;

  // For a tangential discontinuity with finite width, linear interpolation of B
  // across a field reversal or shear reduces |B(s)|^2 inside the layer by
  //   dP_B(s) = ((1 - alpha)*|B1|^2 + alpha*|B2|^2 - |B(s)|^2) / (8*pi) >= 0.
  // Compensate this magnetic pressure dip in the kinetic thermal pressure so
  // total pressure P_tot(s) is constant throughout the finite-width layer.
  if (type == DiscontinuityType::TD) {
    amrex::Real bCur[3];
    eval_B(pos, t, bCur);
    const amrex::Real b1Sq = dot3(b1, b1);
    const amrex::Real b2Sq = dot3(b2, b2);
    const amrex::Real bCurSq = dot3(bCur, bCur);
    const amrex::Real dPB = std::max(
        0.0, ((1.0 - alpha) * b1Sq + alpha * b2Sq - bCurSq) / (2.0 * fourPI));

    // Also account for any non-linearity of adiabatic hybrid P_e(s) if gammaE
    // != 1
    amrex::Real dPPe = 0.0;
    if (useHybridElectronP && te0Code > 0.0 && ne0Code > 0.0) {
      const amrex::Real pe1 = compute_hybrid_pe(state1, useHybridElectronP,
                                                te0Code, ne0Code, gammaE);
      const amrex::Real pe2 = compute_hybrid_pe(state2, useHybridElectronP,
                                                te0Code, ne0Code, gammaE);
      amrex::Real neCur = 0.0;
      for (size_t k = 0; k < state1.size(); ++k) {
        if (state1[k].charge > 0.0) {
          const amrex::Real nk =
              (1.0 - alpha) * state1[k].nDens + alpha * state2[k].nDens;
          neCur += state1[k].charge * nk;
        }
      }
      const amrex::Real peCur =
          (neCur > 0.0)
              ? (ne0Code * te0Code * std::pow(neCur / ne0Code, gammaE))
              : 0.0;
      dPPe = (1.0 - alpha) * pe1 + alpha * pe2 - peCur;
    }

    const amrex::Real dPTot = std::max(0.0, dPB + dPPe);
    if (dPTot > 0.0) {
      // Distribute thermal pressure compensation across kinetic species
      // weighted by local number density (or thermal pressure if non-zero).
      amrex::Real nSum = 0.0;
      for (size_t k = 0; k < state1.size(); ++k) {
        if (isHybridPIC && state1[k].charge < 0.0)
          continue;
        const amrex::Real nk =
            (1.0 - alpha) * state1[k].nDens + alpha * state2[k].nDens;
        if (nk > 0.0)
          nSum += nk;
      }
      if (nSum > 0.0 && nOut > 0.0) {
        pCurS += dPTot * (nOut / nSum);
      }
    }
  }

  if (nOut > 0.0 && s1.mass > 0.0 && pCurS > 0.0) {
    vthOut = std::sqrt(pCurS / (nOut * s1.mass));
  } else {
    vthOut = (1.0 - alpha) * s1.vth + alpha * s2.vth;
  }
}
