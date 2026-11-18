#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>

#include <AMReX_Utility.H>

#include "BitArray.h"
#include "GpuBackend.h"
#include "TestParticles.h"
#include "Utility.h"

using namespace amrex;

namespace {

ParticlesInfo make_test_particles_info(const ParticlesInfo& picInfo) {
  ParticlesInfo info;
  info.nPartPerCell = IntVect(AMREX_D_DECL(2, 2, 2));
  info.isParticleLocationRandom = picInfo.isParticleLocationRandom;
  info.pBCs = picInfo.pBCs;
  info.pBCsSet = picInfo.pBCsSet;
  return info;
}

inline std::uint64_t pack_tp_key(int cpu, int id) noexcept {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(cpu)) << 32) |
         static_cast<std::uint32_t>(id);
}

std::string default_output_string_from_compile_size() {
  switch (ptRecordSize) {
    case 7:
      return "P";
    case 10:
      return "PB";
    case 13:
      return "PBE";
    case 19:
      return "PBED";
    case 22:
      return "PBEG";
    default:
      Abort("Error: unsupported compile-time ptRecordSize = " +
            std::to_string(ptRecordSize));
      return "PBED";
  }
}

void compute_magnetic_drifts(Real mass, Real charge, Real gamma,
                             const Real bp[3], const Real vRec[3],
                             const Real gradB[3][3], Real vGradB[3],
                             Real vCurv[3]) {
  vGradB[0] = 0.0;
  vGradB[1] = 0.0;
  vGradB[2] = 0.0;
  vCurv[0] = 0.0;
  vCurv[1] = 0.0;
  vCurv[2] = 0.0;

  const Real B2 = bp[0] * bp[0] + bp[1] * bp[1] + bp[2] * bp[2];
  if (B2 <= 1e-30 || std::abs(charge) <= 1e-30)
    return;

  const Real invB = 1.0 / std::sqrt(B2);
  const Real bx = bp[0] * invB;
  const Real by = bp[1] * invB;
  const Real bz = bp[2] * invB;

  const Real vpar = vRec[0] * bx + vRec[1] * by + vRec[2] * bz;
  const Real v2 = vRec[0] * vRec[0] + vRec[1] * vRec[1] + vRec[2] * vRec[2];
  const Real vperp2 = std::max(Real(0.0), v2 - vpar * vpar);
  const Real mEff = mass * gamma;

  // Grad |B| vector: grad(|B|)_j = sum_i b_i * d(B_i)/d(x_j)
  const Real gradBx = bx * gradB[0][0] + by * gradB[1][0] + bz * gradB[2][0];
  const Real gradBy = bx * gradB[0][1] + by * gradB[1][1] + bz * gradB[2][1];
  const Real gradBz = bx * gradB[0][2] + by * gradB[1][2] + bz * gradB[2][2];

  // b x grad(|B|)
  const Real b_cross_gradB_x = by * gradBz - bz * gradBy;
  const Real b_cross_gradB_y = bz * gradBx - bx * gradBz;
  const Real b_cross_gradB_z = bx * gradBy - by * gradBx;

  // v_gradB = (mEff * vperp2) / (2 * q * B^2) * (b x grad|B|)
  const Real coefGradB = (mEff * vperp2) / (2.0 * charge * B2);
  vGradB[0] = coefGradB * b_cross_gradB_x;
  vGradB[1] = coefGradB * b_cross_gradB_y;
  vGradB[2] = coefGradB * b_cross_gradB_z;

  // w = (b . grad) B -> w_i = sum_j J_{ij} * b_j
  const Real wx = gradB[0][0] * bx + gradB[0][1] * by + gradB[0][2] * bz;
  const Real wy = gradB[1][0] * bx + gradB[1][1] * by + gradB[1][2] * bz;
  const Real wz = gradB[2][0] * bx + gradB[2][1] * by + gradB[2][2] * bz;

  // b x w
  const Real b_cross_w_x = by * wz - bz * wy;
  const Real b_cross_w_y = bz * wx - bx * wz;
  const Real b_cross_w_z = bx * wy - by * wx;

  // v_curv = (mEff * vpar^2) / (q * B^2) * (b x w)
  const Real coefCurv = (mEff * vpar * vpar) / (charge * B2);
  vCurv[0] = coefCurv * b_cross_w_x;
  vCurv[1] = coefCurv * b_cross_w_y;
  vCurv[2] = coefCurv * b_cross_w_z;
}

inline void interpolate_jacobian_matrix(const Array4<const Real>& jacArr,
                                        const IntVect& loIdx,
                                        const Real coef[2][2][2],
                                        const Dim3& lo, const Dim3& hi,
                                        Real gradB[3][3]) {
  for (int k = lo.z; k <= hi.z; ++k)
    for (int j = lo.y; j <= hi.y; ++j)
      for (int i = lo.x; i <= hi.x; ++i) {
        const IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + i, loIdx[iy_] + j,
                                           loIdx[iz_] + k) };
        const Real c = coef[i][j][k];
        for (int iRow = 0; iRow < 3; ++iRow) {
          for (int jCol = 0; jCol < 3; ++jCol) {
            gradB[iRow][jCol] += jacArr(ijk, iRow * 3 + jCol) * c;
          }
        }
      }
}

void compute_nodal_grad_b(const Array4<const Real>& nodeBArr,
                          const IntVect& loIdx, const RealVect& dShift,
                          const Real* invDx, Real gradB[3][3]) {
  Real b[3][2][2][2];
  const int kMax = (nDim > 2) ? 2 : 1;
  for (int k = 0; k < kMax; ++k)
    for (int j = 0; j < 2; ++j)
      for (int i = 0; i < 2; ++i) {
        const IntVect ijk = { AMREX_D_DECL(loIdx[ix_] + i, loIdx[iy_] + j,
                                           loIdx[iz_] + k) };
        for (int iDim = 0; iDim < nDim3; ++iDim) {
          b[iDim][i][j][k] = nodeBArr(ijk, iDim);
          if (nDim == 2)
            b[iDim][i][j][1] = b[iDim][i][j][0];
        }
      }

  const Real sx = dShift[ix_];
  const Real sy = dShift[iy_];
  const Real sz = (nDim > 2) ? dShift[iz_] : 0.0;

  for (int iDim = 0; iDim < nDim3; ++iDim) {
    gradB[iDim][ix_] = (((b[iDim][1][0][0] - b[iDim][0][0][0]) * (1.0 - sy) +
                         (b[iDim][1][1][0] - b[iDim][0][1][0]) * sy) *
                            (1.0 - sz) +
                        ((b[iDim][1][0][1] - b[iDim][0][0][1]) * (1.0 - sy) +
                         (b[iDim][1][1][1] - b[iDim][0][1][1]) * sy) *
                            sz) *
                       invDx[ix_];

    gradB[iDim][iy_] = (((b[iDim][0][1][0] - b[iDim][0][0][0]) * (1.0 - sx) +
                         (b[iDim][1][1][0] - b[iDim][1][0][0]) * sx) *
                            (1.0 - sz) +
                        ((b[iDim][0][1][1] - b[iDim][0][0][1]) * (1.0 - sx) +
                         (b[iDim][1][1][1] - b[iDim][1][0][1]) * sx) *
                            sz) *
                       invDx[iy_];

#if (AMREX_SPACEDIM == 3)
    gradB[iDim][iz_] = (((b[iDim][0][0][1] - b[iDim][0][0][0]) * (1.0 - sx) +
                         (b[iDim][1][0][1] - b[iDim][1][0][0]) * sx) *
                            (1.0 - sy) +
                        ((b[iDim][0][1][1] - b[iDim][0][1][0]) * (1.0 - sx) +
                         (b[iDim][1][1][1] - b[iDim][1][1][0]) * sx) *
                            sy) *
                       invDx[iz_];
#else
    gradB[iDim][iz_] = 0.0;
#endif
  }
}

inline void reflect_record_at_walls(const BoxBC<ParticleBC::Type>& bc,
                                    const Real* plo, const Real* phi,
                                    Real& xRec, Real& yRec, Real& zRec,
                                    Real& uRec, Real& vRec, Real& wRec) {
  Real* pos[3] = { &xRec, &yRec, &zRec };
  Real* vel[3] = { &uRec, &vRec, &wRec };
  for (int d = 0; d < nDim; ++d) {
    if (bc.lo[d] == ParticleBC::reflect && *pos[d] < plo[d]) {
      *pos[d] = Real(2.0) * plo[d] - *pos[d];
      *vel[d] = -*vel[d];
    } else if (bc.hi[d] == ParticleBC::reflect && *pos[d] > phi[d]) {
      *pos[d] = Real(2.0) * phi[d] - *pos[d];
      *vel[d] = -*vel[d];
    }
  }
}

void append_particle_record(std::vector<Real>& dst, const TPOutputSpec& spec,
                            Real mass, Real charge, bool isRelativistic,
                            Real tRec, Real xRec, Real yRec, Real zRec,
                            Real uRec, Real vRec, Real wRec, const Real* bp,
                            const Real* ep, const Real (*gradB)[3]) {
  dst.push_back(tRec);
  dst.push_back(xRec);
  dst.push_back(yRec);
  dst.push_back(zRec);
  dst.push_back(uRec);
  dst.push_back(vRec);
  dst.push_back(wRec);

  if (spec.groups.empty())
    return;

  const Real vRecArr[3] = { uRec, vRec, wRec };
  const Real v2 = uRec * uRec + vRec * vRec + wRec * wRec;
  Real gamma = 1.0;
  if (isRelativistic) {
    gamma = (v2 < 1.0) ? (Real(1.0) / std::sqrt(Real(1.0) - v2)) : Real(1e10);
  }

  const Real bSafe[3] = { bp ? bp[0] : 0.0, bp ? bp[1] : 0.0,
                          bp ? bp[2] : 0.0 };
  const Real eSafe[3] = { ep ? ep[0] : 0.0, ep ? ep[1] : 0.0,
                          ep ? ep[2] : 0.0 };
  const Real B2 =
      bSafe[0] * bSafe[0] + bSafe[1] * bSafe[1] + bSafe[2] * bSafe[2];

  Real vGradB[3] = { 0.0, 0.0, 0.0 };
  Real vCurv[3] = { 0.0, 0.0, 0.0 };
  bool driftsComputed = false;

  for (TPOutputGroup g : spec.groups) {
    switch (g) {
      case TPOutputGroup::B:
        dst.push_back(bSafe[0]);
        dst.push_back(bSafe[1]);
        dst.push_back(bSafe[2]);
        break;
      case TPOutputGroup::E:
        dst.push_back(eSafe[0]);
        dst.push_back(eSafe[1]);
        dst.push_back(eSafe[2]);
        break;
      case TPOutputGroup::Drifts:
        if (!driftsComputed) {
          if (gradB) {
            compute_magnetic_drifts(mass, charge, gamma, bSafe, vRecArr, gradB,
                                    vGradB, vCurv);
          }
          driftsComputed = true;
        }
        dst.push_back(vGradB[0]);
        dst.push_back(vGradB[1]);
        dst.push_back(vGradB[2]);
        dst.push_back(vCurv[0]);
        dst.push_back(vCurv[1]);
        dst.push_back(vCurv[2]);
        break;
      case TPOutputGroup::GradB:
        for (int iRow = 0; iRow < 3; ++iRow) {
          for (int jCol = 0; jCol < 3; ++jCol) {
            dst.push_back(gradB ? gradB[iRow][jCol] : Real(0.0));
          }
        }
        break;
      case TPOutputGroup::VExB:
        if (B2 > 1e-30) {
          const Real invB2 = Real(1.0) / B2;
          dst.push_back((eSafe[1] * bSafe[2] - eSafe[2] * bSafe[1]) * invB2);
          dst.push_back((eSafe[2] * bSafe[0] - eSafe[0] * bSafe[2]) * invB2);
          dst.push_back((eSafe[0] * bSafe[1] - eSafe[1] * bSafe[0]) * invB2);
        } else {
          dst.push_back(0.0);
          dst.push_back(0.0);
          dst.push_back(0.0);
        }
        break;
      case TPOutputGroup::Mu: {
        Real muVal = 0.0;
        if (B2 > 1e-30) {
          const Real bMag = std::sqrt(B2);
          const Real invB = Real(1.0) / bMag;
          const Real vpar =
              (uRec * bSafe[0] + vRec * bSafe[1] + wRec * bSafe[2]) * invB;
          const Real vperp2 = std::max(Real(0.0), v2 - vpar * vpar);
          muVal = Real(0.5) * mass * gamma * gamma * vperp2 * invB;
        }
        dst.push_back(muVal);
        break;
      }
      case TPOutputGroup::Pitch: {
        Real pitchCos = 0.0;
        if (B2 > 1e-30 && v2 > 1e-30) {
          const Real denom = std::sqrt(v2 * B2);
          const Real raw =
              (uRec * bSafe[0] + vRec * bSafe[1] + wRec * bSafe[2]) / denom;
          pitchCos = std::clamp(raw, Real(-1.0), Real(1.0));
        }
        dst.push_back(pitchCos);
        break;
      }
      case TPOutputGroup::Energy: {
        const Real eCode =
            isRelativistic ? (mass * gamma * gamma * v2 / (gamma + Real(1.0)))
                           : (Real(0.5) * mass * v2);
        dst.push_back(eCode);
        break;
      }
      case TPOutputGroup::Gamma:
        dst.push_back(gamma);
        break;
    }
  }
}

} // namespace

//==========================================================
TPOutputSpec TPOutputSpec::parse(std::string_view specStr) {
  TPOutputSpec out;
  out.colNames = { "time", "x", "y", "z", "vx", "vy", "vz" };
  out.recordSize = 7;
  out.needGradB = false;

  std::string raw = specStr.empty() ? default_output_string_from_compile_size()
                                    : std::string(specStr);

  bool seen[9] = {
    false, false, false, false, false, false, false, false, false
  };

  auto addGroup = [&](TPOutputGroup g) {
    const int idx = static_cast<int>(g);
    if (seen[idx])
      return;
    seen[idx] = true;
    out.groups.push_back(g);

    switch (g) {
      case TPOutputGroup::B:
        out.colNames.insert(out.colNames.end(), { "bx", "by", "bz" });
        out.recordSize += 3;
        break;
      case TPOutputGroup::E:
        out.colNames.insert(out.colNames.end(), { "ex", "ey", "ez" });
        out.recordSize += 3;
        break;
      case TPOutputGroup::Drifts:
        out.colNames.insert(out.colNames.end(),
                            { "v_gradb_x", "v_gradb_y", "v_gradb_z", "v_curv_x",
                              "v_curv_y", "v_curv_z" });
        out.recordSize += 6;
        out.needGradB = true;
        break;
      case TPOutputGroup::GradB:
        out.colNames.insert(out.colNames.end(),
                            { "dbxdx", "dbxdy", "dbxdz", "dbydx", "dbydy",
                              "dbydz", "dbzdx", "dbzdy", "dbzdz" });
        out.recordSize += 9;
        out.needGradB = true;
        break;
      case TPOutputGroup::VExB:
        out.colNames.insert(out.colNames.end(),
                            { "v_exb_x", "v_exb_y", "v_exb_z" });
        out.recordSize += 3;
        break;
      case TPOutputGroup::Mu:
        out.colNames.push_back("mu");
        out.recordSize += 1;
        break;
      case TPOutputGroup::Pitch:
        out.colNames.push_back("pitch");
        out.recordSize += 1;
        break;
      case TPOutputGroup::Energy:
        out.colNames.push_back("energy");
        out.recordSize += 1;
        break;
      case TPOutputGroup::Gamma:
        out.colNames.push_back("gamma");
        out.recordSize += 1;
        break;
    }
  };

  for (char& ch : raw) {
    if (ch == ',' || ch == '+' || ch == '\t' || ch == '\r' || ch == '\n')
      ch = ' ';
  }

  std::istringstream iss(raw);
  std::string rawTok;
  int nTokens = 0;
  while (iss >> rawTok) {
    ++nTokens;
    const std::string tok = normalize_string_token(rawTok);
    if (tok == "p") {
      // Base 7 columns already included.
    } else if (tok == "pb") {
      addGroup(TPOutputGroup::B);
    } else if (tok == "pbe") {
      addGroup(TPOutputGroup::B);
      addGroup(TPOutputGroup::E);
    } else if (tok == "pbed") {
      addGroup(TPOutputGroup::B);
      addGroup(TPOutputGroup::E);
      addGroup(TPOutputGroup::Drifts);
    } else if (tok == "pbeg") {
      addGroup(TPOutputGroup::B);
      addGroup(TPOutputGroup::E);
      addGroup(TPOutputGroup::GradB);
    } else if (tok == "b") {
      addGroup(TPOutputGroup::B);
    } else if (tok == "e") {
      addGroup(TPOutputGroup::E);
    } else if (tok == "drifts") {
      addGroup(TPOutputGroup::Drifts);
    } else if (tok == "gradb") {
      addGroup(TPOutputGroup::GradB);
    } else if (tok == "vexb") {
      addGroup(TPOutputGroup::VExB);
    } else if (tok == "mu") {
      addGroup(TPOutputGroup::Mu);
    } else if (tok == "pitch") {
      addGroup(TPOutputGroup::Pitch);
    } else if (tok == "energy") {
      addGroup(TPOutputGroup::Energy);
    } else if (tok == "gamma") {
      addGroup(TPOutputGroup::Gamma);
    } else {
      Abort("Error [TestParticles]: unknown #TPOUTPUT variable token '" +
            rawTok + "'.");
    }
  }

  if (nTokens == 0 && !specStr.empty()) {
    Abort("Error [TestParticles]: #TPOUTPUT outputVars cannot be empty.");
  }

  return out;
}

//==========================================================
TestParticles::TestParticles(Grid* gridIn, FluidInterface* const fluidIn,
                             TimeCtr* const tcIn, const int speciesID,
                             const Real charge, const Real mass, int id,
                             const ParticlesInfo& picInfo)
    : Particles(gridIn, fluidIn, tcIn, speciesID, charge, mass,
                make_test_particles_info(picInfo), PartMode::PIC) {
  gridID = id;

  gridName = std::string("FLEKS") + std::to_string(gridID);
  printPrefix = gridName + ": ";

  outputDir = component + "/plots/test_particles";

  nInitPart = 0;
  rebuild_col_scales();
}

//==========================================================
void TestParticles::set_IO_units(Real no2outLIn, Real no2outVIn, Real no2outMIn,
                                 Real no2outBIn, Real no2outEIn, Real no2siVIn,
                                 Real no2siBIn, const std::string& ioUnitIn) {
  no2outL = no2outLIn;
  no2outV = no2outVIn;
  no2outM = no2outMIn;
  no2outB = no2outBIn;
  no2outE = no2outEIn;
  no2siV = no2siVIn;
  no2siB = no2siBIn;
  ioUnit = normalize_string_token(ioUnitIn);

  if (ioUnit == "pic") {
    energyScale = 1.0;
    muScale = 1.0;
  } else if (ioUnit == "si") {
    energyScale = cProtonMassSI * no2siV * no2siV;
    muScale = (no2siB > 0.0) ? (energyScale / no2siB) : 1.0;
  } else {
    // Planet / planetary units: energy in [eV], magnetic moment in [eV/nT].
    energyScale = (cProtonMassSI * no2siV * no2siV) / cUnitChargeSI;
    muScale = (no2siB > 0.0) ? (energyScale / (no2siB * 1.0e9)) : 1.0;
  }

  rebuild_col_scales();
}

void TestParticles::rebuild_col_scales() {
  colScales.clear();
  colScales.reserve(outputSpec.recordSize);
  colScales.push_back(1.0);
  colScales.push_back(no2outL);
  colScales.push_back(no2outL);
  colScales.push_back(no2outL);
  colScales.push_back(no2outV);
  colScales.push_back(no2outV);
  colScales.push_back(no2outV);

  const Real no2outG = (no2outL > 0.0) ? (no2outB / no2outL) : Real(1.0);
  for (TPOutputGroup g : outputSpec.groups) {
    switch (g) {
      case TPOutputGroup::B:
        colScales.insert(colScales.end(), 3, no2outB);
        break;
      case TPOutputGroup::E:
        colScales.insert(colScales.end(), 3, no2outE);
        break;
      case TPOutputGroup::Drifts:
        colScales.insert(colScales.end(), 6, no2outV);
        break;
      case TPOutputGroup::GradB:
        colScales.insert(colScales.end(), 9, no2outG);
        break;
      case TPOutputGroup::VExB:
        colScales.insert(colScales.end(), 3, no2outV);
        break;
      case TPOutputGroup::Mu:
        colScales.push_back(muScale);
        break;
      case TPOutputGroup::Pitch:
        colScales.push_back(1.0);
        break;
      case TPOutputGroup::Energy:
        colScales.push_back(energyScale);
        break;
      case TPOutputGroup::Gamma:
        colScales.push_back(1.0);
        break;
    }
  }
}

//==========================================================
void TestParticles::interpolate_record_trajectory(
    Real xp, Real yp, Real zp, Real up, Real vp, Real wp, Real unp1, Real vnp1,
    Real wnp1, Real dtStep, Real dtElse, Real tNowSI, Real& tRec, Real& xRec,
    Real& yRec, Real& zRec, Real& uRec, Real& vRec, Real& wRec) const {
  if (dtSave > 0.0) {
    tRec = tNextSave;
    const Real dtSI = tc->get_dt_si();
    Real theta = 1.0;
    if (dtSI > 0.0) {
      theta = (tNextSave - (tNowSI - dtSI)) / dtSI;
      theta = std::max(Real(0.0), std::min(Real(1.0), theta));
    }
    xRec = xp + theta * unp1 * dtStep;
    yRec = yp + theta * vnp1 * dtStep;
    zRec = zp + theta * wnp1 * dtStep;
    uRec = (1.0 - theta) * up + theta * unp1;
    vRec = (1.0 - theta) * vp + theta * vnp1;
    wRec = (1.0 - theta) * wp + theta * wnp1;
  } else {
    tRec = tNowSI;
    xRec = xp + unp1 * dtElse;
    yRec = yp + vnp1 * dtElse;
    zRec = zp + wnp1 * dtElse;
    uRec = unp1;
    vRec = vnp1;
    wRec = wnp1;
  }
}

void TestParticles::commit_stages(const std::vector<SlotStage>& stages) {
  const std::size_t rSize = static_cast<std::size_t>(outputSpec.recordSize);
  for (const auto& stage : stages) {
    const std::size_t nRec = stage.metas.size();
    for (std::size_t k = 0; k < nRec; ++k) {
      const auto& m = stage.metas[k];
      auto& entry = trajMap_[pack_tp_key(m.cpu, m.id)];
      if (entry.data.empty()) {
        entry.cpu = m.cpu;
        entry.id = m.id;
        entry.qp = m.qp;
      }
      const Real* recPtr = &stage.values[k * rSize];
      entry.data.insert(entry.data.end(), recPtr, recPtr + rSize);
    }
  }
}

//==========================================================
void TestParticles::move_and_save_particles(
    int iLev, const MultiFab& nodeEMF, const MultiFab* nodeE0MF,
    const MultiFab& nodeBMF, const MultiFab* nodeB0MF, Real dt, Real dtNext,
    Real tNowSI, bool doSave, const MultiFab* nodeJacBMF) {
  if (is_neutral()) {
    move_and_save_neutrals(iLev, dt, tNowSI, doSave);
  } else {
    move_and_save_charged_particles(iLev, nodeEMF, nodeE0MF, nodeBMF, nodeB0MF,
                                    dt, dtNext, tNowSI, doSave, nodeJacBMF);
  }
}

//==========================================================
void TestParticles::move_and_save_charged_particles(
    int iLev, const MultiFab& nodeEMF, const MultiFab* nodeE0MF,
    const MultiFab& nodeBMF, const MultiFab* nodeB0MF, Real dt, Real dtNext,
    Real tNowSI, bool doSave, const MultiFab* nodeJacBMF) {
  timing_func("TestParticles::move_charged_particles");

  const Real dtLoc = 0.5 * (dt + dtNext);
  const Real qdto2mc = charge / mass * 0.5 * dt;
  const LevelGeomBox lb = make_level_geom_box(iLev);
  const BoxArray& ba = ParticleBoxArray(iLev);
  const bool hasB0 = (nodeB0MF != nullptr);
  const bool hasE0 = (nodeE0MF != nullptr);
  const bool hasJacB = (nodeJacBMF != nullptr);

  struct TPFabCtx {
    ActiveRegionBox ab;
    Array4<int const> status;
    Array4<Real const> EArr;
    Array4<Real const> E0Arr;
    Array4<Real const> BArr;
    Array4<Real const> B0Arr;
    Array4<Real const> jacBArr;
  };
  const auto fabTable =
      build_local_fab_table(cell_status(iLev), [&](const MFIter& mfi) {
        const Box& bx = cell_status(iLev)[mfi].box();
        return TPFabCtx{
          make_active_region_box(lb, ba[mfi.index()], bx),
          cell_status(iLev)[mfi].const_array(),
          nodeEMF[mfi].const_array(),
          hasE0 ? (*nodeE0MF)[mfi].const_array() : Array4<Real const>{},
          nodeBMF[mfi].const_array(),
          hasB0 ? (*nodeB0MF)[mfi].const_array() : Array4<Real const>{},
          hasJacB ? (*nodeJacBMF)[mfi].const_array() : Array4<Real const>{}
        };
      });

  std::vector<SlotStage> stages(doSave ? omp_slot_count() : 0);
  std::atomic<bool> superluminal{ false };

  parallel_for_cell_tiles(*this, iLev, [&](const auto& item, int /*ord*/) {
    const auto& ctx = fabTable[item.grid];
    auto& particles = item.ptile->GetArrayOfStructs();
    const Dim3 lo = init_dim3(0);
    const Dim3 hi = init_dim3(1);
    const Real ploTile[3] = { lb.plo[0], lb.plo[1], lb.plo[2] };
    const Real invDxTile[3] = { lb.invDx[0], lb.invDx[1], lb.invDx[2] };
    const PusherType pusherTile = pusherType;
    const bool isRelTile = isRelativistic;
    const Real qdto2mcTile = qdto2mc;
    const Real dtLocTile = dtLoc;
    const Real dtTile = dt;
    const bool doSaveTile = doSave;
    const bool needGradB = outputSpec.needGradB;
    SlotStage* stage = doSaveTile ? &stages[omp_slot_id()] : nullptr;
    MoveTallyDelta localTally;

    for (auto& p : particles) {
      if (p.id() < 0)
        continue;

      const Real vIn[3] = { p.rdata(iup_), p.rdata(ivp_), p.rdata(iwp_) };
      const Real xp = p.pos(ix_);
      const Real yp = p.pos(iy_);
      const Real zp = nDim > 2 ? p.pos(iz_) : 0.0;

      IntVect loIdx;
      RealVect dShift;
      Real coef[2][2][2];
      find_node_interpolation(p.pos(), ploTile, invDxTile, loIdx, dShift, coef);

      Real bp[3] = { 0.0, 0.0, 0.0 };
      Real ep[3] = { 0.0, 0.0, 0.0 };
      const Array4<Real const> fields[2] = { ctx.BArr, ctx.EArr };
      Real* values[2] = { bp, ep };
      interpolate_vector_fields(fields, loIdx, coef, lo, hi, values);
      if (hasB0) {
        Real b0p[3] = { 0.0, 0.0, 0.0 };
        interpolate_vector_field(ctx.B0Arr, loIdx, coef, lo, hi, b0p);
        bp[0] += b0p[0];
        bp[1] += b0p[1];
        bp[2] += b0p[2];
      }
      if (hasE0) {
        Real e0p[3] = { 0.0, 0.0, 0.0 };
        interpolate_vector_field(ctx.E0Arr, loIdx, coef, lo, hi, e0p);
        ep[0] += e0p[0];
        ep[1] += e0p[1];
        ep[2] += e0p[2];
      }

      Real vOut[3] = { vIn[0], vIn[1], vIn[2] };
      if (!push_lorentz_velocity(pusherTile, isRelTile, vIn, ep, bp,
                                 qdto2mcTile, vOut)) {
        superluminal.store(true, std::memory_order_relaxed);
      }

      p.rdata(iup_) = vOut[0];
      p.rdata(ivp_) = vOut[1];
      p.rdata(iwp_) = vOut[2];

      p.pos(ix_) = xp + vOut[0] * dtLocTile;
      p.pos(iy_) = yp + vOut[1] * dtLocTile;
      if (nDim > 2)
        p.pos(iz_) = zp + vOut[2] * dtLocTile;

      // Apply passive boundary condition (reflect or delete; never mutates
      // PIC) before recording so reflected coordinates stay inside walls and
      // absorbed particles do not sample unphysical interior/exterior cells.
      if (reflect_or_delete_particle(p, ctx.status, ctx.ab, {}, {},
                                     &localTally)) {
        p.id() = -1;
        continue;
      }

      if (doSaveTile) {
        Real tRec, xRec, yRec, zRec, uRec, vRec, wRec;
        interpolate_record_trajectory(xp, yp, zp, vIn[0], vIn[1], vIn[2],
                                      vOut[0], vOut[1], vOut[2], dtLocTile,
                                      Real(0.5) * dtTile, tNowSI, tRec, xRec,
                                      yRec, zRec, uRec, vRec, wRec);
        reflect_record_at_walls(bc, lb.plo, lb.phi, xRec, yRec, zRec, uRec,
                                vRec, wRec);

        Real gradB[3][3] = { { 0.0 } };
        if (needGradB) {
          if (hasJacB) {
            interpolate_jacobian_matrix(ctx.jacBArr, loIdx, coef, lo, hi,
                                        gradB);
          } else {
            compute_nodal_grad_b(ctx.BArr, loIdx, dShift, invDxTile, gradB);
          }
        }

        stage->metas.push_back({ static_cast<int>(p.cpu()),
                                 static_cast<int>(p.id()), p.rdata(iqp_) });
        append_particle_record(stage->values, outputSpec, mass, charge,
                               isRelTile, tRec, xRec, yRec, zRec, uRec, vRec,
                               wRec, bp, ep, needGradB ? gradB : nullptr);
      }
    }
  });

  if (superluminal.load(std::memory_order_relaxed))
    Abort("Error: particle speed is fast than light!");

  if (doSave)
    commit_stages(stages);
}

//==========================================================
void TestParticles::move_and_save_neutrals(int iLev, Real dt, Real tNowSI,
                                           bool doSave) {
  timing_func("TestParticles::move_neutrals");

  const LevelGeomBox lb = make_level_geom_box(iLev);
  const BoxArray& ba = ParticleBoxArray(iLev);

  struct NeutralFabCtx {
    ActiveRegionBox ab;
    Array4<int const> status;
  };
  const auto fabTable =
      build_local_fab_table(cell_status(iLev), [&](const MFIter& mfi) {
        const Box& bx = cell_status(iLev)[mfi].box();
        return NeutralFabCtx{ make_active_region_box(lb, ba[mfi.index()], bx),
                              cell_status(iLev)[mfi].const_array() };
      });

  std::vector<SlotStage> stages(doSave ? omp_slot_count() : 0);

  parallel_for_cell_tiles(*this, iLev, [&](const auto& item, int /*ord*/) {
    const auto& ctx = fabTable[item.grid];
    auto& particles = item.ptile->GetArrayOfStructs();
    const bool doSaveTile = doSave;
    const bool isRelTile = isRelativistic;
    SlotStage* stage = doSaveTile ? &stages[omp_slot_id()] : nullptr;
    MoveTallyDelta localTally;

    for (auto& p : particles) {
      if (p.id() < 0)
        continue;

      const Real up = p.rdata(iup_);
      const Real vp = p.rdata(ivp_);
      const Real wp = p.rdata(iwp_);
      const Real xp = p.pos(ix_);
      const Real yp = p.pos(iy_);
      const Real zp = nDim > 2 ? p.pos(iz_) : 0.0;

      p.pos(ix_) = xp + up * dt;
      p.pos(iy_) = yp + vp * dt;
      if (nDim > 2)
        p.pos(iz_) = zp + wp * dt;

      if (reflect_or_delete_particle(p, ctx.status, ctx.ab, {}, {},
                                     &localTally)) {
        p.id() = -1;
        continue;
      }

      if (doSaveTile) {
        Real tRec, xRec, yRec, zRec, uRec, vRec, wRec;
        interpolate_record_trajectory(xp, yp, zp, up, vp, wp, up, vp, wp, dt,
                                      dt, tNowSI, tRec, xRec, yRec, zRec, uRec,
                                      vRec, wRec);
        reflect_record_at_walls(bc, lb.plo, lb.phi, xRec, yRec, zRec, uRec,
                                vRec, wRec);
        stage->metas.push_back({ static_cast<int>(p.cpu()),
                                 static_cast<int>(p.id()), p.rdata(iqp_) });
        append_particle_record(stage->values, outputSpec, mass, charge,
                               isRelTile, tRec, xRec, yRec, zRec, uRec, vRec,
                               wRec, nullptr, nullptr, nullptr);
      }
    }
  });

  if (doSave)
    commit_stages(stages);
}

//======================================================================
void TestParticles::read_test_particle_list(
    const Vector<std::string>& listFiles) {
  std::string funcName = "TP::read_test_particle_list";
  timing_func(funcName);

  std::string listName;
  for (int i = 0; i < listFiles.size(); ++i) {
    // Assume the particle list file's name is *N.dat, where N is
    // the species ID (Usually 0 or 1).
    if (listFiles[i].find(std::to_string(speciesID) + ".dat") !=
        std::string::npos) {
      listName = listFiles[i];
      break;
    }
  }

  if (listName.empty())
    return;

  vIDs.clear();

  // Read particle IDs
  std::ifstream source;
  source.open(listName, std::ifstream::in);
  PID id;
  while (source >> id.cpu >> id.id) {
    id.flag = true;
    vIDs.push_back(id);
  }

  std::sort(vIDs.begin(), vIDs.end());
}

//======================================================================
void TestParticles::add_test_particles_from_pic(PicParticles* pts) {
  std::string funcName = "TP::add_test_particles_from_pic";
  timing_func(funcName);

  if (!pts)
    Abort("Error [TestParticles]: #TPINITFROMPIC requires #PIC T.");

  if (vIDs.empty())
    return;

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    const auto& lOther = pts->GetParticles(iLev);

    for (MFIter mfi = pts->MakeMFIter(iLev); mfi.isValid(); ++mfi) {
      auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());

      if (lOther.find(index) == lOther.end())
        continue;

      const auto& tileOther = lOther.at(index);

      if (tileOther.numParticles() == 0)
        continue;

      auto& particles = GetParticles(iLev)[index];

      const auto& aosOther = tileOther.GetArrayOfStructs();

      for (const auto& pOther : aosOther) {
        if (pOther.id() < 0)
          continue;

        PID id;
        id.cpu = pOther.cpu();
        id.id = pOther.id();

        auto it = std::lower_bound(vIDs.begin(), vIDs.end(), id);

        if (it != vIDs.end() && *it == id) {
          auto p = make_particle();
          p.cpu() = pOther.cpu();
          p.id() = pOther.id();
          p.idata(iSupID_) = pOther.idata(iSupID_);
          p.pos(ix_) = pOther.pos(ix_);
          p.pos(iy_) = pOther.pos(iy_);
          if (nDim > 2)
            p.pos(iz_) = pOther.pos(iz_);
          p.rdata(iup_) = pOther.rdata(iup_);
          p.rdata(ivp_) = pOther.rdata(ivp_);
          p.rdata(iwp_) = pOther.rdata(iwp_);
          p.rdata(iqp_) = pOther.rdata(iqp_);

          particles.push_back(p);

          (*it).flag = false;
        }
      }
    }
  }

  { // If a test particle has been generated on one proc, all procs should
    // remove this particle from the test particle list (vIDs). Use a bit array
    // to do mpi_allreduce to improve performance.
    BitArray ba(vIDs.size());
    for (unsigned int i = 0; i < vIDs.size(); ++i) {
      if (vIDs[i].flag) {
        ba.set(i, 1);
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, ba.get(), ba.size_int(), MPI_INT, MPI_BAND,
                  ParallelDescriptor::Communicator());

    for (unsigned int i = 0; i < vIDs.size(); ++i) {
      if (ba.get(i) == 0) {
        vIDs[i].flag = false;
      }
    }
  }

  vIDs.erase(std::remove_if(vIDs.begin(), vIDs.end(),
                            [](const PID& t) { return !t.flag; }),
             vIDs.end());

  redistribute_particles();
}

//======================================================================
void TestParticles::add_test_particles_from_fluid(const Vector<Vel>& tpStates) {
  std::string funcName = "TP::add_test_particles_from_fluid";
  timing_func(funcName);
  Print() << funcName << " : nInitPart = " << nInitPart
          << " : current number = " << TotalNumberOfParticles(true, false)
          << std::endl;

  Vel tpVel;
  for (const auto& state : tpStates) {
    if (state.tag == speciesID) {
      tpVel = state;
    }
  }

  for (int iLev = 0; iLev < n_lev(); ++iLev) {
    const int domHiX = Geom(iLev).Domain().bigEnd(ix_);
    for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
      const auto& status = cell_status(iLev)[mfi].const_array();
      const Box& bx = mfi.validbox();
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int i = lo.x; i <= hi.x; i += nIntervalCell[ix_])
        for (int j = lo.y; j <= hi.y; j += nIntervalCell[iy_])
#if (AMREX_SPACEDIM == 3)
          for (int k = lo.z; k <= hi.z; k += nIntervalCell[iz_])
#else
          for (int k = 0; k <= 0; ++k)
#endif
          {
            const IntVect ijk{ AMREX_D_DECL(i, j, k) };
            const int st = status(ijk);
            // Seed only leaf cells across the AMR hierarchy.
            if (bit::is_refined(st))
              continue;

            const bool isEdge =
                (iLev == 0) ? bit::is_lev_edge(st) : bit::is_domain_edge(st);

            if (iPartRegion == iRegionUniform_ ||
                (iPartRegion == iRegionBoundary_ && isEdge)) {
              add_particles_cell(iLev, mfi, ijk, fi, false, IntVect(0), tpVel);
            } else if (iPartRegion == iRegionUser_) {
              Real xyz[nDim];
              Geom(iLev).CellCenter(ijk, xyz);
              if (tpRegions.is_inside(xyz)) {
                add_particles_cell(iLev, mfi, ijk, fi, false, IntVect(0),
                                   tpVel);
              }
            }
          }

      if (iPartRegion == iRegionSideXp_ && (iLev == 0 || hi.x == domHiX)) {
        const int i = hi.x;
        for (int j = lo.y; j <= hi.y; j += nIntervalCell[iy_])
#if (AMREX_SPACEDIM == 3)
          for (int k = lo.z; k <= hi.z; k += nIntervalCell[iz_]) {
#else
          for (int k = 0; k <= 0; ++k) {
#endif
            const IntVect ijk{ AMREX_D_DECL(i, j, k) };
            const int st = status(ijk);
            if (bit::is_refined(st))
              continue;
            const bool isSideXp =
                (iLev == 0) ? bit::is_lev_edge(st) : bit::is_domain_edge(st);
            if (isSideXp)
              add_particles_cell(iLev, mfi, ijk, fi, false, IntVect(0), tpVel);
          }
      }
    }
  }

  redistribute_particles();
}

//======================================================================
void TestParticles::coalesce_trajectories(
    std::vector<TrajEntry>& sortedEntries) {
  const int nProcs = ParallelDescriptor::NProcs();
  const int myProc = ParallelDescriptor::MyProc();
  const std::size_t rSize = static_cast<std::size_t>(outputSpec.recordSize);

  if (nProcs > 1) {
    std::vector<std::vector<char> > sendBufs(nProcs);
    std::unordered_map<std::uint64_t, TrajEntry> ownedMap;

    auto appendBytes = [](std::vector<char>& buf, const void* ptr,
                          std::size_t n) {
      const char* c = static_cast<const char*>(ptr);
      buf.insert(buf.end(), c, c + n);
    };

    for (auto& kv : trajMap_) {
      TrajEntry& e = kv.second;
      if (e.data.empty())
        continue;
      const int owner = ((e.cpu % nProcs) + nProcs) % nProcs;
      if (owner == myProc) {
        ownedMap.emplace(kv.first, std::move(e));
      } else {
        const int nRec = static_cast<int>(e.data.size() / rSize);
        auto& buf = sendBufs[owner];
        appendBytes(buf, &e.cpu, sizeof(int));
        appendBytes(buf, &e.id, sizeof(int));
        appendBytes(buf, &nRec, sizeof(int));
        appendBytes(buf, &e.qp, sizeof(Real));
        appendBytes(buf, e.data.data(), e.data.size() * sizeof(Real));
      }
    }

    Vector<int> sendCounts(nProcs, 0), recvCounts(nProcs, 0);
    Vector<int> sdispls(nProcs, 0), rdispls(nProcs, 0);
    int totalSend = 0;
    for (int r = 0; r < nProcs; ++r) {
      sendCounts[r] = static_cast<int>(sendBufs[r].size());
      sdispls[r] = totalSend;
      totalSend += sendCounts[r];
    }

    MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT,
                 ParallelDescriptor::Communicator());

    int totalRecv = 0;
    for (int r = 0; r < nProcs; ++r) {
      rdispls[r] = totalRecv;
      totalRecv += recvCounts[r];
    }

    std::vector<char> sendFlat(totalSend);
    for (int r = 0; r < nProcs; ++r) {
      if (sendCounts[r] > 0) {
        std::memcpy(sendFlat.data() + sdispls[r], sendBufs[r].data(),
                    sendCounts[r]);
      }
    }

    std::vector<char> recvFlat(totalRecv);
    MPI_Alltoallv(sendFlat.data(), sendCounts.data(), sdispls.data(), MPI_BYTE,
                  recvFlat.data(), recvCounts.data(), rdispls.data(), MPI_BYTE,
                  ParallelDescriptor::Communicator());

    std::size_t off = 0;
    while (off < recvFlat.size()) {
      int cpu = 0, id = 0, nRec = 0;
      Real qp = 0.0;
      std::memcpy(&cpu, recvFlat.data() + off, sizeof(int));
      off += sizeof(int);
      std::memcpy(&id, recvFlat.data() + off, sizeof(int));
      off += sizeof(int);
      std::memcpy(&nRec, recvFlat.data() + off, sizeof(int));
      off += sizeof(int);
      std::memcpy(&qp, recvFlat.data() + off, sizeof(Real));
      off += sizeof(Real);

      const std::size_t nVals = static_cast<std::size_t>(nRec) * rSize;
      auto& dst = ownedMap[pack_tp_key(cpu, id)];
      if (dst.data.empty()) {
        dst.cpu = cpu;
        dst.id = id;
        dst.qp = qp;
      } else {
        dst.merged = true;
      }
      const std::size_t oldSize = dst.data.size();
      dst.data.resize(oldSize + nVals);
      std::memcpy(dst.data.data() + oldSize, recvFlat.data() + off,
                  nVals * sizeof(Real));
      off += nVals * sizeof(Real);
    }

    trajMap_.swap(ownedMap);
  }

  sortedEntries.clear();
  sortedEntries.reserve(trajMap_.size());

  const bool backwardTime = (tc != nullptr && tc->get_dt_si() < 0.0);
  for (auto& kv : trajMap_) {
    TrajEntry& e = kv.second;
    if (e.data.empty())
      continue;

    if (e.merged) {
      const std::size_t nRec = e.data.size() / rSize;
      std::vector<std::size_t> order(nRec);
      for (std::size_t i = 0; i < nRec; ++i)
        order[i] = i;
      std::stable_sort(order.begin(), order.end(),
                       [&](std::size_t a, std::size_t b) {
                         const Real ta = e.data[a * rSize];
                         const Real tb = e.data[b * rSize];
                         return backwardTime ? (ta > tb) : (ta < tb);
                       });
      std::vector<Real> sortedData(e.data.size());
      for (std::size_t i = 0; i < nRec; ++i) {
        std::memcpy(&sortedData[i * rSize], &e.data[order[i] * rSize],
                    rSize * sizeof(Real));
      }
      e.data.swap(sortedData);
    }

    sortedEntries.push_back(std::move(e));
  }

  trajMap_.clear();

  std::sort(sortedEntries.begin(), sortedEntries.end(),
            [](const TrajEntry& a, const TrajEntry& b) {
              return (a.cpu != b.cpu) ? (a.cpu < b.cpu) : (a.id < b.id);
            });
}

//======================================================================
void TestParticles::pack_particles_and_records(
    const std::vector<TrajEntry>& entries, char* recordBuff, char* listBuff,
    unsigned long long int shift) const {
  constexpr std::size_t listUnitSize =
      2 * sizeof(int) + sizeof(unsigned long long);
  const int rSize = outputSpec.recordSize;
  unsigned long long int nByteCount = 0;
  std::vector<float> recordData(rSize);

  for (std::size_t idx = 0; idx < entries.size(); ++idx) {
    const auto& e = entries[idx];
    const int nRecord = static_cast<int>(e.data.size() / rSize);
    const unsigned long long int nBytePerPart =
        3ULL * sizeof(int) +
        (1ULL + static_cast<unsigned long long int>(rSize) * nRecord) *
            sizeof(float);

    if (listBuff) {
      const int i2[2] = { e.cpu, e.id };
      std::memcpy(listBuff + idx * listUnitSize, i2, 2 * sizeof(int));
      const unsigned long long tmp = nByteCount + shift;
      std::memcpy(listBuff + idx * listUnitSize + 2 * sizeof(int), &tmp,
                  sizeof(unsigned long long));
    }

    if (recordBuff) {
      std::size_t iCountLoc = 0;
      const int i3[3] = { e.cpu, e.id, nRecord };
      std::memcpy(recordBuff + nByteCount + iCountLoc, i3, 3 * sizeof(int));
      iCountLoc += 3 * sizeof(int);

      const float weight = static_cast<float>(e.qp / charge * mass * no2outM);
      std::memcpy(recordBuff + nByteCount + iCountLoc, &weight, sizeof(float));
      iCountLoc += sizeof(float);

      const std::size_t recBytes = sizeof(float) * rSize;
      for (int i = 0; i < nRecord; ++i) {
        const Real* src = &e.data[static_cast<std::size_t>(i) * rSize];
        for (int c = 0; c < rSize; ++c) {
          recordData[c] = static_cast<float>(src[c] * colScales[c]);
        }
        std::memcpy(recordBuff + nByteCount + iCountLoc, recordData.data(),
                    recBytes);
        iCountLoc += recBytes;
      }
    }

    nByteCount += nBytePerPart;
  }
}

//======================================================================
bool TestParticles::write_particles(int cycle) {
  std::string funcName = "TP::write_particles";
  timing_func(funcName);

  std::vector<TrajEntry> entries;
  coalesce_trajectories(entries);

  const int rSize = outputSpec.recordSize;
  int nPartLoc = static_cast<int>(entries.size());
  int nPartAhead = 0;
  gather_accumulate_and_scatter(nPartLoc, nPartAhead);

  unsigned long long int nByteLoc = 0;
  for (const auto& e : entries) {
    const unsigned long long int nRecord = e.data.size() / rSize;
    nByteLoc += 3ULL * sizeof(int) +
                (1ULL + static_cast<unsigned long long int>(rSize) * nRecord) *
                    sizeof(float);
  }
  unsigned long long int nByteAhead = 0;
  gather_accumulate_and_scatter(nByteLoc, nByteAhead);

  constexpr int listUnitSize = 2 * sizeof(int) + sizeof(unsigned long long);
  Vector<char> dataBuffer(nByteLoc);
  Vector<char> partList(static_cast<std::size_t>(nPartLoc) * listUnitSize);

  pack_particles_and_records(entries, dataBuffer.data(), partList.data(),
                             nByteAhead);

  std::stringstream ss;
  ss << "n" << std::setfill('0') << std::setw(8) << cycle;

  UtilCreateDirectory(outputDir, 0755);
  std::string fileNamePartList = outputDir + "/" + gridName +
                                 "_particle_list_species_" +
                                 std::to_string(speciesID) + "_" + ss.str();

  std::string fileNamePartRecord = outputDir + "/" + gridName +
                                   "_particle_species_" +
                                   std::to_string(speciesID) + "_" + ss.str();

  Print() << printPrefix << "Saving test particles into " << fileNamePartRecord
          << std::endl;

  MPI_File recordFile, listFile;
  MPI_Status status;

  MPI_File_open(ParallelDescriptor::Communicator(), fileNamePartList.c_str(),
                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &listFile);
  MPI_File_set_size(listFile, 0);
  MPI_File_write_at_all(listFile,
                        static_cast<MPI_Offset>(nPartAhead) * listUnitSize,
                        partList.data(), partList.size(), MPI_CHAR, &status);
  MPI_File_close(&listFile);

  MPI_File_open(ParallelDescriptor::Communicator(), fileNamePartRecord.c_str(),
                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &recordFile);
  MPI_File_set_size(recordFile, 0);
  MPI_File_write_at_all(recordFile, static_cast<MPI_Offset>(nByteAhead),
                        dataBuffer.data(), dataBuffer.size(), MPI_CHAR,
                        &status);
  MPI_File_close(&recordFile);

  if (ParallelDescriptor::IOProcessor()) {
    std::string headerName = outputDir + "/Header";

    std::ofstream headerFile;
    headerFile.open(headerName.c_str(), std::ios::out | std::ios::trunc);

    if (!headerFile.good())
      FileOpenFailed(headerName);

    headerFile << outputSpec.recordSize << "\n";
    for (std::size_t i = 0; i < outputSpec.colNames.size(); ++i) {
      if (i > 0)
        headerFile << ' ';
      headerFile << outputSpec.colNames[i];
    }
    headerFile << "\n";
  }

  return true;
}

//======================================================================
void TestParticles::reset_record_counter() { trajMap_.clear(); }
