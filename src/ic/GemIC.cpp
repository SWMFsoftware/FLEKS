#include <algorithm>
#include <cmath>

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include "Constants.h"
#include "GemIC.h"
#include "Particles.h"
#include "ReadParam.h"

using namespace amrex;

namespace {
inline amrex::Real sech2(amrex::Real u) {
  if (std::abs(u) > 40.0)
    return 0.0;
  const amrex::Real ch = std::cosh(u);
  return 1.0 / (ch * ch);
}
} // namespace

void GemIC::read_param(ReadParam& param) {
  bool progress = true;
  while (progress) {
    progress = false;
    if (param.read_optional("b0", b0_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("b1", b1_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("b2", b2_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("tp", tp_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("t0", tp_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("t1", t1_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("t2", t2_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("lambda0", lambda0_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("currentSheetWidth", lambda0_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("currentSheetThickness", lambda0_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("apert", apert_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("amplitude", apert_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("bg", bg_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("nb", nb_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("teOverTi", teOverTi_)) {
      progress = true;
      continue;
    }
    std::string sVal;
    if (param.read_optional("useDoubleCurrentSheet", sVal)) {
      useDoubleCurrentSheet_ = (sVal == "T" || sVal == "true" || sVal == "1");
      progress = true;
      continue;
    }
    if (param.read_optional("isAsymmetryReconnection", sVal)) {
      isAsymmetryReconnection_ = (sVal == "T" || sVal == "true" || sVal == "1");
      progress = true;
      continue;
    }
    if (param.read_optional("useGEMReflected", sVal)) {
      useGEMReflected_ = (sVal == "T" || sVal == "true" || sVal == "1");
      progress = true;
      continue;
    }
    if (param.read_optional("useStandardGem", sVal)) {
      useStandardGem_ = (sVal == "T" || sVal == "true" || sVal == "1");
      progress = true;
      continue;
    }
    if (param.read_optional("useUniformPressure", sVal)) {
      useUniformPressure_ = (sVal == "T" || sVal == "true" || sVal == "1");
      progress = true;
      continue;
    }
    if (param.read_optional("useUniformIonPressure", sVal)) {
      useUniformIonPressure_ = (sVal == "T" || sVal == "true" || sVal == "1");
      progress = true;
      continue;
    }
    if (param.read_optional("gaussX", gaussX_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("gaussY", gaussY_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("waveLengthX", waveLengthX_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("waveLengthY", waveLengthY_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("pertType", pertType_)) {
      progress = true;
      continue;
    }
  }
}

void GemIC::init_geometry(amrex::Real Lx, amrex::Real Ly) const {
  Lx_ = Lx;
  Ly_ = Ly;
  Wx_ = (waveLengthX_ > 0.0) ? waveLengthX_ : Lx;
  Wy_ = (waveLengthY_ > 0.0) ? waveLengthY_ : Ly;

  Kx_ = (Wx_ > 0.0) ? (2.0 * dPI / Wx_) : 0.0;
  Ky_ = (Wy_ > 0.0) ? (2.0 * dPI / Wy_) : 0.0;
  gaussXInv_ = (gaussX_ > 0.0) ? (1.0 / gaussX_) : 0.0;
  gaussYInv_ = (gaussY_ > 0.0) ? (1.0 / gaussY_) : 0.0;

  if (isAsymmetryReconnection_) {
    b0Eff_ = 0.5 * (b1_ + b2_);
    tpEff_ = 0.5 * (t1_ + t2_);
    pB_ = (b1_ > b2_) ? (b1_ * b1_) : (b2_ * b2_);
    ySheet_ = 0.25 * Wy_;
    xT_ = 0.25 * Wx_;
    xB_ = -0.25 * Wx_;
  } else if (useDoubleCurrentSheet_) {
    b0Eff_ = b0_;
    tpEff_ = tp_;
    pB_ = b0_ * b0_;
    ySheet_ = 0.25 * Wy_;
    xT_ = 0.25 * Wx_;
    xB_ = -0.25 * Wx_;
  } else if (useGEMReflected_) {
    b0Eff_ = b0_;
    tpEff_ = tp_;
    pB_ = b0_ * b0_;
    ySheet_ = 0.5 * Wy_;
    xT_ = 0.0;
    xB_ = 0.0;
  } else {
    b0Eff_ = b0_;
    tpEff_ = tp_;
    pB_ = b0_ * b0_;
    ySheet_ = 0.0;
    xT_ = 0.0;
    xB_ = 0.0;
  }
  geomInit_ = true;
}

amrex::Real GemIC::eval_Bx0(amrex::Real y) const {
  const amrex::Real l0 = (lambda0_ > 0.0) ? lambda0_ : 1.0;
  if (useDoubleCurrentSheet_) {
    return b0Eff_ * (std::tanh((y + ySheet_) / l0) -
                     std::tanh((y - ySheet_) / l0) - 1.0);
  } else if (isAsymmetryReconnection_) {
    return -b0Eff_ * (std::tanh((y - 0.25 * Wy_) / l0) -
                      std::tanh((y - 0.75 * Wy_) / l0) +
                      std::tanh((y - 1.25 * Wy_) / l0) -
                      std::tanh((y + 0.25 * Wy_) / l0) + 1.0) +
           0.5 * (b1_ - b2_);
  } else if (useGEMReflected_) {
    if (y > 0.0)
      return b0Eff_ * std::tanh((y - ySheet_) / l0);
    else
      return -b0Eff_ * std::tanh((y + ySheet_) / l0);
  }
  return b0Eff_ * std::tanh(y / l0);
}

amrex::Real GemIC::eval_Jz(amrex::Real y) const {
  const amrex::Real l0 = (lambda0_ > 0.0) ? lambda0_ : 1.0;
  const amrex::Real factor = -b0Eff_ / l0;
  if (useDoubleCurrentSheet_) {
    return factor * (sech2((y + ySheet_) / l0) - sech2((y - ySheet_) / l0));
  } else if (isAsymmetryReconnection_) {
    return factor *
           (sech2((y - 0.25 * Wy_) / l0) - sech2((y - 0.75 * Wy_) / l0) +
            sech2((y - 1.25 * Wy_) / l0) - sech2((y + 0.25 * Wy_) / l0));
  } else if (useGEMReflected_) {
    if (y > 0.0)
      return factor * sech2((y - ySheet_) / l0);
    else
      return -factor * sech2((y + ySheet_) / l0);
  }
  return factor * sech2(y / l0);
}

void GemIC::eval_perturbation(amrex::Real x, amrex::Real y, amrex::Real& dbx,
                              amrex::Real& dby) const {
  dbx = 0.0;
  dby = 0.0;
  if (apert_ == 0.0)
    return;

  // Periodic reduction of x into [-Wx/2, Wx/2] for multi-island periodic tiling
  const amrex::Real xRel = (Wx_ > 0.0) ? (x - Wx_ * std::round(x / Wx_)) : x;

  if (useDoubleCurrentSheet_ || isAsymmetryReconnection_) {
    const amrex::Real a1 =
        -apert_ * b0Eff_ *
        std::exp(-((xRel - xT_) * (xRel - xT_)) * (gaussXInv_ * gaussXInv_) -
                 ((y - ySheet_) * (y - ySheet_)) * (gaussYInv_ * gaussYInv_));
    const amrex::Real a2 =
        apert_ * b0Eff_ *
        std::exp(-((xRel - xB_) * (xRel - xB_)) * (gaussXInv_ * gaussXInv_) -
                 ((y + ySheet_) * (y + ySheet_)) * (gaussYInv_ * gaussYInv_));

    dbx =
        a1 * (-2.0 * (y - ySheet_) * (gaussYInv_ * gaussYInv_) *
                  std::cos(Kx_ * (xRel - xT_)) * std::cos(Ky_ * (y - ySheet_)) -
              Ky_ * std::cos(Kx_ * (xRel - xT_)) *
                  std::sin(Ky_ * (y - ySheet_))) +
        a2 * (-2.0 * (y + ySheet_) * (gaussYInv_ * gaussYInv_) *
                  std::cos(Kx_ * (xRel - xB_)) * std::cos(Ky_ * (y + ySheet_)) -
              Ky_ * std::cos(Kx_ * (xRel - xB_)) *
                  std::sin(Ky_ * (y + ySheet_)));

    dby =
        a1 * (2.0 * (xRel - xT_) * (gaussXInv_ * gaussXInv_) *
                  std::cos(Kx_ * (xRel - xT_)) * std::cos(Ky_ * (y - ySheet_)) +
              Kx_ * std::sin(Kx_ * (xRel - xT_)) *
                  std::cos(Ky_ * (y - ySheet_))) +
        a2 * (2.0 * (xRel - xB_) * (gaussXInv_ * gaussXInv_) *
                  std::cos(Kx_ * (xRel - xB_)) * std::cos(Ky_ * (y + ySheet_)) +
              Kx_ * std::sin(Kx_ * (xRel - xB_)) *
                  std::cos(Ky_ * (y + ySheet_)));
  } else if (useGEMReflected_) {
    const amrex::Real kx = 2.0 * dPI / Lx_;
    const amrex::Real ky = 2.0 * dPI / Ly_;
    dbx = -2.0 * apert_ * b0Eff_ * (dPI / Ly_) * std::cos(kx * x) *
          std::sin(ky * (y - ySheet_));
    dby = apert_ * b0Eff_ * (2.0 * dPI / Lx_) * std::sin(kx * x) *
          std::cos(ky * (y - ySheet_));
  } else if (useStandardGem_ && pertType_ != "gaussian") {
    // Classic GEM challenge modal perturbation
    const amrex::Real kx = 2.0 * dPI / Lx_;
    const amrex::Real ky = dPI / Ly_;
    dbx = -apert_ * b0Eff_ * ky * std::cos(kx * x) * std::sin(ky * y);
    dby = apert_ * b0Eff_ * kx * std::sin(kx * x) * std::cos(ky * y);
  } else {
    // Single sheet with Gaussian perturbation
    const amrex::Real a = apert_ * b0Eff_ *
                          std::exp(-(xRel * xRel) * (gaussXInv_ * gaussXInv_) -
                                   (y * y) * (gaussYInv_ * gaussYInv_));
    dbx = a * (-2.0 * y * (gaussYInv_ * gaussYInv_) * std::cos(Kx_ * xRel) *
                   std::cos(Ky_ * y) -
               Ky_ * std::cos(Kx_ * xRel) * std::sin(Ky_ * y));
    dby = a * (2.0 * xRel * (gaussXInv_ * gaussXInv_) * std::cos(Kx_ * xRel) *
                   std::cos(Ky_ * y) +
               Kx_ * std::sin(Kx_ * xRel) * std::cos(Ky_ * y));
  }
}

amrex::Real GemIC::eval_Tp(amrex::Real y) const {
  if (isAsymmetryReconnection_) {
    const amrex::Real l0 = (lambda0_ > 0.0) ? lambda0_ : 1.0;
    return t2_ +
           0.5 * (t1_ - t2_) *
               (std::tanh((y + ySheet_) / l0) - std::tanh((y - ySheet_) / l0));
  }
  return tpEff_;
}

amrex::Real GemIC::eval_density(amrex::Real y, amrex::Real Bx0) const {
  amrex::Real deltaP = 0.5 * (pB_ - Bx0 * Bx0);
  if (deltaP < 0.0)
    deltaP = 0.0;
  const amrex::Real Ptot = nb_ * tpEff_ + deltaP;
  const amrex::Real Tp = eval_Tp(y);
  if (Tp <= 0.0)
    return nb_;
  return Ptot / Tp;
}

void GemIC::set_fields(PicICFields& fields) const {
  const int nLev = fields.n_lev();
  if (nLev == 0)
    return;

  for (int iLev = 0; iLev < nLev; ++iLev) {
    const auto& geom = fields.geom(iLev);
    const auto& prob_lo = geom.ProbLo();
    const auto& prob_hi = geom.ProbHi();
    const auto& dx = geom.CellSize();

    const amrex::Real Lx = prob_hi[0] - prob_lo[0];
    const amrex::Real Ly = prob_hi[1] - prob_lo[1];
    init_geometry(Lx, Ly);

    MultiFab& nodeB = fields.node_B(iLev);
    MultiFab& centerB = fields.center_B(iLev);
    MultiFab& nodeE = fields.node_E(iLev);

    nodeB.setVal(0.0);
    centerB.setVal(0.0);
    nodeE.setVal(0.0);

    const amrex::Real bg = bg_;
    const amrex::Real b0Eff = b0Eff_;
    const bool useUniformPressure = useUniformPressure_;

    // Node-centered B
    for (MFIter mfi(nodeB); mfi.isValid(); ++mfi) {
      FArrayBox& fab = nodeB[mfi];
      const Box& box = mfi.fabbox();
      const Array4<Real>& arrB = fab.array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = prob_lo[0] + dx[0] * i;
        const amrex::Real y = prob_lo[1] + dx[1] * j;
        const amrex::Real bx0 = eval_Bx0(y);
        amrex::Real dbx = 0.0, dby = 0.0;
        eval_perturbation(x, y, dbx, dby);

        amrex::Real bz = b0Eff * bg;
        if (useUniformPressure) {
          const amrex::Real bz2 =
              b0Eff * b0Eff - bx0 * bx0 + (b0Eff * bg) * (b0Eff * bg);
          bz = (bz2 > 0.0) ? std::sqrt(bz2) : 0.0;
        }

        arrB(i, j, k, ix_) = bx0 + dbx;
        arrB(i, j, k, iy_) = dby;
        arrB(i, j, k, iz_) = bz;
      });
    }

    // Cell-centered B
    for (MFIter mfi(centerB); mfi.isValid(); ++mfi) {
      FArrayBox& fab = centerB[mfi];
      const Box& box = mfi.fabbox();
      const Array4<Real>& arrB = fab.array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
        const amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
        const amrex::Real bx0 = eval_Bx0(y);
        amrex::Real dbx = 0.0, dby = 0.0;
        eval_perturbation(x, y, dbx, dby);

        amrex::Real bz = b0Eff * bg;
        if (useUniformPressure) {
          const amrex::Real bz2 =
              b0Eff * b0Eff - bx0 * bx0 + (b0Eff * bg) * (b0Eff * bg);
          bz = (bz2 > 0.0) ? std::sqrt(bz2) : 0.0;
        }

        arrB(i, j, k, ix_) = bx0 + dbx;
        arrB(i, j, k, iy_) = dby;
        arrB(i, j, k, iz_) = bz;
      });
    }
  }

  fields.fill_boundary_E_B();
}

void GemIC::modify_particle_weight(ParticleICState& s) const {
  if (!geomInit_ || nb_ <= 0.0)
    return;
  const amrex::Real bx0 = eval_Bx0(s.y);
  const amrex::Real n = eval_density(s.y, bx0);
  const amrex::Real factor = n / nb_;
  s.q *= (factor > 0.0) ? factor : 0.0;
}

void GemIC::modify_particle_velocity(ParticleICState& s) const {
  if (!geomInit_ || s.charge == 0.0)
    return;
  const amrex::Real bx0 = eval_Bx0(s.y);
  const amrex::Real n = eval_density(s.y, bx0);
  if (n <= 0.0)
    return;

  const amrex::Real jz = eval_Jz(s.y);
  amrex::Real uz = 0.0;

  if (useUniformIonPressure_) {
    if (s.charge > 0.0) {
      uz = 0.0; // ions carry no current
    } else {
      uz = -jz / (s.charge * n); // electrons carry full current
    }
  } else {
    // Distribute current according to species temperatures (diamagnetic drift)
    amrex::Real frac = 0.5;
    if (teOverTi_ >= 0.0) {
      if (s.charge > 0.0) {
        frac = 1.0 / (1.0 + teOverTi_);
      } else {
        frac = teOverTi_ / (1.0 + teOverTi_);
      }
    }
    uz = frac * jz / (s.charge * n);
  }

  s.wBulk += uz;
}

void GemIC::modify_particle_thermal_velocity(ParticleICState& s) const {
  if (!geomInit_ || !isAsymmetryReconnection_)
    return;
  const amrex::Real tp = eval_Tp(s.y);
  if (tpEff_ <= 0.0 || tp <= 0.0)
    return;
  const amrex::Real factor = std::sqrt(tp / tpEff_);
  s.uThermal *= factor;
  s.vThermal *= factor;
  s.wThermal *= factor;
}
