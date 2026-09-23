#include <cmath>

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include "Constants.h"
#include "ForceFreeIC.h"
#include "ReadParam.h"

using namespace amrex;

void ForceFreeIC::read_param(ReadParam& param) {
  bool progress = true;
  while (progress) {
    progress = false;
    if (param.read_optional("lambda", lambda_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("lambda0", lambda_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("L", lambda_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("b0", b0_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("bg", bg_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("perturb", perturb_)) {
      progress = true;
      continue;
    }
    if (param.read_optional("dB", perturb_)) {
      progress = true;
      continue;
    }
    std::string sVal;
    if (param.read_optional("useUniformIonPressure", sVal)) {
      useUniformIonPressure_ = (sVal == "T" || sVal == "true" || sVal == "1");
      progress = true;
      continue;
    }
    if (param.read_optional("teOverTi", teOverTi_)) {
      progress = true;
      continue;
    }
  }
}

void ForceFreeIC::set_fields(PicICFields& fields) const {
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
    const amrex::Real invLx = (Lx > 0.0) ? 1.0 / Lx : 0.0;
    const amrex::Real invLy = (Ly > 0.0) ? 1.0 / Ly : 0.0;
    const amrex::Real xMid = 0.5 * (prob_lo[0] + prob_hi[0]);
    const amrex::Real yMid = 0.5 * (prob_lo[1] + prob_hi[1]);

    const amrex::Real b0 = b0_;
    const amrex::Real bg = bg_;
    const amrex::Real lambda = (lambda_ > 0.0) ? lambda_ : 1.0;
    const amrex::Real perturb = perturb_;

    MultiFab& nodeB = fields.node_B(iLev);
    MultiFab& centerB = fields.center_B(iLev);
    MultiFab& nodeE = fields.node_E(iLev);

    nodeB.setVal(0.0);
    centerB.setVal(0.0);
    nodeE.setVal(0.0);

    // Node-centered B
    for (MFIter mfi(nodeB); mfi.isValid(); ++mfi) {
      FArrayBox& fab = nodeB[mfi];
      const Box& box = mfi.fabbox();
      const Array4<Real>& arrB = fab.array();
      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = prob_lo[0] + dx[0] * i - xMid;
        const amrex::Real y = prob_lo[1] + dx[1] * j - yMid;
        const amrex::Real ch = std::cosh(y / lambda);
        const amrex::Real th = std::tanh(y / lambda);
        const amrex::Real sech = 1.0 / ch;
        const amrex::Real dbx = -perturb * b0 * (0.5 * Lx * invLy) *
                                std::cos(2.0 * dPI * x * invLx) *
                                std::sin(dPI * y * invLy);
        const amrex::Real dby = perturb * b0 *
                                std::sin(2.0 * dPI * x * invLx) *
                                std::cos(dPI * y * invLy);
        const amrex::Real bz = std::sqrt(bg * bg + b0 * b0 * sech * sech);

        arrB(i, j, k, ix_) = b0 * th + dbx;
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
        const amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5) - xMid;
        const amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5) - yMid;
        const amrex::Real ch = std::cosh(y / lambda);
        const amrex::Real th = std::tanh(y / lambda);
        const amrex::Real sech = 1.0 / ch;
        const amrex::Real dbx = -perturb * b0 * (0.5 * Lx * invLy) *
                                std::cos(2.0 * dPI * x * invLx) *
                                std::sin(dPI * y * invLy);
        const amrex::Real dby = perturb * b0 *
                                std::sin(2.0 * dPI * x * invLx) *
                                std::cos(dPI * y * invLy);
        const amrex::Real bz = std::sqrt(bg * bg + b0 * b0 * sech * sech);

        arrB(i, j, k, ix_) = b0 * th + dbx;
        arrB(i, j, k, iy_) = dby;
        arrB(i, j, k, iz_) = bz;
      });
    }
  }

  fields.fill_boundary_E_B();
}

void ForceFreeIC::modify_particle_velocity(ParticleICState& s) const {
  if (s.charge == 0.0)
    return;

  const amrex::Real lambda = (lambda_ > 0.0) ? lambda_ : 1.0;
  const amrex::Real b0 = b0_;
  const amrex::Real bg = bg_;

  const amrex::Real ch = std::cosh(s.y / lambda);
  const amrex::Real th = std::tanh(s.y / lambda);
  const amrex::Real sech = 1.0 / ch;
  const amrex::Real sech2 = sech * sech;
  const amrex::Real bz = std::sqrt(bg * bg + b0 * b0 * sech2);
  if (bz <= 0.0)
    return;

  // Equilibrium current density components J = curl(B):
  // J_x = d(Bz)/dy = -b0^2 / (lambda * Bz) * sech^2(y/lambda) * tanh(y/lambda)
  // J_z = -d(Bx)/dy = -b0 / lambda * sech^2(y/lambda)
  const amrex::Real jx = -(b0 * b0 / (lambda * bz)) * sech2 * th;
  const amrex::Real jz = -(b0 / lambda) * sech2;

  // In code units, background plasma number density is n0 = 1.0
  const amrex::Real n0 = 1.0;

  amrex::Real ux = 0.0;
  amrex::Real uz = 0.0;

  if (useUniformIonPressure_) {
    if (s.charge > 0.0) {
      ux = 0.0;
      uz = 0.0;
    } else {
      ux = jx / (s.charge * n0);
      uz = jz / (s.charge * n0);
    }
  } else {
    amrex::Real frac = 0.5;
    if (teOverTi_ >= 0.0) {
      if (s.charge > 0.0) {
        frac = 1.0 / (1.0 + teOverTi_);
      } else {
        frac = teOverTi_ / (1.0 + teOverTi_);
      }
    }
    ux = frac * jx / (s.charge * n0);
    uz = frac * jz / (s.charge * n0);
  }

  s.uBulk += ux;
  s.wBulk += uz;
}
