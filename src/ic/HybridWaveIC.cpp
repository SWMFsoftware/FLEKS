#include "HybridWaveIC.h"

#include <AMReX.H>
#include <AMReX_MultiFab.H>

// Field seeding, routed through the PicICFields facade:
//   * Bx0 is read from the node-B guide field already deposited by
//     Pic::fill_E_B_fields() from #UNIFORMSTATE.
//   * node-B and center-B are ZEROED and re-written with the wave (this is a
//     replace, matching the kernel which calls setVal(0.0) first).
//   * the transverse perturbation uses  kx = 2*pi/Lx  (waveMode is accepted
//     for compatibility but is NOT used in the wavenumber, exactly like the
//     original kernel).
void HybridWaveIC::set_fields(PicICFields& fields) const {
  const int nLev = fields.n_lev();
  B1_.resize(nLev);
  kx_.resize(nLev);

  for (int iLev = 0; iLev < nLev; ++iLev) {
    // Guide field Bx0 already deposited by fill_E_B_fields() from #UNIFORMSTATE.
    amrex::Real Bx0 = 1.0;
    {
      amrex::MFIter mfi(fields.node_B(iLev));
      if (mfi.isValid()) {
        const amrex::Array4<amrex::Real const>& a =
            fields.node_B(iLev).array(mfi);
        const amrex::Box& b = mfi.validbox();
        Bx0 = a(b.smallEnd(), ix_);
      }
    }

    const amrex::Real B1 = waveFrac_ * Bx0;
    const auto& prob_lo = fields.geom(iLev).ProbLo();
    const auto& dx = fields.geom(iLev).CellSize();
    const amrex::Real Lx = (fields.geom(iLev).ProbHi())[0] - prob_lo[0];
    const amrex::Real kx = (Lx > 0) ? 2.0 * (dPI) / Lx : 0.0;

    B1_[iLev] = B1;
    kx_[iLev] = kx;

    fields.node_B(iLev).setVal(0.0);
    fields.center_B(iLev).setVal(0.0);

    for (amrex::MFIter mfi(fields.node_B(iLev)); mfi.isValid(); ++mfi) {
      amrex::FArrayBox& fab = fields.node_B(iLev)[mfi];
      const amrex::Box& box = mfi.fabbox();
      const amrex::Array4<amrex::Real>& arrB = fab.array();
      amrex::ParallelFor(box, [&](int i, int j, int k) {
        amrex::IntVect ijk = {AMREX_D_DECL(i, j, k)};
        const amrex::Real x = prob_lo[0] + dx[0] * i;
        arrB(ijk, ix_) = Bx0;
        arrB(ijk, iy_) = B1 * std::cos(kx * x);
        arrB(ijk, iz_) = B1 * std::sin(kx * x);
      });
    }

    for (amrex::MFIter mfi(fields.center_B(iLev)); mfi.isValid(); ++mfi) {
      amrex::FArrayBox& fab = fields.center_B(iLev)[mfi];
      const amrex::Box& box = mfi.fabbox();
      const amrex::Array4<amrex::Real>& arrB = fab.array();
      amrex::ParallelFor(box, [&](int i, int j, int k) {
        amrex::IntVect ijk = {AMREX_D_DECL(i, j, k)};
        const amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
        arrB(ijk, ix_) = Bx0;
        arrB(ijk, iy_) = B1 * std::cos(kx * x);
        arrB(ijk, iz_) = B1 * std::sin(kx * x);
      });
    }
  }

  fields.fill_boundary_E_B();
}

// Reproduces Pic::perturb_hybrid_wave_velocities(): for a parallel-propagating
// wave with delta B = (0, B1 cos kx, B1 sin kx) on guide field B0, the Alfven
// limit gives delta U_perp = -v_A * delta B_perp / B0. In code units
// (v_A = B0 = 1):  delta U_y = -B1 cos kx,  delta U_z = -B1 sin kx.
// Applies to every particle routed through here (kinetic ions only in a hybrid
// run).
void HybridWaveIC::modify_particle_velocity(ParticleICState& s) const {
  if (s.iLev < 0 || s.iLev >= static_cast<int>(B1_.size()))
    return;
  const amrex::Real B1 = B1_[s.iLev];
  const amrex::Real kx = kx_[s.iLev];
  s.vBulk += -B1 * std::cos(kx * s.x);
  s.wBulk += -B1 * std::sin(kx * s.x);
}
