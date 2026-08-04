#include "WaveIC.h"

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include <cmath>

#include "Particles.h"
#include "ReadParam.h"

// The pure-EM lightwave test seeds no macroparticles (the electron is a
// background fluid exercised only by the field solver). Zero the per-cell count
// so add_particles_domain() returns early and the PIC cost collapses to the
// field solve. Other wave profiles keep their kinetic ions.
void WaveIC::apply_particle_override(ParticlesInfo& pInfo) const {
  if (profile_ == LightWave)
    pInfo.nPartPerCell = amrex::IntVect::Zero;
}

using namespace amrex;

std::string WaveIC::name() const {
  switch (profile_) {
    case LightWave:        return "lightwave";
    case HybridWave:       return "hybridwave";
    case ConvectionWave:   return "convectionwave";
    case IonAcousticWave:  return "ionacousticwave";
    case Generic:          return "waveic";
  }
  return "waveic";
}

// Fill the preset defaults for the active profile. Called at the start of
// read_param so that any #WAVEIC sub-parameter the user supplies overrides the
// preset (whether the sub-params sit under #TESTCASE, as the legacy tests have
// them, or under a dedicated #WAVEIC block).
void WaveIC::apply_preset() {
  switch (profile_) {
    case LightWave:
      seedE_ = true;
      seedB_ = true;
      oblique_ = true;
      guideField_ = false;
      velKick_ = false;
      seedWeight_ = false;
      dir_[0] = 0.6; dir_[1] = 0.8; dir_[2] = 0.0;
      waveLength_ = 48.0;
      frac_ = 0.0;    // unused for lightwave
      pert_ = 0.0;
      break;
    case HybridWave:
      seedE_ = false;
      seedB_ = true;
      oblique_ = false;
      guideField_ = true;
      velKick_ = true;
      seedWeight_ = false;
      waveMode_ = 1;
      frac_ = 0.02;
      break;
    case ConvectionWave:
      seedE_ = false;
      seedB_ = true;
      oblique_ = false;
      guideField_ = true;
      velKick_ = false;   // the distinguishing flag vs hybridwave
      seedWeight_ = false;
      waveMode_ = 1;
      frac_ = 0.20;
      break;
    case IonAcousticWave:
      seedE_ = false;
      seedB_ = false;
      oblique_ = false;
      guideField_ = false;
      velKick_ = false;
      seedWeight_ = true;
      waveMode_ = 1;
      pert_ = 0.1;
      break;
    case Generic:
      // No behaviour by default: the user MUST supply the relevant #WAVEIC
      // sub-parameters (seedE/seedB/seedWeight, frac/pert, waveMode, ...).
      seedE_ = false;
      seedB_ = false;
      oblique_ = false;
      guideField_ = false;
      velKick_ = false;
      seedWeight_ = false;
      waveMode_ = 1;
      frac_ = 0.0;
      pert_ = 0.0;
      break;
  }
}

void WaveIC::read_param(ReadParam& param) {
  apply_preset();
  // All sub-parameters are optional: legacy tests carry them under #TESTCASE,
  // new tests carry them under #WAVEIC. read_optional returns false (without
  // aborting) when the token is absent or is the start of the next #COMMAND.
  param.read_optional("seedE", seedE_);
  param.read_optional("seedB", seedB_);
  param.read_optional("oblique", oblique_);
  param.read_optional("guideField", guideField_);
  param.read_optional("velKick", velKick_);
  param.read_optional("seedWeight", seedWeight_);
  param.read_optional("dir", dir_[0]);
  param.read_optional("dir", dir_[1]);
  param.read_optional("dir", dir_[2]);
  param.read_optional("waveLength", waveLength_);
  param.read_optional("waveMode", waveMode_);
  param.read_optional("frac", frac_);
  param.read_optional("pert", pert_);
  param.read_optional("anisoTPerpOverTPar", anisoTPerpOverTPar_);
}

void WaveIC::set_fields(PicICFields& fields) const {
  const int nLev = fields.n_lev();
  if (nLev == 0)
    return;

  for (int iLev = 0; iLev < nLev; ++iLev) {
    // Guide field Bx0 already deposited by fill_E_B_fields() from #UNIFORMSTATE.
    amrex::Real Bx0 = 1.0;
    if (guideField_) {
      MFIter mfi(fields.node_B(iLev));
      if (mfi.isValid()) {
        const Array4<Real const>& a = fields.node_B(iLev).array(mfi);
        const Box& b = mfi.validbox();
        Bx0 = a(b.smallEnd(), ix_);
      }
    }
    const amrex::Real B1 = frac_ * Bx0;
    B1_ = B1;

    const auto& prob_lo = fields.geom(iLev).ProbLo();
    const auto& dx = fields.geom(iLev).CellSize();
    const amrex::Real Lx = (fields.geom(iLev).ProbHi())[0] - prob_lo[0];
    Lx_ = Lx;

    // Phase wavenumber K.
    amrex::Real Kx, Ky, Kz;
    if (oblique_) {
      const amrex::Real kMag = 2.0 * dPI / waveLength_;
      Kx = kMag * dir_[0];
      Ky = kMag * dir_[1];
      Kz = kMag * dir_[2];
      kx_ = kMag;
    } else {
      kx_ = (Lx > 0.0) ? 2.0 * dPI * waveMode_ / Lx : 0.0;
      Kx = kx_; Ky = 0.0; Kz = 0.0;
    }

    MultiFab& nodeE = fields.node_E(iLev);
    MultiFab& nodeB = fields.node_B(iLev);
    MultiFab& centerB = fields.center_B(iLev);

    if (oblique_) {
      // LightWave: zero then write E and B on the node grid; B also on the
      // cell-centred grid. Phase uses node coordinates (i, j, k) on the node
      // grid and (i+0.5, j+0.5, k+0.5) on the centre grid, exactly as the
      // former fill_lightwaves(48.0).
      if (seedE_) nodeE.setVal(0.0);
      if (seedB_) { nodeB.setVal(0.0); centerB.setVal(0.0); }

      const amrex::Real n0 = dir_[0], n1 = dir_[1];

      if (seedE_ || seedB_) {
        for (MFIter mfi(nodeE); mfi.isValid(); ++mfi) {
          FArrayBox& fabE = nodeE[mfi];
          FArrayBox& fabB = nodeB[mfi];
          const Box& box = mfi.fabbox();
          const Array4<Real>& arrE = fabE.array();
          const Array4<Real>& arrB = fabB.array();
          ParallelFor(box, [&](int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + dx[0] * i;
            const amrex::Real y = prob_lo[1] + dx[1] * j;
            const amrex::Real z = prob_lo[2] + dx[2] * k;
            const amrex::Real phase = Kx * x + Ky * y + Kz * z;
            const amrex::Real cphi = std::cos(phase);
            const amrex::Real sphi = std::sin(phase);
            if (seedE_) {
              arrE(i, j, k, ix_) = -n1 * sphi;
              arrE(i, j, k, iy_) = n0 * sphi;
              arrE(i, j, k, iz_) = -cphi;
            }
            if (seedB_) {
              arrB(i, j, k, ix_) = -n1 * cphi;
              arrB(i, j, k, iy_) = n0 * cphi;
              arrB(i, j, k, iz_) = sphi;
            }
          });
        }
      }
      if (seedB_) {
        for (MFIter mfi(centerB); mfi.isValid(); ++mfi) {
          FArrayBox& fabcB = centerB[mfi];
          const Box& box = mfi.fabbox();
          const Array4<Real>& arrcB = fabcB.array();
          ParallelFor(box, [&](int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
            const amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
            const amrex::Real z = prob_lo[2] + dx[2] * (k + 0.5);
            const amrex::Real phase = Kx * x + Ky * y + Kz * z;
            const amrex::Real cphi = std::cos(phase);
            const amrex::Real sphi = std::sin(phase);
            arrcB(i, j, k, ix_) = -n1 * cphi;
            arrcB(i, j, k, iy_) = n0 * cphi;
            arrcB(i, j, k, iz_) = sphi;
          });
        }
      }
    } else {
      // x-aligned transverse circularly-polarized wave (hybrid / convection):
      // B = (Bx0, B1 cos kx, B1 sin kx), i.e. a pure transverse perturbation
      // delta B = (0, B1 cos kx, B1 sin kx) on the guide field Bx0, matching
      // the Alfven velocity kick in modify_particle_velocity(). This is the
      // self-consistent right-hand-circularly-polarized whistler/ion-cyclotron
      // eigenmode seed; it must NOT add any longitudinal delta Bx (a non-zero
      // Bx perturbation would make the field non-transverse and break the
      // circular polarization the velocity kick assumes). E is left as
      // deposited by fill_E_B_fields (not zeroed).
      if (seedB_) {
        nodeB.setVal(0.0);
        centerB.setVal(0.0);
        for (MFIter mfi(nodeB); mfi.isValid(); ++mfi) {
          FArrayBox& fab = nodeB[mfi];
          const Box& box = mfi.fabbox();
          const Array4<Real>& arrB = fab.array();
          ParallelFor(box, [&](int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + dx[0] * i;
            const amrex::Real cphi = std::cos(kx_ * x);
            const amrex::Real sphi = std::sin(kx_ * x);
            arrB(i, j, k, ix_) = Bx0;
            arrB(i, j, k, iy_) = B1 * cphi;
            arrB(i, j, k, iz_) = B1 * sphi;
          });
        }
        for (MFIter mfi(centerB); mfi.isValid(); ++mfi) {
          FArrayBox& fab = centerB[mfi];
          const Box& box = mfi.fabbox();
          const Array4<Real>& arrB = fab.array();
          ParallelFor(box, [&](int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
            const amrex::Real cphi = std::cos(kx_ * x);
            const amrex::Real sphi = std::sin(kx_ * x);
            arrB(i, j, k, ix_) = Bx0;
            arrB(i, j, k, iy_) = B1 * cphi;
            arrB(i, j, k, iz_) = B1 * sphi;
          });
        }
      }
    }
  }

  fields.fill_boundary_E_B();
}

void WaveIC::modify_particle_weight(ParticleICState& s) const {
  if (!seedWeight_ || Lx_ <= 0.0)
    return;
  const amrex::Real factor = 1.0 + pert_ * std::sin(kx_ * s.x);
  s.q *= (factor > 0.0) ? factor : 0.0;
}

void WaveIC::modify_particle_velocity(ParticleICState& s) const {
  if (!velKick_)
    return;
  // Matching Alfven ion velocity kick, reproducing the former HybridWave B
  // perturbation delta B = (0, B1 cos kx, B1 sin kx) on guide field B0 with
  // v_A = B0 = 1 in code units:
  //   delta U_y = -B1 cos(kx*x),  delta U_z = -B1 sin(kx*x).
  // The kick is purely transverse (y, z); the x-bulk is untouched.
  const amrex::Real cphi = std::cos(kx_ * s.x);
  const amrex::Real sphi = std::sin(kx_ * s.x);
  s.vBulk -= B1_ * cphi;
  s.wBulk -= B1_ * sphi;
}

// Bi-Maxwellian thermal seeding for the proton-cyclotron anisotropy
// instability.  The #UNIFORMSTATE isotropic sampling draws every component with
// the thermal speed v_th = sqrt(kB*T/m) where T is the #UNIFORMSTATE
// temperature.  Here T is interpreted as the PARALLEL temperature T_par, so:
//   * the component parallel to the guide field is left untouched (already v_th
//     = v_th_par), and
//   * the two perpendicular components are inflated by sqrt(T_perp/T_par).
// The guide field (and hence the parallel direction) is x in all the hybrid
// x-aligned wave presets (kx || B0 || x), so u_thermal is the parallel draw and
// v_thermal / w_thermal are the perpendicular draws.  This gives a bi-Maxwellian
// with T_perp/T_par = anisoTPerpOverTPar_ and beta_par unchanged from the
// #UNIFORMSTATE (since the parallel temperature is unchanged), which drives the
// proton-cyclotron anisotropy instability (PCAI).
void WaveIC::modify_particle_thermal_velocity(ParticleICState& s) const {
  if (anisoTPerpOverTPar_ <= 0.0)
    return;
  const amrex::Real scale = std::sqrt(anisoTPerpOverTPar_);
  // Parallel (x) direction keeps the #UNIFORMSTATE T_par thermal speed.
  s.vThermal *= scale;
  s.wThermal *= scale;
}
