#include <cstdlib>

#include <AMReX_ParReduce.H>

#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::outflow_bc(const MFIter& mfi,
                                                    const IntVect ijkGst,
                                                    const IntVect ijkPhy) {
  const int iLev = 0;

  ParticleTileType& pGst = get_particle_tile(iLev, mfi, ijkGst);

  ParticleTileType& pPhy = get_particle_tile(iLev, mfi, ijkPhy);

  AoS& phyParts = pPhy.GetArrayOfStructs();

  RealVect dxshift;
  for (int i = 0; i < nDim; ++i) {
    dxshift[i] = Geom(iLev).CellSize(i) * (ijkGst[i] - ijkPhy[i]);
  }

  Vector<ParticleType> pList;
  for (const auto& p : phyParts) {
    IntVect iv = Index(p, iLev);
    // Q: Why do we need to check if the physical domain contains the particle?
    // A: Even if tiling with tile_size=1 is used, it seems the ghost cells
    // still share the the same tile with a physical cell. Therefore, we need to
    // make sure a particle in a "physical tile" is actually inside the physical
    // domain.
    if (mfi.validbox().contains(IntVect(iv))) {
      ParticleType pNew = p;
      set_ids(pNew);
      for (int i = 0; i < nDim; ++i) {
        pNew.pos(i) = p.pos(i) + dxshift[i];
      }

      pList.push_back(pNew);
    }
  }

  // Q: Why do not push the new particles into pGst inside previous loop?
  // A: Sometimes, if not always, pPhy and pGst share the same tile. Previous
  // for-loop loops through al particles in pPhy. If we push the new particles
  // into pGst, which is the same as pPhy sometimes, the loop behavior is not
  // well defined.
  for (auto& p : pList) {
    pGst.push_back(p);
  }
}

//==========================================================
static bool particle_bc_no_injection(const ParticleBC::Type type) {
  switch (type) {
    case ParticleBC::inflow:
    case ParticleBC::outflow:
    case ParticleBC::vacuum:
    case ParticleBC::reflect:
    case ParticleBC::absorb:
      return true;
    case ParticleBC::coupled:
    case ParticleBC::periodic:
      return false;
    default:
      return true;
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::inject_particles_at_boundary() {
  timing_func("Pts::inject_particles_at_boundary");

  // Only inject nGstInject layers.
  const int nGstInject = 1;

  // Only launch particles to the base grid boundary cells. The particle moments
  // of the domain edge nodes can be corrected by calling
  // interp_from_coarse_to_fine_for_domain_edge() in order from coarest level to
  // finest level.
  int iLev = 0;

  for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
    const auto& status = cell_status(iLev)[mfi].array();
    const Box& bx = mfi.validbox();
    const IntVect bxLo = bx.smallEnd();
    const IntVect bxHi = bx.bigEnd();

    Box bxGst = bx;
    for (int iDim = 0; iDim < fi->get_fluid_dimension(); iDim++) {
      bxGst.grow(iDim, nGstInject);
    }

    // Host-only kernel: CPU particle allocation into ParticleContainer
    ParallelFor(bxGst, [&](int i, int j, int k) noexcept {
      IntVect ijk = { AMREX_D_DECL(i, j, k) };
      IntVect ijksrc;
      if (do_inject_particles_for_this_cell(bx, status, ijk, ijksrc)) {
        // A ghost cell can lie beyond more than one face (a corner); it is
        // skipped when ANY of those faces must not carry a ghost-cell
        // population -- see particle_bc_no_injection() above.
        bool skip = false;
        for (int d = 0; d < nDim && !skip; ++d) {
          if (ijk[d] < bxLo[d])
            skip = particle_bc_no_injection(
                static_cast<ParticleBC::Type>(bc.face(d, 0)));
          else if (ijk[d] > bxHi[d])
            skip = particle_bc_no_injection(
                static_cast<ParticleBC::Type>(bc.face(d, 1)));
        }
        if (skip)
          return;

        // Seed ghost cell from prescribed #INFLOW state if defined, else fall
        // back to the fluid interface.
        Vel inflowVel;
        if (fi->get_inflow_defined()) {
          const auto* iv = fi->get_inflow_vel(speciesID);
          if (iv) {
            inflowVel.tag = speciesID; // enable the userState override path
            inflowVel.nDens = iv->nDens;
            inflowVel.vth = iv->vth;
            inflowVel.vx = iv->ux;
            inflowVel.vy = iv->uy;
            inflowVel.vz = iv->uz;
          }
        }
        add_particles_cell(iLev, mfi, ijk, fi, true, IntVect(), inflowVel, -1);
      }
    });
  }
}

//==========================================================

namespace {
constexpr Real injPI = 3.14159265358979323846264338328;
constexpr Real sqpi = 1.77245385090551602729816748334;

// Draw a standard normal variate from two uniforms (Box-Muller).
inline Real inj_gaussian(Real r1, Real r2) {
  const Real rr = std::sqrt(-2.0 * std::log(std::max(r1, 1e-300)));
  return rr * std::cos(2.0 * injPI * r2);
}

// Mean inward flux of a drifting Maxwellian through a boundary face, in
// units of n * vtherm (vtherm = sqrt(2) * sigma, sigma = 1-D thermal std):
//   g(vd) = [ exp(-vd^2)/sqrt(pi) + vd * erfc(-vd) ] / 2
// with vd = (inward drift speed) / vtherm. Using erfc(-vd) == 1 + erf(vd)
// avoids catastrophic cancellation for negative outward drift (vd < 0).
inline Real inj_mean_inward_flux(Real vd) {
  return 0.5 * (std::exp(-vd * vd) / sqpi + vd * std::erfc(-vd));
}

// Fast speed sampler for the flux-weighted half-space Maxwellian
//   f(w) \propto w * exp(-(w - vd)^2), w >= 0.
// Precomputes face-level constants and a 64-point LUT once per face to
// eliminate costly 60-step bisections per particle, converging in 2-3
// Newton-Raphson steps.
struct InflowSpeedSampler {
  Real vd{ 0.0 };
  Real e0{ 1.0 };
  Real erfc_mvd{ 1.0 };
  Real twoZ{ 1.0 };
  Real wHi{ 8.0 };

  static constexpr int LUT_SIZE = 64;
  Real lut[LUT_SIZE + 1];

  void init(Real vd_in) {
    vd = vd_in;
    e0 = std::exp(-vd * vd);
    erfc_mvd = std::erfc(-vd);
    twoZ = e0 + vd * sqpi * erfc_mvd;
    wHi = std::max(vd, 0.0) + 8.0;

    lut[0] = 0.0;
    Real wLo = 0.0;
    for (int k = 1; k < LUT_SIZE; ++k) {
      const Real rTarget = static_cast<Real>(k) / LUT_SIZE;
      const Real target = 0.5 * rTarget * twoZ;
      Real lo = wLo;
      Real hi = wHi;
      for (int it = 0; it < 30; ++it) {
        Real mid = 0.5 * (lo + hi);
        Real Fmid = 0.5 * ((e0 - std::exp(-(mid - vd) * (mid - vd))) +
                           vd * sqpi * (erfc_mvd - std::erfc(mid - vd)));
        if (Fmid < target)
          lo = mid;
        else
          hi = mid;
      }
      lut[k] = 0.5 * (lo + hi);
      wLo = lut[k];
    }
    lut[LUT_SIZE] = wHi;
  }

  inline Real draw(Real r) const {
    const Real target = 0.5 * r * twoZ;
    const Real rIdx = r * LUT_SIZE;
    const int idx = std::min(std::max(static_cast<int>(rIdx), 0), LUT_SIZE - 1);
    const Real frac = rIdx - idx;
    Real w = lut[idx] + frac * (lut[idx + 1] - lut[idx]);

    // 2-3 Newton-Raphson polish iterations
    for (int it = 0; it < 3; ++it) {
      const Real diff = w - vd;
      const Real exp_term = std::exp(-diff * diff);
      const Real Fw =
          0.5 * ((e0 - exp_term) + vd * sqpi * (erfc_mvd - std::erfc(diff)));
      const Real dF = w * exp_term;
      if (std::abs(dF) < 1e-300)
        break;
      Real wNew = w - (Fw - target) / dF;
      if (wNew < 0.0 || wNew > wHi)
        wNew = 0.5 * (std::max(static_cast<Real>(0.0), w) + wHi);
      if (std::abs(wNew - w) < 1e-13 * (1.0 + w))
        break;
      w = wNew;
    }
    return w;
  }
};
} // namespace

//==========================================================
template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::inject_flux_at_inflow_faces(Real dt) {
  timing_func("Pts::inject_flux_at_inflow_faces");

  if (dt <= 0)
    return;

  // Any inflow face for this species?
  bool hasInflow = false;
  for (int iDim = 0; iDim < nDim && !hasInflow; ++iDim)
    hasInflow = (bc.lo[iDim] == ParticleBC::inflow) ||
                (bc.hi[iDim] == ParticleBC::inflow);
  if (!hasInflow)
    return;

  // Prescribed upstream state (#INFLOW, code units).  vth is the 1-D
  // thermal std sigma = sqrt(kT/m); the Maxwellian is written with
  // vtherm = sqrt(2)*sigma.
  const auto* inv =
      fi->get_inflow_defined() ? fi->get_inflow_vel(speciesID) : nullptr;
  if (inv == nullptr || inv->nDens <= 0)
    return;

  const Real nDens = inv->nDens;
  const Real sigma = inv->vth;
  const Real vtherm = std::sqrt(2.0) * sigma;
  const Real uIn[3] = { inv->ux, inv->uy, inv->uz };

  // Macroparticle weight identical to the interior cell particles
  // (add_particles_cell convention): qp = qomSign * n * V / nppc.
  const Real q = qomSign * dx[0].product() * nDens / product(nPartPerCell);
  const int nppc = product(nPartPerCell);

  // Inject on the base level only (same policy as
  // inject_particles_at_boundary).
  const int iLev = 0;
  const auto& geom = Geom(iLev);
  const Box& domain = geom.Domain();
  const Real* probLo = geom.ProbLo();
  const Real* probHi = geom.ProbHi();
  const Real* cellSize = geom.CellSize();

  // Cache speed samplers for each inflow face (iDim, side)
  InflowSpeedSampler speedSamplers[3][2];
  bool samplerInit[3][2] = { { false, false },
                             { false, false },
                             { false, false } };

  for (MFIter mfi = MakeMFIter(iLev, false); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.validbox();
    ParticleTileType& particles = get_particle_tile(iLev, mfi);

    for (int iDim = 0; iDim < nDim; ++iDim) {
      const int t1 = (iDim + 1) % nDim; // first transverse grid direction
      const int t2 = (nDim > 2) ? (iDim + 2) % nDim
                                : 0;     // second transverse grid direction
      const int trans1 = (iDim + 1) % 3; // first transverse velocity direction
      const int trans2 = (iDim + 2) % 3; // second transverse velocity direction

      for (int side = 0; side < 2; ++side) {
        const bool isHi = (side == 1);
        const int faceBc = isHi ? bc.hi[iDim] : bc.lo[iDim];
        if (faceBc != ParticleBC::inflow)
          continue;

        // This tile must touch the global domain edge on this face.
        const int domEdge = isHi ? domain.bigEnd(iDim) : domain.smallEnd(iDim);
        const int tileEdge = isHi ? bx.bigEnd(iDim) : bx.smallEnd(iDim);
        if (tileEdge != domEdge)
          continue;

        const Real dxn = cellSize[iDim];
        // Inward-pointing direction in the outward-normal coordinate:
        // hi face -> -1 (inward is -d), lo face -> +1.
        const Real inward = isHi ? -1.0 : 1.0;

        // Inward drift speed = bulk velocity dotted with the inward unit
        // normal (uOut * inward; > 0 when plasma flows into the domain).
        const Real uOut = uIn[iDim];
        const Real vd = (sigma > 0) ? (uOut * inward / vtherm) : 1.0e30;

        if (vtherm > 0 && !samplerInit[iDim][side]) {
          speedSamplers[iDim][side].init(vd);
          samplerInit[iDim][side] = true;
        }

        // Mean influx per boundary-transverse cell per step, in
        // macroparticles (Hybrid-VPIC shock deck, bright() accumulator):
        //   dn = nppc * vtherm * g(vd) * dt / dx
        const Real fluxRate =
            (vtherm > 0) ? nppc * vtherm * inj_mean_inward_flux(vd) * dt / dxn
                         : nppc * (uOut * inward) * dt / dxn;
        if (fluxRate <= 0.0)
          continue;

        // Transverse cell range for this tile's face
        const int lo1 = bx.smallEnd(t1), hi1 = bx.bigEnd(t1);
        const int lo2 = (nDim > 2) ? bx.smallEnd(t2) : 0;
        const int hi2 = (nDim > 2) ? bx.bigEnd(t2) : 0;
        const int numCells = (hi1 - lo1 + 1) * (hi2 - lo2 + 1);

        // Tile-face accumulator: packs (iLev, iDim, side, tileLocalIndex).
        // Accumulating over the tile face eliminates the artificial coherent
        // "pulsing sheets" and replaces 10^4 map lookups per step with 1
        // scalar.
        const int64_t tileKey =
            ((((int64_t)iLev * 3 + iDim) * 2 + side) << 40) |
            (((int64_t)mfi.index()) << 16) |
            (int64_t)(mfi.LocalTileIndex() & 0xFFFF);
        Real& acc = injectFluxAcc[tileKey];
        acc += numCells * fluxRate;
        int nInject = static_cast<int>(acc);
        if (nInject <= 0) {
          if (acc < 0)
            acc = 0;
          continue;
        }
        acc -= nInject;

        const Real span1 = (hi1 - lo1 + 1) * cellSize[t1];
        const Real base1 =
            probLo[t1] + (lo1 - domain.smallEnd(t1)) * cellSize[t1];
        const Real span2 = (nDim > 2) ? (hi2 - lo2 + 1) * cellSize[t2] : 0.0;
        const Real base2 =
            (nDim > 2) ? probLo[t2] + (lo2 - domain.smallEnd(t2)) * cellSize[t2]
                       : 0.0;
        const Real xFace = isHi ? probHi[iDim] : probLo[iDim];

        for (int np = 0; np < nInject; ++np) {
          // Velocity: inward normal speed from speed sampler; transverse
          // components from paired Box-Muller. All 3 velocity components
          // are sampled so 2D3V (AMREX_SPACEDIM=2) has correct out-of-plane
          // velocity vz and thermal pressure.
          Real wIn;
          if (vtherm > 0) {
            wIn = vtherm * speedSamplers[iDim][side].draw(randNum());
          } else {
            wIn = uOut * inward; // cold beam
          }

          // Fractional ingress advancement: particles crossed the face at
          // random times t' in [0, dt], so at dt they have penetrated
          // distance wIn * dt * randNum(). Bound strictly inside the cell.
          const Real pDist =
              (wIn > 0)
                  ? std::min(1.0e-3 * dxn + wIn * dt * randNum(), 0.999 * dxn)
                  : 1.0e-3 * dxn;

          // Position: on the boundary face, nudged into the cell according to
          // sub-step ingress distance; uniform in the transverse directions.
          RealVect pos;
          pos[iDim] = xFace + inward * pDist;
          pos[t1] = base1 + randNum() * span1;
          if (nDim > 2) {
            pos[t2] = base2 + randNum() * span2;
          }

          Real vel[3] = { 0.0, 0.0, 0.0 };
          vel[iDim] = inward * wIn;
          if (sigma > 0) {
            const Real r1 = randNum();
            const Real r2 = randNum();
            const Real R = std::sqrt(-2.0 * std::log(std::max(r1, 1e-300)));
            const Real phi = 2.0 * injPI * r2;
            vel[trans1] = uIn[trans1] + sigma * (R * std::cos(phi));
            vel[trans2] = uIn[trans2] + sigma * (R * std::sin(phi));
          } else {
            vel[trans1] = uIn[trans1];
            vel[trans2] = uIn[trans2];
          }

          auto p = make_particle();
          set_ids(p);
          for (int d = 0; d < nDim; ++d)
            p.pos(d) = pos[d];
          p.rdata(iup_) = vel[ix_];
          p.rdata(ivp_) = vel[iy_];
          p.rdata(iwp_) = vel[iz_];
          p.rdata(iqp_) = q;

          if (NStructInt > iRecordCount_) {
            p.idata(iRecordCount_) = 0;
          }

          particles.push_back(p);
        }
      }
    }
  }
}

template <int NStructReal, int NStructInt>
bool Particles<NStructReal, NStructInt>::do_inject_particles_for_this_cell(
    const Box& bx, const Array4<const int>& status, const IntVect ijk,
    IntVect& ijksrc) {

  // This cell should be a boundary cell at least.
  if (!bit::is_lev_boundary(status(ijk)))
    return false;

  for (int iloop = 1; iloop <= 3; iloop++) {
    // iloop==1: loop through faces;
    // iloop==2: loop through edges;
    // iloop==3: loop through corners;
    for (int di = -1; di <= 1; di++)
      for (int dj = -1; dj <= 1; dj++)
        for (int dk = -1; dk <= 1; dk++) {
          const int sum = std::abs(di) + std::abs(dj) + std::abs(dk);
          if (iloop != sum)
            continue;

          IntVect ijk1 = ijk + IntVect{ AMREX_D_DECL(di, dj, dk) };
          if (!bit::is_lev_boundary(status(ijk1))) {
            // The first neighbor cell that is NOT a boundary cell.
            if (bx.contains(ijk1)) {
              ijksrc = ijk1;
              return true;
            } else {
              return false;
            }
          }
        }
  }
  Abort("do_inject_particles_for_this_cell:something is wrong!");
  return false; // to suppress compilation warning.
}

// Explicit template instantiations.
template class Particles<nPicPartReal, nPicPartInt>;
template class Particles<nPTPartReal, nPTPartInt>;
