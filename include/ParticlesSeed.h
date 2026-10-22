#ifndef _PARTICLES_SEED_H_
#define _PARTICLES_SEED_H_

// Creation of macroparticles from the fluid state, cell by cell, on the CPU
// or the GPU (see Particles::add_particles_cells).
//
// The particles of a cell are drawn from the cell's own random number stream
// (RandNumHD seeded from the cell index, the cycle and the species), in a
// fixed order, by one thread. So a GPU run creates exactly the particles a
// CPU run creates (up to the round-off of the math library), whatever the
// domain decomposition, and the particle IDs follow the order of the cells.

#include <AMReX_GpuAsyncArray.H>
#include <AMReX_Scan.H>

#include "ICHooks.h"
#include "Particles.h"
#include "WaveBC.h"

// Everything the seeding of a cell needs. It is copied to device memory once
// per call (it is too large for a kernel argument).
struct CellSeeder {
  int iLev = 0;
  int speciesID = 0;
  amrex::Real plo[3] = { 0, 0, 0 };
  amrex::Real dx[3] = { 0, 0, 0 };
  amrex::Real vol = 0;

  // Particles per cell: the default, or the target PPC of the cell.
  amrex::IntVect nPPC0;
  bool useTargetPPC = false;
  amrex::Array4<int const> targetPPC;

  // No particle is created in the cells of the absorbing inner body.
  bool useBody = false;
  DomainView domain;

  // Random seed of a cell (see Particles::calc_random_seed).
  int nCell[3] = { 1, 1, 1 };
  int iCycle = 0;

  amrex::Real qomSign = 1, charge = 0, mass = 1, vacuum = 0, dt = 1;
  bool doVacuumLimit = false;
  bool randomLocation = true;
  bool aniso = false;

  // User-defined state instead of the fluid state (test particles, inflow).
  bool userState = false;
  Vel tpVel;

  FluidView fv;
  ICHooks hooks;

  // Velocity kick of the wave boundary in the ghost cells next to a wave face.
  bool useWaveKick = false;
  WaveKickView waveKick;
  amrex::Real tNow = 0;
  bool isWaveFace[6] = { false, false, false, false, false, false };
  int validLo[3] = { 0, 0, 0 };
  int validHi[3] = { 0, 0, 0 };

  AMREX_GPU_HOST_DEVICE long random_seed(const amrex::IntVect& ijk,
                                         const amrex::IntVect& nPPC) const {
    const int nRandom = 7;
    const int nxcg = nCell[ix_];
    const int nycg = nCell[iy_];
    const int nzcg = nDim > 2 ? nCell[iz_] : 1;
    const int i = ijk[0];
    const int j = ijk[1];
    const int k = nDim > 2 ? ijk[2] : 0;
    const long nCellOffset = static_cast<long>(nxcg) * nycg * nzcg * iCycle +
                             static_cast<long>(nycg) * nzcg * i +
                             static_cast<long>(nzcg) * j + k;
    return static_cast<long>(speciesID + 3) * nRandom * product(nPPC) *
           nCellOffset;
  }

  // The integer n-th root of x (n = 2 or 3), robust to the round-off of
  // std::sqrt / std::cbrt.
  AMREX_GPU_HOST_DEVICE static int int_root(int x, int n) {
    if (x <= 0)
      return 0;
    int r = static_cast<int>(n == 3 ? std::cbrt(static_cast<double>(x))
                                    : std::sqrt(static_cast<double>(x)));
    auto pw = [n](long a) { return n == 3 ? a * a * a : a * a; };
    while (r > 0 && pw(r) > x)
      --r;
    while (pw(r + 1) <= x)
      ++r;
    return r;
  }

  // Create the particles of cell ijk; ppc = 0 selects the default number per
  // cell. emit(xyz, u, v, w, q) is called for every particle, in order.
  template <class Emit>
  AMREX_GPU_HOST_DEVICE void operator()(const amrex::IntVect& ijk,
                                        const amrex::IntVect& ppc,
                                        Emit& emit) const {
    amrex::IntVect nPPC = nPPC0;
    if (!(ppc == 0)) {
      nPPC = ppc;
    }

    if (nPPC == 0)
      return;

    // The test uses the cell center, i.e., the same cell-based staircase as
    // the particle absorption and the field mask.
    if (useBody) {
      amrex::Real cellCenter[AMREX_SPACEDIM];
      for (int iDim = 0; iDim < AMREX_SPACEDIM; iDim++)
        cellCenter[iDim] = 0.0;
      for (int iDim = 0; iDim < nDim; iDim++)
        cellCenter[iDim] = (ijk[iDim] + 0.5) * dx[iDim] + plo[iDim];

      if (domain.is_inside_body(cellCenter))
        return;
    }

    if (useTargetPPC) {
      for (int i = 0; i < nDim; ++i)
        nPPC[i] = int_root(targetPPC(ijk), nDim > 2 ? 3 : 2);
    }

    RandNumHD randNum(random_seed(ijk, nPPC));

    const int npcel = product(nPPC);
    const amrex::Real vol2Npcel = qomSign * vol / npcel;

    int icount = 0;
    const int kmax = nDim > 2 ? nPPC[iz_] : 1;

    for (int ii = 0; ii < nPPC[ix_]; ++ii)
      for (int jj = 0; jj < nPPC[iy_]; ++jj)
        for (int kk = 0; kk < kmax; ++kk) {
          amrex::RealVect xyz, xyz0;

          amrex::IntVect ijk0 = { AMREX_D_DECL(ii, jj, kk) };

          for (int iDim = 0; iDim < nDim; iDim++) {
            xyz[iDim] = (ijk0[iDim] + randNum()) * (dx[iDim] / nPPC[iDim]) +
                        ijk[iDim] * dx[iDim] + plo[iDim];

            // If the particle weight is sampled in a random location, the sum
            // of particle mass is NOT the same as the integral of the grid
            // density. It is more convenient for debugging if mass is exactly
            // conserved. For a production run, it makes little difference.
            xyz0[iDim] = (ijk0[iDim] + 0.5) * (dx[iDim] / nPPC[iDim]) +
                         ijk[iDim] * dx[iDim] + plo[iDim];
          }

          if (nDim == 2) {
            // For comparison with the 3D case only.
            randNum();
          }

          if (!randomLocation) {
            xyz = xyz0;
          }

          const amrex::Real zp = (nDim > 2) ? xyz[2] : 0.0;

          const amrex::Real nDens =
              (userState && tpVel.nDens >= 0.0)
                  ? tpVel.nDens
                  : fv.get_number_density(fv.d, xyz0, speciesID);

          if (doVacuumLimit && nDens * dt < vacuum)
            continue;

          amrex::Real q = vol2Npcel * nDens;

          // Per-particle weight modification of the initial condition; runs
          // before the q != 0 guard so an IC can suppress a particle.
          if (hooks.modifiesWeights) {
            ParticleICState pics;
            pics.iLev = iLev;
            pics.iSpec = speciesID;
            pics.iCount = icount;
            pics.nPerCell = npcel;
            pics.x = xyz[0];
            pics.y = xyz[1];
            pics.z = zp;
            pics.q = q;
            hooks.modify_weight(pics);
            q = pics.q;
          }

          if (q == 0)
            continue;

          amrex::Real u, v, w;
          const amrex::Real rand1 = randNum();
          const amrex::Real rand2 = randNum();
          const amrex::Real rand3 = randNum();
          const amrex::Real rand4 = randNum();

          const amrex::Real uth = (userState ? tpVel.vth : -1);

          if (aniso) {
            fv.set_particle_uth_aniso(fv.d, xyz, &u, &v, &w, rand1, rand2,
                                      rand3, rand4, speciesID, uth, uth);
          } else {
            fv.set_particle_uth_iso(fv.d, xyz, &u, &v, &w, rand1, rand2, rand3,
                                    rand4, speciesID, uth);
          }

          // Per-particle thermal-velocity override of the initial condition
          // (e.g. a bi-Maxwellian); runs before the bulk velocity is added.
          if (hooks.modifiesThermalVelocity) {
            ParticleICState pics;
            pics.iLev = iLev;
            pics.iSpec = speciesID;
            pics.iCount = icount;
            pics.nPerCell = npcel;
            pics.x = xyz[0];
            pics.y = xyz[1];
            pics.z = zp;
            pics.uThermal = u;
            pics.vThermal = v;
            pics.wThermal = w;
            pics.rand[0] = rand1;
            pics.rand[1] = rand2;
            pics.rand[2] = rand3;
            pics.rand[3] = rand4;
            hooks.modify_thermal_velocity(pics);
            u = pics.uThermal;
            v = pics.vThermal;
            w = pics.wThermal;
          }

          amrex::Real uBulk =
              userState ? tpVel.vx : fv.get_ux(fv.d, xyz, speciesID);
          amrex::Real vBulk =
              userState ? tpVel.vy : fv.get_uy(fv.d, xyz, speciesID);
          amrex::Real wBulk =
              userState ? tpVel.vz : fv.get_uz(fv.d, xyz, speciesID);

          // Per-particle velocity / weight modification of the initial
          // condition (e.g. a beam); runs after the bulk velocity is computed
          // and before it is added to the thermal velocity.
          if (hooks.modifiesVelocities) {
            ParticleICState pics;
            pics.iLev = iLev;
            pics.iSpec = speciesID;
            pics.iCount = icount;
            pics.nPerCell = npcel;
            pics.x = xyz[0];
            pics.y = xyz[1];
            pics.z = zp;
            pics.uBulk = uBulk;
            pics.vBulk = vBulk;
            pics.wBulk = wBulk;
            pics.qScale = 1.0;
            pics.charge = charge;
            pics.temperature = fv.get_uth_iso(fv.d, xyz, speciesID);
            pics.temperature *= pics.temperature * mass;
            hooks.modify_velocity(pics);
            uBulk = pics.uBulk;
            vBulk = pics.vBulk;
            wBulk = pics.wBulk;
            q *= pics.qScale;
          }

          // Wave boundary: add the velocity kick to particles in a wave cell.
          if (useWaveKick) {
            const bool onWaveX = (isWaveFace[0] && ijk[ix_] < validLo[0]) ||
                                 (isWaveFace[1] && ijk[ix_] > validHi[0]);
            const bool onWaveY = (isWaveFace[2] && ijk[iy_] < validLo[1]) ||
                                 (isWaveFace[3] && ijk[iy_] > validHi[1]);
            const bool onWaveZ =
                (nDim > 2) && ((isWaveFace[4] && ijk[iz_] < validLo[2]) ||
                               (isWaveFace[5] && ijk[iz_] > validHi[2]));
            if (onWaveX || onWaveY || onWaveZ) {
              const amrex::Real ppos[3] = { xyz[0], xyz[1], zp };
              amrex::Real dvx = 0, dvy = 0, dvz = 0;
              waveKick.kick(ppos, tNow, dvx, dvy, dvz);
              uBulk += dvx;
              vBulk += dvy;
              wBulk += dvz;
            }
          }
          u += uBulk;
          v += vBulk;
          w += wBulk;

          // Convert 'density changing rate' to 'density' if necessary.
          emit(xyz, u, v, w, q * dt);
          icount++;
        }
  }
};

// Counts the particles of a cell.
struct SeedCounter {
  int n = 0;
  AMREX_GPU_HOST_DEVICE void operator()(const amrex::RealVect&, amrex::Real,
                                        amrex::Real, amrex::Real, amrex::Real) {
    ++n;
  }
};

// A block of consecutive particle IDs, with the wrap-around of
// Particles::get_next_ids: an id above idMax restarts from 1 with the next
// supID.
struct ParticleIDBlock {
  long id0 = 1;
  int supID0 = 1;
  static constexpr long idMax = 2147483647L;

  AMREX_GPU_HOST_DEVICE void get(long g, int& id, int& supID) const {
    const long raw = id0 + g;
    if (raw <= idMax) {
      id = static_cast<int>(raw);
      supID = supID0;
    } else {
      // The first id past idMax becomes 1, and the counter then restarts
      // from 1 (see get_next_ids), so the ids after the wrap are 1, 1, 2, ...
      const long k = raw - idMax;
      id = static_cast<int>(k > 1 ? k - 1 : 1);
      supID = supID0 + 1;
    }
  }
};

// Writes the particles of a cell, starting at dst, with IDs from the block
// starting at the global index g0.
template <class ParticleType, int NStructInt> struct SeedWriter {
  ParticleType* dst = nullptr;
  ParticleIDBlock ids;
  long g0 = 0;
  int cpu = 0;
  int n = 0;

  AMREX_GPU_HOST_DEVICE void operator()(const amrex::RealVect& xyz,
                                        amrex::Real u, amrex::Real v,
                                        amrex::Real w, amrex::Real q) {
    ParticleType p{};
    int id, supID;
    ids.get(g0 + n, id, supID);
    p.id() = id;
    p.cpu() = cpu;
    p.idata(0) = supID; // iSupID_
    for (int iDim = 0; iDim < nDim; iDim++) {
      p.pos(iDim) = xyz[iDim];
    }
    p.rdata(0) = u; // iup_
    p.rdata(1) = v; // ivp_
    p.rdata(2) = w; // iwp_
    p.rdata(3) = q; // iqp_
    if (NStructInt > 1) {
      // For test particle only: the record counter.
      p.idata(NStructInt > 1 ? 1 : 0) = 0;
    }
    dst[n] = p;
    ++n;
  }
};

template <int NStructReal, int NStructInt>
CellSeeder Particles<NStructReal, NStructInt>::make_cell_seeder(
    int iLev, const amrex::MFIter& mfi, const FluidInterface* interface,
    bool doVacuumLimit, const Vel& tpVel, amrex::Real dt) const {
  CellSeeder s;
  s.iLev = iLev;
  s.speciesID = speciesID;
  for (int d = 0; d < nDim; ++d) {
    s.plo[d] = plo[iLev][d];
    s.dx[d] = dx[iLev][d];
  }
  s.vol = dx[iLev].product();
  s.nPPC0 = nPartPerCell;
  s.useTargetPPC = isTargetPPCDefined && !isFake2D;
  if (s.useTargetPPC)
    s.targetPPC = target_PPC(iLev).const_array(mfi);
  s.useBody = (grid != nullptr && grid->use_body());
  if (grid != nullptr)
    s.domain = grid->domain_view();
  const amrex::IntVect nCell = Geom(iLev).Domain().size();
  for (int d = 0; d < AMREX_SPACEDIM; ++d)
    s.nCell[d] = nCell[d];
  s.iCycle = tc->get_cycle();
  s.qomSign = qomSign;
  s.charge = charge;
  s.mass = mass;
  s.vacuum = vacuum;
  // If dt > 0, the 'density' of the interface is a density changing rate.
  s.dt = (dt <= 0) ? 1 : dt;
  s.doVacuumLimit = doVacuumLimit;
  s.randomLocation = isParticleLocationRandom;
  s.aniso = !is_neutral() && interface->get_UseAnisoP() &&
            (speciesID > 0 || interface->get_useElectronFluid());
  s.userState = (tpVel.tag == speciesID);
  s.tpVel = tpVel;
  s.fv = interface->fluid_view(mfi, iLev);
  if (ic_) {
    s.hooks = ic_->particle_hooks();
    if ((ic_->modifies_weights() || ic_->modifies_velocities() ||
         ic_->modifies_thermal_velocity()) &&
        s.hooks.type == ICHooks::none)
      amrex::Abort("The initial condition '" + ic_->name() +
                   "' modifies the particles but does not provide "
                   "particle_hooks().");
  }
  if (waveBC != nullptr) {
    s.tNow = tc ? tc->get_time() : 0.0;
    s.waveKick = make_wave_kick_view(*waveBC, s.tNow);
    s.useWaveKick = s.waveKick.comps.size() > 0;
    for (int i = 0; i < 6; ++i)
      s.isWaveFace[i] = isWaveFace[i];
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      s.validLo[d] = mfi.validbox().smallEnd(d);
      s.validHi[d] = mfi.validbox().bigEnd(d);
    }
  }
  return s;
}

// Reserve n consecutive particle IDs (see get_next_ids).
template <int NStructReal, int NStructInt>
ParticleIDBlock Particles<NStructReal, NStructInt>::reserve_particle_ids(
    long n) {
  ParticleIDBlock ids;
  ids.id0 = ParticleType::NextID();
  ids.supID0 = supID;
  if (ids.id0 > ParticleIDBlock::idMax) {
    // As in get_next_ids: restart from 1 with the next supID.
    ids.id0 = 1;
    supID++;
    ids.supID0 = supID;
  }
  const long last = ids.id0 + n - 1;
  if (last <= ParticleIDBlock::idMax) {
    ParticleType::NextID(last + 1);
  } else {
    // The block wraps: the IDs after the wrap restart from 1.
    const long k = last - ParticleIDBlock::idMax;
    supID++;
    ParticleType::NextID(k);
  }
  return ids;
}

// Create the particles of the cells of `box` (part of the FAB mfi of level
// iLev, ghost cells allowed) for which rule(ijk, ppc) returns true; rule may
// set ppc, the number of particles per cell (0 = default). This runs on the
// device in a GPU build. With one particle tile per cell (CPU), the particles
// of a cell go to its tile; otherwise all go to the tile of the box, in cell
// order.
template <int NStructReal, int NStructInt>
template <class Rule>
void Particles<NStructReal, NStructInt>::add_particles_cells(
    int iLev, const amrex::MFIter& mfi, const amrex::Box& box, const Rule& rule,
    const FluidInterface* interface, bool doVacuumLimit, const Vel& tpVel,
    amrex::Real dt) {
  using namespace amrex;

  if (box.isEmpty())
    return;

  const CellSeeder seeder =
      make_cell_seeder(iLev, mfi, interface, doVacuumLimit, tpVel, dt);
  Gpu::AsyncArray<CellSeeder> seederStore(&seeder, 1);
  const CellSeeder* const sp = seederStore.data();

  const int nCellBox = static_cast<int>(box.numPts());
  const Dim3 lo = lbound(box);
  const Dim3 len = length(box);
  auto cell_of = [=] AMREX_GPU_HOST_DEVICE(int l) {
    const int k = l / (len.x * len.y);
    const int j = (l - k * len.x * len.y) / len.x;
    const int i = l - (k * len.y + j) * len.x;
    return IntVect(AMREX_D_DECL(lo.x + i, lo.y + j, lo.z + k));
  };

  // 1. Count the particles of each cell.
  Gpu::DeviceVector<int> counts(nCellBox);
  int* const cnt = counts.data();
  ParallelFor(nCellBox, [=] AMREX_GPU_DEVICE(int l) {
    const IntVect ijk = cell_of(l);
    IntVect ppc(0);
    int n = 0;
    if (rule(ijk, ppc)) {
      SeedCounter counter;
      (*sp)(ijk, ppc, counter);
      n = counter.n;
    }
    cnt[l] = n;
  });

  // 2. Position of the first particle of each cell, in cell order.
  Gpu::DeviceVector<int> offsets(nCellBox);
  int* const off = offsets.data();
  const int total = Scan::ExclusiveSum(nCellBox, cnt, off, Scan::retSum);
  if (total == 0)
    return;

  const ParticleIDBlock ids = reserve_particle_ids(total);
  const int cpu = ParallelDescriptor::MyProc();
  using Writer = SeedWriter<ParticleType, NStructInt>;

  // 3. Create the particles.
  if (do_tiling) {
    // One tile per cell: fill the tiles of the cells one after the other.
    // Only a CPU build uses this layout, so it runs on the host.
    for (int l = 0; l < nCellBox; ++l) {
      if (cnt[l] == 0)
        continue;
      const IntVect ijk = cell_of(l);
      IntVect ppc(0);
      rule(ijk, ppc);
      ParticleTileType& tile = get_particle_tile(iLev, mfi, ijk);
      const int nOld = tile.numParticles();
      tile.resize(nOld + cnt[l]);
      Writer writer{ tile.GetArrayOfStructs()().data() + nOld, ids, off[l],
                     cpu };
      seeder(ijk, ppc, writer);
    }
  } else {
    ParticleTileType& tile = get_particle_tile(iLev, mfi, box.smallEnd());
    const int nOld = tile.numParticles();
    tile.resize(nOld + total);
    ParticleType* const dst = tile.GetArrayOfStructs()().data() + nOld;
    ParallelFor(nCellBox, [=] AMREX_GPU_DEVICE(int l) {
      if (cnt[l] == 0)
        return;
      const IntVect ijk = cell_of(l);
      IntVect ppc(0);
      rule(ijk, ppc);
      Writer writer{ dst + off[l], ids, off[l], cpu };
      (*sp)(ijk, ppc, writer);
    });
  }
}

#endif
