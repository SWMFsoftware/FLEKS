#include "Particles.h"

using namespace amrex;

Particles::Particles(const Geometry& geom, const DistributionMapping& dm,
                     const BoxArray& ba, const int speciesIDIn,
                     const Real chargeIn, const Real massIn,
                     const IntVect& nPartPerCellIn)
    : ParticleContainer<4, 0, 0, 0>(geom, dm, ba),
      speciesID(speciesIDIn),
      charge(chargeIn),
      mass(massIn),
      nPartPerCell(nPartPerCellIn) {}

void Particles::add_particles_cell(const MFIter& mfi,
                                   const FluidPicInterface& fluidInterface,
                                   int iBlock, int i, int j, int k, int loi,
                                   int loj, int lok) {
  int ig, jg, kg, nxcg, nycg, nzcg, iCycle, npcel, nRandom = 7;
  // Why +1? for comparison with iPIC3D.-----
  ig = i + 1;
  jg = j + 1;
  kg = k;
  if (fluidInterface.getnDim() > 2)
    kg++; // just for comparison with iPIC3D;
  //----------------------------------------

  nxcg = fluidInterface.getFluidNxc();
  nycg = fluidInterface.getFluidNyc();
  nzcg = fluidInterface.getFluidNzc();
  iCycle = fluidInterface.getCycle();
  npcel = nPartPerCell[ix_] * nPartPerCell[iy_] * nPartPerCell[iz_];
  // What if the seed overflow?
  const long seed =
      (speciesID + 3) * nRandom * npcel *
      (nxcg * nycg * nzcg * iCycle + nycg * nzcg * ig + nzcg * jg + kg);
  randNum.set_seed(seed);
  Print() << " ns = " << speciesID << " iCycle = " << iCycle
          << " npcel = " << npcel << " nxcg =" << nxcg << " nycg = " << nycg
          << " nzcg = " << nzcg << " i = " << i << "  j = " << j << " k = " << k
          << " seed = " << seed << std::endl;

  double x, y, z; // Particle location.

  auto dx = Geom(0).CellSize();
  auto plo = Geom(0).ProbLo();

  const Real invVol = 1.0 / dx[ix_] / dx[iy_] / dx[iz_];

  Print() << "dxz = " << dx[iz_] << std::endl;

  const int lev = 0;
  auto& particles =
      GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

  // loop over particles inside grid cell i, j, k
  for (int ii = 0; ii < nPartPerCell[ix_]; ii++)
    for (int jj = 0; jj < nPartPerCell[iy_]; jj++)
      for (int kk = 0; kk < nPartPerCell[iz_]; kk++) {

        x = (ii + randNum()) * (dx[ix_] / nPartPerCell[ix_]) + i * dx[ix_];
        y = (jj + randNum()) * (dx[iy_] / nPartPerCell[iy_]) + j * dx[iy_];
        z = (kk + randNum()) * (dx[iz_] / nPartPerCell[iz_]) + k * dx[iz_];

        double q =
            (charge / mass / fabs(charge / mass)) *
            (fluidInterface.getPICRhoNum(iBlock, x, y, z, speciesID) / npcel) *
            invVol;
        if (q != 0) {
          double rand;
          double u, v, w;
          double rand1 = randNum();
          double rand2 = randNum();
          double rand3 = randNum();
          double rand4 = randNum();

          if (fluidInterface.getUseAnisoP() &&
              (speciesID > 0 || fluidInterface.get_useElectronFluid())) {
            fluidInterface.setPICAnisoUth(iBlock, x, y, z, &u, &v, &w, rand1,
                                          rand2, rand3, rand4, speciesID);
          } else {
            fluidInterface.setPICIsoUth(iBlock, x, y, z, &u, &v, &w, rand1,
                                        rand2, rand3, rand4, speciesID);
          }
          u += fluidInterface.getPICUx(iBlock, x, y, z, speciesID);
          v += fluidInterface.getPICUy(iBlock, x, y, z, speciesID);
          w += fluidInterface.getPICUz(iBlock, x, y, z, speciesID);

          ParticleType p;
          p.id() = ParticleType::NextID();
          p.cpu() = ParallelDescriptor::MyProc();
          p.pos(ix_) = x + plo[ix_];
          p.pos(iy_) = y + plo[iy_];
          p.pos(iz_) = z + plo[iz_];
          p.rdata(iup_) = u;
          p.rdata(ivp_) = v;
          p.rdata(iwp_) = w;
          p.rdata(iqp_) = q;
          particles.push_back(p);

          Print() << " i = " << i << " j = " << j << " k = " << k
                  << " x = " << p.pos(ix_) << " y = " << p.pos(iy_)  << " z = " << p.pos(iz_) 
                  << " q = " << p.rdata(iqp_) << " u = " << p.rdata(iup_) << " v = " << p.rdata(ivp_)
                  << " w = " << p.rdata(iwp_) << " charge = " << charge
                  << " mass = " << mass << std::endl;

          // Set particle velocity
          // MaxwellianVelocityFromFluidCell(x, y, z, &u, &v, &w);
        }
      }
}

void Particles::add_particles_domain(const FluidPicInterface& fluidInterface) {
  BL_PROFILE("Particles::add_particles");

  const int lev = 0;
  const auto dx = Geom(lev).CellSizeArray();
  const auto plo = Geom(lev).ProbLoArray();

  int iBlock = 0;
  for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
    const Box& tile_box = mfi.validbox();
    const auto lo = amrex::lbound(tile_box);
    const auto hi = amrex::ubound(tile_box);

    for (int i = lo.x; i <= hi.x; ++i)
      for (int j = lo.y; j <= hi.y; ++j)
        for (int k = lo.z; k <= hi.z; ++k) {
          add_particles_cell(mfi, fluidInterface, iBlock, i, j, k, lo.x, lo.y,
                             lo.z);
        }

    Print() << " lo.x = " << lo.x << " lo.y = " << lo.y << " lo.z = " << lo.z
            << " hi.x = " << hi.x << " hi.y = " << hi.y << " hi.z = " << hi.z
            << std::endl;

    iBlock++;
  }

  Print() << " speciesID = " << speciesID << "dx[0] = " << dx[0]
          << "dx[1] = " << dx[1] << "dx[2] = " << dx[2]
          << " plo[0] = " << plo[0] << " plo[1] = " << plo[1]
          << " plo[2] = " << plo[2] << std::endl;
}

void Particles::sum_moments(const MultiFab& Ex, const MultiFab& Ey,
                            const MultiFab& Ez, const MultiFab& Bx,
                            const MultiFab& By, const MultiFab& Bz,
                            MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt) {
  BL_PROFILE("Particles::sum_moments");
}

void Particles::mover(const MultiFab& Ex, const MultiFab& Ey,
                      const MultiFab& Ez, const MultiFab& Bx,
                      const MultiFab& By, const MultiFab& Bz, Real dt) {
  BL_PROFILE("Particles::mover");
}
