#include <cstdlib>

#include <AMReX_ParReduce.H>

#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::limit_weight(Real maxRatio,
                                                      bool seperateVelocity) {
  timing_func("Pts::limit_weight");

  if (maxRatio <= 1)
    return;

  IntVect iv(1);
  if (!(do_tiling && tile_size == iv))
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {

      Vector<ParticleType> newparticles;

      auto& pTile = get_particle_tile(iLev, pti);
      AoS& particles = pti.GetArrayOfStructs();

      // Sort the particles first to make sure the results
      // are the same for different number of processors
      // TODO C++20 support with better syntax:
      // std::ranges::sort(particles, compare_two_parts);
      std::sort(particles.begin(), particles.end(), compare_two_parts);

      Real totalMass = 0;
      Real totalMoment[nDim3] = { 0, 0, 0 };
      for (auto& p : particles) {
        totalMass += fabs(p.rdata(iqp_));
        for (int i = 0; i < nDim3; ++i)
          totalMoment[i] += fabs(p.rdata(iqp_) * p.rdata(iup_ + i));
      }
      Real avg = totalMass / particles.size();

      // Real maxWeight = avg + maxRatio * vars;
      Real maxWeight = avg * maxRatio;

      if (seperateVelocity) {
        Box bx = pti.tilebox();
        set_random_seed(iLev, bx.smallEnd(), IntVect(444));
        Vector<ParticleType*> pold;
        for (size_t ip = 0; ip < particles.size(); ip++) {
          Real qp1 = particles[ip].rdata(iqp_);
          if (fabs(qp1) < maxWeight)
            continue;
          pold.push_back(&(particles[ip]));
        }
        split_particles_by_velocity(pold, newparticles);
      } else {

        const auto lo = lbound(pti.tilebox());
        const auto hi = ubound(pti.tilebox());

        const Real xMin = Geom(iLev).LoEdge(lo.x, ix_) +
                          Geom(iLev).CellSize()[ix_] * 1e-10,
                   xMax = Geom(iLev).HiEdge(hi.x, ix_) -
                          Geom(iLev).CellSize()[ix_] * 1e-10;

        const Real yMin = Geom(iLev).LoEdge(lo.y, iy_) +
                          Geom(iLev).CellSize()[iy_] * 1e-10,
                   yMax = Geom(iLev).HiEdge(hi.y, iy_) -
                          Geom(iLev).CellSize()[iy_] * 1e-10;

        const Real zMin = nDim > 2 ? Geom(iLev).LoEdge(lo.z, iz_) +
                                         Geom(iLev).CellSize()[iz_] * 1e-10
                                   : 0.0,
                   zMax = nDim > 2 ? Geom(iLev).HiEdge(hi.z, iz_) -
                                         Geom(iLev).CellSize()[iz_] * 1e-10
                                   : 0.0;

        if (is_neutral()) {
          Box bx = pti.tilebox();
          set_random_seed(iLev, bx.smallEnd(), IntVect(999));
        }

        for (auto& p : particles) {
          Real qp1 = p.rdata(iqp_);
          if (fabs(qp1) < maxWeight)
            continue;

          Real up1 = p.rdata(iup_);
          Real vp1 = p.rdata(ivp_);
          Real wp1 = p.rdata(iwp_);

          Real xp1 = p.pos(ix_);
          Real yp1 = p.pos(iy_);
          Real zp1 = p.pos(iz_);

          int nNew = is_neutral() ? 7 : 1;

          p.rdata(iqp_) = qp1 / (nNew + 1);

          for (int iNew = 0; iNew < nNew; iNew++) {
            auto pnew = make_particle();
            set_ids(pnew);

            Real xp2 = xp1 + (xMax - xMin) * (randNum() - 0.5);
            Real yp2 = yp1 + (yMax - yMin) * (randNum() - 0.5);
            Real zp2 = zp1 + (zMax - zMin) * (randNum() - 0.5);

            xp2 = std::clamp(xp2, xMin, xMax);
            yp2 = std::clamp(yp2, yMin, yMax);
            zp2 = std::clamp(zp2, zMin, zMax);

            pnew.pos(ix_) = xp2;
            pnew.pos(iy_) = yp2;
            pnew.pos(iz_) = zp2;
            pnew.rdata(iup_) = up1;
            pnew.rdata(ivp_) = vp1;
            pnew.rdata(iwp_) = wp1;

            pnew.rdata(iqp_) = qp1 / (nNew + 1);
            newparticles.push_back(pnew);
          }
        }
      }

      for (auto& p : newparticles) {
        pTile.push_back(p);
      }
    }
  }
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::split_particles_by_velocity(
    Vector<ParticleType*>& plist, Vector<ParticleType>& newparticles) {

  if (plist.size() < 2)
    return;

  const int nCell = 8;
  // Assign the particle IDs to the corresponding velocity space cells.
  Vector<int> phasePartIdx_III[nCell][nCell][nCell];

  // Velocity domain range.
  Real velMin_D[nDim3] = { 1e9, 1e9, 1e9 },
       velMax_D[nDim3] = { -1e9, -1e9, -1e9 };

  for (int pid = 0; pid < plist.size(); pid++) {
    auto& pcl = *plist[pid];
    for (int iDir = 0; iDir < nDim3; iDir++) {
      if (pcl.rdata(iup_ + iDir) > velMax_D[iDir])
        velMax_D[iDir] = pcl.rdata(iup_ + iDir);
      if (pcl.rdata(iup_ + iDir) < velMin_D[iDir])
        velMin_D[iDir] = pcl.rdata(iup_ + iDir);
    }
  }

  Real dvMax = 0;
  for (int iDir = 0; iDir < nDim3; iDir++) {
    Real dv = velMax_D[iDir] - velMin_D[iDir];
    const Real vref = 0.5 * (fabs(velMax_D[iDir]) + fabs(velMin_D[iDir]));
    if (dv < 1e-6 * vref)
      dv = 1e-6 * vref;
    velMax_D[iDir] += 1e-3 * dv;
    velMin_D[iDir] -= 1e-3 * dv;
    dv = velMax_D[iDir] - velMin_D[iDir];

    if (dv > dvMax)
      dvMax = dv;
  }

  Real dvCell = dvMax == 0 ? 1e-9 : dvMax / nCell;
  Real invDv = 1 / dvCell;

  int iCell_D[nDim3];
  for (int pid = 0; pid < plist.size(); pid++) {
    auto& pcl = *plist[pid];
    for (int iDim = 0; iDim < nDim3; iDim++) {
      iCell_D[iDim] = fastfloor((pcl.rdata(iDim) - velMin_D[iDim]) * invDv);
    }

    phasePartIdx_III[iCell_D[ix_]][iCell_D[iy_]][iCell_D[iz_]].push_back(pid);
  }

  Vector<std::array<int, 3> > morton_idx(pow(nCell, 3));

  for (int iu = 0; iu < nCell; iu++)
    for (int iv = 0; iv < nCell; iv++)
      for (int iw = 0; iw < nCell; iw++) {
        morton_idx[encode_morton_3d(iu, iv, iw)] = { iu, iv, iw };
      }

  Vector<ParticleType*> p_morton;

  for (int i = 0; i < morton_idx.size(); ++i) {
    int iu = morton_idx[i][0];
    int iv = morton_idx[i][1];
    int iw = morton_idx[i][2];

    // printf("1 iu = %d iv = %d iw = %d\n", iu, iv, iw);
    for (int ip = 0; ip < phasePartIdx_III[iu][iv][iw].size(); ip++) {
      p_morton.push_back(plist[phasePartIdx_III[iu][iv][iw][ip]]);
    }
  }

  int nPair = floor(p_morton.size() / 2.0);
  for (int ip = 0; ip < nPair * 2; ip += 2) {
    ParticleType& p1 = *p_morton[ip];
    ParticleType& p2 = *p_morton[ip + 1];
    auto p3 = make_particle();
    auto p4 = make_particle();

    Real du = p1.rdata(iup_) - p2.rdata(iup_);
    Real dv = p1.rdata(ivp_) - p2.rdata(ivp_);
    Real dw = p1.rdata(iwp_) - p2.rdata(iwp_);
    Real dspeed = sqrt(du * du + dv * dv + dw * dw);

    if (dspeed / dvCell > 2)
      continue;

    bool doSucceed = split_by_seperate_velocity(p1, p2, p3, p4);
    if (doSucceed) {
      newparticles.push_back(p3);
      newparticles.push_back(p4);
    }
  }
}

template <int NStructReal, int NStructInt>
bool Particles<NStructReal, NStructInt>::split_by_seperate_velocity(
    ParticleType& p1, ParticleType& p2, ParticleType& p3, ParticleType& p4) {
  // AllPrint() << "Old: p1 = " << p1 << std::endl;
  // AllPrint() << "Old: p2 = " << p2 << std::endl;

  Real mt = p1.rdata(iqp_) + p2.rdata(iqp_);
  Real wavg = mt / 4.0;

  // Calculate the average velocity and total energy.
  Real et = 0, uavg2 = 0;
  Real u[nDim3], du1[nDim3], du2[nDim3];
  for (int i = 0; i < nDim3; ++i) {
    u[i] = (p1.rdata(iqp_) * p1.rdata(iup_ + i) +
            p2.rdata(iqp_) * p2.rdata(iup_ + i)) /
           mt;

    const Real u1 = p1.rdata(iup_ + i);
    const Real u2 = p2.rdata(iup_ + i);
    et += 0.5 * p1.rdata(iqp_) * u1 * u1;
    et += 0.5 * p2.rdata(iqp_) * u2 * u2;

    uavg2 += u[i] * u[i];

    // Get the direction of du1.
    du1[i] = p1.rdata(iup_ + i) - p2.rdata(iup_ + i);
  }

  Real du1Amp = l2_norm(du1, nDim3);
  if (du1Amp < 1e-16) {
    // If p1 and p2 have essentially the same velocity, do not split them. But
    // why the velocity difference can be so small?
    // A: with billions of particles, it can happen. I have done some
    // investigation, and it does not look like a bug.
    return false;
  }

  // The amplitude of du1 and du2.
  Real duAmp2 = et / (2 * wavg) - uavg2;
  if (duAmp2 < 0) {
    // Q: Why duAmp2 can be negative?
    // A: Rounding error.
    return false;
  }
  Real duAmp = sqrt(duAmp2);

  // Scale the amplitude of du1
  Real scale = duAmp / du1Amp;
  for (int i = 0; i < nDim3; ++i) {
    du1[i] *= scale;
  }

  {
    // Get the direction of du2
    Real utmp[nDim3];
    const Real r1 = randNum();
    const Real r2 = randNum();
    random_vector(r1, r2, utmp);

    // Correct the amplitide of du2
    for (int i = 0; i < nDim3; ++i) {
      du2[i] = utmp[i] * duAmp;
    }
  }

  // auto p_energy = [](const ParticleType& p) {
  //   Real energy = 0;
  //   for (int i = 0; i < nDim3; ++i) {
  //     energy += 0.5 * p.rdata(iqp_) * pow(p.rdata(iup_ + i), 2);
  //   }
  //   return energy;
  // };

  // Real eold = p_energy(p1) + p_energy(p2);
  // Real mold[3];
  // for (int i = 0; i < nDim3; ++i) {
  //   mold[i] = p1.rdata(iqp_) * p1.rdata(iup_ + i) +
  //             p2.rdata(iqp_) * p2.rdata(iup_ + i);
  // }

  set_ids(p3);
  set_ids(p4);

  p1.rdata(iqp_) = wavg;
  p2.rdata(iqp_) = wavg;
  p3.rdata(iqp_) = wavg;
  p4.rdata(iqp_) = wavg;

  for (int i = 0; i < nDim3; ++i) {
    p1.rdata(iup_ + i) = u[i] + du1[i];
    p2.rdata(iup_ + i) = u[i] - du1[i];

    p3.rdata(iup_ + i) = u[i] + du2[i];
    p4.rdata(iup_ + i) = u[i] - du2[i];
  }

  for (int i = 0; i < nDim; ++i) {
    p3.pos(i) = p1.pos(i);
    p4.pos(i) = p2.pos(i);
  }

  // Real enew = p_energy(p1) + p_energy(p2) + p_energy(p3) + p_energy(p4);
  // Real mnew[3];
  // for (int i = 0; i < nDim3; ++i) {
  //   mnew[i] = p1.rdata(iqp_) * p1.rdata(iup_ + i) +
  //             p2.rdata(iqp_) * p2.rdata(iup_ + i) +
  //             p3.rdata(iqp_) * p3.rdata(iup_ + i) +
  //             p4.rdata(iqp_) * p4.rdata(iup_ + i);
  // }

  // AllPrint() << "eold = " << eold << " enew = " << enew
  //            << " eold - enew = " << eold - enew << " mold - mnew "
  //            << mold[0] - mnew[0] << " " << mold[1] - mnew[1] << " "
  //            << mold[2] - mnew[2] << std::endl;

  // AllPrint() << "New: p1 = " << p1 << std::endl;
  // AllPrint() << "New: p2 = " << p2 << std::endl;
  // AllPrint() << "New: p3 = " << p3 << std::endl;
  // AllPrint() << "New: p4 = " << p4 << std::endl;

  return true;
}
//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::split_new(Real limit,
                                                   bool seperateVelocity) {
  timing_func("Pts::split");

  const int nInitial = product(nPartPerCell);

  IntVect iv = { AMREX_D_DECL(1, 1, 1) };
  if (!(do_tiling && tile_size == iv))
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {

    const Real vol = dx[iLev].product();
    const Real vacuumMass = vacuum * vol;

    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      Real dl = 0.1 * Geom(iLev).CellSize()[ix_] / (nPartPerCell.max());
      int nLowerLimit = nInitial * limit;
      int nGoal = nInitial;

      if (doPreSplitting) {
        const Array4<int const>& status = cell_status(iLev)[pti].array();
        Box bx = pti.tilebox();
        IntVect ibx = bx.smallEnd();
        if (bit::is_refined_neighbour(status(ibx))) {
          nLowerLimit = nLowerLimit * (pow(get_ref_ratio(iLev).max(), nDim));
          nGoal = nGoal * (pow(get_ref_ratio(iLev).max(), nDim));
          dl = dl / (get_ref_ratio(iLev).max());
        }
      }

      Vector<ParticleType> newparticles;

      auto& pTile = get_particle_tile(iLev, pti);
      AoS& particles = pTile.GetArrayOfStructs();

      const int nPartOrig = particles.size();
      if (nPartOrig > nLowerLimit)
        continue;

      const int nSplit =
          nGoal - nPartOrig > nPartOrig ? nPartOrig : nGoal - nPartOrig;
      Real totalMass = 0;
      for (auto& p : particles) {
        // So far, the vacuum limit is designed for OH-PT neutrals only. It is
        // not clear how should it be done for PIC, where electrons and ions
        // have different mass. --Yuxi
        totalMass += qomSign * p.rdata(iqp_);
      }
      if (totalMass < vacuumMass)
        continue;

      // Find the 'heaviest' nNew particles by sorting the weight
      // (charge).-----

      // Sort the particles by the location first to make sure the results
      // are the same for different number of processors
      std::sort(particles.begin(), particles.end(), compare_two_parts);

      // Sort the particles by the weight in decending order.
      std::sort(particles.begin(), particles.end(),
                [](const ParticleType& pl, const ParticleType& pr) {
                  const Real ql = fabs(pl.rdata(iqp_));
                  const Real qr = fabs(pr.rdata(iqp_));
                  if (fabs(ql - qr) > 1e-9 * (ql + qr)) {
                    return ql > qr;
                  }

                  if (fabs(pl.pos(ix_) - pr.pos(ix_)) >
                      1e-9 * (fabs(pl.pos(ix_)) + fabs(pr.pos(ix_)))) {
                    return pl.pos(ix_) > pr.pos(ix_);
                  }
                  return false;
                });
      //----------------------------------------------------------------

      const auto lo = lbound(pti.tilebox());
      const auto hi = ubound(pti.tilebox());

      const Real xMin = Geom(iLev).LoEdge(lo.x, ix_) +
                        Geom(iLev).CellSize()[ix_] * 1e-10,
                 xMax = Geom(iLev).HiEdge(hi.x, ix_) -
                        Geom(iLev).CellSize()[ix_] * 1e-10;

      const Real yMin = Geom(iLev).LoEdge(lo.y, iy_) +
                        Geom(iLev).CellSize()[iy_] * 1e-10,
                 yMax = Geom(iLev).HiEdge(hi.y, iy_) -
                        Geom(iLev).CellSize()[iy_] * 1e-10;

      const Real zMin = nDim > 2 ? Geom(iLev).LoEdge(lo.z, iz_) +
                                       Geom(iLev).CellSize()[iz_] * 1e-10
                                 : 0.0,
                 zMax = nDim > 2 ? Geom(iLev).HiEdge(hi.z, iz_) -
                                       Geom(iLev).CellSize()[iz_] * 1e-10
                                 : 0.0;

      if (is_neutral() || seperateVelocity) {
        Box bx = pti.tilebox();
        set_random_seed(iLev, bx.smallEnd(), IntVect(888));
      }

      if (seperateVelocity) {
        Vector<ParticleType*> pold;
        for (int ip = 0; ip < nSplit; ip++) {
          pold.push_back(&(particles[ip]));
        }
        split_particles_by_velocity(pold, newparticles);
      } else {
        for (int ip = 0; ip < nSplit; ip++) {
          auto& p = particles[ip];
          Real qp1 = p.rdata(iqp_);
          Real xp1 = p.pos(ix_);
          Real yp1 = p.pos(iy_);
          Real zp1 = nDim > 2 ? p.pos(iz_) : 0;
          Real up1 = p.rdata(iup_);
          Real vp1 = p.rdata(ivp_);
          Real wp1 = p.rdata(iwp_);

          const Real u2 = up1 * up1 + vp1 * vp1 + wp1 * wp1;

          Real coef = (u2 < 1e-13) ? 0 : dl / sqrt(u2);
          const Real dpx = coef * up1;
          const Real dpy = coef * vp1;
          const Real dpz = coef * wp1;

          Real xp2 = xp1 + dpx;
          Real yp2 = yp1 + dpy;
          Real zp2 = zp1 + dpz;

          int nNew = is_neutral() ? 7 : 1;

          p.rdata(iqp_) = qp1 / (nNew + 1.0);

          for (int iNew = 0; iNew < nNew; iNew++) {

            if (is_neutral()) {
              xp2 = xp1 + (xMax - xMin) * (randNum() - 0.5);
              yp2 = yp1 + (yMax - yMin) * (randNum() - 0.5);
              zp2 = zp1 + (zMax - zMin) * (randNum() - 0.5);
            } else {
              xp1 -= dpx;
              yp1 -= dpy;
              zp1 -= dpz;

              xp1 = std::clamp(xp1, xMin, xMax);
              yp1 = std::clamp(yp1, yMin, yMax);
              zp1 = std::clamp(zp1, zMin, zMax);
              p.pos(ix_) = xp1;
              p.pos(iy_) = yp1;

              if (nDim > 2)
                p.pos(iz_) = zp1;
            }

            xp2 = std::clamp(xp2, xMin, xMax);
            yp2 = std::clamp(yp2, yMin, yMax);
            zp2 = std::clamp(zp2, zMin, zMax);

            auto pnew = make_particle();
            set_ids(pnew);

            pnew.pos(ix_) = xp2;
            pnew.pos(iy_) = yp2;
            if (nDim > 2)
              pnew.pos(iz_) = zp2;
            pnew.rdata(iup_) = up1;
            pnew.rdata(ivp_) = vp1;
            pnew.rdata(iwp_) = wp1;
            pnew.rdata(iqp_) = qp1 / (nNew + 1.0);
            newparticles.push_back(pnew);
          }
        }
      }

      for (auto& p : newparticles) {
        pTile.push_back(p);
      }
    }
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::split(Real limit,
                                               bool seperateVelocity) {
  timing_func("Pts::split");

  const int nInitial = product(nPartPerCell);

  IntVect iv = { AMREX_D_DECL(1, 1, 1) };
  if (!(do_tiling && tile_size == iv))
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {

    const Real dl = 0.1 * Geom(iLev).CellSize()[ix_] / nPartPerCell.max();

    const int nLowerLimit = nInitial * limit * pow(pLevRatio, iLev);

    const int nGoal = nLowerLimit > nInitial ? nLowerLimit : nInitial;

    const Real vol = dx[iLev].product();

    const Real vacuumMass = vacuum * vol;

    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {

      Vector<ParticleType> newparticles;

      auto& pTile = get_particle_tile(iLev, pti);
      AoS& particles = pTile.GetArrayOfStructs();

      const int nPartOrig = particles.size();

      if (nPartOrig > nLowerLimit)
        continue;

      const int nSplit =
          nGoal - nPartOrig > nPartOrig ? nPartOrig : nGoal - nPartOrig;

      Real totalMass = 0;
      for (auto& p : particles) {
        // So far, the vacuum limit is designed for OH-PT neutrals only. It is
        // not clear how should it be done for PIC, where electrons and ions
        // have different mass. --Yuxi
        totalMass += qomSign * p.rdata(iqp_);
      }
      if (totalMass < vacuumMass)
        continue;

      // Find the 'heaviest' nNew particles by sorting the weight
      // (charge).-----

      // Sort the particles by the location first to make sure the results
      // are the same for different number of processors
      std::sort(particles.begin(), particles.end(), compare_two_parts);

      // Sort the particles by the weight in decending order.
      std::sort(particles.begin(), particles.end(),
                [](const ParticleType& pl, const ParticleType& pr) {
                  const Real ql = fabs(pl.rdata(iqp_));
                  const Real qr = fabs(pr.rdata(iqp_));
                  if (fabs(ql - qr) > 1e-9 * (ql + qr)) {
                    return ql > qr;
                  }

                  if (fabs(pl.pos(ix_) - pr.pos(ix_)) >
                      1e-9 * (fabs(pl.pos(ix_)) + fabs(pr.pos(ix_)))) {
                    return pl.pos(ix_) > pr.pos(ix_);
                  }
                  return false;
                });
      //----------------------------------------------------------------

      const auto lo = lbound(pti.tilebox());
      const auto hi = ubound(pti.tilebox());

      const Real xMin = Geom(iLev).LoEdge(lo.x, ix_) +
                        Geom(iLev).CellSize()[ix_] * 1e-10,
                 xMax = Geom(iLev).HiEdge(hi.x, ix_) -
                        Geom(iLev).CellSize()[ix_] * 1e-10;

      const Real yMin = Geom(iLev).LoEdge(lo.y, iy_) +
                        Geom(iLev).CellSize()[iy_] * 1e-10,
                 yMax = Geom(iLev).HiEdge(hi.y, iy_) -
                        Geom(iLev).CellSize()[iy_] * 1e-10;

      const Real zMin = nDim > 2 ? Geom(iLev).LoEdge(lo.z, iz_) +
                                       Geom(iLev).CellSize()[iz_] * 1e-10
                                 : 0.0,
                 zMax = nDim > 2 ? Geom(iLev).HiEdge(hi.z, iz_) -
                                       Geom(iLev).CellSize()[iz_] * 1e-10
                                 : 0.0;

      if (is_neutral() || seperateVelocity) {
        Box bx = pti.tilebox();
        set_random_seed(iLev, bx.smallEnd(), IntVect(888));
      }

      if (seperateVelocity) {
        Vector<ParticleType*> pold;
        for (int ip = 0; ip < nSplit; ip++) {
          pold.push_back(&(particles[ip]));
        }
        split_particles_by_velocity(pold, newparticles);
      } else {
        for (int ip = 0; ip < nSplit; ip++) {
          auto& p = particles[ip];
          Real qp1 = p.rdata(iqp_);
          Real xp1 = p.pos(ix_);
          Real yp1 = p.pos(iy_);
          Real zp1 = nDim > 2 ? p.pos(iz_) : 0;
          Real up1 = p.rdata(iup_);
          Real vp1 = p.rdata(ivp_);
          Real wp1 = p.rdata(iwp_);

          const Real u2 = up1 * up1 + vp1 * vp1 + wp1 * wp1;

          Real coef = (u2 < 1e-13) ? 0 : dl / sqrt(u2);
          const Real dpx = coef * up1;
          const Real dpy = coef * vp1;
          const Real dpz = coef * wp1;

          Real xp2 = xp1 + dpx;
          Real yp2 = yp1 + dpy;
          Real zp2 = zp1 + dpz;

          int nNew = is_neutral() ? 7 : 1;

          p.rdata(iqp_) = qp1 / (nNew + 1.0);

          for (int iNew = 0; iNew < nNew; iNew++) {

            if (is_neutral()) {
              xp2 = xp1 + (xMax - xMin) * (randNum() - 0.5);
              yp2 = yp1 + (yMax - yMin) * (randNum() - 0.5);
              zp2 = zp1 + (zMax - zMin) * (randNum() - 0.5);
            } else {
              xp1 -= dpx;
              yp1 -= dpy;
              zp1 -= dpz;

              xp1 = std::clamp(xp1, xMin, xMax);
              yp1 = std::clamp(yp1, yMin, yMax);
              zp1 = std::clamp(zp1, zMin, zMax);
              p.pos(ix_) = xp1;
              p.pos(iy_) = yp1;

              if (nDim > 2)
                p.pos(iz_) = zp1;
            }

            xp2 = std::clamp(xp2, xMin, xMax);
            yp2 = std::clamp(yp2, yMin, yMax);
            zp2 = std::clamp(zp2, zMin, zMax);

            auto pnew = make_particle();
            set_ids(pnew);

            pnew.pos(ix_) = xp2;
            pnew.pos(iy_) = yp2;
            if (nDim > 2)
              pnew.pos(iz_) = zp2;
            pnew.rdata(iup_) = up1;
            pnew.rdata(ivp_) = vp1;
            pnew.rdata(iwp_) = wp1;
            pnew.rdata(iqp_) = qp1 / (nNew + 1.0);
            newparticles.push_back(pnew);
          }
        }
      }

      for (auto& p : newparticles) {
        pTile.push_back(p);
      }
    }
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
bool Particles<NStructReal, NStructInt>::merge_particles_accurate(
    int iLev, AoS& particles, Vector<int>& partIdx, Vector<int>& idx_I,
    int nPartCombine, int nPartNew, Vector<Real>& x, Real velNorm) {
  timing_func("Pts::merge_particles_accurate");

  constexpr int nVar = 5;
  constexpr int iq_ = 0, iu_ = 1, iv_ = 2, iw_ = 3, ie_ = 4;
  const Real coefVel = 1, coefPos = 1;

  Vector<Real> ref(nVar, 0);
  Array2D<Real, 0, nVar - 1, 0, nVar> a;
  for (int i = 0; i < nVar; ++i)
    for (int j = 0; j < nVar + 1; ++j) {
      a(i, j) = 0;
    }

  // Find the center of the particles, and sort the particles based
  // on its distance to the 6-D center.
  //----------------------------------------------------------
  Vector<Real> middle(nDim + nDim3, 0);
  for (int pID : partIdx) {
    for (int iDir = ix_; iDir <= iz_; iDir++) {
      if (iDir < nDim)
        middle[iDir] += particles[pID].pos(iDir);
      middle[nDim + iDir] += particles[pID].rdata(iDir);
    }
  }

  for (int i = 0; i < middle.size(); ++i) {
    middle[i] /= partIdx.size();
  }

  auto calc_distance2_to_center = [&, this](int pID) {
    Real dl2 = 0, dvel2 = 0;
    for (int iDir = ix_; iDir <= iz_; iDir++) {

      if (iDir < nDim) {
        Real pos = particles[pID].pos(iDir);
        const Real distance = (pos - middle[iDir]) * invDx[iLev][iDir];
        dl2 += distance * distance;
      }

      Real v = particles[pID].rdata(iDir);
      const Real velocity = (v - middle[nDim + iDir]) * velNorm;
      dvel2 += velocity * velocity;
    }
    return coefPos * dl2 + coefVel * dvel2;
  };

  std::sort(partIdx.begin(), partIdx.end(),
            [calc_distance2_to_center](const int& idl, const int& idr) {
              Real dll = calc_distance2_to_center(idl);
              Real dlr = calc_distance2_to_center(idr);

              if (fabs(dll - dlr) > 1e-9 * (dll + dlr)) {
                return dll < dlr;
              }
              return false;
            });

  /*
      Delete 1 particle out of 6 particles:
      1) Choose two particles that are closest to each other.
      2) Delete the lighter one.
      3) Distribute its weights to another 5 particles to conserve
         mass, momentum and energy.
   */

  idx_I.resize(nPartCombine, 0);
  for (int ip = 0; ip < nPartCombine; ip++) {
    idx_I[ip] = partIdx[ip];
  }

  // Calculate the center of the particles for combination
  for (int i = 0; i < middle.size(); ++i) {
    middle[i] = 0;
  }
  for (int pID : idx_I) {
    for (int iDir = ix_; iDir <= iz_; iDir++) {
      if (iDir < nDim)
        middle[iDir] += particles[pID].pos(iDir);

      middle[nDim + iDir] += particles[pID].rdata(iDir);
    }
  }
  for (int i = 0; i < middle.size(); ++i) {
    middle[i] /= nPartCombine;
  }

  bool doCombine = true;
  for (int pID : idx_I) {
    Real distance = sqrt(calc_distance2_to_center(pID));
    if (distance > mergeThresholdDistance) {
      // printf("Warning: distance=%e\n", distance);
      doCombine = false;
    }
  }

  if (!doCombine)
    return false;

  // Find the pair that is closest to each other in phase space
  int pair1 = 0, pair2 = 0;
  Real dis2Min = 2;
  for (int ip1 = 0; ip1 < nPartCombine - 1; ip1++)
    for (int ip2 = ip1 + 1; ip2 < nPartCombine; ip2++) {

      // Distance between two particles in 6D space.
      Real dl2 = 0, dv2 = 0;
      for (int iDir = 0; iDir < nDim; iDir++) {
        Real dv = velNorm * (particles[idx_I[ip1]].rdata(iDir) -
                             particles[idx_I[ip2]].rdata(iDir));
        dv2 += dv * dv;

        Real dx = invDx[iLev][iDir] * (particles[idx_I[ip1]].pos(iDir) -
                                       particles[idx_I[ip2]].pos(iDir));
        dv2 += dx * dx;
      }

      const Real dis2 = dv2 * coefVel + dl2 * coefPos;

      if (dis2 < dis2Min) {
        dis2Min = dis2;
        pair1 = ip1;
        pair2 = ip2;
      }
    }
  //-------------------------------

  // Delete the lighter one.
  int iPartDel = pair1;
  // Q: Why is it 'l>(1+1e-9)*r' instead of 'l>r'?
  // A: The particle weights can be the same for some cases. 'l>r'
  // may return random results due to the truncation error.
  if (fabs(particles[idx_I[pair1]].rdata(iqp_)) >
      (1 + 1e-9) * fabs(particles[idx_I[pair2]].rdata(iqp_))) {
    iPartDel = pair2;
  }

  std::swap(idx_I[iPartDel], idx_I[nPartCombine - 1]);

  //-----------Solve the new particle weights-------
  for (int ip = 0; ip < nPartCombine; ip++) {
    const Real qp = particles[idx_I[ip]].rdata(iqp_);
    const Real up = particles[idx_I[ip]].rdata(iup_);
    const Real vp = particles[idx_I[ip]].rdata(ivp_);
    const Real wp = particles[idx_I[ip]].rdata(iwp_);
    const Real v2 = up * up + vp * vp + wp * wp;

    if (ip < nVar) {
      a(iq_, ip) = 1;
      a(iu_, ip) = up;
      a(iv_, ip) = vp;
      a(iw_, ip) = wp;
      a(ie_, ip) = v2;
    }

    a(iq_, nVar) += qp;
    a(iu_, nVar) += qp * up;
    a(iv_, nVar) += qp * vp;
    a(iw_, nVar) += qp * wp;
    a(ie_, nVar) += qp * v2;
  }

  const Real csmall = 1e-9;
  const Real tmp = csmall * fabs(1. / a(iq_, nVar));
  for (int i = iq_; i <= ie_; ++i) {
    // Ensure a strictly positive threshold so constraint rows with zero RHS and
    // zero coefficients (e.g., vz = 0 in 1D/2D) are flagged as singular,
    // preventing near-zero pivots or 0/0 NaNs in the weights.
    ref[i] = std::max(fabs(a(i, nVar) * tmp), csmall * tmp);
  }

  x.resize(nVar, 0);
  bool isSolved = linear_solver_Gauss_Elimination<Real, nVar, nVar + 1>(
      nVar, nVar + 1, a, x, ref);

  if (isSolved) {
    // All the particle weights should have the same sign.
    Real qt = x[0];
    for (int ip = 0; ip < nPartNew; ip++) {
      if (qt * x[ip] <= 0) {
        isSolved = false;
        break;
      }
    }
  }

  return isSolved;
}

//==========================================================

template <int NStructReal, int NStructInt>
bool Particles<NStructReal, NStructInt>::merge_particles_fast(
    int iLev, AoS& particles, Vector<int>& partIdx, Vector<int>& idx_I,
    int nPartCombine, int nPartNew, Vector<Real>& x, long seed) {
  timing_func("Pts::merge_particles_fast");

  constexpr int iq_ = 0, iu_ = 1, iv_ = 2, iw_ = 3, ie_ = 4;
  constexpr int nPartNewMax = 16;
  constexpr int nVarMax = nPartNewMax + 5;

  int nVar = nPartNew + 5;

  Vector<Real> ref(nVar, 0);
  Array2D<Real, 0, nVarMax - 1, 0, nVarMax> a;
  for (int i = 0; i < nVar; ++i)
    for (int j = 0; j < nVar + 1; ++j) {
      a(i, j) = 0;
    }

  // Q: Sort the particles by weights in ascending order.
  // But, why is it required here?
  // A: Eliminate randomness.
  std::sort(partIdx.begin(), partIdx.end(),
            [&particles](int idLeft, int idRight) {
              const Real ql = fabs(particles[idLeft].rdata(iqp_));
              const Real qr = fabs(particles[idRight].rdata(iqp_));
              if (fabs(ql - qr) > 1e-9 * (ql + qr)) {
                return ql < qr;
              }

              Real xl = particles[idLeft].pos(ix_);
              Real xr = particles[idRight].pos(ix_);
              if (fabs(xl - xr) > 1e-9 * (fabs(xl) + fabs(xr))) {
                return xl < xr;
              }
              return false;
            });

  if (mergeLight) {
    idx_I.resize(nPartCombine, 0);
    for (int ip = 0; ip < nPartCombine; ip++) {
      idx_I[ip] = partIdx[ip];
    }

    Real plight = 1e99, pheavy = 0;
    for (int ip = 0; ip < nPartCombine; ip++) {
      auto& p = particles[idx_I[ip]];
      Real w = fabs(p.rdata(iqp_));
      if (w < plight)
        plight = w;
      if (w > pheavy)
        pheavy = w;
    }

    if (pheavy / plight > mergePartRatioMax)
      return false;

    randNum.set_seed(seed);
    shuffle_fish_yates(idx_I, randNum);

  } else {
    randNum.set_seed(seed);
    shuffle_fish_yates(partIdx, randNum);

    idx_I.resize(nPartCombine, 0);
    for (int ip = 0; ip < nPartCombine; ip++) {
      idx_I[ip] = partIdx[ip];
    }
  }

  // Sum the moments of all the old particles.
  for (int ip = 0; ip < nPartCombine; ip++) {
    const Real qp = particles[idx_I[ip]].rdata(iqp_);
    const Real up = particles[idx_I[ip]].rdata(iup_);
    const Real vp = particles[idx_I[ip]].rdata(ivp_);
    const Real wp = particles[idx_I[ip]].rdata(iwp_);
    const Real v2 = 0.5 * (up * up + vp * vp + wp * wp);
    a(nPartNew + iq_, nVar) += qp;
    a(nPartNew + iu_, nVar) += qp * up;
    a(nPartNew + iv_, nVar) += qp * vp;
    a(nPartNew + iw_, nVar) += qp * wp;
    a(nPartNew + ie_, nVar) += qp * v2;
  }

  const Real invAvg = 2 * nPartNew / a(nPartNew + iq_, nVar);
  for (int ip = 0; ip < nPartNew; ip++) {
    const Real qp = particles[idx_I[ip]].rdata(iqp_);
    const Real up = particles[idx_I[ip]].rdata(iup_);
    const Real vp = particles[idx_I[ip]].rdata(ivp_);
    const Real wp = particles[idx_I[ip]].rdata(iwp_);
    const Real v2 = 0.5 * (up * up + vp * vp + wp * wp);

    a(ip, nVar) = 2;

    a(ip, ip) = 2. / qp;
    a(ip, nPartNew + iq_) = 1;
    a(ip, nPartNew + iu_) = up;
    a(ip, nPartNew + iv_) = vp;
    a(ip, nPartNew + iw_) = wp;
    a(ip, nPartNew + ie_) = v2;

    a(nPartNew + iq_, ip) = 1;
    a(nPartNew + iu_, ip) = up;
    a(nPartNew + iv_, ip) = vp;
    a(nPartNew + iw_, ip) = wp;
    a(nPartNew + ie_, ip) = v2;
  }

  const Real csmall = 1e-9;
  const Real tmp = csmall * invAvg;
  for (int i = 0; i < nVar; ++i) {
    if (i < nPartNew) {
      ref[i] = tmp;
    } else {
      // Floor at strictly positive: zero-RHS constraint rows (e.g. vz = 0 in
      // 1D/2D) must be flagged as singular rather than accepting near-zero
      // pivots.
      ref[i] = std::max(fabs(a(i, nVar) * tmp * csmall), csmall * tmp);
    }
  }

  x.resize(nVar, 0);
  bool isSolved = linear_solver_Gauss_Elimination<Real, nVarMax, nVarMax + 1>(
      nVar, nVar + 1, a, x, ref);

  if (isSolved) {
    // All the particle weights should have the same sign.
    Real qt = x[0];
    for (int ip = 0; ip < nPartNew; ip++) {
      if (qt * x[ip] <= 0) {
        isSolved = false;
        break;
      }
    }
  }

  if (isSolved) {
    for (int ip = 0; ip < nPartNew; ip++) {
      Real pold = particles[idx_I[ip]].rdata(iqp_);
      Real pnew = x[ip];

      Real c0 = Real(nPartCombine) / nPartNew * mergeRatioMax;
      if (pnew / pold > c0 || pold / pnew > c0) {
        isSolved = false;
        break;
      }
    }
  }

  return isSolved;
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::merge(Real limit) {
  timing_func("Pts::merge");
  IntVect iv = { AMREX_D_DECL(1, 1, 1) };
  if (!(do_tiling && tile_size == iv))
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {

    int nPartGoal = product(nPartPerCell) * limit * pow(pLevRatio, iLev);

    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {

      // It is assumed the tile size is 1x1x1.
      Box bx = pti.tilebox();
      long seed = set_random_seed(iLev, bx.smallEnd(), IntVect(777));

      AoS& particles = pti.GetArrayOfStructs();

      const int nPartOrig = particles.size();

      if (nPartOrig <= nPartGoal)
        continue;

      // The range of the velocity domain:
      // [-r0,r0]*thermal_velocity+bulk_velocity
      const Real r0 = fastMerge ? 2.0 : 1.0;

      // Phase space cell number in one direction.
      // The const 0.5/0.8 is choosen by experiments.
      int nCell = 0;
      if (fastMerge) {
        nCell = r0 * ceil(0.5 * pow(nPartOrig, 1. / nDim3));
      } else {
        nCell = r0 * ceil(0.8 * pow(nPartOrig, 1. / nDim3));
      }

      if (nCell < 3)
        continue;

      // Sort the particles by the location first to make sure the results
      // are the same for different number of processors
      std::sort(particles.begin(), particles.end(), compare_two_parts);

      // One particle may belong to more than one velocity bins, but it can be
      // only merged at most once.
      std::vector<bool> merged;
      merged.resize(nPartOrig, false);

      //----------------------------------------------------------------
      // Estimate the bulk velocity and thermal velocity.
      Real uBulk[nDim3] = { 0, 0, 0 };
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];
        for (int iDir = 0; iDir < 3; iDir++) {
          uBulk[iDir] += pcl.rdata(iDir);
        }
      }

      for (int iDir = 0; iDir < nDim3; iDir++) {
        uBulk[iDir] /= nPartOrig;
      }

      Real thVel = 0, thVel2 = 0;
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];
        for (int iDir = 0; iDir < nDim3; iDir++) {
          thVel2 += pow(pcl.rdata(iDir) - uBulk[iDir], 2);
        }
      }

      thVel2 /= nPartOrig;
      thVel = sqrt(thVel2);

      // The coef 0.5 if choosen by experience.
      const Real velNorm = (thVel < 1e-13) ? 0 : 1.0 / (0.5 * thVel);
      //----------------------------------------------------------------

      //----------------------------------------------------------------
      // Assign the particle IDs to the corresponding velocity space cells.
      Vector<Vector<int> > phasePartIdx_I(nCell * nCell * nCell);
      const auto bin_index = [nCell](int i, int j, int k) {
        return (i * nCell + j) * nCell + k;
      };

      Real dv = (2.0 * r0 * thVel) / nCell;
      Real invDv = (dv < 1e-13) ? 0 : 1.0 / dv;

      // Velocity domain range.
      Real velMin_D[nDim3], velMax_D[nDim3];
      for (int iDir = 0; iDir < nDim3; iDir++) {
        Real dvshift = (randNum() - 0.5) * dv;
        velMin_D[iDir] = -r0 * thVel + uBulk[iDir] + dvshift;
        velMax_D[iDir] = r0 * thVel + uBulk[iDir] + dvshift;
      }

      int iCell_D[nDim3];
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];

        bool isOutside = false;
        for (int iDim = 0; iDim < nDim3; iDim++) {
          if (pcl.rdata(iDim) < velMin_D[iDim] ||
              pcl.rdata(iDim) > velMax_D[iDim])
            isOutside = true;
        }
        if (isOutside)
          continue;

        for (int iDim = 0; iDim < nDim3; iDim++) {
          iCell_D[iDim] = fastfloor((pcl.rdata(iDim) - velMin_D[iDim]) * invDv);
        }

        // One particle may belong to multiple bins when each bin has a buffer
        // region.
        for (int xCell = iCell_D[ix_] - 1; xCell <= iCell_D[ix_] + 1; xCell++)
          for (int yCell = iCell_D[iy_] - 1; yCell <= iCell_D[iy_] + 1; yCell++)
            for (int zCell = iCell_D[iz_] - 1; zCell <= iCell_D[iz_] + 1;
                 zCell++) {

              if (xCell < 0 || xCell >= nCell || yCell < 0 || yCell >= nCell ||
                  zCell < 0 || zCell >= nCell)
                continue;

              Vector<int> cellIdx = { xCell, yCell, zCell };

              Real binMin_D[nDim3], binMax_D[nDim3];

              for (int iDim = 0; iDim < nDim3; iDim++) {
                binMin_D[iDim] =
                    velMin_D[iDim] + (cellIdx[iDim] - velBinBufferSize) * dv;

                binMax_D[iDim] = velMin_D[iDim] +
                                 (cellIdx[iDim] + 1 + velBinBufferSize) * dv;
              }

              bool isInside = true;
              for (int iDim = 0; iDim < nDim3; iDim++) {
                if (pcl.rdata(iDim) < binMin_D[iDim] ||
                    pcl.rdata(iDim) > binMax_D[iDim])
                  isInside = false;
              }

              if (isInside) {
                phasePartIdx_I[bin_index(cellIdx[ix_], cellIdx[iy_],
                                         cellIdx[iz_])]
                    .push_back(pid);
              }
            }
      }
      //----------------------------------------------------------------

      for (int iu = 0; iu < nCell; iu++)
        for (int iv = 0; iv < nCell; iv++)
          for (int iw = 0; iw < nCell; iw++) {
            Vector<int> partIdx;
            auto& pIdx = phasePartIdx_I[bin_index(iu, iv, iw)];
            for (int i = 0; i < pIdx.size(); ++i) {
              int pid = pIdx[i];
              if (!merged[pid]) {
                partIdx.push_back(pid);
              }
            }

            if (partIdx.size() < nPartNew + 1)
              continue;

            Vector<Real> x;
            Vector<int> idx_I;

            int nOld = nPartCombine;
            bool isSolved;
            if (fastMerge) {
              if (nOld > partIdx.size())
                nOld = partIdx.size();

              for (int iTry = 0; iTry < nMergeTry; iTry++) {
                long sd = seed + iu * 777 + iv * 77 + iw + iTry;
                isSolved = merge_particles_fast(iLev, particles, partIdx, idx_I,
                                                nOld, nPartNew, x, sd);
                if (isSolved)
                  break;
              }
            } else {
              isSolved = merge_particles_accurate(
                  iLev, particles, partIdx, idx_I, nOld, nPartNew, x, velNorm);
            }
            if (!isSolved)
              continue;

            //----------------------------------------------

            // Reject merge if any solved weight is non-finite (NaN / Inf).
            bool xIsFinite = true;
            for (int ip = 0; ip < nPartNew; ++ip) {
              if (!std::isfinite(x[ip])) {
                xIsFinite = false;
                break;
              }
            }
            if (!xIsFinite)
              continue;

            // Adjust weight.
            for (int ip = 0; ip < nPartNew; ip++) {
              auto& p = particles[idx_I[ip]];
              p.rdata(iqp_) = x[ip];
              merged[idx_I[ip]] = true;
            }
            // Mark for deletion
            for (int ip = nPartNew; ip < nOld; ip++) {
              particles[idx_I[ip]].id() = -1;
              particles[idx_I[ip]].rdata(iqp_) = 0;
              merged[idx_I[ip]] = true;
            }
          }
    }
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::merge_new(Real limit) {
  timing_func("Pts::merge");
  IntVect iv = { AMREX_D_DECL(1, 1, 1) };
  if (!(do_tiling && tile_size == iv))
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      const auto tppc = target_PPC(iLev)[pti].array();
      const Box& bx = pti.tilebox();
      IntVect ibx = bx.smallEnd();
      int target = tppc(ibx);
      int nPartGoal = target * limit;

      // It is assumed the tile size is 1x1x1.
      long seed = set_random_seed(iLev, bx.smallEnd(), IntVect(777));

      AoS& particles = pti.GetArrayOfStructs();

      const int nPartOrig = particles.size();

      if (nPartOrig <= nPartGoal)
        continue;

      // The range of the velocity domain:
      // [-r0,r0]*thermal_velocity+bulk_velocity
      const Real r0 = fastMerge ? 2.0 : 1.0;

      // Phase space cell number in one direction.
      // The const 0.5/0.8 is choosen by experiments.
      int nCell = 0;
      if (fastMerge) {
        nCell = r0 * ceil(0.5 * pow(nPartOrig, 1. / nDim3));
      } else {
        nCell = r0 * ceil(0.8 * pow(nPartOrig, 1. / nDim3));
      }

      if (nCell < 3)
        continue;

      // Sort the particles by the location first to make sure the results
      // are the same for different number of processors
      std::sort(particles.begin(), particles.end(), compare_two_parts);

      // One particle may belong to more than one velocity bins, but it can be
      // only merged at most once.
      std::vector<bool> merged;
      merged.resize(nPartOrig, false);

      //----------------------------------------------------------------
      // Estimate the bulk velocity and thermal velocity.
      Real uBulk[nDim3] = { 0, 0, 0 };
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];
        for (int iDir = 0; iDir < 3; iDir++) {
          uBulk[iDir] += pcl.rdata(iDir);
        }
      }

      for (int iDir = 0; iDir < nDim3; iDir++) {
        uBulk[iDir] /= nPartOrig;
      }

      Real thVel = 0, thVel2 = 0;
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];
        for (int iDir = 0; iDir < nDim3; iDir++) {
          thVel2 += pow(pcl.rdata(iDir) - uBulk[iDir], 2);
        }
      }

      thVel2 /= nPartOrig;
      thVel = sqrt(thVel2);

      // The coef 0.5 if choosen by experience.
      const Real velNorm = (thVel < 1e-13) ? 0 : 1.0 / (0.5 * thVel);
      //----------------------------------------------------------------

      //----------------------------------------------------------------
      // Assign the particle IDs to the corresponding velocity space cells.
      Vector<Vector<int> > phasePartIdx_I(nCell * nCell * nCell);
      const auto bin_index = [nCell](int i, int j, int k) {
        return (i * nCell + j) * nCell + k;
      };

      Real dv = (2.0 * r0 * thVel) / nCell;
      Real invDv = (dv < 1e-13) ? 0 : 1.0 / dv;

      // Velocity domain range.
      Real velMin_D[nDim3], velMax_D[nDim3];
      for (int iDir = 0; iDir < nDim3; iDir++) {
        Real dvshift = (randNum() - 0.5) * dv;
        velMin_D[iDir] = -r0 * thVel + uBulk[iDir] + dvshift;
        velMax_D[iDir] = r0 * thVel + uBulk[iDir] + dvshift;
      }

      int iCell_D[nDim3];
      for (int pid = 0; pid < nPartOrig; pid++) {
        auto& pcl = particles[pid];

        bool isOutside = false;
        for (int iDim = 0; iDim < nDim3; iDim++) {
          if (pcl.rdata(iDim) < velMin_D[iDim] ||
              pcl.rdata(iDim) > velMax_D[iDim])
            isOutside = true;
        }
        if (isOutside)
          continue;

        for (int iDim = 0; iDim < nDim3; iDim++) {
          iCell_D[iDim] = fastfloor((pcl.rdata(iDim) - velMin_D[iDim]) * invDv);
        }

        // One particle may belong to multiple bins when each bin has a buffer
        // region.
        for (int xCell = iCell_D[ix_] - 1; xCell <= iCell_D[ix_] + 1; xCell++)
          for (int yCell = iCell_D[iy_] - 1; yCell <= iCell_D[iy_] + 1; yCell++)
            for (int zCell = iCell_D[iz_] - 1; zCell <= iCell_D[iz_] + 1;
                 zCell++) {

              if (xCell < 0 || xCell >= nCell || yCell < 0 || yCell >= nCell ||
                  zCell < 0 || zCell >= nCell)
                continue;

              Vector<int> cellIdx = { xCell, yCell, zCell };

              Real binMin_D[nDim3], binMax_D[nDim3];

              for (int iDim = 0; iDim < nDim3; iDim++) {
                binMin_D[iDim] =
                    velMin_D[iDim] + (cellIdx[iDim] - velBinBufferSize) * dv;

                binMax_D[iDim] = velMin_D[iDim] +
                                 (cellIdx[iDim] + 1 + velBinBufferSize) * dv;
              }

              bool isInside = true;
              for (int iDim = 0; iDim < nDim3; iDim++) {
                if (pcl.rdata(iDim) < binMin_D[iDim] ||
                    pcl.rdata(iDim) > binMax_D[iDim])
                  isInside = false;
              }

              if (isInside) {
                phasePartIdx_I[bin_index(cellIdx[ix_], cellIdx[iy_],
                                         cellIdx[iz_])]
                    .push_back(pid);
              }
            }
      }
      //----------------------------------------------------------------

      for (int iu = 0; iu < nCell; iu++)
        for (int iv = 0; iv < nCell; iv++)
          for (int iw = 0; iw < nCell; iw++) {
            Vector<int> partIdx;
            auto& pIdx = phasePartIdx_I[bin_index(iu, iv, iw)];
            for (int i = 0; i < pIdx.size(); ++i) {
              int pid = pIdx[i];
              if (!merged[pid]) {
                partIdx.push_back(pid);
              }
            }

            if (partIdx.size() < nPartNew + 1)
              continue;

            Vector<Real> x;
            Vector<int> idx_I;

            int nOld = nPartCombine;
            bool isSolved;
            if (fastMerge) {
              if (nOld > partIdx.size())
                nOld = partIdx.size();

              for (int iTry = 0; iTry < nMergeTry; iTry++) {
                long sd = seed + iu * 777 + iv * 77 + iw + iTry;
                isSolved = merge_particles_fast(iLev, particles, partIdx, idx_I,
                                                nOld, nPartNew, x, sd);
                if (isSolved)
                  break;
              }
            } else {
              isSolved = merge_particles_accurate(
                  iLev, particles, partIdx, idx_I, nOld, nPartNew, x, velNorm);
            }
            if (!isSolved)
              continue;

            //----------------------------------------------

            // Reject merge if any solved weight is non-finite (NaN / Inf).
            bool xIsFinite = true;
            for (int ip = 0; ip < nPartNew; ++ip) {
              if (!std::isfinite(x[ip])) {
                xIsFinite = false;
                break;
              }
            }
            if (!xIsFinite)
              continue;

            // Adjust weight.
            for (int ip = 0; ip < nPartNew; ip++) {
              auto& p = particles[idx_I[ip]];
              p.rdata(iqp_) = x[ip];
              merged[idx_I[ip]] = true;
            }
            // Mark for deletion
            for (int ip = nPartNew; ip < nOld; ip++) {
              particles[idx_I[ip]].id() = -1;
              particles[idx_I[ip]].rdata(iqp_) = 0;
              merged[idx_I[ip]] = true;
            }
          }
    }
  }
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::limit_weight_new(
    Real maxRatio, bool seperateVelocity) {
  timing_func("Pts::limit_weight");

  if (maxRatio <= 1)
    return;

  IntVect iv(1);
  if (!(do_tiling && tile_size == iv))
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      const auto tppc = target_PPC(iLev)[pti].array();
      const Box& bx = pti.tilebox();
      IntVect ibx = bx.smallEnd();
      int target = tppc(ibx);
      Vector<ParticleType> newparticles;
      auto& pTile = get_particle_tile(iLev, pti);
      AoS& particles = pti.GetArrayOfStructs();
      std::sort(particles.begin(), particles.end(), compare_two_parts);
      Real totalMass = 0;
      for (auto& p : particles) {
        totalMass += fabs(p.rdata(iqp_));
      }
      Real avg = totalMass / target;

      // Real maxWeight = avg + maxRatio * vars;
      Real maxWeight = avg * maxRatio;
      Real dl = 4.0 * Geom(iLev).CellSize()[ix_] / sqrt(tppc(ibx));
      {

        const auto lo = lbound(pti.tilebox());
        const auto hi = ubound(pti.tilebox());

        const Real xMin = Geom(iLev).LoEdge(lo.x, ix_) +
                          Geom(iLev).CellSize()[ix_] * 1e-10,
                   xMax = Geom(iLev).HiEdge(hi.x, ix_) -
                          Geom(iLev).CellSize()[ix_] * 1e-10;

        const Real yMin = Geom(iLev).LoEdge(lo.y, iy_) +
                          Geom(iLev).CellSize()[iy_] * 1e-10,
                   yMax = Geom(iLev).HiEdge(hi.y, iy_) -
                          Geom(iLev).CellSize()[iy_] * 1e-10;

        const Real zMin = nDim > 2 ? Geom(iLev).LoEdge(lo.z, iz_) +
                                         Geom(iLev).CellSize()[iz_] * 1e-10
                                   : 0.0,
                   zMax = nDim > 2 ? Geom(iLev).HiEdge(hi.z, iz_) -
                                         Geom(iLev).CellSize()[iz_] * 1e-10
                                   : 0.0;

        for (auto& p : particles) {
          Real qp1 = p.rdata(iqp_);
          if (fabs(qp1) < maxWeight)
            continue;

          Real xp1 = p.pos(ix_);
          Real yp1 = p.pos(iy_);
          Real zp1 = nDim > 2 ? p.pos(iz_) : 0;
          Real up1 = p.rdata(iup_);
          Real vp1 = p.rdata(ivp_);
          Real wp1 = p.rdata(iwp_);
          const Real u2 = up1 * up1 + vp1 * vp1 + wp1 * wp1;
          Real coef = (u2 < 1e-13) ? 0 : dl / sqrt(u2);
          p.rdata(iqp_) = qp1 / (2.0);
          const Real dpx = coef * up1;
          const Real dpy = coef * vp1;
          const Real dpz = coef * wp1;

          Real xp2 = xp1 + dpx;
          Real yp2 = yp1 + dpy;
          Real zp2 = zp1 + dpz;

          xp1 -= dpx;
          yp1 -= dpy;
          zp1 -= dpz;

          xp1 = std::clamp(xp1, xMin, xMax);
          yp1 = std::clamp(yp1, yMin, yMax);
          zp1 = std::clamp(zp1, zMin, zMax);

          p.pos(ix_) = xp1;
          p.pos(iy_) = yp1;

          if (nDim > 2)
            p.pos(iz_) = zp1;

          xp2 = std::clamp(xp2, xMin, xMax);
          yp2 = std::clamp(yp2, yMin, yMax);
          zp2 = std::clamp(zp2, zMin, zMax);

          auto pnew = make_particle();
          set_ids(pnew);

          pnew.pos(ix_) = xp2;
          pnew.pos(iy_) = yp2;
          if (nDim > 2)
            pnew.pos(iz_) = zp2;
          pnew.rdata(iup_) = up1;
          pnew.rdata(ivp_) = vp1;
          pnew.rdata(iwp_) = wp1;
          pnew.rdata(iqp_) = qp1 / (2.0);
          newparticles.push_back(pnew);
        }
      }
      for (auto& p : newparticles) {
        pTile.push_back(p);
      }
    }
  }
}
// Since Particles is a template, it is necessary to explicitly instantiate
// with template arguments.

// Explicit template instantiations.
template class Particles<nPicPartReal, nPicPartInt>;
template class Particles<nPTPartReal, nPTPartInt>;
