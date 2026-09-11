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
void Particles<NStructReal, NStructInt>::sample_charge_exchange(
    Real* vp, Real* vh, Real* up, Real vth, CrossSection cs) {

  timing_func("Pts::sample_charge_exchange");

  // M (holds normalization costants for both distributions), g(vp) is the
  // maxwellian distribution
  Real sup[3] = { up[0] - 3. * vth, up[1] - 3. * vth, up[2] - 3. * vth };
  Real M = charge_exchange_dis(sup, vh, up, vth, cs) /
           exp(-((sup[0] - up[0]) * (sup[0] - up[0]) +
                 (sup[1] - up[1]) * (sup[1] - up[1]) +
                 (sup[2] - up[2]) * (sup[2] - up[2])) /
               (vth * vth));

  bool accepted = false;
  while (!accepted) {
    {
      Real prob, theta, uth;
      // u = X velocity
      prob = sqrt(-2.0 * log(1.0 - .999999999 * randNum()));
      theta = 2.0 * M_PI * randNum();
      uth = vth / sqrt(2.0);
      vp[0] = uth * prob * cos(theta) + up[0];
      // v = Y velocity
      vp[1] = uth * prob * sin(theta) + up[1];
      // w = Z velocity
      prob = sqrt(-2.0 * log(1.0 - .999999999 * randNum()));
      theta = 2.0 * M_PI * randNum();
      vp[2] = uth * prob * cos(theta) + up[2];
    }

    if (randNum() < charge_exchange_dis(vp, vh, up, vth, cs) /
                        (M * exp(-((vp[0] - up[0]) * (vp[0] - up[0]) +
                                   (vp[1] - up[1]) * (vp[1] - up[1]) +
                                   (vp[2] - up[2]) * (vp[2] - up[2])) /
                                 (vth * vth)))) {
      accepted = true;
    }
  }
}

template <int NStructReal, int NStructInt>
Real Particles<NStructReal, NStructInt>::charge_exchange_dis(Real* vp, Real* vh,
                                                             Real* up, Real vth,
                                                             CrossSection cs) {
  Real dv_D[3], dv2 = 0, dv = 0;
  for (int i = 0; i < 3; ++i) {
    dv_D[i] = vh[i] - vp[i];
    dv2 += dv_D[i] * dv_D[i];
  }

  if (dv2 == 0)
    return 0.0;

  dv = sqrt(dv2);

  Real erel = 0.5 * 1.674E-27 * dv2 * 6.2415E15; // in keV

  Real sigma = 0;
  if (cs == CrossSection::LS) {
    sigma = (4.15 - 0.531 * log(erel)) * (4.15 - 0.531 * log(erel)) *
            pow(1 - exp(-67.3 / erel), 4.5) * 1E-20; // cross section in m^2
  } else if (cs == CrossSection::MT) {
    Real dvcm = dv * 1E2;                             // velocity in cm/s
    sigma = pow(1.6 - 0.0695 * log(dvcm), 2) * 1e-18; // cross section in m^2
  }

  Real dvpup2 = 0;
  for (int i = 0; i < 3; ++i) {
    dvpup2 += pow(vp[i] - up[i], 2);
  }

  return dv * sigma * exp(-dvpup2 / (vth * vth));
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::get_analytic_ion_fluid(
    const RealVect xyz, Real& rhoIon, Real& cs2Ion, Real (&uIon)[nDim3]) {

  // The units of ionOH are assumed to be:
  // r: AU
  // rho: amu/cc
  // T: K
  // U: km/s

  Real r = max(xyz.vectorLength(), 1e-9 * Geom(0).CellSize(0));

  Real rAU = r * fi->get_No2SiL() / cAUSI;

  if (rAU > ionOH.rAnalytic) {
    Abort("Error: rAU > ionOH.rAnalytic");
  }

  if (ionOH.doGetFromOH) {
    Real xSI = xyz[ix_] * fi->get_No2SiL();
    Real ySI = xyz[iy_] * fi->get_No2SiL();
    Real zSI = (nDim > 2 ? xyz[iz_] : 0.0) * fi->get_No2SiL();

    Real temp, ur, b[nDim3];
    OH_get_solar_wind(&xSI, &ySI, &zSI, &rhoIon, &ur, &temp, b);

    // v_th = sqrt(2kT/m); m/s
    cs2Ion = 2 * cBoltzmannSI * temp / cProtonMassSI; // m^2/s^2

    for (int i = 0; i < nDim; ++i) {
      uIon[i] = ur * xyz[i] / r;
    }
  } else {
    Real r0 = 0;
    if (rAU < ionOH.rCutoff) {
      r0 = ionOH.rAnalytic / ionOH.rCutoff;
    } else {
      r0 = ionOH.rAnalytic / rAU;
    }

    rhoIon = ionOH.swRho * pow(r0, 2) * 1e6; // amu/cc -> amu/m^3

    // v_th = sqrt(2kT/m); m/s
    cs2Ion = 2 * cBoltzmannSI * ionOH.swT / cProtonMassSI; // m^2/s^2

    for (int i = 0; i < nDim; ++i) {
      uIon[i] = ionOH.swU * xyz[i] / r * 1e3; // km/s -> m/s
    }
  }

  // AllPrint() << "r = " << r << " xyz = " << xyz[0] << ", " << xyz[1] << ", "
  //            << xyz[2] << " rhoIon = " << rhoIon << ", cs2Ion = " << cs2Ion
  //            << ", uIon = " << uIon[0] << ", " << uIon[1] << ", " << uIon[2]
  //            << std::endl;
}

// Get the iFluid-th ion fluid properties at the location xyz

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::get_ion_fluid(
    FluidInterface* stateOH, PIter& pti, const int iLev, const int iFluid,
    const RealVect xyz, Real& rhoIon, Real& cs2Ion, Real (&uIon)[nDim3]) {

  Real rAU = xyz.vectorLength() * stateOH->get_No2SiL() / cAUSI;

  if (iFluid == 0 && rAU < ionOH.rAnalytic) {
    get_analytic_ion_fluid(xyz, rhoIon, cs2Ion, uIon);
    return;
  }

  // amu/m^3
  rhoIon = stateOH->get_fluid_mass_density(pti, xyz, iFluid, iLev) *
           stateOH->get_No2SiRho() / cProtonMassSI;

  // cs = sqrt(P/n); m/s
  // Assume p = pi + pe = 2pi, so divide by sqrt(2.0).
  Real cs = stateOH->get_fluid_uth(pti, xyz, iFluid, iLev) *
            stateOH->get_No2SiV() / sqrt(2.0);

  // cs2Ion = 2*P/n. The definition of thermal speed in get_uth_iso() is
  // different from the requirement in OH_get_charge_exchange_wrapper().
  // See page 92 of Adam Michael's thesis.
  cs2Ion = 2 * pow(cs, 2);

  uIon[ix_] =
      stateOH->get_fluid_ux(pti, xyz, iFluid, iLev) * stateOH->get_No2SiV();
  uIon[iy_] =
      stateOH->get_fluid_uy(pti, xyz, iFluid, iLev) * stateOH->get_No2SiV();
  uIon[iz_] =
      stateOH->get_fluid_uz(pti, xyz, iFluid, iLev) * stateOH->get_No2SiV();
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::charge_exchange(
    Real dt, FluidInterface* stateOH, FluidInterface* sourcePT2OH,
    SourceInterface* source, bool kineticSource,
    Vector<std::unique_ptr<PicParticles> >& sourceParts, bool doSelectRegion,
    int nppc, Real& maxExchangeRatio) {
  std::string nameFunc = "Pts::charge_exchange";

  timing_func(nameFunc);

  if (dt <= 0)
    return;

  // for (auto& ptr : sourceParts) {
  //   ptr->clearParticles();
  // }

  struct NeuPlasmaPair {
    Real q; // weight
    Real vp[3];
    Real vh[3];
    Real up[3];
    Real xyz[3];
    Real vth;
  };

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      AoS& particles = pti.GetArrayOfStructs();

      if (kineticSource) {
        // Sort the particles by the location first to make sure the results
        // are the same for different number of processors
        std::sort(particles.begin(), particles.end(), compare_two_parts);
      }

      // It is assumed the tile size is 1x1x1.
      Box bx = pti.tilebox();
      auto cellIdx = bx.smallEnd();

      int iRegion = 0;
      if (kineticSource && doSelectRegion) {
        iRegion = stateOH->get_neu_source_region(pti, cellIdx, iLev);
        if (iRegion < 0)
          continue;
      }
      // ParticleTileType
      auto& spTile =
          sourceParts[iRegion]->get_particle_tile(iLev, pti, cellIdx);

      Vector<NeuPlasmaPair> neuPlasmaPairs;

      if (kineticSource) {
        set_random_seed(iLev, cellIdx, IntVect(999));
      }

      for (auto& p : particles) {
        if (p.id() < 0)
          continue;

        RealVect xyz;
        for (int i = 0; i < nDim; ++i) {
          xyz[i] = p.pos(i);
        }

        Real cs2Neu = 0, uNeu[3], rhoNeu;
        Real cs2Ion, uIon[3], rhoIon;
        Real ion2neu[5], neu2ion[5];
        const int iRho_ = 0, iUx_ = 1, iUy_ = 2, iUz_ = 3, iP_ = 4;
        const int iRhoUx_ = iUx_, iRhoUy_ = iUy_, iRhoUz_ = iUz_, iE_ = iP_;

        // amu/m^3
        rhoNeu = qomSign * p.rdata(iqp_) * get_mass() * invVol[iLev] *
                 stateOH->get_No2SiRho() / cProtonMassSI;

        for (int i = 0; i < nDim; ++i) {
          uNeu[i] = p.rdata(iup_ + i) * stateOH->get_No2SiV();
        }

        // A neutral particle interacts with all the ion fluids.
        for (int fluidID = 0; fluidID < stateOH->get_nFluid(); fluidID++) {
          get_ion_fluid(stateOH, pti, iLev, fluidID, xyz, rhoIon, cs2Ion, uIon);

          OH_get_charge_exchange_wrapper(&rhoIon, &cs2Ion, uIon, &rhoNeu,
                                         &cs2Neu, uNeu, ion2neu, neu2ion);

          // The function above returns number density changing rate.
          ion2neu[iRho_] *= cProtonMassSI;
          neu2ion[iRho_] *= cProtonMassSI;

          Real dtSI = dt * stateOH->get_No2SiT();
          // Print() << "rhoion = " << rhoIon << " cs2Ion = " << cs2Ion
          //         << " rhoNeu = " << rhoNeu << " cs2Neu = " << cs2Neu
          //         << " dtSI = " << dtSI << std::endl;
          for (int i = iRho_; i <= iP_; ++i) {
            ion2neu[i] *= dtSI;
            neu2ion[i] *= dtSI;
            // Print() << " i = " << i << " ion2neu = " << ion2neu[i]
            //         << " neu2ion = " << neu2ion[i] << std::endl;
          }

          Real massExchange =
              neu2ion[iRho_] * stateOH->get_Si2NoRho() / invVol[iLev];

          if (massExchange == 0) {
            // It can happen for some special cases. For example, when neutral
            // density is zero.
            continue;
          }

          // Print() << "nden = " << p.rdata(iqp_)
          //         << " massExchange = " << massExchange << std::endl;
          if (p.rdata(iqp_) - massExchange <= 0) {
            // Mark for deletion
            p.id() = -1;

            // Reduce the sources accordingly to conserve total masses.
            const Real ratio = p.rdata(iqp_) / massExchange;
            for (int i = iRho_; i <= iP_; ++i) {
              ion2neu[i] *= ratio;
              neu2ion[i] *= ratio;
            }
          } else {
            // Reduce particle mass due to charge exchange

            Real ratio = massExchange / p.rdata(iqp_);
            if (ratio > maxExchangeRatio)
              maxExchangeRatio = ratio;

            p.rdata(iqp_) = p.rdata(iqp_) - massExchange;
          }

          {
            int iFluidAddTo = 0;
            Real rhoIonAddTo, cs2IonAddTo, uIonAddTo[3];

            const int iSW = 2;
            const int iSheath = 1;
            const int iOutSheath = 0;
            switch (stateOH->get_nFluid()) {
              case 1:
                iFluidAddTo = fluidID;
                rhoIonAddTo = rhoIon;
                break;
              case 2:
                if (iRegion == iSW) {
                  iFluidAddTo = 1; // Pu3
                } else {
                  iFluidAddTo = 0; // background
                }
                break;
              case 3:
                if (iRegion == iSW) {
                  iFluidAddTo = 1; // Pu3
                } else if (iRegion == iSheath || iRegion == iOutSheath) {
                  iFluidAddTo = 2; // Pu2
                } else {
                  iFluidAddTo = 0; // background
                }
                break;
              default:
                Abort("Error: nFluid > 3 is not supported yet.");
            }

            if (stateOH->get_nFluid() > 1) {
              get_ion_fluid(stateOH, pti, iLev, iFluidAddTo, xyz, rhoIonAddTo,
                            cs2IonAddTo, uIonAddTo);
            }

            // Q: Why is (neu2ion-ion2neu) divided by rhoIon?
            // A: What passed between PT and OH is 'source per ion density'
            // instead of source. The ion density will be multiplied back in OH
            // ModUser.f90
            { // Sources for ion fluid: Neu -> Ion
              sourcePT2OH->add_rho_to_loc(neu2ion[iRho_] / rhoIonAddTo, pti,
                                          xyz, iFluidAddTo, iLev);
              sourcePT2OH->add_mx_to_loc(neu2ion[iRhoUx_] / rhoIonAddTo, pti,
                                         xyz, iFluidAddTo, iLev);
              sourcePT2OH->add_my_to_loc(neu2ion[iRhoUy_] / rhoIonAddTo, pti,
                                         xyz, iFluidAddTo, iLev);
              sourcePT2OH->add_mz_to_loc(neu2ion[iRhoUz_] / rhoIonAddTo, pti,
                                         xyz, iFluidAddTo, iLev);
              sourcePT2OH->add_p_to_loc(neu2ion[iP_] / rhoIonAddTo, pti, xyz,
                                        iFluidAddTo, iLev);
            }

            { // Loses for ion fluid: Ion -> Neu
              sourcePT2OH->add_rho_to_loc(-ion2neu[iRho_] / rhoIon, pti, xyz,
                                          fluidID, iLev);
              sourcePT2OH->add_mx_to_loc(-ion2neu[iRhoUx_] / rhoIon, pti, xyz,
                                         fluidID, iLev);
              sourcePT2OH->add_my_to_loc(-ion2neu[iRhoUy_] / rhoIon, pti, xyz,
                                         fluidID, iLev);
              sourcePT2OH->add_mz_to_loc(-ion2neu[iRhoUz_] / rhoIon, pti, xyz,
                                         fluidID, iLev);
              sourcePT2OH->add_p_to_loc(-ion2neu[iP_] / rhoIon, pti, xyz,
                                        fluidID, iLev);
            }
          }

          if (ion2neu[iRho_] > 0) { // Add source to nodes.

            if (kineticSource) {

              NeuPlasmaPair pair;
              pair.q = massExchange;
              for (int i = 0; i < nDim3; ++i) {
                pair.vh[i] = uNeu[i];
                pair.up[i] = uIon[i];
                pair.vth = sqrt(cs2Ion);
              }

              for (int i = 0; i < nDim; ++i) {
                pair.xyz[i] = p.pos(ix_ + i);
              }

              neuPlasmaPairs.push_back(pair);

            } else {

              Real si2no_v[5];
              si2no_v[iRho_] = source->get_Si2NoRho();
              si2no_v[iRhoUx_] = source->get_Si2NoV() * si2no_v[iRho_];
              si2no_v[iRhoUy_] = si2no_v[iRhoUx_];
              si2no_v[iRhoUz_] = si2no_v[iRhoUx_];
              si2no_v[iP_] = source->get_Si2NoP();

              Real m2 = 0;
              for (int i = iRhoUx_; i <= iRhoUz_; ++i) {
                m2 += pow(ion2neu[i], 2);
              }

              // P = (gamma-1)*(E - 0.5*rho*u2)
              ion2neu[iP_] =
                  (gamma0 - 1) * (ion2neu[iE_] - 0.5 * m2 / ion2neu[iRho_]);

              if (ion2neu[iP_] < 0) {
                ion2neu[iP_] = 0;
              }

              // source saves changing rate (density/s...).
              source->add_rho_to_loc(ion2neu[iRho_] * si2no_v[iRho_] / dt, pti,
                                     xyz, speciesID, iLev);
              source->add_mx_to_loc(ion2neu[iRhoUx_] * si2no_v[iRhoUx_] / dt,
                                    pti, xyz, speciesID, iLev);
              source->add_my_to_loc(ion2neu[iRhoUy_] * si2no_v[iRhoUy_] / dt,
                                    pti, xyz, speciesID, iLev);
              source->add_mz_to_loc(ion2neu[iRhoUz_] * si2no_v[iRhoUz_] / dt,
                                    pti, xyz, speciesID, iLev);
              source->add_p_to_loc(ion2neu[iP_] * si2no_v[iP_] / dt, pti, xyz,
                                   speciesID, iLev);
            }
          }
          // p.id() = -1;
        }
      } // for p

      if (kineticSource) {
        // Sample the velocity distribution function.

        Vector<NeuPlasmaPair> newPairs;

        Real wt = 0;

        Vector<Real> weights;
        weights.resize(neuPlasmaPairs.size());
        for (int i = 0; i < neuPlasmaPairs.size(); ++i) {
          weights[i] = neuPlasmaPairs[i].q;
          wt += neuPlasmaPairs[i].q;
        }

        std::vector<int> idx = random_select_weighted_n(weights, nppc, randNum);

        Real wtnew = 0;
        for (int i : idx) {
          newPairs.push_back(neuPlasmaPairs[i]);
          wtnew += neuPlasmaPairs[i].q;
        }
        Real scale = wt / wtnew;

        for (auto& pair : newPairs) {
          sample_charge_exchange(pair.vp, pair.vh, pair.up, pair.vth,
                                 CrossSection::MT);
          for (int i = 0; i < nDim3; ++i) {
            pair.vp[i] *= stateOH->get_Si2NoV();
          }

          auto newp = sourceParts[iRegion]->make_particle();
          newp.rdata(iqp_) = pair.q * scale;

          for (int i = 0; i < nDim3; ++i) {
            newp.rdata(iup_ + i) = pair.vp[i];
          }

          for (int i = 0; i < nDim; ++i) {
            newp.pos(ix_ + i) = pair.xyz[i];
          }

          spTile.push_back(newp);
        }
      }
    } // for pti
  }

  ParallelDescriptor::ReduceRealMax(maxExchangeRatio);

  Print() << "maxExchangeRatio = " << maxExchangeRatio << std::endl;
  if (maxExchangeRatio > 0.2) {
    Print() << "Warning: maybe the charge exchange rate is too high within one "
               "time step! Reducing the time step will help to slow down the "
               "charge exchange."
            << std::endl;
  }
}

//==========================================================
// Apply chemical loss (recombination, etc.) by proportionally
// reducing particle weights.  For each cell, reads the ion loss
// rate from source->nodeLossFluid and the existing ion mass density
// from fi, then reduces every particle's weight by the fraction
//   fraction = min(lossRate * dt / rhoExisting, 1.0).
// In hybrid PIC mode, only ion species (charge > 0) have particles;
// in full PIC mode, electron particle weights are also reduced proportionally.

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::apply_loss(
    const SourceInterface* source, Real dt) {
  std::string nameFunc = "Pts::apply_loss";
  timing_func(nameFunc);

  if (!source || !source->use_loss_source())
    return;

  if (speciesID < 0 || speciesID >= fi->get_nS())
    return;

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    if (NumberOfParticlesAtLevel(iLev, true, true) == 0)
      continue;

    if (!source->has_loss_array(iLev))
      continue;

#ifdef AMREX_USE_GPU
    const auto plo_geom = Geom(iLev).ProbLoArray();
    const auto inv_dx_geom = Geom(iLev).InvCellSizeArray();
    const int spID = speciesID;
    const int iqp = iqp_;
    const int iRho = fi->get_iRho(speciesID);
    const int ndim = nDim;

    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      AoS& particles = pti.GetArrayOfStructs();
      const int np = particles.numParticles();
      if (np == 0)
        continue;

      auto p_ptr = particles.data();
      const auto rhoArr = fi->get_node_fluid(iLev)[pti].array();
      const auto lossArr = source->get_node_loss_fluid(iLev)[pti].array();

      amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(int ip) noexcept {
        auto& p = p_ptr[ip];
        if (p.id() < 0)
          return;

        int i = static_cast<int>(
            std::floor((p.pos(0) - plo_geom[0]) * inv_dx_geom[0]));
        int j = (ndim > 1) ? static_cast<int>(std::floor(
                                 (p.pos(1) - plo_geom[1]) * inv_dx_geom[1]))
                           : 0;
        int k = (ndim > 2) ? static_cast<int>(std::floor(
                                 (p.pos(2) - plo_geom[2]) * inv_dx_geom[2]))
                           : 0;

        Real rhoExisting = rhoArr(i, j, k, iRho);
        if (rhoExisting <= 0.0)
          return;

        Real lossRate = lossArr(i, j, k, spID);
        if (lossRate <= 0.0)
          return;

        Real fraction = lossRate * dt / rhoExisting;
        if (fraction > 1.0)
          fraction = 1.0;
        if (fraction <= 0.0)
          return;

        p.rdata(iqp) *= (1.0 - fraction);
        if (fraction >= 1.0) {
          p.id() = -1;
        }
      });
    }
#else
    const Real* plo = Geom(iLev).ProbLo();
    const Real* inv_dx = Geom(iLev).InvCellSize();

    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      AoS& particles = pti.GetArrayOfStructs();
      const int np = particles.numParticles();
      if (np == 0)
        continue;

      for (int ip = 0; ip < np; ++ip) {
        ParticleType& p = particles[ip];
        if (p.id() < 0)
          continue;

        // Find the cell index for this particle.
        IntVect ijk;
        for (int iDim = 0; iDim < nDim; ++iDim) {
          ijk[iDim] = static_cast<int>(
              std::floor((p.pos(iDim) - plo[iDim]) * inv_dx[iDim]));
        }

        // Existing mass density (normalized) from the plasma state.
        Real rhoExisting =
            fi->get_fluid_mass_density(pti, ijk, speciesID, iLev);
        if (rhoExisting <= 0.0)
          continue;

        // Loss rate (normalized mass-density rate) from nodeLossFluid.
        Real lossRate = source->get_loss_value(pti, ijk, speciesID, iLev);
        if (lossRate <= 0.0)
          continue;

        // Fraction of mass to remove in this step.
        Real fraction = lossRate * dt / rhoExisting;
        if (fraction > 1.0)
          fraction = 1.0;
        if (fraction <= 0.0)
          continue;

        // Reduce particle weight proportionally.  The sign is preserved
        // (ions have positive weight, electrons negative).
        p.rdata(iqp_) *= (1.0 - fraction);
        if (fraction >= 1.0) {
          p.id() = -1;
        }
      }
    }
#endif
  }

  // Remove particles whose weight has been driven to (near) zero.
  redistribute_particles();
}

// Explicit template instantiations.
template class Particles<nPicPartReal, nPicPartInt>;
template class Particles<nPTPartReal, nPTPartInt>;
