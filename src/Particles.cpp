#include <cstdlib>

#include <AMReX_ParReduce.H>

#include "InitialCondition.h"
#include "Morton.h"
#include "Particles.h"
#include "SWMFInterface.h"
#include "Timer.h"
#include "Utility.h"

using namespace amrex;

namespace {

ParticlesInfo make_io_particles_info() {
  ParticlesInfo info;
  info.nPartPerCell = IntVect(AMREX_D_DECL(-1, -1, -1));
  return info;
}

} // namespace

//==========================================================
template <int NStructReal, int NStructInt>
Particles<NStructReal, NStructInt>::Particles(
    Grid* gridIn, FluidInterface* const fluidIn, TimeCtr* const tcIn,
    const int speciesIDIn, const Real chargeIn, const Real massIn,
    const ParticlesInfo& pInfo, const PartMode pModeIn,
    const InitialCondition* icIn)
    : AmrParticleContainer<NStructReal, NStructInt>(gridIn),
      grid(gridIn),
      fi(fluidIn),
      tc(tcIn),
      pMode(pModeIn),
      speciesID(speciesIDIn),
      charge(chargeIn),
      mass(massIn),
      nPartPerCell(pInfo.nPartPerCell),
      ic_(icIn) {

  isParticleLocationRandom = pInfo.isParticleLocationRandom;
  isPPVconstant = pInfo.isPPVconstant;
  doPreSplitting = pInfo.doPreSplitting;
  fastMerge = pInfo.fastMerge;
  mergeLight = pInfo.mergeLight;
  nPartCombine = pInfo.nPartCombine;
  nPartNew = pInfo.nPartNew;
  nMergeTry = pInfo.nMergeTry;
  mergeThresholdDistance = pInfo.mergeThresholdDistance;
  velBinBufferSize = pInfo.velBinBufferSize;
  mergeRatioMax = pInfo.mergeRatioMax;
  pLevRatio = pInfo.pLevRatio;
  mergePartRatioMax = pInfo.mergePartRatioMax;
  if (fi)
    vacuum = pInfo.vacuumIO * cProtonMassSI * 1e6 * fi->get_Si2NoRho();
  else
    vacuum = pInfo.vacuumIO;
  ionOH = pInfo.ionOH;
  bc = pInfo.particle_bc(speciesID);
  supID = pInfo.initial_sup_id(speciesID);
  do_tiling = true;

  qom = charge / mass;
  qomSign = qom >= 0 ? 1 : -1;

  absorbTallies.assign(18, 0.0);

  plo.resize(n_lev_max());
  phi.resize(n_lev_max());
  dx.resize(n_lev_max());
  invDx.resize(n_lev_max());
  invVol.resize(n_lev_max());

  for (int iLev = 0; iLev < n_lev_max(); iLev++) {
    for (int i = 0; i < nDim; ++i) {
      tile_size[i] = 1;
      plo[iLev][i] = Geom(iLev).ProbLo(i);
      phi[iLev][i] = Geom(iLev).ProbHi(i);
      dx[iLev][i] = Geom(iLev).CellSize(i);
      invDx[iLev][i] = Geom(iLev).InvCellSize(i);
    }
    invVol[iLev] = invDx[iLev].product();
  }

  isFake2D = (nDim == 3) &&
             (Geom(0).Domain().bigEnd(iz_) == Geom(0).Domain().smallEnd(iz_));

  // The following line is used to avoid an MPI bug (feature?) on Frontera. It
  // should be removed after the bug being fixed.
  SetUseUnlink(false);
}

IOParticles::IOParticles(Particles& other, Grid* gridIn, Real no2outL,
                         Real no2outV, Real no2outM, RealBox IORange)
    : Particles(gridIn, nullptr, nullptr, other.get_speciesID(),
                other.get_charge(), other.get_mass(), make_io_particles_info(),
                other.part_mode()) {

  no2outM *= qomSign * get_mass();

  const bool doLimit = IORange.ok();

  for (int iLev = 0; iLev < n_lev(); iLev++) {
    const auto& plevelOther = other.GetParticles(iLev);
    auto& plevel = GetParticles(iLev);
    for (MFIter mfi = other.MakeMFIter(iLev); mfi.isValid(); ++mfi) {
      auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());

      if (plevelOther.find(index) == plevelOther.end())
        continue;

      const auto& tileOther = plevelOther.at(index);

      if (tileOther.numParticles() == 0)
        continue;

      const AoS& aosOther = tileOther.GetArrayOfStructs();

      const Box& bx = other.cell_status(iLev)[mfi].box();
      const Array4<int const>& status = other.cell_status(iLev)[mfi].array();

      const IntVect lowCorner = bx.smallEnd();
      const IntVect highCorner = bx.bigEnd();

      for (auto p : aosOther) {
        if (other.is_outside_active_region(p, status, lowCorner, highCorner,
                                           iLev)) {
          // redistribute_particles() may fail if the ghost cell particles'
          // IDs are not -1 (marked for deletion);
          p.id() = -1;
        }

        for (int iDim = 0; iDim < nDim; iDim++) {
          p.pos(ix_ + iDim) = no2outL * p.pos(ix_ + iDim);
        }

        if (doLimit && !IORange.contains(RealVect(
                           AMREX_D_DECL(p.pos(ix_), p.pos(iy_), p.pos(iz_)))))
          continue;

        for (int iDim = 0; iDim < nDim3; iDim++) {
          p.rdata(iup_ + iDim) = no2outV * p.rdata(iup_ + iDim);
        }
        p.rdata(iqp_) = no2outM * p.rdata(iqp_);

        plevel[index].push_back(p);
      }
    }
  }
  redistribute_particles();
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::select_particle(
    Vector<std::array<int, 3> >& selectParticleIn) {

  timing_func("Pts::select_particle");

  int numParticlesFoundLocal = 0, numParticlesFoundTotal = 0;

  // output files
  std::string filename = "select_particle_out_sp" + std::to_string(speciesID) +
                         "_pe" + std::to_string(ParallelDescriptor::MyProc()) +
                         ".dat";
  std::ofstream outFile;
  outFile.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
  outFile.precision(12);

  // loop through particles
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      AoS& particles = pti.GetArrayOfStructs();
      for (auto& p : particles) {
        if (p.id() < 0)
          continue;
        if (p.idata(iSupID_) < 0)
          continue;

        for (auto& currentTargetParticle : selectParticleIn) {
          if (p.cpu() == currentTargetParticle[0] &&
              p.idata(iSupID_) == currentTargetParticle[1] &&
              p.id() == currentTargetParticle[2]) {
            numParticlesFoundLocal++;
            outFile << p.cpu() << " " << p.idata(iSupID_) << " " << p.id()
                    << " " << p.pos(ix_) << " " << p.pos(iy_) << " "
                    << p.pos(iz_) << " " << p.rdata(iup_) << " "
                    << p.rdata(ivp_) << " " << p.rdata(iwp_) << "\n";
          }
        }
      }
    }
  }
  outFile.close();
  MPI_Reduce(&numParticlesFoundLocal, &numParticlesFoundTotal, 1, MPI_INT,
             MPI_SUM, ParallelDescriptor::IOProcessorNumber(),
             ParallelDescriptor::Communicator());
  Print() << "select particle finished... " << numParticlesFoundTotal
          << "particles found..." << std::endl;
  amrex::Abort("Abort: select particle finished!");
}

//==========================================================

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::calculate_particle_quality(
    amrex::Vector<amrex::MultiFab>& quality) {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    quality[iLev].setVal(0.0);
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      Box bx = pti.tilebox();
      IntVect ibx = bx.smallEnd();
      const auto tppc = target_PPC(iLev)[pti].array();
      const auto qArr = quality[iLev][pti].array();
      auto& pTile = get_particle_tile(iLev, pti);
      AoS& particles = pTile.GetArrayOfStructs();
      Real totalMass = 0;
      for (auto& p : particles) {
        totalMass += fabs(p.rdata(iqp_));
      }
      Real perfectaverage = totalMass / tppc(ibx);

      for (auto& p : particles) {
        for (int aa = 0; aa <= 8; aa++) {
          if (fabs(p.rdata(iqp_)) > pow(2, aa + 1) * perfectaverage) {
            qArr(ibx, aa) += 1.0;
          }
        }
        for (int aa = 9; aa <= 17; aa++) {
          if (fabs(p.rdata(iqp_)) < perfectaverage / pow(2, aa - 8)) {
            qArr(ibx, aa) += 1.0;
          }
        }
      }
    }
  }
}

template <int NStructReal, int NStructInt>
void Particles<NStructReal, NStructInt>::add_velocity_perturbation(Real ampY,
                                                                   Real ampZ,
                                                                   Real kx) {
  for (int iLev = 0; iLev < n_lev(); iLev++) {
    for (PIter pti(*this, iLev); pti.isValid(); ++pti) {
      AoS& particles = pti.GetArrayOfStructs();
      for (auto& p : particles) {
        if (p.id() < 0)
          continue;
        Real x = p.pos(ix_);
        p.rdata(ivp_) += ampY * std::cos(kx * x);
        p.rdata(iwp_) += ampZ * std::sin(kx * x);
      }
    }
  }
}

// Explicit template instantiations.
template class Particles<nPicPartReal, nPicPartInt>;
template class Particles<nPTPartReal, nPTPartInt>;
