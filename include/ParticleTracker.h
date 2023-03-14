#ifndef _PARTICLETRACKER_H_
#define _PARTICLETRACKER_H_

#include <AMReX_Vector.H>

#include "Grid.h"
#include "Pic.h"
#include "TestParticles.h"

class ParticleTracker : public Grid {
public:
  ParticleTracker(amrex::Geometry const &gm, amrex::AmrInfo const &amrInfo,
                  int nGst, std::shared_ptr<FluidInterface> &fluidIn,
                  std::shared_ptr<TimeCtr> &tcIn, int id)
      : Grid(gm, amrInfo, nGst, id), tc(tcIn), fi(fluidIn) {
    nSpecies = fi->get_nS();
  };

  ~ParticleTracker() {
    if (!isGridInitialized)
      return;

    bool doSave = savectr->is_time_to(true);
    for (auto &tps : parts) {
      if (doSave) {
        tps->write_particles(tc->get_cycle());
      }
    }
  };

  void init(std::shared_ptr<FluidInterface> &fluidIn,
            std::shared_ptr<TimeCtr> &tcIn, int id = 0);

  void post_process_param();

  void regrid(const amrex::BoxArray &activeRegionBAIn,
              const amrex::BoxArray &centerBAIn,
              const amrex::DistributionMapping &dmIn, Pic &pic);

  void update_field(Pic &pic);
  void update_cell_status(Pic &pic);
  void set_ic(Pic &pic);
  void update(Pic &pic);

  void complete_parameters();

  void save_restart_data();
  void save_restart_header(std::ofstream &headerFile);
  void read_restart();
  void read_param(const std::string &command, ReadParam &param);
  void write_log(bool doForce = false, bool doCreateFile = false);

private:
  bool usePT = false;

  std::shared_ptr<FluidInterface> fi;
  std::shared_ptr<TimeCtr> tc;

  int nSpecies;
  amrex::Vector<std::unique_ptr<TestParticles> > parts;
  amrex::MultiFab nodeE;

  amrex::Vector<amrex::MultiFab> nodeB;

  amrex::Vector<unsigned long int> initPartNumber;

  std::unique_ptr<PlotCtr> savectr;
  int dnSave = 1;

  amrex::IntVect nTPPerCell = { 1, 1, 1 };
  amrex::IntVect nTPIntervalCell = { 1, 1, 1 };

  std::string sPartRegion;

  std::string sIOUnit = "planet";

  bool isRelativistic = false;

  amrex::Vector<std::string> listFiles;
  bool doInitFromPIC = false;

  amrex::Vector<Vel> tpStates;

  std::string logFile;
};

#endif
