#ifndef _PARTICLETRACKER_H_
#define _PARTICLETRACKER_H_

#include <AMReX_Vector.H>

#include "Pic.h"
#include "TestParticles.h"

class ParticleTracker {
public:
  ParticleTracker(){};
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
            std::shared_ptr<TimeCtr> &tcIn, int domainIDIn = 0);

  void post_process_param();

  void set_geom(int nGstIn, const amrex::Geometry &geomIn);

  void regrid(const amrex::BoxArray &ptRegionBAIn,
              const amrex::BoxArray &centerBAIn,
              const amrex::DistributionMapping &dmIn, Pic &pic);

  bool is_grid_empty() const { return isGridEmpty; }

  void update_field(Pic &pic);
  void update_cell_status(Pic &pic);
  void set_ic(Pic &pic);
  void update(Pic &pic);

  void complete_parameters();

  void save_restart_data();
  void save_restart_header(std::ofstream &headerFile);
  void read_restart();
  void read_param(const std::string &command, ReadParam &readParam);

private:
  bool usePT = false;

  std::string printPrefix;
  std::string domainName;
  int domainID;

  std::shared_ptr<FluidInterface> fluidInterface;
  std::shared_ptr<TimeCtr> tc;

  int nSpecies;
  amrex::Vector<std::unique_ptr<TestParticles> > parts;
  amrex::MultiFab nodeE, nodeB;

  amrex::Vector<unsigned long int> initPartNumber;

  int nGst;
  amrex::DistributionMapping dm;
  amrex::Geometry gm;
  amrex::BoxArray nodeBA, centerBA;

  // A collection of boxes to describe the PT domain. The boxes have been
  // combined if possible. It covers the same region as centerBA, but usually
  // contains less boxes.
  amrex::BoxArray ptRegionBA;

  amrex::iMultiFab cellStatus;

  std::unique_ptr<PlotCtr> savectr;
  int dnSave = 1;

  bool isGridInitialized = false;

  bool isGridEmpty = false;

  amrex::IntVect nTPPerCell = { 1, 1, 1 };
  amrex::IntVect nTPIntervalCell = { 1, 1, 1 };

  std::string sPartRegion;

  std::string sIOUnit = "planet";

  bool isRelativistic = false;

  amrex::Vector<std::string> listFiles;
  bool doInitFromPIC = false;

  amrex::Vector<Vel> tpStates;
};

#endif
