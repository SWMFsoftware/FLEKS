#ifndef _PARTICLETRACKER_H_
#define _PARTICLETRACKER_H_

#include <AMReX_Vector.H>

#include "Grid.h"
#include "Pic.h"
#include "TestParticles.h"

class ParticleTracker : public Grid {
public:
  ParticleTracker(const amrex::RealBox &rb, const amrex::Vector<int> &nCell,
                  int coord = 0, int nLevel = 0, const int *isPer = nullptr)
      : Grid(rb, nCell, coord, nLevel, isPer){};

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
};

#endif
