#ifndef _PARTICLETRACKER_H_
#define _PARTICLETRACKER_H_

#include <AMReX_Vector.H>

#include "Grid.h"
#include "Pic.h"
#include "TestParticles.h"

class ParticleTracker : public Grid {
public:
  ParticleTracker(amrex::Geometry const &gm, amrex::AmrInfo const &amrInfo,
                  int nGst, FluidInterface *fluidIn, TimeCtr *tcIn, int id)
      : Grid(gm, amrInfo, nGst, id, "pt"), tc(tcIn), fi(fluidIn) {
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

  void post_process_param();

  void pre_regrid() override;
  void post_regrid() override;

  void update_field(Pic &pic);
  void set_ic(Pic &pic);
  void update(Pic &pic, bool doReport = false);

  void complete_parameters();

  void save_restart_data();
  void save_restart_header(std::ofstream &headerFile);
  void read_restart();
  void read_param(const std::string &command, ReadParam &param);
  void write_log(bool doForce = false, bool doCreateFile = false);

private:
  TimeCtr *tc = nullptr;
  FluidInterface *fi = nullptr;

  int nSpecies;
  amrex::Vector<std::unique_ptr<TestParticles> > parts;
  amrex::Vector<amrex::MultiFab> nodeE;

  amrex::Vector<amrex::MultiFab> nodeB;

  amrex::Vector<unsigned long int> initPartNumber;

  std::unique_ptr<PlotCtr> savectr;
  int dnSave = 1;

  amrex::IntVect nTPPerCell = { AMREX_D_DECL(1, 1, 1) };
  amrex::IntVect nTPIntervalCell = { AMREX_D_DECL(1, 1, 1) };

  std::string sPartRegion;

  std::string sIOUnit = "planet";

  bool isRelativistic = false;

  amrex::Vector<std::string> listFiles;
  bool doInitFromPIC = false;

  amrex::Vector<Vel> tpStates;

  std::string logFile;
};

#endif
