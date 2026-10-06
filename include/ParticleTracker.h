#ifndef _PARTICLETRACKER_H_
#define _PARTICLETRACKER_H_

#include <AMReX_Vector.H>

#include "Grid.h"
#include "Particles.h"
#include "Pic.h"
#include "TestParticles.h"

class ParticleTracker {
public:
  ParticleTracker(Grid &gridIn, FluidInterface *fluidIn, TimeCtr *tcIn, int id,
                  ParticleTrackerInfo &info, const DomainParameters & /*dp*/)
      : grid(gridIn),
        tc(tcIn),
        fi(fluidIn),
        pInfo(&info),
        gridID(id),
        nGst(gridIn.get_n_ghost()),
        isGridEmpty(gridIn.is_grid_empty_ref()),
        isNewGrid(gridIn.is_new_grid_ref()),
        cGrids(gridIn.box_arrays()),
        nGrids(gridIn.node_box_arrays()) {
    gridName = std::string("FLEKS") + std::to_string(gridID);
    printPrefix = gridName + " pt: ";
  }

  ~ParticleTracker();

  int n_lev() const { return grid.n_lev(); }
  int n_lev_max() const { return grid.n_lev_max(); }
  const amrex::Geometry &Geom(int iLev) const { return grid.Geom(iLev); }
  const amrex::DistributionMapping &DistributionMap(int iLev) const {
    return grid.get_dmap(iLev);
  }
  bool is_grid_empty() const { return grid.is_grid_empty(); }
  bool is_new_grid() const { return grid.is_new_grid(); }
  void is_new_grid(bool in) { grid.is_new_grid(in); }

  void post_process_param();

  void pre_regrid();
  void post_regrid();

  void update_field(Pic &pic, bool needJacobian = false);
  void set_ic(Pic &pic);
  void update(Pic &pic, bool doReport = false);

  void complete_parameters();

  void save_restart_data();
  void save_restart_header(std::ofstream &headerFile);
  void read_restart();
  void write_log(bool doForce = false, bool doCreateFile = false);

  void set_tp_init_shapes(amrex::Vector<std::shared_ptr<Shape> > &shapes);

private:
  Grid &grid;
  TimeCtr *tc = nullptr;
  FluidInterface *fi = nullptr;

  std::string tag = "pt";
  std::string gridName;
  std::string printPrefix;
  int gridID;
  int nGst;
  const bool &isGridEmpty;
  const bool &isNewGrid;
  const amrex::Vector<amrex::BoxArray> &cGrids;
  const amrex::Vector<amrex::BoxArray> &nGrids;

  // Parameter container populated by Domain during read_param and resolved in
  // ParticleTrackerInfo::post_process_param (after fi is fully processed).
  ParticleTrackerInfo *pInfo = nullptr;

  amrex::Vector<std::unique_ptr<TestParticles> > parts;
  amrex::Vector<amrex::MultiFab> nodeE;
  amrex::Vector<amrex::MultiFab> nodeB;

  // Nodal magnetic field Jacobian (9 components) computed via
  // jacobian_center_to_node.
  amrex::Vector<amrex::MultiFab> nodeJacB;

  // Scratch holding the total cell-centered field B1 + B0 when an intrinsic
  // field is configured; the Jacobian is the gradient of the total field.
  amrex::Vector<amrex::MultiFab> centerBtotal;

  std::unique_ptr<PlotCtr> savectr;

  std::string logFile;
  std::ofstream ptLogStream;

  // Test Particle initialization regions (set from the domain #REGION blocks).
  amrex::Vector<std::shared_ptr<Shape> > tpShapes;
};

#endif
