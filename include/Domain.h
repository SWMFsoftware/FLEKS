#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <memory>
#include <set>
#include <string>
#include <vector>

#include <AMReX_AmrMesh.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_RealBox.H>
#include <AMReX_Vector.H>

#include "DomainParameters.h"
#include "GridInfo.h"
#include "OHInterface.h"
#include "ParticleTracker.h"
#include "Pic.h"
#include "ReadParam.h"
#include "RefineRegions.h"
#include "SourceInterface.h"

enum class ParameterOwner { Pic, ParticleTracker, Source, FluidInterface };

class Domain {
private:
  // Inputs used to construct Grid; runtime geometry belongs to Grid.
  static constexpr int nGst = 2;
  static constexpr int coord = 0;
  amrex::Vector<int> nCell = { 1, 1, 1 };
  amrex::IntVect maxBlockSize = amrex::IntVect(8);
  amrex::IntVect periodicity = amrex::IntVect(0);
  amrex::RealBox domainBox;
  int refinementRatio = 2;
  bool isFake2D = false;

  void set_periodicity(int iDir, bool isPeriodic) {
    periodicity[iDir] = isPeriodic ? 1 : 0;
  }

  GridInfo gridInfo;
  int iGrid = 1;
  int iDecomp = 1;
  int gridID = 1;
  std::string printPrefix;
  std::string gridName;
  bool initialRegridPending = true;
  bool doSplitLevs = false;
  amrex::Vector<std::shared_ptr<Shape> > shapeList;
  std::set<std::string> shapeNames;

  void add_shape(const std::shared_ptr<Shape> &shape) {
    const std::string name = shape->get_name();
    if (!shapeNames.insert(name).second)
      amrex::Abort("Duplicate #REGION name: " + name);
    shapeList.push_back(shape);
  }
  DomainParameters domainParameters;

  ReadParam readParam;
  std::vector<std::string> sourceParameterCommands;

  ParticleTrackerInfo ptInfo;

  bool isTCInitialized = false;
  bool readingSessionUpdate = false;

  bool find_parameter_owner(const std::string &command,
                            ParameterOwner &owner) const;
  void prepare_grid_info(const amrex::Vector<double> &info,
                         amrex::Geometry &geometry, amrex::AmrInfo &amrInfo);

public:
  std::unique_ptr<TimeCtr> tc;

  // Single simulation grid owned by Domain. Declared before consumers so
  // it is destroyed after them.
  std::unique_ptr<Grid> grid;

  std::unique_ptr<Pic> pic;
  std::unique_ptr<ParticleTracker> pt;

  // fi stores EM field and/or the fluid quantities that are used to generate
  // particles.
  std::unique_ptr<FluidInterface> fi;

  std::unique_ptr<SourceInterface> source;

  // So far, stateOH is only used for OH-PT coupling. It stores the plasma
  // variables obtained from OH.
  std::unique_ptr<OHInterface> stateOH;
  // The Neu->Ion source for PT->OH coupling.
  std::unique_ptr<OHInterface> sourcePT2OH;
  //--------------------------------

  int couplerMarker;

public:
  Domain() { tc = std::make_unique<TimeCtr>(); }

  ~Domain() = default;

  void update();

  void calc_refine_region();

  //--------------Initialization begin-------------------------------
  void init(const double time = 0, const int iDomain = 1,
            const std::string &paramString = std::string(),
            const amrex::Vector<int> &paramInt = amrex::Vector<int>(),
            const amrex::Vector<double> &paramRegion = amrex::Vector<double>(),
            const amrex::Vector<double> &paramComm = amrex::Vector<double>(),
            const bool isStandalone = false);
  void update_param(const std::string &paramString);
  void set_ic();
  //----------------Initialization end-------------------------------

  //------------Coupler related begin--------------
  void set_state_var(double *data, int *index, std::vector<std::string> &names);

  int get_grid_nodes_number();

  void get_grid(double *pos_DI);

  void find_mpi_rank_for_points(const int nPoint, const double *const xyz_I,
                                int *const rank_I);

  void get_fluid_state_for_points(const int nDim, const int nPoint,
                                  const double *const xyz_I,
                                  double *const data_I, const int nVar);
  void get_source_for_points(const int nDim, const int nPoint,
                             const double *const xyz_I, double *const data_I,
                             const int nVar);
  //------------Coupler related end--------------

  //--------------- IO begin--------------------------------
  void read_param();

  // Read the commands that decide whether and how grid, fi, pic and pt are
  // constructed; the caller rolls back so read_param() reads the rest.
  void read_construction_parameters(ReadParam &param);
  void validate_configuration() const;
  void save_restart(std::string restartOutDir);
  void save_restart_header();
  void save_restart_data();
  void read_restart();
  void write_plots(bool doForce = false);
  //--------------- IO end----------------------------------

  //-------------- grid begin-------------------------------
  // Preparing grid information for Grid/AmrCore initialization.
  void regrid();
  void receive_grid_info(int *status = nullptr);
  //-------------- grid end---------------------------------

  void load_balance();

  // void make_data();
  void init_time_ctr();

  bool receive_ic_only() { return domainParameters.receiveICOnly; };

  int get_iGrid() const { return iGrid; }
  int get_iDecomp() const { return iDecomp; }
  int get_dim() const {
    return grid ? grid->get_dim() : ((isFake2D || nDim == 2) ? 2 : nDim);
  }
};

#endif
