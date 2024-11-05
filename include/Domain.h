#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "DomainGrid.h"
#include "OHInterface.h"
#include "ParticleTracker.h"
#include "Pic.h"
#include "ReadParam.h"
#include "SourceInterface.h"

class Domain : public DomainGrid {
private:
  bool doRestart = false;

  bool doRestartPT = false;

  bool doRestartFIOnly = false;

  ReadParam param;

  // Number of files per AMREX output.
  int nFileField = 64, nFileParticle = 256;

  bool initFromSWMF = true;

  bool receiveICOnly = false;

  bool usePT = false;

public:
  // Q: Why are pic and pt defined as pointers?
  // A: Pic and particleTracker are derived from AmrCore, whose initialization
  // requires the informaion of the domain, such as the domain range and cell
  // size. Such information is not known untill Domain received information from
  // GM and read the input parameters. pic and pt are used as pointers so that
  // their initialization is defered until the grid information is obtained.
  std::unique_ptr<Pic> pic;
  std::unique_ptr<ParticleTracker> pt;

  // fi stores EM field and/or the fluid quantities that are used to generate
  // particles.
  std::unique_ptr<FluidInterface> fi;

  std::unique_ptr<SourceInterface> source;
  bool useFluidSource = false;

  // So far, stateOH is only used for OH-PT coupling. It stores the plasma
  // variables obtained from OH.
  std::unique_ptr<OHInterface> stateOH;
  // The Neu->Ion source for PT->OH coupling.
  std::unique_ptr<OHInterface> sourcePT2OH;

  std::unique_ptr<TimeCtr> tc;
  //--------------------------------

  int couplerMarker;

public:
  Domain() { tc = std::make_unique<TimeCtr>(); }

  ~Domain() = default;

  void update();

  //--------------Initialization begin-------------------------------
  void init(const double time = 0, const int iDomain = 1,
            const std::string &paramString = std::string(),
            const amrex::Vector<int> &paramInt = amrex::Vector<int>(),
            const amrex::Vector<double> &paramRegion = amrex::Vector<double>(),
            const amrex::Vector<double> &paramComm = amrex::Vector<double>());
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
  void read_param(const bool readGridInfo);
  void save_restart(std::string restartOutDir);
  void save_restart_header();
  void save_restart_data();
  void read_restart();
  void write_plots(bool doForce = false);
  //--------------- IO end----------------------------------

  //-------------- grid begin-------------------------------
  // Preparing grid information for Grid/AmrCore initialization.
  void prepare_grid_info(const amrex::Vector<double> &info);
  void regrid();
  void receive_grid_info(int *status = nullptr);
  //-------------- grid end---------------------------------

  void load_balance();

  // void make_data();
  void init_time_ctr();

  bool receive_ic_only() { return receiveICOnly; };
};

#endif
