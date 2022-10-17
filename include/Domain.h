#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "DomainGrid.h"
#include "ParticleTracker.h"
#include "Pic.h"
#include "ReadParam.h"

class Domain : public DomainGrid {
private:
  bool doRestart = false;

  ReadParam param;

  // Number of files per AMREX output.
  int nFileField = 64, nFileParticle = 256;

public:
  // Q: Why are pic and pt defined as pointers?
  // A: Pic and particleTracker are derived from AmrCore, whose initialization
  // requires the informaion of the domain, such as the domain range and cell
  // size. Such information is not known untill Domain received information from
  // GM and read the input parameters. pic and pt are used as pointers so that
  // their initialization is defered until the grid information is obtained.
  std::unique_ptr<Pic> pic;
  std::unique_ptr<ParticleTracker> pt;

  // Conceptually, both the Domain class and the Pic class may use the
  // following classes, so they are handled by shared pointers.
  std::shared_ptr<FluidInterface> fi;
  std::shared_ptr<TimeCtr> tc;

  int couplerMarker;

public:
  Domain() { tc = std::make_shared<TimeCtr>(); }

  ~Domain() = default;

  void update();

  //--------------Initialization begin-------------------------------
  void init(double time, const std::string &paramString, int *paramInt,
            double *gridDim, double *paramReal, int iDomain = 1);
  void update_param(const std::string &paramString);
  void set_ic();
  //----------------Initialization end-------------------------------

  //------------Coupler related begin--------------
  void set_state_var(double *data, int *index);

  int get_grid_nodes_number();

  void get_grid(double *pos_DI);

  void find_mpi_rank_for_points(const int nPoint, const double *const xyz_I,
                                int *const rank_I);

  void get_fluid_state_for_points(const int nDim, const int nPoint,
                                  const double *const xyz_I,
                                  double *const data_I, const int nVar);
  //------------Coupler related end--------------

  //--------------- IO begin--------------------------------
  void read_param(const bool readGridInfo);
  void save_restart();
  void save_restart_header();
  void save_restart_data();
  void read_restart();
  void write_plots(bool doForce = false);
  //--------------- IO end----------------------------------

  //-------------- grid begin-------------------------------
  // Preparing grid information for Grid/AmrCore initialization.
  void prepare_grid_info(const double *const info);
  void regrid();
  void receive_grid_info(int *status);
  //-------------- grid end---------------------------------

  // void make_data();
  void init_time_ctr();
};

#endif
