#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "DomainGrid.h"
#include "Pic.h"

class Domain : public DomainGrid {
private:
  bool doRestart;

public:
  Pic pic;

  std::shared_ptr<FluidInterface> fluidInterface;

  std::shared_ptr<TimeCtr> tc;

public:
  Domain() {
    fluidInterface = std::make_shared<FluidInterface>();
    tc = std::make_unique<TimeCtr>();
    doRestart = false;
  }

  ~Domain() = default;

  void update();

  //--------------Initialization begin-------------------------------
  void init(amrex::Real timeIn, const std::string &paramString, int *paramInt,
            double *gridDim, double *paramReal, int iDomain = 1);
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
  void read_param();
  void save_restart();
  void save_restart_header();
  void save_restart_data();
  void read_restart();
  void write_plots(bool doForce = false);
  //--------------- IO end----------------------------------

  void make_grid();
  void make_data();
  void init_time_ctr();
};

#endif