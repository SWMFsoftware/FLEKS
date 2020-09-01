#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "DomainGrid.h"
#include "Pic.h"
#include "ParticleTracker.h"

class Domain : public DomainGrid {
private:
  bool doRestart;

  bool usePT; 
public:
  Pic pic;
  ParticleTracker pt; 

  //  Conceptually, both the Domain class and the Pic class may use the
  //  following classes, so they are handled by shared ponters.
  std::shared_ptr<FluidInterface> fluidInterface;
  std::shared_ptr<TimeCtr> tc;

  int couplerMarker;   

public:
  Domain() {
    fluidInterface = std::make_shared<FluidInterface>();
    tc = std::make_shared<TimeCtr>();
    doRestart = false;
    usePT = false; 
  }

  ~Domain() = default;

  void update();

  //--------------Initialization begin-------------------------------
  void init(double time, const std::string &paramString, int *paramInt, double *gridDim,
            double *paramReal, int iDomain = 1);
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

  //-------------- grid begin-------------------------------
  void make_grid();
  void regrid();
  void receive_grid_info(int *status);
  //-------------- grid end---------------------------------

  // void make_data();
  void init_time_ctr();
};

#endif
