#include "Domain.h"

//------------------------------------------------------------------------
void Domain::update() {
    pic.update();
  };

//------------------------------------------------------------------------
void Domain::init(amrex::Real timeIn, const std::string &paramString,
                  int *paramInt, double *gridDim, double *paramReal,
                  int iDomain) {

  pic.init(timeIn, paramString, paramInt, gridDim, paramReal, fluidInterface, tc);
};

//------------------------------------------------------------------------
void Domain::set_ic() {
  pic.set_ic();
};

//------------------------------------------------------------------------
void Domain::set_state_var(double *data, int *index) {
  pic.set_state_var(data, index);
};

//------------------------------------------------------------------------
int Domain::get_grid_nodes_number() {
  return pic.get_grid_nodes_number();
};

//------------------------------------------------------------------------
void Domain::get_grid(double *pos_DI) {
  pic.get_grid(pos_DI);
};

//------------------------------------------------------------------------
void Domain::find_mpi_rank_for_points(const int nPoint,
                                      const double *const xyz_I,
                                      int *const rank_I) {
  pic.find_mpi_rank_for_points(nPoint, xyz_I, rank_I);
};

//------------------------------------------------------------------------
void Domain::get_fluid_state_for_points(const int nDim, const int nPoint,
                                        const double *const xyz_I,
                                        double *const data_I, const int nVar) {
  pic.get_fluid_state_for_points(nDim, nPoint, xyz_I, data_I, nVar);
};

//------------------------------------------------------------------------
void Domain::save_restart() {
  pic.save_restart();
};
