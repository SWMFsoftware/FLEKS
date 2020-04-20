#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "SWMFDomains.h"
#include "SWMFInterface.h"

#include "ReadParam.h"

std::vector<Domain *> FLEKSs;
// Domain *FLEKSs;

// // The string contains the param.in
std::string paramString;

bool isFirstSession = true;
bool isInitialized = false;

double timenow = 0;

int pic_init_mpi_(MPI_Fint *iComm, signed int *iProc, signed int *nProc) {
  // fortran communicator tranlated to C comunicator
  MPI_Comm c_iComm = MPI_Comm_f2c(*iComm);

  amrex::Initialize(c_iComm);

  return 0;
}

int pic_init_(double *inittime) {
  timenow = *inittime;
  return 0;
}

int pic_read_param_(char *paramIn, int *nlines, int *ncharline, int *iProc) {
  if (!isFirstSession)
    return 0;

  paramString = char_to_string(paramIn, (*nlines) * (*ncharline), *ncharline);

  isFirstSession = false;
  return 0;
}

int pic_from_gm_init_(int *paramint, double *paramreal, char *NameVar) {
  if (isInitialized)
    return 0;

  // // number of dimensions in GM
  int nDim = paramint[0];

  const int nDomain = paramint[1];

  Domain *tmp = new Domain();
  FLEKSs.push_back(tmp);

  int nParamRegion = 21;
  for (int i = 0; i < nDomain; i++) {
    FLEKSs[0]->init(timenow, paramString, paramint,
                    &paramreal[i * nParamRegion],
                    &paramreal[nDomain * nParamRegion], i);
  }

  isInitialized = true;

  return 0;
}

int pic_finalize_init_() {
  FLEKSs[0]->set_ic();
  return 0;
}

int pic_run_(double *time) {
  double timenow = *time;

  FLEKSs[0]->update();

  *time = (double)(FLEKSs[0]->tc->get_time_si());

  return 0;
}

int pic_save_restart_() {
  FLEKSs[0]->save_restart();
  return 0;
}

int pic_get_ngridpoints_(int *nPoint) {
  *nPoint = FLEKSs[0]->get_grid_nodes_number();
  return 0;
}

int pic_get_grid_(double *Pos_DI, int *n) {
  FLEKSs[0]->get_grid(Pos_DI);
  return 0;
}

int pic_set_state_var_(double *Data_VI, int *iPoint_I) {
  FLEKSs[0]->set_state_var(Data_VI, iPoint_I);
  return 0;
}

int pic_get_state_var_(int *nDim, int *nPoint, double *Xyz_I, double *data_I,
                       int *nVar) {
  FLEKSs[0]->get_fluid_state_for_points(*nDim, *nPoint, Xyz_I, data_I, *nVar);
  return 0;
}

int pic_find_points_(int *nPoint, double *Xyz_I, int *iProc_I) {
  FLEKSs[0]->find_mpi_rank_for_points(*nPoint, Xyz_I, iProc_I);
  return 0;
}

int pic_set_dt_(double *DtSi) {
  FLEKSs[0]->tc->set_dt_si(*DtSi);
  return 0;
}

int pic_cal_dt_(double *dtOut) {
  *dtOut = FLEKSs[0]->tc->get_dt_si();
  return 0;
}

int pic_get_grid_info_(int *iGrid, int *iDecomp) {
  (*iGrid) = FLEKSs[0]->get_iGrid();
  (*iDecomp) = FLEKSs[0]->get_iDecomp();
  return 0;
}

int pic_end_() {
  {
    // Saving plots before exiting.
    FLEKSs[0]->write_plots(true);

    for (auto& FLEKS : FLEKSs) {
      delete FLEKS;
    }
    FLEKSs.clear();

    // BL_PROFILE_VAR_STOP(pmain);
    // The curly bracket here is necessary!!! It ensures the destructor is
    // called and finished before the amrex::Finalize().
  }

  amrex::Finalize();

  return 0;
}

int pic_set_grid_info_(int *nInt, int *status) {
  FLEKSs[0]->receive_grid_info(status);
  FLEKSs[0]->regrid();

  return 0;
}
