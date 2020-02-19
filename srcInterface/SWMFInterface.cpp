#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "SWMFDomains.h"
#include "SWMFInterface.h"

#include "ReadParam.h"

using namespace std;

Domain *MPICs;

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

  MPICs = new Domain();

  int nParamRegion = 21;
  for (int i = 0; i < nDomain; i++) {
    MPICs->init(timenow, paramString, paramint, &paramreal[i * nParamRegion],
                &paramreal[nDomain * nParamRegion], i);
  }

  isInitialized = true;

  return 0;
}

int pic_finalize_init_() {
  MPICs->set_ic();
  return 0;
}

int pic_run_(double *time) {
  BL_PROFILE("pic_run");

  double timenow = *time;

  MPICs->update();

  *time = (double)(MPICs->tc->get_time_si());

  return 0;
}

int pic_save_restart_() {
  MPICs->save_restart();
  return 0;
}

int pic_get_ngridpoints_(int *nPoint) {
  *nPoint = MPICs->get_grid_nodes_number();
  return 0;
}

int pic_get_grid_(double *Pos_DI, int *n) {
  MPICs->get_grid(Pos_DI);
  return 0;
}

int pic_set_state_var_(double *Data_VI, int *iPoint_I) {
  MPICs->set_state_var(Data_VI, iPoint_I);
  return 0;
}

int pic_get_state_var_(int *nDim, int *nPoint, double *Xyz_I, double *data_I,
                       int *nVar) {
  MPICs->get_fluid_state_for_points(*nDim, *nPoint, Xyz_I, data_I, *nVar);
  return 0;
}

int pic_find_points_(int *nPoint, double *Xyz_I, int *iProc_I) {
  MPICs->find_mpi_rank_for_points(*nPoint, Xyz_I, iProc_I);
  return 0;
}

int pic_set_dt_(double *DtSi) {
  MPICs->tc->set_dt_si(*DtSi);
  return 0;
}

int pic_cal_dt_(double *dtOut) {
  *dtOut = MPICs->tc->get_dt_si();
  return 0;
}

int pic_get_grid_info_(int *iGrid, int *iDecomp) {
  (*iGrid) = MPICs->get_iGrid();
  (*iDecomp) = MPICs->get_iDecomp();
  return 0;
}

int pic_end_() {
  {
    // Saving plots before exiting.
    MPICs->write_plots(true);
    delete MPICs;

    //BL_PROFILE_VAR_STOP(pmain);
    // The curly bracket here is necessary!!! It ensures the destructor is
    // called and finished before the amrex::Finalize().
  }

  amrex::Finalize();

  return 0;
}


int pic_set_grid_info_(int *nInt, int *status){
  for(int i = 0; i<(*nInt); i++){
    std::cout<<"status = "<<status[i]<<std::endl;
  }
  MPICs->receive_grid_info(status);
  MPICs->regrid();

  return 0; 
}
