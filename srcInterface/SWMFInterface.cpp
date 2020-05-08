#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "SWMFDomains.h"
#include "SWMFInterface.h"

#include "ReadParam.h"

Domains FLEKSs;

// // The string contains the param.in
std::string paramString;

bool isFirstSession = true;
bool isInitialized = false;

int pic_init_mpi_(MPI_Fint *iComm, signed int *iProc, signed int *nProc) {
  // fortran communicator tranlated to C comunicator
  MPI_Comm c_iComm = MPI_Comm_f2c(*iComm);

  amrex::Initialize(c_iComm);

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
  for (int iDomain = 0; iDomain < nDomain; iDomain++)
    FLEKSs.add_new_domain();

  int nParamRegion = 21;
  for (int i = 0; i < FLEKSs.size(); i++) {
    FLEKSs(i).init(paramString, paramint, &paramreal[i * nParamRegion],
                   &paramreal[nDomain * nParamRegion], i);
  }

  isInitialized = true;

  return 0;
}

int pic_finalize_init_() {
  for (int i = 0; i < FLEKSs.size(); i++){
    FLEKSs.select(i); 
    FLEKSs(i).set_ic();
  }
  return 0;
}

int pic_run_(double *time) {

  for (int i = 0; i < FLEKSs.size(); i++) {
    FLEKSs.select(i);
    FLEKSs(i).update();
    *time = (double)(FLEKSs(i).tc->get_time_si());
  }

  return 0;
}

int pic_save_restart_() {
  for (int i = 0; i < FLEKSs.size(); i++)
    FLEKSs(i).save_restart();
  return 0;
}

int pic_get_ngridpoints_(int *nPoint) {
  *nPoint = 0; 
  for (int i = 0; i < FLEKSs.size(); i++){
    FLEKSs(i).couplerMarker = (*nPoint); 
    *nPoint += FLEKSs(i).get_grid_nodes_number();
  }
  return 0;
}

int pic_get_grid_(double *Pos_DI, int *n) {
  for (int i = 0; i < FLEKSs.size(); i++){
    int idx = FLEKSs(i).couplerMarker*FLEKSs(i).fluidInterface->getnDim();
    FLEKSs(i).get_grid(&Pos_DI[idx]);
  }
  return 0;
}

int pic_set_state_var_(double *Data_VI, int *iPoint_I) {
  for (int i = 0; i < FLEKSs.size(); i++){
    int idx = FLEKSs(i).couplerMarker;
    FLEKSs(i).set_state_var(Data_VI, &iPoint_I[idx]);
  }
  return 0;
}

int pic_get_state_var_(int *nDim, int *nPoint, double *Xyz_I, double *data_I,
                       int *nVar) {
  for (int i = 0; i < FLEKSs.size(); i++)
    FLEKSs(i).get_fluid_state_for_points(*nDim, *nPoint, Xyz_I, data_I, *nVar);
  return 0;
}

int pic_find_points_(int *nPoint, double *Xyz_I, int *iProc_I) {
  for (int i = 0; i < FLEKSs.size(); i++)
    FLEKSs(i).find_mpi_rank_for_points(*nPoint, Xyz_I, iProc_I);
  return 0;
}

int pic_set_dt_(double *DtSi) {
  for (int i = 0; i < FLEKSs.size(); i++)
    FLEKSs(i).tc->set_dt_si(*DtSi);
  return 0;
}

int pic_cal_dt_(double *dtOut) {
  for (int i = 0; i < FLEKSs.size(); i++)
    *dtOut = FLEKSs(i).tc->get_dt_si();
  return 0;
}

int pic_get_grid_info_(int *iGrid, int *iDecomp) {
  (*iGrid) = 0;
  (*iDecomp) = 0;
  for (int i = 0; i < FLEKSs.size(); i++) {
    (*iGrid) += FLEKSs(i).get_iGrid();
    (*iDecomp) += FLEKSs(i).get_iDecomp();
  }
  return 0;
}

int pic_end_() {
  {
    // Saving plots before exiting.
    for (int i = 0; i < FLEKSs.size(); i++){      
      FLEKSs.select(i); 
      FLEKSs(i).write_plots(true);
    }

    FLEKSs.clear();

    // BL_PROFILE_VAR_STOP(pmain);
    // The curly bracket here is necessary!!! It ensures the destructor is
    // called and finished before the amrex::Finalize().
  }

  amrex::Finalize();

  return 0;
}

int pic_set_grid_info_(int *nInt, int *accumulatedSize, int *status) {
  for (int i = 0; i < FLEKSs.size(); i++) {
    int idxStart = 0;
    if (i > 0)
      idxStart = accumulatedSize[i-1];
    printf("pic_set_grid_info i = %d\n",i);
    FLEKSs(i).receive_grid_info(&status[idxStart]);
    FLEKSs(i).regrid();
  }

  return 0;
}
