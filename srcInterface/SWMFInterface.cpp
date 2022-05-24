#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "ReadParam.h"
#include "SWMFDomains.h"
#include "SWMFInterface.h"

Domains fleksDomains;

// The string contains the input string.
std::string paramString;
bool isFirstSession = true;
bool isInitialized = false;
double timeNow = -1;

//==========================================================
int pic_init_mpi_(MPI_Fint *iComm, signed int *iProc, signed int *nProc) {
  // fortran communicator tranlated to C comunicator
  MPI_Comm c_iComm = MPI_Comm_f2c(*iComm);

  amrex::Initialize(c_iComm);

  return 0;
}

//==========================================================
int pic_init_(double *time) {
  timeNow = (*time);
  return 0;
}

//==========================================================
int pic_read_param_(char *paramIn, int *nlines, int *ncharline, int *iProc) {

  paramString.clear();
  paramString = char_to_string(paramIn, (*nlines) * (*ncharline), *ncharline);

  if (!isFirstSession) {
    for (int i = 0; i < fleksDomains.size(); i++) {
      fleksDomains.select(i);
      fleksDomains(i).update_param(paramString);
    }
  }

  isFirstSession = false;
  return 0;
}

//==========================================================
int pic_from_gm_init_(int *paramint, double *paramreal, char *NameVar) {
  if (isInitialized)
    return 0;

  const int nDomain = paramint[1];
  for (int iDomain = 0; iDomain < nDomain; iDomain++)
    fleksDomains.add_new_domain();

  int nParamRegion = 22;
  for (int i = 0; i < fleksDomains.size(); i++) {
    fleksDomains.select(i);
    fleksDomains(i).init(timeNow, paramString, paramint,
                         &paramreal[i * nParamRegion],
                         &paramreal[nDomain * nParamRegion], i);
  }

  isInitialized = true;

  return 0;
}

//==========================================================
int pic_finalize_init_() {
  for (int i = 0; i < fleksDomains.size(); i++) {
    fleksDomains.select(i);
    fleksDomains(i).set_ic();
  }
  return 0;
}

//==========================================================
int pic_run_(double *time) {
  // For multiple FLEKS domains/grids, some domains may need to run a few steps
  // to make sure the time difference among the domains is within one step:
  // max(t_i)-min(t_i)<min(dt_i)

  bool mayNeedUpdata = true;
  double tMax = 0;

  int iLoop = 0;

  while (mayNeedUpdata) {
    mayNeedUpdata = false;

    for (int i = 0; i < fleksDomains.size(); i++) {
      double domainTime = (double)(fleksDomains(i).tc->get_time_si());
      double domainDt = (double)(fleksDomains(i).tc->get_dt_si());

      if (iLoop == 0 || domainTime + domainDt <= tMax + 1e-10 * domainDt) {
        // Update at least once, and a domain should be updated to a time that
        // is close to but smaller or equal to tMax.
        fleksDomains.select(i);
        fleksDomains(i).update();
        mayNeedUpdata = true;
      }

      if (iLoop == 0) {
        domainTime = (double)(fleksDomains(i).tc->get_time_si());
        // Find out tMax for the first while loop.
        if (i == 0) {
          tMax = domainTime;
        } else {
          if (tMax < domainTime)
            tMax = (domainTime);
        }
      }
    }
    iLoop++;
  }

  // Not all domains have reached tMax, but the difference should be within one
  // step.
  (*time) = tMax;
  return 0;
}

//==========================================================
int pic_save_restart_() {
  for (int i = 0; i < fleksDomains.size(); i++)
    fleksDomains(i).save_restart();
  return 0;
}

//==========================================================
int pic_get_ngridpoints_(int *nPoint) {
  *nPoint = 0;
  for (int i = 0; i < fleksDomains.size(); i++) {
    fleksDomains(i).couplerMarker = (*nPoint);
    *nPoint += fleksDomains(i).get_grid_nodes_number();
  }
  return 0;
}

//==========================================================
int pic_get_grid_(double *Pos_DI, int *n) {
  for (int i = 0; i < fleksDomains.size(); i++) {
    int idx = fleksDomains(i).couplerMarker *
              fleksDomains(i).fluidInterface->get_fluid_dimension();
    fleksDomains(i).get_grid(&Pos_DI[idx]);
  }
  return 0;
}

//==========================================================
int pic_set_state_var_(double *Data_VI, int *iPoint_I) {
  for (int i = 0; i < fleksDomains.size(); i++) {
    int idx = fleksDomains(i).couplerMarker;
    fleksDomains(i).set_state_var(Data_VI, &iPoint_I[idx]);
  }
  return 0;
}

//==========================================================
int pic_get_state_var_(int *nDim, int *nPoint, double *Xyz_I, double *data_I,
                       int *nVar) {
  for (int i = 0; i < fleksDomains.size(); i++)
    fleksDomains(i).get_fluid_state_for_points(*nDim, *nPoint, Xyz_I, data_I,
                                               *nVar);
  return 0;
}

//==========================================================
int pic_find_points_(int *nPoint, double *Xyz_I, int *iProc_I) {
  for (int i = 0; i < fleksDomains.size(); i++)
    fleksDomains(i).find_mpi_rank_for_points(*nPoint, Xyz_I, iProc_I);
  return 0;
}

//==========================================================
int pic_set_dt_(double *DtSi) {
  for (int i = 0; i < fleksDomains.size(); i++)
    fleksDomains(i).tc->set_dt_si(*DtSi);
  return 0;
}

//==========================================================
int pic_cal_dt_(double *dtOut) {
  for (int i = 0; i < fleksDomains.size(); i++)
    *dtOut = fleksDomains(i).tc->get_dt_si();
  return 0;
}

//==========================================================
int pic_get_grid_info_(int *iGrid, int *iDecomp) {
  (*iGrid) = 0;
  (*iDecomp) = 0;
  for (int i = 0; i < fleksDomains.size(); i++) {
    (*iGrid) += fleksDomains(i).get_iGrid();
    (*iDecomp) += fleksDomains(i).get_iDecomp();
  }
  return 0;
}

//==========================================================
int pic_end_() {
  {
    // Saving plots before exiting.
    for (int i = 0; i < fleksDomains.size(); i++) {
      fleksDomains.select(i);
      fleksDomains(i).write_plots(true);
    }

    fleksDomains.clear();

    // Q: Why the curly bracket is here?
    // A: It ensures the destructor is called and finished before the
    // amrex::Finalize().
  }

  amrex::Finalize();

  return 0;
}

//==========================================================
int pic_set_grid_info_(int *nInt, int *accumulatedSize, int *status) {
  for (int i = 0; i < fleksDomains.size(); i++) {
    int idxStart = 0;
    if (i > 0)
      idxStart = accumulatedSize[i - 1];
    fleksDomains(i).receive_grid_info(&status[idxStart]);
    fleksDomains(i).regrid();
  }

  return 0;
}
