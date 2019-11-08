#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include <AMReX.H>
#include <AMReX_Print.H>

#include "ipic3d_interface.h"
#include "multi_ipic3d_domain.h"

// Files in SWMF/share/Library/src/
//#include "Timing_c.h"
#include "ReadParam.h"

// using namespace iPic3D;
using namespace std;

Domain *MPICs;

// int nIPIC;
// int iIPIC;
// int *iSimCycle;

// // pointer to all IPIC3D
// iPic3D::c_Solver **SimRun;

// // The string contains the param.in
std::string paramString;

// // first time for region
// bool firstcall[nmaxIPIC];
bool isFirstSession = true;
bool isInitialized = false;

// // store start time
// double starttime[nmaxIPIC];
double timenow = 0;

// // Store variables asosiated with comunicating the grid to and from GM
// // When we recive data form BATSRUS we need to split it up and give it to
// // the respective SimRun object. This to arrays will be used for this purpus.
// int nGridPntSim[nmaxIPIC]; // number of grid pints *ndim used for each sim
// box int nShiftGridPntSim[nmaxIPIC]; // fist index for each sims grid ponts

int ipic3d_init_mpi_(MPI_Fint *iComm, signed int *iProc, signed int *nProc) {
  // fortran communicator tranlated to C comunicator
  MPI_Comm c_iComm = MPI_Comm_f2c(*iComm);

  amrex::Initialize(c_iComm);
  amrex::Print() << "init_amrex is called!" << std::endl;

  // At this time we do not have a proper timstep
  // so we only do steping internaly iSimCycle

  // if (MPI_SUCCESS != MPI_Comm_dup(c_iComm, &ipic3dComm)) {
  //   cout << "IPIC3D_init_mpi :: Can not copy  MPI comunicator, arborting!"
  //        << endl;
  //   cout.flush();
  //   abort();
  // }

  // myrank = *iProc;
  // numProc = *nProc;
  // nIPIC = 0;
  // SimRun = new iPic3D::c_Solver *[nmaxIPIC];
  // iSimCycle = new int[nmaxIPIC];

  // isInitialized = false;
  // for (int i = 0; i < nmaxIPIC; i++) {
  //   firstcall[i] = true;
  //   starttime[i] = 0.0;
  //   SimRun[i] = NULL;
  // }

  // myProc = *iProc;

  // MPIdata::init(ipic3dComm, &myrank, &numProc);

  return (0);
}

int ipic3d_init_(double *inittime) {
  //  Called by  PC_init_session
  //  Not used as this time as
  //  most of the grid are set up when
  //  we read PARAM.in, and the rest
  //  when IPIC3D recive data for the
  //  first time.
  timenow = *inittime;
  return (0);
}

int ipic3d_read_param_(char *paramIn, int *nlines, int *ncharline, int *iProc) {
  // So far, iPIC3D does not support multi sessions.
  if (!isFirstSession)
    return 0;
  std::string nameFunc = "ipic3d_read_param";
  amrex::Print() << nameFunc << std::endl;
  // // convert character array to string stream object
  // myProc = *iProc;

  paramString = char_to_string(paramIn, (*nlines) * (*ncharline), *ncharline);

  // std::cout<<"paramString = "<<paramString<<std::endl;

  isFirstSession = false;
  return (0);
}

int ipic3d_from_gm_init_(int *paramint, double *paramreal, char *NameVar) {
  std::cout << "isInitialized = " << (isInitialized ? 'T' : 'F') << std::endl;
  if (isInitialized)
    return 0;
  std::stringstream *ss;
  ss = new std::stringstream;
  (*ss) << NameVar;

  // int firstIPIC;
  std::string nameFunc = "ipic3d_read_param";
  amrex::Print() << nameFunc << std::endl;

  // timing_start(nameFunc);

  // // number of dimensions in GM
  int nDim = paramint[0];

  // firstIPIC = 0;
  const int nDomain = paramint[1];
  // for (int i = 0; i < nDomain; i++) {
  //   Domain* ptr = new Domain();
  //   MPICs->push_back(ptr);
  // }
  MPICs = new Domain();

  int nParamRegion = 21;
  for (int i = 0; i < nDomain; i++) {
    MPICs->init(timenow, paramString, paramint, &paramreal[i * nParamRegion],
                &paramreal[nDomain * nParamRegion], i);
  }

  // char **dummy = NULL; // dummy argument

  // // the number of PIC fluids ns = nFluid + nSpecies - 1
  // int ns = paramint[3] + paramint[4] - 1;
  // for (int i = firstIPIC; i < nIPIC; i++) {
  //   starttime[i] = timenow;
  //   SimRun[i] = new iPic3D::c_Solver;
  //   int nParamRegion = 21; // See GM_couple_pc.f90
  //   SimRun[i]->Init(0, dummy, timenow, paramString, i, paramint,
  //                   &paramreal[(i - firstIPIC) * nParamRegion],
  //                   &paramreal[(nIPIC - firstIPIC) * nParamRegion], ss,
  //                   true);
  //   SimRun[i]->SetCycle(0);
  // }
  // timing_stop(nameFunc);
  isInitialized = true;

  delete ss;
  return (0);
}

int ipic3d_finalize_init_() {
  std::string nameFunc = "ipic3d_finalize_init";
  amrex::Print() << nameFunc << std::endl;

  MPICs->set_ic();

  // This function is called by the coupler the first time
  // it want to couple from GM -> PC. At this point we should have all
  // information to finnish the initialization of SimRun[i]

  // string nameFunc = "PC: ipic3d_finalize_init_";
  // timing_start(nameFunc);

  // string tmp;

  // char **dummy = NULL; // dummy argument

  // // now we should have all the the infomation
  // for (int i = 0; i < nIPIC; i++) {
  //   SimRun[i]->Init(0, dummy, timenow, tmp);
  //   SimRun[i]->CalculateMoments();
  //   iSimCycle[i] = SimRun[i]->FirstCycle();
  //   SimRun[i]->WriteOutput(iSimCycle[i], true);
  // }
  // timing_stop(nameFunc);
  return (0);
}

int ipic3d_run_(double *time) {
  // bool b_err = false;

  std::string nameFunc = "ipic3d_run";
  amrex::Print() << nameFunc << std::endl;

  double timenow = *time;
  // string nameFunc = "PC: IPIC3D_run";
  // timing_start(nameFunc);

  // for (int i = 0; i < nIPIC; i++) {
  //   iIPIC = i;
  //   if (SimRun[i]->get_myrank() == 0)
  //     cout << " ======= Cycle " << iSimCycle[i] << ", dt=" <<
  //     SimRun[i]->getDt()
  //          << " ======= " << endl;

  //   SimRun[i]->SetCycle(iSimCycle[i]);

  //   timing_start("PC: SyncWithFluid");
  //   SimRun[i]->SyncWithFluid(iSimCycle[i]);
  //   timing_stop("PC: SyncWithFluid");

  //   if (iSimCycle[i] == 0) {
  //     // iSimCycle[i]==0 is true for the first iteration when the PIC code
  //     // stars from scratch (NOT from restart files.). When CalculateMoments
  //     // is called during the finalization stage, dt is still unknow, and
  //     // the mass matrix, which is used for ECSIM, is not correct. So
  //     // CalculateMoments should be called again to calculate mass matrix
  //     // for this case.
  //     timing_start("PC: GatherMoments");
  //     SimRun[i]->CalculateMoments();
  //     timing_stop("PC: GatherMoments");
  //   }

  //   timing_start("PC: CalculateField");
  //   SimRun[i]->CalculateField(iSimCycle[i]);
  //   timing_stop("PC: CalculateField");

  //   timing_start("PC: ParticlesMover");
  //   b_err = SimRun[i]->ParticlesMover();
  //   timing_stop("PC: ParticlesMover");

  //   timing_start("PC: CalculateBField");
  //   if (!b_err)
  //     SimRun[i]->CalculateB();
  //   timing_stop("PC: CalculateBField");

  //   timing_start("PC: GatherMoments");
  //   if (!b_err)
  //     SimRun[i]->CalculateMoments(true);
  //   timing_stop("PC: GatherMoments");

  //   if (b_err) {
  //     iSimCycle[i] = SimRun[i]->LastCycle() + 1;
  //   }

  //   SimRun[i]->updateSItime();

  //   timing_start("PC: WriteOutput");
  //   SimRun[i]->WriteOutput(iSimCycle[i] + 1);
  //   timing_stop("PC: WriteOutput");

  //   iSimCycle[i]++;
  // }

  // // All simulations are in the same time zone
  MPICs->update();
  *time = (double)(MPICs->tc.get_time_si());
  // timing_stop(nameFunc);
  return (0);
}

int ipic3d_save_restart_() {
  // string nameFunc = "PC: ipic3d_save_restart";
  // timing_start(nameFunc);

  // for (int i = 0; i < nIPIC; i++)
  //   SimRun[i]->WriteRestart(iSimCycle[i] - 1);

  // timing_stop(nameFunc);
  return (0);
}

int ipic3d_get_ngridpoints_(int *nPoint) {

  // *nPoint = 0;
  // for (int i = 0; i < nIPIC; i++) {
  //   SimRun[i]->GetNgridPnt(nGridPntSim[i]);
  //   nShiftGridPntSim[i] = *nPoint;
  //   *nPoint += nGridPntSim[i];
  // }
  // // Fortran operates on an 2D array [ndim,nPoint]
  *nPoint = MPICs->get_grid_nodes_number();
  std::cout << " nPoint = " << (*nPoint) << std::endl;
  return (0);
}

int ipic3d_get_grid_(double *Pos_DI, int *n) {
  // for (int i = 0; i < nIPIC; i++)
  //   SimRun[i]->GetGridPnt(&Pos_DI[nShiftGridPntSim[i]]);
  MPICs->get_grid(Pos_DI);

  return (0);
}

int ipic3d_set_state_var_(double *Data_VI, int *iPoint_I) {
  std::string nameFunc = "PC: ipic3d_set_state_var_";
  amrex::Print() << nameFunc << std::endl;

  // timing_start(nameFunc);

  // // WARNING  iPoint_I is a reindexing of state var array.

  // for (int i = 0; i < nIPIC; i++) {
  //   SimRun[i]->setStateVar(Data_VI, &iPoint_I[(nShiftGridPntSim[i] / nDim)]);
  // }
  // timing_stop(nameFunc);
  MPICs->set_state_var(Data_VI, iPoint_I);
  return (0);
}

int ipic3d_get_state_var_(int *nDim, int *nPoint, double *Xyz_I, double *data_I,
                          int *nVar) {
  std::string nameFunc = "PC: IPIC3D_get_state_var";
  amrex::Print() << nameFunc << std::endl;

  for (int i = 0; i < (*nPoint) * (*nVar); i++) {
    data_I[i] = -1;
  }

  // amrex::Abort();
  // timing_start(nameFunc);

  // for (int i = 0; i < nIPIC; i++) {
  //   SimRun[i]->getStateVar(*nDim, *nPoint, Xyz_I, data_I, *nVar);
  // }
  // timing_stop(nameFunc);
  return (0);
}

int ipic3d_find_points_(int *nPoint, double *Xyz_I, int *iProc_I) {
  std::string nameFunc = "PC: ipic3d_find_points_";
  amrex::Print() << nameFunc << " begin" << std::endl;

  for (int i = 0; i < (*nPoint); i++) {
    iProc_I[i] = 0;
  }

  // amrex::Abort();
  // timing_start(nameFunc);

  // for (int i = 0; i < nIPIC; i++) {
  //   SimRun[i]->findProcForPoint(*nPoint, Xyz_I, iProc_I);
  // }
  // timing_stop(nameFunc);
  return (0);
}

int ipic3d_set_dt_(double *DtSi) {
  //  for (int i = 0; i < nIPIC; i++)
  std::string nameFunc = "PC: IPIC3D_set_dt";
  amrex::Print() << nameFunc << std::endl;

  MPICs->tc.set_dt_si(*DtSi);

  return (0);
}

int ipic3d_cal_dt_(double *dtOut) {
  // Each PIC domain should be treated seperately in the future. Now, all the
  // domains use the same time step.

  // double dt = 1e10, dt0;

  // for (int i = 0; i < nIPIC; i++) {
  //   dt0 = SimRun[i]->calSIDt();
  //   if (dt0 < dt)
  //     dt = dt0;
  // }

  std::string nameFunc = "PC: IPIC3D_cal_dt";
  amrex::Print() << nameFunc << std::endl;

  *dtOut = MPICs->tc.get_dt_si();
  return (0);
}

int ipic3d_end_() {  
  { // Saving plots before exiting.
    MPICs->tc.write_plots(true);
    delete MPICs;
  }

  amrex::Print() << "finalize_amrex is called!" << std::endl;
  amrex::Finalize();

  return 0;
}
