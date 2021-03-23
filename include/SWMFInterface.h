#ifndef _SWMFINTERFACE_H_
#define _SWMFINTERFACE_H_

#include <mpi.h>
#include <sstream>

extern "C" {
int pic_init_mpi_(MPI_Fint *iComm, signed int *iProc, signed int *nProc);
int pic_init_(double *timeNow);
int pic_finalize_init_();
int pic_run_(double *time);
int pic_save_restart_();
int pic_end_();
int pic_read_param_(char *param, int *nlines, int *ncharline, int *iProc);
int pic_from_gm_init_(int *paramint, double *paramreal, char *NameVar);
int pic_get_ngridpoints_(int *nPoint);
int pic_get_grid_(double *Pos_DI, int *n);
int pic_set_state_var_(double *Data_VI, int *iPoint_I);
int pic_cal_dt_(double *dt);
int pic_find_points_(int *nPoint, double *Xyz_DI, int *iProc_I);
int pic_get_state_var_(int *nDim, int *nPoint, double *Xyz_I, double *data_I,
                       int *nVar);
int pic_get_grid_info_(int *iGrid, int *iDecomp);
int pic_set_grid_info_(int *nInt, int *accumulatedSize, int *status);
}
#endif
