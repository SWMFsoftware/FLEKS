#ifndef _SWMFINTERFACE_H_
#define _SWMFINTERFACE_H_

#include <mpi.h>
#include <sstream>

extern "C" {
void OH_get_charge_exchange_wrapper(double *rhoIon, double *cs2Ion,
                                    double uIon_D[3], double *rhoNeu,
                                    double *cs2Neu, double uNeu_D[3],
                                    double sourceIon_V[5],
                                    double sourceNeu_V[5]);
}


//   {
//     double rhoIon, cs2Ion, uIon_D[3], rhoNeu, cs2Neu, uNeu_D[3], sourceIon_V[5],
//         sourceNeu_V[5];
//     OH_get_charge_exchange_wrapper(&rhoIon, &cs2Ion, uIon_D, &rhoNeu, &cs2Neu,
//                                    uNeu_D, sourceIon_V, sourceNeu_V);
//   }


#endif
