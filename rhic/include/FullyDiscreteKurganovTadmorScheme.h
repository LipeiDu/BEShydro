//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef FULLYDISCRETEKURGANOVTADMORSCHEME_H_
#define FULLYDISCRETEKURGANOVTADMORSCHEME_H_

#include "../include/DynamicalVariables.h"

void rungeKutta2(PRECISION t, PRECISION dt, CONSERVED_VARIABLES * __restrict__ q, CONSERVED_VARIABLES * __restrict__ Q, void * latticeParams, void * hydroParams);

#endif /* FULLYDISCRETEKURGANOVTADMORSCHEME_H_ */
