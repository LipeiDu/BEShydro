//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef FLUXLIMITER_H_
#define FLUXLIMITER_H_

#include "../include/DynamicalVariables.h"

#define THETA 1.8

PRECISION approximateDerivative(PRECISION x, PRECISION y, PRECISION z);

#endif /* FLUXLIMITER_H_ */
