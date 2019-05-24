//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef SPECTRALRADIUS_H_
#define SPECTRALRADIUS_H_

#include "../include/DynamicalVariables.h"

PRECISION spectralRadiusX(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);
PRECISION spectralRadiusY(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);
PRECISION spectralRadiusZ(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);

#endif /* SPECTRALRADIUS_H_ */
