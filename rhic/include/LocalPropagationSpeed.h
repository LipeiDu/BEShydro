//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef LOCALPROPAGATIONSPEED_H_
#define LOCALPROPAGATIONSPEED_H_

#include "../include/DynamicalVariables.h"

PRECISION localPropagationSpeed(PRECISION utr, PRECISION uxr, PRECISION uyr, PRECISION unr, PRECISION utl, PRECISION uxl, PRECISION uyl, PRECISION unl, PRECISION (*spectralRadius)(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un));

#endif /* LOCALPROPAGATIONSPEED_H_ */
