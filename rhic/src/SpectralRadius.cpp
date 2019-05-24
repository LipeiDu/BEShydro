//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <math.h>

#include "../include/SpectralRadius.h"
#include "../include/PrimaryVariables.h"
#include "../include/DynamicalVariables.h"

PRECISION spectralRadiusX(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return fabs(ux/ut);
}

PRECISION spectralRadiusY(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return fabs(uy/ut);
}

PRECISION spectralRadiusZ(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return fabs(un/ut);
}
