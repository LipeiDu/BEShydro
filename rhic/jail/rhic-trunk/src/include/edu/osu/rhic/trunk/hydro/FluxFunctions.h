/*
 * FluxFunctions.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef FLUXFUNCTIONS_H_
#define FLUXFUNCTIONS_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

PRECISION Fx(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);
PRECISION Fy(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);
PRECISION Fz(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);

#endif /* FLUXFUNCTIONS_H_ */
