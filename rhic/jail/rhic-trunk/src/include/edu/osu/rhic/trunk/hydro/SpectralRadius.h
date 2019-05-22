/*
 * SpectralRadius.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef SPECTRALRADIUS_H_
#define SPECTRALRADIUS_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

PRECISION spectralRadiusX(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);
PRECISION spectralRadiusY(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);
PRECISION spectralRadiusZ(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un);

#endif /* SPECTRALRADIUS_H_ */
