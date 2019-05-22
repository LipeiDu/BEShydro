/*
 * SpectralRadius.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include <math.h>

#include "edu/osu/rhic/trunk/hydro/SpectralRadius.h"
#include "edu/osu/rhic/trunk/hydro/EnergyMomentumTensor.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

PRECISION spectralRadiusX(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return fabs(ux/ut);
}

PRECISION spectralRadiusY(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return fabs(uy/ut);
}

PRECISION spectralRadiusZ(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return fabs(un/ut);
}
