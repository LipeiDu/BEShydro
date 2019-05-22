/*
 * LocalPropagationSpeed.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <stdlib.h>
#include <math.h>

#include "edu/osu/rhic/core/muscl/LocalPropagationSpeed.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

// maximal local speed at the cell boundaries x_{j\pm 1/2}
PRECISION localPropagationSpeed(PRECISION utr, PRECISION uxr, PRECISION uyr, PRECISION unr,
		PRECISION utl, PRECISION uxl, PRECISION uyl, PRECISION unl,
		PRECISION (*spectralRadius)(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un)
) {
	PRECISION rhoLeftMovingWave = spectralRadius(utl,uxl,uyl,unl);
	PRECISION rhoRightMovingWave = spectralRadius(utr,uxr,uyr,unr);
	PRECISION a = fmax(rhoLeftMovingWave, rhoRightMovingWave);
	return a;
}
