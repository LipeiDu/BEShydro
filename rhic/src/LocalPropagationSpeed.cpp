/*
 * LocalPropagationSpeed.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <stdlib.h>
#include <math.h>

#include "../include/LocalPropagationSpeed.h"
#include "../include/DynamicalVariables.h"

// maximal local speed at the cell boundaries x_{j\pm 1/2}
PRECISION localPropagationSpeed(PRECISION utr, PRECISION uxr, PRECISION uyr, PRECISION unr, PRECISION utl, PRECISION uxl, PRECISION uyl, PRECISION unl, PRECISION (*spectralRadius)(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un))
{
	PRECISION rhoLeftMovingWave = spectralRadius(utl,uxl,uyl,unl);
	PRECISION rhoRightMovingWave = spectralRadius(utr,uxr,uyr,unr);
	PRECISION a = fmax(rhoLeftMovingWave, rhoRightMovingWave);
	return a;
}
