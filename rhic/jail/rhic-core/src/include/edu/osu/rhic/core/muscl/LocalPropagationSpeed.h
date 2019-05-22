/*
 * LocalPropagationSpeed.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef LOCALPROPAGATIONSPEED_H_
#define LOCALPROPAGATIONSPEED_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

PRECISION localPropagationSpeed(PRECISION utr, PRECISION uxr, PRECISION uyr, PRECISION unr,
		PRECISION utl, PRECISION uxl, PRECISION uyl, PRECISION unl,
		PRECISION (*spectralRadius)(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un)
);

#endif /* LOCALPROPAGATIONSPEED_H_ */
