/*
 * EquationOfState.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

#define CONFORMAL_EOS

// ideal gas of massless quarks and gluons
//#define EOS_FACTOR 15.6269 // Nc=3, Nf=3
#define EOS_FACTOR 13.8997 // Nc=3, Nf=2.5


void getEquationOfStateTable();//Lipei
void testEOS();

PRECISION baryonDiffusionConstant(PRECISION T, PRECISION mub);

PRECISION equilibriumPressure(PRECISION e, PRECISION rhob);
PRECISION equilibriumPressure(PRECISION e);

PRECISION speedOfSoundSquared(PRECISION e, PRECISION rhob);

PRECISION effectiveTemperature(PRECISION e, PRECISION rhob);
PRECISION effectiveTemperature(PRECISION e);

PRECISION chemicalPotentialOverT(PRECISION e, PRECISION rhob);

PRECISION dPdRhob(PRECISION e, PRECISION rhob);

PRECISION equilibriumEnergyDensity(PRECISION T);

#endif /* EQUATIONOFSTATE_H_ */
