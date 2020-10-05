//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "../include/DynamicalVariables.h"

#define CONFORMAL_EOS

// ideal gas of massless quarks and gluons
#define EOS_FACTOR 13.8997 // Nc=3, Nf=2.5

#define EOS_ALPHA 0.2 //0.5 //


void getEquationOfStateTable();

void getPrimaryVariablesCombo(PRECISION e, PRECISION rhob, PRECISION * const __restrict__ PrimaryVariables);

PRECISION InferredPrimaryVariable(PRECISION e, PRECISION rhob, PRECISION e_start, PRECISION d_e, int nrhob, PRECISION d_rhob, int index_start, const PRECISION * const __restrict__ EOS_Variable);

// EoS with baryon

PRECISION equilibriumPressure(PRECISION e, PRECISION rhob);

PRECISION effectiveTemperature(PRECISION e, PRECISION rhob);

PRECISION chemicalPotentialOverT(PRECISION e, PRECISION rhob);

PRECISION speedOfSoundSquared(PRECISION e, PRECISION rhob);

PRECISION equilibriumEntropy(PRECISION e, PRECISION rhob, PRECISION p, PRECISION T, PRECISION alphaB);

PRECISION dPdRhob(PRECISION e, PRECISION rhob);

PRECISION dPdE(PRECISION e, PRECISION rhob);

PRECISION dPdT(PRECISION e, PRECISION rhob);

// Wuppertal-Budapest EoS without baryon

PRECISION dpdeWB(PRECISION e);

PRECISION equilibriumPressureWB(PRECISION e);

PRECISION effectiveTemperatureWB(PRECISION e);

PRECISION equilibriumEnergyDensityWB(PRECISION T);

#endif /* EQUATIONOFSTATE_H_ */
