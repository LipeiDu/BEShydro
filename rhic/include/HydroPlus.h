//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

// Equation indices from PRD 98 (2018) 036006

#ifndef HydroPlus_h
#define HydroPlus_h

#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroParameters.h"


PRECISION correlationLength(PRECISION T, PRECISION muB);

PRECISION corrLen(PRECISION T, PRECISION muB, PRECISION T_C, PRECISION MU_C);

PRECISION equilibriumPhiQ(PRECISION e, PRECISION rhob, PRECISION T, PRECISION muB, PRECISION s, PRECISION Q);

PRECISION relaxationCoefficientPhi(PRECISION rhob, PRECISION s, PRECISION T, PRECISION corrL2);

PRECISION relaxationCoefficientPhiQ(PRECISION gammaPhi, PRECISION corrL2, PRECISION Q);

PRECISION dlnXidrhob(PRECISION e0, PRECISION rhob0);

PRECISION dlnXide(PRECISION e0, PRECISION rhob0);

void getCorrelationLengthTable();

void setInitialConditionSlowModes(void * latticeParams, void * hydroParams);

void getPressurePlusFromSlowModes(PRECISION * const __restrict__ deltaVariables, PRECISION * const __restrict__ pPlus, const PRECISION * const __restrict__ equiPhiQ, const PRECISION * const __restrict__ PhiQ, PRECISION e, PRECISION rhob, PRECISION p, PRECISION T, PRECISION alphaB, PRECISION s);

#endif /* HydroPlus_h */
