//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef TransportCoefficients_h
#define TransportCoefficients_h

#include <stdio.h>

//shear transport coefficients

const PRECISION delta_pipi = 1.33333;
const PRECISION tau_pipi = 1.42857;
const PRECISION delta_PiPi = 0.666667;
const PRECISION lambda_piPi = 1.2;
const PRECISION tau_piw = 1.0; //coupling between shear and vorticity

//bulk transport coefficients

PRECISION bulkViscosityToEntropyDensity(PRECISION T);

// baryon transport coefficients

const PRECISION Cb = 0.4;//4.0;
const PRECISION delta_nn = 1.0;//0.;//
const PRECISION lambda_nn = 0.6;//0.;//
const PRECISION tau_nw = 1.0;//0.;// coupling between baryon diffusion and vorticity

// functions

void getBaryonDiffusionCoefficientTable();

void baryonDiffusionCoefficient(PRECISION T, PRECISION muB, PRECISION * const __restrict__ diffusionCoeff);

PRECISION baryonDiffusionCoefficientKinetic(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p);

PRECISION baryonDiffusionCoefficientTest(PRECISION T, PRECISION rhob, PRECISION alphaB);

PRECISION baryonDiffusionCoefficientAdscft(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p, PRECISION seq);

PRECISION baryonDiffusionCoefficientHydroPlus(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p, PRECISION seq);


#endif /* TransportCoefficients_hpp */
