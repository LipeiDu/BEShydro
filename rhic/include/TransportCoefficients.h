//
//  TransportCoefficients.hpp
//  
//
//  Created by Lipei Du on 3/10/19.
//

#ifndef TransportCoefficients_h
#define TransportCoefficients_h

#include <stdio.h>

//Transport coefficients

const PRECISION delta_pipi = 1.33333;
const PRECISION tau_pipi = 0.0;//1.42857;
const PRECISION delta_PiPi = 0.666667;
const PRECISION lambda_piPi = 1.2;
const PRECISION tau_piw = 0.0;//1.0; coupling between shear and vorticity

PRECISION bulkViscosityToEntropyDensity(PRECISION T);

// baryon diffusion coefficients

const PRECISION Cb = 4.0;

void getBaryonDiffusionCoefficientTable();

void baryonDiffusionCoefficient(PRECISION T, PRECISION muB, PRECISION * const __restrict__ diffusionCoeff);

PRECISION baryonDiffusionCoefficientKinetic(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p);

PRECISION baryonDiffusionCoefficientAdscft(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p, PRECISION seq);

PRECISION baryonDiffusionCoefficientHydroPlus(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p, PRECISION seq);


#endif /* TransportCoefficients_hpp */
