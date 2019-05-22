/*
 * PRIMARYVARIABLES.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef PRIMARYVARIABLES_H_
#define PRIMARYVARIABLES_H_

#include "../include/DynamicalVariables.h"

void getInferredVariables(PRECISION t, const PRECISION * const __restrict__ q, PRECISION ePrev, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, PRECISION utPrev, PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un, PRECISION rhobPrev, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ T, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ seq, PRECISION *  const __restrict__ equiPhiQ);

void setInferredVariablesKernel(const CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, const FLUID_VELOCITY * const __restrict__ uPrev, FLUID_VELOCITY * const __restrict__ u, PRECISION t, void * latticeParams, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ);

PRECISION Ttt(PRECISION e, PRECISION p, PRECISION ut, PRECISION pitt);
PRECISION Ttx(PRECISION e, PRECISION p, PRECISION ut, PRECISION ux, PRECISION pitx);
PRECISION Tty(PRECISION e, PRECISION p, PRECISION ut, PRECISION uy, PRECISION pity);
PRECISION Ttn(PRECISION e, PRECISION p, PRECISION ut, PRECISION un, PRECISION pitn);
PRECISION Txx(PRECISION e, PRECISION p, PRECISION ux, PRECISION pixx);
PRECISION Txy(PRECISION e, PRECISION p, PRECISION ux, PRECISION uy, PRECISION pixy);
PRECISION Txn(PRECISION e, PRECISION p, PRECISION ux, PRECISION un, PRECISION pixn);
PRECISION Tyy(PRECISION e, PRECISION p, PRECISION uy, PRECISION piyy);
PRECISION Tyn(PRECISION e, PRECISION p, PRECISION uy, PRECISION un, PRECISION piyn);
PRECISION Tnn(PRECISION e, PRECISION p, PRECISION un, PRECISION pinn, PRECISION t);

PRECISION Nbt(PRECISION rhob, PRECISION ut, PRECISION nbt);//Lipei

#endif /* PRIMARYVARIABLES_H_ */
