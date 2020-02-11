//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef SOURCETERMS_H_
#define SOURCETERMS_H_

#include "../include/DynamicalVariables.h"

void loadSourceTermsX(const PRECISION * const __restrict__ I, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION d_dx);

void loadSourceTermsY(const PRECISION * const __restrict__ J, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION d_dy);

void loadSourceTermsZ(const PRECISION * const __restrict__ K, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION t, PRECISION d_dz);

void loadSourceTerms2(const PRECISION * const __restrict__ Q, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp, PRECISION t, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ pvec, int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_dt, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz, const DYNAMICAL_SOURCES * const __restrict__ Source, const PRECISION * const __restrict__ rhobvec, PRECISION rhobp, const PRECISION * const __restrict__ alphaBvec, PRECISION alphaBp, const PRECISION * const __restrict__ Tvec, PRECISION Tp, PRECISION seq, const SLOW_MODES *  const __restrict__ eqPhiQ, void * hydroParams);

#endif /* SOURCETERMS_H_ */
