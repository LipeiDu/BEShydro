/*
 * FiniteDifference.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef FINITEDIFFERENCE_H_
#define FINITEDIFFERENCE_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

PRECISION finiteDifferenceX(const PRECISION * const var, int i,int j,int k, int NX, int NY, int NZ, double dx);
PRECISION finiteDifferenceY(const PRECISION * const var, int i,int j,int k, int NX, int NY, int NZ, double dx);
PRECISION finiteDifferenceZ(const PRECISION * const var, int i,int j,int k, int NX, int NY, int NZ, double dz);

#endif /* FINITEDIFFERENCE_H_ */
