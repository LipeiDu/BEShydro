/*
 * FiniteDifference.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include "edu/osu/rhic/core/util/FiniteDifference.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

PRECISION finiteDifferenceX(const PRECISION * const var, int i,int j,int k, int NX, int NY, int NZ, double dx) {
	return 0.5 * (var[i+1 + NX * (j + NY * k)]-var[i-1 + NX * (j + NY * k)])/dx;
}

PRECISION finiteDifferenceY(const PRECISION * const var, int i,int j,int k, int NX, int NY, int NZ, double dx) {
	return 0.5 * (var[i+ NX * (j+1 + NY * k)]-var[i+ NX * (j-1 + NY * k)])/dx;
}
 
PRECISION finiteDifferenceZ(const PRECISION * const var, int i,int j,int k, int NX, int NY, int NZ, double dz) {
	return 0.5 * (var[i+ NX * (j + NY * (k+1))]-var[i+ NX * (j + NY * (k-1))])/dz;
}
