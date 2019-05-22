/*
 * SemiDiscreteKurganovTadmorScheme.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef SEMIDISCRETEKURGANOVTADMORSCHEME_H_
#define SEMIDISCRETEKURGANOVTADMORSCHEME_H_

#include "../include/DynamicalVariables.h"

void flux(const PRECISION * const __restrict__ data, PRECISION * const __restrict__ result,
		PRECISION (* const rightHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
		PRECISION (* const leftHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
		PRECISION (* const spectralRadius)(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un),
		PRECISION (* const fluxFunction)(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un),
		PRECISION t, PRECISION ePrev, PRECISION rhobPrev, PRECISION utPrev
);

#endif /* SEMIDISCRETEKURGANOVTADMORSCHEME_H_ */
