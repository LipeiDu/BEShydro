/*
 * HalfSiteExtrapolation.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef HALFSITEEXTRAPOLATION_H_
#define HALFSITEEXTRAPOLATION_H_

#include "../include/DynamicalVariables.h"

PRECISION rightHalfCellExtrapolationForward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);
PRECISION rightHalfCellExtrapolationBackwards(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);
PRECISION leftHalfCellExtrapolationForward(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);
PRECISION leftHalfCellExtrapolationBackwards(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp);

#endif /* HALFSITEEXTRAPOLATION_H_ */
