/*
 * GlauberModel.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef GLAUBERMODEL_H_
#define GLAUBERMODEL_H_

double woodsSaxonDistribution(double r, double A);
void energyDensityTransverseProfileAA(double * const __restrict__ energyDensityTransverse, int nx, int ny, double dx, double dy, void * initCondParams, double * const __restrict__ TA, double * const __restrict__ TB//thickness TA&TB and n1&n2 by Lipei
);


#endif /* GLAUBERMODEL_H_ */
