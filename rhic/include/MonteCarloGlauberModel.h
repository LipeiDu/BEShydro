//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef MONTECARLOGLAUBERMODEL_H_
#define MONTECARLOGLAUBERMODEL_H_

void monteCarloGlauberEnergyDensityTransverseProfile(double * const __restrict__ energyDensityTransverse, int nx, int ny, double dx, double dy, void * initCondParams, double * const __restrict__ TA, double * const __restrict__ TB, int * const __restrict__ n1, int * const __restrict__ n2);

#endif /* MONTECARLOGLAUBERMODEL_H_ */
