//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef GLAUBERMODEL_H_
#define GLAUBERMODEL_H_

double woodsSaxonDistribution(double r, double A);
void energyDensityTransverseProfileAA(double * const __restrict__ energyDensityTransverse, int nx, int ny, double dx, double dy, void * initCondParams, double * const __restrict__ TA, double * const __restrict__ TB, int * const __restrict__ wn);


#endif /* GLAUBERMODEL_H_ */
