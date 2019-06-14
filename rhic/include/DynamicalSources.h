//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef DYNAMICALSOURCE_H_
#define DYNAMICALSOURCE_H_

void zeroSource(void * latticeParams, void * initCondParams);

void readInSource(int n, void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory);

void setDynamicalSources(void * latticeParams, void * initCondParams, double *dp_dtau, double *pos); //function to set dynamical source terms for stress tensor and baryon current

#endif /* DYNAMICALSOURCE_H_ */
