//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef FILEIO_H_
#define FILEIO_H_

#include "../include/DynamicalVariables.h"

void output(const PRECISION * const var, double t, const char *pathToOutDir, const char *name, void * latticeParams);

#endif /* FILEIO_H_ */
