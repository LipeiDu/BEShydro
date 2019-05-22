/*
 * FileIO.h
 *
 *  Created on: Oct 24, 2015
 *      Author: bazow
 */

#ifndef FILEIO_H_
#define FILEIO_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

void output(const PRECISION * const var, double t, const char *pathToOutDir, const char *name, void * latticeParams);

#endif /* FILEIO_H_ */
