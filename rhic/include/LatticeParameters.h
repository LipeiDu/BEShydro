/*
 * LatticeParameters.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef LATTICEPARAMETERS_H_
#define LATTICEPARAMETERS_H_

#include <libconfig.h>

#define N_GHOST_CELLS_M 2
#define N_GHOST_CELLS_P 2
#define N_GHOST_CELLS 4

struct LatticeParameters
{
	int numLatticePointsX;
	int numLatticePointsY;
	int numLatticePointsRapidity;
	int numComputationalLatticePointsX;
	int numComputationalLatticePointsY;
	int numComputationalLatticePointsRapidity;
	int numProperTimePoints;

	double latticeSpacingX;
	double latticeSpacingY;
	double latticeSpacingRapidity;
	double latticeSpacingProperTime;
};

void loadLatticeParameters(config_t *cfg, const char* configDirectory, void * params);

#endif /* LATTICEPARAMETERS_H_ */
