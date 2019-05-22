/*
 * InitialConditionParameters.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef INITIALCONDITIONPARAMETERS_H_
#define INITIALCONDITIONPARAMETERS_H_

#include <libconfig.h>

struct InitialConditionParameters
{
	int initialConditionType;
	int numberOfNucleonsPerNuclei;
    int sourceType;
    int numberOfSourceFiles;
    
    double initialBaryonDensity;
	double initialEnergyDensity;
	double scatteringCrossSectionNN;
	double impactParameter;
	double fractionOfBinaryCollisions;

	// longitudinal energy density profile parameters
	double rapidityVariance; // \sigma^{2}_{\eta}
	double rapidityMean; // flat region around \ets_s = 0
    
    double bRapidityVariance1;
    double bRapidityVariance2;
    double bRapidityMean;
};

void loadInitialConditionParameters(config_t *cfg, const char* configDirectory, void * params);

#endif /* INITIALCONDITIONPARAMETERS_H_ */
