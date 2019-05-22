/*
 * InitialConditionParameters.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/harness/util/Properties.h"

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

// longitudinal baryon density profile parameters
double bRapidityVariance1;
double bRapidityVariance2;
double bRapidityMean;

void loadInitialConditionParameters(config_t *cfg, const char* configDirectory, void * params) {

	char fname[255];
	sprintf(fname, "%s/%s", configDirectory, "ic.properties");
	if (!config_read_file(cfg, fname)) {
		fprintf(stderr, "No configuration file  %s found for initial condition parameters - %s.\n", fname, config_error_text(cfg));
		fprintf(stderr, "Using default initial condition configuration parameters.\n");
	}

	getIntegerProperty(cfg, "initialConditionType", &initialConditionType, 2);
	getIntegerProperty(cfg, "numberOfNucleonsPerNuclei", &numberOfNucleonsPerNuclei, 208);

    getIntegerProperty(cfg, "numberOfSourceFiles", &numberOfSourceFiles, 30);
    getIntegerProperty(cfg, "sourceType", &sourceType, 0);
    getDoubleProperty(cfg, "initialBaryonDensity", &initialBaryonDensity, 0.0);
    getDoubleProperty(cfg, "bRapidityVariance1", &bRapidityVariance1, 0.2);
    getDoubleProperty(cfg, "bRapidityVariance2", &bRapidityVariance2, 2.0);
    getDoubleProperty(cfg, "bRapidityMean", &bRapidityMean, 2.0);

	getDoubleProperty(cfg, "initialEnergyDensity", &initialEnergyDensity, 1.0);
	getDoubleProperty(cfg, "scatteringCrossSectionNN", &scatteringCrossSectionNN, 62);
	getDoubleProperty(cfg, "impactParameter", &impactParameter, 7);
	getDoubleProperty(cfg, "fractionOfBinaryCollisions", &fractionOfBinaryCollisions, 0.5);
	getDoubleProperty(cfg, "rapidityVariance", &rapidityVariance, 0.5);
	getDoubleProperty(cfg, "rapidityMean", &rapidityMean, 0.5);

	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) params;
	initCond->initialConditionType = initialConditionType;
	initCond->numberOfNucleonsPerNuclei = numberOfNucleonsPerNuclei;
    initCond->numberOfSourceFiles = numberOfSourceFiles;
    initCond->sourceType = sourceType;
    initCond->initialBaryonDensity = initialBaryonDensity;
	initCond->initialEnergyDensity = initialEnergyDensity;
	initCond->scatteringCrossSectionNN = scatteringCrossSectionNN;
	initCond->impactParameter = impactParameter;
	initCond->fractionOfBinaryCollisions = fractionOfBinaryCollisions;
	initCond->rapidityVariance = rapidityVariance;
	initCond->rapidityMean = rapidityMean;
    initCond->bRapidityVariance1 = bRapidityVariance1;
    initCond->bRapidityVariance2 = bRapidityVariance2;
    initCond->bRapidityMean = bRapidityMean;
}
