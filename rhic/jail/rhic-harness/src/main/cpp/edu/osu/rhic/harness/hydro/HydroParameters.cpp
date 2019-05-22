/*
 * HydroParameters.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include "edu/osu/rhic/harness/hydro/HydroParameters.h"
#include "edu/osu/rhic/harness/util/Properties.h"

double initialProperTimePoint;
double shearViscosityToEntropyDensity;
double freezeoutTemperatureGeV;
int initializePimunuNavierStokes;

void loadHydroParameters(config_t *cfg, const char* configDirectory, void * params) {
	// Read the file
	char fname[255];
	sprintf(fname, "%s/%s", configDirectory, "hydro.properties");
	if (!config_read_file(cfg, fname)) {
		fprintf(stderr, "No configuration file  %s found for hydrodynamic parameters - %s.\n", fname, config_error_text(cfg));
		fprintf(stderr, "Using default hydrodynamic configuration parameters.\n");
	}

	getDoubleProperty(cfg, "initialProperTimePoint", &initialProperTimePoint, 0.1);
	getDoubleProperty(cfg, "shearViscosityToEntropyDensity", &shearViscosityToEntropyDensity, 0.0795775);
	getDoubleProperty(cfg, "freezeoutTemperatureGeV", &freezeoutTemperatureGeV, 0.155);

	getIntegerProperty(cfg, "initializePimunuNavierStokes", &initializePimunuNavierStokes, 1);

	struct HydroParameters * hydro = (struct HydroParameters *) params;
	hydro->initialProperTimePoint = initialProperTimePoint;
	hydro->shearViscosityToEntropyDensity = shearViscosityToEntropyDensity;
	hydro->freezeoutTemperatureGeV = freezeoutTemperatureGeV;
	hydro->initializePimunuNavierStokes = initializePimunuNavierStokes;
}
