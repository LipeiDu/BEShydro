//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include "../include/HydroParameters.h"
#include "../include/Properties.h"

double initialProperTimePoint;
double shearViscosityToEntropyDensity;
double freezeoutEnergyDensityGeV;
int initializePimunuNavierStokes;
int initializePiNavierStokes;

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
	getDoubleProperty(cfg, "freezeoutEnergyDensityGeV", &freezeoutEnergyDensityGeV, 0.155);

	getIntegerProperty(cfg, "initializePimunuNavierStokes", &initializePimunuNavierStokes, 1);
    getIntegerProperty(cfg, "initializePiNavierStokes", &initializePiNavierStokes, 1);

	struct HydroParameters * hydro = (struct HydroParameters *) params;
	hydro->initialProperTimePoint = initialProperTimePoint;
	hydro->shearViscosityToEntropyDensity = shearViscosityToEntropyDensity;
	hydro->freezeoutEnergyDensityGeV = freezeoutEnergyDensityGeV;
	hydro->initializePimunuNavierStokes = initializePimunuNavierStokes;
    hydro->initializePiNavierStokes = initializePiNavierStokes;
}
