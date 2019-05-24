//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include "gtest/gtest.h"
#include <libconfig.h>
#include<unistd.h>

#include "../include/HydroParameters.h"

TEST(loadHydroParameters, HyrdoParametersFromConfFile) {
	struct HydroParameters params;
	config_t config;
	config_init(&config);

	char *rootDirectory = NULL;
	size_t size;
	char pathToConfigFile[255];
	rootDirectory = getcwd(rootDirectory,size);
	sprintf(pathToConfigFile, "%s/rhic/rhic-harness/src/test/resources", rootDirectory);
	loadHydroParameters(&config, pathToConfigFile, &params);
	config_destroy(&config);
	EXPECT_EQ(0.5, params.initialProperTimePoint);
	EXPECT_EQ(0.2, params.shearViscosityToEntropyDensity);
}

TEST(loadHydroParameters, DefaultHyrdoParameters) {
	struct HydroParameters params;
	config_t config;
	config_init(&config);
	loadHydroParameters(&config, "", &params);
	config_destroy(&config);
	EXPECT_EQ(0.1, params.initialProperTimePoint);
	EXPECT_EQ(0.0795775, params.shearViscosityToEntropyDensity);
}
