/*
 * InitialConditionParameterTest.cpp
 *
 *  Created on: Oct 29, 2015
 *      Author: bazow
 */

#include "gtest/gtest.h"
#include <libconfig.h>
#include<unistd.h>

#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"

TEST(loadInitialConditionParameters, InitialConditionParametersFromConfFile) {
	struct InitialConditionParameters params;
	config_t config;
	config_init(&config);

	char *rootDirectory = NULL;
	size_t size;
	char pathToConfigFile[255];
	rootDirectory = getcwd(rootDirectory,size);
	sprintf(pathToConfigFile, "%s/rhic/rhic-harness/src/test/resources", rootDirectory);
	loadInitialConditionParameters(&config, pathToConfigFile, &params);
	config_destroy(&config);
	EXPECT_EQ(0, params.initialConditionType);
	EXPECT_EQ(63, params.numberOfNucleonsPerNuclei);
	EXPECT_EQ(0.6, params.initialEnergyDensity);
	EXPECT_EQ(52.0, params.scatteringCrossSectionNN);
	EXPECT_EQ(20.0, params.impactParameter);
	EXPECT_EQ(0.1, params.fractionOfBinaryCollisions);
	EXPECT_EQ(0.2, params.rapidityVariance);
	EXPECT_EQ(0.3, params.rapidityMean);
}

TEST(loadInitialConditionParameters, DefaultInitialConditionParameters) {
	struct InitialConditionParameters params;
	config_t config;
	config_init(&config);
	loadInitialConditionParameters(&config, "", &params);
	config_destroy(&config);
	EXPECT_EQ(2, params.initialConditionType);
	EXPECT_EQ(208, params.numberOfNucleonsPerNuclei);
	EXPECT_EQ(1, params.initialEnergyDensity);
	EXPECT_EQ(62, params.scatteringCrossSectionNN);
	EXPECT_EQ(7, params.impactParameter);
	EXPECT_EQ(0.5, params.fractionOfBinaryCollisions);
	EXPECT_EQ(0.5, params.rapidityVariance);
	EXPECT_EQ(0.5, params.rapidityMean);
}
