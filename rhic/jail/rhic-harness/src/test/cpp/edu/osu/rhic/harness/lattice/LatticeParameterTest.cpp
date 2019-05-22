/*
 * LatticeParameterTest.cpp
 *
 *  Created on: Oct 29, 2015
 *      Author: bazow
 */

#include "gtest/gtest.h"
#include <libconfig.h>
#include <unistd.h>

#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"

TEST(loadLatticeParameters, LatticeParametersFromConfFile) {
	struct LatticeParameters params;
	config_t config;
	config_init(&config);

	char *rootDirectory = NULL;
	size_t size;
	char pathToConfigFile[255];
	rootDirectory = getcwd(rootDirectory,size);
	sprintf(pathToConfigFile, "%s/rhic/rhic-harness/src/test/resources", rootDirectory);
	loadLatticeParameters(&config, pathToConfigFile, &params);
	config_destroy(&config);
	EXPECT_EQ(10, params.numLatticePointsX);
	EXPECT_EQ(15, params.numLatticePointsY);
	EXPECT_EQ(20, params.numLatticePointsRapidity);
	EXPECT_EQ(133, params.numProperTimePoints);

	EXPECT_EQ(0.1, params.latticeSpacingX);
	EXPECT_EQ(0.2, params.latticeSpacingY);
	EXPECT_EQ(1, params.latticeSpacingRapidity);
	EXPECT_EQ(0.5, params.latticeSpacingProperTime);
}

TEST(loadLatticeParameters, DefaultLatticeParameters) {
	struct LatticeParameters params;
	config_t config;
	config_init(&config);
	loadLatticeParameters(&config, "", &params);
	config_destroy(&config);
	EXPECT_EQ(128, params.numLatticePointsX);
	EXPECT_EQ(128, params.numLatticePointsY);
	EXPECT_EQ(64, params.numLatticePointsRapidity);
	EXPECT_EQ(10, params.numProperTimePoints);

	EXPECT_EQ(0.08, params.latticeSpacingX);
	EXPECT_EQ(0.08, params.latticeSpacingY);
	EXPECT_EQ(0.3, params.latticeSpacingRapidity);
	EXPECT_EQ(0.01, params.latticeSpacingProperTime);
}

