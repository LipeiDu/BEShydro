/*
 ============================================================================
 Name        : Run.c
 Author      : Dennis Bazow
 Version     :
 Copyright   :
 Description : Run viscous hydrodynamic simulation of a relativistic heavy ion collision
 ============================================================================
 */

#include <stdlib.h>
#include <stdio.h> // for printf
#include <sys/time.h> // for timing
#include <unistd.h>		// for current working directory
#include <libconfig.h>

//#include "gtest/gtest.h" // for unit testing

#include "edu/osu/rhic/harness/cli/CommandLineArguments.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"
#include "edu/osu/rhic/harness/hydro/HydroPlugin.h"

const char *version = "";
const char *address = "";
/*
int runTest(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
*/
void runHydro(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, const char *outputDir) {
	run(latticeParams, initCondParams, hydroParams, rootDirectory, outputDir);
}

int main(int argc, char **argv) {

	struct CommandLineArguments cli;
	struct LatticeParameters latticeParams;
	struct InitialConditionParameters initCondParams;
	struct HydroParameters hydroParams;

	loadCommandLineArguments(argc, argv, &cli, version, address);

	char *rootDirectory = NULL;
	size_t size;
	rootDirectory = getcwd(rootDirectory,size);

	// Print argument values
	printf("configDirectory = %s\n", cli.configDirectory);
	printf("outputDirectory = %s\n", cli.outputDirectory);
	if (cli.runHydro)
		printf("runHydro = True\n");
	else
		printf("runHydro = False\n");
	if (cli.runTest)
		printf("runTest = True\n");
	else
		printf("runTest = False\n");

	//=========================================
	// Set parameters from configuration files
	//=========================================
	config_t latticeConfig, initCondConfig, hydroConfig;

	// Set lattice parameters from configuration file
	config_init(&latticeConfig);
	loadLatticeParameters(&latticeConfig, cli.configDirectory, &latticeParams);
	config_destroy(&latticeConfig);
	// Set initial condition parameters from configuration file
	config_init(&initCondConfig);
	loadInitialConditionParameters(&initCondConfig, cli.configDirectory, &initCondParams);
	config_destroy(&initCondConfig);
	// Set hydrodynamic parameters from configuration file
	config_init(&hydroConfig);
	loadHydroParameters(&hydroConfig, cli.configDirectory, &hydroParams);
	config_destroy (&hydroConfig);

	//=========================================
	// Run tests
	//=========================================
	/*
	if (cli.runTest) {
		int status = runTest(argc, argv);		
		printf("Done tests.\n");
	}
	*/
	//=========================================
	// Run hydro
	//=========================================
	if (cli.runHydro) {
		runHydro(&latticeParams, &initCondParams, &hydroParams, rootDirectory, cli.outputDirectory);
		printf("Done hydro.\n");
	}

	// TODO: Probably should free host memory here since the freezeout plugin will need
	// to access the energy density, pressure, and fluid velocity.

	return 0;
}
