//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <stdlib.h>
#include <stdio.h> // for printf
#include <sys/time.h> // for timing
#include <unistd.h> // for current working directory
#include <libconfig.h>

#include "../include/BEShydroLOGO.h"
#include "../include/CommandLineArguments.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroParameters.h"
#include "../include/HydroPlugin.h"

#ifdef _OPENMP
#include <omp.h>
#endif

const char *version = "";
const char *address = "";

void runHydro(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, const char *outputDir) {
	run(latticeParams, initCondParams, hydroParams, rootDirectory, outputDir);
}

int main(int argc, char **argv) {
    
#ifdef _OPENMP
    int tid = omp_get_max_threads();
    printf("Get %d threads...\n", tid);

    //clock
    double sec = 0.0;
    sec = omp_get_wtime();
#endif

	struct CommandLineArguments cli;
	struct LatticeParameters latticeParams;
	struct InitialConditionParameters initCondParams;
	struct HydroParameters hydroParams;

	loadCommandLineArguments(argc, argv, &cli, version, address);

	char *rootDirectory = NULL;
	size_t size;
	rootDirectory = getcwd(rootDirectory,size);
    
    // show logo and copyright information
    displayLogo();
    displayCopyright();

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
	// Run hydro
	//=========================================
	if (cli.runHydro) {
		runHydro(&latticeParams, &initCondParams, &hydroParams, rootDirectory, cli.outputDirectory);
		printf("Done hydro.\n");
	}
    
#ifdef _OPENMP
    sec = omp_get_wtime() - sec;

    printf("Hydro took %f seconds\n", sec);
#endif

	return 0;
}
