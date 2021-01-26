//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

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
    double bNorm;
    
    double energyCutOff;
    double baryonCutOff;
};

void loadInitialConditionParameters(config_t *cfg, const char* configDirectory, void * params);

#endif /* INITIALCONDITIONPARAMETERS_H_ */
