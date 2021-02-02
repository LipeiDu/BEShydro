//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef HYDROPARAMETERS_H_
#define HYDROPARAMETERS_H_

#include <libconfig.h>

#include "../include/DynamicalVariables.h"

struct HydroParameters
{
	double initialProperTimePoint;
	double shearViscosityToEntropyDensity;
	double freezeoutEnergyDensityGeV;
	int initializePimunuNavierStokes;
    int initializePiNavierStokes;
    int kappaType;
    int gradientType;
    int criticalSlowingDown;
    double cB;
    double Tc;
    double muc;
    int criticalRelaxationTime;
    double correlationLengthMax;
};

void loadHydroParameters(config_t *cfg, const char* configDirectory, void * params);

#endif /* HYDROPARAMETERS_H_ */
