//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include "../include/GlauberModel.h"
#include "../include/InitialConditionParameters.h"

#include <gsl/gsl_integration.h>

struct nuclearThicknessFunctionParams {
	double x;
	double y;
	double A;
};

double woodsSaxonDistribution(double r, double A) {
	double n0 = 0.17; // central density in fm^(-3)
	double Rn = 1.12*pow(A,1./3.) - 0.86*pow(A,-1./3.); // radius in fm
	double d = 0.54; // thickness in fm
	return n0/(1+exp((r-Rn)/d));
}

double nuclearThicknessFunctionIntegrand(double s, void * params) {
	struct nuclearThicknessFunctionParams * ta_params = (struct nuclearThicknessFunctionParams *) params;
	double x = ta_params->x;
	double y = ta_params->y;
	double A = ta_params->A;
	double z = (1.-s)/s;
	double s2 = s*s;
	return woodsSaxonDistribution(sqrt(x * x + y * y + z * z), A)/s2;
}

double nuclearThicknessFunction(double x, double y, double A) {
	double result, error;
	double epsabs = 0.;
	double epsrel = 1.e-8;
	double eps = 1.e-16;
	int n = 1000;
	struct nuclearThicknessFunctionParams params = { x, y, A };
	gsl_function F;
	F.function = &nuclearThicknessFunctionIntegrand;
	F.params = &params;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(n);
	gsl_integration_qags(&F, 0., 1., epsabs, epsrel, n, w, &result, &error);
	gsl_integration_workspace_free(w);
	return 2*result;
}

double binaryCollisionPairs(double x, double y, double TAplus, double TBminus, double snn) {
	return snn * TAplus * TBminus;
}

double woundedNucleons(double x, double y, double TAminus, double TAplus, double TBminus, double TBplus, double A, double B, double snn) {
	double fermiToMilliBarns = 0.1;
	double part, res;
	part = 1-snn*TBminus*fermiToMilliBarns/B;
	res = TAplus*(1-pow(1-part,B));
	part = 1-snn*TAplus*fermiToMilliBarns/A;
	res += TBminus*(1-pow(1-part,A));
	return res;
}

//************************************************************************************\
 * Mixture of wounded nucleon and binary collision energy density profile
//************************************************************************************/
void energyDensityTransverseProfileAA(double * const __restrict__ energyDensityTransverse, int nx, int ny, double dx, double dy, void * initCondParams, double * const __restrict__ TA, double * const __restrict__ TB, int * const __restrict__ wn) {
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	double A = initCond->numberOfNucleonsPerNuclei;
	double B = A;
	double b = initCond->impactParameter;
	double snn = initCond->scatteringCrossSectionNN;
	double alpha = initCond->fractionOfBinaryCollisions;

	// Normalization factors
	double TAminusNorm = nuclearThicknessFunction(0,0,A);
	double TAplusNorm = nuclearThicknessFunction(0,0,A);
	double TBminusNorm = TAminusNorm;
	double TBplusNorm = TAplusNorm;
    
	double nbcNorm = 1./binaryCollisionPairs(0,0,TAplusNorm,TBminusNorm,snn);
    *wn = woundedNucleons(0,0,TAminusNorm,TAplusNorm,TBminusNorm,TBplusNorm,A,B,snn);
    double wnNorm = 1./ (*wn);

	for(int i = 0; i < nx; ++i) {
		double x = (i - (nx-1)/2.)*dx;
		for(int j = 0; j < ny; ++j) {
			double y = (j - (ny-1)/2.)*dy;
			double TAminus = nuclearThicknessFunction(x-b/2,y,A);
			double TAplus  = nuclearThicknessFunction(x+b/2,y,A);
			double TBminus = TAminus;
			double TBplus  = TAplus;
            
            // Thickness function to construct baryon density;
            TA[i+j*nx] = TAplus;
            TB[i+j*nx] = TAminus;
			// Binary collision energy density profile
			double ed = nbcNorm * binaryCollisionPairs(x,y,TAplus,TBminus,snn);
			ed *= alpha;
			// Wounded nucleon energy density profile
			ed += (1-alpha) * wnNorm * woundedNucleons(x,y,TAminus,TAplus,TBminus,TBplus,A,B,snn);
			energyDensityTransverse[i+j*nx] = ed;
		}
	}
}
