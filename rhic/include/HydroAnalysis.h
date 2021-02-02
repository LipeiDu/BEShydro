//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef HydroAnalysis_h
#define HydroAnalysis_h

void testEOS();
void testBaryCoeff();
void testHydroPlus();
void testCorreLength(double xi_max, double muc);

void outputAnalysisa(int n, double t, FILE *fpan1, FILE *fpan2, FILE *fpan3, void * latticeParams, void * hydroParams);
void outputHydroPlus(double t, const char *pathToOutDir, void * latticeParams);
void outputBaryonCP(double t, const char *pathToOutDir, void * latticeParams);

#endif /* HydroAnalysis_h */
