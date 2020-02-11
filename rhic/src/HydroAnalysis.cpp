//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <string>
#include <iomanip>

#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/PrimaryVariables.h"
#include "../include/EquationOfState.h"
#include "../include/HydroAnalysis.h"
#include "../include/HydroPlus.h"
#include "../include/TransportCoefficients.h"
#include "../include/FluxLimiter.h"

#define HBARC 0.197326938
#define zFREQ 25
#define tFREQ 50

using namespace std;

void outputHydroPlus(double t, const char *pathToOutDir, void * latticeParams) {
    
    FILE *fpgamma0, *fpgamma1, *fpgamma2, *fpxi;
    FILE *fpdp, *fpds, *fpdt, *fpdmu;
    char fname[255];
    
    sprintf(fname, "%s/gammaQ0_%.3f.dat", pathToOutDir, t);
    fpgamma0=fopen(fname, "w");
    
    sprintf(fname, "%s/gammaQ1_%.3f.dat", pathToOutDir, t);
    fpgamma1=fopen(fname, "w");
    
    sprintf(fname, "%s/gammaQ2_%.3f.dat", pathToOutDir, t);
    fpgamma2=fopen(fname, "w");
    
    sprintf(fname, "%s/xi_%.3f.dat", pathToOutDir, t);
    fpxi=fopen(fname, "w");
    
    sprintf(fname, "%s/deltaP_%.3f.dat", pathToOutDir, t);
    fpdp=fopen(fname, "w");
    
    sprintf(fname, "%s/deltaS_%.3f.dat", pathToOutDir, t);
    fpds=fopen(fname, "w");
    
    sprintf(fname, "%s/deltaBeta_%.3f.dat", pathToOutDir, t);
    fpdt=fopen(fname, "w");
    
    sprintf(fname, "%s/deltaAlphaB_%.3f.dat", pathToOutDir, t);
    fpdmu=fopen(fname, "w");
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double x,y,z;
    
    int i,j,k;
    int s;
    
    double Q = Qvec[1];
    
    for(k = 2; k < nz+2; ++k) {
        z = (k-2 - (nz-1)/2.)*dz;
        
        for(j = 2; j < ny+2; ++j) {
            y = (j-2 - (ny-1)/2.)*dy;
            
            for(i = 2; i < nx+2; ++i) {
                x = (i-2 - (nx-1)/2.)*dx;
                
                s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                
                double Ts = T[s];
                double alphaBs = alphaB[s];
                double seqs = seq[s];
                double rhobs = rhob[s];
                
                // output correlation length and relaxation rate
                double corrL = correlationLength(Ts, Ts*alphaBs);
                double corrL2 = corrL * corrL;
                
                double gammaPhi = relaxationCoefficientPhi(rhobs, seqs, Ts, corrL2);
                double gammaQ0 = relaxationCoefficientPhiQ(gammaPhi, corrL2, Qvec[0]);
                double gammaQ1 = relaxationCoefficientPhiQ(gammaPhi, corrL2, Qvec[1]);
                double gammaQ2 = relaxationCoefficientPhiQ(gammaPhi, corrL2, Qvec[2]);
                
                fprintf(fpgamma0, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,gammaQ0);
                fprintf(fpgamma1, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,gammaQ1);
                fprintf(fpgamma2, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,gammaQ2);
                fprintf(fpxi, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,corrL);
                
                // output contributions from slow modes to p, T, s and alphaB
                PRECISION deltaVariables[4], equiPhiQ[NUMBER_SLOW_MODES], PhiQ[NUMBER_SLOW_MODES], pPlus;
                
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    equiPhiQ[n] = eqPhiQ->phiQ[n][s];
                    PhiQ[n] = q->phiQ[n][s];
                }
                
                getPressurePlusFromSlowModes(deltaVariables, &pPlus, equiPhiQ, PhiQ, e[s], rhobs, p[s], Ts, alphaBs, seqs);
                
                fprintf(fpdmu, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,deltaVariables[0]);
                fprintf(fpdt, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,deltaVariables[1]);
                fprintf(fpdp, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,deltaVariables[2]);
                fprintf(fpds, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,deltaVariables[3]);
            }
        }
    }
    
    fclose(fpgamma0);
    fclose(fpgamma1);
    fclose(fpgamma2);
    fclose(fpxi);
    fclose(fpdp);
    fclose(fpds);
    fclose(fpdt);
    fclose(fpdmu);
}


void outputBaryonCP(double t, const char *pathToOutDir, void * latticeParams)
{
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    
    int d_ncx = lattice->numLatticePointsX;
    int d_ncy = lattice->numLatticePointsY;
    int d_ncz = lattice->numLatticePointsRapidity;
    
    double d_dx = lattice->latticeSpacingX;
    double d_dy = lattice->latticeSpacingY;
    double d_dz = lattice->latticeSpacingRapidity;
    double d_dt = lattice->latticeSpacingProperTime;
    
    FILE *fpn;
    char fname[255];
    
    int nt = t/d_dt;
    if(nt % tFREQ == 0){
        for(int i = 2; i < d_ncx+2; ++i) {
            for(int j = 2; j < d_ncy+2; ++j) {
                
                // only print out the center of the transverse plane
                //if(i== (d_ncx+3)/2 && j== (d_ncy+3)/2){
                    
                    for(int k = 2; k < d_ncz+2; ++k) {
                        
                        // print out info at some rapidities, with zPREQ
                        if((k-2 - (d_ncz-1)/2.) >= 0 && (k-2) % zFREQ == 0){
                            
                            double x = (i-2 - (d_ncx-1)/2.) * d_dx;
                            double y = (j-2 - (d_ncy-1)/2.) * d_dy;
                            double z = (k-2 - (d_ncz-1)/2.) * d_dz;
                            
                            int s = columnMajorLinearIndex(i, j, k, d_ncx+4, d_ncy+4);
                                                    
                            PRECISION es = e[s];
                            PRECISION ps  = p[s];
                            PRECISION Ts  = T[s];
                            PRECISION rhobs = rhob[s];
                            PRECISION alphaBs  = alphaB[s];
                            PRECISION muBs = Ts * alphaBs;
                            PRECISION corrL = correlationLength(Ts, muBs);
                                                        
                            sprintf(fname, "%s/BaryonCP_%.3f.dat", pathToOutDir, z);
                            fpn = fopen(fname, "a+");
                            
                            fprintf(fpn, "%.3f\t%.3f\t%.3f\t%.3f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",t,x,y,z,muBs,Ts,es,rhobs,corrL);
                            
                            fclose(fpn);
                        }
                    }
                //}
            }
        }
    }
}

void outputAnalysisa(int n, double t, FILE *fpan, void * latticeParams)
{
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    int dnz = floor(nz/2/zFREQ) + 1;
    double zv[dnz], Tv[dnz], muBv[dnz];
    
    if((n-1) % tFREQ == 0){
        
        int m = 0;
        
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    
                    if((k-2 - (nz-1)/2) >= 0 && (k-2 - (nz-1)/2) % zFREQ == 0){
                        
                        if(i == 2 && j == 2){
                            double z = (k-2 - (nz-1)/2.) * dz;
                            
                            int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                            
                            zv[m] = z;
                            Tv[m] = T[s];
                            muBv[m] = T[s] * alphaB[s];
                                                    
                            m++;
                        }
                    }
                }
            }
        }
        
        fprintf(fpan, "%.8f\t",t);
        for(int l = 0; l < dnz; ++l){
            fprintf(fpan, "%.8f\t%.8f\t%.8f\t",zv[l],muBv[l],Tv[l]);
        }
        fprintf(fpan, "\n");
    }
}

// To test the interpolatin function to see it reproduce the EOS table
void testEOS(){

    char EOStable1[] = "output/EOS_dPdT1.dat";
    ofstream eos_table1(EOStable1);
    
    char EOStable2[] = "output/EOS_dPdT2.dat";
    ofstream eos_table2(EOStable2);

    for(int i = 0; i < 13; ++i) {
        for(int j = 0; j < 500; ++j){

            PRECISION etest = (0.0 + i * 0.0003)/HBARC;
            PRECISION rhobtest = (0.0 + j * 0.00001)/HBARC;

            eos_table1 << setprecision(6) << setw(18) << etest*HBARC << setprecision(6) << setw(18) << rhobtest
                       << setprecision(6) << setw(18) << dPdT(etest, rhobtest) << endl;

        }
    }
    
    for(int i = 0; i < 20; ++i) {
        for(int j = 0; j < 300; ++j){

            PRECISION etest = (0.0036 + i * 0.0006)/HBARC;
            PRECISION rhobtest = (0.0 + j * 0.00005)/HBARC;

            eos_table2 << setprecision(6) << setw(18) << etest*HBARC << setprecision(6) << setw(18) << rhobtest
                       << setprecision(6) << setw(18) << dPdT(etest, rhobtest) << endl;

        }
    }
    
    eos_table1.close();
    eos_table2.close();
    printf("EOS table is reproduced.\n");
}

void testCorreLength(){

    char EOStable1[] = "output/correL.dat";
    ofstream eos_table1(EOStable1);

    for(int i = 0; i < 101; ++i) {
        for(int j = 0; j < 86; ++j){

            PRECISION Ttest = (0.07 + i * 0.002)/HBARC;
            PRECISION muBtest = (0.11 + j * 0.002)/HBARC;

            eos_table1 << setprecision(6) << setw(18) << Ttest*HBARC << setprecision(6) << setw(18) << muBtest*HBARC
                       << setprecision(6) << setw(18) << correlationLength(Ttest, muBtest) << endl;
        }
    }
    
    eos_table1.close();
    printf("CorreL table is reproduced.\n");
}

void testBaryCoeff(){
    
    char EOStable1[] = "output/baryonCoefficients.dat";
    ofstream eos_table1(EOStable1);
    
    
    for(int i = 0; i < 101; ++i) {
        for(int j = 0; j < 101; ++j){
            
            PRECISION diffusionCoeff[2];

            PRECISION Ttest = (50.0 + i * 7.0)*0.001/HBARC;
            PRECISION muBtest = (1.0 + j * 7.0)*0.001/HBARC;
            
            baryonDiffusionCoefficient(Ttest, muBtest, diffusionCoeff);

            eos_table1
            << setprecision(6) << setw(18) << Ttest*1000*HBARC << setprecision(6) << setw(18) << muBtest*1000*HBARC
            << setprecision(6) << setw(18) << diffusionCoeff[0]/Ttest/Ttest << setprecision(6) << setw(18) << diffusionCoeff[1]*Ttest << endl;

        }
    }
    
    eos_table1.close();
    printf("Baryon coeff table is reproduced.\n");
}


void testHydroPlus()
{

    // make derivatives of correlation length tables
    
    FILE *filedxi;
    filedxi = fopen ("output/dxi_table.dat","w");
    
    float de = 0.02;
    float dn = 0.005;
    
    for(int i = 0; i < 301; ++i) {
        for(int j = 0; j < 1; ++j) {
            
            float es = i * de;
            float rhobs = j * dn;
            
            float dlogxide = dlnXide(es,rhobs);
            float dlogxidn = dlnXidrhob(es,rhobs);
            
            fprintf(filedxi, "%.3f\t%.3f\t%.3f\t%.3f\n",es,rhobs,dlogxide,dlogxidn);
        }
    }
    
    fclose(filedxi);
     
    // correlation length as a function of (e, rhob)

    char xitable[] = "output/xi_enb_table.dat";
    ofstream xifile(xitable);
    for(int i = 0; i < 900; ++i) {
        for(int j = 0; j < 180; ++j) {
            
            float e = i * 0.02;
            float rhob = j * 0.005;
            
            PRECISION PrimaryVariables[3];
            
            getPrimaryVariablesCombo(e, rhob, PrimaryVariables);

            float Teq = PrimaryVariables[1];
            float alphaBeq = PrimaryVariables[2];
            float muB = alphaBeq*Teq;
            
            float xi = correlationLength(Teq, muB);
            float logxi = log(xi);
            
            if(logxi<0.0) printf("Teq=%f,\t muB=%f,\t xi=%f,\t lnxi=%f\n",Teq,muB,xi,logxi);
            
            xifile
            << setprecision(5) << setw(10) << e//*HBARC
            << setprecision(5) << setw(10) << rhob
            << setprecision(6) << setw(18) << Teq*HBARC
            << setprecision(6) << setw(18) << muB*HBARC
            << setprecision(6) << setw(18) << logxi
            << endl;
        }
    }
    xifile.close();
}


/*PRECISION baryonDiffusionConstant(PRECISION T, PRECISION mub){
 PRECISION T0 = T*HBARC*1000;
 PRECISION mub0 = mub*HBARC*1000;
 if((100<=T0)&&(T0<=450)){
 if((0<=mub0)&&(mub0<=400))
 return InferredPrimaryVariable(mub0, T0-100, 0, 5, 71, 5, 0, 0, EOState->sigmaB)/HBARC/1000 * T;
 else
 return InferredPrimaryVariable(400, T0-100, 0, 5, 71, 5, 0, 0, EOState->sigmaB)/HBARC/1000 * T;
 }else if(T0<100)
 {
 if((0<=mub0)&&(mub0<=400))
 return InferredPrimaryVariable(mub0, 0, 0, 5, 71, 5, 0, 0, EOState->sigmaB)/HBARC/1000 * T;
 else
 return 0.0543361/HBARC/1000 * T;
 }else
 {
 if((0<=mub0)&&(mub0<=400))
 return InferredPrimaryVariable(mub0, 350, 0, 5, 71, 5, 0, 0, EOState->sigmaB)/HBARC/1000 * 2.28; //2.28 [1/fm^4] = 450 MeV
 else
 return 22.5093/HBARC/1000 * 2.28;
 }
 }*/

