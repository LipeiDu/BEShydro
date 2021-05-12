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
#define zFREQ 1
#define xFREQ 5
#define tFREQ 5

using namespace std;

PRECISION BaryonDiffusionNS(PRECISION t, int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_dt, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz, const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ alphaBpvec, const PRECISION * const __restrict__ alphaBvec){
    
    PRECISION facX = 1/d_dx/2;
    PRECISION facY = 1/d_dy/2;
    PRECISION facZ = 1/d_dz/2;
    
    int stride = d_ncx * d_ncy;
    
    
    PRECISION *utvec = u->ut;
    PRECISION *uxvec = u->ux; 
    PRECISION *uyvec = u->uy;
    PRECISION *unvec = u->un;

    PRECISION ut = utvec[s];
    PRECISION ux = uxvec[s];
    PRECISION uy = uyvec[s];
    PRECISION un = unvec[s];
    
    PRECISION alphaBp = alphaBpvec[s];
    PRECISION alphaBs = alphaBvec[s];
    
    // derivatives of muB/T
    PRECISION dtalphaB = (alphaBs - alphaBp) / d_dt;
    PRECISION dxalphaB = (*(alphaBvec + s + 1) - *(alphaBvec + s - 1)) * facX;
    PRECISION dyalphaB = (*(alphaBvec + s + d_ncx) - *(alphaBvec + s - d_ncx)) * facY;
    PRECISION dnalphaB = (*(alphaBvec + s + stride) - *(alphaBvec + s - stride)) * facZ;
    
    // gradient of muB/T
    PRECISION ukdk_alphaB = ut * dtalphaB + ux * dxalphaB + uy * dyalphaB + un * dnalphaB;
    PRECISION Nablat_alphaB =  dtalphaB - ut * ukdk_alphaB;
    PRECISION Nablax_alphaB = -dxalphaB - ux * ukdk_alphaB;
    PRECISION Nablay_alphaB = -dyalphaB - uy * ukdk_alphaB;
    PRECISION Nablan_alphaB = -1/pow(t,2)*dnalphaB - un * ukdk_alphaB;
    
    return Nablan_alphaB;
}


PRECISION BaryonDiffusionKnInvRe(PRECISION t, int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_dt, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz, PRECISION Cb, const FLUID_VELOCITY * const __restrict__ u, const FLUID_VELOCITY * const __restrict__ up, const PRECISION * const __restrict__ rhobvec, PRECISION * const __restrict__ Kn, PRECISION * const __restrict__ InvRe){
    
    PRECISION facX = 1/d_dx/2;
    PRECISION facY = 1/d_dy/2;
    PRECISION facZ = 1/d_dz/2;
    
    int stride = d_ncx * d_ncy;
    
    PRECISION t2 = t*t;
    
    // expansion rate
    
    PRECISION *utpvec = up->ut;
    PRECISION utp = utpvec[s];
    
    PRECISION *utvec = u->ut;
    PRECISION *uxvec = u->ux; 
    PRECISION *uyvec = u->uy;
    PRECISION *unvec = u->un;

    PRECISION ut = utvec[s];
    PRECISION ux = uxvec[s];
    PRECISION uy = uyvec[s];
    PRECISION un = unvec[s];
    
    PRECISION dtut = (ut - utp) / d_dt;
    PRECISION dxux = (*(uxvec + s + 1) - *(uxvec + s - 1)) * facX;
    PRECISION dyuy = (*(uyvec + s + d_ncx) - *(uyvec + s - d_ncx)) * facY;
    PRECISION dnun = (*(unvec + s + stride) - *(unvec + s - stride)) * facZ;
    
    PRECISION theta = ut / t + dtut + dxux + dyuy + dnun;
    
    // magnitude of baryon diffusion
#ifdef VMU
    PRECISION nbt = q->nbt[s];
    PRECISION nbx = q->nbx[s];
    PRECISION nby = q->nby[s];
    PRECISION nbn = q->nbn[s];
#else
    PRECISION nbt = 0.;
    PRECISION nbx = 0.;
    PRECISION nby = 0.;
    PRECISION nbn = 0.;
#endif
    
    PRECISION normBaryon = fabs(rhob[s]);
    
    PRECISION nbnb = nbt*nbt - nbx*nbx - nby*nby - nbn*nbn*t2;
    PRECISION snbnb = sqrt(fabs(nbnb));
    if(isnan(snbnb)) printf("found snbnb Nan\n");
    
    // relaxation time
    PRECISION taun = Cb/T[s];
    
    // Knudsen number and inverse Renolds number
    *Kn = taun * theta;
    *InvRe = snbnb / (normBaryon + 1.e-2);
}


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

void outputAnalysis3D(int n, double t, FILE *fpan1, FILE *fpan2, FILE *fpan3, void * latticeParams, void * hydroParams)
{
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    int ncx = lattice->numComputationalLatticePointsX;
    int ncy = lattice->numComputationalLatticePointsY;
    int ncz = lattice->numComputationalLatticePointsRapidity;
    
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    double dt = lattice->latticeSpacingProperTime;
    
    double xi_max = hydro->correlationLengthMax;
    double muc = hydro->muc;
    PRECISION Cb = (PRECISION)(hydro->cB);
    
    double zs, xs, Ts, muBs, corrLs, kappaBs, tauns, qns, seqs, rhobs;
    double nablaEtaAlphaB;
    
    if((n-1) % tFREQ == 0){// time steps
        
        for(int k = 2; k < nz+2; ++k) {
            if((k-2 - (nz-1)/2) >= 0 && (k-2 - (nz-1)/2) % zFREQ == 0){
                zs = (k-2 - (nz-1)/2.) * dz;// eta steps
                
                int j = (ny+1)/2;// y=0
                
                    for(int i = 2; i < nx+2; ++i) {
                        if((i-2 - (nx-1)/2) >= 0 && (i-2 - (nx-1)/2) % xFREQ == 0){
                            xs = (i-2 - (nx-1)/2.) * dx;// x steps

                            int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                            
                            Ts = T[s];
                            muBs = T[s] * alphaB[s];
                            
                            corrLs = corrLen(Ts, muBs, xi_max, muc);
                            
                            kappaBs = baryonDiffusionCoefficientKinetic(Cb, T[s], rhob[s], alphaB[s], e[s], p[s]);
                            
                            tauns = Cb/Ts;
#ifdef VMU                            
                            qns = q->nbn[s];
#else
                            qns = 0.0;
#endif                            
                            seqs = seq[s];
                                
                            rhobs = rhob[s];
                                                     
                            fprintf(fpan1, "%.4f\t%.4f\t%.4f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", t, zs, xs, muBs, Ts, corrLs, kappaBs, tauns, qns, seqs, rhobs);

#ifdef VMU                            
                            nablaEtaAlphaB = BaryonDiffusionNS(t, s, ncx, ncy, ncz, dt, dx, dy, dz, u, alphaBp, alphaB);
#else
                            nablaEtaAlphaB = 0.0;
#endif
                            fprintf(fpan2, "%.4f\t%.4f\t%.4f\t%.8f\n", t, zs, xs, nablaEtaAlphaB);
                            
                            // Kn and Inverse Re
                            
                            PRECISION Kn, InvRe;
                            
                            BaryonDiffusionKnInvRe(t, s, ncx, ncy, ncz, dt, dx, dy, dz, Cb, u, up, rhob, &Kn, &InvRe);
                            
                            if((i-2 - (nx-1)/2) == 0){
                                
                                fprintf(fpan3, "%.4f\t%.4f\t%.8f\t%.8f\n", t, zs, Kn, InvRe);
                                
                            }
                    }
                }
            }
        }
    }
}


void outputAnalysisa(int n, double t, FILE *fpan1, FILE *fpan2, FILE *fpan3, void * latticeParams, void * hydroParams)
{
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double xi_max = hydro->correlationLengthMax;
    double muc = hydro->muc;
    PRECISION Cb = (PRECISION)(hydro->cB);
    
    int dnz = floor(nz/2/zFREQ) + 1;
    double zv[dnz], Tv[dnz], muBv[dnz], corrLv[dnz], kappaB[dnz], tau_n[dnz];
    
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
                            
                            corrLv[m] = corrLen(Tv[m], muBv[m], xi_max, muc);
                            
                            kappaB[m] = baryonDiffusionCoefficientKinetic(Cb, T[s], rhob[s], alphaB[s], e[s], p[s]);
                            
                            tau_n[m] = Cb/T[s];
                                                    
                            m++;
                        }
                    }
                }
            }
        }
        
        fprintf(fpan1, "%.8f\t",t);
        fprintf(fpan3, "%.8f\t",t);
        for(int l = 0; l < dnz; ++l){
            fprintf(fpan1, "%.8f\t",zv[l]);
            fprintf(fpan2, "%.8f\t%.8f\t",muBv[l],Tv[l]);
            fprintf(fpan3, "%.8f\t%.8f\t%.8f\t",corrLv[l], kappaB[l], tau_n[l]);
        }
        fprintf(fpan1, "\n");
        fprintf(fpan2, "\n");
        fprintf(fpan3, "\n");
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

void testCorreLength(double xi_max, double muc){

    char EOStable1[] = "output/correL.dat";
    ofstream eos_table1(EOStable1);

    for(int i = 0; i < 201; ++i) {
        for(int j = 0; j < 125; ++j){

            PRECISION Ttest = (0.144 + i * 0.00005)/HBARC;
            PRECISION muBtest = (0.233 + j * 0.00025)/HBARC;

            eos_table1 << setprecision(6) << setw(18) << Ttest*HBARC << setprecision(6) << setw(18) << muBtest*HBARC
                      // << setprecision(6) << setw(18) << correlationLength(Ttest, muBtest) << endl;
                        << setprecision(6) << setw(18) << corrLen(Ttest, muBtest, xi_max, muc) << endl;
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

