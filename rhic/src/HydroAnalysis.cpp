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
#define zFREQ 10
#define tFREQ 10

using namespace std;

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
    
    PRECISION facX = 1/d_dx/2;
    PRECISION facY = 1/d_dy/2;
    PRECISION facZ = 1/d_dz/2;
    
    int stride = (d_ncx+4) * (d_ncy+4);
    
    FILE *fpn;
    char fname[255];
    
    for(int i = 2; i < d_ncx+2; ++i) {
        for(int j = 2; j < d_ncy+2; ++j) {
            
            // only print out the center of the transverse plane
            if(i== (d_ncx+3)/2 && j== (d_ncy+3)/2){
                
                for(int k = 2; k < d_ncz+2; ++k) {
                    
                    // print out info at some rapidities, with zPREQ
                    if((k-2) >= 0 && (k-2) % zFREQ == 0){
                        
                        double x = (i-2 - (d_ncx-1)/2.) * d_dx;
                        double y = (j-2 - (d_ncy-1)/2.) * d_dy;
                        double z = (k-2 - (d_ncz-1)/2.) * d_dz;
                        
                        int s = columnMajorLinearIndex(i, j, k, d_ncx+4, d_ncy+4);
                        
                        PRECISION *utvec = u->ut;
                        PRECISION *uxvec = u->ux;
                        PRECISION *uyvec = u->uy;
                        PRECISION *unvec = u->un;
                        
                        PRECISION ps  = p[s];
                        PRECISION ut = utvec[s];
                        PRECISION ux = uxvec[s];
                        PRECISION uy = uyvec[s];
                        PRECISION un = unvec[s];
                        
                        //*********************************************************\
                        //* PART I
                        //*********************************************************/
                        
                        PRECISION es = e[s];
                        PRECISION Ts  = T[s];
                        PRECISION rhobs = rhob[s];
                        PRECISION alphaBs  = alphaB[s];
                        
                        PRECISION Tps  = Tp[s];
                        PRECISION rhobps = rhobp[s];
                        PRECISION alphaBps  = alphaBp[s];
                        
                        // derivatives of rhob
                        PRECISION dtrhob = (rhobs - rhobps) / d_dt;
                        PRECISION dxrhob = (rhob[s+1]- rhob[s - 1]) * facX;
                        PRECISION dyrhob = (rhob[s + d_ncx] - rhob[s - d_ncx]) * facY;
                        PRECISION dnrhob = (rhob[s + stride] - rhob[s - stride]) * facZ;
                        
                        // derivatives of T
                        PRECISION dtT = (Ts - Tps) / d_dt;
                        PRECISION dxT = (*(T + s + 1) - *(T + s - 1)) * facX;
                        PRECISION dyT = (*(T + s + d_ncx) - *(T + s - d_ncx)) * facY;
                        PRECISION dnT = (*(T + s + stride) - *(T + s - stride)) * facZ;
                        
                        // derivatives of muB/T
                        /*PRECISION dtalphaB = (alphaBs - alphaBps) / d_dt;
                         PRECISION dxalphaB = (*(alphaB + s + 1) - *(alphaB + s - 1)) * facX;
                         PRECISION dyalphaB = (*(alphaB + s + d_ncx) - *(alphaB + s - d_ncx)) * facY;
                         PRECISION dnalphaB = (*(alphaB + s + stride) - *(alphaB + s - stride)) * facZ;*/
                        
                        PRECISION alphabxp1 = alphaB[s+1];
                        PRECISION alphabxp2 = alphaB[s+2];
                        PRECISION alphabxm1 = alphaB[s-1];
                        PRECISION alphabxm2 = alphaB[s-2];
                        PRECISION alphabyp1 = alphaB[s+d_ncx];
                        PRECISION alphabyp2 = alphaB[s+2*d_ncx];
                        PRECISION alphabym1 = alphaB[s-d_ncx];
                        PRECISION alphabym2 = alphaB[s-2*d_ncx];
                        PRECISION alphabnp1 = alphaB[s+stride];
                        PRECISION alphabnp2 = alphaB[s+2*stride];
                        PRECISION alphabnm1 = alphaB[s-stride];
                        PRECISION alphabnm2 = alphaB[s-2*stride];
                        
                        PRECISION dtalphaB = (alphaBs - alphaBps) / d_dt;
                        PRECISION dxalphaB = approximateDerivative(alphabxm1,alphaBs,alphabxp1) * facX * 2;
                        PRECISION dyalphaB = approximateDerivative(alphabym1,alphaBs,alphabyp1) * facY * 2;
                        PRECISION dnalphaB = approximateDerivative(alphabnm1,alphaBs,alphabnp1) * facZ * 2;
                        
                        // gradient of rhob
                        PRECISION ukdk_rhob = ut * dtrhob + ux * dxrhob + uy * dyrhob + un * dnrhob;
                        PRECISION Nablat_rhob =  dtrhob - ut * ukdk_rhob;
                        PRECISION Nablax_rhob = -dxrhob - ux * ukdk_rhob;
                        PRECISION Nablay_rhob = -dyrhob - uy * ukdk_rhob;
                        PRECISION Nablan_rhob = -1/pow(t,2)*dnrhob - un * ukdk_rhob;
                        
                        // gradient of T
                        PRECISION ukdk_T = ut * dtT + ux * dxT + uy * dyT + un * dnT;
                        PRECISION Nablat_T =  dtT - ut * ukdk_T;
                        PRECISION Nablax_T = -dxT - ux * ukdk_T;
                        PRECISION Nablay_T = -dyT - uy * ukdk_T;
                        PRECISION Nablan_T = -1/pow(t,2)*dnT - un * ukdk_T;
                        
                        // gradient of muB/T
                        PRECISION ukdk_alphaB = ut * dtalphaB + ux * dxalphaB + uy * dyalphaB + un * dnalphaB;
                        PRECISION Nablat_alphaB =  dtalphaB - ut * ukdk_alphaB;
                        PRECISION Nablax_alphaB = -dxalphaB - ux * ukdk_alphaB;
                        PRECISION Nablay_alphaB = -dyalphaB - uy * ukdk_alphaB;
                        PRECISION Nablan_alphaB = -1/pow(t,2)*dnalphaB - un * ukdk_alphaB;
                        
                        //printf("dnrhob = %f\t dnT = %f\t dnalphaB = %f\t rhob[s + stride] = %f\t rhob[s - stride] = %f\t strie = %d\n",dnrhob,dnT,dnalphaB,rhob[s+stride],rhob[s - stride],stride);
                        
                        //*********************************************************\
                        //* PART II
                        //*********************************************************/
                        
                        PRECISION muBs = Ts * alphaBs;
                        PRECISION corrL = correlationLength(Ts, muBs);
                        
                        // relaxation time
                        PRECISION tau_n0 = Cb/Ts;
                        
                        // baryon diffusion coefficients
                        PRECISION kappaB0 = 0.0;
                        PRECISION DB0 = 0.0;
                        
                        PRECISION diffusionCoeff[2];
                        baryonDiffusionCoefficient(Ts, muBs, diffusionCoeff);
                        
                        PRECISION kappaB = diffusionCoeff[0];//baryonDiffusionCoefficientKinetic(Ts, rhobs, alphaBs, es, ps);;//
                        PRECISION DB = diffusionCoeff[1];
                        
                        // critical scaling
                        //#ifdef CRITICAL
                        //                        PRECISION tau_n = tau_n0 * corrL;
                        //                        PRECISION kappaB = kappaB0 * corrL;
                        //                        PRECISION DB = DB0 / corrL;
                        //#endif
                        
                        // other transport coefficients
                        //PRECISION delta_nn = tau_n;
                        //PRECISION lambda_nn = 0.60 * tau_n;
                        
                        // kappaB * gradient of muB/T term
                        PRECISION NBI0t1 = kappaB * Nablat_alphaB;
                        PRECISION NBI0x1 = kappaB * Nablax_alphaB;
                        PRECISION NBI0y1 = kappaB * Nablay_alphaB;
                        PRECISION NBI0n1 = kappaB * Nablan_alphaB;
                        
                        PRECISION TInv = 1 / Ts;
                        PRECISION facG1 = kappaB / rhobs * TInv;
                        PRECISION facG2 = dPdT(es, rhobs) - (es + ps) * TInv;
                        PRECISION facG = facG1 * facG2;
                        
                        PRECISION NBI0t2 = DB * Nablat_rhob + facG *  Nablat_T;
                        PRECISION NBI0x2 = DB * Nablax_rhob + facG *  Nablax_T;
                        PRECISION NBI0y2 = DB * Nablay_rhob + facG *  Nablay_T;
                        PRECISION NBI0n2 = DB * Nablan_rhob + facG *  Nablan_T;
                        
                        int nt = t/d_dt;
                        if(nt % tFREQ == 0){
                            sprintf(fname, "%s/BaryonCP_%.3f.dat", pathToOutDir, z);
                            fpn = fopen(fname, "a+");
                            
                            //fprintf(fpn, "%.3f\t%.3f\t%.3f\t%.3f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t %.8f\t%.8f\t%.8f\t%.8f\t %.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",
                            //t,x,y,z,muBs,Ts,corrL,DB0,DB,kappaB0,kappaB,NBI0t1,NBI0x1,NBI0y1,NBI0n1,NBI0t2,NBI0x2,NBI0y2,NBI0n2,es,rhobs);//17
                            
                            fprintf(fpn, "%.3f\t%.3f\t%.3f\t%.3f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t %.8f\t%.8f\t%.8f\t%.8f\t %.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",t,x,y,z,muBs,Ts,corrL,dtrhob,dxrhob,dyrhob,dnrhob,dtalphaB,dxalphaB,dyalphaB,dnalphaB,dtT,dxT,dyT,dnT,es,rhobs);//17
                            
                            fclose(fpn);
                        }
                    }
                }
            }
        }
    }
}


void outputAnalysis(double t, FILE *outputDir, void * latticeParams)
{
    /*FILE *fp;
    char fname[255];
    sprintf(fname, "%s/AnalysisData.dat", outputDir);
    fp=fopen(fname, "a+");*/
    
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
    
    /*double v2t,becc,eecc,v2t1,v2t2;
    v2t = 0;
    v2t1 = 0;
    v2t2 = 0;
    becc = 0;
    eecc = 0;
    double bymx = 0;
    double bxy = 0;
    double bypx = 0;
    double eymx = 0;
    double exy = 0;
    double eypx = 0;*/
    
    //PRECISION phiQ;
    
    //k=(nz+3)/2;
    //j=(ny+3)/2;
    //i=(nx+3)/2;
    for(i = 2; i < nx+2; ++i) {
        for(j = 2; j < ny+2; ++j) {
            
            PRECISION mub0,t0,mub1,t1,mub2,t2,mub3,t3,mub4,t4;
            x = (i-2 - (nx-1)/2.) * dx;
            y = (j-2 - (ny-1)/2.) * dy;
            
            if(x==0&&y==0){
                
            for(k = 2; k < nz+2; ++k) {
                
                z = (k-2 - (nz-1)/2.) * dz;
                
                s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                
                if(z==0.0){
                    t0 = T[s];
                    mub0 = t0*alphaB[s];
                }
                if(z==0.5){
                    t1 = T[s];
                    mub1 = t1*alphaB[s];
                }
                if(z==1.0){
                    t2 = T[s];
                    mub2 = t2*alphaB[s];
                }
                if(z==1.5){
                    //printf("here!\n");
                    t3 = T[s];
                    mub3 = t3*alphaB[s];
                }
                if(fabs(z-1.7)<1.e-3){
                    //printf("here\n");
                    t4 = T[s];
                    mub4 = t4*alphaB[s];
                }
                
                
                
                
                    
                
                //if(x==0.0&&y==0.0&&t==1.5)
                //    fprintf(fp, "%.3f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",z,e[s],rhob[s],seq[s],alphaB[s],T[s]);
                
                //if(x==0&&y==0)
                //fprintf(fp, "%.3f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",t,e[s],q->pinn[s],q->Pi[s],rhob[s],q->nbn[s],T[s],p[s]);
                
                //if(x==0&&y==0)
                    //fprintf(fp, "%.3f\t%.8f\t%.8f\t%.8f\n",t,e[s],q->pinn[s],p[s]);
                
                //double tt=Ttt(e[s],p[s],u->ut[s],q->pitt[s]);
                //double tx=Ttx(e[s],p[s],u->ut[s],u->ux[s],q->pitx[s]);
                //double ty=Tty(e[s],p[s],u->ut[s],u->uy[s],q->pity[s]);
                //double tn=Ttn(e[s],p[s],u->ut[s],u->un[s],q->pitn[s]);
/*#ifndef PIMUNU
                double pixx=0;
                double piyy=0;
#else
                double pixx=q->pixx[s];
                double piyy=q->piyy[s];
#endif
                double xx=Txx(e[s],p[s],u->ux[s],pixx);
                //double xy=Txy(e[s],p[s],u->ux[s],u->uy[s],q->pixy[s]);
                //double xn=Txn(e[s],p[s],u->ux[s],u->un[s],q->pixn[s]);
                double yy=Tyy(e[s],p[s],u->uy[s],piyy);
                //double yn=Tyn(e[s],p[s],u->uy[s],u->un[s],q->piyn[s]);
                //double nn=Tnn(e[s],p[s],u->un[s],q->pinn[s],t);
                
                bymx = bymx + (y*y - x*x)*rhob[s];
                bxy = bxy + x*y*rhob[s];
                bypx = bypx + (y*y + x*x)*rhob[s];
                eymx = eymx + (y*y - x*x)*e[s];
                exy = exy + x*y*e[s];
                eypx = eypx + (y*y + x*x)*e[s];
                v2t1 = v2t1 + (xx-yy);
                v2t2 = v2t2 + (xx+yy);*/
                
                //if(x==0&&y==0)
               // phiQ = q->phiQ[0][s];
                
                /*PRECISION eIn = e[s];
                PRECISION rhobIn = rhob[s];
                PRECISION pIn = p[s];
                PRECISION TIn = T[s];
                PRECISION alphaBIn = alphaB[s];
                PRECISION equiPhiQ[NUMBER_SLOW_MODES];
                PRECISION PhiQ[NUMBER_SLOW_MODES];
                
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    
                    equiPhiQ[n] = eqPhiQ->phiQ[n][s];
                    PhiQ[n] = q->phiQ[n][s];//for test
                }
                
                PRECISION p,T,alphaB;

                getPrimaryVariablesFromSlowModes(&p, &T, &alphaB, equiPhiQ, PhiQ, eIn, rhobIn, pIn, TIn, alphaBIn);*/
            }
                fprintf(outputDir, "%.3f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",t,mub0,t0,mub1,t1,mub2,t2,mub3,t3,mub4,t4);
            }
            
            
        }
    }
    
    /*becc = sqrt(bymx*bymx + 4*bxy*bxy)/bypx;
    eecc = sqrt(eymx*eymx + 4*exy*exy)/eypx;
    v2t = v2t1/v2t2;
    
    fprintf(fp, "%.3f\t%.8f\t%.8f\t%.8f\n",t,v2t,becc,eecc);*/
    
    //fprintf(fp, "%.3f\t%.8f\n",t,phiQ);
    //if(t>20)
    //exit(0);
    
    //fclose(fp);
}

// To test the interpolatin function to see it reproduce the EOS table
void testEOS(){
    
    //char EOStable[] = "output/sigmaB_test.dat";
    //ofstream eos_table(EOStable);
    char EOStable1[] = "output/EOS_t_cpuvh.dat";
    ofstream eos_table1(EOStable1);
    //char EOStable2[] = "output/eos_dpdrhob_test5.dat";
    //ofstream eos_table2(EOStable2);
    //char EOStable3[] = "output/eos_p_test5.dat";
    //ofstream eos_table3(EOStable3);
    //char EOStable4[] = "output/eos_cs2_test5.dat";
    //ofstream eos_table4(EOStable4);
    
    for(int i = 0; i < 41; ++i) {
        for(int j = 0; j < 250; ++j){
            //PRECISION ttest=(0.5+j*5)/HBARC/1000.0;
            ///PRECISION mubtest=i*5/HBARC/1000.0;
            PRECISION etest = (0.045+i*0.01)/HBARC;
            PRECISION rhobtest = (0.0+j*0.006)/HBARC;
            //printf("ttest=%lf, mubtest=%lf.\n",ttest,mubtest);
            //eos_table4  << setprecision(6) << setw(18) << etest*HBARC << setprecision(6) << setw(18) << rhobtest
            //           << setprecision(6) << setw(18) << speedOfSoundSquared(etest, rhobtest) << endl;
            //eos_table3 << setprecision(6) << setw(18) << etest*HBARC << setprecision(6) << setw(18) << rhobtest
            //           << setprecision(6) << setw(18) << equilibriumPressure(etest, rhobtest)*HBARC << endl;
            //eos_table2 << setprecision(6) << setw(18) << etest*HBARC << setprecision(6) << setw(18) << rhobtest
            //           << setprecision(6) << setw(18) << dPdRhob(etest,rhobtest)*HBARC << endl;
            eos_table1 << setprecision(6) << setw(18) << etest*HBARC << setprecision(6) << setw(18) << rhobtest
                       << setprecision(6) << setw(18) << chemicalPotentialOverT(etest, rhobtest) << endl;
            //eos_table  << setprecision(6) << setw(18) << ttest << setprecision(6) << setw(18) << mubtest
            //<< setprecision(6) << setw(18) << baryonDiffusionConstant(ttest, mubtest) << endl;
        }
    }
    
    //eos_table.close();
    eos_table1.close();
    //eos_table2.close();
    //eos_table3.close();
    //eos_table4.close();
    printf("EOS table is reproduced.\n");
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

void testJetCoeff(){
    
    char EOStable1[] = "output/jetCoefficients.dat";
    ofstream eos_table1(EOStable1);
    
    PRECISION qhat;
    
    for(int i = 0; i < 75; ++i) {
        for(int j = 0; j < 90; ++j){
            
            PRECISION Ttest = (90.0 + i * 5.0)*0.001/HBARC;
            PRECISION muBtest = (0.0 + j * 5.0)*0.001/HBARC;
            
            qhat = qHat(Ttest, muBtest);

            eos_table1
            << setprecision(6) << setw(18) << Ttest*1000*HBARC << setprecision(6) << setw(18) << muBtest*1000*HBARC
            << setprecision(6) << setw(18) << qhat*HBARC*HBARC << endl;

        }
    }
    
    eos_table1.close();
    printf("Jet coeff table is reproduced.\n");
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

