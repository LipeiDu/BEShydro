//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

// Equation indices from PRD 98 (2018) 036006

#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <istream>
#include <fstream>
#include <stdio.h>
#include <cassert>
#include <string>
#include <iomanip>

#include "../include/InitialConditions.h"
#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroParameters.h"
#include "../include/EquationOfState.h"
#include "../include/HydroPlus.h"

using namespace std;

//#define GaussianXI

#define dQvec 0.7 // difference between the elements of Q vector of slow modes
#define Q0 0.5
#define HBARC 0.197326938

#define xi0 1.0
#define xi02 1.0 // correlation length squared
#define xiw 0.0192 // gaussian xi width

#define Cr 4606.67 // LambdaT = Cr * T^2, a normalization factor to match the setup in the paper, Gamma0*T(t0,r=0)*h^2/(alphaB*g)/2

/**************************************************************************************************************************************************/
/* Conductivity, heat capacity, and some derivatives
/**************************************************************************************************************************************************/

// heat conductivity
PRECISION lambdaT(PRECISION T){
    return Cr * T * T;
}

// heat capacity density
PRECISION Cp(PRECISION s, PRECISION rhob, PRECISION corrL2){
    return (s * s / rhob) * (corrL2 / xi02);
}

// derivative of log(phi0) with respect to energy
PRECISION dlnPhi0de(PRECISION T, PRECISION s, PRECISION dlnXi_de){
    return 2/(s*T) + 2 * dlnXi_de;
}

// derivative of log(phi0) with respect to baryon density
PRECISION dlnPhi0drhob(PRECISION alphaB, PRECISION rhob, PRECISION s, PRECISION dlnXi_drhob){
    return -2 * alphaB / s - 3 / rhob + 2 * dlnXi_drhob;
}

/**************************************************************************************************************************************************/
/* Conductivity, heat capacity, correlation length and relaxation coefficients of different slow modes, and their derivatives
/**************************************************************************************************************************************************/

// universal function
PRECISION f2(PRECISION x){
    return 1.0 / (1.0 + x*x); // Eq. (93)
}

// relaxation coefficents of fluctuations, without (Q*xi)f2(Q*xi), just 2*lambdaT/(Cp*xi^2).
PRECISION relaxationCoefficientPhi(PRECISION rhob, PRECISION s, PRECISION T, PRECISION corrL2)
{
    PRECISION lambdat = lambdaT(T);
    PRECISION cp = Cp(s, rhob, corrL2);
    
    return 2 * lambdat/(cp * corrL2);
}

// relaxation coefficents of fluctuations, only work for f2 defined above.
PRECISION relaxationCoefficientPhiQ(PRECISION gammaPhi, PRECISION corrL2, PRECISION Q)
{
    PRECISION Q2 = Q * Q;
    PRECISION qL2 = corrL2 * Q2;
    PRECISION qL4 = qL2 * qL2;
    
    return gammaPhi * (qL2 + qL4);
}


/**************************************************************************************************************************************************/
/* slow modes out and at equilibrium with different Q
/**************************************************************************************************************************************************/

// slow modes with zero Q
PRECISION equilibriumPhi0(PRECISION rhob, PRECISION s, PRECISION corrL2)
{
    return Cp(s, rhob, corrL2) / (rhob * rhob); // slow modes at equilibrium with Q = 0, Eq.(90)
}

// slow modes with nonzero Q
PRECISION equilibriumPhiQ(PRECISION e, PRECISION rhob, PRECISION T, PRECISION muB, PRECISION s, PRECISION Q)
{
    PRECISION corrL = correlationLength(T, muB);
    PRECISION corrL2 = corrL * corrL;
    PRECISION qL = Q * corrL;
    
    return equilibriumPhi0(rhob, s, corrL2) * f2(qL); // Magnitude of mode Q at Equilibrium, Eq. (89)
}

// initialization of slow modes
void setInitialConditionSlowModes(void * latticeParams, void * hydroParams)
{
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
#ifdef HydroPlus
    printf("HYDRO+ is ON, number of slow modes is %d, Q0 is %f, dQ is %f...\n",NUMBER_SLOW_MODES, Q0, dQvec);
    
    // initialization of the Q vector
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        Qvec[n] = Q0 + n * dQvec;
    }

    // initialization of slow mdoes at/out of equilibrium
    for(int k = 2; k < nz+2; ++k) {
        for(int j = 2; j < ny+2; ++j) {
            for(int i = 2; i < nx+2; ++i) {
                
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                
                PRECISION es = e[s];
                PRECISION rhobs = rhob[s];
                PRECISION Ts = T[s];
                PRECISION muBs = Ts * alphaB[s];
                PRECISION seqs = seq[s];
                
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    
                    PRECISION equiPhiQ = equilibriumPhiQ(es, rhobs, Ts, muBs, seqs, Qvec[n]);
                    
                    eqPhiQ->phiQ[n][s] = equiPhiQ;
                    q->phiQ[n][s] = equiPhiQ;
                    
                }
            }
        }
    }
#endif
}

/**************************************************************************************************************************************************/
/* contributions from slow modes to inferred variables for a specific cell
/**************************************************************************************************************************************************/

// note: the integrands of alpha and beta only work for f2 defined above.
// this function takes e/p/rhob/T/alphaB and slow modes PhiQ/eqPhiQ, then returns variables with contributions from slow modes, including p/T/alphaB
void getPressurePlusFromSlowModes(PRECISION * const __restrict__ deltaVariables, PRECISION * const __restrict__ pPlus, const PRECISION * const __restrict__ equiPhiQ, const PRECISION * const __restrict__ PhiQ, PRECISION e, PRECISION rhob, PRECISION p, PRECISION T, PRECISION alphaB, PRECISION s)
{

    PRECISION muB = alphaB * T;
    PRECISION corrL = correlationLength(T, muB); // correlation length
    PRECISION corrL2 = corrL * corrL;
    
    PRECISION heatC = Cp(s, rhob, corrL2); // heat capacity
    
    // derivaties
    PRECISION dlnXi_de = dlnXide(e, rhob);
    PRECISION dlnXi_drhob = dlnXidrhob(e, rhob);
    PRECISION dlnPhi0_de = dlnPhi0de(T, s, dlnXi_de);
    PRECISION dlnPhi0_drhob = dlnPhi0drhob(alphaB, rhob, s, dlnXi_drhob);
    
    PRECISION entropy = 0.0;
    PRECISION alpha = 0.0;
    PRECISION beta = 0.0;
    
    // dQ/(2*pi)^2
    PRECISION facQ = dQvec/(4 * M_PI * M_PI);
    
    // contributions from slow modes to alpha, beta and entropy
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        
        // ln(Phi/eqPhi) and (Phi/eqPhi-1)
        PRECISION phiRatio = PhiQ[n] / (equiPhiQ[n] + 1e-15);
        PRECISION phiRatioLog = log(phiRatio);
        PRECISION phiRatioOne = phiRatio - 1;

        // (Q*xi)f2(Q*xi), Riemann sum for the intergal below
        PRECISION Q = Qvec[n] + 0.5 * dQvec;
        PRECISION qL = Q * corrL;
        PRECISION qL2= qL * qL;
        PRECISION qLf2 = qL / (1 + qL2);
        
        // Q^2
        PRECISION Q2 = Q * Q;
        // Q^2*(Phi/eqPhi-1)
        PRECISION QphiRatioOne = Q2 * phiRatioOne;
        // Q^2*ln(Phi/eqPhi)
        PRECISION QphiRatioLog = Q2 * log(phiRatio);
        
        // delta beta, Eq.(106)
        PRECISION intBeta = QphiRatioOne * (dlnPhi0_de - 2 * qLf2 * dlnXi_de);
        beta += intBeta;
        
        // delta alpha, Eq.(106)
        PRECISION intAlpha = QphiRatioOne * (dlnPhi0_drhob - 2 * qLf2 * dlnXi_drhob);
        alpha += intAlpha;
        
        // delta entropy, Eq.(85)
        PRECISION intEntropy = QphiRatioLog - QphiRatioOne;
        entropy += intEntropy;
    }
    
    // contributions from slow modes to entropy, inverse temperature, chemical potential over temperature and pressure
    PRECISION deltaS = facQ * entropy;
    PRECISION deltaAlphaB = - facQ * alpha;
    PRECISION deltaBeta = facQ * beta;
    
    // variables(+) with contribution from slow modes
    T = 1 / (1/T + deltaBeta);
    
    PRECISION deltaP = T * (deltaS - (e + p) * deltaBeta + rhob * deltaAlphaB);
    *pPlus = p;// + deltaP;
    
    deltaVariables[0] = deltaAlphaB;
    deltaVariables[1] = deltaBeta;
    deltaVariables[2] = deltaP;
    deltaVariables[3] = deltaS;
}

/**************************************************************************************************************************************************/
/* dlnxide, dlnxidn, correlation length table
/**************************************************************************************************************************************************/

PRECISION dlnXide(PRECISION e0, PRECISION rhob0){
    PRECISION ep, em, xip, xim, Tp, Tm, muBp, muBm;
    PRECISION delta_e = 0.02;
    
    ep = e0 + 0.1 * delta_e;
    em = e0 - 0.1 * delta_e;
    
    if(em >= 0 && ep <= 18.0)
    {
        PRECISION mPrimaryVariables[3], pPrimaryVariables[3];
        
        getPrimaryVariablesCombo(em, rhob0, mPrimaryVariables);
        getPrimaryVariablesCombo(ep, rhob0, pPrimaryVariables);
        
        Tm = mPrimaryVariables[1];
        muBm = Tm * mPrimaryVariables[2];
        xim = correlationLength(Tm, muBm);
        
        Tp = pPrimaryVariables[1];
        muBp = Tp * pPrimaryVariables[2];
        xip = correlationLength(Tp, muBp);
        
        return (log(xip)-log(xim))/(2 * 0.1 * delta_e);
    }
    else
        return 0.0;
}

PRECISION dlnXidrhob(PRECISION e0, PRECISION rhob0){
    PRECISION np, nm, xip, xim, Tp, Tm, muBp, muBm;
    PRECISION delta_n = 0.005;
    
    np = rhob0 + delta_n;
    nm = rhob0 - delta_n;
    
    if(nm >= 0 && np <= 0.9)
    {
        PRECISION mPrimaryVariables[3], pPrimaryVariables[3];
        
        getPrimaryVariablesCombo(e0, nm, mPrimaryVariables);
        getPrimaryVariablesCombo(e0, np, pPrimaryVariables);
        
        Tm = mPrimaryVariables[1];
        muBm = Tm * mPrimaryVariables[2];
        xim = correlationLength(Tm, muBm);
        
        Tp = pPrimaryVariables[1];
        muBp = Tp * pPrimaryVariables[2];
        xip = correlationLength(Tp, muBp);
        
        return (log(xip)-log(xim))/(2 * delta_n);
    }
    else
        return 0.0;
}

// Read in the table of correlation length
void getCorrelationLengthTable(){
#ifdef CRITICAL
    printf("CRITICAL is ON...\n");

    FILE *filexi;
    PRECISION x, y;
    
    filexi = fopen ("input/profiles/ximut.dat","r");
    if(filexi==NULL){
        printf("ximut.dat was not opened...\n");
        exit(-1);
    }
    else
    {
        fseek(filexi,0L,SEEK_SET);
        for(int i = 0; i < 8686; ++i){
            fscanf(filexi,"%lf %lf %lf", & x, & y, & xieq[i]);
        }
    }
    fclose(filexi);
#endif
}

// bilinear interpolation of correlation length table
PRECISION correlationLength(PRECISION T, PRECISION muB){
#ifdef CRITICAL
#ifndef GaussianXI
    PRECISION T0 = T*HBARC;
    PRECISION muB0 = muB*HBARC;
    
    if((0.07<=T0)&&(T0<=0.27)&&(0.11<=muB0)&&(muB0<=0.28))
        return InferredPrimaryVariable(T0, muB0-0.11, 0.07, 0.002, 86, 0.002, 0, xieq);
    else
        return 1.0;
#else
    PRECISION T0 = T*HBARC;
    return xi0 + 2.0 * exp(-(T0-0.16)*(T0-0.16)/(2*xiw*xiw));
    
    //long double A = pow(tanh((T0-0.16)/(0.4*0.16)),2);
    //long double B = 1/81;
    //double xi = 1.0 / pow((A*(1-B) + B),0.25);
    //return xi;
#endif
#else
    return 1.0;
#endif
}

PRECISION corrLen(PRECISION T, PRECISION muB, PRECISION xi_max, PRECISION MU_C){

    // T and mu in GeV
    PRECISION T_GeV = T*HBARC;
    PRECISION muB_GeV = muB*HBARC;
    
    // some parameters in the parametrization
    
//     PRECISION pi = 3.14159265359;
    
//     PRECISION delta_T = 0.02;
//     PRECISION delta_mu = 0.065;
    
//     PRECISION Alpha = 4.6 * (pi / 180); // angle between h and r axes
    
    
    PRECISION T0 = 0.155;
    PRECISION k2 = 0.0149;
    PRECISION Alpha = atan(2*k2*MU_C/T0);
    
    PRECISION T_C = T0 - k2 * T0 * (MU_C/T0)*(MU_C/T0);
    
//     PRECISION A = 0.8;
//     PRECISION fac = pow((sin(Alpha)/sin(0.0524)),1.5);
    
//     PRECISION delta_T = T_C * sin(Alpha)/cos(Alpha) * fac;
//     PRECISION delta_mu = A * T_C * cos(Alpha) * fac;
    
    PRECISION delta_T = 0.021;
    PRECISION delta_mu = 0.066;
    
    PRECISION xi_min = 1.0;

    // critical exponent
    PRECISION nu = 2./3.;
    PRECISION betaDeltaInv = 0.6;

    // rotation by angle Alpha
    PRECISION Tp  = (muB_GeV - MU_C) * sin(Alpha) + (T_GeV - T_C) * cos(Alpha);
    PRECISION mup = (muB_GeV - MU_C) * cos(Alpha) - (T_GeV - T_C) * sin(Alpha);

    // terms in the parameters
    PRECISION xi_ratio = xi_min / xi_max;
    PRECISION xi_ratio_nu = pow(xi_ratio, 2./nu);

    PRECISION mu2 = pow((mup / delta_mu), 2);
    PRECISION T2  = pow(fabs(Tp / delta_T), 2*betaDeltaInv);

    PRECISION Tanh_term = tanh(mu2 + T2);
    PRECISION Br_term = Tanh_term * (1 - xi_ratio_nu) + xi_ratio_nu;

    return xi_min / pow(Br_term, nu/2);
    
}
