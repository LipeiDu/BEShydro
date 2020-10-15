//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <math.h>
#include <cmath>
#include <iostream>
#include <istream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <string>
#include <iomanip>

#include "../include/DynamicalVariables.h"
#include "../include/EquationOfState.h"

using namespace std;

#define HBARC 0.197326938

/**************************************************************************************************************************************************/
/* Equation of State tables, X(e, rhob)
/* (1) conformal see e.g. 1204.6710
/* (2) nonconformal PRC 98 (2018) 034916
/**************************************************************************************************************************************************/

// Read in the table of the Equation of State
void getEquationOfStateTable(){
#ifdef EOS_with_baryon
    
    printf("The EOS with baryon density is used...\n");
    
    FILE *eosfilet, *eosfilep, *eosfilembt, *eosfileprhob;
    PRECISION energy;
    PRECISION baryon;
    
    // temperature table
#ifndef CONFORMAL_EOS
    eosfilet = fopen ("eos/eos_t.dat","r");
#else
    eosfilet = fopen ("eos/conformal_eos_t.dat","r");
#endif
    if(eosfilet==NULL){
        printf("The EOS file eos_t.dat was not opened...\n");
        exit(-1);
    }
    else
    {
        fseek(eosfilet,0L,SEEK_SET);
        for(int i = 0; i < 188580; ++i){
            fscanf(eosfilet,"%lf %lf %lf", & energy, & baryon, & EOState->Temperature[i]);
        }
    }
    fclose(eosfilet);
    
    // pressure table, for conformal EoS, p = e/3
#ifndef CONFORMAL_EOS
    eosfilep = fopen ("eos/eos_p.dat","r");

    if(eosfilep==NULL){
        printf("The EOS file eos_p.dat was not opened...\n");
        exit(-1);
    }
    else
    {
        fseek(eosfilep,0L,SEEK_SET);
        for(int i = 0; i < 188580; ++i){
            fscanf(eosfilep,"%lf %lf %lf", & energy, & baryon, & EOState->Pressure[i]);
        }
    }
    fclose(eosfilep);
#endif
    
    // baryon chemical potential over temperature table
#ifndef CONFORMAL_EOS
    eosfilembt = fopen ("eos/eos_mbovert.dat","r");
#else
    eosfilembt = fopen ("eos/conformal_eos_mbovert.dat","r");
#endif
    if(eosfilembt==NULL){
        printf("The EOS file eos_mbovert.dat was not opened...\n");
        exit(-1);
    }
    else
    {
        fseek(eosfilembt,0L,SEEK_SET);
        for(int i = 0; i < 188580; ++i){
            fscanf(eosfilembt,"%lf %lf %lf", & energy, & baryon, & EOState->alphab[i]);
        }
    }
    fclose(eosfilembt);
 
    // dPdRhob table, for conformal EoS, dpdrhob = 0
#ifndef CONFORMAL_EOS
    eosfileprhob = fopen ("eos/eos_dpdrhob.dat","r");

    if(eosfileprhob==NULL){
        printf("The EOS file eos_dpdrhob.dat was not opened...\n");
        exit(-1);
    }
    else
    {
        fseek(eosfileprhob,0L,SEEK_SET);
        for(int i = 0; i < 188580; ++i){
            fscanf(eosfileprhob,"%lf %lf %lf", & energy, & baryon, & EOState->dpdrhob[i]);
        }
    }
    fclose(eosfileprhob);
#endif
    
    printf("Equation of State table is read in.\n");
#endif
}

/**************************************************************************************************************************************************/
/* Functions for bilinear interpolation
/**************************************************************************************************************************************************/

// column major index for 2d array
int columnIndex(int i, int j, int nrhob){
    return j + nrhob * i;
}

// Bilinear interplotion function
inline PRECISION biLinearInterpolation(PRECISION x, PRECISION y, PRECISION q11, PRECISION q12, PRECISION q21, PRECISION q22, PRECISION x1, PRECISION x2, PRECISION y1, PRECISION y2){
    PRECISION x2x1, y2y1, x2x, y2y, yy1, xx1;
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x  = x2 - x;
    y2y  = y2 - y;
    yy1  = y - y1;
    xx1  = x - x1;
    return 1.0 / (x2x1 * y2y1) * (q11 * x2x * y2y + q21 * xx1 * y2y + q12 * x2x * yy1 + q22 * xx1 * yy1);
}

// If (e, rhob) lies in the nontivial zone, this function can do 2D interpolation to give the inferred value for (e, rhob)
PRECISION InferredPrimaryVariable(PRECISION e, PRECISION rhob, PRECISION e_start, PRECISION d_e, int nrhob, PRECISION d_rhob, int index_start, const PRECISION * const __restrict__ EOS_Variable){
    
    int m, n;
    int s11, s12, s21, s22;
    PRECISION e1, e2, rhob1, rhob2;
    PRECISION Q11, Q12, Q21, Q22;
    
    m = floor((e-e_start)/d_e);
    n = floor(rhob/d_rhob);
    
    e1 = e_start + d_e * m;
    e2 = e_start + d_e * (m+1);
    
    rhob1 = d_rhob * n;
    rhob2 = d_rhob * (n+1);
    
    s11 = index_start + columnIndex(m,n,nrhob);
    s12 = index_start + columnIndex(m,n+1,nrhob);
    s21 = index_start + columnIndex(m+1,n,nrhob);
    s22 = index_start + columnIndex(m+1,n+1,nrhob);
    
    Q11 = EOS_Variable[s11];
    Q12 = EOS_Variable[s12];
    Q21 = EOS_Variable[s21];
    Q22 = EOS_Variable[s22];
    
    return biLinearInterpolation(e, rhob, Q11, Q12, Q21, Q22, e1, e2, rhob1, rhob2);
}


/**************************************************************************************************************************************************/
/* The following two functions only defined for a specific EoS table, By Lipei Jan 18, 2018
/**************************************************************************************************************************************************/

// (e, rhob) can lie in the triavial zone, the value will be trivially equal to the value at the edge of the nontrivial zone.
PRECISION primaryVariablesEOS(PRECISION e, PRECISION rhob, const PRECISION * const __restrict__ EOS_Variable){
    
    PRECISION e0 = e*HBARC;
    PRECISION rhob0 = rhob;
    
    if((0<=e0) && (e0<0.0036))
    {
        if((0<=rhob0) && (rhob0<=0.00498))//zone 1
            return InferredPrimaryVariable(e0, rhob0, 0.0, 0.0003, 500, 0.00001, 0, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 0.00499, 0.0, 0.0003, 500, 0.00001, 0, EOS_Variable)/HBARC;
    }
    else if((0.0036<=e0) && (e0<0.015))
    {
        if((0<=rhob0) && (rhob0<=0.0149))//zone 2
            return InferredPrimaryVariable(e0, rhob0, 0.0036, 0.0006, 300, 0.00005, 6500, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 0.01495, 0.0036, 0.0006, 300, 0.00005, 6500, EOS_Variable)/HBARC;
    }
    else if((0.015<=e0) && (e0<0.045))
    {
        if((0<=rhob0) && (rhob0<=0.0445))//zone 3
            return InferredPrimaryVariable(e0, rhob0, 0.015, 0.001, 180, 0.00025, 12500, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 0.04475, 0.015, 0.001, 180, 0.00025, 12500, EOS_Variable)/HBARC;
    }
    else if((0.045<=e0) && (e0<0.455))
    {
        if((0<=rhob0) && (rhob0<=0.496))//zone 4
            return InferredPrimaryVariable(e0, rhob0, 0.045, 0.01, 250, 0.002, 18080, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 0.498, 0.045, 0.01, 250, 0.002, 18080, EOS_Variable)/HBARC;
    }
    else if((0.455<=e0) && (e0<20.355))
    {
        if((0<=rhob0) && (rhob0<=3.48))//zone 5
            return InferredPrimaryVariable(e0, rhob0, 0.455, 0.1, 350, 0.01, 28580, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 3.49, 0.455, 0.1, 350, 0.01, 28580, EOS_Variable)/HBARC;
    }
    else if((20.355<=e0)&(e0<219.355))
    {
        if((0<=rhob0) && (rhob0<=12.4))//zone 6
            return InferredPrimaryVariable(e0, rhob0, 20.355, 1, 250, 0.05, 98580, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 12.45, 20.355, 1, 250, 0.05, 98580, EOS_Variable)/HBARC;
    }
    else
    {
        if((0<=rhob0) && (rhob0<=39.6))//zone 7
            return InferredPrimaryVariable(e0, rhob0, 219.355, 10, 200, 0.2, 148580, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 39.8, 219.355, 10, 200, 0.2, 148580, EOS_Variable)/HBARC;
    }
}

PRECISION primaryVariablesConformalEOS(PRECISION e, PRECISION rhob, const PRECISION * const __restrict__ EOS_Variable){
    
    PRECISION e0 = e*HBARC;
    PRECISION rhob0 = rhob;
    
    if((0<=e0) && (e0<0.0036))
    {
        if((0<=rhob0) && (rhob0<=0.0249))//zone 1
            return InferredPrimaryVariable(e0, rhob0, 0.0, 0.0003, 500, 0.00005, 0, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 0.02495, 0.0, 0.0003, 500, 0.00005, 0, EOS_Variable)/HBARC;
    }
    else if((0.0036<=e0) && (e0<0.015))
    {
        if((0<=rhob0) && (rhob0<=0.0740))//zone 2
            return InferredPrimaryVariable(e0, rhob0, 0.0036, 0.0006, 300, 0.00025, 6500, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 0.07475, 0.0036, 0.0006, 300, 0.00025, 6500, EOS_Variable)/HBARC;
    }
    else if((0.015<=e0) && (e0<0.045))
    {
        if((0<=rhob0) && (rhob0<=0.445))//zone 3
            return InferredPrimaryVariable(e0, rhob0, 0.015, 0.001, 180, 0.0025, 12500, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 0.4475, 0.015, 0.001, 180, 0.0025, 12500, EOS_Variable)/HBARC;
    }
    else if((0.045<=e0) && (e0<0.455))
    {
        if((0<=rhob0) && (rhob0<=1.488))//zone 4
            return InferredPrimaryVariable(e0, rhob0, 0.045, 0.01, 250, 0.006, 18080, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 1.494, 0.045, 0.01, 250, 0.006, 18080, EOS_Variable)/HBARC;
    }
    else if((0.455<=e0) && (e0<20.355))
    {
        if((0<=rhob0) && (rhob0<=6.96))//zone 5
            return InferredPrimaryVariable(e0, rhob0, 0.455, 0.1, 350, 0.02, 28580, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 6.98, 0.455, 0.1, 350, 0.02, 28580, EOS_Variable)/HBARC;
    }
    else if((20.355<=e0)&(e0<219.355))
    {
        if((0<=rhob0) && (rhob0<=24.8))//zone 6
            return InferredPrimaryVariable(e0, rhob0, 20.355, 1, 250, 0.1, 98580, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 24.9, 20.355, 1, 250, 0.1, 98580, EOS_Variable)/HBARC;
    }
    else
    {
        if((0<=rhob0) && (rhob0<=39.6))//zone 7
            return InferredPrimaryVariable(e0, rhob0, 219.355, 10, 200, 0.2, 148580, EOS_Variable)/HBARC;
        else
            return InferredPrimaryVariable(e0, 39.8, 219.355, 10, 200, 0.2, 148580, EOS_Variable)/HBARC;
    }
}

/**************************************************************************************************************************************************/
/* Equation of State, functions X(e, rhob)
/* 1. without rhob: (1) conformal (2) Wuppertal-Budapest
/* 2. with rhob: (1) conformal (2) PRC 98 (2018) 034916
/**************************************************************************************************************************************************/

PRECISION equilibriumEntropy(PRECISION e, PRECISION rhob, PRECISION p, PRECISION T, PRECISION alphaB){
    PRECISION s = (e + p) / T - alphaB * rhob;
    if(s<0.0){
        //printf("Warning: s=%f\t e=%f\t rhob=%f, set to be zero!\n",s,e,rhob);
        s = 1.e-10;
    }
    return s;
}

PRECISION chemicalPotentialOverT(PRECISION e, PRECISION rhob){
#ifndef EOS_with_baryon
    return 0.0;
#else
#ifndef CONFORMAL_EOS
    return primaryVariablesEOS(e, rhob, EOState->alphab)*HBARC;
#else
    return primaryVariablesConformalEOS(e, rhob, EOState->alphab)*HBARC;
    
    //return EOS_ALPHA;
#endif
#endif
}

PRECISION dPdRhob(PRECISION e, PRECISION rhob){
#ifndef EOS_with_baryon
    return 0.0;
#else
#ifndef CONFORMAL_EOS
    return primaryVariablesEOS(e, rhob, EOState->dpdrhob);
#else
    return 0.0;
#endif
#endif
}

PRECISION dPdE(PRECISION e, PRECISION rhob) {
#ifndef EOS_with_baryon
#ifndef CONFORMAL_EOS
    return dpdeWB(e);
#else
    return 1/3;
#endif
#else
#ifndef CONFORMAL_EOS
    PRECISION e0 = e*HBARC;
    PRECISION delta_e = 0;
    if((0<=e0) && (e0<0.0036))
        delta_e = 0.0003;
    else if((0.0036<=e0) && (e0<0.015))
        delta_e = 0.0006;
    else if((0.015<=e0) && (e0<0.045))
        delta_e = 0.001;
    else if((0.045<=e0) && (e0<0.455))
        delta_e = 0.01;
    else if((0.455<=e0) && (e0<20.355))
        delta_e = 0.1;
    else if((20.355<=e0)&(e0<219.355))
        delta_e = 1;
    else
        delta_e = 10;
    
    PRECISION ep, em, pp, pm, p0;
    ep = (e0 + 0.1 * delta_e)/HBARC;
    em = (e0 - 0.1 * delta_e)/HBARC;
    
    pp = primaryVariablesEOS(ep, rhob, EOState->Pressure);
    if(em>=0){
        pm = primaryVariablesEOS(em, rhob, EOState->Pressure);
        return (pp-pm)/(2*0.1 * delta_e/HBARC);
    }else{
        p0 = primaryVariablesEOS(e0, rhob, EOState->Pressure);
        return (pp-p0)/(0.1*delta_e)/HBARC;
    }
#else
    return 1/3;
#endif
#endif
}

PRECISION equilibriumPressure(PRECISION e, PRECISION rhob) {
#ifndef EOS_with_baryon
#ifndef CONFORMAL_EOS
    return equilibriumPressureWB(e);
#else
    return e/3;
#endif
#else
#ifndef CONFORMAL_EOS
    return primaryVariablesEOS(e, fabs(rhob), EOState->Pressure);
#else
    return e/3;
#endif
#endif
}

PRECISION effectiveTemperature(PRECISION e, PRECISION rhob) {
#ifndef EOS_with_baryon
#ifndef CONFORMAL_EOS
    return effectiveTemperatureWB(e);
#else
    return powf(e/EOS_FACTOR, 0.25);
#endif
#else
#ifndef CONFORMAL_EOS
    return primaryVariablesEOS(e, fabs(rhob), EOState->Temperature);
#else
    return primaryVariablesConformalEOS(e, fabs(rhob), EOState->Temperature);
    
    //return powf(e/EOS_FACTOR, 0.25);
#endif
#endif
}

PRECISION speedOfSoundSquared(PRECISION e, PRECISION rhob) {
#ifndef EOS_with_baryon
#ifndef CONFORMAL_EOS
    return dpdeWB(e);
#else
    return 1/3;
#endif
#else
#ifndef CONFORMAL_EOS
    PRECISION p = equilibriumPressure(e, rhob);
    PRECISION dp_drhob = dPdRhob(e, rhob);
    PRECISION dp_de = dPdE(e, rhob);
    return dp_de + rhob / (e + p) * dp_drhob;
#else
    return 1/3;
#endif
#endif
}

PRECISION dPdT(PRECISION e, PRECISION rhob) {
    PRECISION ep = 1.1 * e;
    PRECISION em = 0.9 * e;
    
    PRECISION pPrimaryVariables[3], mPrimaryVariables[3];
    
    getPrimaryVariablesCombo(ep, rhob, pPrimaryVariables);
    getPrimaryVariablesCombo(em, rhob, mPrimaryVariables);
    
    PRECISION dp = pPrimaryVariables[0] - mPrimaryVariables[0];
    PRECISION dT = pPrimaryVariables[1] - mPrimaryVariables[1];
    
    return dp/dT;
}

// a function returns temperature, pressure and chemical potential over temperature altogether
void getPrimaryVariablesCombo(PRECISION e, PRECISION rhob, PRECISION * const __restrict__ PrimaryVariables){
#ifndef EOS_with_baryon
#ifndef CONFORMAL_EOS
    PrimaryVariables[0] = equilibriumPressureWB(e);
    PrimaryVariables[1] = effectiveTemperatureWB(e);
    PrimaryVariables[2] = 0.0;
#else
    PrimaryVariables[0] = e/3;
    PrimaryVariables[1] = powf(e/EOS_FACTOR, 0.25);
    PrimaryVariables[2] = 0.0;
#endif
#else
#ifndef CONFORMAL_EOS
    PRECISION e0 = e*HBARC;
    PRECISION rhob0 = rhob;
    
    if((0<=e0) && (e0<0.0036))
    {
        if((0<=rhob0) && (rhob0<=0.00498)){//zone 1
            PrimaryVariables[0] = InferredPrimaryVariable(e0, rhob0, 0.0, 0.0003, 500, 0.00001, 0, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 0.0, 0.0003, 500, 0.00001, 0, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 0.0, 0.0003, 500, 0.00001, 0, EOState->alphab);
        }
        else{
            PrimaryVariables[0] = InferredPrimaryVariable(e0, 0.00499, 0.0, 0.0003, 500, 0.00001, 0, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 0.00499, 0.0, 0.0003, 500, 0.00001, 0, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 0.00499, 0.0, 0.0003, 500, 0.00001, 0, EOState->alphab);
        }
    }
    else if((0.0036<=e0) && (e0<0.015))
    {
        if((0<=rhob0) && (rhob0<=0.0149)){//zone 2
            PrimaryVariables[0] = InferredPrimaryVariable(e0, rhob0, 0.0036, 0.0006, 300, 0.00005, 6500, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 0.0036, 0.0006, 300, 0.00005, 6500, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 0.0036, 0.0006, 300, 0.00005, 6500, EOState->alphab);
        }
        else{
            PrimaryVariables[0] = InferredPrimaryVariable(e0, 0.01495, 0.0036, 0.0006, 300, 0.00005, 6500, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 0.01495, 0.0036, 0.0006, 300, 0.00005, 6500, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 0.01495, 0.0036, 0.0006, 300, 0.00005, 6500, EOState->alphab);
        }
    }
    else if((0.015<=e0) && (e0<0.045))
    {
        if((0<=rhob0) && (rhob0<=0.0445)){//zone 3
            PrimaryVariables[0] = InferredPrimaryVariable(e0, rhob0, 0.015, 0.001, 180, 0.00025, 12500, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 0.015, 0.001, 180, 0.00025, 12500, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 0.015, 0.001, 180, 0.00025, 12500, EOState->alphab);
        }
        else{
            PrimaryVariables[0] = InferredPrimaryVariable(e0, 0.04475, 0.015, 0.001, 180, 0.00025, 12500, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 0.04475, 0.015, 0.001, 180, 0.00025, 12500, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 0.04475, 0.015, 0.001, 180, 0.00025, 12500, EOState->alphab);
        }
    }
    else if((0.045<=e0) && (e0<0.455))
    {
        if((0<=rhob0) && (rhob0<=0.496)){//zone 4
            PrimaryVariables[0] = InferredPrimaryVariable(e0, rhob0, 0.045, 0.01, 250, 0.002, 18080, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 0.045, 0.01, 250, 0.002, 18080, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 0.045, 0.01, 250, 0.002, 18080, EOState->alphab);
        }
        else{
            PrimaryVariables[0] = InferredPrimaryVariable(e0, 0.498, 0.045, 0.01, 250, 0.002, 18080, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 0.498, 0.045, 0.01, 250, 0.002, 18080, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 0.498, 0.045, 0.01, 250, 0.002, 18080, EOState->alphab);
        }
    }
    else if((0.455<=e0) && (e0<20.355))
    {
        if((0<=rhob0) && (rhob0<=3.48)){//zone 5
            PrimaryVariables[0] = InferredPrimaryVariable(e0, rhob0, 0.455, 0.1, 350, 0.01, 28580, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 0.455, 0.1, 350, 0.01, 28580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 0.455, 0.1, 350, 0.01, 28580, EOState->alphab);
        }
        else{
            PrimaryVariables[0] = InferredPrimaryVariable(e0, 3.49, 0.455, 0.1, 350, 0.01, 28580, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 3.49, 0.455, 0.1, 350, 0.01, 28580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 3.49, 0.455, 0.1, 350, 0.01, 28580, EOState->alphab);
        }
    }
    else if((20.355<=e0)&(e0<219.355))
    {
        if((0<=rhob0) && (rhob0<=12.4)){//zone 6
            PrimaryVariables[0] = InferredPrimaryVariable(e0, rhob0, 20.355, 1, 250, 0.05, 98580, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 20.355, 1, 250, 0.05, 98580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 20.355, 1, 250, 0.05, 98580, EOState->alphab);
        }
        else{
            PrimaryVariables[0] = InferredPrimaryVariable(e0, 12.45, 20.355, 1, 250, 0.05, 98580, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 12.45, 20.355, 1, 250, 0.05, 98580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 12.45, 20.355, 1, 250, 0.05, 98580, EOState->alphab);
        }
    }
    else
    {
        if((0<=rhob0) && (rhob0<=39.6)){//zone 7
            PrimaryVariables[0] = InferredPrimaryVariable(e0, rhob0, 219.355, 10, 200, 0.2, 148580, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 219.355, 10, 200, 0.2, 148580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 219.355, 10, 200, 0.2, 148580, EOState->alphab);
        }
        else{
            PrimaryVariables[0] = InferredPrimaryVariable(e0, 39.8, 219.355, 10, 200, 0.2, 148580, EOState->Pressure)/HBARC;
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 39.8, 219.355, 10, 200, 0.2, 148580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 39.8, 219.355, 10, 200, 0.2, 148580, EOState->alphab);
        }
    }
#else
    PRECISION e0 = e*HBARC;
    PRECISION rhob0 = rhob;
    
    PrimaryVariables[0] = e/3.;
    /*PrimaryVariables[1] = powf(e/EOS_FACTOR, 0.25);
    PrimaryVariables[2] = EOS_ALPHA;*/
    
    if((0<=e0) && (e0<0.0036))
    {
        if((0<=rhob0) && (rhob0<=0.0249)){//zone 1
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 0.0, 0.0003, 500, 0.00005, 0, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 0.0, 0.0003, 500, 0.00005, 0, EOState->alphab);
        }
        else{
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 0.02495, 0.0, 0.0003, 500, 0.00005, 0, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 0.02495, 0.0, 0.0003, 500, 0.00005, 0, EOState->alphab);
        }
    }
    else if((0.0036<=e0) && (e0<0.015))
    {
        if((0<=rhob0) && (rhob0<=0.0740)){//zone 2
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 0.0036, 0.0006, 300, 0.00025, 6500, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 0.0036, 0.0006, 300, 0.00025, 6500, EOState->alphab);
        }
        else{
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 0.07475, 0.0036, 0.0006, 300, 0.00025, 6500, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 0.07475, 0.0036, 0.0006, 300, 0.00025, 6500, EOState->alphab);
        }
    }
    else if((0.015<=e0) && (e0<0.045))
    {
        if((0<=rhob0) && (rhob0<=0.445)){//zone 3
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 0.015, 0.001, 180, 0.0025, 12500, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 0.015, 0.001, 180, 0.0025, 12500, EOState->alphab);
        }
        else{
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 0.4475, 0.015, 0.001, 180, 0.0025, 12500, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 0.4475, 0.015, 0.001, 180, 0.0025, 12500, EOState->alphab);
        }
    }
    else if((0.045<=e0) && (e0<0.455))
    {
        if((0<=rhob0) && (rhob0<=1.488)){//zone 4
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 0.045, 0.01, 250, 0.006, 18080, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 0.045, 0.01, 250, 0.006, 18080, EOState->alphab);
        }
        else{
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 1.494, 0.045, 0.01, 250, 0.006, 18080, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 1.494, 0.045, 0.01, 250, 0.006, 18080, EOState->alphab);
        }
    }
    else if((0.455<=e0) && (e0<20.355))
    {
        if((0<=rhob0) && (rhob0<=6.96)){//zone 5
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 0.455, 0.1, 350, 0.02, 28580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 0.455, 0.1, 350, 0.02, 28580, EOState->alphab);
        }
        else{
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 6.98, 0.455, 0.1, 350, 0.02, 28580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 6.98, 0.455, 0.1, 350, 0.02, 28580, EOState->alphab);
        }
    }
    else if((20.355<=e0)&(e0<219.355))
    {
        if((0<=rhob0) && (rhob0<=24.8)){//zone 6
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 20.355, 1, 250, 0.1, 98580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 20.355, 1, 250, 0.1, 98580, EOState->alphab);
        }
        else{
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 24.9, 20.355, 1, 250, 0.1, 98580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 24.9, 20.355, 1, 250, 0.1, 98580, EOState->alphab);
        }
    }
    else
    {
        if((0<=rhob0) && (rhob0<=39.6)){//zone 7
            PrimaryVariables[1] = InferredPrimaryVariable(e0, rhob0, 219.355, 10, 200, 0.2, 148580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, rhob0, 219.355, 10, 200, 0.2, 148580, EOState->alphab);
        }
        else{
            PrimaryVariables[1] = InferredPrimaryVariable(e0, 39.8, 219.355, 10, 200, 0.2, 148580, EOState->Temperature)/HBARC;
            PrimaryVariables[2] = InferredPrimaryVariable(e0, 39.8, 219.355, 10, 200, 0.2, 148580, EOState->alphab);
        }
    }
#endif
#endif
}

/**************************************************************************************************************************************************/
/* Equation of state from the Wuppertal-Budapest collaboration, functions X(e)
/**************************************************************************************************************************************************/

PRECISION equilibriumPressureWB(PRECISION e) {

    double e1 = (double)e;
    double e2 = e*e;
    double e3 = e2*e;
    double e4 = e3*e;
    double e5 = e4*e;
    double e6 = e5*e;
    double e7 = e6*e;
    double e8 = e7*e;
    double e9 = e8*e;
    double e10 = e9*e;
    double e11 = e10*e;
    double e12 = e11*e;
    
    double a0 = -0.25181736420168666;
    double a1 = 9737.845799644809;
    double a2 = 1.077580993288114e6;
    double a3 = 3.1729694865420084e6;
    double a4 = 1.6357487344679043e6;
    double a5 = 334334.4309240126;
    double a6 = 41913.439282708554;
    double a7 = 6340.448389300905;
    double a8 = 141.5073484468774;
    double a9 = 0.7158279081255019;
    double a10 = 0.0009417586777847889;
    double a11 = 3.1188455176941583e-7;
    double a12 = 1.9531729608963267e-11;
    PRECISION a = (PRECISION)fma(a12,e12,fma(a11,e11,fma(a10,e10,fma(a9,e9,fma(a8,e8,fma(a7,e7,fma(a6,e6,fma(a5,e5,fma(a4,e4,fma(a3,e3,fma(a2,e2,fma(a1,e1,a0))))))))))));
    
    double b0 = 45829.44617893836;
    double b1 = 4.0574329080826794e6;
    double b2 = 2.0931169138134286e7;
    double b3 = 1.3512402226067686e7;
    double b4 = 1.7851642641834426e6;
    double b5 = 278581.2989342773;
    double b6 = 26452.34905933697;
    double b7 = 499.04919730607065;
    double b8 = 2.3405487982094204;
    double b9 = 0.002962497695527404;
    double b10 = 9.601103399348206e-7;
    double b11 = 5.928138360995685e-11;
    double b12 = 3.2581066229887368e-18;
    PRECISION b = (PRECISION)fma(b12,e12,fma(b11,e11,fma(b10,e10,fma(b9,e9,fma(b8,e8,fma(b7,e7,fma(b6,e6,fma(b5,e5,fma(b4,e4,fma(b3,e3,fma(b2,e2,fma(b1,e1,b0))))))))))));
    
    return a/b;

}

PRECISION effectiveTemperatureWB(PRECISION e) {

    double e1 = (double) e;
    double e2 = e * e1;
    double e3 = e2 * e1;
    double e4 = e3 * e1;
    double e5 = e4 * e1;
    double e6 = e5 * e1;
    double e7 = e6 * e1;
    double e8 = e7 * e1;
    double e9 = e8 * e1;
    double e10 = e9 * e1;
    double e11 = e10 * e1;
    
    return (1.510073201405604e-29 + 8.014062800678687e-18 * e
            + 2.4954778310451065e-10 * e2 + 0.000063810382643387 * e3
            + 0.4873490574161924 * e4 + 207.48582344326206 * e5
            + 6686.07424325115 * e6 + 14109.766109389702 * e7
            + 1471.6180520527757 * e8 + 14.055788949565482 * e9
            + 0.015421252394182246 * e10 + 1.5780479034557783e-6 * e11)
    / (7.558667139355393e-28 + 1.3686372302041508e-16 * e
       + 2.998130743142826e-9 * e2 + 0.0005036835870305458 * e3
       + 2.316902328874072 * e4 + 578.0778724946719 * e5
       + 11179.193315394154 * e6 + 17965.67607192861 * e7
       + 1051.0730543534657 * e8 + 5.916312075925817 * e9
       + 0.003778342768228011 * e10 + 1.8472801679382593e-7 * e11);

}

PRECISION equilibriumEnergyDensityWB(PRECISION T) {

	double T1 = (double) T;
	double T2 = T1 * T1;
	double T3 = T2 * T1;
	double T4 = T3 * T1;
	double T5 = T4 * T1;
	double T6 = T5 * T1;
	double T7 = T6 * T1;
	double T8 = T7 * T1;
	double T9 = T8 * T1;
	double T10 = T9 * T1;
	double T11 = T10 * T1;
	double T12 = T11 * T1;
	double T13 = T12 * T1;
	double T14 = T13 * T1;
	double T15 = T14 * T1;
	double T16 = T15 * T1;
	double T17 = T16 * T1;
	double T18 = T17 * T1;
	double T19 = T18 * T1;
	double T20 = T19 * T1;
	double T21 = T20 * T1;
	double T22 = T21 * T1;
	double T23 = T22 * T1;
    
	return (-0.011958188410851651 + 119.89423098138208 * T
			- 3156.9475699248055 * T2 + 32732.86844374939 * T3
			- 187899.8994764422 * T4 + 712537.3610845465 * T5
			- 1.557049803609345e6 * T6 + 1.4852519861308339e6 * T7
			+ 532132.6079941876 * T8 - 1.963099445042592e6 * T9
			- 4484.44579242679 * T10 + 1.7984228830058286e6 * T11
			+ 119345.25619517374 * T12 - 1.3499773937058165e6 * T13
			- 207838.4995663606 * T14 + 654970.2138652403 * T15
			- 78643.00334616247 * T16 + 40274.00078068926 * T17
			+ 422619.58977657766 * T18 - 409688.07836393174 * T19
			- 62005.75915066359 * T20 + 46788.14270090656 * T21
			+ 40784.330477857235 * T22 - 12589.47744840392 * T23)
			/ (31630.074365558292 - 127100.88940643385 * T
					+ 173528.1225422275 * T2 - 39403.297956865215 * T3
					- 85582.57873541754 * T4 + 9320.560804233442 * T5
					+ 50882.74198960172 * T6 + 20335.926473421183 * T7
					- 14897.725710713818 * T8 - 23836.484117457 * T9
					- 13726.013896090335 * T10 + 4517.908673107615 * T11
					+ 18056.19917986404 * T12 + 14954.82860467155 * T13
					+ 2569.623976952738 * T14 - 9304.046211514986 * T15
					- 15606.429173842751 * T16 + 8383.710735812094 * T17
					+ 1591.3177623932843 * T18 - 678.748230997762 * T19
					- 33.58687934953277 * T20 + 3.2520554133126285 * T21
					- 0.19647288043440464 * T22 + 0.005443394551264717 * T23);

}

PRECISION dpdeWB(PRECISION e) {

    double e1 = (double) e;
    double e2 = e * e1;
    double e3 = e2 * e1;
    double e4 = e3 * e1;
    double e5 = e4 * e1;
    double e6 = e5 * e1;
    double e7 = e6 * e1;
    double e8 = e7 * e1;
    double e9 = e8 * e1;
    double e10 = e9 * e1;
    double e11 = e10 * e1;
    double e12 = e11 * e1;
    double e13 = e12 * e1;
    
    return (5.191934309650155e-32 + 4.123605749683891e-23 * e
            + 3.1955868410879504e-16 * e2 + 1.4170364808063119e-10 * e3
            + 6.087136671592452e-6 * e4 + 0.02969737949090831 * e5
            + 15.382615282179595 * e6 + 460.6487249985994 * e7
            + 1612.4245252438795 * e8 + 275.0492627924299 * e9
            + 58.60283714484669 * e10 + 6.504847576502024 * e11
            + 0.03009027913262399 * e12 + 8.189430244031285e-6 * e13)
    / (1.4637868900982493e-30 + 6.716598285341542e-22 * e
       + 3.5477700458515908e-15 * e2 + 1.1225580509306008e-9 * e3
       + 0.00003551782901018317 * e4 + 0.13653226327408863 * e5
       + 60.85769171450653 * e6 + 1800.5461219450308 * e7
       + 15190.225535036281 * e8 + 590.2572000057821 * e9
       + 293.99144775704605 * e10 + 21.461303090563028 * e11
       + 0.09301685073435291 * e12 + 0.000024810902623582917 * e13);

}

//****************************************************************************\
//* Parameterization based on the Equation of state from the Wuppertal-Budapest collaboration
//* Tref 1.01355
//* h0 0.1396
//* h1 (-0.1800)
//* h2 0.0350
//* alpha 0.01
//* nf = 2+1+1
//* f0 5.59
//* f1 7.34
//* f2 (-5.60)
//* g1 1.42
//* g2 0.5
//****************************************************************************/
