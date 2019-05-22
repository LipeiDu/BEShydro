/*
 * EnergyMomentumTensor.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <math.h> // for math functions

#include "edu/osu/rhic/trunk/hydro/EnergyMomentumTensor.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"

#include "edu/osu/rhic/trunk/hydro/FullyDiscreteKurganovTadmorScheme.h" // for const params
#include "edu/osu/rhic/trunk/eos/EquationOfState.h"
 
#define MAX_ITERS 10000000

PRECISION velocityFromConservedVariables(PRECISION ePrev, PRECISION M0, PRECISION M, PRECISION Pi, PRECISION rhobPrev, PRECISION delta_nbt, PRECISION vPrev) {
    PRECISION e0 = ePrev;
    PRECISION rhob0 = rhobPrev;
    PRECISION v0 = vPrev;
    PRECISION alpha = 1;
    
    for(int j = 0; j < MAX_ITERS; ++j) {
        PRECISION p    = equilibriumPressure(e0, rhob0);
        PRECISION cs2  = speedOfSoundSquared(e0, rhob0);
        PRECISION dp_drhob = dPdRhob(e0, rhob0);
        
        PRECISION cst2 = p/e0;
        PRECISION dcst2_de = 1/e0 * cs2 - p/e0/e0;
        PRECISION dcst2_drhob = 1/e0 * dp_drhob;
        PRECISION dcst2_dv = -dcst2_de * M - dcst2_drhob * delta_nbt * v0/sqrt(1-v0*v0);
        
        PRECISION A = M0*(1+cst2)+Pi;
        PRECISION B = A*A-4*cst2*M*M;
        
        PRECISION f = A*v0 + sqrt(B)*v0 - 2*M;
        PRECISION fp = A + v0*M0*dcst2_dv + sqrt(B) + v0*(A*M0*dcst2_dv - 2*M*M*dcst2_dv)/sqrt(B);
        PRECISION v = v0 - alpha*f/fp;

        if(fabs(v - v0) <=  0.001 * fabs(v)) return v;
        
        v0 = v;
        e0 = M0 - v0*M;
        rhob0 = delta_nbt*sqrt(1-v0*v0);
        
        if(j == floor(0.4 * MAX_ITERS)) alpha = 0.9;//In case the root solver oscillate back and forth forever
    }
    
    printf("v0 = Maxiter.\t vPrev=%5e,\t ePrev=%5e,\t M0=%5e,\t M =%5e,\t Pi=%5e,\t rhobPrev=%5e,\t d_nbt=%5e.\n",vPrev,ePrev,M0,M,Pi,rhobPrev,delta_nbt);
    
    return v0;
}

PRECISION utauFromConservedVariables(PRECISION ePrev, PRECISION M0, PRECISION M, PRECISION Pi, PRECISION rhobPrev, PRECISION delta_nbt, PRECISION utPrev) {
    PRECISION e0 = ePrev;
    PRECISION rhob0 = rhobPrev;
    PRECISION u0 = utPrev;
    PRECISION alpha = 1;
    
    for(int j = 0; j < MAX_ITERS; ++j) {
        PRECISION p    = equilibriumPressure(e0, rhob0);
        PRECISION cs2  = speedOfSoundSquared(e0, rhob0);
        PRECISION dp_drhob = dPdRhob(e0, rhob0);
        
        PRECISION cst2 = p/e0;
        PRECISION dcst2_de = 1/e0 * cs2 - p/e0/e0;
        PRECISION dcst2_drhob = 1/e0 * dp_drhob;
        PRECISION groot = sqrt(u0*u0-1);
        
        PRECISION dcst2_du = -dcst2_de * M/(u0*u0*groot) - dcst2_drhob * delta_nbt/(u0*u0);
        
        PRECISION A = M0*(1+cst2)+Pi;
        PRECISION B = A*A-4*cst2*M*M;
        
        PRECISION f = A + sqrt(B) - 2*M*u0/groot;
        PRECISION fp = (M0 + (A*M0-2*M*M)/sqrt(B))*dcst2_du + 2*M/(groot*groot*groot);
        PRECISION u = u0 - alpha*f/fp;
        
        if(fabs(u - u0) <=  0.001 * fabs(u)) return u;
        
        u0 = u;
        e0 = M0 - M*sqrt(1-1/(u0*u0));
        rhob0 = delta_nbt/u0;
        
        if(j == floor(0.4 * MAX_ITERS)) alpha = 0.9;//In case the root solver oscillate back and forth forever
    }
    
    printf("u0 = Maxiter.\t uPrev=%5e,\t ePrev=%5e,\t M0=%5e,\t M =%5e,\t Pi=%5e,\t rhobPrev=%5e,\t d_nbt=%5e.\n",utPrev,ePrev,M0,M,Pi,rhobPrev,delta_nbt);
    
    return u0;
}

PRECISION energyDensityFromConservedVariables(PRECISION ePrev, PRECISION M0, PRECISION M, PRECISION Pi, PRECISION rhobPrev) {
#ifndef CONFORMAL_EOS
	PRECISION e0 = ePrev;	//initial guess for energy density
    PRECISION rhob0 = rhobPrev;
	for(int j = 0; j < MAX_ITERS; ++j) {
		PRECISION p = equilibriumPressure(e0, rhob0);
		PRECISION cs2 = speedOfSoundSquared(e0, rhob0);
		PRECISION cst2 = p/e0;

		PRECISION A = M0*(1-cst2)+Pi;
		PRECISION B = M0*(M0+Pi)-M;
		PRECISION H = sqrtf(fabs(A*A+4*cst2*B));
		PRECISION D = (A-H)/(2*cst2);

		PRECISION f = e0 + D;
		PRECISION fp = 1 - ((cs2 - cst2)*(B + D*H - ((cs2 - cst2)*cst2*D*M0)/e0))/(cst2*e0*H);

		PRECISION e = e0 - f/fp;
		if(fabs(e - e0) <=  0.001 * fabs(e)) return e;
		e0 = e;
	}
	printf("ev0 = Maxi_inter.\t ePrev=%.3f,\t M0=%.3f,\t M=%.3f,\t Pi=%.3f \n",ePrev,M0,M,Pi);
	return e0;
#else
	return fabs(sqrtf(fabs(4 * M0 * M0 - 3 * M)) - M0);
#endif
}

void getInferredVariables(PRECISION t, const PRECISION * const __restrict__ q, PRECISION ePrev,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, PRECISION utPrev,
PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un,
PRECISION rhobPrev, PRECISION * const __restrict__ rhob) {
	PRECISION ttt = q[0];
	PRECISION ttx = q[1];
	PRECISION tty = q[2];
	PRECISION ttn = q[3];
#ifdef PIMUNU
	PRECISION pitt = q[4];
	PRECISION pitx = q[5];
	PRECISION pity = q[6];
	PRECISION pitn = q[7];
#else
	PRECISION pitt = 0;
	PRECISION pitx = 0;
	PRECISION pity = 0;
	PRECISION pitn = 0;
#endif
#ifdef PI
	PRECISION Pi = q[14];
#else
	PRECISION Pi = 0;
#endif
#ifdef NBMU
    PRECISION Nbt = q[NUMBER_CONSERVED_VARIABLES];
#else
    PRECISION Nbt = 0;
#endif
#ifdef VMU
    PRECISION nbt = q[NUMBER_CONSERVED_VARIABLES+1];
    PRECISION nbx = q[NUMBER_CONSERVED_VARIABLES+2];
    PRECISION nby = q[NUMBER_CONSERVED_VARIABLES+3];
    PRECISION nbn = q[NUMBER_CONSERVED_VARIABLES+4];
#else
    PRECISION nbt = 0;
    PRECISION nbx = 0;
    PRECISION nby = 0;
    PRECISION nbn = 0;
#endif
    
	PRECISION M0 = ttt - pitt;
	PRECISION M1 = ttx - pitx;
	PRECISION M2 = tty - pity;
	PRECISION M3 = ttn - pitn;
	PRECISION M  = M1 * M1 + M2 * M2 + t * t * M3 * M3;
    

    //===================================================================
    // When baryon density is not involved in the root solver
    //===================================================================
    
#ifndef RootSolver_with_Baryon
    
#ifdef Pi
	if ((M0 * M0 - M + M0 * Pi) < 0)
		Pi = M / M0 - M0;
#endif

	if (ePrev <= 0.1) {
        *e = M0 - M / M0;
	}
    else {
		*e = energyDensityFromConservedVariables(ePrev, M0, M, Pi, rhobPrev);
    }
    
	if (isnan(*e)) {
		printf("e is nan. \n M0=%.3f,\t M1=%.3f,\t M2=%.3f,\t M3=%.3f\n", M0, M1, M2, M3);
        printf("ttt=%.3f,\t ttx=%.3f,\t tty=%.3f,\t ttn=%.3f\n", ttt, ttx, tty, ttn);
        printf("pitt=%.3f,\t pitx=%.3f,\t pity=%.3f,\t pitn=%.3f\n", pitt, pitx, pity, pitn);
	}
    
    if (isinf(*e)) {
        printf("e is inf. \n M0=%.3f,\t M1=%.3f,\t M2=%.3f,\t M3=%.3f\n", M0, M1, M2, M3);
        printf("ttt=%.3f,\t ttx=%.3f,\t tty=%.3f,\t ttn=%.3f\n", ttt, ttx, tty, ttn);
        printf("pitt=%.3f,\t pitx=%.3f,\t pity=%.3f,\t pitn=%.3f\n", pitt, pitx, pity, pitn);
    }
    
	if (*e < 1.e-7) {
		*e = 1.e-7;
		*p = 1.e-7;
    }else{
        *p = equilibriumPressure(*e);
    }

	PRECISION P = *p + Pi;
	PRECISION E = 1/(*e + P);
	*ut = sqrtf(fabs((M0 + P) * E));
	PRECISION E2 = E/(*ut);
	*ux = M1 * E2;
	*uy = M2 * E2;
	*un = M3 * E2;

    //baryon by Lipei; baryon density is not evolving by default
    *rhob = rhobPrev;
    
    //===================================================================
    // When baryon number is involved in the root solver
    //===================================================================
#else
    
    PRECISION Ms = sqrt(M);

#ifdef Pi
    if(M0 + Pi - Ms < 0){
        Pi = Ms - M0;
    }
#endif
    
    PRECISION v0 = 0;
    PRECISION u0 = 0;
    PRECISION delta_nbt = Nbt - nbt;
    PRECISION vPrev = sqrt(1-1/(utPrev*utPrev));
    
    if(ePrev <= 0.1)
        v0 = Ms/M0;
    else
        v0 = velocityFromConservedVariables(ePrev, M0, Ms, Pi, rhobPrev, delta_nbt, vPrev);
    
    if (isnan(v0)) {
        printf("v0 = nan.\t vPrev=%5e,\t ePrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t rhobPrev=%5e,\t d_nbt=%5e.\n",vPrev,ePrev,M0,Ms,Pi,rhobPrev,delta_nbt);
    }
    
    if(v0<0.563624&&v0>=0){
        *e    = M0 - v0 * Ms;
        *rhob = delta_nbt * sqrt(1 - v0*v0);

        if (*e < 3.e-4)
        {
            *e = 3.e-4;
            *p = 3.e-4;
        }else{
            *p = equilibriumPressure(*e, *rhob);
        }
        
        if (*rhob<1.e-4 && *rhob >=0) *rhob = 1.e-4;
        else if (*rhob<0 && *rhob > -1.e-4) *rhob = -1.e-4;
    
        PRECISION P  = *p + Pi;
        PRECISION v1 = M1/(M0 + P);
        PRECISION v2 = M2/(M0 + P);
        PRECISION v3 = M3/(M0 + P);
        u0  = 1/sqrt(1 - v0*v0);
        *ut = u0;
        *ux = u0 * v1;
        *uy = u0 * v2;
        *un = u0 * v3;
    }
    else{
        if(ePrev <= 0.1)
            u0 = 1/sqrt(1-(Ms/M0)*(Ms/M0));
        else
            u0 = utauFromConservedVariables(ePrev, M0, Ms, Pi, rhobPrev, delta_nbt, utPrev);
        
        if (isnan(u0)) {
            printf("u0 = nan.\t utPrev=%5e,\t ePrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t rhobPrev=%5e,\t d_nbt=%5e.\n",utPrev,ePrev,M0,Ms,Pi,rhobPrev,delta_nbt);
        }
        
        *e    = M0 - Ms * sqrt(1 - 1/(u0*u0));
        *rhob = delta_nbt/u0;
        
        if (*e < 3.e-4)
        {
            *e = 3.e-4;
            *p = 3.e-4;
        }else{
            *p = equilibriumPressure(*e, *rhob);
        }
        
        if (*rhob<1.e-4 && *rhob >=0) *rhob = 1.e-4;
        else if (*rhob<0 && *rhob > -1.e-4) *rhob = -1.e-4;
        
        PRECISION P  = *p + Pi;
        *ut = u0;
        *ux = u0 * M1/(M0 + P);
        *uy = u0 * M2/(M0 + P);
        *un = u0 * M3/(M0 + P);
    }
#endif
}


void setInferredVariablesKernel(const CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, const FLUID_VELOCITY * const __restrict__ uPrev, FLUID_VELOCITY * const __restrict__ u, PRECISION t, void * latticeParams, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ muB, PRECISION * const __restrict__ T) {
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    
    int ncx,ncy,ncz;
    ncx = lattice->numComputationalLatticePointsX;
    ncy = lattice->numComputationalLatticePointsY;
    ncz = lattice->numComputationalLatticePointsRapidity;
    
    for(int i = 2; i < ncx-2; ++i) {
        for(int j = 2; j < ncy-2; ++j) {
            for(int k = 2; k < ncz-2; ++k) {
                
                int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                
                PRECISION q_s[ALL_NUMBER_CONSERVED_VARIABLES],_e,_p,ut,ux,uy,un,_rhob;//rhob by Lipei
                
                q_s[0] = q->ttt[s];
                q_s[1] = q->ttx[s];
                q_s[2] = q->tty[s];
                q_s[3] = q->ttn[s];
#ifdef PIMUNU
                q_s[4] = q->pitt[s];
                q_s[5] = q->pitx[s];
                q_s[6] = q->pity[s];
                q_s[7] = q->pitn[s];
                q_s[8] = q->pixx[s];
                q_s[9] = q->pixy[s];
                q_s[10] = q->pixn[s];
                q_s[11] = q->piyy[s];
                q_s[12] = q->piyn[s];
                q_s[13] = q->pinn[s];
#endif
#ifdef PI
                q_s[14] = q->Pi[s];
#endif
#ifdef NBMU
                q_s[NUMBER_CONSERVED_VARIABLES] = q->Nbt[s];
#endif
#ifdef VMU
                q_s[NUMBER_CONSERVED_VARIABLES+1] = q->nbt[s];
                q_s[NUMBER_CONSERVED_VARIABLES+2] = q->nbx[s];
                q_s[NUMBER_CONSERVED_VARIABLES+3] = q->nby[s];
                q_s[NUMBER_CONSERVED_VARIABLES+4] = q->nbn[s];
#endif
            
                PRECISION ut_s = uPrev->ut[s];
                
                getInferredVariables(t,q_s,e[s],&_e,&_p,ut_s,&ut,&ux,&uy,&un,rhob[s],&_rhob);
                
                e[s] = _e;
                p[s] = _p;
                rhob[s]  = _rhob;
                u->ut[s] = ut;
                u->ux[s] = ux;
                u->uy[s] = uy;
                u->un[s] = un;
                
                T[s] = effectiveTemperature(_e, _rhob);
                if (T[s] < 1.e-7) T[s] = 1.e-7;
                
#ifdef NBMU
                muB[s] = chemicalPotentialOverT(_e, _rhob);
                if (muB[s]>=0 && muB[s] < 1.e-7) muB[s] = 1.e-7;
                else if (muB[s]<=0 && muB[s] > -1.e-7)  muB[s] = -1.e-7;

                term2[s] = muB[s]*T[s];
#endif
            }
        }
    }
}


//===================================================================
// Components of T^{\mu\nu} and Nb^{\mu} in (\tau,x,y,\eta_s)-coordinates
//===================================================================
PRECISION Ttt(PRECISION e, PRECISION p, PRECISION ut, PRECISION pitt) {
	return (e+p)*ut*ut-p+pitt;
}

PRECISION Ttx(PRECISION e, PRECISION p, PRECISION ut, PRECISION ux, PRECISION pitx) {
	return (e+p)*ut*ux+pitx;
}
 
PRECISION Tty(PRECISION e, PRECISION p, PRECISION ut, PRECISION uy, PRECISION pity) {
	return (e+p)*ut*uy+pity;
}

PRECISION Ttn(PRECISION e, PRECISION p, PRECISION ut, PRECISION un, PRECISION pitn) {
	return (e+p)*ut*un+pitn;
}
 
PRECISION Txx(PRECISION e, PRECISION p, PRECISION ux, PRECISION pixx) {
	return (e+p)*ux*ux+p+pixx;
}

PRECISION Txy(PRECISION e, PRECISION p, PRECISION ux, PRECISION uy, PRECISION pixy) {
	return (e+p)*ux*uy+pixy;
}

PRECISION Txn(PRECISION e, PRECISION p, PRECISION ux, PRECISION un, PRECISION pixn) {
	return (e+p)*ux*un+pixn;
}

PRECISION Tyy(PRECISION e, PRECISION p, PRECISION uy, PRECISION piyy) {
	return (e+p)*uy*uy+p+piyy;
}

PRECISION Tyn(PRECISION e, PRECISION p, PRECISION uy, PRECISION un, PRECISION piyn) {
	return (e+p)*uy*un+piyn;
}
 
PRECISION Tnn(PRECISION e, PRECISION p, PRECISION un, PRECISION pinn, PRECISION t) {
	return (e+p)*un*un+p/t/t+pinn;
}

PRECISION Nbt(PRECISION rhob, PRECISION ut, PRECISION nbt){
    return rhob*ut+nbt;
}

