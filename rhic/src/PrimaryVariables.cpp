//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <math.h>
#include <stdlib.h>


#include "../include/PrimaryVariables.h"
#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/FullyDiscreteKurganovTadmorScheme.h"
#include "../include/EquationOfState.h"
#include "../include/HydroPlus.h"

#define MAX_ITERS 1000

#define Method_Newton

/**************************************************************************************************************************************************/
/* calculate v, u^tau or energy density, Newton method
/**************************************************************************************************************************************************/

// designed for the case without baryon evolution
PRECISION energyDensityFromConservedVariables(PRECISION ePrev, PRECISION M0, PRECISION M, PRECISION Pi, PRECISION rhobPrev) {
#ifndef CONFORMAL_EOS
    PRECISION e0 = ePrev;    //initial guess for energy density
    PRECISION rhob0 = rhobPrev;
    for(int j = 0; j < MAX_ITERS; ++j) {
        PRECISION p = equilibriumPressure(e0, rhobPrev);
        PRECISION cs2 = dPdE(e0, rhobPrev);
        PRECISION cst2 = p/e0;
        
        PRECISION A = M0*(1-cst2)+Pi;
        PRECISION B = M0*(M0+Pi)-M;
        PRECISION H = sqrtf(A*A+4*cst2*B);
        PRECISION D = (A-H)/(2*cst2);
        
        PRECISION f = e0 + D;
        PRECISION fp = 1 - ((cs2 - cst2)*(B + D*H - ((cs2 - cst2)*cst2*D*M0)/e0))/(cst2*e0*H);
        
        PRECISION e = e0 - f/fp;
        if(fabs(e - e0) <=  0.001 * fabs(e)) return e;
        e0 = e;
    }
    printf("ev0 = Maxiter.\t ePrev=%.3f,\t M0=%.3f,\t M=%.3f,\t Pi=%.3f \n",ePrev,M0,M,Pi);
    return e0;
#else
    return fabs(sqrtf(4 * M0 * M0 - 3 * M) - M0);
#endif
}

PRECISION energyDensityFromConservedVariables2(PRECISION ePrev, PRECISION M0, PRECISION M, PRECISION Pi, PRECISION rhobPrev) {
#ifndef CONFORMAL_EOS
    PRECISION e0 = ePrev;    //initial guess for energy density
    PRECISION rhob0 = rhobPrev;
    for(int j = 0; j < MAX_ITERS; ++j) {
        PRECISION p = equilibriumPressure(e0, rhobPrev);
        PRECISION cs2 = dPdE(e0, rhobPrev);
        
        PRECISION A = 1 / (M0 + p + Pi);
        PRECISION A2 = A * A;
        PRECISION B = M * M;
        
        PRECISION f = e0 - M0 + B * A;
        PRECISION fp = 1 - B * A2 * cs2;
        
        PRECISION e = e0 - f/fp;
        if(fabs(e - e0) <=  0.001 * fabs(e)) return e;
        e0 = e;
    }
    printf("ev0 = Maxiter.\t ePrev=%.3f,\t M0=%.3f,\t M=%.3f,\t Pi=%.3f \n",ePrev,M0,M,Pi);
    return e0;
#else
    return fabs(sqrtf(4 * M0 * M0 - 3 * M) - M0);
#endif
}

// designed for the case with baryon evolution, but no slow modes
PRECISION velocityFromConservedVariables(PRECISION M0, PRECISION Ms, PRECISION Pi, PRECISION N0, PRECISION vPrev) {

    PRECISION v0 = vPrev;
    PRECISION e0 = M0 - v0*Ms;
    PRECISION rhob0 = N0*sqrt(1-v0*v0);
    PRECISION alpha = 1;
    
    for(int j = 0; j < MAX_ITERS; ++j) {
        PRECISION p    = equilibriumPressure(e0, rhob0);
        PRECISION dp_de  = dPdE(e0, rhob0);
        PRECISION dp_drhob = dPdRhob(e0, rhob0);
        
        PRECISION f = v0 - Ms / (M0 + p + Pi);
        PRECISION fp = 1 - Ms /(M0 + p + Pi)/(M0 + p + Pi) * ( Ms * dp_de + N0 * dp_drhob * v0 / sqrt(1-v0*v0));
        PRECISION v = v0 - alpha*f/fp;

        if(fabs(v - v0) <=  0.001 * fabs(v)) return v;
        
        v0 = v;
        e0 = M0 - v0*Ms;
        rhob0 = N0*sqrt(1-v0*v0);
        
        if(j == floor(0.4 * MAX_ITERS)) alpha = 0.9;//In case the root solver oscillates back and forth forever
    }
    
    printf("v0 = Maxiter.\t vPrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t N0=%5e.\n",vPrev,M0,Ms,Pi,N0);
    
    return v0;
}

// designed for the case with baryon evolution, but no slow modes
PRECISION utauFromConservedVariables(PRECISION M0, PRECISION Ms, PRECISION Pi, PRECISION N0, PRECISION utPrev) {

    PRECISION u0 = utPrev;
    PRECISION e0 = M0 - Ms * sqrt(1 - 1 / (u0 * u0));
    PRECISION rhob0 = N0 / u0;
    PRECISION alpha = 1;
    
    for(int j = 0; j < MAX_ITERS; ++j) {
        PRECISION p    = equilibriumPressure(e0, rhob0);
        PRECISION dp_de  = dPdE(e0, rhob0);
        PRECISION dp_drhob = dPdRhob(e0, rhob0);
        
        PRECISION u02 = u0 * u0;
        PRECISION v0 = sqrt(1 - 1 / u02);
        PRECISION de_du0 = - Ms / (u02 * u0 * v0 + 1.e-10);
        PRECISION drhob_du0 = - N0 / (u02 + 1.e-10);
        PRECISION dp_du0 = dp_drhob * drhob_du0 + dp_de * de_du0;
        
        PRECISION epsum = e0 + p + Pi;
        PRECISION mpsum = M0 + p + Pi;
        
        PRECISION f = u0 - sqrt(mpsum / epsum);
        PRECISION fp  = 1 - 0.5 * dp_du0 * (1 - mpsum / epsum) / sqrt(mpsum * epsum);
        
        PRECISION u = u0 - alpha*f/fp;
        
        if(fabs(u - u0) <=  0.001 * fabs(u)) return u;
        
        u0 = u;
        e0 = M0 - Ms * sqrt(1 - 1 / (u0 * u0));
        rhob0 = N0 / u0;
        
        if(j == floor(0.4 * MAX_ITERS)) alpha = 0.9;//In case the root solver oscillates back and forth forever
    }
    
    printf("u0 = Maxiter.\t vPrev=%5e,\t M0=%5e,\t Ms =%5e,\t Pi=%5e,\t N0=%5e.\n",utPrev,M0,Ms,Pi,N0);
    
    return u0;
}

/**************************************************************************************************************************************************/
/* calculate v, u^tau, non-Newton method, without slow modes
/**************************************************************************************************************************************************/

// designed for the case with baryon evolution, but no slow modes, not Newton method
PRECISION velocityIterationFromConservedVariables(PRECISION M0, PRECISION Ms, PRECISION Pi, PRECISION N0, PRECISION vPrev) {
    
    PRECISION vRelativeError = 10.0;
    
    PRECISION v0 = vPrev;
    
    PRECISION e0 = M0 - vPrev * Ms;
    PRECISION rhob0 = N0 * sqrt(1 - vPrev * vPrev);
    
    for(int j = 0; j < MAX_ITERS; ++j) {
        
        PRECISION p0 = equilibriumPressure(e0, rhob0);
        
        v0 = Ms/(M0 + p0 + Pi);
        
        vRelativeError = fabs(v0 - vPrev)/(v0 + 1.e-10);
        
        if (vRelativeError < 1.e-6) return v0;
        
        vPrev = v0;
        e0 = M0 - vPrev * Ms;
        rhob0 = N0 * sqrt(1 - vPrev * vPrev);
    }
    
    printf("v0 = Maxiter.\t vPrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t N0=%5e.\n",vPrev,M0,Ms,Pi,N0);
    
    return v0;
}

// designed for the case with baryon evolution, but no slow modes, not Newton method
PRECISION utauIterationFromConservedVariables(PRECISION M0, PRECISION Ms, PRECISION Pi, PRECISION N0, PRECISION utPrev) {
    
    PRECISION utRelativeError = 10.0;
    
    PRECISION ut0 = utPrev;
    
    PRECISION e0 = M0 - Ms * sqrt(1 - 1 / (utPrev * utPrev));
    PRECISION rhob0 = N0 / utPrev;
    
    for(int j = 0; j < MAX_ITERS; ++j) {
        
        PRECISION p0 = equilibriumPressure(e0, rhob0);
        
        ut0 = sqrt((M0 + p0 + Pi)/(e0 + p0 + Pi));
        
        utRelativeError = fabs(ut0 - utPrev)/(ut0 + 1.e-10);
        
        if (utRelativeError < 1.e-4) return ut0;
        
        utPrev = ut0;
        e0 = M0 - Ms * sqrt(1 - 1 / (utPrev * utPrev));
        rhob0 = N0 / utPrev;
    }
    
    printf("u0 = Maxiter.\t uPrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t N0=%5e.\n",utPrev,M0,Ms,Pi,N0);
    
    return ut0;
}

/**************************************************************************************************************************************************/
/* calculate v, u^tau, non Newton, with slow modes
/**************************************************************************************************************************************************/

// designed for the case with baryon evolution and slow modes
void InferredVariablesVelocityIterationHydroPlus(PRECISION * const __restrict__ e, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ p, PRECISION * const __restrict__ T, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ seq, PRECISION * const __restrict__ v, PRECISION * const __restrict__ equiPhiQ, PRECISION M0, PRECISION Ms, PRECISION Pi, PRECISION N0, PRECISION vPrev, const PRECISION * const __restrict__ PhiQ)
{
    PRECISION vRelativeError = 10.0;
    
    PRECISION e0, rhob0, pPlus, peq, Teq, alphaBeq, muBeq, Seq, v0, equiPhiQ0[NUMBER_SLOW_MODES];
    PRECISION deltaVariables[4];
    
    for(int j = 0; j < MAX_ITERS; ++j)
    {
       
        e0 = M0 - vPrev * Ms;
        rhob0 = N0 * sqrt(1 - vPrev * vPrev);
        
        // without contributions from the slow modes, from EOS
        PRECISION PrimaryVariables[3];
        
        getPrimaryVariablesCombo(e0, rhob0, PrimaryVariables);
        
        peq = PrimaryVariables[0];
        Teq = PrimaryVariables[1];
        alphaBeq = PrimaryVariables[2];
        muBeq = Teq * alphaBeq;
        Seq = equilibriumEntropy(e0, rhob0, peq, Teq, alphaBeq);
        
        // include contributions from slow modes into p, T and alphaB
        for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
        {
            equiPhiQ0[n] = equilibriumPhiQ(e0, rhob0, Teq, muBeq, Seq, Qvec[n]);
        }
        
        pPlus = peq; // initialize p(+) as equilibrium value
        
        //getPressurePlusFromSlowModes(deltaVariables, &pPlus, equiPhiQ0, PhiQ, e0, rhob0, peq, Teq, alphaBeq, Seq);
        
        // recalculate the flow velocity, with p(+) from Hydro+
        v0 = Ms/(M0 + pPlus + Pi);
        
        vRelativeError = fabs(v0 - vPrev)/(v0 + 1.e-10);
        
        if (vRelativeError < 1.e-6)
            break;
        
        vPrev = v0;
        
        if(j == MAX_ITERS-1) printf("v = Maxiter.\t vPrev=%5e,\t M0=%5e,\t Ms =%5e,\t Pi=%5e,\t N0=%5e.\n",vPrev,M0,Ms,Pi,N0);
    }
    
    // after iteration, update v, e, rhob, p, T, alphaB and eqPhiQ
    *e = e0;
    *rhob = rhob0;
    *v = v0;
    
    *T = Teq;
    *alphaB = alphaBeq;
    *seq = Seq;
    
    *p = pPlus;

    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
    {
        equiPhiQ[n] = equiPhiQ0[n];
    }
}

// designed for the case with baryon evolution and slow modes
void InferredVariablesUtauIterationHydroPlus(PRECISION * const __restrict__ e, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ p, PRECISION * const __restrict__ T, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ seq, PRECISION * const __restrict__ ut, PRECISION * const __restrict__ equiPhiQ, PRECISION M0, PRECISION Ms, PRECISION Pi, PRECISION N0, PRECISION utPrev, const PRECISION * const __restrict__ PhiQ)
{

    PRECISION utRelativeError = 10.0;
    
    PRECISION e0, rhob0, pPlus, peq, Teq, alphaBeq, muBeq, Seq, ut0, equiPhiQ0[NUMBER_SLOW_MODES];
    
    for(int j = 0; j < MAX_ITERS; ++j)
    {
        e0 = M0 - Ms * sqrt(1 - 1 / (utPrev * utPrev));
        rhob0 = N0/utPrev;
        
        // without contributions from the slow modes, from EOS
        PRECISION PrimaryVariables[3];
        
        getPrimaryVariablesCombo(e0, rhob0, PrimaryVariables);
        
        peq = PrimaryVariables[0];
        Teq = PrimaryVariables[1];
        alphaBeq = PrimaryVariables[2];
        muBeq = Teq * alphaBeq;
        Seq = equilibriumEntropy(e0, rhob0, peq, Teq, alphaBeq);
        
        // include contributions from slow modes into p, T and alphaB
        for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
        {
            equiPhiQ0[n] = equilibriumPhiQ(e0, rhob0, Teq, muBeq, Seq, Qvec[n]);
        }
        
        pPlus = peq; // initialize p(+) as equilibrium value
        
        //getPressurePlusFromSlowModes(&pPlus, equiPhiQ0, PhiQ, e0, rhob0, peq, Teq, alphaBeq, Seq);
        
        // recalculate the flow velocity, with p(+) from Hydro+
        ut0 = sqrt((M0 + pPlus + Pi)/(e0 + pPlus + Pi));
        
        utRelativeError = fabs(ut0 - utPrev)/(ut0 + 1.e-10);
        
        if (utRelativeError < 1.e-4)
            break;
        
        utPrev = ut0;
        
        if(j == MAX_ITERS-1) printf("ut = Maxiter.\t vPrev=%5e,\t M0=%5e,\t Ms =%5e,\t Pi=%5e,\t N0=%5e.\n",utPrev,M0,Ms,Pi,N0);
    }
    
    
    // after iteration, update ut, e, rhob, p, T, alphaB and eqPhiQ
    *e = e0;
    *rhob = rhob0;
    *ut = ut0;
    
    *T = Teq;
    *alphaB = alphaBeq;
    *seq = Seq;
    
    *p = pPlus;
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
    {
        equiPhiQ[n] = equiPhiQ0[n];
    }
}

/**************************************************************************************************************************************************/
/* calculate primary variables u^\mu, e, rhob and pressure etc., from T^\tau^mu N^\tau and dissipative components
/**************************************************************************************************************************************************/

void getInferredVariables(PRECISION t, const PRECISION * const __restrict__ q, PRECISION ePrev, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, PRECISION utPrev, PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un, PRECISION rhobPrev, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ T, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ seq, PRECISION * const __restrict__ equiPhiQ) {
    
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
	PRECISION Pi = q[NUMBER_CONSERVED_VARIABLES_NO_BULK];
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
#ifdef HydroPlus
    PRECISION PhiQ[NUMBER_SLOW_MODES];
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
    {
        PhiQ[n]= q[ALL_NUMBER_CONSERVED_VARIABLES+n];
    }
#else
    PRECISION PhiQ[NUMBER_SLOW_MODES];
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
    {
        PhiQ[n]= 0;
    }
#endif
    
	PRECISION M0 = ttt - pitt;
	PRECISION M1 = ttx - pitx;
	PRECISION M2 = tty - pity;
	PRECISION M3 = ttn - pitn;
    PRECISION N0 = Nbt - nbt;
	PRECISION M  = M1 * M1 + M2 * M2 + t * t * M3 * M3;
    PRECISION Ms = sqrt(M);

    //**************************************************************************
    // SECTION I: baryon evolution is off
    //**************************************************************************
    
#ifndef RootSolver_with_Baryon // Root Solver without Baryon
    
    *e = energyDensityFromConservedVariables2(ePrev, M0, M, Pi, rhobPrev);

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
    
    *p = equilibriumPressure(*e, *rhob);

	PRECISION P = *p + Pi;
	PRECISION E = 1/(*e + P);
    
	*ut = sqrtf((M0 + P) * E);
	PRECISION E2 = E/(*ut);
	*ux = M1 * E2;
	*uy = M2 * E2;
	*un = M3 * E2;

    //baryon density is not evolving by default
    *rhob = rhobPrev;
    
    *T = effectiveTemperature(*e, *rhob);
    
    *seq = equilibriumEntropy(*e, 0.0, *p, *T, 0.0);
    
    //**************************************************************************
    // SECTION II: baryon evolution is on, but no Hydro+
    //**************************************************************************
    
#else // RootSolver_with_Baryon

#ifndef HydroPlus // HydroPlus

    //**************************** Newton method ****************************
#ifdef Method_Newton
    
    PRECISION v0 = 0;
    PRECISION u0 = 0;

    PRECISION vPrev = sqrt(1-1/(utPrev*utPrev));
    
    v0 = velocityFromConservedVariables(M0, Ms, Pi, N0, vPrev);
    
    if (isnan(v0)) {
        printf("v0 = nan.\t vPrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t rhobPrev=%5e,\t N0=%5e.\n",vPrev,M0,Ms,Pi,rhobPrev,N0);
    }
    
    if(v0<0.563624&&v0>=0)
    {
        *e    = M0 - v0 * Ms;
        *rhob = N0 * sqrt(1 - v0*v0);

        *p = equilibriumPressure(*e, *rhob);
    
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
    else
    {

        u0 = utauFromConservedVariables(M0, Ms, Pi, N0, utPrev);
        
        if (isnan(u0)) {
            printf("u0 = nan.\t utPrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t rhobPrev=%5e,\t N0=%5e.\n",utPrev,M0,Ms,Pi,rhobPrev,N0);
        }
        
        *e    = M0 - Ms * sqrt(1 - 1/(u0*u0));
        *rhob = N0/u0;
        
        *p = equilibriumPressure(*e, *rhob);
        
        PRECISION P  = *p + Pi;
        *ut = u0;
        *ux = u0 * M1/(M0 + P);
        *uy = u0 * M2/(M0 + P);
        *un = u0 * M3/(M0 + P);
    }
    
    if (isnan(*e)||isinf(*e)) {
        printf("Error: e is nan or inf. Quit...\n");
        exit (EXIT_FAILURE);
    }
    
//     //what should we do when solution for flow velocity is too large?
// 	//MUSIC has revert_grid in this case...
// 	if ( *ut > 1.0e2 )
// 	{
// 		printf("\n found ut = %f in getInferredVariables\n", *ut);
// 		printf("M0 = %.9f , e = %.9f, p = %.9f, Pi = %.9f\n", M0, *e, *p, Pi);
// 		printf("***REGULATING FLOW VELOCITY FOR THIS CELL! (ut -> 1.0, ui -> 0.0)*** \n");
// 		*ut = 1.0;
// 		*ux = 0.0;
// 		*uy = 0.0;
// 		*un = 0.0;
// 	}
    
    // code stablization
      if (*rhob<=1.e-4){
          *rhob = 1.e-4;
      }
    
    *T = effectiveTemperature(*e, *rhob);
    *alphaB = chemicalPotentialOverT(*e, *rhob);
    
    // code stablization
    
    if(isnan(*alphaB)){
        printf("alphaB is nan\n");
        *rhob = 1.e-3;
        *alphaB = chemicalPotentialOverT(*e, *rhob);
    }

    *seq = equilibriumEntropy(*e, *rhob, *p, *T, *alphaB);
    
    //**************************** non-Newton method ****************************
#else
    
    PRECISION v0 = 0;
    PRECISION u0 = 0;
    
    PRECISION vPrev = sqrt(1-1/(utPrev*utPrev));

    v0 = velocityIterationFromConservedVariables(M0, Ms, Pi, N0, vPrev);
    
    if (isnan(v0)) {
        printf("v0 = nan.\t vPrev=%5e,\t ePrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t rhobPrev=%5e,\t d_nbt=%5e.\n",vPrev,ePrev,M0,Ms,Pi,rhobPrev,N0);
    }
    
    if(v0<0.563624&&v0>=0){
        
        *e    = M0 - v0 * Ms;
        *rhob = N0 * sqrt(1 - v0*v0);
        
        *p = equilibriumPressure(*e, *rhob);
        
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
        
        u0 = utauIterationFromConservedVariables(M0, Ms, Pi, N0, utPrev);
        
        if (isnan(u0)) {
            printf("u0 = nan.\t utPrev=%5e,\t ePrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t rhobPrev=%5e,\t d_nbt=%5e.\n",utPrev,ePrev,M0,Ms,Pi,rhobPrev,N0);
        }
        
        *e    = M0 - Ms * sqrt(1 - 1/(u0*u0));
        *rhob = N0/u0;
        
        *p = equilibriumPressure(*e, *rhob);
        
        PRECISION P  = *p + Pi;
        *ut = u0;
        *ux = u0 * M1/(M0 + P);
        *uy = u0 * M2/(M0 + P);
        *un = u0 * M3/(M0 + P);
    }
    
    *T = effectiveTemperature(*e, *rhob);
    *alphaB = chemicalPotentialOverT(*e, *rhob);
    *seq = equilibriumEntropy(*e, *rhob, *p, *T, *alphaB);
    
#endif
    
#else // HydroPlus
    
    //**************************************************************************
    // SECTION III: baryon evolution is on, slow modes are evolved
    //**************************************************************************
    // Note: the root sloving routine below, which is not Newton's method, is adapted from a version of MUSIC
    
    // STEP I: setup

    PRECISION vPrev = sqrt(1-1/(utPrev*utPrev));
    
    PRECISION v = 0;
    
    // STEP II: update v or utau, and then calculate e, rhob and p with contributions from slow modes
    
    InferredVariablesVelocityIterationHydroPlus(e, rhob, p, T, alphaB, seq, &v, equiPhiQ, M0, Ms, Pi, N0, vPrev, PhiQ);
    
    if (isnan(v))
    {
        printf("v = nan.\t vPrev=%5e,\t ePrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t rhobPrev=%5e,\t N0=%5e.\n",vPrev,ePrev,M0,Ms,Pi,rhobPrev,N0);
    }
    
    if(v < 0.563624 && v >= 0){
        *ut  = 1/sqrt(1 - v * v);
    }
    else
    {
        InferredVariablesUtauIterationHydroPlus(e, rhob, p, T, alphaB, seq, ut, equiPhiQ, M0, Ms, Pi, N0, utPrev, PhiQ);
        if (isnan(*ut))
        {
            printf("ut = nan.\t utPrev=%5e,\t ePrev=%5e,\t M0=%5e,\t Ms=%5e,\t Pi=%5e,\t rhobPrev=%5e,\t N0=%5e.\n",utPrev,ePrev,M0,Ms,Pi,rhobPrev,N0);
        }
    }
    
    // STEP III: with updated pressure(+), update ux, uy, un
    
    PRECISION P  = *p + Pi;
    PRECISION ufac = (*ut)/(M0 + P);
    
    *ux = ufac * M1;
    *uy = ufac * M2;
    *un = ufac * M3;
    
#endif // HydroPlus
#endif // RootSolver_with_Baryon
}


/**************************************************************************************************************************************************/
/* kernel for setting primary variables
/**************************************************************************************************************************************************/

void setInferredVariablesKernel(const CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, const FLUID_VELOCITY * const __restrict__ uPrev, FLUID_VELOCITY * const __restrict__ u, PRECISION t, void * latticeParams, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ) {
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    
    int ncx,ncy,ncz;
    ncx = lattice->numComputationalLatticePointsX;
    ncy = lattice->numComputationalLatticePointsY;
    ncz = lattice->numComputationalLatticePointsRapidity;
    
#pragma omp parallel for collapse(3)
#ifdef TILE
#pragma unroll_and_jam
#endif
    
    for(int k = 2; k < ncz-2; ++k) {
        for(int j = 2; j < ncy-2; ++j) {
            for(int i = 2; i < ncx-2; ++i) {
                
                int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                
                //**************************************************************************
                // updated Tmunu, Nbmu, shear and bulk pressure, slow modes
                //**************************************************************************
                
                PRECISION q_s[NUMBER_ALL_EVOLVING_VARIABLES];
                
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
                q_s[NUMBER_CONSERVED_VARIABLES_NO_BULK] = q->Pi[s];
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
#ifdef HydroPlus
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
                {
                    q_s[ALL_NUMBER_CONSERVED_VARIABLES+n] = q->phiQ[n][s];
                }
#endif
            
                // old flow velocity
                PRECISION ut_s = uPrev->ut[s];
                
                //**************************************************************************
                // _e, _rhob, _p, _alphaB and _equiPhiQ etc will have updated values.
                //**************************************************************************
                
                PRECISION _e,_rhob,_p,_T,_alphaB,_ut,_ux,_uy,_un,_seq;
                PRECISION _equiPhiQ[NUMBER_SLOW_MODES];
                
                getInferredVariables(t,q_s,e[s],&_e,&_p,ut_s,&_ut,&_ux,&_uy,&_un,rhob[s],&_rhob,&_T,&_alphaB,&_seq,_equiPhiQ);
                
                //**************************************************************************
                // pass the updated value to external arries
                //**************************************************************************
                // NOTES: T, alphaB and seq have no contributions from slow modes even when hydro+ is on, but pressure will have.
#ifdef VMU                
                if (rhob[s]<=1.e-2){

                    q->nbt[s] = 1.e-6;
                    q->nbn[s] = 1.e-6;
                        
                }
#endif                
                e[s] = _e;
                rhob[s]  = _rhob;
                T[s] = _T;
                seq[s] = _seq;
                
                p[s] = _p;

                u->ut[s] = _ut;
                u->ux[s] = _ux;
                u->uy[s] = _uy;
                u->un[s] = _un;

#ifdef NBMU
                alphaB[s] = _alphaB;
#endif
#ifdef HydroPlus
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
                {
                    eqPhiQ->phiQ[n][s] = _equiPhiQ[n];
                }
#endif
            }
        }
    }
}


/**************************************************************************************************************************************************/
/* Components of T^{\mu\nu} and Nb^{\mu} in (\tau,x,y,\eta_s)-coordinates
/**************************************************************************************************************************************************/

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

