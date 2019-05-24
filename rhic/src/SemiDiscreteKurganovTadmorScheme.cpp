//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <stdio.h> // for printf
#include <math.h> 

#include "../include/SemiDiscreteKurganovTadmorScheme.h"
#include "../include/DynamicalVariables.h"
#include "../include/PrimaryVariables.h"
#include "../include/LocalPropagationSpeed.h"
 
void flux(const PRECISION * const __restrict__ data, PRECISION * const __restrict__ result,
		PRECISION (* const rightHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
		PRECISION (* const leftHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
		PRECISION (* const spectralRadius)(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un),
		PRECISION (* const fluxFunction)(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un),
		PRECISION t, PRECISION ePrev, PRECISION rhobPrev, PRECISION utPrev
) {
	// left and right cells
	PRECISION qR[NUMBER_ALL_EVOLVING_VARIABLES], qL[NUMBER_ALL_EVOLVING_VARIABLES];

	// left and right extrapolated values of the conserved variables
	int ptr = 0;
	PRECISION qmm, qm, q, qp, qpp;
	for (unsigned int n = 0; n < NUMBER_ALL_EVOLVING_VARIABLES; ++n) {
		qmm = *(data+ptr);
		qm 	= *(data+ptr+1);
		q 	= *(data+ptr+2);
		qp 	= *(data+ptr+3);
		qpp = *(data+ptr+4);
		ptr+=5;
		qR[n]	= rightHalfCellExtrapolation(qmm, qm, q, qp, qpp);
		qL[n]	= leftHalfCellExtrapolation(qmm, qm, q, qp, qpp);
	}

	// left and right extrapolated values of the primary variables
	PRECISION eR,pR,utR,uxR,uyR,unR,rhobR,TR,alphaBR,sR,equiPhiQR[NUMBER_SLOW_MODES];
	getInferredVariables(t,qR,ePrev,&eR,&pR,utPrev,&utR,&uxR,&uyR,&unR,rhobPrev,&rhobR,&TR,&alphaBR,&sR,equiPhiQR);
	PRECISION eL,pL,utL,uxL,uyL,unL,rhobL,TL,alphaBL,sL,equiPhiQL[NUMBER_SLOW_MODES];
	getInferredVariables(t,qL,ePrev,&eL,&pL,utPrev,&utL,&uxL,&uyL,&unL,rhobPrev,&rhobL,&TL,&alphaBL,&sL,equiPhiQL);

	PRECISION a,qR_n,qL_n,FqR,FqL,res;
	a = localPropagationSpeed(utR,uxR,uyR,unR,utL,uxL,uyL,unL,spectralRadius);
	for (unsigned int n = 0; n < NUMBER_ALL_EVOLVING_VARIABLES; ++n) {
		qR_n = qR[n];
		qL_n = qL[n];
		FqR = fluxFunction(qR_n, utR, uxR, uyR, unR);
		FqL = fluxFunction(qL_n, utL, uxL, uyL, unL);
		res = FqR + FqL - a * (qR_n - qL_n);
		res /= 2;
		result[n] = res; 
	}
}
