/*
 * DynamicalVariables.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <stdlib.h>
#include <math.h>

#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/PrimaryVariables.h"
#include "../include/FullyDiscreteKurganovTadmorScheme.h"
#include "../include/EquationOfState.h"
#include "../include/HydroPlus.h"

/**************************************************************************************************************************************************/
/* instances
/**************************************************************************************************************************************************/

CONSERVED_VARIABLES *q, *Q, *qS;

FLUID_VELOCITY *u, *up, *uS;

PRECISION *e, *p, *seq;

PRECISION *rhob, *rhobp, *rhobS;

PRECISION *alphaB, *alphaBp, *alphaBS;

PRECISION *T, *Tp, *TS;

EQUATION_OF_STATE *EOState;

BARYON_DIFFUSION_COEFF *BaryDiffCoeff;

DYNAMICAL_SOURCE *Source;

PRECISION *Qvec; // Q vectors of slow modes

SLOW_MODES *eqPhiQ; // Slow modes at equilibrium

PRECISION *xieq; // correlation length


/**************************************************************************************************************************************************/
/* column compact index
/**************************************************************************************************************************************************/

int columnMajorLinearIndex(int i, int j, int k, int nx, int ny) {
	return i + nx * (j + ny * k);
}


/**************************************************************************************************************************************************/
/* allocate memory for all external variables
/**************************************************************************************************************************************************/

void allocateHostMemory(int len) {
	size_t bytes = sizeof(PRECISION);

	//=======================================================================
	// Primary variables
	//=======================================================================
    
    // energy density, pressure, baryon density, entropy density
	e = (PRECISION *)calloc(len, bytes);
	p = (PRECISION *)calloc(len, bytes);
    seq = (PRECISION *)calloc(len, bytes);
    
    rhob = (PRECISION *)calloc(len, bytes);
    rhobp = (PRECISION *)calloc(len, bytes);
    rhobS = (PRECISION *)calloc(len, bytes);
    
    // baryon chemical potential
    alphaB = (PRECISION *)calloc(len, bytes);
    alphaBp= (PRECISION *)calloc(len, bytes);
    alphaBS= (PRECISION *)calloc(len, bytes);
    
    // temperature
    T = (PRECISION *)calloc(len, bytes);
    Tp= (PRECISION *)calloc(len, bytes);
    TS= (PRECISION *)calloc(len, bytes);
    
	// fluid velocity at current time step
	u = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	u->ut = (PRECISION *)calloc(len,bytes);
	u->ux = (PRECISION *)calloc(len,bytes);
	u->uy = (PRECISION *)calloc(len,bytes);
	u->un = (PRECISION *)calloc(len,bytes);
    
	// fluid velocity at previous time step
	up = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	up->ut = (PRECISION *)calloc(len,bytes);
	up->ux = (PRECISION *)calloc(len,bytes);
	up->uy = (PRECISION *)calloc(len,bytes);
	up->un = (PRECISION *)calloc(len,bytes);
    
	// fluid velocity at intermediate time step
	uS = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	uS->ut = (PRECISION *)calloc(len,bytes);
	uS->ux = (PRECISION *)calloc(len,bytes);
	uS->uy = (PRECISION *)calloc(len,bytes);
	uS->un = (PRECISION *)calloc(len,bytes);

	//=======================================================================
	// Conserved variables
	//=======================================================================
    
    // current variables at the n time step
	q = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	q->ttt = (PRECISION *)calloc(len, bytes);
	q->ttx = (PRECISION *)calloc(len, bytes);
	q->tty = (PRECISION *)calloc(len, bytes);
	q->ttn = (PRECISION *)calloc(len, bytes);
#ifdef PIMUNU
	q->pitt = (PRECISION *)calloc(len, bytes);
	q->pitx = (PRECISION *)calloc(len, bytes);
	q->pity = (PRECISION *)calloc(len, bytes);
	q->pitn = (PRECISION *)calloc(len, bytes);
	q->pixx = (PRECISION *)calloc(len, bytes);
	q->pixy = (PRECISION *)calloc(len, bytes);
	q->pixn = (PRECISION *)calloc(len, bytes);
	q->piyy = (PRECISION *)calloc(len, bytes);
	q->piyn = (PRECISION *)calloc(len, bytes);
	q->pinn = (PRECISION *)calloc(len, bytes);
#endif
#ifdef PI
	q->Pi = (PRECISION *)calloc(len, bytes);
#endif
#ifdef NBMU
    q->Nbt = (PRECISION *)calloc(len, bytes);
#endif
#ifdef VMU
    q->nbt = (PRECISION *)calloc(len, bytes);
    q->nbx = (PRECISION *)calloc(len, bytes);
    q->nby = (PRECISION *)calloc(len, bytes);
    q->nbn = (PRECISION *)calloc(len, bytes);
#endif
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        q->phiQ[n]  = (PRECISION *)calloc(len, bytes);
    }

	// upated variables at the n+1 time step
	Q = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	Q->ttt = (PRECISION *)calloc(len, bytes);
	Q->ttx = (PRECISION *)calloc(len, bytes);
	Q->tty = (PRECISION *)calloc(len, bytes);
	Q->ttn = (PRECISION *)calloc(len, bytes);
#ifdef PIMUNU
	Q->pitt = (PRECISION *)calloc(len, bytes);
	Q->pitx = (PRECISION *)calloc(len, bytes);
	Q->pity = (PRECISION *)calloc(len, bytes);
	Q->pitn = (PRECISION *)calloc(len, bytes);
	Q->pixx = (PRECISION *)calloc(len, bytes);
	Q->pixy = (PRECISION *)calloc(len, bytes);
	Q->pixn = (PRECISION *)calloc(len, bytes);
	Q->piyy = (PRECISION *)calloc(len, bytes);
	Q->piyn = (PRECISION *)calloc(len, bytes);
	Q->pinn = (PRECISION *)calloc(len, bytes);
#endif
#ifdef PI
	Q->Pi = (PRECISION *)calloc(len, bytes);
#endif
#ifdef NBMU
    Q->Nbt = (PRECISION *)calloc(len, bytes);
#endif
#ifdef VMU
    Q->nbt = (PRECISION *)calloc(len, bytes);
    Q->nbx = (PRECISION *)calloc(len, bytes);
    Q->nby = (PRECISION *)calloc(len, bytes);
    Q->nbn = (PRECISION *)calloc(len, bytes);
#endif
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        Q->phiQ[n]  = (PRECISION *)calloc(len, bytes);
    }

	// updated variables at the intermediate time step
	qS = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	qS->ttt = (PRECISION *)calloc(len, bytes);
	qS->ttx = (PRECISION *)calloc(len, bytes);
	qS->tty = (PRECISION *)calloc(len, bytes);
	qS->ttn = (PRECISION *)calloc(len, bytes);
#ifdef PIMUNU
	qS->pitt = (PRECISION *)calloc(len, bytes);
	qS->pitx = (PRECISION *)calloc(len, bytes);
	qS->pity = (PRECISION *)calloc(len, bytes);
	qS->pitn = (PRECISION *)calloc(len, bytes);
	qS->pixx = (PRECISION *)calloc(len, bytes);
	qS->pixy = (PRECISION *)calloc(len, bytes);
	qS->pixn = (PRECISION *)calloc(len, bytes);
	qS->piyy = (PRECISION *)calloc(len, bytes);
	qS->piyn = (PRECISION *)calloc(len, bytes);
	qS->pinn = (PRECISION *)calloc(len, bytes);
#endif
#ifdef PI
	qS->Pi = (PRECISION *)calloc(len, bytes);
#endif
#ifdef NBMU
    qS->Nbt = (PRECISION *)calloc(len, bytes);
#endif
#ifdef VMU
    qS->nbt = (PRECISION *)calloc(len, bytes);
    qS->nbx = (PRECISION *)calloc(len, bytes);
    qS->nby = (PRECISION *)calloc(len, bytes);
    qS->nbn = (PRECISION *)calloc(len, bytes);
#endif
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        qS->phiQ[n]  = (PRECISION *)calloc(len, bytes);
    }
    
    //=======================================================================
    // dynamical sources and Equation of State
    //=======================================================================
    
    // dynamical source terms
    Source = (DYNAMICAL_SOURCE *)calloc(1, sizeof(DYNAMICAL_SOURCE));
    Source->sourcet = (PRECISION *)calloc(len, bytes);
    Source->sourcex = (PRECISION *)calloc(len, bytes);
    Source->sourcey = (PRECISION *)calloc(len, bytes);
    Source->sourcen = (PRECISION *)calloc(len, bytes);
    Source->sourceb = (PRECISION *)calloc(len, bytes);
    
    // equation of state table
    EOState = (EQUATION_OF_STATE *)calloc(1, sizeof(EQUATION_OF_STATE));
    EOState->Pressure      = (PRECISION *)calloc(188580, bytes);
    EOState->Temperature   = (PRECISION *)calloc(188580, bytes);
    EOState->alphab      = (PRECISION *)calloc(188580, bytes);
    EOState->dpdrhob       = (PRECISION *)calloc(188580, bytes);
    
    // baryon diffusion coefficients
#ifdef VMU
    BaryDiffCoeff = (BARYON_DIFFUSION_COEFF *)calloc(1, sizeof(BARYON_DIFFUSION_COEFF));
    //BaryDiffCoeff->sigmaB = (PRECISION *)calloc(5751, bytes);
    BaryDiffCoeff->sigmaB = (PRECISION *)calloc(128721, bytes);
    BaryDiffCoeff->DB = (PRECISION *)calloc(128721, bytes);
#endif
    
    //=======================================================================
    // Slow modes at equilibrium
    //=======================================================================
    
#ifdef HydroPlus
    Qvec = (PRECISION *)calloc(NUMBER_SLOW_MODES, bytes);
    
    eqPhiQ  = (SLOW_MODES *)calloc(1, sizeof(SLOW_MODES));
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        eqPhiQ->phiQ[n]  = (PRECISION *)calloc(len, bytes);
    }
#endif
#ifdef CRITICAL
    xieq = (PRECISION *)calloc(9396, bytes);
#endif
}


/**************************************************************************************************************************************************/
/* set values for variables, e.g. T^tau^mu and N^tau calculated from flow velocity and energy density etc. at the beginning of the event
/**************************************************************************************************************************************************/

void setConservedVariables(double t, void * latticeParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	int ncx = lattice->numComputationalLatticePointsX;
	int ncy = lattice->numComputationalLatticePointsY;

	for (int k = N_GHOST_CELLS_M; k < nz+N_GHOST_CELLS_M; ++k) {
		for (int j = N_GHOST_CELLS_M; j < ny+N_GHOST_CELLS_M; ++j) {
			for (int i = N_GHOST_CELLS_M; i < nx+N_GHOST_CELLS_M; ++i) {
                
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);

                //=======================================================================
                // u, e, p, pi and Pi, rhob, nb initialized in initial condition
                //=======================================================================
                
				PRECISION ux_s = u->ux[s];
				PRECISION uy_s = u->uy[s];
				PRECISION un_s = u->un[s];
				PRECISION ut_s = u->ut[s];
                
				PRECISION e_s = e[s];
				PRECISION p_s = p[s];

				PRECISION pitt_s = 0;
				PRECISION pitx_s = 0;
				PRECISION pity_s = 0;
				PRECISION pitn_s = 0;
#ifdef PIMUNU
				pitt_s = q->pitt[s];
				pitx_s = q->pitx[s];
				pity_s = q->pity[s];
				pitn_s = q->pitn[s];
#endif
				PRECISION Pi_s = 0;
#ifdef PI
				Pi_s = q->Pi[s];
#endif

                //=======================================================================
                // initialize Tmunu, Nbmu, Temperature and potential
                //=======================================================================
                
				q->ttt[s] = Ttt(e_s, p_s+Pi_s, ut_s, pitt_s);
				q->ttx[s] = Ttx(e_s, p_s+Pi_s, ut_s, ux_s, pitx_s);
				q->tty[s] = Tty(e_s, p_s+Pi_s, ut_s, uy_s, pity_s);
				q->ttn[s] = Ttn(e_s, p_s+Pi_s, ut_s, un_s, pitn_s);

                PRECISION rhob_s = rhob[s];
                
                PRECISION nbt_s = 0;
                T[s] = effectiveTemperature(e_s, rhob_s);
                
#ifdef VMU
                nbt_s = q->nbt[s];
#endif
#ifdef NBMU
                q->Nbt[s] = Nbt(rhob_s, ut_s, nbt_s);

                alphaB[s] = chemicalPotentialOverT(e_s, rhob_s);

                seq[s] = equilibriumEntropy(e[s], rhob[s], p[s], T[s], alphaB[s]);
#endif
                
                //=======================================================================
                // initialize quantities at the previous step
                //=======================================================================
                
                up->ux[s] = ux_s;
                up->uy[s] = uy_s;
                up->un[s] = un_s;
                up->ut[s] = ut_s;
                
                rhobp[s] = rhob[s];
                
                Tp[s] = T[s];
                alphaBp[s] = alphaB[s];
			}
		}
	}
}


/**************************************************************************************************************************************************/
/* set values for ghost cells
/**************************************************************************************************************************************************/

void setGhostCells(CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, void * latticeParams, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ)
{
	setGhostCellsKernelI(q,e,p,u,latticeParams,rhob,alphaB,T,seq,eqPhiQ);
	setGhostCellsKernelJ(q,e,p,u,latticeParams,rhob,alphaB,T,seq,eqPhiQ);
	setGhostCellsKernelK(q,e,p,u,latticeParams,rhob,alphaB,T,seq,eqPhiQ);
}

void setGhostCellVars(CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, int s, int sBC, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ)
{
    e[s] = e[sBC];
    p[s] = p[sBC];
    u->ut[s] = u->ut[sBC];
    u->ux[s] = u->ux[sBC];
    u->uy[s] = u->uy[sBC];
    u->un[s] = u->un[sBC];
    q->ttt[s] = q->ttt[sBC];
    q->ttx[s] = q->ttx[sBC];
    q->tty[s] = q->tty[sBC];
    q->ttn[s] = q->ttn[sBC];
#ifdef PIMUNU
    q->pitt[s] = q->pitt[sBC];
    q->pitx[s] = q->pitx[sBC];
    q->pity[s] = q->pity[sBC];
    q->pitn[s] = q->pitn[sBC];
    q->pixx[s] = q->pixx[sBC];
    q->pixy[s] = q->pixy[sBC];
    q->pixn[s] = q->pixn[sBC];
    q->piyy[s] = q->piyy[sBC];
    q->piyn[s] = q->piyn[sBC];
    q->pinn[s] = q->pinn[sBC];
#endif
#ifdef PI
    q->Pi[s] = q->Pi[sBC];
#endif
    
    seq[s] = seq[sBC];
    rhob[s] = rhob[sBC];
    T[s] = T[sBC];
#ifdef NBMU
    q->Nbt[s] = q->Nbt[sBC];
    alphaB[s] = alphaB[sBC];
#endif
#ifdef VMU
    q->nbt[s] = q->nbt[sBC];
    q->nbx[s] = q->nbx[sBC];
    q->nby[s] = q->nby[sBC];
    q->nbn[s] = q->nbn[sBC];
#endif
    
#ifdef HydroPlus
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        eqPhiQ->phiQ[n][s] = eqPhiQ->phiQ[n][sBC];
        q->phiQ[n][s] = q->phiQ[n][sBC];
    }
#endif
}

void setGhostCellsKernelI(CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, void * latticeParams, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ)
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int nx,ncx,ncy,ncz;
	nx = lattice->numLatticePointsX;
	ncx = lattice->numComputationalLatticePointsX;
	ncy = lattice->numComputationalLatticePointsY;
	ncz = lattice->numComputationalLatticePointsRapidity;

	int iBC,s,sBC;
	for(int j = 2; j < ncy; ++j) {
		for(int k = 2; k < ncz; ++k) {
			iBC = 2;
			for (int i = 0; i <= 1; ++i) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				sBC = columnMajorLinearIndex(iBC, j, k, ncx, ncy);
				setGhostCellVars(q,e,p,u,s,sBC,rhob,alphaB,T,seq,eqPhiQ);
			}
			iBC = nx + 1;
			for (int i = nx + 2; i <= nx + 3; ++i) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				sBC = columnMajorLinearIndex(iBC, j, k, ncx, ncy);
				setGhostCellVars(q,e,p,u,s,sBC,rhob,alphaB,T,seq,eqPhiQ);
			}
		}
	}
}

void setGhostCellsKernelJ(CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, void * latticeParams, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ)
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int ny,ncx,ncy,ncz;
	ny = lattice->numLatticePointsY;
	ncx = lattice->numComputationalLatticePointsX;
	ncy = lattice->numComputationalLatticePointsY;
	ncz = lattice->numComputationalLatticePointsRapidity;

	int jBC,s,sBC;
	for(int i = 2; i < ncx; ++i) {
		for(int k = 2; k < ncz; ++k) {
			jBC = 2;
			for (int j = 0; j <= 1; ++j) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				sBC = columnMajorLinearIndex(i, jBC, k, ncx, ncy);
				setGhostCellVars(q,e,p,u,s,sBC,rhob,alphaB,T,seq,eqPhiQ);
			}
			jBC = ny + 1;
			for (int j = ny + 2; j <= ny + 3; ++j) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				sBC = columnMajorLinearIndex(i, jBC, k, ncx, ncy);
				setGhostCellVars(q,e,p,u,s,sBC,rhob,alphaB,T,seq,eqPhiQ);
			}
		}
	}
}

void setGhostCellsKernelK(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq,
SLOW_MODES *  const __restrict__ eqPhiQ)
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int nz,ncx,ncy;
	nz = lattice->numLatticePointsRapidity;
	ncx = lattice->numComputationalLatticePointsX;
	ncy = lattice->numComputationalLatticePointsY;

	int kBC,s,sBC;
	for(int i = 2; i < ncx; ++i) {
		for(int j = 2; j < ncy; ++j) {
			kBC = 2;
			for (int k = 0; k <= 1; ++k) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				sBC = columnMajorLinearIndex(i, j, kBC, ncx, ncy);
				setGhostCellVars(q,e,p,u,s,sBC,rhob,alphaB,T,seq,eqPhiQ);
			}
			kBC = nz + 1;
			for (int k = nz + 2; k <= nz + 3; ++k) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				sBC = columnMajorLinearIndex(i, j, kBC, ncx, ncy);
				setGhostCellVars(q,e,p,u,s,sBC,rhob,alphaB,T,seq,eqPhiQ);
			}
		}
	}
}


/**************************************************************************************************************************************************/
/* give the current value of variables to be the ones of previous step
/**************************************************************************************************************************************************/

void swap(CONSERVED_VARIABLES **arr1, CONSERVED_VARIABLES **arr2) {
    CONSERVED_VARIABLES *tmp = *arr1;
    *arr1 = *arr2;
    *arr2 = tmp;
}

void setCurrentConservedVariables() {
	swap(&q, &Q);
}

void swapFluidVelocity(FLUID_VELOCITY **arr1, FLUID_VELOCITY **arr2) {
	FLUID_VELOCITY *tmp = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;
}

void swapPrimaryVariables(PRECISION **arr1, PRECISION **arr2) {
    PRECISION *tmp = *arr1;
    *arr1 = *arr2;
    *arr2 = tmp;
}


/**************************************************************************************************************************************************/
/* free host memory
/**************************************************************************************************************************************************/

void freeHostMemory() {
    free(Source->sourcet);
    free(Source->sourcex);
    free(Source->sourcey);
    free(Source->sourcen);
    free(Source->sourceb);

    free(EOState->Pressure);
    free(EOState->Temperature);
    free(EOState->alphab);
    free(EOState->dpdrhob);
#ifdef VMU
    free(BaryDiffCoeff->sigmaB);
    free(BaryDiffCoeff->DB);
#endif

    free(rhob);
    free(rhobp);
    free(rhobS);
    free(alphaB);
    free(alphaBp);
    free(alphaBS);
    free(T);
    free(Tp);
    free(TS);
#ifdef NBMU
    free(q->Nbt);
#endif
#ifdef VMU
    free(q->nbt);
    free(q->nbx);
    free(q->nby);
    free(q->nbn);
#endif

	free(e);
	free(p);
    free(seq);
	free(u->ut);
	free(u->ux);
	free(u->uy);
	free(u->un);
    
    free(Qvec);
    
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
        free(eqPhiQ->phiQ[n]);
    }
    
    free(eqPhiQ);
#ifdef CRITICAL
    free(xieq);
#endif

	free(q->ttt);
	free(q->ttx);
	free(q->tty);
	free(q->ttn);
#ifdef PIMUNU
	free(q->pitt);
	free(q->pitx);
	free(q->pity);
	free(q->pitn);
	free(q->pixx);
	free(q->pixy);
	free(q->pixn);
	free(q->piyy);
	free(q->piyn);
	free(q->pinn);
#endif
#ifdef PI
	free(q->Pi);
#endif
	free(q);
}
