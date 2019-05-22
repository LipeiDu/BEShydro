/*
 * DynamicalVariables.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <stdlib.h>
#include <math.h>

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/trunk/hydro/EnergyMomentumTensor.h"
#include "edu/osu/rhic/trunk/hydro/FullyDiscreteKurganovTadmorScheme.h" // for ghost cells
#include "edu/osu/rhic/trunk/eos/EquationOfState.h"//Lipei

CONSERVED_VARIABLES *q, *Q, *qS;

FLUID_VELOCITY *u, *up, *uS;

PRECISION *e, *p, *rhob;

PRECISION *muB, *muBp, *muBS;

PRECISION *T, *Tp, *TS;

EQUATION_OF_STATE *EOState;//Lipei

DYNAMICAL_SOURCE *Source;//Lipei

PRECISION *termX;
PRECISION *termY;
PRECISION *termZ;
PRECISION *term2;//Lipei

int columnMajorLinearIndex(int i, int j, int k, int nx, int ny) {
	return i + nx * (j + ny * k);
}

void allocateHostMemory(int len) {
	size_t bytes = sizeof(PRECISION);

    // dynamical source terms
    Source = (DYNAMICAL_SOURCE *)calloc(1, sizeof(DYNAMICAL_SOURCE));
    Source->sourcet = (PRECISION *)calloc(len, bytes);//Lipei
    Source->sourcex = (PRECISION *)calloc(len, bytes);//Lipei
    Source->sourcey = (PRECISION *)calloc(len, bytes);//Lipei
    Source->sourcen = (PRECISION *)calloc(len, bytes);//Lipei
    Source->sourceb = (PRECISION *)calloc(len, bytes);//Lipei
    // baryon density
    rhob = (PRECISION *)calloc(len, bytes);//Lipei
    // baryon chemical potential
    muB = (PRECISION *)calloc(len, bytes);//Lipei
    muBp= (PRECISION *)calloc(len, bytes);//Lipei
    muBS= (PRECISION *)calloc(len, bytes);//Lipei
    // temperature
    T = (PRECISION *)calloc(len, bytes);//Lipei
    Tp= (PRECISION *)calloc(len, bytes);//Lipei
    TS= (PRECISION *)calloc(len, bytes);//Lipei
    
    termX = (PRECISION *)calloc(len, bytes);//Lipei
    termY = (PRECISION *)calloc(len, bytes);//Lipei
    termZ = (PRECISION *)calloc(len, bytes);//Lipei
    term2 = (PRECISION *)calloc(len, bytes);//Lipei

    // equation of state table
    EOState = (EQUATION_OF_STATE *)calloc(1, sizeof(EQUATION_OF_STATE));
    EOState->Pressure      = (PRECISION *)calloc(188580, bytes);//Lipei
    EOState->Temperature   = (PRECISION *)calloc(188580, bytes);//Lipei
    EOState->Mubovert      = (PRECISION *)calloc(188580, bytes);//Lipei
    EOState->dpdrhob       = (PRECISION *)calloc(188580, bytes);//Lipei
    
    EOState->sigmaB = (PRECISION *)calloc(5751, bytes);//Lipei
	//=======================================================
	// Primary variables
	//=======================================================
	e = (PRECISION *)calloc(len, bytes);
	p = (PRECISION *)calloc(len, bytes);
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

	//=======================================================
	// Conserved variables
	//=======================================================
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
}

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

				q->ttt[s] = Ttt(e_s, p_s+Pi_s, ut_s, pitt_s);
				q->ttx[s] = Ttx(e_s, p_s+Pi_s, ut_s, ux_s, pitx_s);
				q->tty[s] = Tty(e_s, p_s+Pi_s, ut_s, uy_s, pity_s);
				q->ttn[s] = Ttn(e_s, p_s+Pi_s, ut_s, un_s, pitn_s);

                PRECISION rhob_s = rhob[s];
                PRECISION nbt_s = 0;
                T[s] = effectiveTemperature(e_s, rhob_s);
                if (T[s] < 1.e-7)
                    T[s] = 1.e-7;
                
                Tp[s] = T[s];
#ifdef VMU
                nbt_s = q->nbt[s];
#endif
#ifdef NBMU
                q->Nbt[s] = Nbt(rhob_s, ut_s, nbt_s);

                muB[s] = chemicalPotentialOverT(e_s, rhob_s);
                
                if (muB[s]>=0 && muB[s] < 1.e-7) muB[s] = 1.e-7;
                else if (muB[s]<=0 && muB[s] > -1.e-7)  muB[s] = -1.e-7;

                muBp[s] = muB[s];
                term2[s] = muB[s]*T[s];
#endif
			}
		}
	}
}

void setGhostCells(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ muB, PRECISION * const __restrict__ T//by Lipei
) {
	setGhostCellsKernelI(q,e,p,u,latticeParams,rhob,muB,T);
	setGhostCellsKernelJ(q,e,p,u,latticeParams,rhob,muB,T);
	setGhostCellsKernelK(q,e,p,u,latticeParams,rhob,muB,T);
}

void setGhostCellVars(CONSERVED_VARIABLES * const __restrict__ q,
                      PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
                      FLUID_VELOCITY * const __restrict__ u,
                      int s, int sBC,
                      PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ muB, PRECISION * const __restrict__ T//by Lipei
) {
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
    // set \pi^\mu\nu ghost cells if evolved
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
    // set \Pi ghost cells if evolved
#ifdef PI
    q->Pi[s] = q->Pi[sBC];
#endif
    
    rhob[s] = rhob[sBC];
    T[s] = T[sBC];
#ifdef NBMU
    q->Nbt[s] = q->Nbt[sBC];
    muB[s] = muB[sBC];
#endif
#ifdef VMU
    q->nbt[s] = q->nbt[sBC];
    q->nbx[s] = q->nbx[sBC];
    q->nby[s] = q->nby[sBC];
    q->nbn[s] = q->nbn[sBC];
#endif
}

void setGhostCellsKernelI(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ muB, PRECISION * const __restrict__ T//by Lipei
) {
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
				setGhostCellVars(q,e,p,u,s,sBC,rhob,muB,T);
			}
			iBC = nx + 1;
			for (int i = nx + 2; i <= nx + 3; ++i) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				sBC = columnMajorLinearIndex(iBC, j, k, ncx, ncy);
				setGhostCellVars(q,e,p,u,s,sBC,rhob,muB,T);
			}
		}
	}
}

void setGhostCellsKernelJ(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ muB, PRECISION * const __restrict__ T//by Lipei
) {
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
				setGhostCellVars(q,e,p,u,s,sBC,rhob,muB,T);
			}
			jBC = ny + 1;
			for (int j = ny + 2; j <= ny + 3; ++j) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				sBC = columnMajorLinearIndex(i, jBC, k, ncx, ncy);
				setGhostCellVars(q,e,p,u,s,sBC,rhob,muB,T);
			}
		}
	}
}

void setGhostCellsKernelK(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ muB, PRECISION * const __restrict__ T//by Lipei
) {
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
				setGhostCellVars(q,e,p,u,s,sBC,rhob,muB,T);
			}
			kBC = nz + 1;
			for (int k = nz + 2; k <= nz + 3; ++k) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				sBC = columnMajorLinearIndex(i, j, kBC, ncx, ncy);
				setGhostCellVars(q,e,p,u,s,sBC,rhob,muB,T);
			}
		}
	}
}

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

//swap energy density or baryon density; Lipei
void swapPrimaryVariables(PRECISION **arr1, PRECISION **arr2) {
    PRECISION *tmp = *arr1;
    *arr1 = *arr2;
    *arr2 = tmp;
}

void freeHostMemory() {
    free(Source->sourcet);//Lipei
    free(Source->sourcex);//Lipei
    free(Source->sourcey);//Lipei
    free(Source->sourcen);//Lipei
    free(Source->sourceb);//Lipei
    
    free(termX);//Lipei
    free(termY);//Lipei
    free(termZ);//Lipei
    free(term2);//Lipei
    

    //free(EOState->ChemicalPotential);//Lipei
    free(EOState->Pressure);//Lipei
    free(EOState->Temperature);//Lipei
    free(EOState->Mubovert);//Lipei
    free(EOState->dpdrhob);
    
    free(EOState->sigmaB);//Lipei

    free(rhob);//Lipei
    free(muB);
    free(muBp);
    free(muBS);
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
	free(u->ut);
	free(u->ux);
	free(u->uy);
	free(u->un);

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
