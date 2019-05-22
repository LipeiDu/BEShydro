/*
 * DynamicalVariables.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_

/*********************************************************/
//Main switch//

#define PIMUNU
//#define PI

#define NBMU
#define VMU

#ifdef NBMU
#define RootSolver_with_Baryon//Lipei
#define EOS_with_baryon
#endif

/*********************************************************/
//Conservation laws//

#define NUMBER_CONSERVATION_LAWS 4

#ifndef NBMU
#define NUMBER_BARYON_COMPONENTS 0
#else
#define NUMBER_BARYON_COMPONENTS 1
#endif

#define ALL_NUMBER_CONSERVATION_LAWS (NUMBER_CONSERVATION_LAWS+NUMBER_BARYON_COMPONENTS)

/*********************************************************/
//Pressures//

#ifndef PI
#define NUMBER_PI_COMPONENTS 0
#else
#define NUMBER_PI_COMPONENTS 1
#endif

#ifndef PIMUNU
#define NUMBER_PROPAGATED_PIMUNU_COMPONENTS 0
#else
#define NUMBER_PROPAGATED_PIMUNU_COMPONENTS 10
#endif

/*********************************************************/
//Dissipative currents//

#ifndef VMU
#define NUMBER_PROPAGATED_VMU_COMPONENTS 0
#else
#define NUMBER_PROPAGATED_VMU_COMPONENTS 4
#endif

#define NUMBER_DISSIPATIVE_CURRENTS (NUMBER_PI_COMPONENTS+NUMBER_PROPAGATED_PIMUNU_COMPONENTS)

#define ALL_NUMBER_DISSIPATIVE_CURRENTS (NUMBER_PI_COMPONENTS+NUMBER_PROPAGATED_VMU_COMPONENTS+NUMBER_PROPAGATED_PIMUNU_COMPONENTS)

#if ALL_NUMBER_DISSIPATIVE_CURRENTS==0
#define IDEAL
#endif

#define NUMBER_CONSERVED_VARIABLES (NUMBER_CONSERVATION_LAWS+NUMBER_DISSIPATIVE_CURRENTS)//excluding baryon part

#define ALL_NUMBER_CONSERVED_VARIABLES (ALL_NUMBER_CONSERVATION_LAWS+ALL_NUMBER_DISSIPATIVE_CURRENTS)//baryon part is defined seperately

/*********************************************************/

#define PRECISION double

typedef struct
{
	PRECISION *ttt;
	PRECISION *ttx;
	PRECISION *tty;
	PRECISION *ttn;
#ifdef PIMUNU
	PRECISION *pitt;
	PRECISION *pitx;
	PRECISION *pity;
	PRECISION *pitn;
	PRECISION *pixx;
	PRECISION *pixy;
	PRECISION *pixn;
	PRECISION *piyy;
	PRECISION *piyn;
	PRECISION *pinn;
#endif
#ifdef PI
	PRECISION *Pi;
#endif
#ifdef NBMU
    PRECISION *Nbt;
#endif
#ifdef VMU
    PRECISION *nbt;
    PRECISION *nbx;
    PRECISION *nby;
    PRECISION *nbn;
#endif
} CONSERVED_VARIABLES;

typedef struct
{
	PRECISION *ut;
	PRECISION *ux;
	PRECISION *uy;
	PRECISION *un;
} FLUID_VELOCITY;

typedef struct
{
    PRECISION *sourcet;
    PRECISION *sourcex;
    PRECISION *sourcey;
    PRECISION *sourcen;
    PRECISION *sourceb;
} DYNAMICAL_SOURCE;

typedef struct
{
    PRECISION *Pressure;
    PRECISION *Temperature;
    PRECISION *Mubovert;
    PRECISION *dpdrhob;
    PRECISION *sigmaB;
} EQUATION_OF_STATE;

extern PRECISION *termX;
extern PRECISION *termY;
extern PRECISION *termZ;
extern PRECISION *term2;
//test Lipei


extern CONSERVED_VARIABLES *q,*Q,*qS;
extern FLUID_VELOCITY *u,*up,*uS;
extern DYNAMICAL_SOURCE *Source;
extern EQUATION_OF_STATE *EOState;

extern PRECISION *e, *p, *rhob;
extern PRECISION *muB, *muBp, *muBS;
extern PRECISION *T, *Tp, *TS;

int columnMajorLinearIndex(int i, int j, int k, int nx, int ny);

void allocateHostMemory(int len);

void setConservedVariables(double t, void * latticeParams);
void setCurrentConservedVariables();
void swapFluidVelocity(FLUID_VELOCITY **arr1, FLUID_VELOCITY **arr2) ;
void swapPrimaryVariables(PRECISION **arr1, PRECISION **arr2);

void setGhostCells(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ muB, PRECISION * const __restrict__ T
);

void setGhostCellsKernelI(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ muB, PRECISION * const __restrict__ T
);

void setGhostCellsKernelJ(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ muB, PRECISION * const __restrict__ T
);

void setGhostCellsKernelK(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams,
PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ muB, PRECISION * const __restrict__ T
);

void freeHostMemory();

#endif /* DYNAMICALVARIABLES_H_ */
