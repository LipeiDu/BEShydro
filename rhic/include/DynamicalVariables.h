//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#ifndef DYNAMICALVARIABLES_H_
#define DYNAMICALVARIABLES_H_

/**************************************************************************************************************************************************/
/* components: energy-momentun tensor, shear stress, bulk pressure, net baryon current, baryon diffusion current, slow modes
/**************************************************************************************************************************************************/

/*********************************************************/
//Main switch//

#define PIMUNU
//#define PI

//#define NBMU
//#define VMU

//#define RootSolver_with_Baryon
//#define EOS_with_baryon

//#define DYNAMICAL_SOURCE

//#define HydroPlus
//#define CRITICAL

/**************************************************************************************************************************************************/
/*********************************************************/
//Conservation laws//

#define NUMBER_CONSERVATION_LAWS 4

#ifndef NBMU
#define NUMBER_BARYON_COMPONENTS 0
#else
#define NUMBER_BARYON_COMPONENTS 1
#endif

/*********************************************************/
//Dissipative components//

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

#ifndef VMU
#define NUMBER_PROPAGATED_VMU_COMPONENTS 0
#else
#define NUMBER_PROPAGATED_VMU_COMPONENTS 4
#endif

/*********************************************************/
//HydroPlus extra modes//

#ifndef HydroPlus
#define NUMBER_SLOW_MODES 0
#else
#define NUMBER_SLOW_MODES 3
#endif

/**************************************************************************************************************************************************/
/*********************************************************/
//Conservation variables//

#define ALL_NUMBER_CONSERVATION_LAWS (NUMBER_CONSERVATION_LAWS+NUMBER_BARYON_COMPONENTS)// T^\tau\mu + N^\tau

#define NUMBER_DISSIPATIVE_CURRENTS (NUMBER_PI_COMPONENTS+NUMBER_PROPAGATED_PIMUNU_COMPONENTS)// \pi^\tau\mu + \Pi

#define ALL_NUMBER_DISSIPATIVE_CURRENTS (NUMBER_PI_COMPONENTS+NUMBER_PROPAGATED_VMU_COMPONENTS+NUMBER_PROPAGATED_PIMUNU_COMPONENTS)// \pi^\tau\mu + \Pi + V^\mu

#define NUMBER_CONSERVED_VARIABLES_NO_BULK (NUMBER_CONSERVATION_LAWS+NUMBER_PROPAGATED_PIMUNU_COMPONENTS)// T^\tau\mu + \pi^\tau\mu

#define NUMBER_CONSERVED_VARIABLES (NUMBER_CONSERVATION_LAWS+NUMBER_DISSIPATIVE_CURRENTS)// T^\tau\mu + \pi^\tau\mu + \Pi

#define ALL_NUMBER_CONSERVED_VARIABLES (ALL_NUMBER_CONSERVATION_LAWS+ALL_NUMBER_DISSIPATIVE_CURRENTS)// T^\tau\mu + N^\tau + \pi^\tau\mu + \Pi + V^\mu

#define NUMBER_ALL_EVOLVING_VARIABLES (ALL_NUMBER_CONSERVED_VARIABLES+NUMBER_SLOW_MODES)// T^\tau\mu + N^\tau + \pi^\tau\mu + \Pi + V^\mu + slow modes


/**************************************************************************************************************************************************/
/*********************************************************/
//Ideal or DNMR hydro//

#if ALL_NUMBER_DISSIPATIVE_CURRENTS==0
#define IDEAL
#endif

/*********************************************************/
//Precision control//

#define PRECISION double


/**************************************************************************************************************************************************/
/* conserved variables: energy momentum tensor, shear and bulk pressures, net baryon current, baryon diffusion current, slow modes
/**************************************************************************************************************************************************/

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
    PRECISION *phiQ[NUMBER_SLOW_MODES];// slow modes for Hydro+
} CONSERVED_VARIABLES;


/**************************************************************************************************************************************************/
/* flow velocity, Equation of State table, slow modes, dynamical source terms
/**************************************************************************************************************************************************/

// flow velocity
typedef struct
{
	PRECISION *ut;
	PRECISION *ux;
	PRECISION *uy;
	PRECISION *un;
} FLUID_VELOCITY;

// equation of state
typedef struct
{
    PRECISION *Pressure;
    PRECISION *Temperature;
    PRECISION *alphab;
    PRECISION *dpdrhob;
    PRECISION *dPdT;
    PRECISION *Chib;
} EQUATION_OF_STATE;

// baryon diffusion coefficients
typedef struct
{
    PRECISION *DB;
    PRECISION *sigmaB;
} BARYON_DIFFUSION_COEFF;

// slow modes
typedef struct
{
    PRECISION *phiQ[NUMBER_SLOW_MODES];
} SLOW_MODES;

// dynamical sources
typedef struct
{
    PRECISION *sourcet;
    PRECISION *sourcex;
    PRECISION *sourcey;
    PRECISION *sourcen;
    PRECISION *sourceb;
} DYNAMICAL_SOURCES;


/**************************************************************************************************************************************************/
/* instances
/**************************************************************************************************************************************************/

extern SLOW_MODES *eqPhiQ;  // Slow modes at equilibrium: updated, previous, intermediate values
extern PRECISION *Qvec; // Q vectors of slow modes
extern PRECISION *xieq; // correlation length

extern CONSERVED_VARIABLES *q,*Q,*qS;
extern FLUID_VELOCITY *u,*up,*uS;

extern PRECISION *e, *p, *seq;
extern PRECISION *rhob, *rhobp, *rhobS;
extern PRECISION *alphaB, *alphaBp, *alphaBS;
extern PRECISION *T, *Tp, *TS;

extern DYNAMICAL_SOURCES *Source;

extern EQUATION_OF_STATE *EOState;

extern BARYON_DIFFUSION_COEFF *BaryDiffCoeff;

/**************************************************************************************************************************************************/
/* initialization and update energy momentum tensor, net baryon current and slow modes
/**************************************************************************************************************************************************/

void setConservedVariables(double t, void * latticeParams);
void setCurrentConservedVariables();

/**************************************************************************************************************************************************/
/* update values for the "previous step
/**************************************************************************************************************************************************/

void swapFluidVelocity(FLUID_VELOCITY **arr1, FLUID_VELOCITY **arr2) ;
void swapPrimaryVariables(PRECISION **arr1, PRECISION **arr2);

/**************************************************************************************************************************************************/
/* Ghost cells
/**************************************************************************************************************************************************/

void setGhostCells(CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, void * latticeParams, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ);

void setGhostCellsKernelI(CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, void * latticeParams, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ);

void setGhostCellsKernelJ(CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, void * latticeParams, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ);

void setGhostCellsKernelK(CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u, void * latticeParams, PRECISION * const __restrict__ rhob, PRECISION * const __restrict__ alphaB, PRECISION * const __restrict__ T, PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ);

/**************************************************************************************************************************************************/
/* column index and memory
/**************************************************************************************************************************************************/

int columnMajorLinearIndex(int i, int j, int k, int nx, int ny);

void allocateHostMemory(int len);

void freeHostMemory();

#endif /* DYNAMICALVARIABLES_H_ */
