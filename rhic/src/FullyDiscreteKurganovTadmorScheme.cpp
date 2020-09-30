//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../include/FullyDiscreteKurganovTadmorScheme.h"
#include "../include/LatticeParameters.h"
#include "../include/DynamicalVariables.h"
#include "../include/SemiDiscreteKurganovTadmorScheme.h"
#include "../include/HalfSiteExtrapolation.h"
#include "../include/FluxFunctions.h"
#include "../include/SpectralRadius.h"
#include "../include/SourceTerms.h"
#include "../include/PrimaryVariables.h"
#include "../include/HydroParameters.h"
#include "../include/HydroPlus.h"
#include "../include/RegulationScheme.h"

#ifdef _OPENMP
#include <omp.h>
#endif


/**************************************************************************************************************************************************/
/* store nearest and next-to-nearest neighbors of all conserved variables at each grid to a 5 times longer array
/**************************************************************************************************************************************************/

void setNeighborCellsJK2(const PRECISION * const __restrict__ in, PRECISION * const __restrict__ out,
int s, int ptr, int smm, int sm, int sp, int spp)
{
	PRECISION data_ns= in[s];
	*(out + ptr    ) = in[smm];
	*(out + ptr + 1) = in[sm];
	*(out + ptr + 2) = data_ns;
	*(out + ptr + 3) = in[sp];
	*(out + ptr + 4) = in[spp];
}

/**************************************************************************************************************************************************/
/* source terms J times delta \tau for all components, excluding terms involving gradients of shear, bulk and baryon diffusion in J
/* contribution from t gradient of slow modes to tau component of baryon diffusion current included
/**************************************************************************************************************************************************/

void eulerStepKernelSource(PRECISION t, const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars, const PRECISION * const __restrict__ e, const PRECISION * const __restrict__ p, const FLUID_VELOCITY * const __restrict__ u, const FLUID_VELOCITY * const __restrict__ up, int ncx, int ncy, int ncz, PRECISION dt, PRECISION dx, PRECISION dy, PRECISION dz, const PRECISION * const __restrict__ rhob, const PRECISION * const __restrict__ rhobp, const PRECISION * const __restrict__ alphaB, const PRECISION * const __restrict__ alphaBp, const PRECISION * const __restrict__ T, const PRECISION * const __restrict__ Tp, const PRECISION * const __restrict__ seq, SLOW_MODES *  const __restrict__ eqPhiQ, void * hydroParams)
{
    
#pragma omp parallel for collapse(3)
#ifdef TILE
#pragma unroll_and_jam
#endif
    
	for(int k = 2; k < ncz-2; ++k) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int i = 2; i < ncx-2; ++i) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                
                //=====================================================================================
                // I: current variables
                //=====================================================================================
                
				PRECISION Q[NUMBER_ALL_EVOLVING_VARIABLES];
				
				Q[0] = currrentVars->ttt[s];
				Q[1] = currrentVars->ttx[s];
				Q[2] = currrentVars->tty[s];
				Q[3] = currrentVars->ttn[s];
#ifdef PIMUNU
				Q[4] = currrentVars->pitt[s];
				Q[5] = currrentVars->pitx[s];
				Q[6] = currrentVars->pity[s];
				Q[7] = currrentVars->pitn[s];
				Q[8] = currrentVars->pixx[s];
				Q[9] = currrentVars->pixy[s];
				Q[10] = currrentVars->pixn[s];
				Q[11] = currrentVars->piyy[s];
				Q[12] = currrentVars->piyn[s];
				Q[13] = currrentVars->pinn[s];
#endif
#ifdef PI
				Q[NUMBER_CONSERVED_VARIABLES_NO_BULK] = currrentVars->Pi[s];
#endif
#ifdef NBMU
                Q[NUMBER_CONSERVED_VARIABLES] = currrentVars->Nbt[s];
#endif
#ifdef VMU
                Q[NUMBER_CONSERVED_VARIABLES+1] = currrentVars->nbt[s];
                Q[NUMBER_CONSERVED_VARIABLES+2] = currrentVars->nbx[s];
                Q[NUMBER_CONSERVED_VARIABLES+3] = currrentVars->nby[s];
                Q[NUMBER_CONSERVED_VARIABLES+4] = currrentVars->nbn[s];
#endif
#ifdef HydroPlus
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    Q[ALL_NUMBER_CONSERVED_VARIABLES+n] = currrentVars->phiQ[n][s];
                }
#endif

                //=====================================================================================
                // II: source terms
                //=====================================================================================
                
                PRECISION S[NUMBER_ALL_EVOLVING_VARIABLES];
				loadSourceTerms2(Q, S, u, up->ut[s], up->ux[s], up->uy[s], up->un[s], t, e, p, s, ncx, ncy, ncz, dt, dx, dy, dz, Source, rhob, rhobp[s], alphaB, alphaBp[s], T, Tp[s], seq[s], eqPhiQ, hydroParams);
                
				PRECISION result[NUMBER_ALL_EVOLVING_VARIABLES];
                
#ifdef SIMD
#pragma omp simd
#endif
				for (unsigned int n = 0; n < NUMBER_ALL_EVOLVING_VARIABLES; ++n) {
					*(result+n) = *(Q+n) + dt * ( *(S+n) );
				}
                
                //=====================================================================================
                // III: update variables
                //=====================================================================================
                
				updatedVars->ttt[s] = result[0];
				updatedVars->ttx[s] = result[1];
				updatedVars->tty[s] = result[2];
				updatedVars->ttn[s] = result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] = result[4];
				updatedVars->pitx[s] = result[5];
				updatedVars->pity[s] = result[6];
				updatedVars->pitn[s] = result[7];
				updatedVars->pixx[s] = result[8];
				updatedVars->pixy[s] = result[9];
				updatedVars->pixn[s] = result[10];
				updatedVars->piyy[s] = result[11];
				updatedVars->piyn[s] = result[12];
				updatedVars->pinn[s] = result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] = result[NUMBER_CONSERVED_VARIABLES_NO_BULK];
#endif
#ifdef NBMU
                updatedVars->Nbt[s] = result[NUMBER_CONSERVED_VARIABLES];
#endif
#ifdef VMU
                updatedVars->nbt[s] = result[NUMBER_CONSERVED_VARIABLES+1];
                updatedVars->nbx[s] = result[NUMBER_CONSERVED_VARIABLES+2];
                updatedVars->nby[s] = result[NUMBER_CONSERVED_VARIABLES+3];
                updatedVars->nbn[s] = result[NUMBER_CONSERVED_VARIABLES+4];
#endif
#ifdef HydroPlus
                for (unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n) {
                    updatedVars->phiQ[n][s] = result[ALL_NUMBER_CONSERVED_VARIABLES+n];
                }
#endif
			}
		}
	}
}


/**************************************************************************************************************************************************/
/* (flux deltaHX/deltaX + J) times delta \tau for T^\tau^mu and N^tau, in J terms involving x gradient of shear, bulk and baryon diffusion
/* contribution from x gradient of slow modes to x component of baryon diffusion current included
/**************************************************************************************************************************************************/

void eulerStepKernelX(PRECISION t, const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars, const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ e, int ncx, int ncy, int ncz, PRECISION dt, PRECISION dx, const PRECISION * const __restrict__ rhob)
{
    
#pragma omp parallel for collapse(3)
#ifdef TILE
#pragma unroll_and_jam
#endif
    
    for(int k = 2; k < ncz-2; ++k) {
        for(int j = 2; j < ncy-2; ++j) {
            for(int i = 2; i < ncx-2; ++i) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                
                //=====================================================================================
                // I: current variables
                //=====================================================================================
                
				PRECISION I[5 * NUMBER_ALL_EVOLVING_VARIABLES];

				// calculate neighbor cell indices;
				int sim  = s-1;
				int simm = sim-1;
				int sip  = s+1;
				int sipp = sip+1;

				int ptr=0;
				setNeighborCellsJK2(currrentVars->ttt,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->tty,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
#ifdef PIMUNU
				setNeighborCellsJK2(currrentVars->pitt,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pity,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixy,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyy,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pinn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
#endif
#ifdef PI
				setNeighborCellsJK2(currrentVars->Pi,I,s,ptr,simm,sim,sip,sipp);  ptr+=5;
#endif
#ifdef NBMU
                setNeighborCellsJK2(currrentVars->Nbt,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
#endif
#ifdef VMU
                setNeighborCellsJK2(currrentVars->nbt,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbx,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nby,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbn,I,s,ptr,simm,sim,sip,sipp); ptr+=5;
#endif
#ifdef HydroPlus
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    setNeighborCellsJK2(currrentVars->phiQ[n],I,s,ptr,simm,sim,sip,sipp); ptr+=5;
                }
#endif
                
                //=====================================================================================
                // II: source terms
                //=====================================================================================
                
				PRECISION result[NUMBER_ALL_EVOLVING_VARIABLES];
                
                PRECISION Hforward[NUMBER_ALL_EVOLVING_VARIABLES];
				flux(I, Hforward, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusX, &Fx, t, e[s], rhob[s], u->ut[s]);
                
                PRECISION Hbackwards[NUMBER_ALL_EVOLVING_VARIABLES];
				flux(I, Hbackwards, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusX, &Fx, t, e[s], rhob[s], u->ut[s]);
                
                PRECISION Hx[NUMBER_ALL_EVOLVING_VARIABLES];
                loadSourceTermsX(I, Hx, u, s, dx);
                
#ifdef SIMD
#pragma omp simd
#endif
				for (unsigned int n = 0; n < NUMBER_ALL_EVOLVING_VARIABLES; ++n) {
                    *(result+n) = (*(Hbackwards+n) - *(Hforward+n))/dx + *(Hx+n);
                    *(result+n) *= dt;
				}
                
                
                //=====================================================================================
                // III: update variables
                //=====================================================================================
                
				updatedVars->ttt[s] += result[0];
				updatedVars->ttx[s] += result[1];
				updatedVars->tty[s] += result[2];
				updatedVars->ttn[s] += result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] += result[4];
				updatedVars->pitx[s] += result[5];
				updatedVars->pity[s] += result[6];
				updatedVars->pitn[s] += result[7];
				updatedVars->pixx[s] += result[8];
				updatedVars->pixy[s] += result[9];
				updatedVars->pixn[s] += result[10];
				updatedVars->piyy[s] += result[11];
				updatedVars->piyn[s] += result[12];
				updatedVars->pinn[s] += result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] += result[NUMBER_CONSERVED_VARIABLES_NO_BULK];
#endif
#ifdef NBMU
                updatedVars->Nbt[s] += result[NUMBER_CONSERVED_VARIABLES];
#endif
#ifdef VMU
                updatedVars->nbt[s] += result[NUMBER_CONSERVED_VARIABLES+1];
                updatedVars->nbx[s] += result[NUMBER_CONSERVED_VARIABLES+2];
                updatedVars->nby[s] += result[NUMBER_CONSERVED_VARIABLES+3];
                updatedVars->nbn[s] += result[NUMBER_CONSERVED_VARIABLES+4];
#endif
#ifdef HydroPlus
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    updatedVars->phiQ[n][s] += result[ALL_NUMBER_CONSERVED_VARIABLES+n];
                }
#endif
			}
		}
	}
}


/**************************************************************************************************************************************************/
/* (flux deltaHY/deltaY + J) times delta \tau for T^\tau^mu and N^tau, in J terms involving y gradient of shear, bulk and baryon diffusion
/* contribution from y gradient of slow modes to y component of baryon diffusion current included
/**************************************************************************************************************************************************/

void eulerStepKernelY(PRECISION t, const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars, const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ e, int ncx, int ncy, int ncz, PRECISION dt, PRECISION dy, const PRECISION * const __restrict__ rhob)
{
    
#pragma omp parallel for collapse(3)
#ifdef TILE
#pragma unroll_and_jam
#endif
    
    for(int k = 2; k < ncz-2; ++k) {
        for(int j = 2; j < ncy-2; ++j) {
            for(int i = 2; i < ncx-2; ++i) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                
                //=====================================================================================
                // I: current variables
                //=====================================================================================
                
				PRECISION J[5* NUMBER_ALL_EVOLVING_VARIABLES];
				PRECISION H[NUMBER_ALL_EVOLVING_VARIABLES];

				// calculate neighbor cell indices;
				int sjm = s-ncx;
				int sjmm = sjm-ncx;
				int sjp = s+ncx;
				int sjpp = sjp+ncx;

				int ptr=0;
				setNeighborCellsJK2(currrentVars->ttt,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->tty,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
#ifdef PIMUNU
				setNeighborCellsJK2(currrentVars->pitt,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pity,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixy,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyy,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pinn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
#endif
#ifdef PI
				setNeighborCellsJK2(currrentVars->Pi,J,s,ptr,sjmm,sjm,sjp,sjpp);  ptr+=5;
#endif
#ifdef NBMU
                setNeighborCellsJK2(currrentVars->Nbt,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
#endif
#ifdef VMU
                setNeighborCellsJK2(currrentVars->nbt,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbx,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nby,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbn,J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
#endif
#ifdef HydroPlus
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    setNeighborCellsJK2(currrentVars->phiQ[n],J,s,ptr,sjmm,sjm,sjp,sjpp); ptr+=5;
                }
#endif
                
                //=====================================================================================
                // II: source terms
                //=====================================================================================
                
				PRECISION result[NUMBER_ALL_EVOLVING_VARIABLES];
                
                PRECISION Hforward[NUMBER_ALL_EVOLVING_VARIABLES];
				flux(J, Hforward, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusY, &Fy, t, e[s], rhob[s], u->ut[s]);

                PRECISION Hbackwards[NUMBER_ALL_EVOLVING_VARIABLES];
				flux(J, Hbackwards, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusY, &Fy, t, e[s], rhob[s], u->ut[s]);

                PRECISION Hy[NUMBER_ALL_EVOLVING_VARIABLES];
				loadSourceTermsY(J, Hy, u, s, dy);

#ifdef SIMD
#pragma omp simd
#endif
                for (unsigned int n = 0; n < NUMBER_ALL_EVOLVING_VARIABLES; ++n) {
                    *(result+n) = (*(Hbackwards+n) - *(Hforward+n))/dy + *(Hy+n);
                    *(result+n) *= dt;
                }

                //=====================================================================================
                // III: update variables
                //=====================================================================================
                
				updatedVars->ttt[s] += result[0];
				updatedVars->ttx[s] += result[1];
				updatedVars->tty[s] += result[2];
				updatedVars->ttn[s] += result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] += result[4];
				updatedVars->pitx[s] += result[5];
				updatedVars->pity[s] += result[6];
				updatedVars->pitn[s] += result[7];
				updatedVars->pixx[s] += result[8];
				updatedVars->pixy[s] += result[9];
				updatedVars->pixn[s] += result[10];
				updatedVars->piyy[s] += result[11];
				updatedVars->piyn[s] += result[12];
				updatedVars->pinn[s] += result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] += result[NUMBER_CONSERVED_VARIABLES_NO_BULK];
#endif
#ifdef NBMU
                updatedVars->Nbt[s] += result[NUMBER_CONSERVED_VARIABLES];
#endif
#ifdef VMU
                updatedVars->nbt[s] += result[NUMBER_CONSERVED_VARIABLES+1];
                updatedVars->nbx[s] += result[NUMBER_CONSERVED_VARIABLES+2];
                updatedVars->nby[s] += result[NUMBER_CONSERVED_VARIABLES+3];
                updatedVars->nbn[s] += result[NUMBER_CONSERVED_VARIABLES+4];
#endif
#ifdef HydroPlus
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    updatedVars->phiQ[n][s] += result[ALL_NUMBER_CONSERVED_VARIABLES+n];
                }
#endif
			}
		}
	}
}


/**************************************************************************************************************************************************/
/* (flux deltaHZ/deltaZ + J) times delta \tau for T^\tau^mu and N^tau, in J terms involving eta_s gradient of shear, bulk and baryon diffusion
/* contribution from y gradient of slow modes to y component of baryon diffusion current included
/**************************************************************************************************************************************************/

void eulerStepKernelZ(PRECISION t, const CONSERVED_VARIABLES * const __restrict__ currrentVars, CONSERVED_VARIABLES * const __restrict__ updatedVars, const FLUID_VELOCITY * const __restrict__ u, const PRECISION * const __restrict__ e, int ncx, int ncy, int ncz, PRECISION dt, PRECISION dz, const PRECISION * const __restrict__ rhob)
{
    
#pragma omp parallel for collapse(3)
#ifdef TILE
#pragma unroll_and_jam
#endif
    
    for(int k = 2; k < ncz-2; ++k) {
        for(int j = 2; j < ncy-2; ++j) {
            for(int i = 2; i < ncx-2; ++i) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                
                //=====================================================================================
                // I: current variables
                //=====================================================================================
                
				PRECISION K[5 * NUMBER_ALL_EVOLVING_VARIABLES];
				PRECISION H[NUMBER_ALL_EVOLVING_VARIABLES];

				// calculate neighbor cell indices;
				int stride = ncx * ncy;
				int skm = s-stride;
				int skmm = skm-stride;
				int skp = s+stride;
				int skpp = skp+stride;

				int ptr=0;
				setNeighborCellsJK2(currrentVars->ttt,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->tty,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->ttn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#ifdef PIMUNU
				setNeighborCellsJK2(currrentVars->pitt,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pity,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pitn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixy,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pixn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyy,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->piyn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
				setNeighborCellsJK2(currrentVars->pinn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#endif
#ifdef PI
				setNeighborCellsJK2(currrentVars->Pi,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#endif
#ifdef NBMU
                setNeighborCellsJK2(currrentVars->Nbt,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#endif
#ifdef VMU
                setNeighborCellsJK2(currrentVars->nbt,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbx,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nby,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
                setNeighborCellsJK2(currrentVars->nbn,K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
#endif
#ifdef HydroPlus
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    setNeighborCellsJK2(currrentVars->phiQ[n],K,s,ptr,skmm,skm,skp,skpp); ptr+=5;
                }
#endif
                
                //=====================================================================================
                // II: source terms
                //=====================================================================================
                
				PRECISION result[NUMBER_ALL_EVOLVING_VARIABLES];
                
                PRECISION Hforward[NUMBER_ALL_EVOLVING_VARIABLES];
				flux(K, Hforward, &rightHalfCellExtrapolationForward, &leftHalfCellExtrapolationForward, &spectralRadiusZ, &Fz, t, e[s], rhob[s], u->ut[s]);

                PRECISION Hbackwards[NUMBER_ALL_EVOLVING_VARIABLES];
				flux(K, Hbackwards, &rightHalfCellExtrapolationBackwards, &leftHalfCellExtrapolationBackwards, &spectralRadiusZ, &Fz, t, e[s], rhob[s], u->ut[s]);

                PRECISION Hz[NUMBER_ALL_EVOLVING_VARIABLES];
				loadSourceTermsZ(K, Hz, u, s, t, dz);

#ifdef SIMD
#pragma omp simd
#endif
                for (unsigned int n = 0; n < NUMBER_ALL_EVOLVING_VARIABLES; ++n) {
                    *(result+n) = (*(Hbackwards+n) - *(Hforward+n))/dz + *(Hz+n);
                    *(result+n) *= dt;
                }
                
                //=====================================================================================
                // III: update variables
                //=====================================================================================
                
				updatedVars->ttt[s] += result[0];
				updatedVars->ttx[s] += result[1];
				updatedVars->tty[s] += result[2];
				updatedVars->ttn[s] += result[3];
#ifdef PIMUNU
				updatedVars->pitt[s] += result[4];
				updatedVars->pitx[s] += result[5];
				updatedVars->pity[s] += result[6];
				updatedVars->pitn[s] += result[7];
				updatedVars->pixx[s] += result[8];
				updatedVars->pixy[s] += result[9];
				updatedVars->pixn[s] += result[10];
				updatedVars->piyy[s] += result[11];
				updatedVars->piyn[s] += result[12];
				updatedVars->pinn[s] += result[13];
#endif
#ifdef PI
				updatedVars->Pi[s] += result[NUMBER_CONSERVED_VARIABLES_NO_BULK];
#endif
#ifdef NBMU
                updatedVars->Nbt[s] += result[NUMBER_CONSERVED_VARIABLES];
#endif
#ifdef VMU
                updatedVars->nbt[s] += result[NUMBER_CONSERVED_VARIABLES+1];
                updatedVars->nbx[s] += result[NUMBER_CONSERVED_VARIABLES+2];
                updatedVars->nby[s] += result[NUMBER_CONSERVED_VARIABLES+3];
                updatedVars->nbn[s] += result[NUMBER_CONSERVED_VARIABLES+4];
#endif
#ifdef HydroPlus
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    updatedVars->phiQ[n][s] += result[ALL_NUMBER_CONSERVED_VARIABLES+n];
                }
#endif
			}
		}
	}
}


/**************************************************************************************************************************************************/
/* two-step explicit Runge-Kutta for time integration
/**************************************************************************************************************************************************/

void convexCombinationEulerStepKernel(const CONSERVED_VARIABLES * const __restrict__ q, CONSERVED_VARIABLES * const __restrict__ Q, int ncx, int ncy, int ncz)
{
    
#pragma omp parallel for collapse(3)
#ifdef TILE
#pragma unroll_and_jam
#endif
    for(int k = 2; k < ncz-2; ++k) {
        for(int j = 2; j < ncy-2; ++j) {
            for(int i = 2; i < ncx-2; ++i) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
				Q->ttt[s] += q->ttt[s];
				Q->ttt[s] /= 2;
				Q->ttx[s] += q->ttx[s];
				Q->ttx[s] /= 2;
				Q->tty[s] += q->tty[s];
				Q->tty[s] /= 2;
				Q->ttn[s] += q->ttn[s];
				Q->ttn[s] /= 2;
#ifdef PIMUNU
				Q->pitt[s] += q->pitt[s];
				Q->pitt[s] /= 2;
				Q->pitx[s] += q->pitx[s];
				Q->pitx[s] /= 2;
				Q->pity[s] += q->pity[s];
				Q->pity[s] /= 2;
				Q->pitn[s] += q->pitn[s];
				Q->pitn[s] /= 2;
				Q->pixx[s] += q->pixx[s];
				Q->pixx[s] /= 2;
				Q->pixy[s] += q->pixy[s];
				Q->pixy[s] /= 2;
				Q->pixn[s] += q->pixn[s];
				Q->pixn[s] /= 2;
				Q->piyy[s] += q->piyy[s];
				Q->piyy[s] /= 2;
				Q->piyn[s] += q->piyn[s];
				Q->piyn[s] /= 2;
				Q->pinn[s] += q->pinn[s];
				Q->pinn[s] /= 2;
#endif
#ifdef PI
                Q->Pi[s]   += q->Pi[s];
                Q->Pi[s]   /= 2;
#endif
#ifdef NBMU
                Q->Nbt[s]  += q->Nbt[s];
                Q->Nbt[s]  /= 2;
#endif
#ifdef VMU
                Q->nbt[s]  += q->nbt[s];
                Q->nbt[s]  /= 2;
                Q->nbx[s]  += q->nbx[s];
                Q->nbx[s]  /= 2;
                Q->nby[s]  += q->nby[s];
                Q->nby[s]  /= 2;
                Q->nbn[s]  += q->nbn[s];
                Q->nbn[s]  /= 2;
#endif
#ifdef HydroPlus
                for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n){
                    Q->phiQ[n][s] += q->phiQ[n][s];
                    Q->phiQ[n][s] /= 2;
                }
#endif
            }
		}
	}
}


/**************************************************************************************************************************************************/
/* Kurganov-Tadmor algorithm for integration of hyperbolic equations
/**************************************************************************************************************************************************/

void rungeKutta2(PRECISION t, PRECISION dt, CONSERVED_VARIABLES * __restrict__ q, CONSERVED_VARIABLES * __restrict__ Q, void * latticeParams, void * hydroParams)
{
    
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
    
	int ncx = lattice->numComputationalLatticePointsX;
	int ncy = lattice->numComputationalLatticePointsY;
	int ncz = lattice->numComputationalLatticePointsRapidity;

	PRECISION dx = (PRECISION)(lattice->latticeSpacingX);
	PRECISION dy = (PRECISION)(lattice->latticeSpacingY);
	PRECISION dz = (PRECISION)(lattice->latticeSpacingRapidity);

	//===================================================
	// Euler STEP 1: predicted step
	//===================================================
    
    // Q: previous value, q: current value, qS: to be updated value. Q is added for time gradient of slow modes.
	eulerStepKernelSource(t, q, qS, e, p, u, up, ncx, ncy, ncz, dt, dx, dy, dz, rhob, rhobp, alphaB, alphaBp, T, Tp, seq, eqPhiQ, hydroParams);
    eulerStepKernelX(t, q, qS, u, e, ncx, ncy, ncz, dt, dx, rhob);
	eulerStepKernelY(t, q, qS, u, e, ncx, ncy, ncz, dt, dy, rhob);
    eulerStepKernelZ(t, q, qS, u, e, ncx, ncy, ncz, dt, dz, rhob);
    
	t+=dt;

    // calculate e, p and T etc. from the updated T^tau^mu and shear etc.
	setInferredVariablesKernel(qS, e, p, u, uS, t, latticeParams, rhobS, alphaBS, TS, seq, eqPhiQ);

#ifndef IDEAL
	regulateDissipativeCurrents(t, qS, e, p, rhobS, uS, ncx, ncy, ncz);
#endif

	setGhostCells(qS, e, p, uS, latticeParams, rhobS, alphaBS, TS, seq, eqPhiQ);

    
	//===================================================
	// Euler STEP 2: corrected step
	//===================================================
    
    // q: previous value, qS: current value, Q: to be updated value. q is added for time gradient slow modes
	eulerStepKernelSource(t, qS, Q, e, p, uS, u, ncx, ncy, ncz, dt, dx, dy, dz, rhobS, rhob, alphaBS, alphaB, TS, T, seq, eqPhiQ, hydroParams);
	eulerStepKernelX(t, qS, Q, uS, e, ncx, ncy, ncz, dt, dx, rhobS);
	eulerStepKernelY(t, qS, Q, uS, e, ncx, ncy, ncz, dt, dy, rhobS);
	eulerStepKernelZ(t, qS, Q, uS, e, ncx, ncy, ncz, dt, dz, rhobS);

	convexCombinationEulerStepKernel(q, Q, ncx, ncy, ncz);

    // u, alphaB etc will store the final updated value after setInferredVariablesKernel() and, in the next step, become the "current value".
	swapFluidVelocity(&up, &u);
    swapPrimaryVariables(&alphaBp, &alphaB);
    swapPrimaryVariables(&Tp, &T);
    swapPrimaryVariables(&rhobp, &rhob);
    
    // calculate e, p and T etc. from the updated T^tau^mu and shear etc. Q, e, p, u, rhob etc will store the final updated values
	setInferredVariablesKernel(Q, e, p, uS, u, t, latticeParams, rhob, alphaB, T, seq, eqPhiQ);
    
#ifndef IDEAL
	regulateDissipativeCurrents(t, Q, e, p, rhob, u, ncx, ncy, ncz);
#endif
    
	setGhostCells(Q, e, p, u, latticeParams, rhob, alphaB, T, seq, eqPhiQ);
    
}
