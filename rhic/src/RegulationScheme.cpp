//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../include/DynamicalVariables.h"

/**************************************************************************************************************************************************/
/* regulation of dissipative currents
/**************************************************************************************************************************************************/

// regulation method 1
const PRECISION xi0 = 0.1;
const PRECISION rhomax = 1.0;

// regulation method 2
const PRECISION chi0 = 10.0;
const PRECISION epsilon0 = 0.1;
const PRECISION reg_width = 0.01;
const PRECISION rmax_pi = 1.0;
const PRECISION rmax_Pi = 1.0;
const PRECISION rmax_n = 1.0;

#define Regulation_OLD

void regulateDissipativeCurrents(PRECISION t, const CONSERVED_VARIABLES * const __restrict__ currrentVars, const PRECISION * const __restrict__ e, const PRECISION * const __restrict__ p, const PRECISION * const __restrict__ rhob, const FLUID_VELOCITY * const __restrict__ u, int ncx, int ncy, int ncz) {
 
#pragma omp parallel for collapse(3)
#ifdef TILE
#pragma unroll_and_jam
#endif
    for(int k = 2; k < ncz-2; ++k) {
        for(int j = 2; j < ncy-2; ++j) {
            for(int i = 2; i < ncx-2; ++i) {
                int s = columnMajorLinearIndex(i, j, k, ncx, ncy);
#ifdef PIMUNU
                PRECISION pitt = currrentVars->pitt[s];
                PRECISION pitx = currrentVars->pitx[s];
                PRECISION pity = currrentVars->pity[s];
                PRECISION pitn = currrentVars->pitn[s];
                PRECISION pixx = currrentVars->pixx[s];
                PRECISION pixy = currrentVars->pixy[s];
                PRECISION pixn = currrentVars->pixn[s];
                PRECISION piyy = currrentVars->piyy[s];
                PRECISION piyn = currrentVars->piyn[s];
                PRECISION pinn = currrentVars->pinn[s];
#else
                PRECISION pitt = 0;
                PRECISION pitx = 0;
                PRECISION pity = 0;
                PRECISION pitn = 0;
                PRECISION pixx = 0;
                PRECISION pixy = 0;
                PRECISION pixn = 0;
                PRECISION piyy = 0;
                PRECISION piyn = 0;
                PRECISION pinn = 0;
#endif
#ifdef PI
                PRECISION Pi = currrentVars->Pi[s];
#else
                PRECISION Pi = 0;
#endif
#ifdef VMU
                PRECISION nbt = currrentVars->nbt[s];
                PRECISION nbx = currrentVars->nbx[s];
                PRECISION nby = currrentVars->nby[s];
                PRECISION nbn = currrentVars->nbn[s];
#else
                PRECISION nbt = 0;
                PRECISION nbx = 0;
                PRECISION nby = 0;
                PRECISION nbn = 0;
#endif
                
                PRECISION ut = u->ut[s];
                PRECISION ux = u->ux[s];
                PRECISION uy = u->uy[s];
                PRECISION un = u->un[s];
                
                PRECISION t2 = t*t;
                PRECISION normShear = sqrtf(e[s]*e[s]+3*p[s]*p[s]);
                PRECISION normBulk = normShear;
                PRECISION normBaryon = fabs(rhob[s]);
                
                PRECISION pipi = pitt*pitt-2*pitx*pitx-2*pity*pity+pixx*pixx+2*pixy*pixy+piyy*piyy-2*pitn*pitn*t2+2*pixn*pixn*t2+2*piyn*piyn*t2+pinn*pinn*t2*t2;
                PRECISION spipi = sqrt(fabs(pipi));
                if(isnan(spipi)) printf("found spipi Nan\n");
                
                PRECISION nbnb = nbt*nbt - nbx*nbx - nby*nby - nbn*nbn*t2;
                PRECISION snbnb = sqrt(fabs(nbnb));
                if(isnan(snbnb)) printf("found snbnb Nan\n");

                PRECISION facpi = 1.0;
                PRECISION facPi = 1.0;
                PRECISION facn  = 1.0;
                
                const double hbarc = 0.197326938;
                
                //===================================================
                // regulation factors
                //===================================================
#ifdef Regulation_OLD
                //***************************************************
                // method 1
#ifdef PIMUNU
                PRECISION pimumu = pitt - pixx - piyy - pinn*t2; // trace of shear stress
                PRECISION piu0 = -(pitn*t2*un) + pitt*ut - pitx*ux - pity*uy;
                PRECISION piu1 = -(pixn*t2*un) + pitx*ut - pixx*ux - pixy*uy;
                PRECISION piu2 = -(piyn*t2*un) + pity*ut - pixy*ux - piyy*uy;
                PRECISION piu3 = -(pinn*t2*un) + pitn*ut - pixn*ux - piyn*uy;
                
                PRECISION a1 = spipi/rhomax/normShear;
                PRECISION a2 = pimumu/xi0/rhomax/spipi;
                PRECISION a3 = piu0/xi0/rhomax/spipi;
                PRECISION a4 = piu1/xi0/rhomax/spipi;
                PRECISION a5 = piu2/xi0/rhomax/spipi;
                PRECISION a6 = piu3/xi0/rhomax/spipi;
                
                PRECISION a12 = fmax(a1,a2);
                PRECISION a34 = fmax(a3,a4);
                PRECISION a56 = fmax(a5,a6);
                PRECISION a3456 = fmax(a34,a56);
                PRECISION rho_pi = fmax(a12,a3456);
                
                if(fabs(rho_pi)>1.e-7) facpi = tanh(rho_pi)/rho_pi;
                if(isnan(facpi))   printf("found facpi Nan\n");
#endif
#ifdef PI
                PRECISION rho_Pi = 1.732*fabs(Pi)/rhomax/normBulk;

                if(fabs(rho_Pi)>1.e-7) facPi = tanh(rho_Pi)/rho_Pi;
                if(isnan(facPi))  printf("found facPi Nan\n");
#endif
#ifdef VMU
                PRECISION nbu = nbt*ut - nbx*ux - nby*uy - nbn*t2*un;
                PRECISION b1 = snbnb/rhomax/normBaryon;
                PRECISION b2 = nbu/xi0/rhomax/snbnb;
                PRECISION rho_n = fmax(b1,b2);

                if(fabs(rho_n)>1.e-7) facn = tanh(rho_n)/rho_n;
                if(isnan(facn)) printf("found facn Nan\n");
#endif
                //***************************************************
                // method 2
#else
                PRECISION factorShear = chi0 * (1./(exp(-(e[s]*hbarc - epsilon0)/reg_width) + 1.) - 1./(exp(epsilon0/reg_width) + 1.));
                PRECISION factorBulk = factorShear;
                PRECISION factorBaryon = factorShear;
                
                PRECISION r_pi = 1/factorShear * spipi/normShear;
                PRECISION r_Pi = 1/factorBulk * 1.732*fabs(Pi)/normBulk;
                PRECISION r_n = 1/factorBaryon * snbnb/normBaryon;
                
                if(r_pi>rmax_pi) facpi = rmax_pi/r_pi;
                if(isnan(facpi))   printf("found facpi Nan\n");
                
                if(r_Pi>rmax_Pi) facpi = rmax_Pi/r_Pi;
                if(isnan(facPi))   printf("found facPi Nan\n");
                
                if(r_n>rmax_n) facn = rmax_n/r_n;
                if(isnan(facn))   printf("found facn Nan\n");
#endif
                
                //===================================================
                // regulating dissipative components
                //===================================================
                
#ifdef PIMUNU
                currrentVars->pitt[s] *= facpi;
                currrentVars->pitx[s] *= facpi;
                currrentVars->pity[s] *= facpi;
                currrentVars->pitn[s] *= facpi;
                currrentVars->pixx[s] *= facpi;
                currrentVars->pixy[s] *= facpi;
                currrentVars->pixn[s] *= facpi;
                currrentVars->piyy[s] *= facpi;
                currrentVars->piyn[s] *= facpi;
                currrentVars->pinn[s] *= facpi;
#endif
#ifdef PI
                currrentVars->Pi[s] *= facPi;
#endif
#ifdef VMU
                currrentVars->nbt[s] *= facn;
                currrentVars->nbx[s] *= facn;
                currrentVars->nby[s] *= facn;
                currrentVars->nbn[s] *= facn;
#endif
            }
        }
    }
}
