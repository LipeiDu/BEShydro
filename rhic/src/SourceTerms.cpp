//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <stdio.h> // for printf
#include <math.h> // for math functions
#include <stdlib.h>

#include "../include/SourceTerms.h"
#include "../include/PrimaryVariables.h"
#include "../include/DynamicalVariables.h"
#include "../include/FullyDiscreteKurganovTadmorScheme.h" // for const params
#include "../include/EquationOfState.h" // for bulk terms
#include "../include/DynamicalSources.h"
#include "../include/FluxLimiter.h"
#include "../include/TransportCoefficients.h"
#include "../include/HydroPlus.h"

//#define USE_CARTESIAN_COORDINATES
//#define MINMOD_FOR_U_AND_P
//#define MINMOD_FOR_BARYON

#define HBARC 0.197326938

/**************************************************************************************************************************************************/
/* calculate source terms for dissipative compnents, e.g. shear, bulk and baryon diffusion etc.
/**************************************************************************************************************************************************/

void setDissipativeSourceTerms(PRECISION * const __restrict__ pimunuRHS, PRECISION * const __restrict__  piRHS, PRECISION * const __restrict__ nbmuRHS, PRECISION * const __restrict__ phiQRHS, PRECISION nbt, PRECISION nbx, PRECISION nby, PRECISION nbn, PRECISION rhob, PRECISION Nablat_rhob, PRECISION Nablax_rhob, PRECISION Nablay_rhob, PRECISION Nablan_rhob, PRECISION alphaB, PRECISION Nablat_alphaB, PRECISION Nablax_alphaB, PRECISION Nablay_alphaB, PRECISION Nablan_alphaB, PRECISION T, PRECISION Nablat_T, PRECISION Nablax_T, PRECISION Nablay_T, PRECISION Nablan_T, PRECISION seq, PRECISION t, PRECISION e, PRECISION p, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un, PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp, PRECISION pitt, PRECISION pitx, PRECISION pity, PRECISION pitn, PRECISION pixx, PRECISION pixy, PRECISION pixn, PRECISION piyy, PRECISION piyn, PRECISION pinn, PRECISION Pi, PRECISION dxut, PRECISION dyut, PRECISION dnut, PRECISION dxux, PRECISION dyux, PRECISION dnux, PRECISION dxuy, PRECISION dyuy, PRECISION dnuy, PRECISION dxun, PRECISION dyun, PRECISION dnun, PRECISION dkvk, PRECISION d_dt, const PRECISION * const __restrict__ PhiQ, const PRECISION * const __restrict__ equiPhiQ, void * hydroParams)
{
    
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    
    PRECISION etabar = (PRECISION)(hydro->shearViscosityToEntropyDensity);
    PRECISION Cb = (PRECISION)(hydro->cB);
    int kappaType = hydro->kappaType;
    int gradientType = hydro->gradientType;
    int criticalSlowingDown = hydro->criticalSlowingDown;
    
    /*******************************************************************************/
    /* Sec. I basic elements
    /*******************************************************************************/
    
    PRECISION t2 = t * t;
    PRECISION t3 = t * t2;
    
	//*********************************************************\
	//* transport coefficients
	//*********************************************************/
    
    // shear stress tensor
    PRECISION taupiInv = T / 5  / etabar;//0.6 * rhob;// //d_etabar is eta/(e+p)
    PRECISION beta_pi = (e + p) / 5;//0.8 * T * rhob;//
    
	// bulk transport coefficients
	PRECISION cs2 = speedOfSoundSquared(e, rhob);
	PRECISION a = 1.0/3.0 - cs2;
	PRECISION a2 = a*a;
	PRECISION beta_Pi = 15*a2*(e+p);
	PRECISION lambda_Pipi = 8*a/5;
    PRECISION zetabar = bulkViscosityToEntropyDensity(T);
    PRECISION tauPiInv = 15*a2*T/zetabar;
    
    // baryon diffusion
    PRECISION tau_n = Cb/T;//9.0/4.0/rhob;//
    PRECISION taunInv = 1/tau_n;
    PRECISION kappaB = 0.0;
    
    switch (kappaType) {
        case 0: {
            kappaB = baryonDiffusionCoefficientKinetic(Cb, T, rhob, alphaB, e, p);
            break;
        }
        case 1: {
            PRECISION diffusionCoeff[2];
            baryonDiffusionCoefficient(T, T*alphaB, diffusionCoeff);
            
            kappaB = diffusionCoeff[0];
            break;
        }
        case 2: {
            kappaB = baryonDiffusionCoefficientAdscft(T, rhob, alphaB, e, p, seq);
            break;
        }
        case 3: {
            kappaB = baryonDiffusionCoefficientHydroPlus(T, rhob, alphaB, e, p, seq);
            break;
        }
        case 4: {
            kappaB = baryonDiffusionCoefficientTest(T, rhob, alphaB);
            break;
        }
        default: {
            exit(-1);
        }
    }
    
    //*********************************************************\
    //* flow velocity derivatives
    //*********************************************************/

	PRECISION ut2 = ut * ut;
	PRECISION un2 = un * un;
    
	// time derivatives of u
	PRECISION dtut = (ut - utp) / d_dt;
	PRECISION dtux = (ux - uxp) / d_dt;
	PRECISION dtuy = (uy - uyp) / d_dt;
	PRECISION dtun = (un - unp) / d_dt;

	// covariant derivatives

	PRECISION Dut = ut*dtut + ux*dxut + uy*dyut + un*dnut + t*un*un;
	PRECISION DuxUpper = ut*dtux + ux*dxux + uy*dyux + un*dnux;
	PRECISION Dux = -DuxUpper;
	PRECISION DuyUpper = ut*dtuy + ux*dxuy + uy*dyuy + un*dnuy;
	PRECISION Duy = -DuyUpper;
	PRECISION DunUpper = ut*dtun + ux*dxun + uy*dyun + un*dnun + 2*ut*un/t;
	PRECISION Dun = -t2*DunUpper;

	PRECISION dut = Dut - t*un*un;
	PRECISION dux = ut*dtux + ux*dxux + uy*dyux + un*dnux;
	PRECISION duy = ut*dtuy + ux*dxuy + uy*dyuy + un*dnuy;
	PRECISION dun = ut*dtun + ux*dxun + uy*dyun + un*dnun;

	//*********************************************************\
	//* expansion rate
	//*********************************************************/
	PRECISION theta = ut / t + dtut + dxux + dyuy + dnun;

	//*********************************************************\
	//* shear tensor
	//*********************************************************/
	PRECISION stt = -t * ut * un2 + (dtut - ut * dut) + (ut2 - 1) * theta / 3;
	PRECISION stx = -(t * un2 * ux) / 2 + (dtux - dxut) / 2 - (ux * dut + ut * dux) / 2 + ut * ux * theta / 3;
	PRECISION sty = -(t * un2 * uy) / 2 + (dtuy - dyut) / 2 - (uy * dut + ut * duy) / 2 + ut * uy * theta / 3;
	PRECISION stn = -un * (2 * ut2 + t2 * un2) / (2 * t) + (dtun - dnut / t2) / 2 - (un * dut + ut * dun) / 2 + ut * un * theta / 3;
	PRECISION sxx = -(dxux + ux * dux) + (1 + ux*ux) * theta / 3;
	PRECISION sxy = -(dxuy + dyux) / 2 - (uy * dux + ux * duy) / 2	+ ux * uy * theta / 3;
	PRECISION sxn = -ut * ux * un / t - (dxun + dnux / t2) / 2 - (un * dux + ux * dun) / 2 + ux * un * theta / 3;
	PRECISION syy = -(dyuy + uy * duy) + (1 + uy*uy) * theta / 3;
	PRECISION syn = -ut * uy * un / t - (dyun + dnuy / t2) / 2 - (un * duy + uy * dun) / 2 + uy * un * theta / 3;
	PRECISION snn = -ut * (1 + 2 * t2 * un2) / t3 - dnun / t2 - un * dun + (1 / t2 + un2) * theta / 3;

	//*********************************************************\
	//* vorticity tensor
	//*********************************************************/
    PRECISION wtx = (-dtux  -  dxut  +  ut * dux  -  ux * dut  -  t * un2 * ux) / 2.0;
    PRECISION wty = (-dtuy  -  dyut  +  ut * duy  -  uy * dut  -  t * un2 * uy) / 2.0;
    PRECISION wtn = (-t2 * dtun  -  2.0 * t * un  - dnut  -  ut * Dun  - t2 * un * dut  -  t3 * un * un2) / 2.0;
    PRECISION wxy = (dxuy  -  dyux  +  ux * duy  -  uy * dux) / 2.0;
    PRECISION wxn = (t2 * dxun  -  dnux  -  ux * Dun  -  t2 * un * dux) / 2.0;
    PRECISION wyn = (t2 * dyun  -  dnuy  -  uy * Dun  -  t2 * un * duy) / 2.0;
	// anti-symmetric vorticity components
	PRECISION wxt = wtx;
	PRECISION wyt = wty;
	PRECISION wnt = wtn / t2;
	PRECISION wyx = -wxy;
	PRECISION wnx = -wxn / t2;
	PRECISION wny = -wyn / t2;
    
    //*********************************************************\
    //* pi^\mu\nu \sigma_\mu\nu
    //*********************************************************/
    PRECISION ux2 = ux*ux;
    PRECISION uy2 = uy*uy;
    PRECISION ps  = pitt*stt-2*pitx*stx-2*pity*sty+pixx*sxx+2*pixy*sxy+piyy*syy - 2*pitn*stn*t2+2*pixn*sxn*t2+2*piyn*syn*t2+pinn*snn*t2*t2;

    /*******************************************************************************/
    /* Sec. II source terms on RHS
    /*******************************************************************************/
    
#ifdef PIMUNU
	//*********************************************************\
	//* I1
	//*********************************************************/
	PRECISION I1tt = 2 * ut * (pitt * Dut + pitx * Dux + pity * Duy + pitn * Dun);
	PRECISION I1tx = (pitt * ux + pitx * ut) * Dut + (pitx * ux + pixx * ut) * Dux + (pity * ux + pixy * ut) * Duy + (pitn * ux + pixn * ut) * Dun;
	PRECISION I1ty = (pitt * uy + pity * ut) * Dut + (pitx * uy + pixy * ut) * Dux + (pity * uy + piyy * ut) * Duy + (pitn * uy + piyn * ut) * Dun;
	PRECISION I1tn = (pitt * un + pitn * ut) * Dut + (pitx * un + pixn * ut) * Dux + (pity * un + piyn * ut) * Duy + (pitn * un + pinn * ut) * Dun;
	PRECISION I1xx = 2 * ux * (pitx * Dut + pixx * Dux + pixy * Duy + pixn * Dun);
	PRECISION I1xy = (pitx * uy + pity * ux) * Dut + (pixx * uy + pixy * ux) * Dux + (pixy * uy + piyy * ux) * Duy + (pixn * uy + piyn * ux) * Dun;
	PRECISION I1xn = (pitx * un + pitn * ux) * Dut + (pixx * un + pixn * ux) * Dux + (pixy * un + piyn * ux) * Duy + (pixn * un + pinn * ux) * Dun;
	PRECISION I1yy = 2 * uy * (pity * Dut + pixy * Dux + piyy * Duy + piyn * Dun);
	PRECISION I1yn = (pity * un + pitn * uy) * Dut + (pixy * un + pixn * uy) * Dux + (piyy * un + piyn * uy) * Duy + (piyn * un + pinn * uy) * Dun;
	PRECISION I1nn = 2 * un * (pitn * Dut + pixn * Dux + piyn * Duy + pinn * Dun);

	//*********************************************************\
	//* I2
	//*********************************************************/
	PRECISION I2tt = theta * pitt;
	PRECISION I2tx = theta * pitx;
	PRECISION I2ty = theta * pity;
	PRECISION I2tn = theta * pitn;
	PRECISION I2xx = theta * pixx;
	PRECISION I2xy = theta * pixy;
	PRECISION I2xn = theta * pixn;
	PRECISION I2yy = theta * piyy;
	PRECISION I2yn = theta * piyn;
	PRECISION I2nn = theta * pinn;

	//*********************************************************\
	//* I3
	//*********************************************************/
	PRECISION I3tt = 2 * (pitx * wtx + pity * wty + pitn * wtn);
	PRECISION I3tx = pitt * wxt + pity * wxy + pitn * wxn + pixx * wtx + pixy * wty + pixn * wtn;
	PRECISION I3ty = pitt * wyt + pitx * wyx + pitn * wyn + pixy * wtx + piyy * wty + piyn * wtn;
	PRECISION I3tn = pitt * wnt + pitx * wnx + pity * wny + pixn * wtx + piyn * wty + pinn * wtn;
	PRECISION I3xx = 2 * (pitx * wxt + pixy * wxy + pixn * wxn);
	PRECISION I3xy = pitx * wyt + pity * wxt + pixx * wyx + piyy * wxy + pixn * wyn + piyn * wxn;
	PRECISION I3xn = pitx * wnt + pitn * wxt + pixx * wnx + pixy * wny + piyn * wxy + pinn * wxn;
	PRECISION I3yy = 2 * (pity * wyt + pixy * wyx + piyn * wyn);
	PRECISION I3yn = pity * wnt + pitn * wyt + pixy * wnx + pixn * wyx + piyy * wny + pinn * wyn;
	PRECISION I3nn = 2 * (pitn * wnt + pixn * wnx + piyn * wny);

	//*********************************************************\
	//* I4
	//*********************************************************/

	PRECISION I4tt = (pitt * stt - pitx * stx  - pity * sty - t2 * pitn * stn) - (1 - ut2) * ps / 3;
	PRECISION I4tx = (pitt * stx + pitx * stt) / 2 - (pitx * sxx + pixx * stx) / 2 - (pity * sxy + pixy * sty) / 2
                    - t2 * (pitn * sxn + pixn * stn) / 2 + (ut * ux) * ps / 3;
	PRECISION I4ty = (pitt * sty + pity * stt) / 2 - (pitx * sxy + pixy * stx) / 2 - (pity * syy + piyy * sty) / 2
                    - t2 * (pitn * syn + piyn * stn) / 2 + (ut * uy) * ps / 3;
	PRECISION I4tn = (pitt * stn + pitn * stt) / 2 - (pitx * sxn + pixn * stx) / 2 - (pity * syn + piyn * sty) / 2
                    - t2 * (pitn * snn + pinn * stn) / 2 + (ut * un) * ps / 3;
	PRECISION I4xx = (pitx * stx - pixx * sxx  - pixy * sxy - t2 * pixn * sxn) + (1 + ux2) * ps / 3;
	PRECISION I4xy = (pitx * sty + pity * stx) / 2 - (pixx * sxy + pixy * sxx) / 2 - (pixy * syy + piyy * sxy) / 2
                    - t2 * (pixn * syn + piyn * sxn) / 2 + (ux * uy) * ps / 3;
	PRECISION I4xn = (pitx * stn + pitn * stx) / 2 - (pixx * sxn + pixn * sxx) / 2 - (pixy * syn + piyn * sxy) / 2
                    - t2 * (pixn * snn + pinn * sxn) / 2 + (ux * un) * ps / 3;
	PRECISION I4yy = (pity * sty - pixy * sxy  - piyy * syy - t2 * piyn * syn) + (1 + uy2) * ps / 3;
	PRECISION I4yn = (pity * stn + pitn * sty) / 2 - (pixy * sxn + pixn * sxy) / 2 - (piyy * syn + piyn * syy) / 2
                    - t2 * (piyn * snn + pinn * syn) / 2 + (uy * un) * ps / 3;
	PRECISION I4nn = (pitn * stn - pixn * sxn  - piyn * syn - t2 * pinn * snn) + (1 / t2 + un2) * ps / 3;

	//*********************************************************\
	//* I_pi
	//*********************************************************/
	PRECISION Itt = I1tt + delta_pipi * I2tt - tau_piw * I3tt + tau_pipi * I4tt - lambda_piPi * Pi * stt;
	PRECISION Itx = I1tx + delta_pipi * I2tx - tau_piw * I3tx + tau_pipi * I4tx - lambda_piPi * Pi * stx;
	PRECISION Ity = I1ty + delta_pipi * I2ty - tau_piw * I3ty + tau_pipi * I4ty - lambda_piPi * Pi * sty;
	PRECISION Itn = I1tn + delta_pipi * I2tn - tau_piw * I3tn + tau_pipi * I4tn - lambda_piPi * Pi * stn;
	PRECISION Ixx = I1xx + delta_pipi * I2xx - tau_piw * I3xx + tau_pipi * I4xx - lambda_piPi * Pi * sxx;
	PRECISION Ixy = I1xy + delta_pipi * I2xy - tau_piw * I3xy + tau_pipi * I4xy - lambda_piPi * Pi * sxy;
	PRECISION Ixn = I1xn + delta_pipi * I2xn - tau_piw * I3xn + tau_pipi * I4xn - lambda_piPi * Pi * sxn;
	PRECISION Iyy = I1yy + delta_pipi * I2yy - tau_piw * I3yy + tau_pipi * I4yy - lambda_piPi * Pi * syy;
	PRECISION Iyn = I1yn + delta_pipi * I2yn - tau_piw * I3yn + tau_pipi * I4yn - lambda_piPi * Pi * syn;
	PRECISION Inn = I1nn + delta_pipi * I2nn - tau_piw * I3nn + tau_pipi * I4nn - lambda_piPi * Pi * snn;

	//*********************************************************\
	//* shear stress tensor source terms, i.e. terms on RHS
	//*********************************************************/
	PRECISION dpitt = 2 * beta_pi * stt - pitt * taupiInv - Itt - 2 * un * t * pitn;
	PRECISION dpitx = 2 * beta_pi * stx - pitx * taupiInv - Itx - un * t * pixn;
	PRECISION dpity = 2 * beta_pi * sty - pity * taupiInv - Ity - un * t * piyn;
	PRECISION dpitn = 2 * beta_pi * stn - pitn * taupiInv - Itn - un * t * pinn - (ut * pitn + un * pitt) / t;
	PRECISION dpixx = 2 * beta_pi * sxx - pixx * taupiInv - Ixx;
	PRECISION dpixy = 2 * beta_pi * sxy - pixy * taupiInv - Ixy;
	PRECISION dpixn = 2 * beta_pi * sxn - pixn * taupiInv - Ixn - (ut * pixn + un * pitx) / t;
	PRECISION dpiyy = 2 * beta_pi * syy - piyy * taupiInv - Iyy;
	PRECISION dpiyn = 2 * beta_pi * syn - piyn * taupiInv - Iyn - (ut * piyn + un * pity) / t;
	PRECISION dpinn = 2 * beta_pi * snn - pinn * taupiInv - Inn - 2 * (ut * pinn + un * pitn) / t;
#endif
#ifdef PI
    //*********************************************************\
    //* I_Pi
    //*********************************************************/
    PRECISION I_Pi = delta_PiPi*Pi*theta - lambda_Pipi*ps;
    
	//*********************************************************\
	//* bulk viscous pressure source terms, i.e. terms on RHS
	//*********************************************************/
	PRECISION dPi = -beta_Pi*theta - Pi*tauPiInv - I_Pi;
#endif

    // correlation length returns 1.0 if CRITICAL is not defined.
    PRECISION muB = T * alphaB;
    PRECISION corrL = correlationLength(T, muB);

#ifdef VMU
    
    //*********************************************************\
    //* I1, I2, I3, I4
    //*********************************************************/
    
    PRECISION nbD  = nbt * Dut + nbx * Dux + nby * Duy + nbn * Dun;
    PRECISION I1t = ut * nbD;
    PRECISION I1x = ux * nbD;
    PRECISION I1y = uy * nbD;
    PRECISION I1n = un * nbD;
    
    PRECISION I2t = nbt * theta;
    PRECISION I2x = nbx * theta;
    PRECISION I2y = nby * theta;
    PRECISION I2n = nbn * theta;
    
    PRECISION I3t = -nbx * wxt - nby * wyt - t2 * nbn * wnt;
    PRECISION I3x = -nbt * wtx + nby * wyx + t2 * nbn * wnx;
    PRECISION I3y = -nbt * wty + nbx * wxy + t2 * nbn * wny;
    PRECISION I3n = 1/t2 * (-nbt * wtn + nbx * wxn + nby * wyn);
    
    PRECISION I4t = stt * nbt - stx * nbx - sty * nby - t2 * stn * nbn;
    PRECISION I4x = stx * nbt - sxx * nbx - sxy * nby - t2 * sxn * nbn;
    PRECISION I4y = sty * nbt - sxy * nbx - syy * nby - t2 * syn * nbn;
    PRECISION I4n = stn * nbt - sxn * nbx - syn * nby - t2 * snn * nbn;
    
    //*********************************************************\
    //* In
    //*********************************************************/
    
    PRECISION It = I1t + delta_nn * I2t + tau_nw * I3t + lambda_nn * I4t;
    PRECISION Ix = I1x + delta_nn * I2x + tau_nw * I3x + lambda_nn * I4x;
    PRECISION Iy = I1y + delta_nn * I2y + tau_nw * I3y + lambda_nn * I4y;
    PRECISION In = I1n + delta_nn * I2n + tau_nw * I3n + lambda_nn * I4n;
    
    //*********************************************************\
    //* baryon diffusion source terms, i.e. terms on RHS
    //*********************************************************/
    
    PRECISION chib = chiB(e, rhob+1.e-5);
    
    if(criticalSlowingDown){
        chib = chib * corrL * corrL;
        kappaB = kappaB * corrL;
        taunInv = taunInv / corrL;
    }
    
    switch (gradientType) {
        case 0: { // gradient of mu/T
            break;
        }
        case 1: { // gradient of rhob + T
            PRECISION facn = 1/(T+1.e-5)/(chib/pow(HBARC,2)+1.e-5);
            PRECISION facT = 1/(T+1.e-5)/(rhob+1.e-5) * (dPdT(e, rhob) - (e+p)/T);
            
            Nablat_alphaB = facn * Nablat_rhob + facT * Nablat_T;
            Nablax_alphaB = facn * Nablax_rhob + facT * Nablax_T;
            Nablay_alphaB = facn * Nablay_rhob + facT * Nablay_T;
            Nablan_alphaB = facn * Nablan_rhob + facT * Nablan_T;
            break;
        }
        case 2: { // gradient of rhob
            PRECISION facn = 1/(T+1.e-5)/(chib/pow(HBARC,2)+1.e-5);
            
            Nablat_alphaB = facn * Nablat_rhob;
            Nablax_alphaB = facn * Nablax_rhob;
            Nablay_alphaB = facn * Nablay_rhob;
            Nablan_alphaB = facn * Nablan_rhob;
            break;
        }
        case 3: { // gradeint of T
            PRECISION facT = 1/(T+1.e-5)/(rhob+1.e-5) * (dPdT(e, rhob) - (e+p)/T);
            
            Nablat_alphaB = facT * Nablat_T;
            Nablax_alphaB = facT * Nablax_T;
            Nablay_alphaB = facT * Nablay_T;
            Nablan_alphaB = facT * Nablan_T;
            break;
        }
        default: {
            exit(-1);
        }
    }
    
    PRECISION dnt = taunInv * kappaB * Nablat_alphaB - taunInv * nbt - It - t * un * nbn;
    PRECISION dnx = taunInv * kappaB * Nablax_alphaB - taunInv * nbx - Ix;
    PRECISION dny = taunInv * kappaB * Nablay_alphaB - taunInv * nby - Iy;
    PRECISION dnn = taunInv * kappaB * Nablan_alphaB - taunInv * nbn - In - 1/t * (ut * nbn + un * nbt) ;
    
#endif
#ifdef HydroPlus
    //*********************************************************\
    //* for the slow modes from Hydro+
    //*********************************************************/

    PRECISION corrL2 = corrL * corrL;
    PRECISION gammaPhi = relaxationCoefficientPhi(rhob, seq, T, corrL2);
#endif
    
    
    /*******************************************************************************/
    /* Sec. III load the sources
    /*******************************************************************************/
    
	//*********************************************************\
	//* time derivative of the dissipative quantities
	//*********************************************************/
#ifdef PIMUNU
	pimunuRHS[0] = dpitt / ut + pitt * dkvk;
	pimunuRHS[1] = dpitx / ut + pitx * dkvk;
	pimunuRHS[2] = dpity / ut + pity * dkvk;
	pimunuRHS[3] = dpitn / ut + pitn * dkvk;
	pimunuRHS[4] = dpixx / ut + pixx * dkvk;
	pimunuRHS[5] = dpixy / ut + pixy * dkvk;
	pimunuRHS[6] = dpixn / ut + pixn * dkvk;
	pimunuRHS[7] = dpiyy / ut + piyy * dkvk;
	pimunuRHS[8] = dpiyn / ut + piyn * dkvk;
	pimunuRHS[9] = dpinn / ut + pinn * dkvk;
#endif
#ifdef PI
	*piRHS = dPi / ut + Pi * dkvk;
#endif
#ifdef VMU
    nbmuRHS[0] = dnt / ut + nbt * dkvk;
    nbmuRHS[1] = dnx / ut + nbx * dkvk;
    nbmuRHS[2] = dny / ut + nby * dkvk;
    nbmuRHS[3] = dnn / ut + nbn * dkvk;
#endif
#ifdef HydroPlus
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
    {
        PRECISION gammaQ = relaxationCoefficientPhiQ(gammaPhi, corrL2, Qvec[n]);
        phiQRHS[n] = - 1/ut * gammaQ * (PhiQ[n] - equiPhiQ[n]) + PhiQ[n] * dkvk;
    }
#endif
}


/**************************************************************************************************************************************************/
/* calculate source terms for T^tau^mu and N^tau, only include source terms directly dependent on x gradient of shear, bulk and baryon diffusion
/**************************************************************************************************************************************************/

void loadSourceTermsX(const PRECISION * const __restrict__ I, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s,
PRECISION d_dx
) {
    for (unsigned int n = 0; n < NUMBER_ALL_EVOLVING_VARIABLES; ++n) {
        S[n] = 0.0;
    }
    
#ifndef IDEAL
	//=========================================================
	// spatial derivatives of the conserved variables \pi^{\mu\nu}
	//=========================================================
	PRECISION facX = 1/d_dx/2;
	int ptr = 20; // 5 * n (with n = 4 corresponding to pitt)
#ifdef PIMUNU
	PRECISION dxpitt = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpitx = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpity = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpitn = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpixx = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpixy = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;	ptr+=5;
	PRECISION dxpixn = (*(I + ptr + 3) - *(I + ptr + 1)) *facX; ptr+=20;
#endif
#ifdef PI
    PRECISION dxPi = (*(I + ptr + 3) - *(I + ptr + 1)) *facX; ptr+=5;
#endif
    //=========================================================
    // spatial derivatives of the conserved variables \nb^{\mu}
    //=========================================================
#ifdef VMU
    ptr+=5; // increased by 5 to start with nbt
    PRECISION dxnbt = (*(I + ptr + 3) - *(I + ptr + 1)) *facX; ptr+=5;
    PRECISION dxnbx = (*(I + ptr + 3) - *(I + ptr + 1)) *facX;
#endif

	PRECISION ut = u->ut[s];
	PRECISION ux = u->ux[s];

	//=========================================================
	// set dx terms in the source terms of T^\tau\mu
	//=========================================================
	PRECISION vx = ux / ut;
#ifdef PIMUNU
#ifndef PI
	S[0] = dxpitt*vx - dxpitx;
	S[1] = dxpitx*vx - dxpixx;
#else
	S[0] = dxpitt*vx - dxpitx - vx*dxPi;
	S[1] = dxpitx*vx - dxpixx - dxPi;
#endif
	S[2] = dxpity*vx - dxpixy;
	S[3] = dxpitn*vx - dxpixn;
#else
    S[0] = 0.0;
    S[1] = 0.0;
    S[2] = 0.0;
    S[3] = 0.0;
#endif
    //=========================================================
    // set dx terms in the source terms of baryon current
    //=========================================================
#ifdef NBMU
#ifdef VMU
    S[NUMBER_CONSERVED_VARIABLES] = dxnbt*vx - dxnbx;
#else
    S[NUMBER_CONSERVED_VARIABLES] = 0.0;
#endif
#endif
#endif
}


/**************************************************************************************************************************************************/
/* calculate source terms for T^tau^mu and N^tau, only include source terms directly dependent on y gradient of shear, bulk and baryon diffusion
/**************************************************************************************************************************************************/

void loadSourceTermsY(const PRECISION * const __restrict__ J, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s,
PRECISION d_dy
) {
    for (unsigned int n = 0; n < NUMBER_ALL_EVOLVING_VARIABLES; ++n) {
        S[n] = 0.0;
    }
    
#ifndef IDEAL
	//=========================================================
	// spatial derivatives of the conserved variables \pi^{\mu\nu}
	//=========================================================
	PRECISION facY = 1/d_dy/2;
	int ptr = 20; // 5 * n (with n = 4 corresponding to pitt)
#ifdef PIMUNU
	PRECISION dypitt = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=5;
	PRECISION dypitx = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=5;
	PRECISION dypity = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=5;
	PRECISION dypitn = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=10;
	PRECISION dypixy = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=10;
	PRECISION dypiyy = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;	ptr+=5;
	PRECISION dypiyn = (*(J + ptr + 3) - *(J + ptr + 1)) *facY; ptr+=10;
#endif
#ifdef PI
    PRECISION dyPi = (*(J + ptr + 3) - *(J + ptr + 1)) *facY; ptr+=5;
#endif
    //=========================================================
    // spatial derivatives of the conserved variables \nb^{\mu}
    //=========================================================
#ifdef VMU
    ptr+=5;
    PRECISION dynbt = (*(J + ptr + 3) - *(J + ptr + 1)) *facY; ptr+=10;
    PRECISION dynby = (*(J + ptr + 3) - *(J + ptr + 1)) *facY;
#endif
    

	PRECISION ut = u->ut[s];
	PRECISION uy = u->uy[s];

	//=========================================================
	// set dy terms in the source terms
	//=========================================================
	PRECISION vy = uy / ut;
#ifdef PIMUNU
#ifndef PI
	S[0] = dypitt*vy - dypity;
	S[2] = dypity*vy - dypiyy;
#else
	S[0] = dypitt*vy - dypity - vy*dyPi;
	S[2] = dypity*vy - dypiyy - dyPi;
#endif
	S[1] = dypitx*vy - dypixy;
	S[3] = dypitn*vy - dypiyn;
#else
    S[0] = 0.0;
    S[1] = 0.0;
    S[2] = 0.0;
    S[3] = 0.0;
#endif
    //=========================================================
    // set dy terms in the source terms of baryon current
    //=========================================================
#ifdef NBMU
#ifdef VMU
    S[NUMBER_CONSERVED_VARIABLES] = dynbt*vy - dynby;
#else
    S[NUMBER_CONSERVED_VARIABLES] = 0.0;
#endif
#endif
#endif
}


/**************************************************************************************************************************************************/
/* calculate source terms for T^tau^mu and N^tau, only include source terms directly dependent on eta_s gradient of shear, bulk and baryon diffusion
/**************************************************************************************************************************************************/

void loadSourceTermsZ(const PRECISION * const __restrict__ K, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, int s, PRECISION t,
PRECISION d_dz
) {
    for (unsigned int n = 0; n < NUMBER_ALL_EVOLVING_VARIABLES; ++n) {
        S[n] = 0.0;
    }
    
#ifndef IDEAL
	//=========================================================
	// spatial derivatives of the conserved variables \pi^{\mu\nu}
	//=========================================================
	PRECISION facZ = 1/d_dz/2;
	int ptr = 20; // 5 * n (with n = 4 corresponding to pitt)
#ifdef PIMUNU
	PRECISION dnpitt = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=5;
	PRECISION dnpitx = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=5;
	PRECISION dnpity = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=5;
	PRECISION dnpitn = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=15;
	PRECISION dnpixn = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=10;
	PRECISION dnpiyn = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;	ptr+=5;
	PRECISION dnpinn = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ; ptr+=5;
#endif
#ifdef PI
    PRECISION dnPi = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ; ptr+=5;
#endif
    //=========================================================
    // spatial derivatives of the conserved variables \nb^{\mu}
    //=========================================================
#ifdef VMU
    ptr+=5;
    PRECISION dnnbt = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ; ptr+=15;
    PRECISION dnnbn = (*(K + ptr + 3) - *(K + ptr + 1)) *facZ;
#endif
    
    
	PRECISION ut = u->ut[s];
	PRECISION un = u->un[s];

	//=========================================================
	// set dn terms in the source terms
	//=========================================================
	PRECISION vn = un / ut;
#ifdef PIMUNU
#ifndef PI
	S[0] = dnpitt*vn - dnpitn;
	S[3] = dnpitn*vn - dnpinn;
#else
	S[0] = dnpitt*vn - dnpitn - vn*dnPi;
	S[3] = dnpitn*vn - dnpinn - dnPi/pow(t,2);
#endif
	S[1] = dnpitx*vn - dnpixn;
	S[2] = dnpity*vn - dnpiyn;
#else
    S[0] = 0.0;
    S[1] = 0.0;
    S[2] = 0.0;
    S[3] = 0.0;
#endif
    //=========================================================
    // set dn terms in the source terms of baryon current
    //=========================================================
#ifdef NBMU
#ifdef VMU
    S[NUMBER_CONSERVED_VARIABLES] = dnnbt*vn - dnnbn;
#else
    S[NUMBER_CONSERVED_VARIABLES] = 0.0;
#endif
#endif
#endif
}


/**************************************************************************************************************************************************/
/* load source terms for all components, but exclude source terms directly dependent of gradients of shear, bulk and baryon diffusion
/**************************************************************************************************************************************************/

void loadSourceTerms2(const PRECISION * const __restrict__ Q, PRECISION * const __restrict__ S, const FLUID_VELOCITY * const __restrict__ u, PRECISION utp, PRECISION uxp, PRECISION uyp, PRECISION unp, PRECISION t, const PRECISION * const __restrict__ evec, const PRECISION * const __restrict__ pvec, int s, int d_ncx, int d_ncy, int d_ncz, PRECISION d_dt, PRECISION d_dx, PRECISION d_dy, PRECISION d_dz, const DYNAMICAL_SOURCES * const __restrict__ Source, const PRECISION * const __restrict__ rhobvec, PRECISION rhobp, const PRECISION * const __restrict__ alphaBvec, PRECISION alphaBp, const PRECISION * const __restrict__ Tvec, PRECISION Tp, PRECISION seq, const SLOW_MODES * const __restrict__ eqPhiQ, void * hydroParams)
{
	//=========================================================
	// conserved variables
	//=========================================================
	PRECISION ttt = Q[0];
	PRECISION ttx = Q[1];
	PRECISION tty = Q[2];
	PRECISION ttn = Q[3];
#ifdef PIMUNU
	PRECISION pitt = Q[4];
	PRECISION pitx = Q[5];
	PRECISION pity = Q[6];
	PRECISION pitn = Q[7];
	PRECISION pixx = Q[8];
	PRECISION pixy = Q[9];
	PRECISION pixn = Q[10];
	PRECISION piyy = Q[11];
	PRECISION piyn = Q[12];
	PRECISION pinn = Q[13];
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
	PRECISION Pi = Q[NUMBER_CONSERVED_VARIABLES_NO_BULK];
#else
	PRECISION Pi = 0;
#endif
#ifdef NBMU
    PRECISION Nbt = Q[NUMBER_CONSERVED_VARIABLES];
#else
    PRECISION Nbt = 0;
#endif
#ifdef VMU
    PRECISION nbt = Q[NUMBER_CONSERVED_VARIABLES+1];
    PRECISION nbx = Q[NUMBER_CONSERVED_VARIABLES+2];
    PRECISION nby = Q[NUMBER_CONSERVED_VARIABLES+3];
    PRECISION nbn = Q[NUMBER_CONSERVED_VARIABLES+4];
#else
    PRECISION nbt = 0;
    PRECISION nbx = 0;
    PRECISION nby = 0;
    PRECISION nbn = 0;
#endif
#ifdef HydroPlus
    PRECISION PhiQ[NUMBER_SLOW_MODES];
    PRECISION equiPhiQ[NUMBER_SLOW_MODES];
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
    {
        PhiQ[n] = Q[ALL_NUMBER_CONSERVED_VARIABLES+n];
        equiPhiQ[n] = eqPhiQ->phiQ[n][s];
    }
#else
    PRECISION PhiQ[NUMBER_SLOW_MODES];
    PRECISION equiPhiQ[NUMBER_SLOW_MODES];
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n)
    {
        PhiQ[n] = 0;
        equiPhiQ[n] = 0;
    }
#endif

	//=========================================================
	// primary variables
	//=========================================================
	PRECISION *utvec = u->ut;
	PRECISION *uxvec = u->ux;
	PRECISION *uyvec = u->uy;
	PRECISION *unvec = u->un;

	PRECISION p  = pvec[s];
	PRECISION ut = utvec[s];
	PRECISION ux = uxvec[s];
	PRECISION uy = uyvec[s];
	PRECISION un = unvec[s];

	//=========================================================
	// spatial derivatives of primary variables
	//=========================================================
	PRECISION facX = 1/d_dx/2;
	PRECISION facY = 1/d_dy/2;
	PRECISION facZ = 1/d_dz/2;
    
    int stride = d_ncx * d_ncy;
    
#ifndef MINMOD_FOR_U_AND_P
    
    // dx of u^{\mu} components
    PRECISION dxut = (*(utvec + s + 1) - *(utvec + s - 1)) * facX;
    PRECISION dxux = (*(uxvec + s + 1) - *(uxvec + s - 1)) * facX;
    PRECISION dxuy = (*(uyvec + s + 1) - *(uyvec + s - 1)) * facX;
    PRECISION dxun = (*(unvec + s + 1) - *(unvec + s - 1)) * facX;
    
    // dy of u^{\mu} components
    PRECISION dyut = (*(utvec + s + d_ncx) - *(utvec + s - d_ncx)) * facY;
    PRECISION dyux = (*(uxvec + s + d_ncx) - *(uxvec + s - d_ncx)) * facY;
    PRECISION dyuy = (*(uyvec + s + d_ncx) - *(uyvec + s - d_ncx)) * facY;
    PRECISION dyun = (*(unvec + s + d_ncx) - *(unvec + s - d_ncx)) * facY;
    
    // dn of u^{\mu} components
    PRECISION dnut = (*(utvec + s + stride) - *(utvec + s - stride)) * facZ;
    PRECISION dnux = (*(uxvec + s + stride) - *(uxvec + s - stride)) * facZ;
    PRECISION dnuy = (*(uyvec + s + stride) - *(uyvec + s - stride)) * facZ;
    PRECISION dnun = (*(unvec + s + stride) - *(unvec + s - stride)) * facZ;
    
    // pressure
    PRECISION dxp = (*(pvec + s + 1) - *(pvec + s - 1)) * facX;
    PRECISION dyp = (*(pvec + s + d_ncx) - *(pvec + s - d_ncx)) * facY;
    PRECISION dnp = (*(pvec + s + stride) - *(pvec + s - stride)) * facZ;
    
#else
    // dx of u^{\mu} components
    PRECISION utxp1 = utvec[s+1];
    PRECISION utxm1 = utvec[s-1];
    PRECISION uxxp1 = uxvec[s+1];
    PRECISION uxxm1 = uxvec[s-1];
    PRECISION uyxp1 = uyvec[s+1];
    PRECISION uyxm1 = uyvec[s-1];
    PRECISION unxp1 = unvec[s+1];
    PRECISION unxm1 = unvec[s-1];
    
    PRECISION dxut = approximateDerivative(utxm1,ut,utxp1) * facX * 2;
    PRECISION dxux = approximateDerivative(uxxm1,ux,uxxp1) * facX * 2;
    PRECISION dxuy = approximateDerivative(uyxm1,uy,uyxp1) * facX * 2;
    PRECISION dxun = approximateDerivative(unxm1,un,unxp1) * facX * 2;
    
	// dy of u^{\mu} components
    PRECISION utyp1 = utvec[s+d_ncx];
    PRECISION utym1 = utvec[s-d_ncx];
    PRECISION uxyp1 = uxvec[s+d_ncx];
    PRECISION uxym1 = uxvec[s-d_ncx];
    PRECISION uyyp1 = uyvec[s+d_ncx];
    PRECISION uyym1 = uyvec[s-d_ncx];
    PRECISION unyp1 = unvec[s+d_ncx];
    PRECISION unym1 = unvec[s-d_ncx];
    
    PRECISION dyut = approximateDerivative(utym1,ut,utyp1) * facY * 2;
    PRECISION dyux = approximateDerivative(uxym1,ux,uxyp1) * facY * 2;
    PRECISION dyuy = approximateDerivative(uyym1,uy,uyyp1) * facY * 2;
    PRECISION dyun = approximateDerivative(unym1,un,unyp1) * facY * 2;
    
	// dn of u^{\mu} components
    PRECISION utnp1 = utvec[s+stride];
    PRECISION utnm1 = utvec[s-stride];
    PRECISION uxnp1 = uxvec[s+stride];
    PRECISION uxnm1 = uxvec[s-stride];
    PRECISION uynp1 = uyvec[s+stride];
    PRECISION uynm1 = uyvec[s-stride];
    PRECISION unnp1 = unvec[s+stride];
    PRECISION unnm1 = unvec[s-stride];
    
    PRECISION dnut = approximateDerivative(utnm1,ut,utnp1) * facZ * 2;
    PRECISION dnux = approximateDerivative(uxnm1,ux,uxnp1) * facZ * 2;
    PRECISION dnuy = approximateDerivative(uynm1,uy,uynp1) * facZ * 2;
    PRECISION dnun = approximateDerivative(unnm1,un,unnp1) * facZ * 2;
    
	// pressure
    PRECISION pc = pvec[s];
    PRECISION pxp1 = pvec[s+1];
    PRECISION pxm1 = pvec[s-1];
    PRECISION pyp1 = pvec[s+d_ncx];
    PRECISION pym1 = pvec[s-d_ncx];
    PRECISION pnp1 = pvec[s+stride];
    PRECISION pnm1 = pvec[s-stride];
    
    PRECISION dxp = approximateDerivative(pxm1,pc,pxp1) * facX * 2;
    PRECISION dyp = approximateDerivative(pym1,pc,pyp1) * facY * 2;
    PRECISION dnp = approximateDerivative(pnm1,pc,pnp1) * facZ * 2;
#endif
    
    //=========================================================
    //Calculate the gradient of chemical potential to temperature ratio
    //=========================================================
#ifdef VMU
    PRECISION es = evec[s];
    PRECISION Ts  = Tvec[s];
    PRECISION rhobs = rhobvec[s];
    PRECISION alphaBs  = alphaBvec[s];
    
#ifndef MINMOD_FOR_BARYON
    // derivatives of rhob
    PRECISION dtrhob = (rhobs - rhobp) / d_dt;
    PRECISION dxrhob = (*(rhobvec + s + 1) - *(rhobvec + s - 1)) * facX;
    PRECISION dyrhob = (*(rhobvec + s + d_ncx) - *(rhobvec + s - d_ncx)) * facY;
    PRECISION dnrhob = (*(rhobvec + s + stride) - *(rhobvec + s - stride)) * facZ;
    
    // derivatives of T
    PRECISION dtT = (Ts - Tp) / d_dt;
    PRECISION dxT = (*(Tvec + s + 1) - *(Tvec + s - 1)) * facX;
    PRECISION dyT = (*(Tvec + s + d_ncx) - *(Tvec + s - d_ncx)) * facY;
    PRECISION dnT = (*(Tvec + s + stride) - *(Tvec + s - stride)) * facZ;
    
    // derivatives of muB/T
    PRECISION dtalphaB = (alphaBs - alphaBp) / d_dt;
    PRECISION dxalphaB = (*(alphaBvec + s + 1) - *(alphaBvec + s - 1)) * facX;
    PRECISION dyalphaB = (*(alphaBvec + s + d_ncx) - *(alphaBvec + s - d_ncx)) * facY;
    PRECISION dnalphaB = (*(alphaBvec + s + stride) - *(alphaBvec + s - stride)) * facZ;

#else
    // derivatives of rhob
    PRECISION rhobxp1 = rhobvec[s+1];
    PRECISION rhobxp2 = rhobvec[s+2];
    PRECISION rhobxm1 = rhobvec[s-1];
    PRECISION rhobxm2 = rhobvec[s-2];
    PRECISION rhobyp1 = rhobvec[s+d_ncx];
    PRECISION rhobyp2 = rhobvec[s+2*d_ncx];
    PRECISION rhobym1 = rhobvec[s-d_ncx];
    PRECISION rhobym2 = rhobvec[s-2*d_ncx];
    PRECISION rhobnp1 = rhobvec[s+stride];
    PRECISION rhobnp2 = rhobvec[s+2*stride];
    PRECISION rhobnm1 = rhobvec[s-stride];
    PRECISION rhobnm2 = rhobvec[s-2*stride];
    
    PRECISION dtrhob = (rhobs - rhobp) / d_dt;
    PRECISION dxrhob = approximateDerivative(rhobxm1,rhobs,rhobxp1) * facX * 2;
    PRECISION dyrhob = approximateDerivative(rhobym1,rhobs,rhobyp1) * facY * 2;
    PRECISION dnrhob = approximateDerivative(rhobnm1,rhobs,rhobnp1) * facZ * 2;
    
    // derivatives of T
    PRECISION Txp1 = Tvec[s+1];
    PRECISION Txp2 = Tvec[s+2];
    PRECISION Txm1 = Tvec[s-1];
    PRECISION Txm2 = Tvec[s-2];
    PRECISION Typ1 = Tvec[s+d_ncx];
    PRECISION Typ2 = Tvec[s+2*d_ncx];
    PRECISION Tym1 = Tvec[s-d_ncx];
    PRECISION Tym2 = Tvec[s-2*d_ncx];
    PRECISION Tnp1 = Tvec[s+stride];
    PRECISION Tnp2 = Tvec[s+2*stride];
    PRECISION Tnm1 = Tvec[s-stride];
    PRECISION Tnm2 = Tvec[s-2*stride];
    
    PRECISION dtT = (Ts - Tp) / d_dt;
    PRECISION dxT = approximateDerivative(Txm1,Ts,Txp1) * facX * 2;
    PRECISION dyT = approximateDerivative(Tym1,Ts,Typ1) * facY * 2;
    PRECISION dnT = approximateDerivative(Tnm1,Ts,Tnp1) * facZ * 2;
    
    // derivatives of muB/T
    PRECISION alphabxp1 = alphaBvec[s+1];
    PRECISION alphabxp2 = alphaBvec[s+2];
    PRECISION alphabxm1 = alphaBvec[s-1];
    PRECISION alphabxm2 = alphaBvec[s-2];
    PRECISION alphabyp1 = alphaBvec[s+d_ncx];
    PRECISION alphabyp2 = alphaBvec[s+2*d_ncx];
    PRECISION alphabym1 = alphaBvec[s-d_ncx];
    PRECISION alphabym2 = alphaBvec[s-2*d_ncx];
    PRECISION alphabnp1 = alphaBvec[s+stride];
    PRECISION alphabnp2 = alphaBvec[s+2*stride];
    PRECISION alphabnm1 = alphaBvec[s-stride];
    PRECISION alphabnm2 = alphaBvec[s-2*stride];
    
    PRECISION dtalphaB = (alphaBs - alphaBp) / d_dt;
    PRECISION dxalphaB = approximateDerivative(alphabxm1,alphaBs,alphabxp1) * facX * 2;
    PRECISION dyalphaB = approximateDerivative(alphabym1,alphaBs,alphabyp1) * facY * 2;
    PRECISION dnalphaB = approximateDerivative(alphabnm1,alphaBs,alphabnp1) * facZ * 2;
#endif
    
    // gradient of rhob
    PRECISION ukdk_rhob = ut * dtrhob + ux * dxrhob + uy * dyrhob + un * dnrhob;
    PRECISION Nablat_rhob =  dtrhob - ut * ukdk_rhob;
    PRECISION Nablax_rhob = -dxrhob - ux * ukdk_rhob;
    PRECISION Nablay_rhob = -dyrhob - uy * ukdk_rhob;
    PRECISION Nablan_rhob = -1/pow(t,2)*dnrhob - un * ukdk_rhob;
    
    // gradient of T
    PRECISION ukdk_T = ut * dtT + ux * dxT + uy * dyT + un * dnT;
    PRECISION Nablat_T =  dtT - ut * ukdk_T;
    PRECISION Nablax_T = -dxT - ux * ukdk_T;
    PRECISION Nablay_T = -dyT - uy * ukdk_T;
    PRECISION Nablan_T = -1/pow(t,2)*dnT - un * ukdk_T;
    
    // gradient of muB/T
    PRECISION ukdk_alphaB = ut * dtalphaB + ux * dxalphaB + uy * dyalphaB + un * dnalphaB;
    PRECISION Nablat_alphaB =  dtalphaB - ut * ukdk_alphaB;
    PRECISION Nablax_alphaB = -dxalphaB - ux * ukdk_alphaB;
    PRECISION Nablay_alphaB = -dyalphaB - uy * ukdk_alphaB;
    PRECISION Nablan_alphaB = -1/pow(t,2)*dnalphaB - un * ukdk_alphaB;
    
#else
    PRECISION es = evec[s];
    PRECISION Ts  = Tvec[s];
    PRECISION rhobs = rhobvec[s];
    PRECISION alphaBs  = alphaBvec[s];
    
    PRECISION Nablat_rhob = 0;
    PRECISION Nablax_rhob = 0;
    PRECISION Nablay_rhob = 0;
    PRECISION Nablan_rhob = 0;
    
    PRECISION Nablat_T = 0;
    PRECISION Nablax_T = 0;
    PRECISION Nablay_T = 0;
    PRECISION Nablan_T = 0;
    
    PRECISION Nablat_alphaB = 0;
    PRECISION Nablax_alphaB = 0;
    PRECISION Nablay_alphaB = 0;
    PRECISION Nablan_alphaB = 0;
#endif

	//=========================================================
	// T^{\mu\nu} source terms
	//=========================================================
	PRECISION tnn = Tnn(evec[s],p+Pi,un,pinn,t);
	PRECISION vx = ux/ut;
	PRECISION vy = uy/ut;
	PRECISION vn = un/ut;
	PRECISION dxvx = (dxux - vx * dxut)/ ut;
	PRECISION dyvy = (dyuy - vy * dyut)/ ut;
	PRECISION dnvn = (dnun - vn * dnut)/ ut;
    PRECISION dkvk = dxvx + dyvy + dnvn;
    
    // added dynamical source terms
#ifdef DYNAMICAL_SOURCE
	S[0] = Source->sourcet[s] - (ttt / t + t * tnn) + dkvk*(pitt-p-Pi) - vx*dxp - vy*dyp - vn*dnp;
	S[1] = Source->sourcex[s] - ttx/t -dxp + dkvk*pitx;
	S[2] = Source->sourcey[s] - tty/t -dyp + dkvk*pity;
	S[3] = Source->sourcen[s] - 3*ttn/t -dnp/pow(t,2) + dkvk*pitn;
#else
    S[0] = - (ttt / t + t * tnn) + dkvk*(pitt-p-Pi) - vx*dxp - vy*dyp - vn*dnp;
    S[1] = - ttx/t -dxp + dkvk*pitx;
    S[2] = - tty/t -dyp + dkvk*pity;
    S[3] = - 3*ttn/t -dnp/pow(t,2) + dkvk*pitn;
#endif
#ifdef USE_CARTESIAN_COORDINATES
	S[0] = + dkvk*(pitt-p-Pi) - vx*dxp - vy*dyp - vn*dnp;
	S[1] = - dxp + dkvk*pitx;
	S[2] = - dyp + dkvk*pity;
	S[3] = - dnp + dkvk*pitn;
#endif

    //=========================================================
    // N^{\mu} source terms
    //=========================================================
#ifdef NBMU
#ifdef DYNAMICAL_SOURCE
    S[NUMBER_CONSERVED_VARIABLES] = Source->sourceb[s] - Nbt/t + dkvk*nbt;
#else
    S[NUMBER_CONSERVED_VARIABLES] = - Nbt/t + dkvk*nbt;
#endif
#ifdef USE_CARTESIAN_COORDINATES
    S[NUMBER_CONSERVED_VARIABLES] = + dkvk*nbt;
#endif
#endif

	//=========================================================
	// Dissipative components source terms
	//=========================================================

	PRECISION pimunuRHS[NUMBER_PROPAGATED_PIMUNU_COMPONENTS];
    PRECISION piRHS;
    PRECISION nbmuRHS[NUMBER_PROPAGATED_VMU_COMPONENTS];
    PRECISION phiQRHS[NUMBER_SLOW_MODES];
    
	setDissipativeSourceTerms(pimunuRHS, &piRHS, nbmuRHS, phiQRHS, nbt, nbx, nby, nbn, rhobs, Nablat_rhob, Nablax_rhob, Nablay_rhob, Nablan_rhob, alphaBs, Nablat_alphaB, Nablax_alphaB, Nablay_alphaB, Nablan_alphaB, Ts, Nablat_T, Nablax_T, Nablay_T, Nablan_T, seq, t, es, p, ut, ux, uy, un, utp, uxp, uyp, unp, pitt, pitx, pity, pitn, pixx, pixy, pixn, piyy, piyn, pinn, Pi, dxut, dyut, dnut, dxux, dyux, dnux, dxuy, dyuy, dnuy, dxun, dyun, dnun, dkvk, d_dt, PhiQ, equiPhiQ, hydroParams);
    
#ifdef PIMUNU
    for(unsigned int n = 0; n < NUMBER_PROPAGATED_PIMUNU_COMPONENTS; ++n) S[n+4] = pimunuRHS[n];// for shear
#endif
#ifdef PI
    S[NUMBER_CONSERVED_VARIABLES_NO_BULK] = piRHS;// for bulk
#endif
#ifdef VMU
    for(unsigned int n = 0; n < NUMBER_PROPAGATED_VMU_COMPONENTS; ++n) S[n+1+NUMBER_CONSERVED_VARIABLES] = nbmuRHS[n];// for baryon diffusion current
#endif
#ifdef HydroPlus
    for(unsigned int n = 0; n < NUMBER_SLOW_MODES; ++n) S[ALL_NUMBER_CONSERVED_VARIABLES+n] = phiQRHS[n];// for slow modes
#endif
}
