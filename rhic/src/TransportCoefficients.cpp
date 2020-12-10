//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <stdlib.h>
#include <math.h>

#include "../include/EquationOfState.h"
#include "../include/TransportCoefficients.h"
#include "../include/HydroPlus.h"

#define HBARC 0.197326938


/**************************************************************************************************************************************************/
/* bulk viscosity as a function of temperature
/**************************************************************************************************************************************************/

// G. Denicol et. al. PRC80 (2009) 064901. parameters for the analytic parameterization of the bulk viscosity \zeta/S
#define A_1 -13.77
#define A_2 27.55
#define A_3 13.45

#define LAMBDA_1 0.9
#define LAMBDA_2 0.25
#define LAMBDA_3 0.9
#define LAMBDA_4 0.22

#define SIGMA_1 0.025
#define SIGMA_2 0.13
#define SIGMA_3 0.0025
#define SIGMA_4 0.022

PRECISION bulkViscosityToEntropyDensity(PRECISION T) {
    PRECISION x = T/1.01355;
    if(x > 1.05)
        return LAMBDA_1*exp(-(x-1)/SIGMA_1) + LAMBDA_2*exp(-(x-1)/SIGMA_2)+0.001;
    else if(x < 0.995)
        return LAMBDA_3*exp((x-1)/SIGMA_3)+ LAMBDA_4*exp((x-1)/SIGMA_4)+0.03;
    else
        return A_1*x*x + A_2*x - A_3;
}

/**************************************************************************************************************************************************/
/* baryon diffusion coefficient of the medium from kinetic theory
/**************************************************************************************************************************************************/

// G. S. Denicol, C. Gale, S. Jeon, A. Monnai, B. Schenke, and C. Shen, (2018), arXiv:1804.10557 [nucl-th].
PRECISION baryonDiffusionCoefficientKinetic(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p){
    PRECISION HyCotangent = 1/tanh(alphaB);
    if(isnan(HyCotangent)) printf("kappaB is nan. e=%4e,\t rhob=%4e,\t mub_over_T=%4e,\t T=%4e. \n",e,rhob, alphaB, T);
    
    return Cb/T * rhob * (0.3333333 * HyCotangent - rhob*T/(e+p));
}

// compare to MUSIC
PRECISION baryonDiffusionCoefficientTest(PRECISION T, PRECISION rhob, PRECISION alphaB){
    return 0.2 * rhob / (alphaB * T);
}

// D. T. Son and A. O. Starinets, JHEP 03, 052 (2006).
PRECISION baryonDiffusionCoefficientAdscft(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p, PRECISION seq){
    PRECISION fac = rhob*T/(e+p);
    PRECISION mub = (alphaB * T);
    
    return fac * fac * T * seq * 2.0 * M_PI / (mub * mub);
}

// M. Stephanov and Y. Yin, Phys. Rev., D98(3):036006, 2018.
PRECISION baryonDiffusionCoefficientHydroPlus(PRECISION T, PRECISION rhob, PRECISION alphaB, PRECISION e, PRECISION p, PRECISION seq){
    PRECISION fac = rhob*T/(e+p);
    
    return fac * fac * T * seq / rhob / (1.2 * M_PI);
}


/**************************************************************************************************************************************************/
/* Baryon diffusion coefficients
/**************************************************************************************************************************************************/

// Baryon diffusion table from PRL 115 (2015) 202301
void getBaryonDiffusionCoefficientTable(){
#ifdef VMU
    FILE *filebarycoeff;
    PRECISION temp;
    PRECISION chem;
    PRECISION chiB;
    
    // in the table: T [MeV] , muB [MeV] , chi2B / T^2 , sigmaB / T , T * DB
    filebarycoeff = fopen ("input/coefficients/BaryonDiffData.dat","r");
    if(filebarycoeff==NULL){
        printf("The file BaryonDiffData.dat was not opened...\n");
        exit(-1);
    }
    else
    {
        fseek(filebarycoeff,0L,SEEK_SET);
        for(int i = 0; i < 128721; ++i){
            fscanf(filebarycoeff,"%lf %lf %lf %lf %lf", & temp, & chem, & chiB, & BaryDiffCoeff->sigmaB[i], & BaryDiffCoeff->DB[i]);
        }
    }
    fclose(filebarycoeff);
    
    printf("Baryon diffusion coefficients table is read in.\n");
#endif
}

void baryonDiffusionCoefficient(PRECISION T, PRECISION muB, PRECISION * const __restrict__ diffusionCoeff){
    
    PRECISION T0 = T*HBARC*1000; // MeV
    PRECISION muB0 = fabs(muB*HBARC*1000); // MeV
    
    PRECISION sigmaBT, DBT; // sigmaB/T and D_B*T
    
    diffusionCoeff[0] = 0.0;
    diffusionCoeff[1] = 0.0;
    
#ifdef VMU
    if((100<=T0)&&(T0<419)){
        if((0<=muB0)&&(muB0<400)){
            sigmaBT = InferredPrimaryVariable(T0, muB0, 100., 1., 401, 1., 0, BaryDiffCoeff->sigmaB);
            DBT = InferredPrimaryVariable(T0, muB0, 100., 1., 401, 1., 0, BaryDiffCoeff->DB);
        }else{
            sigmaBT = InferredPrimaryVariable(T0, 400., 100., 1., 401, 1., 0, BaryDiffCoeff->sigmaB);
            DBT = InferredPrimaryVariable(T0, 400., 100., 1., 401, 1., 0, BaryDiffCoeff->DB);
        }
    }else if(T0<100)
    {
        if((0<=muB0)&&(muB0<400)){
            sigmaBT = InferredPrimaryVariable(100., muB0, 100., 1., 401, 1., 0, BaryDiffCoeff->sigmaB);
            DBT = InferredPrimaryVariable(100., muB0, 100., 1., 401, 1., 0, BaryDiffCoeff->DB);
        }else{ // using values at T=100 and muB=400
            sigmaBT = 0.0005433609062953698;
            DBT = 0.034289064162551376;
        }
    }else
    {
        if((0<=muB0)&&(muB0<400)){
            sigmaBT = InferredPrimaryVariable(419., muB0, 100., 1., 401, 1., 0, BaryDiffCoeff->sigmaB);
            DBT = InferredPrimaryVariable(419., muB0, 100., 1., 401, 1., 0, BaryDiffCoeff->DB);
        }else{ // using values at T=420 and muB=400
            sigmaBT = 0.049794598243803105;
            DBT = 0.1547150351237801;
        }
    }
    
    diffusionCoeff[0] = sigmaBT * T * T; // kappaB in [fm^-2]
    diffusionCoeff[1] = DBT / T; // DB in [fm]
#endif
}

/**************************************************************************************************************************************************/
/* jet-medium paramemters
/**************************************************************************************************************************************************/

// conformal A. Ficnar et al 1311.6160, fac taken from 1507.06556
PRECISION eHatCFT(PRECISION T, PRECISION x){
    
    PRECISION fac = 4.387;//pi/2*sqrt(lambda_t)
    
    return fac * T * T * (1 + M_PI * T * x) * (1 + M_PI * T * x);
}

// conformal arXiv:hep-ph/0605178, fac taken from 1507.06556
PRECISION qHatCFT(PRECISION T){
    
    PRECISION fac = 21.02;//4.6;//
    
    return fac * T * T * T;
}

// tables from 1507.06556
void getqhatTable(){
#ifdef JET
    FILE *fileqhat;
    PRECISION temp;
    PRECISION chem;
    
    // in the table: T [MeV] , muB [MeV] , qhat/qhat_CFT
    fileqhat = fopen ("input/coefficients/qhat.dat","r");
    if(fileqhat==NULL){
        printf("The file qhat.dat was not opened...\n");
        exit(-1);
    }
    else
    {
        fseek(fileqhat,0L,SEEK_SET);
        for(int i = 0; i < 5525; ++i){
            fscanf(fileqhat,"%lf %lf %lf", & temp, & chem, & JetCoeff->qhat[i]);
        }
    }
    fclose(fileqhat);
    
    printf("qhat table is read in.\n");
#endif
}

PRECISION qHat(PRECISION T, PRECISION muB){
    
    PRECISION T0 = T*HBARC*1000; // MeV
    PRECISION muB0 = fabs(muB*HBARC*1000); // MeV
    
    PRECISION qHAT; // qhat/qhat_CFT
    
    if((100<=T0)&&(T0<419)){
        if((0<=muB0)&&(muB0<419)){
            qHAT = InferredPrimaryVariable(T0, muB0, 100., 5., 85, 5., 0, JetCoeff->qhat);
            //DBT = InferredPrimaryVariable(T0, muB0, 100., 1., 401, 1., 0, BaryDiffCoeff->DB);
        }else{
            qHAT = InferredPrimaryVariable(T0, 419, 100., 5., 85, 5., 0, JetCoeff->qhat);
            //DBT = InferredPrimaryVariable(T0, 400., 100., 1., 401, 1., 0, BaryDiffCoeff->DB);
        }
    }else if(T0<100)
    {
        if((0<=muB0)&&(muB0<419)){
            qHAT = InferredPrimaryVariable(100., muB0, 100., 5., 85, 5., 0, JetCoeff->qhat);
            //DBT = InferredPrimaryVariable(100., muB0, 100., 1., 401, 1., 0, BaryDiffCoeff->DB);
        }else{ // using values at T=100 and muB=400
            qHAT = 1.0214966454406247;
            //DBT = 0.034289064162551376;
        }
    }else
    {
        if((0<=muB0)&&(muB0<419)){
            qHAT = InferredPrimaryVariable(419., muB0, 100., 5., 85, 5., 0, JetCoeff->qhat);
            //DBT = InferredPrimaryVariable(419., muB0, 100., 1., 401, 1., 0, BaryDiffCoeff->DB);
        }else{ // using values at T=420 and muB=400
            qHAT = 1.3084476621169798;
            //DBT = 0.1547150351237801;
        }
    }
    
    // qHAT is the ratio between qhat_nonconformal/qhat_conformal
    
    return qHAT * qHatCFT(T); // in fm^-3 since T is in fm^-1 and so is qHatCFT
}

