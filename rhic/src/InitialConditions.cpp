//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <math.h> // for math functions
#include <stdio.h> // for printf
#include <stdlib.h> //TEMP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>

#include "../include/InitialConditions.h"
#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/GlauberModel.h"
#include "../include/MonteCarloGlauberModel.h"
#include "../include/HydroParameters.h"
#include "../include/EquationOfState.h"

#include "../include/HydroPlus.h"
#include "../include/SourceTerms.h"

#define THETA_FUNCTION(X) ((double)X < (double)0 ? (double)0 : (double)1)
#define REGULATE 1 // 1 to regulate flow in dilute regions
#define GAMMA_MAX 10.0
#define HBARC 0.197326938
//#define Avg_MCGlauber

using namespace std;

/**************************************************************************************************************************************************/
/* Read in all initial profiles from a single file
/**************************************************************************************************************************************************/

//this reads all hydro variables from a single file; this way we do not need to fetch the coordinates many times
//note that the file must contain values for all dissipative currents, even if they are zero !!!
void setInitialTmunuFromFile(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    float x, y, z, e_in, p_in, ut_in, ux_in, uy_in, un_in;
    float pitt_in, pitx_in, pity_in, pitn_in, pixx_in, pixy_in, pixn_in, piyy_in, piyn_in, pinn_in;
    float Pi_in;

    FILE *fileIn;
    char fname[255];
    
    sprintf(fname, "%s/%s", rootDirectory, "/input/Tmunu.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open Tmunu.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &x, &y, &z, &e_in, &p_in, &ut_in, &ux_in, &uy_in, &un_in, &pitt_in, &pitx_in, &pity_in, &pitn_in, &pixx_in, &pixy_in, &pixn_in, &piyy_in, &piyn_in, &pinn_in, &Pi_in);
                    
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    
                    e[s] =  (PRECISION) e_in;
                    p[s] = p_in;
                    
                    u->ut[s] = ut_in;
                    u->ux[s] = ux_in;
                    u->uy[s] = uy_in;
                    u->un[s] = un_in;

#ifdef PIMUNU
                    q->pitt[s] = pitt_in;
                    q->pitx[s] = pitx_in;
                    q->pity[s] = pity_in;
                    q->pitn[s] = pitn_in;
                    q->pixx[s] = pixx_in;
                    q->pixy[s] = pixy_in;
                    q->pixn[s] = pixn_in;
                    q->piyy[s] = piyy_in;
                    q->piyn[s] = piyn_in;
                    q->pinn[s] = pinn_in;
#endif
#ifdef PI
                    q->Pi[s] = Pi_in;
#endif
                }
            }
        }
    }
    fclose(fileIn);
}


/**************************************************************************************************************************************************/
/* Read in all initial profiles from seperate files
/**************************************************************************************************************************************************/

//this function reads a separate file for every hydrodynamic variable
void setInitialTmunuFromFiles(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    float x, y, z, value;
    FILE *fileIn;
    char fname[255];
    
    //energy density
    sprintf(fname, "%s/%s", rootDirectory, "/input/e.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open e.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    e[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pressure
    sprintf(fname, "%s/%s", rootDirectory, "/input/p.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open p.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    p[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //ut
    sprintf(fname, "%s/%s", rootDirectory, "/input/ut.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open ut.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->ut[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //ux
    sprintf(fname, "%s/%s", rootDirectory, "/input/ux.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open ux.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->ux[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //uy
    sprintf(fname, "%s/%s", rootDirectory, "/input/uy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open uy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->uy[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //un
    sprintf(fname, "%s/%s", rootDirectory, "/input/un.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open un.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    u->un[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
#ifdef PIMUNU
    //pitt
    sprintf(fname, "%s/%s", rootDirectory, "/input/pitt.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pitt.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pitt[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pitx
    sprintf(fname, "%s/%s", rootDirectory, "/input/pitx.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pitx.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pitx[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pity
    sprintf(fname, "%s/%s", rootDirectory, "/input/pity.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pity.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pity[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pitn
    sprintf(fname, "%s/%s", rootDirectory, "/input/pitn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pitn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pitn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pixx
    sprintf(fname, "%s/%s", rootDirectory, "/input/pixx.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pixx.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pixx[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pixy
    sprintf(fname, "%s/%s", rootDirectory, "/input/pixy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pixy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pixy[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pixn
    sprintf(fname, "%s/%s", rootDirectory, "/input/pixn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pixn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pixn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //piyy
    sprintf(fname, "%s/%s", rootDirectory, "/input/piyy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open piyy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->piyy[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //piyn
    sprintf(fname, "%s/%s", rootDirectory, "/input/piyn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open piyn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->piyn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    
    //pinn
    sprintf(fname, "%s/%s", rootDirectory, "/input/pinn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pinn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->pinn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
#endif
#ifdef PI
    //bulk
    sprintf(fname, "%s/%s", rootDirectory, "/input/bulk.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open bulk.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    q->Pi[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
#endif
}

/**************************************************************************************************************************************************/
/* Set initial flow profile
/*	- u^\mu = (1, 0, 0, 0)
/* 	- No transverse flow (ux = uy = 0)
/*	- Longitudinal scaling flow (u_z = z/t, i.e. un = 0)
/**************************************************************************************************************************************************/
void setFluidVelocityInitialCondition(void * latticeParams, void * hydroParams) {
    
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double t0 = hydro->initialProperTimePoint;

	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				PRECISION ux = 0;
				PRECISION uy = 0;
				PRECISION un = 0;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = sqrt(1+ux*ux+uy*uy+t0*t0*un*un);
			}
		}
	}
}

/**************************************************************************************************************************************************/
/* Set initial baryon diffusion current to be 0
/**************************************************************************************************************************************************/
void setNbmuInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
    
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
#ifdef VMU
    printf("Initialize \\nb^\\mu to be zero.\n");
    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                q->nbt[s] = 0.0;
                q->nbx[s] = 0.0;
                q->nby[s] = 0.0;
                q->nbn[s] = 0.0;
            }
        }
    }
#endif
}

/**************************************************************************************************************************************************/
/* Set initial shear-stress tensor \pi^\mu\nu
/*	- Navier-Stokes value, i.e. \pi^\mu\nu = 2 * (\epsilon + P) / T * \eta/S * \sigma^\mu\nu
/* 	- No initial pressure anisotropies (\pi^\mu\nu = 0)
/**************************************************************************************************************************************************/
void setPimunuNavierStokesInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
    
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	PRECISION dx = (PRECISION)(lattice->latticeSpacingX);
	PRECISION dz = (PRECISION)(lattice->latticeSpacingRapidity);

	PRECISION etabar = (PRECISION)(hydro->shearViscosityToEntropyDensity);
	PRECISION t = hydro->initialProperTimePoint;

	PRECISION e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				PRECISION T = effectiveTemperatureWB(e[s]);
				if (T == 0) T = 1.e-3;
				PRECISION pinn = -4.0/(3*t*t*t)*etabar*(e[s]+p[s])/T;
#ifdef PIMUNU
				q->pitt[s] = 0;
				q->pitx[s] = 0;
				q->pity[s] = 0;
				q->pitn[s] = 0;
				q->pixx[s] = -t*t*pinn/2;
				q->pixy[s] = 0;
				q->pixn[s] = 0;
				q->piyy[s] = -t*t*pinn/2;
				q->piyn[s] = 0;
				q->pinn[s] = pinn;
#endif
			}
		}
	}
}

void setPiNavierStokesInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    PRECISION dx = (PRECISION)(lattice->latticeSpacingX);
    PRECISION dz = (PRECISION)(lattice->latticeSpacingRapidity);
    
    PRECISION etabar = (PRECISION)(hydro->shearViscosityToEntropyDensity);
    PRECISION t = hydro->initialProperTimePoint;
    
    PRECISION e0 = initCond->initialEnergyDensity;
    
    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                PRECISION T = effectiveTemperatureWB(e[s]);
                if (T == 0) T = 1.e-3;

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
                    PRECISION x = T/1.01355;
                    PRECISION zetabar = A_1*x*x + A_2*x - A_3;
                    if(x > 1.05)
                        zetabar = LAMBDA_1*exp(-(x-1)/SIGMA_1) + LAMBDA_2*exp(-(x-1)/SIGMA_2)+0.001;
                    else if(x < 0.995)
                        zetabar = LAMBDA_3*exp((x-1)/SIGMA_3)+ LAMBDA_4*exp((x-1)/SIGMA_4)+0.03;
#ifdef PI
                    q->Pi[s] = -zetabar*(e[s]+p[s])/T/t;
#endif
            }
        }
    }
}

/**************************************************************************************************************************************************/
/* set initial shear stress tensor and bulk pressure to be 0
/**************************************************************************************************************************************************/
void setPimunuInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
	int initializePimunuNavierStokes = hydro->initializePimunuNavierStokes;
    int initializePiNavierStokes = hydro->initializePiNavierStokes;

#ifdef PIMUNU
	if (initializePimunuNavierStokes==1) {
		printf("Initialize \\pi^\\mu\\nu to its asymptotic Navier-Stokes value.\n");
		setPimunuNavierStokesInitialCondition(latticeParams, initCondParams, hydroParams);
	}
	else {
		printf("Initialize \\pi^\\mu\\nu to zero.\n");
		struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
		int nx = lattice->numLatticePointsX;
		int ny = lattice->numLatticePointsY;
		int nz = lattice->numLatticePointsRapidity;
		for(int i = 2; i < nx+2; ++i) {
			for(int j = 2; j < ny+2; ++j) {
				for(int k = 2; k < nz+2; ++k) {
					int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
			  		q->pitt[s] = 0;
			  		q->pitx[s] = 0;
			  		q->pity[s] = 0;
			  		q->pitn[s] = 0;
			  		q->pixx[s] = 0;
			  		q->pixy[s] = 0;
			  		q->pixn[s] = 0;
			  		q->piyy[s] = 0;
			  		q->piyn[s] = 0;
			  		q->pinn[s] = 0;
				}
			}
		}
	}
#endif
    
#ifdef PI
    if (initializePiNavierStokes==1) {
        printf("Initialize \\Pi to its asymptotic Navier-Stokes value.\n");
        setPiNavierStokesInitialCondition(latticeParams, initCondParams, hydroParams);
        return;
    }
    else {
        printf("Initialize \\Pi to zero.\n");
        struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
        int nx = lattice->numLatticePointsX;
        int ny = lattice->numLatticePointsY;
        int nz = lattice->numLatticePointsRapidity;
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    
                    q->Pi[s] = 0;
                }
            }
        }
        return;
    }
#endif
}

/**************************************************************************************************************************************************/
/* Longitudinal initial energy density distribution
/**************************************************************************************************************************************************/
void longitudinalEnergyDensityDistribution(double * const __restrict__ eL, void * latticeParams, void * initCondParams) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    
    int nz = lattice->numLatticePointsRapidity;
    double dz = lattice->latticeSpacingRapidity;
    
    double etaFlat = initCond->rapidityMean;
    double etaVariance = initCond->rapidityVariance;
    
    for(int k = 0; k < nz; ++k) {
        double eta = (k - (nz-1)/2)*dz;
        double etaScaled = fabs(eta) - etaFlat/2;
        double arg = -etaScaled * etaScaled / (etaVariance * etaVariance) / 2 * THETA_FUNCTION(etaScaled);
        eL[k] = 2 * exp(arg);
    }
}

/**************************************************************************************************************************************************/
/* Longitudinal initial baryon density distribution
/**************************************************************************************************************************************************/
void longitudinalBaryonDensityDistribution(double * const __restrict__ rhoLa, double * const __restrict__ rhoLb, void * latticeParams, void * initCondParams) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    
    int nz = lattice->numLatticePointsRapidity;
    double dz = lattice->latticeSpacingRapidity;
    
    double etaVariance1 = initCond->bRapidityVariance1;
    double etaVariance2 = initCond->bRapidityVariance2;
    double etaMean = initCond->bRapidityMean;

    for(int k = 0; k < nz; ++k) {
        double eta = (k - (nz-1)/2)*dz;
        rhoLa[k] = (exp(-(eta-etaMean)*(eta-etaMean)/(2*etaVariance1*etaVariance1))*THETA_FUNCTION(eta-etaMean+1.e-2) + exp(-(eta-etaMean)*(eta-etaMean)/(2*etaVariance2*etaVariance2))*THETA_FUNCTION(etaMean-eta-1.e-2));
        rhoLb[k] = (exp(-(-eta-etaMean)*(-eta-etaMean)/(2*etaVariance1*etaVariance1))*THETA_FUNCTION(-eta-etaMean+1.e-2) + exp(-(-eta-etaMean)*(-eta-etaMean)/(2*etaVariance2*etaVariance2))*THETA_FUNCTION(etaMean+eta-1.e-2));
    }
}


/**************************************************************************************************************************************************/
/* Initial conditions for hydro with dynamical sources
/**************************************************************************************************************************************************/
void setICfromSource(void * latticeParams, void * initCondParams, void * hydroParams, const char * rootDirectory){
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    
    int sourceType = initCond->sourceType;
    
    double initialEnergyDensity = initCond->initialEnergyDensity;
    double initialBaryonDensity = initCond->initialBaryonDensity;
    
    int ncx = lattice->numComputationalLatticePointsX;
    int ncy = lattice->numComputationalLatticePointsY;
    int ncz = lattice->numComputationalLatticePointsRapidity;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    double t0 = hydro->initialProperTimePoint;
    
    int DIM_X = lattice->numLatticePointsX;
    int DIM_Y = lattice->numLatticePointsY;
    int DIM_ETA = lattice->numLatticePointsRapidity;
    int DIM = ncx*ncy*ncz;
    
    double DX = lattice->latticeSpacingX;
    double DY = lattice->latticeSpacingY;
    double DETA = lattice->latticeSpacingRapidity;
    double TAU = hydro->initialProperTimePoint;
    
    switch(sourceType){
            
            //==================================================
            // Initialize to vacuum
            //==================================================
            
            case 0: {
                
                printf("initialized to be vaccum \n");
                
                double ed = initialEnergyDensity;
                double rhobd = initialBaryonDensity;
                double pd = equilibriumPressure(ed, rhobd);
                
                //--------------------------------------------------
                // Initialize energy, baryon and pressure density, also
                // shear, bulk, flow velocity and baryon diffusion current
                //--------------------------------------------------
                
                printf("Initialize \\pi^\\mu\\nu to zero.\n");
                printf("Initialize \\nb^\\mu to zero.\n");

                for(int i = 2; i < nx+2; ++i) {
                    for(int j = 2; j < ny+2; ++j) {
                        for(int k = 2; k < nz+2; ++k) {
                            int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                            
                            e[s] = (PRECISION) ed;
                            rhob[s] = rhobd;
                            p[s] = pd;
                            
                            PRECISION ux = 0;
                            PRECISION uy = 0;
                            PRECISION un = 0;
                            
                            u->ux[s] = 0;
                            u->uy[s] = 0;
                            u->un[s] = 0;
                            u->ut[s] = sqrt(1+ux*ux+uy*uy+t0*t0*un*un);
                            
#ifdef PIMUNU
                            q->pitt[s] = 0;
                            q->pitx[s] = 0;
                            q->pity[s] = 0;
                            q->pitn[s] = 0;
                            q->pixx[s] = 0;
                            q->pixy[s] = 0;
                            q->pixn[s] = 0;
                            q->piyy[s] = 0;
                            q->piyn[s] = 0;
                            q->pinn[s] = 0;
#endif
#ifdef PI
                            q->Pi[s] = 0;
#endif
#ifdef VMU
                            q->nbt[s] = 0;
                            q->nbx[s] = 0;
                            q->nby[s] = 0;
                            q->nbn[s] = 0;
#endif
                        }
                    }
                }
                
                return;
            }
            
            //==================================================
            // Initialize to distributions in the Landau Frame
            //==================================================
            
            case 1: {
                
                printf("initialized with distributions in the Landau frame \n");
                
                //--------------------------------------------------
                // Solve eigenvalue problem in the Landau Frame
                //--------------------------------------------------
                
                double ttt_in, ttx_in, tty_in, ttn_in, txx_in, txy_in, txn_in, tyy_in, tyn_in, tnn_in;
                double stressTensor[10][DIM];
                
                FILE *fileIn;
                char fname[255];
                
                sprintf(fname, "%s/%s", rootDirectory, "source/initialTmunu.dat");
                
                fileIn = fopen(fname, "r");
                if (fileIn == NULL)
                {
                    printf("Couldn't open initialTmunu.dat!\n");
                }
                else
                {
                    for(int i = 2; i < DIM_X+2; ++i) {
                        for(int j = 2; j < DIM_Y+2; ++j) {
                            for(int k = 2; k < DIM_ETA+2; ++k) {
                                fscanf(fileIn, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &ttt_in, &ttx_in, &tty_in, &ttn_in, &txx_in, &txy_in, &txn_in, &tyy_in, &tyn_in, &tnn_in);
                                int is = columnMajorLinearIndex(i, j, k, DIM_X+4, DIM_Y+4);
                                stressTensor[0][is] = ttt_in;
                                stressTensor[1][is] = ttx_in;
                                stressTensor[2][is] = tty_in;
                                stressTensor[3][is] = ttn_in;
                                stressTensor[4][is] = txx_in;
                                stressTensor[5][is] = txy_in;
                                stressTensor[6][is] = txn_in;
                                stressTensor[7][is] = tyy_in;
                                stressTensor[8][is] = tyn_in;
                                stressTensor[9][is] = tnn_in;
                            }
                        }
                    }
                }
                
                float tolerance = 1.0e-3; //set energy density to tolerance if it is less than tolerance and if REGULATE is true
                
                for (int is = 0; is < DIM; is++)
                {
                    gsl_matrix *Tmunu; //T^(mu,nu) with two contravariant indices; we need to lower an index
                    //using the metric to find the eigenvectors of T^(mu)_(nu) with one contravariant and one contravariant index
                    Tmunu = gsl_matrix_alloc(4,4);
                    gsl_matrix *gmunu;
                    gmunu = gsl_matrix_alloc(4,4);
                    gsl_matrix_complex *eigen_vectors;
                    eigen_vectors = gsl_matrix_complex_alloc(4,4);
                    gsl_vector_complex *eigen_values;
                    eigen_values = gsl_vector_complex_alloc(4);
                    //set the values of the energy momentum tensor
                    gsl_matrix_set(Tmunu, 0, 0, stressTensor[0][is]); //tau,tau
                    gsl_matrix_set(Tmunu, 0, 1, stressTensor[1][is]); //tau,x
                    gsl_matrix_set(Tmunu, 0, 2, stressTensor[2][is]); //tau,y
                    gsl_matrix_set(Tmunu, 0, 3, stressTensor[3][is]); //tau,eta
                    gsl_matrix_set(Tmunu, 1, 1, stressTensor[4][is]); //x,x
                    gsl_matrix_set(Tmunu, 1, 2, stressTensor[5][is]); //x,y
                    gsl_matrix_set(Tmunu, 1, 3, stressTensor[6][is]); //x,eta
                    gsl_matrix_set(Tmunu, 2, 2, stressTensor[7][is]); //y,y
                    gsl_matrix_set(Tmunu, 2, 3, stressTensor[8][is]); //y,eta
                    gsl_matrix_set(Tmunu, 3, 3, stressTensor[9][is]); //eta,eta
                    gsl_matrix_set(Tmunu, 1, 0, stressTensor[1][is]); //x,tau
                    gsl_matrix_set(Tmunu, 2, 0, stressTensor[2][is]); //y,tau
                    gsl_matrix_set(Tmunu, 3, 0, stressTensor[3][is]); //eta,tau
                    gsl_matrix_set(Tmunu, 2, 1, stressTensor[5][is]); //y,x
                    gsl_matrix_set(Tmunu, 3, 1, stressTensor[6][is]); //eta,x
                    gsl_matrix_set(Tmunu, 3, 2, stressTensor[8][is]); //eta,y
                    
                    //set the values of the "metric"; not really the metric, but the numerical constants
                    //which are multiplied by the elements of T^(mu,nu) to get the values of T^(mu)_(nu)
                    //note factors of TAU appropriate for milne coordinates g_(mu.nu) = diag(1,-1,-1,-TAU^2)
                    gsl_matrix_set(gmunu, 0, 0, 1.0); //tau,tau
                    gsl_matrix_set(gmunu, 0, 1, -1.0); //tau,x
                    gsl_matrix_set(gmunu, 0, 2, -1.0); //tau,y
                    gsl_matrix_set(gmunu, 0, 3, -1.0*TAU*TAU); //tau,eta
                    gsl_matrix_set(gmunu, 1, 0, 1.0); //x,tau
                    gsl_matrix_set(gmunu, 1, 1, -1.0); //x,x
                    gsl_matrix_set(gmunu, 1, 2, -1.0); //x,y
                    gsl_matrix_set(gmunu, 1, 3, -1.0*TAU*TAU); //x,eta
                    gsl_matrix_set(gmunu, 2, 0, 1.0); //y,tau
                    gsl_matrix_set(gmunu, 2, 1, -1.0); //y,x
                    gsl_matrix_set(gmunu, 2, 2, -1.0); //y,y
                    gsl_matrix_set(gmunu, 2, 3, -1.0*TAU*TAU); //y,eta
                    gsl_matrix_set(gmunu, 3, 0, 1.0); //eta,tau
                    gsl_matrix_set(gmunu, 3, 1, -1.0); //eta,x
                    gsl_matrix_set(gmunu, 3, 2, -1.0); //eta,y
                    gsl_matrix_set(gmunu, 3, 3, -1.0*TAU*TAU); //eta,eta
                    //lower one index of the stress tensor; save it to the same matrix to save memory
                    gsl_matrix_mul_elements(Tmunu, gmunu); //result stored in Tmunu !this multiplies element-wise, not ordinary matrix multiplication!
                    gsl_eigen_nonsymmv_workspace *eigen_workspace;
                    eigen_workspace = gsl_eigen_nonsymmv_alloc(4);
                    gsl_eigen_nonsymmv(Tmunu, eigen_values, eigen_vectors, eigen_workspace);
                    gsl_eigen_nonsymmv_free(eigen_workspace);
                    
                    //***does this have a solution for energy density and flow at every point?
                    int eigenvalue_exists = 0;
                    for (int i = 0; i < 4; i++)
                    {
                        gsl_complex eigenvalue = gsl_vector_complex_get(eigen_values, i);
                        
                        if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
                        {
                            gsl_complex v0 = gsl_matrix_complex_get(eigen_vectors, 0 , i);
                            gsl_complex v1 = gsl_matrix_complex_get(eigen_vectors, 1 , i);
                            gsl_complex v2 = gsl_matrix_complex_get(eigen_vectors, 2 , i);
                            gsl_complex v3 = gsl_matrix_complex_get(eigen_vectors, 3 , i);
                            
                            if (GSL_IMAG(v0) == 0 && (2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0 - (GSL_REAL(v3) * GSL_REAL(v3) * (TAU * TAU - 1.0) )) > 0) //choose timelike eigenvector
                            {
                                double minkowskiLength = GSL_REAL(v0)*GSL_REAL(v0) - (GSL_REAL(v1)*GSL_REAL(v1) + GSL_REAL(v2)*GSL_REAL(v2) + TAU*TAU*GSL_REAL(v3)*GSL_REAL(v3));
                                double factor = 1.0 / sqrt(minkowskiLength);
                                
                                if (GSL_REAL(v0) < 0) factor=-factor;
                                
                                //ignore eigenvectors with gamma >~ 60
                                if ( (GSL_REAL(v0) * factor) < GAMMA_MAX)
                                {
                                    eigenvalue_exists = 1;
                                    e[is] = GSL_REAL(eigenvalue);
                                    u->ut[is] = GSL_REAL(v0) * factor;
                                    u->ux[is] = GSL_REAL(v1) * factor;
                                    u->uy[is] = GSL_REAL(v2) * factor;
                                    u->un[is] = GSL_REAL(v3) * factor;
                                }
                                
                            } // if (GSL_IMAG(v0) == 0 && (2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0 - (GSL_REAL(v3) * GSL_REAL(v3) * (TAU * TAU - 1.0) )) > 0) //choose timelike eigenvector
                        } // if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
                    } //for (int i = 0; i < 4; ...)
                    
                    if (eigenvalue_exists == 0)
                    {
                        //in dilute regions where we can't find a timelike eigenvector, set e = 0, u^t = 1, u^x=u^y=u^n=0
                        e[is] = 0.0;
                        u->ut[is] = 1.0;
                        u->ux[is] = 0.0;
                        u->uy[is] = 0.0;
                        u->un[is] = 0.0;
                    }
                } // for (int is; is < DIM; ...)
                
                
                //--------------------------------------------------
                // Initialize flow velocity, energy, pressure and baryon densities
                //--------------------------------------------------
                
                //now regulate the flow velocity by a smoothing procedure. Flow can be too large in dilute regions, cause hydro to crash...
                //this method doesnt yield a smooth profile
                
                //try scaling the flow velocity by a smooth profile which goes to zero after some finite radius
                if (REGULATE)
                {
                    printf("Regulating flow velocity profile in dilute regions \n");
                    for(int i = 2; i < DIM_X+2; ++i) {
                        double x = (i-2 - (DIM_X-1)/2.)*DX;
                        for(int j = 2; j < DIM_Y+2; ++j) {
                            double y = (j-2 - (DIM_Y-1)/2.)*DY;
                            for(int k = 2; k < DIM_ETA+2; ++k) {
                                
                                int is = columnMajorLinearIndex(i, j, k, DIM_X+4, DIM_Y+4);
                                
                                double r = sqrt(x*x + y*y);
                                //printf("r=%f\n",r);
                                
                                float R_WIDTH = 0.6;
                                float R_FLAT = 4.5;
                                float arg = (-1.0) * (r - R_FLAT) * (r - R_FLAT) / (2.0 * R_WIDTH * R_WIDTH);
                                arg = arg * THETA_FUNCTION(r - R_FLAT);
                                
                                u->ux[is] = u->ux[is] * exp(arg);
                                u->uy[is] = u->uy[is] * exp(arg);
                                u->un[is] = 0.0;
                                
                                u->ut[is] = sqrt( 1 + u->ux[is]*u->ux[is] + u->uy[is]*u->uy[is] + TAU*TAU*u->un[is]*u->un[is]);
                            }
                        }
                    }
                }
                
                for(int i = 2; i < DIM_X+2; ++i) {
                    for(int j = 2; j < DIM_Y+2; ++j) {
                        for(int k = 2; k < DIM_ETA+2; ++k) {
                            int is = columnMajorLinearIndex(i, j, k, DIM_X+4, DIM_Y+4);
                            
                            u->ux[is] = 0.0;
                            u->uy[is] = 0.0;
                            u->un[is] = 0.0;
                            u->ut[is] = sqrt( 1 + u->ux[is]*u->ux[is] + u->uy[is]*u->uy[is] + TAU*TAU*u->un[is]*u->un[is]);
                            
                            e[is] = e[is] + initialEnergyDensity;

                            p[is] = equilibriumPressureWB(e[is]);
                        }
                    }
                }
                fclose(fileIn);
                
                //--------------------------------------------------
                // Initialize shear, bulk, and baryon diffusion current
                //--------------------------------------------------
                
                printf("Initialize \\pi^\\mu\\nu to zero.\n");
                printf("Initialize \\nb^\\mu to zero.\n");
                
                int nx = lattice->numLatticePointsX;
                int ny = lattice->numLatticePointsY;
                int nz = lattice->numLatticePointsRapidity;
                for(int i = 2; i < nx+2; ++i) {
                    for(int j = 2; j < ny+2; ++j) {
                        for(int k = 2; k < nz+2; ++k) {
                            int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
#ifdef PIMUNU
                            q->pitt[s] = 0;
                            q->pitx[s] = 0;
                            q->pity[s] = 0;
                            q->pitn[s] = 0;
                            q->pixx[s] = 0;
                            q->pixy[s] = 0;
                            q->pixn[s] = 0;
                            q->piyy[s] = 0;
                            q->piyn[s] = 0;
                            q->pinn[s] = 0;
#endif
#ifdef PI
                            q->Pi[s] = 0;
#endif
#ifdef VMU
                            q->nbt[s] = 0;
                            q->nbx[s] = 0;
                            q->nby[s] = 0;
                            q->nbn[s] = 0;
#endif
                        }
                    }
                }
                return;
            }
    }
}

/**************************************************************************************************************************************************/
/* Constant initial energy density distribution
/**************************************************************************************************************************************************/
void setConstantDensityInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
    
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double t0 = hydro->initialProperTimePoint;
    
    double e0 = initCond->initialEnergyDensity;
    double rhob0 = initCond->initialBaryonDensity;

    double eL[nz];
    double rhoLa[nz], rhoLb[nz];
    
    longitudinalEnergyDensityDistribution(eL, latticeParams, initCondParams);
    longitudinalBaryonDensityDistribution(rhoLa, rhoLb, latticeParams, initCondParams);

    double eT = 5.5;
    
    //char rhobb[] = "output/kappaB.dat";
    //ofstream baryonden(rhobb);
    
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                
                e[s] = (PRECISION) e0 / t0 * (eL[k-2] * eT + 1.e-5);
                rhob[s] = (PRECISION) 1 / t0 * (0.332452 * (rhoLa[k-2] + rhoLb[k-2]) * eT + 1.e-5);// normalization factor 0.33
                p[s] = equilibriumPressure(e[s], rhob[s]);
                
/*#ifdef VMU
                double T = effectiveTemperature(e[s], rhob[s]);
                double alphaB = chemicalPotentialOverT(e[s], rhob[s]);
                double seq = equilibriumEntropy(e[s], rhob[s], p[s], T, alphaB);
                
                double kappaKinetic = baryonDiffusionCoefficient(T, rhob[s], alphaB, e[s], p[s]);
                double kappaHolography = baryonDiffusionConstant(T, alphaB*T)*T;
                double kappaAdscft = criticalBaryonDiffusionCoefficientAdscft(T, rhob[s], alphaB, e[s], p[s], seq);
                double kappaPlus = criticalBaryonDiffusionCoefficientPlus(T, rhob[s], alphaB, e[s], p[s], seq);
                
                baryonden
                << setprecision(5) << setw(10) << (i-2 - (nx-1)/2.0) * dx
                << setprecision(5) << setw(10) << (j-2 - (ny-1)/2.0) * dy
                << setprecision(5) << setw(10) << (k-2 - (nz-1)/2.0) * dz
                << setprecision(6) << setw(18) << kappaKinetic
                << setprecision(6) << setw(18) << kappaHolography
                << setprecision(6) << setw(18) << kappaAdscft
                << setprecision(6) << setw(18) << kappaPlus
                << endl;
#endif*/
			}
		}
	}
    //baryonden.close();
}

/**************************************************************************************************************************************************/
/* Bjorken expansion - Constant initial energy density distribution
/**************************************************************************************************************************************************/
void setBjorkenExpansionInitialCondition(void * latticeParams, void * initCondParams) {
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    double initialEnergyDensity = initCond->initialEnergyDensity;
    double initialBaryonDensity = initCond->initialBaryonDensity;
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double ed = initialEnergyDensity;
    double rhobd = initialBaryonDensity;
    
    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

                e[s] = (PRECISION) ed + 1.e-3;
                rhob[s] = (PRECISION) rhobd + 1.e-5;
                p[s] = equilibriumPressureWB(e[s]);
            }
        }
    }
}

/**************************************************************************************************************************************************/
/* Continuous optical glauber Glauber initial energy density distribution
/**************************************************************************************************************************************************/
void setGlauberInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
    double rhob0 = initCond->initialBaryonDensity;

	double eT[nx*ny], eL[nz];
    double rhoLa[nz], rhoLb[nz];
    double Ta[nx*ny], Tb[nx*ny];
    int wn; // number of wounded nucleons, deposisted baryon number
    
    energyDensityTransverseProfileAA(eT, nx, ny, dx, dy, initCondParams, Ta, Tb, &wn);
    longitudinalBaryonDensityDistribution(rhoLa, rhoLb, latticeParams, initCondParams);
    longitudinalEnergyDensityDistribution(eL, latticeParams, initCondParams);
    
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			double transverseDensity = eT[i-2+(j-2)*nx];
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				double energyDensityLongitudinal = eL[k-2];

				e[s] =  (PRECISION) e0 * (transverseDensity * energyDensityLongitudinal + 1.e-5);
                rhob[s] = (PRECISION) rhob0 * wn * ((rhoLa[k-2] + rhoLb[k-2]) * transverseDensity + 1.e-5);
                p[s] = equilibriumPressure(e[s], rhob[s]);
			}
		}
	}

}


/**************************************************************************************************************************************************/
/* Monte carlo Glauber initial energy density distribution
/**************************************************************************************************************************************************/
void setMCGlauberInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
    double rhob0 = initCond->initialBaryonDensity;
    
    double t0 = hydro->initialProperTimePoint;

#ifndef Avg_MCGlauber
    double eT[nx*ny], eL[nz];
    double rhoLa[nz], rhoLb[nz];
    double Ta[nx*ny], Tb[nx*ny];
    int n1, n2;
    
    monteCarloGlauberEnergyDensityTransverseProfile(eT, nx, ny, dx, dy, initCondParams, Ta, Tb, &n1, &n2);
    longitudinalEnergyDensityDistribution(eL, latticeParams, initCondParams);
    longitudinalBaryonDensityDistribution(rhoLa, rhoLb, latticeParams, initCondParams);
    
    
    /*char rhobb[] = "output/kappachun.dat";
    ofstream baryonden1(rhobb);
    
    char rhobc[] = "output/kappabrazil.dat";
    ofstream baryonden2(rhobc);*/
    
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                
                double eta = (k-2 - (nz-1)/2)*dz;
                
                e[s] = e0 / t0 * (eL[k-2] * (Ta[i-2 + (j-2)*nx] + Tb[i-2 + (j-2)*nx]) + 1.e-5);
                rhob[s] = rhob0 / t0 * ((n1 * rhoLa[k-2] * Ta[i-2+(j-2)*nx] + n2 * rhoLb[k-2] * Tb[i-2+(j-2)*nx]) + 1.e-5);
                p[s] = equilibriumPressure(e[s], rhob[s]);
                
                //double T = effectiveTemperature(e[s], rhob[s]);
                //double alphaB = chemicalPotentialOverT(e[s], rhob[s]);

                /*double kappa1 = 0.4/T * rhob[s] * (0.3333333*1/tanh(alphaB) - rhob[s]*T/(e[s]+p[s]));
                double kappa2 = baryonDiffusionConstant(T, alphaB*T)*T;
                
                baryonden1 << setprecision(5) << setw(10) << (i-1 - (nx-1)/2) * dx <<setprecision(5)  << setw(10)  << (j-1 - (ny-1)/2) * dy
                << setprecision(5) << setw(10) << (k-1 - (nz-1)/2) * dz << setprecision(6) << setw(18) << kappa1 << endl;
                
                baryonden2 << setprecision(5) << setw(10) << (i-1 - (nx-1)/2) * dx <<setprecision(5)  << setw(10)  << (j-1 - (ny-1)/2) * dy
                << setprecision(5) << setw(10) << (k-1 - (nz-1)/2) * dz << setprecision(6) << setw(18) << kappa2 << endl;*/
			}
		}
	}
    
    /*baryonden1.close();
    baryonden2.close();*/
#else
    double eT[201*201], eL[nz];
    double rhoLa[nz], rhoLb[nz];
    
    FILE *file;
    char fname[255];
    
    sprintf(fname, "%s/%s", rootDirectory, "/input/profiles/avgMCGlauberAu.dat");
    file = fopen(fname, "r");
    
    if (file == NULL)
    {
        printf("Couldn't open avgMCGlauberAu.dat!\n");
    }
    else
    {
        double x, y, z, ed;
        
        for(int i = 2; i < 201+2; ++i) {
            for(int j = 2; j < 201+2; ++j) {
                fscanf(file,"%lf\t%lf\t%lf\t%lf\n", &x,&y,&z,&ed);
                eT[i-2+(j-2)*201] = (PRECISION) ed;
            }
        }
    }
    
    longitudinalEnergyDensityDistribution(eL, latticeParams, initCondParams);
    longitudinalBaryonDensityDistribution(rhoLa, rhoLb, latticeParams, initCondParams);
    
    /*double es = 1.0 / t0 * (eL[50-2] * eT[100-2 + (100-2)*nx] + 1.e-5);
    double rhobs = 1 / t0 * (0.332452 * (rhoLa[50-2] + rhoLb[50-2]) * eT[100-2+(100-2)*nx] + 1.e-5);
    double seqs = 3.15 / t0 * (eL[50-2] * eT[100-2 + (100-2)*nx] + 1.e-5);
    printf("es=%f\t rhobs=%f\t seqs=%f\n",es,rhobs,seqs);
    printf("eTs=%f\t eLs=%f\t rhobLs=%f\n",eT[100-2 + (100-2)*nx],eL[50-2],(rhoLa[50-2] + rhoLb[50-2]));
    
    PRECISION PrimaryVariables[3];
    getPrimaryVariablesCombo(es, rhobs, PrimaryVariables);
    
    double peq = PrimaryVariables[0];
    double Teq = PrimaryVariables[1];
    double alphaBeq = PrimaryVariables[2];
    double Seq = equilibriumEntropy(es, rhobs, peq, Teq, alphaBeq);

    printf("Seq=%f \n",Seq);*/
    
    /*PRECISION PrimaryVariables[3];
    
    for(int i = 0; i < 100; ++i){
        
        e0 = es * (3.05+i*0.01);
        getPrimaryVariablesCombo(e0, rhobs, PrimaryVariables);
        
        double peq = PrimaryVariables[0];
        double Teq = PrimaryVariables[1];
        double alphaBeq = PrimaryVariables[2];
        double Seq = equilibriumEntropy(e0, rhobs, peq, Teq, alphaBeq);
        double ratio = fabs(seqs - Seq);
        printf("i=%d\t e0=%f\t Seq=%f \n",i,e0,Seq);
        
        if(ratio<0.1) {printf("i=%d",i);exit(1);}
    }*/
    
    
    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                
                e[s] = (double) e0 / t0 * (eL[k-2] * eT[i-2 + (j-2)*nx] + 1.e-5);
                rhob[s] = (double) 1 / t0 * (0.332452 * (rhoLa[k-2] + rhoLb[k-2]) * eT[i-2+(j-2)*nx] + 1.e-5);// normalization factor 0.33
                p[s] = equilibriumPressure(e[s], rhob[s]);
            }
        }
    }
    
    
    /*double Ta[nx*ny], Tb[nx*ny];
    int n1, n2;
    
    int nev = 1000;
    double fac = 1/ (double) nev;
    for(int n = 0; n < nev; ++n) {
        
        printf("nev=%d\n",n);
        
        monteCarloGlauberEnergyDensityTransverseProfile(eT, nx, ny, dx, dy, initCondParams, Ta, Tb, &n1, &n2);
        
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    
                    e[s] += (double) fac * (eT[i-2 + (j-2)*nx] + 1.e-5);
                }
            }
        }
    }
    
    char avgmc[] = "output/avgMCGlauber.dat";
    ofstream avgic(avgmc);
    
    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                 
                 avgic << setprecision(5) << setw(10) << (i-1 - (nx-1)/2) * dx <<setprecision(5)  << setw(10)  << (j-1 - (ny-1)/2) * dy
                 << setprecision(5) << setw(10) << (k-1 - (nz-1)/2) * dz << setprecision(6) << setw(18) << e[s] << endl;
            }
        }
    }
    
    avgic.close();
    
    exit(1);*/
#endif
}

/**************************************************************************************************************************************************/
/* Initial conditions for the Gubser ideal hydro test
/*		- set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
/**************************************************************************************************************************************************/
void setIdealGubserInitialCondition(void * latticeParams, void * initCondParams, const char *rootDirectory) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

    double x,y,ed,u1,u2,rhod;
    
    FILE *file;
    char fname[255];

    sprintf(fname, "%s/%s", rootDirectory, "/tests/Gubser-test/Gubser_InitialProfile_ideal.dat");
    file = fopen(fname, "r");
    
    if (file == NULL)
    {
        printf("Couldn't open Gubser_InitialProfile_ideal.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                int status = fscanf(file,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &x,&y,&ed,&rhod,&u1,&u2);
                for(int k = 2; k < nz+2; ++k) {
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    
                    e[s] = (PRECISION) ed;
                    p[s] = e[s]/3;
                    u->ux[s] = u1;
                    u->uy[s] = u2;
                    u->un[s] = 0;
                    u->ut[s] = sqrt(1 + u1*u1 + u2*u2);
                    rhob[s] = (PRECISION) rhod;
                }
            }
        }
    }
}

/**************************************************************************************************************************************************/
/* Initial conditions for the Gubser viscous hydro test
/*		- set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
/**************************************************************************************************************************************************/
void setISGubserInitialCondition(void * latticeParams, const char *rootDirectory) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double x,y,ed,u1,u2,pitt,pitx,pity,pixx,pixy,piyy,pinn,rhod,nb;

	FILE *file;
	char fname[255];
    sprintf(fname, "%s/%s", rootDirectory, "/tests/Gubser-test/Gubser_InitialProfile_IS_Baryon.dat");
	file = fopen(fname, "r");
    
    if (file == NULL)
    {
        printf("Couldn't open Gubser_InitialProfile_IS_Baryon.dat!\n");
    }
    else
    {
        double pitn=0;
        double pixn=0;
        double piyn=0;
        
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                fscanf(file,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &x,&y,&ed,&u1,&u2,&pixx,&piyy,&pixy,&pitt,&pitx,&pity,&pinn,&rhod,&nb);
                
                for(int k = 2; k < nz+2; ++k) {
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    
                    e[s] = (PRECISION) ed;
                    p[s] = e[s]/3;
                    u->ux[s] = u1;
                    u->uy[s] = u2;
                    u->un[s] = 0;
                    u->ut[s] = sqrt(1 + u1*u1 + u2*u2);
#ifdef PIMUNU
                    q->pitt[s] = (PRECISION) pitt;
                    q->pitx[s] = (PRECISION) pitx;
                    q->pity[s] = (PRECISION) pity;
                    q->pitn[s] = (PRECISION) pitn;
                    q->pixx[s] = (PRECISION) pixx;
                    q->pixy[s] = (PRECISION) pixy;
                    q->pixn[s] = (PRECISION) pixn;
                    q->piyy[s] = (PRECISION) piyy;
                    q->piyn[s] = (PRECISION) piyn;
                    q->pinn[s] = (PRECISION) pinn;
#endif
                    rhob[s] = (PRECISION) rhod;
#ifdef VMU
                    q->nbt[s] = 0;
                    q->nbx[s] = 0;
                    q->nby[s] = 0;
                    q->nbn[s] = nb;
#endif
                }
            }
        }
    }
}

/**************************************************************************************************************************************************/
/* Initial conditions to compare to MUSIC
/*        - set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
/**************************************************************************************************************************************************/
void setMusicInitialCondition(void * latticeParams, const char *rootDirectory) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double eta,ed,rhobd;
    
    FILE *file;
    char fname[255];
    sprintf(fname, "%s/%s", rootDirectory, "/tests/MUSIC-test/musictest.dat");
    file = fopen(fname, "r");
    
    if (file == NULL)
    {
        printf("Couldn't open musictest.dat!\n");
    }
    else
    {
        for(int i = 2; i < 3; ++i) {
            for(int j = 2; j < 3; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                    
                    int status = fscanf(file,"%le\t%le\t%le\n",&eta,&ed,&rhobd);
                    
                    e[s] = (PRECISION) ed/HBARC;
                    rhob[s] = (PRECISION) rhobd;
                    p[s] = equilibriumPressure(e[s], rhob[s]);
                }
            }
        }
    }


    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                int s1= columnMajorLinearIndex(2, 2, k, nx+4, ny+4);

                if(i !=2 || j!=2){
                    e[s] = e[s1];
                    rhob[s] = rhob[s1];
                    p[s] = p[s1];
                }
                
                u->ux[s] = 0;
                u->uy[s] = 0;
                u->un[s] = 0;
                u->ut[s] = 1.0;
#ifdef VMU
                q->nbt[s] = 0;
                q->nbx[s] = 0;
                q->nby[s] = 0;
                q->nbn[s] = 0;
#endif
            }
        }
    }

}

/**************************************************************************************************************************************************/
/* Initial conditions for the relativistic Sod shock-tube test
/*		- set energy density, pressure, fluid velocity u^\mu
/**************************************************************************************************************************************************/
void setSodShockTubeInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
#ifndef NBMU
#ifdef CONFORMAL_EOS
				if(x > 0) 	e[s] = (PRECISION) (0.00778147);
				else 			e[s] = (PRECISION) (0.124503);

				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;
#endif
#else
                if(x > 0)
                {
                    p[s] = (PRECISION) (0.0625);
                    rhob[s] = (PRECISION) (0.125);
                }
                else
                {
                    p[s] = (PRECISION) (1.0);
                    rhob[s] = (PRECISION) (1.0);
                }
                
                e[s] = 3.0*p[s];
                u->ux[s] = 0;
                u->uy[s] = 0;
                u->un[s] = 0;
                u->ut[s] = 1;

#endif
			}
		}
	}
}

/**************************************************************************************************************************************************/
/* Initial conditions for the relativistic Sod shock-tube test
/*		- set energy density, pressure, fluid velocity u^\mu
/**************************************************************************************************************************************************/
void set2dSodShockTubeInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                
				if(y > x) 	e[s] = (PRECISION) (0.00778147);
				else 			e[s] = (PRECISION) (0.124503);
                
				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;
			}
		}
	}
}

/**************************************************************************************************************************************************/
/* Initial conditions for the implosion in a box test
/*		- set energy density, pressure, fluid velocity u^\mu
/**************************************************************************************************************************************************/
void setImplosionBoxInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
    double initialBaryonDensity = initCond->initialBaryonDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

				e[s] = (PRECISION) (0.00778147);
				if (sqrt(x*x+y*y)<=2.5) e[s] = (PRECISION) (5.24503);

				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;

                rhob[s] = (PRECISION) initialBaryonDensity;
			}
		}
	}
}

/**************************************************************************************************************************************************/
/* Initial conditions for the implosion in a box test
/*		- set energy density, pressure, fluid velocity u^\mu
/**************************************************************************************************************************************************/
void setRayleighTaylorInstibilityInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double Lx = ( (nx-1)/2.)*dx;
	double Ly = ( (ny-1)/2.)*dy;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

				double gasGamma = 1.4;
				double gravity = 0.1;

				double rhoTop = 2.0;
				double rhoBot = 1.0;
				double pr0 = 0.01;
				double pert = 0.01;

				double yloc = 0.5 + pert*cos(M_PI*x);
				double pr;
                
				if(y > yloc)
                    pr = rhoTop*gravity*(1-y);
				else
                    pr = rhoTop*gravity*(1-yloc) + rhoBot*gravity*(yloc-y);
                
				e[s] = pr;
				p[s] = e[s]/3;
                
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;
			}
		}
	}
}

/**************************************************************************************************************************************************/
/* Initial conditions for the implosion in a box test
/*		- set energy density, pressure, fluid velocity u^\mu
/**************************************************************************************************************************************************/
void setGaussianPulseInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double Lx = ( (nx-1)/2.)*dx;
	double Ly = ( (ny-1)/2.)*dy;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

				double xc = 0.5;
				double yc = 0.5;
			    double beta = 50.0;
			    double pr = 1.0 + 1e-1*exp(-beta*((x-xc)*(x-xc)+(y-yc)*(y-yc)));

				e[s] = (PRECISION) (pr/(1.4-1.));

				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = (1+cos(2*M_PI*x/Lx))*(1+cos(2*M_PI*y/Ly));
				u->un[s] = 0;
				u->ut[s] = sqrt(1+u->uy[s]*u->uy[s]);
			}
		}
	}
}

/**************************************************************************************************************************************************/
/* Gaussian initial conditon
/**************************************************************************************************************************************************/
void GaussianProfile(void * latticeParams, void * initCondParams) {
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double e0 = initCond->initialEnergyDensity;
    
    double SIG0 = 2;
    
    for(int i = 2; i < nx+2; ++i) {
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                double x = (i - ((double)nx-1.)/2.)*dx;
                double y = (j - ((double)ny-1.)/2.)*dy;
                
                e[s] += 50*exp(-x*x/2/SIG0/SIG0-y*y/2/SIG0/SIG0);
                
                p[s] = e[s]/3;
                u->ux[s] = 0;
                u->uy[s] = 0;
                u->un[s] = 0;
                u->ut[s] = 1;
            }
        }
    }
}

/**************************************************************************************************************************************************/
/* Initial conditions for the sound propagation test
/*        - set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
/**************************************************************************************************************************************************/

void setSoundPropagationInitialCondition(void * latticeParams, void * initCondParams) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    
    double initialEnergyDensity = initCond->initialEnergyDensity;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double T0 = 3.05; // -> 0.6 GeV
    double e0 = initialEnergyDensity*pow(T0, 4);
    
    double cs = 0.57735;
    double de = e0/100.;
    double p0=e0/3;
    double lambda = (nx-1)*dx/2.;
    
    for(int i = 2; i < nx+2; ++i) {
        double x = (i-2 - (nx-1)/2.)*dx;
        double vx = cs*de/(e0+p0)*sin(2*M_PI*x/lambda);
        double ed = e0 + de*sin(2*M_PI*x/lambda);
        // periodic boundary conditions
        if (i==2) {
            vx = cs*de/(e0+p0)*sin(2*M_PI*abs(x)/lambda);
            ed = e0 + de*sin(2*M_PI*abs(x)/lambda);
        }
        
        double u0 = 1/sqrt(1-vx*vx);
        
        for(int j = 2; j < ny+2; ++j) {
            for(int k = 2; k < nz+2; ++k) {
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                e[s] = ed;
                p[s] = equilibriumPressure(e[s], 0.0);
                
                u->ux[s] = u0*vx;
                u->uy[s] = 0;
                u->un[s] = 0;
                u->ut[s] = u0;
#ifdef PIMUNU
                q->pitt[s] = 0;
                q->pitx[s] = 0;
                q->pity[s] = 0;
                q->pitn[s] = 0;
                q->pixx[s] = 0;
                q->pixy[s] = 0;
                q->pixn[s] = 0;
                q->piyy[s] = 0;
                q->piyn[s] = 0;
                q->pinn[s] = 0;
#endif
            }
        }
    }
}

/**************************************************************************************************************************************************/
/* Initial conditions to use.
/*	Set the energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny.
/* 	0 - constant energy density
/*		1 - Isreal-Stewart hydrodynamic Gubser flow test
/*		2 - Continous optical Glauber
/*		3 - Ideal hydrodynamic Gubser flow test
/*		4 - Monte carlo Glauber
/*		5 - Relativistic Sod shock-tube test
/**************************************************************************************************************************************************/
void setInitialConditions(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	int initialConditionType = initCond->initialConditionType;
	printf("Setting initial conditions: ");
	switch (initialConditionType) {
		case 0: {
			printf("Constant energy density.\n");
			setConstantDensityInitialCondition(latticeParams, initCondParams, hydroParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
            setNbmuInitialCondition(latticeParams, initCondParams, hydroParams);
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 1: {
			printf("Isreal-Stewart hydrodynamic Gubser flow test.\n");
			setISGubserInitialCondition(latticeParams, rootDirectory);
			return;
		}
		case 2: {
			printf("Continous optical Glauber.\n");
			setGlauberInitialCondition(latticeParams, initCondParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
            setNbmuInitialCondition(latticeParams, initCondParams, hydroParams);
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 3: {
			printf("Ideal hydrodynamic Gubser flow test.\n");
			setIdealGubserInitialCondition(latticeParams, initCondParams, rootDirectory);
			return;
		}
		case 4: {
			printf("Monte carlo Glauber.\n");
			setMCGlauberInitialCondition(latticeParams, initCondParams, hydroParams, rootDirectory);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
            setNbmuInitialCondition(latticeParams, initCondParams, hydroParams);
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 5: {
			printf("Relativistic Sod shock-tube test.\n");
			setSodShockTubeInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 6: {
			printf("Implosion in a box test.\n");
			setImplosionBoxInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 7: {
			printf("Rayleigh-Taylor instability test.\n");
			setRayleighTaylorInstibilityInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 8: {
			printf("Implosion in a box test.\n");
			setGaussianPulseInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 9: {
			printf("Relativistic 2d Sod shock-tube test.\n");
			set2dSodShockTubeInitialCondition(latticeParams, initCondParams);
			return;
		}
        case 10: {
            printf("Reading initial T ^mu nu from input/e.dat , input/p.dat , etc... \n");
            setInitialTmunuFromFiles(latticeParams, initCondParams, hydroParams, rootDirectory);
            return;
        }
        case 11: {
            printf("Reading initial T ^mu nu from /input/Tmunu.dat \n");
            setInitialTmunuFromFile(latticeParams, initCondParams, hydroParams, rootDirectory);
            return;
        }
        case 12:{
            printf("Gaussian Profile \n");
            GaussianProfile(latticeParams, initCondParams);
            return;
        }
        case 13:{
            printf("Hydro with dynamical sources...\n");
            setICfromSource(latticeParams, initCondParams, hydroParams, rootDirectory);
            return;
        }
        case 14:{
            printf("Test against MUSIC...\n");
            setMusicInitialCondition(latticeParams, rootDirectory);
            return;
        }
        case 15:{
            printf("Bjorken expansion...\n");
            setBjorkenExpansionInitialCondition(latticeParams, initCondParams);
            setFluidVelocityInitialCondition(latticeParams, hydroParams);
            setNbmuInitialCondition(latticeParams, initCondParams, hydroParams);
            setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
            return;
        }
		default: {
			printf("Initial condition type not defined. Exiting ...\n");
			exit(-1);
		}
	}
}
