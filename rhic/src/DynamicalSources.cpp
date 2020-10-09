//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <math.h> // for math functions
#include <cmath>
#include <stdio.h> // for printf
#include <stdlib.h> //TEMP
#include <iostream>
#include <istream>
#include <fstream>
#include <cassert>
#include <string>
#include <iomanip>

#include "../include/DynamicalVariables.h"
#include "../include/DynamicalSources.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroParameters.h"

using namespace std;

/**************************************************************************************************************************************************/
/* Initialize the dynamical source terms
/**************************************************************************************************************************************************/

void readInSource(int n, void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory)
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

    FILE *sourcefile;
    char fname[255];
    sprintf(fname, "%s/%s%d.dat", rootDirectory, "input/dynamical-source/Sources",n);
    sourcefile = fopen(fname, "r");
    
    double x, y, z;

    if(sourcefile==NULL){
        printf("The source file could not be opened...\n");
        exit(-1);
    }
    else
    {
      fseek(sourcefile,0L,SEEK_SET);
        
      for(int i = 2; i < nx+2; ++i){
         for(int j = 2; j < ny+2; ++j){
             for(int k = 2; k < nz+2; ++k){
               int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
               fscanf(sourcefile,"%le %le %le %le %le %le %le %le\n", &x, &y, &z, & Source->sourcet[s], & Source->sourcex[s], & Source->sourcey[s], & Source->sourcen[s], & Source->sourceb[s]);
             }
          }
       }
    }

    fclose(sourcefile);
}

/**************************************************************************************************************************************************/
/* When dynamical source ends, zero it
/**************************************************************************************************************************************************/

void zeroSource(void * latticeParams, void * initCondParams)
{
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    
    for(int k = 2; k < nz+2; ++k){
        for(int j = 2; j < ny+2; ++j){
            for(int i = 2; i < nx+2; ++i){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                Source->sourcet[s] = 0.0;
                Source->sourcex[s] = 0.0;
                Source->sourcey[s] = 0.0;
                Source->sourcen[s] = 0.0;
                Source->sourceb[s] = 0.0;
            }//k
        }//j
    }//i
}

/**************************************************************************************************************************************************/
/* Dynamical source terms from the jet traversing the medium
/**************************************************************************************************************************************************/

void setDynamicalSources(void * latticeParams, void * initCondParams, double *dp_dtau, double *pos) //dp_dtau is the jet energy loss, pos is the jet position
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	int ncx = lattice->numComputationalLatticePointsX;
	int ncy = lattice->numComputationalLatticePointsY;
	int ncz = lattice->numComputationalLatticePointsRapidity;
	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	PRECISION dummy = 0.0;

    double xmin = (-1.0) * ((double)(nx - 1) / 2.0) * dx;
	double ymin = (-1.0) * ((double)(ny - 1) / 2.0) * dy;
	double zmin = (-1.0) * ((double)(nz - 1) / 2.0) * dz;

	//construct an array of the gaussian smeared jet position
	double smearedPosition[ncx * ncy * ncz];
	double width = 0.2; //width of gaussian smearing
    double wid = 2 * width * width;
    double fac = 1 / width / sqrt(2 * M_PI);
    
    //printf("source term: parton has position {%f, %f, %f, %f}\n", pos[0],pos[1],pos[2],pos[3]);

	for(int k = 2; k < nz+2; ++k){
        for(int j = 2; j < ny+2; ++j){
            for(int i = 2; i < nx+2; ++i){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                double x = (double) (i-2) * dx + xmin;
                double y = (double) (j-2) * dy + ymin;
                double z = (double) (k-2) * dz + zmin;
                
                //smearedPosition[s] = exp((-1.0)*(pos[1] - x) * (pos[1] - x) / width) * exp((-1.0)*(pos[2] - y) * (pos[2] - y) / width) * exp((-1.0)*(pos[3] - z) * (pos[3] - z) / width);
                
                smearedPosition[s] = fac * exp((-1.0)*(pos[1] - x) * (pos[1] - x) / wid) * exp((-1.0)*(pos[2] - y) * (pos[2] - y) / wid);
                
            }//k
        }//j
    }//i

    //now multiply the smeared position by energy loss corresponding to vector components
	for(int k = 2; k < nz+2; ++k){
        for(int j = 2; j < ny+2; ++j){
            for(int i = 2; i < nx+2; ++i){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                Source->sourcet[s] += -dp_dtau[0] * smearedPosition[s];
                Source->sourcex[s] += -dp_dtau[1] * smearedPosition[s];
                Source->sourcey[s] += -dp_dtau[2] * smearedPosition[s];
                Source->sourcen[s] += -dp_dtau[3] * smearedPosition[s];
                Source->sourceb[s] += dummy;
            }//k
        }//j
    }//i
}

