//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <math.h> // for math functions
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

#include <H5Cpp.h>
#include <H5File.h>

using namespace std;

#define HDF5

/**************************************************************************************************************************************************/
/* ghost cells for the dynamical source terms
/**************************************************************************************************************************************************/

void setGhostCellVarsSource(DYNAMICAL_SOURCES * Source, int s, int sBC)
{
    Source->sourcet[s] = Source->sourcet[sBC];
    Source->sourcex[s] = Source->sourcex[sBC];
    Source->sourcey[s] = Source->sourcey[sBC];
    Source->sourcen[s] = Source->sourcen[sBC];
    Source->sourceb[s] = Source->sourceb[sBC];
}

void setGhostCellsSource(DYNAMICAL_SOURCES * Source, void * latticeParams)
{
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

    int nx,ny,nz,ncx,ncy,ncz;
    
    nx = lattice->numLatticePointsX;
    ny = lattice->numLatticePointsY;
    nz = lattice->numLatticePointsRapidity;
    
    ncx = lattice->numComputationalLatticePointsX;
    ncy = lattice->numComputationalLatticePointsY;
    ncz = lattice->numComputationalLatticePointsRapidity;

    int iBC,jBC,kBC;
    int s,sBC;
    
    for(int j = 2; j < ncy; ++j) {
        for(int k = 2; k < ncz; ++k) {
            iBC = 2;
            for (int i = 0; i <= 1; ++i) {
                s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                sBC = columnMajorLinearIndex(iBC, j, k, ncx, ncy);
                setGhostCellVarsSource(Source,s,sBC);
            }
            iBC = nx + 1;
            for (int i = nx + 2; i <= nx + 3; ++i) {
                s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                sBC = columnMajorLinearIndex(iBC, j, k, ncx, ncy);
                setGhostCellVarsSource(Source,s,sBC);
            }
        }
    }
    
    for(int i = 2; i < ncx; ++i) {
        for(int k = 2; k < ncz; ++k) {
            jBC = 2;
            for (int j = 0; j <= 1; ++j) {
                s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                sBC = columnMajorLinearIndex(i, jBC, k, ncx, ncy);
                setGhostCellVarsSource(Source,s,sBC);
            }
            jBC = ny + 1;
            for (int j = ny + 2; j <= ny + 3; ++j) {
                s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                sBC = columnMajorLinearIndex(i, jBC, k, ncx, ncy);
                setGhostCellVarsSource(Source,s,sBC);
            }
        }
    }
    
    for(int i = 2; i < ncx; ++i) {
        for(int j = 2; j < ncy; ++j) {
            kBC = 2;
            for (int k = 0; k <= 1; ++k) {
                s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                sBC = columnMajorLinearIndex(i, j, kBC, ncx, ncy);
                setGhostCellVarsSource(Source,s,sBC);
            }
            kBC = nz + 1;
            for (int k = nz + 2; k <= nz + 3; ++k) {
                s = columnMajorLinearIndex(i, j, k, ncx, ncy);
                sBC = columnMajorLinearIndex(i, j, kBC, ncx, ncy);
                setGhostCellVarsSource(Source,s,sBC);
            }
        }
    }
}

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

#ifndef HDF5
    FILE *sourcefile;
    char fname[255];
    //sprintf(fname, "%s/%s%d.dat", rootDirectory, "input/dynamical-source/Sources",n);
    sprintf(fname, "%s/%s%d.dat", rootDirectory, "input/DynamicalSources/Sources",n);
    sourcefile = fopen(fname, "r");
    
    double x, y, z;

    if(sourcefile==NULL){
        printf("The source file could not be opened...\n");
        exit(-1);
    }
    else
    {
      fseek(sourcefile,0L,SEEK_SET);
        
      for(int k = 2; k < nz+2; ++k){
         for(int j = 2; j < ny+2; ++j){
             for(int i = 2; i < nx+2; ++i){
               int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
               fscanf(sourcefile,"%le %le %le %le %le %le %le %le\n", &x, &y, &z, & Source->sourcet[s], & Source->sourcex[s], & Source->sourcey[s], & Source->sourcen[s], & Source->sourceb[s]);
             }
          }
       }
    }
    
    printf("The source file Sources%d.dat was read in...\n",n);

    fclose(sourcefile);
    
#else
    
	char fname[255];
	//sprintf(fname, "%s/%s%d.h5", rootDirectory, "../part2s/sourceTerms/output/Sources",n);
	sprintf(fname, "%s/%s%d.h5", rootDirectory, "input/DynamicalSources/Sources",n);
	//sourcefile = fopen(fname, "r");

	H5::H5File file( fname, H5F_ACC_RDONLY );
	H5::DataSet dataset = file.openDataSet( "data" );
	H5::DataSpace dataspace = dataset.getSpace();
	int rank = dataspace.getSimpleExtentNdims();
	hsize_t dims_out[4];
	int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
	//H5::DataSpace mspace(RANK, dims);

	//read the source terms (Sb, St, Sx, Sy, Sn) into a single array
	float *Sall;
    int nElements = nx * ny * nz;
	Sall = (float *)calloc( 5*nElements, sizeof(float) );
	dataset.read( Sall, H5::PredType::NATIVE_FLOAT, dataspace);

	//split this array into corresponding source terms
	for(int k = 2; k < nz+2; ++k){
		for(int j = 2; j < ny+2; ++j){
			for(int i = 2; i < nx+2; ++i){
                
                //this index runs over grid with ghost cells
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				//this index runs over grid without ghost cells
				int s_m = columnMajorLinearIndex(i-2, j-2, k-2, nx, ny);
				
				Source->sourcet[s] = (PRECISION) Sall[5*s_m];
				Source->sourcex[s] = (PRECISION) Sall[5*s_m + 1];
				Source->sourcey[s] = (PRECISION) Sall[5*s_m + 2];
				if (nz > 1) Source->sourcen[s] = (PRECISION) Sall[5*s_m + 3];
                
                Source->sourceb[s] = (PRECISION) Sall[5*s_m + 4];
                
                Source->sourcex[s] = 0.;
                Source->sourcey[s] = 0.;
                Source->sourcen[s] = 0.;
                
			} //for(int k = 2; k < nz+2; ++k)
		} // for(int j = 2; j < ny+2; ++j)
	} // for(int i = 2; i < nx+2; ++i)
    
    printf("The source file Sources%d.h5 was read in...\n",n);
    
#endif
    
    // ghost cells for dynamical sources
    setGhostCellsSource(Source, latticeParams);
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
    
    // ghost cells for dynamical sources
    setGhostCellsSource(Source, latticeParams);
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
	double width = 0.1; //width of gaussian smearing

	for(int i = 2; i < nx+2; ++i){
        for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                double x = (double)i * dx + xmin;
                double y = (double)j * dy + ymin;
                double z = (double)k * dz + zmin;
                smearedPosition[s] = exp((-1.0)*(pos[1] - x) * (pos[1] - x) / width)
                                   * exp((-1.0)*(pos[2] - y) * (pos[2] - y) / width)
                                   * exp((-1.0)*(pos[3] - z) * (pos[3] - z) / width);
            }//k
        }//j
    }//i

    //now multiply the smeared position by energy loss corresponding to vector components
	for(int i = 2; i < nx+2; ++i){
        for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
                Source->sourcet[s] = -dp_dtau[0] * smearedPosition[s];
                Source->sourcex[s] = -dp_dtau[1] * smearedPosition[s];
                Source->sourcey[s] = -dp_dtau[2] * smearedPosition[s];
                Source->sourcen[s] = -dp_dtau[3] * smearedPosition[s];
                Source->sourceb[s] = dummy;
            }//k
        }//j
    }//i
}

