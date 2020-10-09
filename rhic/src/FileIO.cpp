//**********************************************************************************//
//  BEShydro: A (3+1)-dimensional diffusive relativistic hydrodynamic code          //
//                                                                                  //
//          By Dennis Bazow, Lipei Du, Derek Everett and Ulrich Heinz               //
//**********************************************************************************//

#include <stdlib.h>
#include <stdio.h>

#include "../include/FileIO.h"
#include "../include/LatticeParameters.h"
#include "../include/DynamicalVariables.h"


void output(const PRECISION * const var, double t, const char *pathToOutDir, const char *name, void * latticeParams) {
	FILE *fp;
	char fname[255];
	sprintf(fname, "%s/%s_%.3f.dat", pathToOutDir, name, t);
	fp=fopen(fname, "w");

	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double x,y,z;

	int i,j,k;
	int s;

	for(k = 2; k < nz+2; ++k) {
		z = (k-2 - (nz-1)/2.)*dz;
		for(j = 2; j < ny+2; ++j) {
			y = (j-2 - (ny-1)/2.)*dy;
			for(i = 2; i < nx+2; ++i) {
				x = (i-2 - (nx-1)/2.)*dx;
				s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

				fprintf(fp, "%.3f\t%.3f\t%.3f\t%.8f\n",x,y,z,var[s]);
			}
		}
	}

	fclose(fp);
}

void outputPhaseDiagram(const PRECISION * const var1, const PRECISION * const var2, double t, const char *pathToOutDir, const char *name, void * latticeParams){
    FILE *fp;
    char fname[255];
    sprintf(fname, "%s/%s_%.3f.dat", pathToOutDir, name, t);
    fp=fopen(fname, "w");
    
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;
    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    
    double x,y,z;
    
    int i,j,k;
    int s;
    
    for(k = 2; k < nz+2; ++k) {
        z = (k-2 - (nz-1)/2.)*dz;
        for(j = 2; j < ny+2; ++j) {
            y = (j-2 - (ny-1)/2.)*dy;
            for(i = 2; i < nx+2; ++i) {
                x = (i-2 - (nx-1)/2.)*dx;
                s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);

                fprintf(fp, "%.3f\t%.3f\t%.3f\t%.8f\t%.8f\n",x,y,z,var1[s]*var2[s],var2[s]);
            }
        }
    }
    
    fclose(fp);
}
